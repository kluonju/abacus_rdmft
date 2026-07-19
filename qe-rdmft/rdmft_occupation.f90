!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_occupation
  !------------------------------------------------------------------
  !! Box constraint and electron-number equality for the occupation
  !! sub-problem of the RDMFT solver.
  !!
  !! Mirror of \texttt{rdmft\_occupation.h} in the ABACUS reference.
  !! The two main routines are:
  !!
  !! * \texttt{rdmft\_proximal\_project} -- L2 proximal projector onto
  !!   ``{ 0 <= n <= 1, sum_k w_k sum_i n_ik = N_e }`` solved via the
  !!   shifted clipping form ``y_ki(lam) = clip(x - lam * w_k, 0, 1)``
  !!   and bisection on the monotone weighted sum.
  !!
  !! * \texttt{rdmft\_proximal\_project\_occ} -- dispatches to the
  !!   global projector or to two spin-resolved projectors when
  !!   \texttt{tot\_magnetization} fixes \texttt{nelup}/\texttt{neldw}.
  !!
  !! * \texttt{rdmft\_projected\_grad\_inf} -- the Bertsekas projected
  !!   gradient norm
  !!     ``g_proj = (n - P(n - tau g)) / tau``
  !!   used as the stopping criterion of the inner projected-gradient
  !!   loop (see ABACUS rdmft_derivation.md \S6.2).
  !
  USE kinds, ONLY : DP
  USE rdmft_module, ONLY : rdmft_fix_magnetization, rdmft_n_target, &
                           rdmft_n_target_up, rdmft_n_target_down
  USE lsda_mod, ONLY : lsda, isk
  USE mp,        ONLY : mp_sum, mp_max, mp_min
  USE mp_pools,  ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_proximal_project, rdmft_proximal_project_occ, &
            rdmft_proximal_project_uniform, rdmft_proximal_project_uniform_occ, &
            rdmft_weighted_sum, rdmft_weighted_sum_ispin, &
            rdmft_magnetization_electrons, rdmft_constraint_violation_occ, &
            rdmft_occ_wdot, rdmft_occ_w_norm2, rdmft_occ_wflat_from_wk, &
            rdmft_pg_kkt_residual, rdmft_pg_map_grad_inf, &
            rdmft_projected_grad_inf, rdmft_projected_grad_inf_occ
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_weighted_sum(occ, wk, nbnd, nks)
    !---------------------------------------------------------------
    !! \(\sum_k w_k \sum_i n_{ik}\).
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: occ(nbnd, nks), wk(nks)
    INTEGER :: ib, ik
    REAL(DP) :: s
    !
    s = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          s = s + wk(ik) * occ(ib, ik)
       ENDDO
    ENDDO
    CALL mp_sum(s, inter_pool_comm)
    rdmft_weighted_sum = s
    !
  END FUNCTION rdmft_weighted_sum
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_weighted_sum_ispin(occ, wk, nbnd, nks, ispin_filter)
    !---------------------------------------------------------------
    !! \(\sum_{k,\sigma} w_k \sum_i n_{ik}\) with \(\sigma=\) \texttt{ispin\_filter}.
    INTEGER, INTENT(IN) :: nbnd, nks, ispin_filter
    REAL(DP), INTENT(IN) :: occ(nbnd, nks), wk(nks)
    INTEGER :: ib, ik
    REAL(DP) :: s
    !
    s = 0.0_DP
    DO ik = 1, nks
       IF (isk(ik) /= ispin_filter) CYCLE
       DO ib = 1, nbnd
          s = s + wk(ik) * occ(ib, ik)
       ENDDO
    ENDDO
    CALL mp_sum(s, inter_pool_comm)
    rdmft_weighted_sum_ispin = s
    !
  END FUNCTION rdmft_weighted_sum_ispin
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_magnetization_electrons(occ, wk, nbnd, nks)
    !---------------------------------------------------------------
    !! \(N_\uparrow - N_\downarrow\) in electron units (same as PWscf
    !! \texttt{tot\_magnetization} = \texttt{nelup} - \texttt{neldw}).
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: occ(nbnd, nks), wk(nks)
    !
    rdmft_magnetization_electrons = &
         rdmft_weighted_sum_ispin(occ, wk, nbnd, nks, 1) - &
         rdmft_weighted_sum_ispin(occ, wk, nbnd, nks, 2)
    !
  END FUNCTION rdmft_magnetization_electrons
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_occ_wdot(a, b, wk, nbnd, nks)
    !---------------------------------------------------------------
    !! Weighted inner product \(\langle a, b\rangle_W = \sum_k w_k
    !! \sum_i a_{ik} b_{ik}\) (BZ weights in the numerator).
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: a(nbnd, nks), b(nbnd, nks), wk(nks)
    INTEGER :: ib, ik
    REAL(DP) :: s
    !
    s = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          s = s + wk(ik) * a(ib, ik) * b(ib, ik)
       ENDDO
    ENDDO
    CALL mp_sum(s, inter_pool_comm)
    rdmft_occ_wdot = s
    !
  END FUNCTION rdmft_occ_wdot
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_occ_w_norm2(vec, wk, nbnd, nks)
    !---------------------------------------------------------------
    !! Squared weighted norm \(\|v\|_W^2 = \sum_k w_k \sum_i v_{ik}^2\).
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: vec(nbnd, nks), wk(nks)
    INTEGER :: ib, ik
    REAL(DP) :: s
    !
    s = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          s = s + wk(ik) * vec(ib, ik) * vec(ib, ik)
       ENDDO
    ENDDO
    CALL mp_sum(s, inter_pool_comm)
    rdmft_occ_w_norm2 = s
    !
  END FUNCTION rdmft_occ_w_norm2
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_occ_wflat_from_wk(wk, nbnd, nks, wflat)
    !---------------------------------------------------------------
    !! Flatten ``w_k`` to match ``RESHAPE(occ, [nbnd*nks])`` (``ib`` fastest).
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: wk(nks)
    REAL(DP), INTENT(OUT) :: wflat(:)
    INTEGER :: ib, ik, i
    !
    i = 0
    DO ik = 1, nks
       DO ib = 1, nbnd
          i = i + 1
          wflat(i) = wk(ik)
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_occ_wflat_from_wk
  !
  !-----------------------------------------------------------------
  SUBROUTINE eval_weighted_sum(x, shift, wk, nbnd, nks, ispin_filter, use_spin_filter, ssum)
    !---------------------------------------------------------------
    !! Compute \(\sum_k w_k \sum_i \mathrm{clip}(x_{ik} - \lambda w_k,
    !! 0, 1)\) for a single trial \(\lambda\) (= \texttt{shift}).
    !! When \texttt{use\_spin\_filter}, only k-rows with
    !! \texttt{isk(ik)=ispin\_filter} contribute.
    INTEGER, INTENT(IN) :: nbnd, nks, ispin_filter
    REAL(DP), INTENT(IN) :: x(nbnd, nks), wk(nks), shift
    LOGICAL, INTENT(IN) :: use_spin_filter
    REAL(DP), INTENT(OUT) :: ssum
    INTEGER :: ib, ik
    REAL(DP) :: y, sub
    !
    ssum = 0.0_DP
    DO ik = 1, nks
       IF (use_spin_filter .AND. isk(ik) /= ispin_filter) CYCLE
       sub = shift * wk(ik)
       DO ib = 1, nbnd
          y = MAX(0.0_DP, MIN(1.0_DP, x(ib, ik) - sub))
          ssum = ssum + wk(ik) * y
       ENDDO
    ENDDO
    CALL mp_sum(ssum, inter_pool_comm)
    !
  END SUBROUTINE eval_weighted_sum
  !
  !-----------------------------------------------------------------
  SUBROUTINE eval_weighted_sum_uniform(x, mu, wk, nbnd, nks, ispin_filter, use_spin_filter, ssum)
    !---------------------------------------------------------------
    !! \(\sum_k w_k \sum_i \mathrm{clip}(x_{ik} - \mu, 0, 1)\) (uniform shift).
    INTEGER, INTENT(IN) :: nbnd, nks, ispin_filter
    REAL(DP), INTENT(IN) :: x(nbnd, nks), wk(nks), mu
    LOGICAL, INTENT(IN) :: use_spin_filter
    REAL(DP), INTENT(OUT) :: ssum
    INTEGER :: ib, ik
    REAL(DP) :: y
    !
    ssum = 0.0_DP
    DO ik = 1, nks
       IF (use_spin_filter .AND. isk(ik) /= ispin_filter) CYCLE
       DO ib = 1, nbnd
          y = MAX(0.0_DP, MIN(1.0_DP, x(ib, ik) - mu))
          ssum = ssum + wk(ik) * y
       ENDDO
    ENDDO
    CALL mp_sum(ssum, inter_pool_comm)
    !
  END SUBROUTINE eval_weighted_sum_uniform
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_constraint_violation_occ(occ, wk, nbnd, nks)
    !---------------------------------------------------------------
    !! Scalar equality residual for the occupation constraint set.
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: occ(nbnd, nks), wk(nks)
    !
    IF (rdmft_fix_magnetization .AND. lsda) THEN
       rdmft_constraint_violation_occ = &
            ABS(rdmft_weighted_sum_ispin(occ, wk, nbnd, nks, 1) - rdmft_n_target_up) + &
            ABS(rdmft_weighted_sum_ispin(occ, wk, nbnd, nks, 2) - rdmft_n_target_down)
    ELSE
       rdmft_constraint_violation_occ = &
            ABS(rdmft_weighted_sum(occ, wk, nbnd, nks) - rdmft_n_target)
    ENDIF
    !
  END FUNCTION rdmft_constraint_violation_occ
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_proximal_project(occ, wk, nbnd, nks, target_ne)
    !---------------------------------------------------------------
    !! In-place L2 proximal projection of \texttt{occ} onto
    !! \({ 0 <= n <= 1, \sum_k w_k \sum_i n_{ik} = \texttt{target\_ne} }\).
    !!
    !! The dual variable \(\lambda\) is solved by bisection on the
    !! piecewise-linear, non-increasing function
    !! \(F(\lambda) = \sum_k w_k \sum_i \mathrm{clip}(x_{ik} - \lambda w_k,
    !! 0, 1)\).
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: wk(nks), target_ne
    REAL(DP), INTENT(INOUT) :: occ(nbnd, nks)
    !
    CALL rdmft_proximal_project_core(occ, wk, nbnd, nks, target_ne, &
                                     0, .FALSE.)
    !
  END SUBROUTINE rdmft_proximal_project
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_proximal_project_ispin(occ, wk, nbnd, nks, ispin_filter, target_ne)
    !---------------------------------------------------------------
    !! Like \texttt{rdmft\_proximal\_project}, but only k-rows with
    !! \texttt{isk(ik)=ispin\_filter} are updated; other rows are left
    !! unchanged.
    INTEGER, INTENT(IN) :: nbnd, nks, ispin_filter
    REAL(DP), INTENT(IN) :: wk(nks), target_ne
    REAL(DP), INTENT(INOUT) :: occ(nbnd, nks)
    !
    CALL rdmft_proximal_project_core(occ, wk, nbnd, nks, target_ne, &
                                     ispin_filter, .TRUE.)
    !
  END SUBROUTINE rdmft_proximal_project_ispin
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_proximal_project_occ(occ, wk, nbnd, nks)
    !---------------------------------------------------------------
    !! Euclidean L2 occupation projector ``P_w`` (SPG, joint init, cleanup).
    !! SPG uses ``P_w`` not uniform-shift ``P_u`` so that projection,
    !! chord slopes, and spectral steplength share the Euclidean metric
    !! required by BMR SPG2 (Lemma 2.1; see doc §5.4.1).
    !! Enforces \texttt{nelup}/\texttt{neldw} when \texttt{rdmft\_fix\_magnetization},
    !! otherwise the total \texttt{nelec} constraint only.
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: wk(nks)
    REAL(DP), INTENT(INOUT) :: occ(nbnd, nks)
    !
    IF (rdmft_fix_magnetization .AND. lsda) THEN
       CALL rdmft_proximal_project_ispin(occ, wk, nbnd, nks, 1, rdmft_n_target_up)
       CALL rdmft_proximal_project_ispin(occ, wk, nbnd, nks, 2, rdmft_n_target_down)
    ELSE
       CALL rdmft_proximal_project(occ, wk, nbnd, nks, rdmft_n_target)
    ENDIF
    !
  END SUBROUTINE rdmft_proximal_project_occ
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_proximal_project_core(occ, wk, nbnd, nks, target_ne, &
                                        ispin_filter, use_spin_filter)
    !---------------------------------------------------------------
    !! Core bisection for the Euclidean L2 projector ``P_w`` onto
    !! \(\{0\le n\le 1,\; \sum_k w_k\sum_i n_{ik}=N_e\}\):
    !! \(y_{ik}=\mathrm{clip}(x_{ik}-\lambda w_k,0,1)\) with \(\lambda\)
    !! chosen so \(\sum_k w_k\sum_i y_{ik}=\texttt{target\_ne}\).
    INTEGER, INTENT(IN) :: nbnd, nks, ispin_filter
    REAL(DP), INTENT(IN) :: wk(nks), target_ne
    LOGICAL, INTENT(IN) :: use_spin_filter
    REAL(DP), INTENT(INOUT) :: occ(nbnd, nks)
    !
    REAL(DP) :: x(nbnd, nks)
    REAL(DP) :: lam_lo, lam_hi, lam_mid, sum_lo, sum_hi, sum_mid, s0
    REAL(DP) :: shift, sub
    REAL(DP), PARAMETER :: tol_sum = 1.0e-12_DP, expand_max = 1.0e12_DP
    INTEGER :: ib, ik, it, out_of_box
    LOGICAL :: already_box
    !
    x = occ
    out_of_box = 0
    DO ik = 1, nks
       IF (use_spin_filter .AND. isk(ik) /= ispin_filter) CYCLE
       DO ib = 1, nbnd
          IF (x(ib, ik) < -1.0e-15_DP .OR. x(ib, ik) > 1.0_DP + 1.0e-15_DP) THEN
             out_of_box = 1
          ENDIF
       ENDDO
    ENDDO
    CALL mp_sum(out_of_box, inter_pool_comm)
    already_box = (out_of_box == 0)
    !
    CALL eval_weighted_sum(x, 0.0_DP, wk, nbnd, nks, ispin_filter, use_spin_filter, s0)
    IF (already_box .AND. ABS(s0 - target_ne) < tol_sum) THEN
       DO ik = 1, nks
          IF (use_spin_filter .AND. isk(ik) /= ispin_filter) CYCLE
          DO ib = 1, nbnd
             occ(ib, ik) = MAX(0.0_DP, MIN(1.0_DP, x(ib, ik)))
          ENDDO
       ENDDO
       RETURN
    ENDIF
    !
    IF (s0 > target_ne) THEN
       lam_lo = 0.0_DP
       lam_hi = 1.0_DP
       sum_lo = s0
       CALL eval_weighted_sum(x, lam_hi, wk, nbnd, nks, ispin_filter, use_spin_filter, sum_hi)
       DO WHILE (sum_hi > target_ne .AND. lam_hi < expand_max)
          lam_hi = lam_hi * 2.0_DP
          CALL eval_weighted_sum(x, lam_hi, wk, nbnd, nks, ispin_filter, use_spin_filter, sum_hi)
       ENDDO
    ELSE
       lam_hi = 0.0_DP
       lam_lo = -1.0_DP
       sum_hi = s0
       CALL eval_weighted_sum(x, lam_lo, wk, nbnd, nks, ispin_filter, use_spin_filter, sum_lo)
       DO WHILE (sum_lo < target_ne .AND. ABS(lam_lo) < expand_max)
          lam_lo = lam_lo * 2.0_DP
          CALL eval_weighted_sum(x, lam_lo, wk, nbnd, nks, ispin_filter, use_spin_filter, sum_lo)
       ENDDO
    ENDIF
    !
    IF (.NOT. (sum_lo >= target_ne .AND. sum_hi <= target_ne)) THEN
       IF (ABS(target_ne - sum_lo) < ABS(target_ne - sum_hi)) THEN
          shift = lam_lo
       ELSE
          shift = lam_hi
       ENDIF
    ELSE
       DO it = 1, 100
          lam_mid = 0.5_DP * (lam_lo + lam_hi)
          CALL eval_weighted_sum(x, lam_mid, wk, nbnd, nks, ispin_filter, use_spin_filter, sum_mid)
          IF (ABS(sum_mid - target_ne) < tol_sum) THEN
             lam_lo = lam_mid
             lam_hi = lam_mid
             EXIT
          ENDIF
          IF (sum_mid > target_ne) THEN
             lam_lo = lam_mid
          ELSE
             lam_hi = lam_mid
          ENDIF
       ENDDO
       shift = 0.5_DP * (lam_lo + lam_hi)
    ENDIF
    !
    DO ik = 1, nks
       IF (use_spin_filter .AND. isk(ik) /= ispin_filter) CYCLE
       sub = shift * wk(ik)
       DO ib = 1, nbnd
          occ(ib, ik) = MAX(0.0_DP, MIN(1.0_DP, x(ib, ik) - sub))
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_proximal_project_core
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_proximal_project_uniform(occ, wk, nbnd, nks, target_ne)
    !---------------------------------------------------------------
    !! In-place uniform-shift proximal projection:
    !! \(y_{ik} = \mathrm{clip}(x_{ik} - \mu, 0, 1)\) with \(\mu\) chosen so
    !! \(\sum_k w_k \sum_i y_{ik} = \texttt{target\_ne}\).
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: wk(nks), target_ne
    REAL(DP), INTENT(INOUT) :: occ(nbnd, nks)
    !
    CALL rdmft_proximal_project_uniform_core(occ, wk, nbnd, nks, target_ne, &
                                            0, .FALSE.)
    !
  END SUBROUTINE rdmft_proximal_project_uniform
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_proximal_project_uniform_ispin(occ, wk, nbnd, nks, ispin_filter, target_ne)
    INTEGER, INTENT(IN) :: nbnd, nks, ispin_filter
    REAL(DP), INTENT(IN) :: wk(nks), target_ne
    REAL(DP), INTENT(INOUT) :: occ(nbnd, nks)
    !
    CALL rdmft_proximal_project_uniform_core(occ, wk, nbnd, nks, target_ne, &
                                            ispin_filter, .TRUE.)
    !
  END SUBROUTINE rdmft_proximal_project_uniform_ispin
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_proximal_project_uniform_occ(occ, wk, nbnd, nks)
    !---------------------------------------------------------------
    !! Uniform proximal projector for SPG (spin-split when fixed \(M_z\)).
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: wk(nks)
    REAL(DP), INTENT(INOUT) :: occ(nbnd, nks)
    !
    IF (rdmft_fix_magnetization .AND. lsda) THEN
       CALL rdmft_proximal_project_uniform_ispin(occ, wk, nbnd, nks, 1, rdmft_n_target_up)
       CALL rdmft_proximal_project_uniform_ispin(occ, wk, nbnd, nks, 2, rdmft_n_target_down)
    ELSE
       CALL rdmft_proximal_project_uniform(occ, wk, nbnd, nks, rdmft_n_target)
    ENDIF
    !
  END SUBROUTINE rdmft_proximal_project_uniform_occ
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_proximal_project_uniform_core(occ, wk, nbnd, nks, target_ne, &
                                                 ispin_filter, use_spin_filter)
    !---------------------------------------------------------------
    !! Core bisection for uniform-shift proximal projection (no pre-clip).
    INTEGER, INTENT(IN) :: nbnd, nks, ispin_filter
    REAL(DP), INTENT(IN) :: wk(nks), target_ne
    LOGICAL, INTENT(IN) :: use_spin_filter
    REAL(DP), INTENT(INOUT) :: occ(nbnd, nks)
    !
    REAL(DP) :: x(nbnd, nks)
    REAL(DP) :: mu_lo, mu_hi, mu_mid, sum_lo, sum_hi, sum_mid, s0, mu
    REAL(DP), PARAMETER :: tol_sum = 1.0e-12_DP, expand_max = 1.0e12_DP
    INTEGER :: ib, ik, it
    !
    x = occ
    CALL eval_weighted_sum_uniform(x, 0.0_DP, wk, nbnd, nks, ispin_filter, use_spin_filter, s0)
    IF (ABS(s0 - target_ne) < tol_sum) THEN
       DO ik = 1, nks
          IF (use_spin_filter .AND. isk(ik) /= ispin_filter) CYCLE
          DO ib = 1, nbnd
             occ(ib, ik) = MAX(0.0_DP, MIN(1.0_DP, x(ib, ik)))
          ENDDO
       ENDDO
       RETURN
    ENDIF
    !
    IF (s0 > target_ne) THEN
       mu_lo = 0.0_DP
       mu_hi = 1.0_DP
       sum_lo = s0
       CALL eval_weighted_sum_uniform(x, mu_hi, wk, nbnd, nks, ispin_filter, use_spin_filter, sum_hi)
       DO WHILE (sum_hi > target_ne .AND. mu_hi < expand_max)
          mu_hi = mu_hi * 2.0_DP
          CALL eval_weighted_sum_uniform(x, mu_hi, wk, nbnd, nks, ispin_filter, use_spin_filter, sum_hi)
       ENDDO
    ELSE
       mu_hi = 0.0_DP
       mu_lo = -1.0_DP
       sum_hi = s0
       CALL eval_weighted_sum_uniform(x, mu_lo, wk, nbnd, nks, ispin_filter, use_spin_filter, sum_lo)
       DO WHILE (sum_lo < target_ne .AND. ABS(mu_lo) < expand_max)
          mu_lo = mu_lo * 2.0_DP
          CALL eval_weighted_sum_uniform(x, mu_lo, wk, nbnd, nks, ispin_filter, use_spin_filter, sum_lo)
       ENDDO
    ENDIF
    !
    IF (.NOT. (sum_lo >= target_ne .AND. sum_hi <= target_ne)) THEN
       IF (ABS(target_ne - sum_lo) < ABS(target_ne - sum_hi)) THEN
          mu = mu_lo
       ELSE
          mu = mu_hi
       ENDIF
    ELSE
       DO it = 1, 100
          mu_mid = 0.5_DP * (mu_lo + mu_hi)
          CALL eval_weighted_sum_uniform(x, mu_mid, wk, nbnd, nks, ispin_filter, use_spin_filter, sum_mid)
          IF (ABS(sum_mid - target_ne) < tol_sum) THEN
             mu_lo = mu_mid
             mu_hi = mu_mid
             EXIT
          ENDIF
          IF (sum_mid > target_ne) THEN
             mu_lo = mu_mid
          ELSE
             mu_hi = mu_mid
          ENDIF
       ENDDO
       mu = 0.5_DP * (mu_lo + mu_hi)
    ENDIF
    !
    DO ik = 1, nks
       IF (use_spin_filter .AND. isk(ik) /= ispin_filter) CYCLE
       DO ib = 1, nbnd
          occ(ib, ik) = MAX(0.0_DP, MIN(1.0_DP, x(ib, ik) - mu))
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_proximal_project_uniform_core
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_pg_kkt_residual(n, grad_n, wk, nbnd, nks, boundary_tol, resid)
    !---------------------------------------------------------------
    !! KKT residual for box \([0,1]\) + \(\sum_k w_k n_{ik} = N_e\):
    !! \(\max(0, V - W)\) with \(V=\max (\partial E/\partial n)\) on free
    !! \(\cup\{n=1\}\), \(W=\min (\partial E/\partial n)\) on free
    !! \(\cup\{n=0\}\), using physical \(g_{ik}=\texttt{grad\_n}_{ik}\)
    !! (weighted-metric proximal geometry; ``wk`` kept for call-site
    !! compatibility).
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: n(nbnd, nks), grad_n(nbnd, nks), wk(nks)
    REAL(DP), INTENT(IN), OPTIONAL :: boundary_tol
    REAL(DP), INTENT(OUT) :: resid
    REAL(DP) :: btol, V, W, e
    INTEGER :: ib, ik
    LOGICAL :: at_upper, at_lower
    !
    btol = 1.0e-10_DP
    IF (PRESENT(boundary_tol)) btol = boundary_tol
    V = -HUGE(1.0_DP)
    W = HUGE(1.0_DP)
    DO ik = 1, nks
       DO ib = 1, nbnd
          e = grad_n(ib, ik)
          at_upper = (n(ib, ik) >= 1.0_DP - btol)
          at_lower = (n(ib, ik) <= btol)
          IF (at_upper) V = MAX(V, e)
          IF (at_lower) W = MIN(W, e)
          IF (.NOT. at_upper .AND. .NOT. at_lower) THEN
             V = MAX(V, e)
             W = MIN(W, e)
          ENDIF
       ENDDO
    ENDDO
    ! Pool reduce: V/W are global maxima/minima over all k-points.
    CALL mp_max(V, inter_pool_comm)
    CALL mp_min(W, inter_pool_comm)
    IF (V < -0.5_DP * HUGE(1.0_DP) .OR. W > 0.5_DP * HUGE(1.0_DP)) THEN
       resid = 0.0_DP
    ELSE
       resid = MAX(0.0_DP, V - W)
    ENDIF
    !
  END SUBROUTINE rdmft_pg_kkt_residual
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_pg_map_grad_inf(n, grad_n, wk, nbnd, nks, alpha, gp_inf)
    !---------------------------------------------------------------
    !! \(\|n - P_w(n - \alpha\,\partial E/\partial n)\|_\infty\) (SPG2 map)
    !! with Euclidean L2 proximal projection ``P_w`` (Section 5.4).
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: n(nbnd, nks), grad_n(nbnd, nks), wk(nks)
    REAL(DP), INTENT(IN) :: alpha
    REAL(DP), INTENT(OUT) :: gp_inf
    REAL(DP) :: trial(nbnd, nks)
    INTEGER :: ib, ik
    !
    trial = n
    DO ik = 1, nks
       DO ib = 1, nbnd
          trial(ib, ik) = n(ib, ik) - alpha * grad_n(ib, ik)
       ENDDO
    ENDDO
    CALL rdmft_proximal_project_occ(trial, wk, nbnd, nks)
    gp_inf = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          gp_inf = MAX(gp_inf, ABS(n(ib, ik) - trial(ib, ik)))
       ENDDO
    ENDDO
    CALL mp_max(gp_inf, inter_pool_comm)
    !
  END SUBROUTINE rdmft_pg_map_grad_inf
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_projected_grad_inf(occ, grad, wk, nbnd, nks, &
                                       target_ne, tau, gp_inf)
    !---------------------------------------------------------------
    !! Compute the L-infinity norm of the Bertsekas projected-gradient
    !! map ``g_proj = (n - P(n - tau g)) / tau`` for a single
    !! total-electron constraint.
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: occ(nbnd, nks), grad(nbnd, nks), wk(nks)
    REAL(DP), INTENT(IN) :: target_ne, tau
    REAL(DP), INTENT(OUT) :: gp_inf
    !
    REAL(DP) :: trial(nbnd, nks)
    REAL(DP) :: tau_safe
    INTEGER :: ib, ik
    !
    tau_safe = MAX(tau, 1.0e-30_DP)
    trial = occ - tau_safe * grad
    CALL rdmft_proximal_project(trial, wk, nbnd, nks, target_ne)
    gp_inf = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          gp_inf = MAX(gp_inf, ABS((occ(ib, ik) - trial(ib, ik)) / tau_safe))
       ENDDO
    ENDDO
    CALL mp_max(gp_inf, inter_pool_comm)
    !
  END SUBROUTINE rdmft_projected_grad_inf
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_projected_grad_inf_occ(occ, grad, wk, nbnd, nks, tau, gp_inf)
    !---------------------------------------------------------------
    !! Projected-gradient norm using the same constraint set as
    !! \texttt{rdmft\_proximal\_project\_occ}.
    INTEGER, INTENT(IN) :: nbnd, nks
    REAL(DP), INTENT(IN) :: occ(nbnd, nks), grad(nbnd, nks), wk(nks)
    REAL(DP), INTENT(IN) :: tau
    REAL(DP), INTENT(OUT) :: gp_inf
    !
    REAL(DP) :: trial(nbnd, nks)
    REAL(DP) :: tau_safe
    INTEGER :: ib, ik
    !
    tau_safe = MAX(tau, 1.0e-30_DP)
    trial = occ - tau_safe * grad
    CALL rdmft_proximal_project_occ(trial, wk, nbnd, nks)
    gp_inf = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          gp_inf = MAX(gp_inf, ABS((occ(ib, ik) - trial(ib, ik)) / tau_safe))
       ENDDO
    ENDDO
    CALL mp_max(gp_inf, inter_pool_comm)
    !
  END SUBROUTINE rdmft_projected_grad_inf_occ
  !
END MODULE rdmft_occupation
