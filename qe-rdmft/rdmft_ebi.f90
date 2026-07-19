!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_ebi
  !------------------------------------------------------------------
  !! Explicit-by-implicit (EBI) occupation inner loop.
  !!
  !! Implements the EBI parameterisation of Yao, Zhang, Fang and Su
  !! (J. Phys. Chem. A \textbf{2022}, \textbf{126}, 5654--5662;
  !! \doi{10.1021/acs.jpca.2c02345}) for natural-orbital occupations
  !! in RDMFT:
  !! \[
  !!   n_{ik} = \tfrac{1}{2}\bigl(\mathrm{erf}(x_{ik}+\mu)+1\bigr),
  !! \]
  !! with unconstrained parameters \(x_{ik}\) and a scalar shift
  !! \(\mu\) (one per spin channel when \texttt{tot\_magnetization}
  !! is fixed) solved implicitly so that
  !! \(\sum_k w_k\sum_i n_{ik}=N_e\) at every step.  The analytic
  !! gradient \(\partial E/\partial x\) follows the implicit-function
  !! chain rule (Eqs.~21--22 in the 2022 paper; same structure as the
  !! logistic ``sigma-shift'' parameterisation in the ABACUS reference):
  !! \[
  !!   \frac{\partial E}{\partial x_{ik}} = s(u_{ik})
  !!     \left( g_{ik} - w_k R \right),\quad
  !!   R \equiv \frac{\sum_{jk'} w_{k'} g_{jk'} s(u_{jk'})}
  !!                 {\sum_{jk'} w_{k'} s(u_{jk'})},
  !! \]
  !! with \(s(u) = \pi^{-1/2}\,e^{-u^2}\) and \(u_{ik}=x_{ik}+\mu\).
  !!
  !! The inner optimiser is **first-order steepest descent on \(x\)
  !! (EBI@GD)**, with a quadratic-interpolating Armijo line search
  !! along the straight segment \(x(\alpha)=x_0+\alpha d\).  This is
  !! the combination Yao \emph{et al.} (2022) found to give the
  !! "lowest converged energies for different types of systems, with
  !! the lowest computational scaling'' across HF, Müller, Power, and
  !! various transition states and large molecules: their abstract
  !! concludes that EBI@GD "consistently provides the lowest
  !! converged energies'' and that EBI@NM (Newton) suffers from local
  !! minima while preconditioner-free EBI@CG / EBI@L-BFGS are
  !! ineffective without an additional Hessian-based preconditioner
  !! (cf. the 2024 follow-up arXiv:2402.03532 by the same group, which
  !! introduces a diagonal-Hessian preconditioner for joint
  !! optimisation, not provided here).
  !!
  !! Because \(n(x(\alpha))\) satisfies the electron-count constraint
  !! by construction, no proximal projection is required during the
  !! line search.
  !
  USE kinds, ONLY : DP
  USE mp,        ONLY : mp_sum, mp_max, mp_min
  USE mp_pools,  ONLY : inter_pool_comm, npool
  USE lsda_mod,  ONLY : isk, lsda
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_ebi_occ_block
  !
  REAL(DP), PARAMETER :: ebi_logit_eps = 1.0e-12_DP
  REAL(DP), PARAMETER :: ebi_inv_sqrtpi = 0.5641895835477562869480794515607725858440506293289988108620429630535684860605_DP
  !! \(1/\sqrt{\pi}\); derivative of \((\mathrm{erf}(u)+1)/2\) w.r.t.\
  !! \(u\) is \(\exp(-u^2)/\sqrt{\pi}\).
  !
  ! Shared state for the line-search evaluator.
  REAL(DP), ALLOCATABLE, SAVE :: ebi_x0(:,:), ebi_dir(:,:), ebi_trial(:,:)
  REAL(DP), ALLOCATABLE, SAVE :: ebi_grad(:,:), ebi_wk(:)
  INTEGER, SAVE :: ebi_nbnd = 0, ebi_nks = 0
  REAL(DP), SAVE :: ebi_phi0 = 0.0_DP, ebi_mu_up = 0.0_DP, ebi_mu_dw = 0.0_DP
  LOGICAL, SAVE :: ebi_use_wolfe_grad = .FALSE.
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE ebi_pool_bcast_dp(x)
    !! Pool-agree a scalar line-search / control value.
    !! Barzilai--Borwein seeds are formed from pool-local packed
    !! vectors (different \texttt{nks} per pool when \texttt{npool>1}),
    !! so take the minimum step across pools to keep a single shared
    !! trial (\texttt{mp\_min}, not \texttt{mp\_bcast}, avoids an
    !! extra Bcast/Allreduce ordering mismatch).
    REAL(DP), INTENT(INOUT) :: x
    CALL mp_min(x, inter_pool_comm)
  END SUBROUTINE ebi_pool_bcast_dp
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION ebi_occ_from_arg(arg)
    !! Map \(u=x+\mu\) to \(n=(\mathrm{erf}(u)+1)/2\in[0,1]\).
    REAL(DP), INTENT(IN) :: arg
    ebi_occ_from_arg = 0.5_DP * (ERF(arg) + 1.0_DP)
  END FUNCTION ebi_occ_from_arg
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION ebi_occ_prime(arg)
    !! \(\mathrm{d}n/\mathrm{d}x\) at \(u=\) \texttt{arg} (same w.r.t.\
    !! \(\mu\) because \(n\) depends on \(x+\mu\)).
    REAL(DP), INTENT(IN) :: arg
    ebi_occ_prime = ebi_inv_sqrtpi * EXP(-arg * arg)
  END FUNCTION ebi_occ_prime
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION ebi_erfinv(y)
    !! Inverse error function: \texttt{ERF(ebi\_erfinv(y)) = y} for
    !! \(|y|<1\).  Used to map occupations to EBI parameters
    !! \(x=\mathrm{erf}^{-1}(2n-1)\).
    REAL(DP), INTENT(IN) :: y
    REAL(DP) :: x, dy, dfdx, yclip
    INTEGER :: it
    REAL(DP), PARAMETER :: tol = 1.0e-14_DP
    !
    yclip = MAX(-1.0_DP + ebi_logit_eps, MIN(1.0_DP - ebi_logit_eps, y))
    IF (ABS(yclip) < 1.0e-30_DP) THEN
       ebi_erfinv = 0.0_DP
       RETURN
    ENDIF
    x = SIGN(SQRT(-LOG(0.5_DP * (1.0_DP - ABS(yclip)))), yclip)
    DO it = 1, 20
       dy = ERF(x) - yclip
       IF (ABS(dy) < tol) EXIT
       dfdx = 2.0_DP * ebi_inv_sqrtpi * EXP(-x * x)
       x = x - dy / dfdx
    ENDDO
    ebi_erfinv = x
  END FUNCTION ebi_erfinv
  !
  !-----------------------------------------------------------------
  SUBROUTINE ebi_weighted_occ_sum(x, mu, wk, nbnd, nks, ispin, use_filter, ssum)
    !! \(\sum_k w_k\sum_i n(x_{ik}+\mu)\) over the selected k-set.
    USE rdmft_occupation, ONLY : rdmft_weighted_sum, rdmft_weighted_sum_ispin
    REAL(DP), INTENT(IN)  :: x(nbnd, nks), wk(nks), mu
    INTEGER,  INTENT(IN)  :: nbnd, nks, ispin
    LOGICAL,  INTENT(IN)  :: use_filter
    REAL(DP), INTENT(OUT) :: ssum
    REAL(DP) :: occ(nbnd, nks)
    INTEGER :: ib, ik
    REAL(DP) :: u
    !
    DO ik = 1, nks
       DO ib = 1, nbnd
          u = x(ib, ik) + mu
          occ(ib, ik) = ebi_occ_from_arg(u)
       ENDDO
    ENDDO
    IF (use_filter) THEN
       ssum = rdmft_weighted_sum_ispin(occ, wk, nbnd, nks, ispin)
    ELSE
       ssum = rdmft_weighted_sum(occ, wk, nbnd, nks)
    ENDIF
    !
  END SUBROUTINE ebi_weighted_occ_sum
  !
  !-----------------------------------------------------------------
  SUBROUTINE ebi_weighted_occ_sum_global(xg, mu, wkg, isk_g, nbnd, nkstot, &
                                         ispin, use_filter, ssum)
    !! Global EBI occupation sum on a pool-collected k mesh (no
    !! \texttt{inter\_pool\_comm} reduction inside the bisection loop).
    REAL(DP), INTENT(IN) :: xg(nbnd, nkstot), wkg(nkstot), isk_g(nkstot)
    REAL(DP), INTENT(IN) :: mu
    INTEGER,  INTENT(IN) :: nbnd, nkstot, ispin
    LOGICAL,  INTENT(IN) :: use_filter
    REAL(DP), INTENT(OUT) :: ssum
    INTEGER :: ib, ik
    REAL(DP) :: u
    !
    ssum = 0.0_DP
    DO ik = 1, nkstot
       IF (use_filter .AND. NINT(isk_g(ik)) /= ispin) CYCLE
       DO ib = 1, nbnd
          u = xg(ib, ik) + mu
          ssum = ssum + wkg(ik) * ebi_occ_from_arg(u)
       ENDDO
    ENDDO
    !
  END SUBROUTINE ebi_weighted_occ_sum_global
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ebi_solve_mu_global(xg, wkg, isk_g, nbnd, nkstot, &
                                       target_ne, ispin, use_filter, mu)
    !! Scalar \(\mu\) solve on a pool-collected parameter mesh.
    USE io_global, ONLY : stdout, ionode
    REAL(DP), INTENT(IN)  :: xg(nbnd, nkstot), wkg(nkstot), isk_g(nkstot)
    INTEGER,  INTENT(IN)  :: nbnd, nkstot, ispin
    REAL(DP), INTENT(IN)  :: target_ne
    LOGICAL,  INTENT(IN)  :: use_filter
    REAL(DP), INTENT(OUT) :: mu
    REAL(DP), PARAMETER :: tol = 1.0e-12_DP, mu_bound = 1.0e6_DP
    INTEGER, PARAMETER :: max_bisect = 100
    REAL(DP) :: s0, s_lo, s_hi, s_mid, mu_lo, mu_hi, mu_mid
    INTEGER :: it
    !
    CALL ebi_weighted_occ_sum_global(xg, 0.0_DP, wkg, isk_g, nbnd, nkstot, &
         ispin, use_filter, s0)
    IF (ABS(s0 - target_ne) < tol) THEN
       mu = 0.0_DP
       RETURN
    ENDIF
    IF (s0 < target_ne) THEN
       mu_lo = 0.0_DP
       mu_hi = 1.0_DP
       CALL ebi_weighted_occ_sum_global(xg, mu_hi, wkg, isk_g, nbnd, nkstot, &
            ispin, use_filter, s_hi)
       DO WHILE (s_hi < target_ne .AND. mu_hi < mu_bound)
          mu_hi = 2.0_DP * mu_hi
          CALL ebi_weighted_occ_sum_global(xg, mu_hi, wkg, isk_g, nbnd, nkstot, &
               ispin, use_filter, s_hi)
       ENDDO
    ELSE
       mu_lo = -1.0_DP
       mu_hi = 0.0_DP
       CALL ebi_weighted_occ_sum_global(xg, mu_lo, wkg, isk_g, nbnd, nkstot, &
            ispin, use_filter, s_lo)
       DO WHILE (s_lo > target_ne .AND. mu_lo > -mu_bound)
          mu_lo = 2.0_DP * mu_lo
          CALL ebi_weighted_occ_sum_global(xg, mu_lo, wkg, isk_g, nbnd, nkstot, &
               ispin, use_filter, s_lo)
       ENDDO
    ENDIF
    DO it = 1, max_bisect
       mu_mid = 0.5_DP * (mu_lo + mu_hi)
       CALL ebi_weighted_occ_sum_global(xg, mu_mid, wkg, isk_g, nbnd, nkstot, &
            ispin, use_filter, s_mid)
       IF (ABS(s_mid - target_ne) < tol) THEN
          mu = mu_mid
          RETURN
       ENDIF
       IF (s_mid < target_ne) THEN
          mu_lo = mu_mid
       ELSE
          mu_hi = mu_mid
       ENDIF
    ENDDO
    mu = 0.5_DP * (mu_lo + mu_hi)
    IF (ionode) WRITE(stdout, '(7X,A)') &
         'ebi block: WARNING rdmft_ebi_solve_mu did not fully converge; using last mu.'
  END SUBROUTINE rdmft_ebi_solve_mu_global
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ebi_solve_mu(x, wk, nbnd, nks, target_ne, ispin, use_filter, mu)
    !! Bisection on \(\mu\) so the weighted EBI occupation sum matches
    !! \texttt{target\_ne}.
    !!
    !! With \texttt{npool>1}, gather \(x\) and \(w_k\) on the full
    !! \texttt{nkstot} mesh first and run the bracket/bisection loop
    !! without per-trial \texttt{inter\_pool\_comm} reductions.  Every
    !! pool then executes the same scalar logic and obtains the same
    !! \(\mu\), avoiding MPI ordering deadlocks when \texttt{nproc\_pool=1}.
    USE io_global, ONLY : stdout, ionode
    USE klist,     ONLY : nkstot
    !
    REAL(DP), INTENT(IN)  :: x(nbnd, nks), wk(nks)
    INTEGER,  INTENT(IN)  :: nbnd, nks, ispin
    REAL(DP), INTENT(IN)  :: target_ne
    LOGICAL,  INTENT(IN)  :: use_filter
    REAL(DP), INTENT(OUT) :: mu
    !
    REAL(DP), PARAMETER :: tol = 1.0e-12_DP, mu_bound = 1.0e6_DP
    INTEGER, PARAMETER :: max_bisect = 100
    REAL(DP) :: s0, s_lo, s_hi, s_mid, mu_lo, mu_hi, mu_mid
    INTEGER :: it
    REAL(DP), ALLOCATABLE :: xg(:,:), wkg(:), isk_loc(:,:), isk_g(:,:), wk_loc(:,:), wk_g(:,:)
    INTEGER :: ik
    !
    IF (npool > 1) THEN
       ALLOCATE(xg(nbnd, nkstot), wkg(nkstot), isk_loc(1, nks), isk_g(1, nkstot))
       ALLOCATE(wk_loc(1, nks), wk_g(1, nkstot))
       DO ik = 1, nks
          isk_loc(1, ik) = REAL(isk(ik), DP)
          wk_loc(1, ik)  = wk(ik)
       ENDDO
       CALL poolcollect(nbnd, nks, x, nkstot, xg)
       CALL poolcollect(1, nks, isk_loc, nkstot, isk_g)
       CALL poolcollect(1, nks, wk_loc, nkstot, wk_g)
       DO ik = 1, nkstot
          wkg(ik) = wk_g(1, ik)
       ENDDO
       CALL rdmft_ebi_solve_mu_global(xg, wkg, isk_g(1, :), nbnd, nkstot, &
            target_ne, ispin, use_filter, mu)
       DEALLOCATE(xg, wkg, isk_loc, isk_g, wk_loc, wk_g)
       RETURN
    ENDIF
    !
    CALL ebi_weighted_occ_sum(x, 0.0_DP, wk, nbnd, nks, ispin, use_filter, s0)
    IF (ABS(s0 - target_ne) < tol) THEN
       mu = 0.0_DP
       RETURN
    ENDIF
    !
    IF (s0 < target_ne) THEN
       mu_lo = 0.0_DP
       mu_hi = 1.0_DP
       CALL ebi_weighted_occ_sum(x, mu_hi, wk, nbnd, nks, ispin, use_filter, s_hi)
       DO WHILE (s_hi < target_ne .AND. mu_hi < mu_bound)
          mu_hi = 2.0_DP * mu_hi
          CALL ebi_weighted_occ_sum(x, mu_hi, wk, nbnd, nks, ispin, use_filter, s_hi)
       ENDDO
    ELSE
       mu_lo = -1.0_DP
       mu_hi = 0.0_DP
       CALL ebi_weighted_occ_sum(x, mu_lo, wk, nbnd, nks, ispin, use_filter, s_lo)
       DO WHILE (s_lo > target_ne .AND. mu_lo > -mu_bound)
          mu_lo = 2.0_DP * mu_lo
          CALL ebi_weighted_occ_sum(x, mu_lo, wk, nbnd, nks, ispin, use_filter, s_lo)
       ENDDO
    ENDIF
    !
    DO it = 1, max_bisect
       mu_mid = 0.5_DP * (mu_lo + mu_hi)
       CALL ebi_weighted_occ_sum(x, mu_mid, wk, nbnd, nks, ispin, use_filter, s_mid)
       IF (ABS(s_mid - target_ne) < tol) THEN
          mu = mu_mid
          RETURN
       ENDIF
       IF (s_mid < target_ne) THEN
          mu_lo = mu_mid
       ELSE
          mu_hi = mu_mid
       ENDIF
    ENDDO
    mu = 0.5_DP * (mu_lo + mu_hi)
    IF (ionode) WRITE(stdout, '(7X,A)') &
         'ebi block: WARNING rdmft_ebi_solve_mu did not fully converge; using last mu.'
    !
  END SUBROUTINE rdmft_ebi_solve_mu
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ebi_params_to_occ(x, mu_up, mu_dw, wk, nbnd, nks, occ)
    !! Map EBI parameters to occupations using per-spin shifts when
    !! \texttt{rdmft\_fix\_magnetization} is active.
    USE rdmft_module, ONLY : rdmft_fix_magnetization
    !
    REAL(DP), INTENT(IN)  :: x(nbnd, nks), wk(nks), mu_up, mu_dw
    INTEGER,  INTENT(IN)  :: nbnd, nks
    REAL(DP), INTENT(OUT) :: occ(nbnd, nks)
    INTEGER :: ib, ik
    REAL(DP) :: mu
    LOGICAL :: split_spin
    !
    split_spin = (rdmft_fix_magnetization .AND. lsda)
    DO ik = 1, nks
       IF (split_spin) THEN
          mu = MERGE(mu_up, mu_dw, isk(ik) == 1)
       ELSE
          mu = mu_up
       ENDIF
       DO ib = 1, nbnd
          occ(ib, ik) = ebi_occ_from_arg(x(ib, ik) + mu)
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_ebi_params_to_occ
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ebi_occ_to_params(occ, x)
    !! Map \(n=(\mathrm{erf}(x)+1)/2\) back to \(x=\mathrm{erf}^{-1}(2n-1)\).
    REAL(DP), INTENT(IN)  :: occ(:,:)
    REAL(DP), INTENT(OUT) :: x(:,:)
    INTEGER :: ib, ik
    REAL(DP) :: n_clip
    !
    DO ik = 1, SIZE(occ, 2)
       DO ib = 1, SIZE(occ, 1)
          n_clip = MIN(MAX(ebi_logit_eps, occ(ib, ik)), 1.0_DP - ebi_logit_eps)
          x(ib, ik) = ebi_erfinv(2.0_DP * n_clip - 1.0_DP)
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_ebi_occ_to_params
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ebi_fold_mu(x, wk, nbnd, nks, mu_up, mu_dw)
    !! Canonicalise parameters by folding the implicit shift into
    !! \(x\leftarrow x+\mu\) so the next \(\mu\) solve starts near zero.
    USE rdmft_module, ONLY : rdmft_fix_magnetization, rdmft_n_target, &
                             rdmft_n_target_up, rdmft_n_target_down
    USE lsda_mod,     ONLY : lsda, isk
    !
    REAL(DP), INTENT(INOUT) :: x(nbnd, nks)
    REAL(DP), INTENT(IN)    :: wk(nks)
    INTEGER,  INTENT(IN)    :: nbnd, nks
    REAL(DP), INTENT(OUT)   :: mu_up, mu_dw
    LOGICAL :: split_spin
    INTEGER :: ik
    !
    split_spin = (rdmft_fix_magnetization .AND. lsda)
    IF (split_spin) THEN
       CALL rdmft_ebi_solve_mu(x, wk, nbnd, nks, rdmft_n_target_up, 1, .TRUE., mu_up)
       CALL rdmft_ebi_solve_mu(x, wk, nbnd, nks, rdmft_n_target_down, 2, .TRUE., mu_dw)
       DO ik = 1, nks
          IF (isk(ik) == 1) THEN
             x(:, ik) = x(:, ik) + mu_up
          ELSE
             x(:, ik) = x(:, ik) + mu_dw
          ENDIF
       ENDDO
       mu_up = 0.0_DP
       mu_dw = 0.0_DP
    ELSE
       CALL rdmft_ebi_solve_mu(x, wk, nbnd, nks, rdmft_n_target, 0, .FALSE., mu_up)
       x = x + mu_up
       mu_up = 0.0_DP
       mu_dw = 0.0_DP
    ENDIF
    !
  END SUBROUTINE rdmft_ebi_fold_mu
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ebi_transform_gradient(grad_n, x, mu_up, mu_dw, wk, nbnd, nks, grad_x)
    !! Implicit chain rule for \(\partial E/\partial x\) (EBI / sigma-shift).
    USE rdmft_module, ONLY : rdmft_fix_magnetization
    !
    REAL(DP), INTENT(IN)  :: grad_n(nbnd, nks), x(nbnd, nks), wk(nks)
    REAL(DP), INTENT(IN)  :: mu_up, mu_dw
    INTEGER,  INTENT(IN)  :: nbnd, nks
    REAL(DP), INTENT(OUT) :: grad_x(nbnd, nks)
    !
    REAL(DP), PARAMETER :: sum_s_floor = 1.0e-30_DP
    REAL(DP) :: sum_s, sum_hsp, ratio, mu, sp, wk_ik
    INTEGER :: ib, ik
    LOGICAL :: split_spin
    !
    split_spin = (rdmft_fix_magnetization .AND. lsda)
    grad_x = 0.0_DP
    !
    IF (.NOT. split_spin) THEN
       sum_s = 0.0_DP
       sum_hsp = 0.0_DP
       DO ik = 1, nks
          wk_ik = wk(ik)
          DO ib = 1, nbnd
             sp = ebi_occ_prime(x(ib, ik) + mu_up)
             sum_s = sum_s + wk_ik * sp
             sum_hsp = sum_hsp + grad_n(ib, ik) * sp
          ENDDO
       ENDDO
       CALL mp_sum(sum_s, inter_pool_comm)
       CALL mp_sum(sum_hsp, inter_pool_comm)
       IF (ABS(sum_s) < sum_s_floor) RETURN
       ratio = sum_hsp / sum_s
       DO ik = 1, nks
          wk_ik = wk(ik)
          DO ib = 1, nbnd
             sp = ebi_occ_prime(x(ib, ik) + mu_up)
             grad_x(ib, ik) = sp * (grad_n(ib, ik) - wk_ik * ratio)
          ENDDO
       ENDDO
       RETURN
    ENDIF
    !
    ! Spin-resolved: separate ratios per channel.
    sum_s = 0.0_DP
    sum_hsp = 0.0_DP
    DO ik = 1, nks
       IF (isk(ik) /= 1) CYCLE
       wk_ik = wk(ik)
       DO ib = 1, nbnd
          sp = ebi_occ_prime(x(ib, ik) + mu_up)
          sum_s = sum_s + wk_ik * sp
          sum_hsp = sum_hsp + grad_n(ib, ik) * sp
       ENDDO
    ENDDO
    CALL mp_sum(sum_s, inter_pool_comm)
    CALL mp_sum(sum_hsp, inter_pool_comm)
    IF (ABS(sum_s) >= sum_s_floor) THEN
       ratio = sum_hsp / sum_s
       DO ik = 1, nks
          IF (isk(ik) /= 1) CYCLE
          wk_ik = wk(ik)
          DO ib = 1, nbnd
             sp = ebi_occ_prime(x(ib, ik) + mu_up)
             grad_x(ib, ik) = sp * (grad_n(ib, ik) - wk_ik * ratio)
          ENDDO
       ENDDO
    ENDIF
    !
    sum_s = 0.0_DP
    sum_hsp = 0.0_DP
    DO ik = 1, nks
       IF (isk(ik) /= 2) CYCLE
       wk_ik = wk(ik)
       DO ib = 1, nbnd
          sp = ebi_occ_prime(x(ib, ik) + mu_dw)
          sum_s = sum_s + wk_ik * sp
          sum_hsp = sum_hsp + grad_n(ib, ik) * sp
       ENDDO
    ENDDO
    CALL mp_sum(sum_s, inter_pool_comm)
    CALL mp_sum(sum_hsp, inter_pool_comm)
    IF (ABS(sum_s) >= sum_s_floor) THEN
       ratio = sum_hsp / sum_s
       DO ik = 1, nks
          IF (isk(ik) /= 2) CYCLE
          wk_ik = wk(ik)
          DO ib = 1, nbnd
             sp = ebi_occ_prime(x(ib, ik) + mu_dw)
             grad_x(ib, ik) = sp * (grad_n(ib, ik) - wk_ik * ratio)
          ENDDO
       ENDDO
    ENDIF
    !
  END SUBROUTINE rdmft_ebi_transform_gradient
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ebi_sync_n_from_x(x, wk, nbnd, nks)
    !! Solve \(\mu\), map \(x\to n\), and refresh EXX occupation hooks.
    USE rdmft_module, ONLY : rdmft_n, rdmft_n_target, rdmft_n_target_up, &
                             rdmft_n_target_down, rdmft_fix_magnetization
    USE rdmft_energy, ONLY : rdmft_set_wg_from_n, rdmft_refresh_x_occupation
    USE lsda_mod,     ONLY : lsda
    !
    REAL(DP), INTENT(IN) :: x(nbnd, nks), wk(nks)
    INTEGER,  INTENT(IN) :: nbnd, nks
    REAL(DP) :: mu_up, mu_dw
    !
    IF (rdmft_fix_magnetization .AND. lsda) THEN
       CALL rdmft_ebi_solve_mu(x, wk, nbnd, nks, rdmft_n_target_up, 1, .TRUE., mu_up)
       CALL rdmft_ebi_solve_mu(x, wk, nbnd, nks, rdmft_n_target_down, 2, .TRUE., mu_dw)
    ELSE
       CALL rdmft_ebi_solve_mu(x, wk, nbnd, nks, rdmft_n_target, 0, .FALSE., mu_up)
       mu_dw = 0.0_DP
    ENDIF
    CALL rdmft_ebi_params_to_occ(x, mu_up, mu_dw, wk, nbnd, nks, rdmft_n)
    CALL rdmft_set_wg_from_n(1)
    CALL rdmft_refresh_x_occupation()
    !
  END SUBROUTINE rdmft_ebi_sync_n_from_x
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ebi_ls_eval(alpha, f, g, ierr)
    !! Line-search evaluator: \(x(\alpha)=x_0+\alpha d\), implicit
    !! \(\mu\), then total energy and directional derivative
    !! \(\phi'(\alpha)=\langle\nabla_x E(x(\alpha)), d\rangle\).
    USE rdmft_module, ONLY : rdmft_fix_magnetization, rdmft_n_target, &
                             rdmft_n_target_up, rdmft_n_target_down, &
                             rdmft_n, rdmft_occ_tol
    USE rdmft_energy, ONLY : rdmft_grad_n, rdmft_total_energy
    USE rdmft_occupation, ONLY : rdmft_constraint_violation_occ
    USE lsda_mod,     ONLY : lsda
    !
    REAL(DP), INTENT(IN)  :: alpha
    REAL(DP), INTENT(OUT) :: f, g
    INTEGER, INTENT(OUT)  :: ierr
    REAL(DP) :: phi_prime, mu_up, mu_dw, grad_x(ebi_nbnd, ebi_nks), cviol
    INTEGER  :: ib, ik
    !
    ierr = 0
    DO ik = 1, ebi_nks
       DO ib = 1, ebi_nbnd
          ebi_trial(ib, ik) = ebi_x0(ib, ik) + alpha * ebi_dir(ib, ik)
       ENDDO
    ENDDO
    CALL rdmft_ebi_sync_n_from_x(ebi_trial, ebi_wk, ebi_nbnd, ebi_nks)
    cviol = rdmft_constraint_violation_occ(rdmft_n, ebi_wk, ebi_nbnd, ebi_nks)
    IF (cviol > MAX(10.0_DP * rdmft_occ_tol, 1.0e-10_DP)) THEN
       ierr = 1
       f = HUGE(1.0_DP)
       g = 0.0_DP
       RETURN
    ENDIF
    !
    IF (ebi_use_wolfe_grad) THEN
       CALL rdmft_grad_n(ebi_grad)
       IF (rdmft_fix_magnetization .AND. lsda) THEN
          CALL rdmft_ebi_solve_mu(ebi_trial, ebi_wk, ebi_nbnd, ebi_nks, &
               rdmft_n_target_up, 1, .TRUE., mu_up)
          CALL rdmft_ebi_solve_mu(ebi_trial, ebi_wk, ebi_nbnd, ebi_nks, &
               rdmft_n_target_down, 2, .TRUE., mu_dw)
       ELSE
          CALL rdmft_ebi_solve_mu(ebi_trial, ebi_wk, ebi_nbnd, ebi_nks, &
               rdmft_n_target, 0, .FALSE., mu_up)
          mu_dw = 0.0_DP
       ENDIF
       CALL rdmft_ebi_transform_gradient(ebi_grad, ebi_trial, mu_up, mu_dw, &
            ebi_wk, ebi_nbnd, ebi_nks, grad_x)
       phi_prime = 0.0_DP
       DO ik = 1, ebi_nks
          DO ib = 1, ebi_nbnd
             phi_prime = phi_prime + grad_x(ib, ik) * ebi_dir(ib, ik)
          ENDDO
       ENDDO
       CALL mp_sum(phi_prime, inter_pool_comm)
    ELSE
       phi_prime = ebi_phi0
    ENDIF
    !
    CALL rdmft_total_energy(f)
    g = phi_prime
    !
  END SUBROUTINE rdmft_ebi_ls_eval
  !
  SUBROUTINE rdmft_ebi_ls_post_reject()
    !! Restore the EBI line origin after a rejected trial.
    USE rdmft_module, ONLY : rdmft_invalidate_xi_cache
    !
    CALL rdmft_ebi_sync_n_from_x(ebi_x0, ebi_wk, ebi_nbnd, ebi_nks)
    CALL rdmft_invalidate_xi_cache()
    !
  END SUBROUTINE rdmft_ebi_ls_post_reject
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ebi_occ_block(etot_now, occ_converged, max_inner, log_tag)
    !---------------------------------------------------------------
    !! EBI@GD occupation inner block with energy line search.
    USE io_global,  ONLY : stdout, ionode
    USE wvfct,      ONLY : nbnd
    USE klist,      ONLY : nks, wk
    USE lsda_mod,   ONLY : lsda
    USE rdmft_module
    USE rdmft_energy, ONLY : rdmft_grad_n, rdmft_total_energy
    USE rdmft_occupation, ONLY : rdmft_weighted_sum, &
                                  rdmft_constraint_violation_occ, &
                                  rdmft_pg_kkt_residual
    USE rdmft_linesearch, ONLY : rdmft_zhang_hager_state, rdmft_zhang_hager_init, &
                                  rdmft_zhang_hager_update, rdmft_zhang_hager_ref, &
                                  rdmft_ls_result, rdmft_strong_wolfe_ls, &
                                  rdmft_ls_is_wolfe, rdmft_ls_is_nm, &
                                  rdmft_resolve_occ_ls_type, rdmft_armijo_ls, &
                                  rdmft_occ_ls_alpha0, rdmft_bb_state, rdmft_bb_init, &
                                  rdmft_bb_record, rdmft_quad_ls_history, rdmft_quad_ls_reset
    USE rdmft_inner_log,   ONLY : rdmft_inner_log_record
    USE rdmft_occ_eps_log, ONLY : rdmft_occ_eps_log_record
    USE mp,         ONLY : mp_sum
    USE mp_pools,   ONLY : inter_pool_comm
    !
    REAL(DP), INTENT(INOUT) :: etot_now
    LOGICAL, INTENT(OUT)    :: occ_converged
    INTEGER, INTENT(IN)     :: max_inner
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: log_tag
    !
    REAL(DP), ALLOCATABLE :: x(:,:), x_save(:,:), grad_n(:,:), grad_x(:,:), dir(:,:)
    REAL(DP), ALLOCATABLE :: n_save(:,:), x_flat(:), grad_x_flat(:)
    REAL(DP) :: mu_up, mu_dw, phi0, E0, alpha0, alpha_acc, sum_dn
    REAL(DP) :: ref_energy, ne_now, cviol, gnorm, kkt_resid, kkt_post
    INTEGER  :: inner, ik, ib
    LOGICAL  :: ls_ok, use_wolfe, use_nm, grad_ready
    CHARACTER(LEN=24) :: ls_eff
    CHARACTER(LEN=32) :: tag
    TYPE(rdmft_ls_result) :: lsres
    TYPE(rdmft_zhang_hager_state) :: zh_occ
    TYPE(rdmft_bb_state) :: bb
    TYPE(rdmft_quad_ls_history) :: qhist
    INTEGER :: ndim
    !
    occ_converged = .TRUE.
    IF (max_inner <= 0) RETURN
    !
    tag = 'occ'
    IF (PRESENT(log_tag)) tag = TRIM(log_tag)
    !
    ALLOCATE(x(nbnd, nks), x_save(nbnd, nks), grad_n(nbnd, nks), &
         grad_x(nbnd, nks), dir(nbnd, nks), n_save(nbnd, nks))
    ndim = nbnd * nks
    ALLOCATE(x_flat(ndim), grad_x_flat(ndim))
    ebi_nbnd = nbnd
    ebi_nks  = nks
    IF (.NOT. ALLOCATED(ebi_x0))    ALLOCATE(ebi_x0(nbnd, nks))
    IF (.NOT. ALLOCATED(ebi_dir))   ALLOCATE(ebi_dir(nbnd, nks))
    IF (.NOT. ALLOCATED(ebi_trial)) ALLOCATE(ebi_trial(nbnd, nks))
    IF (.NOT. ALLOCATED(ebi_grad))  ALLOCATE(ebi_grad(nbnd, nks))
    IF (.NOT. ALLOCATED(ebi_wk))    ALLOCATE(ebi_wk(nks))
    ebi_wk = wk
    !
    CALL rdmft_resolve_occ_ls_type(rdmft_occ_ls_type, ls_eff)
    use_wolfe = rdmft_ls_is_wolfe(ls_eff)
    use_nm    = rdmft_ls_is_nm(ls_eff)
    ebi_use_wolfe_grad = use_wolfe
    CALL rdmft_zhang_hager_init(zh_occ, etot_now, rdmft_zhang_hager_eta)
    CALL rdmft_bb_init(bb, ndim, rdmft_bb_alpha_min, rdmft_bb_alpha_max)
    CALL rdmft_quad_ls_reset(qhist)
    !
    IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,A,A,I0,A)') &
         TRIM(tag)//' block (EBI@GD / ', TRIM(ls_eff), ', max ', max_inner, ' inner iters)'
    !
    CALL rdmft_ebi_occ_to_params(rdmft_n, x)
    CALL rdmft_ebi_fold_mu(x, wk, nbnd, nks, mu_up, mu_dw)
    CALL rdmft_ebi_sync_n_from_x(x, wk, nbnd, nks)
    CALL rdmft_total_energy(etot_now)
    !
    occ_converged = .FALSE.
    grad_ready = .FALSE.
    DO inner = 1, max_inner
       ! Reuse the post-step gradient computed at the end of the
       ! previous inner iteration (same n, orbitals, and density) --
       ! saves one full occupation-gradient evaluation per inner
       ! iteration.  The skip decision is uniform across MPI ranks.
       IF (grad_ready) THEN
          grad_ready = .FALSE.
          kkt_resid = kkt_post
       ELSE
          CALL rdmft_grad_n(grad_n)
          CALL rdmft_pg_kkt_residual(rdmft_n, grad_n, wk, nbnd, nks, resid=kkt_resid)
       ENDIF
       IF (rdmft_verbose >= 2) THEN
          ne_now = rdmft_weighted_sum(rdmft_n, wk, nbnd, nks)
          IF (ionode) WRITE(stdout, &
               '(9X,A,I3,A,1PE10.2,0P,A,F12.6)') &
               '['//TRIM(tag)//' inner ', inner, '] KKT_resid=', kkt_resid, &
               ' Ne=', ne_now
       ENDIF
       !
       IF (rdmft_fix_magnetization .AND. lsda) THEN
          CALL rdmft_ebi_solve_mu(x, wk, nbnd, nks, rdmft_n_target_up, 1, .TRUE., mu_up)
          CALL rdmft_ebi_solve_mu(x, wk, nbnd, nks, rdmft_n_target_down, 2, .TRUE., mu_dw)
       ELSE
          CALL rdmft_ebi_solve_mu(x, wk, nbnd, nks, rdmft_n_target, 0, .FALSE., mu_up)
          mu_dw = 0.0_DP
       ENDIF
       CALL rdmft_ebi_transform_gradient(grad_n, x, mu_up, mu_dw, wk, nbnd, nks, &
            grad_x)
       !
       ! Steepest descent in x-space (EBI@GD).
       dir = -grad_x
       gnorm = 0.0_DP
       DO ik = 1, nks
          DO ib = 1, nbnd
             gnorm = gnorm + grad_x(ib, ik) * grad_x(ib, ik)
          ENDDO
       ENDDO
       CALL mp_sum(gnorm, inter_pool_comm)
       gnorm = SQRT(gnorm)
       ! EBI uses erf(x+mu) with mu enforcing Ne; box KKT (rdmft_pg_kkt_residual)
       ! is not a meaningful convergence test here — use ||grad_x|| like GD uses phi0.
       IF ((rdmft_occ_grad_tol > 0.0_DP .AND. gnorm <= rdmft_occ_grad_tol) &
            .OR. gnorm < 1.0e-30_DP) THEN
          occ_converged = .TRUE.
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, &
               '(7X,A,1PE10.2,A,1PE10.2,A)') &
               TRIM(tag)//' block converged: |grad_x|=', gnorm, &
               ' (stationary).  KKT=', kkt_resid
          EXIT
       ENDIF
       !
       x_save  = x
       n_save  = rdmft_n
       ebi_x0  = x
       ebi_dir = dir
       ebi_mu_up = mu_up
       ebi_mu_dw = mu_dw
       E0 = etot_now
       !
       x_flat = RESHAPE(x, [ndim])
       grad_x_flat = RESHAPE(grad_x, [ndim])
       alpha0 = rdmft_occ_ls_alpha0(rdmft_occ_ls_init_step, x_flat, grad_x_flat, &
            rdmft_occ_ls_stepsize, bb, qhist, rdmft_bb_alpha_min, rdmft_bb_alpha_max)
       IF (alpha0 <= 0.0_DP) alpha0 = MAX(rdmft_occ_ls_stepsize, rdmft_bgd_tau)
       CALL ebi_pool_bcast_dp(alpha0)
       CALL rdmft_bb_record(bb, x_flat, grad_x_flat)
       !
       ! phi'(0) = <grad_x, dir> = -||grad_x||^2 <= 0.
       phi0 = 0.0_DP
       DO ik = 1, nks
          DO ib = 1, nbnd
             phi0 = phi0 + grad_x(ib, ik) * dir(ib, ik)
          ENDDO
       ENDDO
       CALL mp_sum(phi0, inter_pool_comm)
       ebi_phi0 = phi0
       !
       IF (phi0 >= 0.0_DP) THEN
          occ_converged = .TRUE.
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A)') &
               TRIM(tag)//' block converged: phi0=', phi0, ' >= 0 (no descent).'
          EXIT
       ENDIF
       !
       IF (rdmft_occ_tol > 0.0_DP .AND. ABS(phi0) <= rdmft_occ_tol) THEN
          occ_converged = .TRUE.
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A)') &
               TRIM(tag)//' block converged: |phi0|=', ABS(phi0), &
               ' <= ', rdmft_occ_tol, ' (rdmft_occ_tol)'
          EXIT
       ENDIF
       !
       ls_ok = .FALSE.
       alpha_acc = 0.0_DP
       !
       IF (use_wolfe) THEN
          IF (use_nm) THEN
             ref_energy = rdmft_zhang_hager_ref(zh_occ)
             CALL rdmft_strong_wolfe_ls(rdmft_ebi_ls_eval, E0, phi0, alpha0, &
                  rdmft_line_search_c1, rdmft_line_search_c2, &
                  rdmft_line_search_max_iter, rdmft_line_search_max_zoom, lsres, &
                  f_ref=ref_energy)
          ELSE
             CALL rdmft_strong_wolfe_ls(rdmft_ebi_ls_eval, E0, phi0, alpha0, &
                  rdmft_line_search_c1, rdmft_line_search_c2, &
                  rdmft_line_search_max_iter, rdmft_line_search_max_zoom, lsres)
          ENDIF
          IF (lsres%success .AND. lsres%step > 0.0_DP) THEN
             alpha_acc = lsres%step
             DO ik = 1, nks
                DO ib = 1, nbnd
                   x(ib, ik) = ebi_x0(ib, ik) + alpha_acc * ebi_dir(ib, ik)
                ENDDO
             ENDDO
             CALL rdmft_ebi_sync_n_from_x(x, wk, nbnd, nks)
             CALL rdmft_total_energy(etot_now)
             ls_ok = .TRUE.
          ELSE IF (rdmft_verbose >= 1 .AND. ionode) THEN
             WRITE(stdout, '(7X,A)') &
                  TRIM(tag)//' block: strong Wolfe failed; Armijo fallback.'
          ENDIF
       ENDIF
       !
       IF (.NOT. ls_ok) THEN
          alpha0 = rdmft_occ_ls_alpha0(rdmft_occ_ls_init_step, x_flat, grad_x_flat, &
               rdmft_occ_ls_stepsize, bb, qhist, rdmft_bb_alpha_min, rdmft_bb_alpha_max)
          IF (alpha0 <= 0.0_DP) alpha0 = MAX(rdmft_occ_ls_stepsize, rdmft_bgd_tau)
          CALL ebi_pool_bcast_dp(alpha0)
          ref_energy = E0
          IF (use_nm) ref_energy = rdmft_zhang_hager_ref(zh_occ)
          CALL rdmft_armijo_ls(rdmft_ebi_ls_eval, E0, phi0, alpha0, &
               rdmft_line_search_c1, rdmft_line_search_rho, &
               rdmft_line_search_max_iter, lsres, rdmft_ebi_ls_post_reject, &
               f_ref=ref_energy, alpha_floor=rdmft_bb_alpha_min, &
               use_polynomial=rdmft_line_search_polynomial, log_tag=tag)
          IF (lsres%success .AND. lsres%step > 0.0_DP) THEN
             alpha_acc = lsres%step
             etot_now = lsres%f_new
             x = ebi_trial
             ls_ok = .TRUE.
          ENDIF
       ENDIF
       !
       IF (.NOT. ls_ok) THEN
          x = x_save
          CALL rdmft_ebi_sync_n_from_x(x, wk, nbnd, nks)
          CALL rdmft_invalidate_xi_cache()
          IF (rdmft_occ_tol > 0.0_DP .AND. ABS(phi0) <= rdmft_occ_tol) THEN
             occ_converged = .TRUE.
             IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A)') &
                  TRIM(tag)//' block converged: |phi0|=', ABS(phi0), &
                  ' ~ 0 (Armijo exhausted at stationary point).'
          ELSE
             CALL rdmft_total_energy(etot_now)
             IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A)') &
                  TRIM(tag)//' block: line search failed; keeping current x.'
          ENDIF
          EXIT
       ENDIF
       !
       CALL rdmft_ebi_fold_mu(x, wk, nbnd, nks, mu_up, mu_dw)
       CALL rdmft_ebi_sync_n_from_x(x, wk, nbnd, nks)
       IF (use_nm) CALL rdmft_zhang_hager_update(zh_occ, etot_now)
       !
       sum_dn = SUM(ABS(rdmft_n - n_save))
       CALL mp_sum(sum_dn, inter_pool_comm)
       ne_now = rdmft_weighted_sum(rdmft_n, wk, nbnd, nks)
       cviol = rdmft_constraint_violation_occ(rdmft_n, wk, nbnd, nks)
       CALL rdmft_grad_n(grad_n)
       CALL rdmft_pg_kkt_residual(rdmft_n, grad_n, wk, nbnd, nks, resid=kkt_post)
       ! grad_n / kkt_post are now current for the accepted iterate; the
       ! next inner iteration reuses them (grad_ready at the loop head).
       grad_ready = .TRUE.
       IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, &
            '(7X,A,I4,A,F18.10,A,1PE10.2,A,1PE10.2,A,1PE10.2,A,1PE10.2,0P,A,F12.6)') &
            TRIM(tag)//' inner ', inner, ' E=', etot_now, ' KKT=', kkt_post, &
            " phi0=", phi0, ' alpha=', alpha_acc, ' sum|dn|=', sum_dn, ' Ne=', ne_now
       CALL rdmft_inner_log_record('occ', inner, 0, etot_now, phi0, &
            alpha_acc, sum_dn, ne_now, kkt_resid=kkt_post)
       !
       IF (sum_dn < 0.0_DP .OR. cviol > MAX(10.0_DP * rdmft_occ_tol, 1.0e-10_DP)) THEN
          x = x_save
          CALL rdmft_ebi_sync_n_from_x(x, wk, nbnd, nks)
          CALL rdmft_invalidate_xi_cache()
          CALL rdmft_total_energy(etot_now)
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A)') &
               TRIM(tag)//' block: rejecting invalid step (sum|dn|=', sum_dn, &
               ', Ne violation=', cviol, '); keeping current x.'
          EXIT
       ENDIF
       !
       CALL rdmft_occ_eps_log_record(inner, grad_n)
       !
       IF (rdmft_occ_tol > 0.0_DP .AND. sum_dn < rdmft_occ_tol) THEN
          occ_converged = .TRUE.
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A)') &
               TRIM(tag)//' block converged: sum|dn|=', sum_dn, &
               ' < ', rdmft_occ_tol, ' (rdmft_occ_tol)'
          EXIT
       ENDIF
    ENDDO
    !
    DEALLOCATE(x, x_save, grad_n, grad_x, dir, n_save, x_flat, grad_x_flat)
    !
  END SUBROUTINE rdmft_ebi_occ_block
  !
END MODULE rdmft_ebi
