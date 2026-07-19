!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_bgd
  !------------------------------------------------------------------
  !! ELK-style box-aware occupation gradient descent (``bgd``) inner loop.
  !!
  !! This is a faithful port of ELK's \texttt{rdmvaryn.f90}
  !! occupation-number update (the RDMFT ``task 300`` occupation
  !! minimisation, \texttt{rdmminn} / \texttt{rdmvaryn}), with one
  !! addition: an explicit energy **line search** along the ELK search
  !! direction.  ELK's original routine takes a single feasibility-only
  !! backtracking step from \texttt{taurdmn} (with a hard-coded
  !! geometric factor of ``0.75``); here we keep ELK's search
  !! direction unchanged and walk the (straight, feasible) segment
  !! ``n(tau) = n_0 + tau * gamma`` looking for a sufficient energy
  !! decrease (Armijo) using the shared occupation line-search driver
  !! (:func:`rdmft_armijo_ls` in :mod:`rdmft_linesearch`) with
  !! :var:`rdmft_line_search_rho` and optional quadratic refinement.
  !!
  !! ## ELK search direction (see \texttt{rdmvaryn.f90})
  !!
    !! Let \(g_{ik}=\partial E/\partial n_{ik}\) be the raw occupation
    !! gradient from :func:`rdmft_grad_n` (includes \(w_k\)).  ELK
    !! stores ``dedn``\(_{ik}=-g_{ik}/w_k\).  With
    !! \(t_{ik}=\texttt{dedn}_{ik}-\kappa\), the box-aware direction is
    !! \[
    !!   \gamma_{ik} =
    !!   \begin{cases}
    !!     t_{ik}\,(n_{\max}-n_{ik}), & t_{ik} > 0,\\
    !!     t_{ik}\, n_{ik},          & t_{ik} \le 0,
    !!   \end{cases}
    !! \]
  !! where the scalar \(\kappa\) (a chemical potential) is chosen so
  !! that the move is charge-neutral, \(\sum_k w_k\sum_i\gamma_{ik}=0\).
  !! ELK then normalises \(\gamma\) when its weighted square norm
  !! exceeds 1 and finds a step \(\tau\le\texttt{taurdmn}\) keeping
  !! every \(n_{ik}+\tau\gamma_{ik}\in[0,n_{\max}]\).
  !!
  !! The same construction makes \(\gamma\) a descent direction for the
  !! total energy: with \(p_{ik}=w_k\times(\text{box factor})\ge 0\),
    !! \(\phi'(0)=\sum_{ik}(\partial E/\partial n_{ik})\,\gamma_{ik}
    !!   = -P\,\mathrm{Var}_p(\texttt{dedn})\le 0\), vanishing exactly at a
  !! KKT point (all active per-state gradients equal).  Because
  !! \(\sum_k w_k\gamma_{ik}=0\), the electron number is conserved
  !! **exactly** along the whole line (no projection needed), and the
  !! box constraint holds for \(\tau\in[0,\tau_{\max}]\); the line
  !! search is therefore projection-free on that interval.
  !
  USE kinds, ONLY : DP
  USE mp,        ONLY : mp_sum, mp_min
  USE mp_pools,  ONLY : inter_pool_comm
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_bgd_occ_block, rdmft_bgd_build_direction, &
            rdmft_bgd_solve_kappa
  !
  REAL(DP), PARAMETER :: bgd_occmax = 1.0_DP
  !! Maximum occupation per (band, k, spin) in the PWscf RDMFT
  !! convention (\(n_{ik}\in[0,1]\)).
  !
  ! Shared state for the line-search evaluator (the straight segment
  ! n(tau) = bgd_n0 + tau * bgd_dir).
  REAL(DP), ALLOCATABLE, SAVE :: bgd_n0(:,:), bgd_dir(:,:), bgd_trial(:,:)
  REAL(DP), ALLOCATABLE, SAVE :: bgd_grad(:,:)
  INTEGER, SAVE :: bgd_nbnd = 0, bgd_nks = 0
  REAL(DP), SAVE :: bgd_tau_feas = 0.0_DP, bgd_phi0 = 0.0_DP
  LOGICAL, SAVE :: bgd_use_wolfe_grad = .FALSE.
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_bgd_solve_kappa(dedn, n, wk, nbnd, nks, ispin, use_filter, kappa)
    !---------------------------------------------------------------
    !! Find the chemical-potential shift \(\kappa\) so that
    !! \(\sum_k w_k\sum_i\gamma_{ik}=0\) over the selected set of
    !! k-points (all of them, or only those with
    !! \texttt{isk(ik)=ispin} when \texttt{use\_filter}).  Direct port
    !! of the bracketing loop in ELK \texttt{rdmvaryn.f90}.
    USE io_global, ONLY : stdout, ionode
    USE lsda_mod,  ONLY : isk
    !
    REAL(DP), INTENT(IN)  :: dedn(nbnd, nks), n(nbnd, nks), wk(nks)
    INTEGER,  INTENT(IN)  :: nbnd, nks, ispin
    LOGICAL,  INTENT(IN)  :: use_filter
    REAL(DP), INTENT(OUT) :: kappa
    !
    INTEGER, PARAMETER :: maxit = 10000
    REAL(DP), PARAMETER :: eps = 1.0e-12_DP
    REAL(DP) :: gs, gsp, dgs, dkapa, sumsq, t1, gam
    INTEGER  :: it, ib, ik
    !
    gsp   = 0.0_DP
    kappa = 0.0_DP
    dkapa = 0.1_DP
    DO it = 1, maxit
       gs    = 0.0_DP
       sumsq = 0.0_DP
       DO ik = 1, nks
          IF (use_filter .AND. isk(ik) /= ispin) CYCLE
          DO ib = 1, nbnd
             t1 = dedn(ib, ik) - kappa
             IF (t1 > 0.0_DP) THEN
                gam = t1 * (bgd_occmax - n(ib, ik))
             ELSE
                gam = t1 * n(ib, ik)
             ENDIF
             gs    = gs    + wk(ik) * gam
             sumsq = sumsq + wk(ik) * gam * gam
          ENDDO
       ENDDO
       ! Occupations / weights are distributed across k-point pools:
       ! reduce so every rank sees the SAME (gs, sumsq) and therefore
       ! takes an identical sequence of kappa updates (otherwise the
       ! per-pool kappa -- and the collective line search built on top
       ! of gamma -- would diverge between pools).
       CALL mp_sum(gs,    inter_pool_comm)
       CALL mp_sum(sumsq, inter_pool_comm)
       sumsq = SQRT(sumsq)
       sumsq = MAX(sumsq, 1.0_DP)
       t1 = ABS(gs) / sumsq
       IF (t1 < eps) RETURN
       IF (it >= 2) THEN
          dgs = gs - gsp
          IF (gs * dgs > 0.0_DP) dkapa = -dkapa
          IF (gs * gsp < 0.0_DP) THEN
             dkapa = 0.5_DP * dkapa
          ELSE
             dkapa = 1.1_DP * dkapa
          ENDIF
       ENDIF
       gsp = gs
       kappa = kappa + dkapa
    ENDDO
    ! ELK aborts here; we keep the last kappa and warn instead so the
    ! line search can still try the (slightly charge-non-neutral)
    ! direction -- the proximal handling downstream tolerates it.
    IF (ionode) WRITE(stdout, '(7X,A)') &
         'bgd block: WARNING rdmft_bgd_solve_kappa did not converge; using last offset.'
    !
  END SUBROUTINE rdmft_bgd_solve_kappa
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_bgd_build_direction(n, grad_n, wk, nbnd, nks, gamma, tau_feas, &
                                      gamma_active)
    !---------------------------------------------------------------
    !! Build the ELK reduced-gradient direction \(\gamma\) and the
    !! feasibility step bound \(\tau_{\max}\).  Mirrors ELK
    !! ``rdmvaryn.f90`` (``dedn`` assembly, \(\kappa\) bracketing,
    !! weighted \(\gamma\), weighted-norm scaling).
    !!
    !! ``grad_n`` is the physical occupation gradient
    !! \(\partial E/\partial n_{ik}\) from :func:`rdmft_grad_n` (includes
    !! \(w_k\)).  ELK ``dedn`` uses the unweighted per-state quantity
    !! ``dedn = -(\partial E/\partial n)/w_k``.
    USE rdmft_module, ONLY : rdmft_fix_magnetization, rdmft_verbose
    USE rdmft_energy, ONLY : rdmft_dedn_from_grad_n
    USE io_global,    ONLY : ionode, stdout
    USE lsda_mod,     ONLY : lsda, isk
    !
    REAL(DP), INTENT(IN)  :: n(nbnd, nks), grad_n(nbnd, nks), wk(nks)
    INTEGER,  INTENT(IN)  :: nbnd, nks
    REAL(DP), INTENT(OUT) :: gamma(nbnd, nks)
    REAL(DP), INTENT(OUT) :: tau_feas
    LOGICAL,  INTENT(OUT) :: gamma_active
    !
    REAL(DP) :: dedn(nbnd, nks), dedn_print(nbnd, nks)
    REAL(DP) :: kappa_up, kappa_dw, kappa_all, kappa, t1, gam, sumsq, cap, gs_check
    INTEGER  :: ib, ik
    LOGICAL  :: split_spin
    !
    ! ELK ``rdmdedn`` internal array uses ``-dF/dn``; ``rdmwritededn`` prints
    ! ``-dedn_internal``.  Build the printed table first, then flip sign
    ! for the ``rdmvaryn`` bracketing loop.
    CALL rdmft_dedn_from_grad_n(grad_n, wk, nbnd, nks, dedn_print)
    dedn = -dedn_print
    !
    split_spin = (rdmft_fix_magnetization .AND. lsda)
    ! ELK ``rdmvaryn`` uses one global ``kapa`` for all spins; when
    ! ``rdmft_fix_magnetization`` is set we solve separate ``kappa`` per
    ! spin channel so each spin block is charge-neutral independently.
    IF (split_spin) THEN
       CALL rdmft_bgd_solve_kappa(dedn, n, wk, nbnd, nks, 1, .TRUE., kappa_up)
       CALL rdmft_bgd_solve_kappa(dedn, n, wk, nbnd, nks, 2, .TRUE., kappa_dw)
    ELSE
       CALL rdmft_bgd_solve_kappa(dedn, n, wk, nbnd, nks, 0, .FALSE., kappa_all)
    ENDIF
    !
    DO ik = 1, nks
       IF (split_spin) THEN
          kappa = MERGE(kappa_up, kappa_dw, isk(ik) == 1)
       ELSE
          kappa = kappa_all
       ENDIF
       DO ib = 1, nbnd
          t1 = dedn(ib, ik) - kappa
          IF (t1 > 0.0_DP) THEN
             gam = t1 * (bgd_occmax - n(ib, ik))
          ELSE
             gam = t1 * n(ib, ik)
          ENDIF
          gamma(ib, ik) = gam
       ENDDO
    ENDDO
    !
    ! Normalise gamma when its weighted square norm exceeds 1 (ELK).
    ! Scaling by a constant preserves the charge-neutral zero sum.
    sumsq = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          sumsq = sumsq + wk(ik) * gamma(ib, ik) * gamma(ib, ik)
       ENDDO
    ENDDO
    CALL mp_sum(sumsq, inter_pool_comm)
    IF (sumsq > 1.0_DP) THEN
       t1 = 1.0_DP / SQRT(sumsq)
       gamma = t1 * gamma
    ENDIF
    !
    ! Largest feasible step tau keeping every n + tau*gamma in [0, occmax].
    tau_feas = HUGE(1.0_DP)
    DO ik = 1, nks
       DO ib = 1, nbnd
          IF (gamma(ib, ik) > 1.0e-30_DP) THEN
             cap = (bgd_occmax - n(ib, ik)) / gamma(ib, ik)
             tau_feas = MIN(tau_feas, MAX(cap, 0.0_DP))
          ELSE IF (gamma(ib, ik) < -1.0e-30_DP) THEN
             cap = -n(ib, ik) / gamma(ib, ik)
             tau_feas = MIN(tau_feas, MAX(cap, 0.0_DP))
          ENDIF
       ENDDO
    ENDDO
    CALL mp_min(tau_feas, inter_pool_comm)
    !
    ! Charge neutrality of the weighted direction (ELK ``gs`` at
    ! convergence).  Should be ~0 after ``rdmft_bgd_solve_kappa``.
    gs_check = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          gs_check = gs_check + wk(ik) * gamma(ib, ik)
       ENDDO
    ENDDO
    CALL mp_sum(gs_check, inter_pool_comm)
    IF (rdmft_verbose >= 2 .AND. ionode) WRITE(stdout, &
         '(11X,A,1PE12.4)') 'ELK weighted gamma sum w_k*gamma_k =', gs_check
    !
    gamma_active = (sumsq > 1.0e-30_DP .AND. tau_feas > 1.0e-30_DP &
                    .AND. tau_feas < HUGE(1.0_DP))
    !
  END SUBROUTINE rdmft_bgd_build_direction
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_bgd_ls_eval(alpha, f, g, ierr)
    !---------------------------------------------------------------
    !! Line-search evaluator for the ELK ``bgd`` direction.  Walks the
    !! straight, feasible segment ``n(alpha) = bgd_n0 + alpha * bgd_dir``
    !! for ``alpha in [0, bgd_tau_feas]``.  Returns the total energy
    !! ``f`` and the directional derivative ``g = <grad(n), bgd_dir>``
    !! (Nocedal-Wright line-search contract).
    USE rdmft_module, ONLY : rdmft_n, rdmft_occ_tol
    USE rdmft_energy, ONLY : rdmft_grad_n, rdmft_total_energy
    USE rdmft_occupation, ONLY : rdmft_constraint_violation_occ
    USE klist, ONLY : wk
    !
    REAL(DP), INTENT(IN)  :: alpha
    REAL(DP), INTENT(OUT) :: f, g
    INTEGER, INTENT(OUT)  :: ierr
    REAL(DP) :: phi_prime, cviol
    INTEGER  :: ib, ik
    !
    ierr = 0
    ! Beyond the feasibility bound the box constraint would be
    ! violated; reject (the Wolfe driver treats ierr/=0 as a hard
    ! bracket boundary, Armijo never extrapolates past it).
    IF (alpha > bgd_tau_feas * (1.0_DP + 1.0e-12_DP) + 1.0e-30_DP) THEN
       ierr = 1
       f = HUGE(1.0_DP)
       g = 0.0_DP
       RETURN
    ENDIF
    !
    DO ik = 1, bgd_nks
       DO ib = 1, bgd_nbnd
          bgd_trial(ib, ik) = bgd_n0(ib, ik) + alpha * bgd_dir(ib, ik)
       ENDDO
    ENDDO
    rdmft_n = bgd_trial
    cviol = rdmft_constraint_violation_occ(bgd_trial, wk, bgd_nbnd, bgd_nks)
    IF (cviol > MAX(10.0_DP * rdmft_occ_tol, 1.0e-10_DP)) THEN
       ierr = 1
       f = HUGE(1.0_DP)
       g = 0.0_DP
       RETURN
    ENDIF
    !
    IF (bgd_use_wolfe_grad) THEN
       CALL rdmft_grad_n(bgd_grad)
       phi_prime = 0.0_DP
       DO ik = 1, bgd_nks
          DO ib = 1, bgd_nbnd
             phi_prime = phi_prime + bgd_grad(ib, ik) * bgd_dir(ib, ik)
          ENDDO
       ENDDO
       CALL mp_sum(phi_prime, inter_pool_comm)
    ELSE
       phi_prime = bgd_phi0
    ENDIF
    !
    CALL rdmft_total_energy(f)
    g = phi_prime
    !
  END SUBROUTINE rdmft_bgd_ls_eval
  !
  SUBROUTINE rdmft_bgd_ls_post_reject()
    !! Restore the line origin and EXX hooks after a rejected trial.
    USE rdmft_module, ONLY : rdmft_n, rdmft_invalidate_xi_cache
    USE rdmft_energy, ONLY : rdmft_set_wg_from_n, rdmft_refresh_x_occupation
    !
    rdmft_n = bgd_n0
    CALL rdmft_set_wg_from_n(1)
    CALL rdmft_refresh_x_occupation()
    CALL rdmft_invalidate_xi_cache()
    !
  END SUBROUTINE rdmft_bgd_ls_post_reject
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_bgd_occ_block(etot_now, occ_converged, max_inner, log_tag)
    !---------------------------------------------------------------
    !! ELK-style occupation gradient-descent inner block with line
    !! search.  Drop-in replacement for
    !! :subroutine:`rdmft_spg_occ_block` when
    !! \texttt{rdmft\_occ\_optimizer = 'bgd'}.
    USE io_global,  ONLY : stdout, ionode
    USE wvfct,      ONLY : nbnd
    USE klist,      ONLY : nks, wk
    USE rdmft_module
    USE rdmft_energy, ONLY : rdmft_grad_n, rdmft_total_energy, &
                              rdmft_set_wg_from_n, rdmft_refresh_x_occupation, &
                              rdmft_build_occ_grad_dir, rdmft_compute_elk_occ_scale
    USE rdmft_occupation, ONLY : rdmft_weighted_sum, &
                                  rdmft_constraint_violation_occ, &
                                  rdmft_pg_kkt_residual
    USE rdmft_linesearch, ONLY : rdmft_zhang_hager_state, rdmft_zhang_hager_init, &
                                  rdmft_zhang_hager_update, rdmft_zhang_hager_ref, &
                                  rdmft_ls_result, rdmft_strong_wolfe_ls, &
                                  rdmft_ls_is_wolfe, rdmft_ls_is_nm, &
                                  rdmft_resolve_occ_ls_type, rdmft_armijo_ls, &
                                  rdmft_bb_state, rdmft_bb_init, rdmft_bb_record, &
                                  rdmft_quad_ls_history, rdmft_quad_ls_reset, &
                                  rdmft_occ_ls_alpha0
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
    REAL(DP), ALLOCATABLE :: grad_n(:,:), grad_dir(:,:), gamma(:,:), n_save(:,:)
    REAL(DP), ALLOCATABLE :: elk_scale(:,:), occ_flat(:), grad_flat(:)
    REAL(DP) :: phi0, E0, alpha0, alpha_acc, tau_feas, sum_dn
    REAL(DP) :: ref_energy, ne_now, cviol, alpha_floor, kkt_resid, kkt_post
    REAL(DP) :: alpha_bb, alpha_cap
    INTEGER  :: inner, ndim, ik, ib
    LOGICAL  :: ls_ok, gamma_active, use_wolfe, use_nm, recompute_elk
    LOGICAL  :: grad_ready
    CHARACTER(LEN=24) :: ls_eff
    CHARACTER(LEN=32) :: tag
    TYPE(rdmft_ls_result) :: lsres
    TYPE(rdmft_zhang_hager_state) :: zh_occ
    TYPE(rdmft_bb_state) :: bb
    TYPE(rdmft_quad_ls_history) :: qhist
    !
    occ_converged = .TRUE.
    IF (max_inner <= 0) RETURN
    !
    tag = 'occ'
    IF (PRESENT(log_tag)) tag = TRIM(log_tag)
    !
    ALLOCATE(grad_n(nbnd, nks), grad_dir(nbnd, nks), gamma(nbnd, nks), &
         n_save(nbnd, nks), elk_scale(nbnd, nks))
    bgd_nbnd = nbnd
    bgd_nks  = nks
    IF (.NOT. ALLOCATED(bgd_n0))    ALLOCATE(bgd_n0(nbnd, nks))
    IF (.NOT. ALLOCATED(bgd_dir))   ALLOCATE(bgd_dir(nbnd, nks))
    IF (.NOT. ALLOCATED(bgd_trial)) ALLOCATE(bgd_trial(nbnd, nks))
    IF (.NOT. ALLOCATED(bgd_grad))  ALLOCATE(bgd_grad(nbnd, nks))
    !
    CALL rdmft_resolve_occ_ls_type(rdmft_occ_ls_type, ls_eff)
    use_wolfe = rdmft_ls_is_wolfe(ls_eff)
    use_nm    = rdmft_ls_is_nm(ls_eff)
    bgd_use_wolfe_grad = use_wolfe
    CALL rdmft_zhang_hager_init(zh_occ, etot_now, rdmft_zhang_hager_eta)
    ndim = nbnd * nks
    ALLOCATE(occ_flat(ndim), grad_flat(ndim))
    CALL rdmft_bb_init(bb, ndim, rdmft_bb_alpha_min, rdmft_bb_alpha_max)
    CALL rdmft_quad_ls_reset(qhist)
    !
    IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,A,A,I0,A)') &
         TRIM(tag)//' block (bgd / ', TRIM(ls_eff), ', max ', max_inner, ' inner iters)'
    !
    recompute_elk = .TRUE.
    occ_converged = .FALSE.
    grad_ready = .FALSE.
    DO inner = 1, max_inner
       ! Reuse the post-step gradient computed at the end of the
       ! previous inner iteration (same n, orbitals, and density) --
       ! this saves one full occupation-gradient evaluation per inner
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
       IF (rdmft_occ_precond) THEN
          IF (recompute_elk) THEN
             CALL rdmft_compute_elk_occ_scale(elk_scale)
             recompute_elk = .FALSE.
          ENDIF
       ENDIF
       CALL rdmft_build_occ_grad_dir(grad_n, elk_scale, grad_dir)
       !
       ! ELK reduced-gradient direction + feasibility bound.
       CALL rdmft_bgd_build_direction(rdmft_n, grad_n, wk, nbnd, nks, gamma, &
            tau_feas, gamma_active)
       IF (.NOT. gamma_active) THEN
          occ_converged = .TRUE.
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, &
               '(7X,A,A,1PE10.2,A)') &
               TRIM(tag)//' block converged: |gamma| ~ 0 (stationary).  KKT=', &
               kkt_resid
          EXIT
       ENDIF
       !
       n_save  = rdmft_n
       bgd_n0   = rdmft_n
       bgd_dir  = gamma
       bgd_tau_feas = tau_feas
       E0 = etot_now
       !
       ! Energy line-search slope phi'(0) = <grad_n, gamma>.
       phi0 = 0.0_DP
       DO ik = 1, nks
          DO ib = 1, nbnd
             phi0 = phi0 + grad_n(ib, ik) * gamma(ib, ik)
          ENDDO
       ENDDO
       CALL mp_sum(phi0, inter_pool_comm)
       bgd_phi0 = phi0
       !
       IF (phi0 >= 0.0_DP) THEN
          occ_converged = .TRUE.
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A)') &
               TRIM(tag)//' block converged: phi0=', phi0, ' >= 0 (no descent).'
          EXIT
       ENDIF
       !
       ! Near a KKT point the directional derivative vanishes; do not
       ! enter a line search that backtracks to tau ~ 1e-9 and accepts
       ! spurious energy drops from stale EXX / density state.
       IF (rdmft_occ_tol > 0.0_DP .AND. ABS(phi0) <= rdmft_occ_tol) THEN
          occ_converged = .TRUE.
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A)') &
               TRIM(tag)//' block converged: |phi0|=', ABS(phi0), &
               ' <= ', rdmft_occ_tol, ' (rdmft_occ_tol)'
          EXIT
       ENDIF
       !
       ! Barzilai-Borwein (or fixed / quadratic) seed, capped by the
       ! ELK feasibility bound and ``rdmft_bgd_tau``.
       alpha_cap = MIN(rdmft_bgd_tau, tau_feas)
       IF (alpha_cap <= 0.0_DP) alpha_cap = tau_feas
       occ_flat = RESHAPE(rdmft_n, [ndim])
       grad_flat = RESHAPE(grad_dir, [ndim])
       alpha_bb = rdmft_occ_ls_alpha0(rdmft_occ_ls_init_step, occ_flat, grad_flat, &
            rdmft_occ_ls_stepsize, bb, qhist, rdmft_bb_alpha_min, rdmft_bb_alpha_max)
       IF (alpha_bb <= 0.0_DP) alpha_bb = rdmft_occ_ls_stepsize
       CALL rdmft_bb_record(bb, occ_flat, grad_flat)
       alpha0 = MIN(alpha_cap, alpha_bb)
       IF (alpha0 <= 0.0_DP) alpha0 = MIN(alpha_cap, rdmft_occ_ls_stepsize)
       IF (TRIM(rdmft_occ_ls_init_step) == 'fixed') &
            alpha0 = MIN(alpha_cap, rdmft_occ_ls_stepsize)
       !
       ls_ok = .FALSE.
       alpha_acc = 0.0_DP
       !
       IF (use_wolfe) THEN
          IF (use_nm) THEN
             ref_energy = rdmft_zhang_hager_ref(zh_occ)
             CALL rdmft_strong_wolfe_ls(rdmft_bgd_ls_eval, E0, phi0, alpha0, &
                  rdmft_line_search_c1, rdmft_line_search_c2, &
                  rdmft_line_search_max_iter, rdmft_line_search_max_zoom, lsres, &
                  f_ref=ref_energy)
          ELSE
             CALL rdmft_strong_wolfe_ls(rdmft_bgd_ls_eval, E0, phi0, alpha0, &
                  rdmft_line_search_c1, rdmft_line_search_c2, &
                  rdmft_line_search_max_iter, rdmft_line_search_max_zoom, lsres)
          ENDIF
          IF (lsres%success .AND. lsres%step > 0.0_DP) THEN
             alpha_acc = lsres%step
             ! Reconstruct the accepted iterate so rdmft_n / density are
             ! consistent with the returned energy (the Wolfe driver may
             ! have left rdmft_n at a different trial).
             DO ik = 1, nks
                DO ib = 1, nbnd
                   rdmft_n(ib, ik) = bgd_n0(ib, ik) + alpha_acc * bgd_dir(ib, ik)
                ENDDO
             ENDDO
             CALL rdmft_total_energy(etot_now)
             ls_ok = .TRUE.
          ELSE IF (rdmft_verbose >= 1 .AND. ionode) THEN
             WRITE(stdout, '(7X,A)') &
                  TRIM(tag)//' block: strong Wolfe failed; Armijo fallback.'
          ENDIF
       ENDIF
       !
       IF (.NOT. ls_ok) THEN
          ! Shared monotone Armijo backtracking on the straight feasible line.
          alpha_floor = MAX(1.0e-6_DP, 1.0e-3_DP * alpha_cap)
          ref_energy = E0
          IF (use_nm) ref_energy = rdmft_zhang_hager_ref(zh_occ)
          CALL rdmft_armijo_ls(rdmft_bgd_ls_eval, E0, phi0, alpha0, &
               rdmft_line_search_c1, rdmft_line_search_rho, &
               rdmft_line_search_max_iter, lsres, rdmft_bgd_ls_post_reject, &
               f_ref=ref_energy, alpha_floor=alpha_floor, &
               use_polynomial=rdmft_line_search_polynomial, log_tag=tag)
          IF (lsres%success .AND. lsres%step > 0.0_DP) THEN
             alpha_acc = lsres%step
             etot_now = lsres%f_new
             ls_ok = .TRUE.
          ENDIF
       ENDIF
       !
       IF (.NOT. ls_ok) THEN
          rdmft_n = n_save
          CALL rdmft_set_wg_from_n(1)
          CALL rdmft_refresh_x_occupation()
          CALL rdmft_invalidate_xi_cache()
          IF (rdmft_occ_tol > 0.0_DP .AND. ABS(phi0) <= rdmft_occ_tol) THEN
             occ_converged = .TRUE.
             IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A)') &
                  TRIM(tag)//' block converged: |phi0|=', ABS(phi0), &
                  ' ~ 0 (Armijo exhausted at stationary point).'
          ELSE
             CALL rdmft_total_energy(etot_now)
             IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A)') &
                  TRIM(tag)//' block: line search failed; keeping current n.'
          ENDIF
          EXIT
       ENDIF
       !
       IF (use_nm) CALL rdmft_zhang_hager_update(zh_occ, etot_now)
       ! ELK occupation preconditioner ``elk_scale`` is a diagonal
       ! conditioning heuristic built from O(nbnd*nks) full ``rdmft_grad_n``
       ! probes at n=0.5.  Within this occupation block the orbitals are
       ! frozen, so it is kept from the first inner iteration and NOT
       ! rebuilt after every accepted step (recompute_elk left .FALSE.),
       ! saving nbnd*nks gradient evaluations per inner iteration.  It is
       ! refreshed once at the next block entry after the orbital block
       ! moves the wavefunctions.  Off by default (rdmft_occ_precond=.FALSE.).
       !
       sum_dn = SUM(ABS(rdmft_n - n_save))
       CALL mp_sum(sum_dn, inter_pool_comm)
       ! Collective (mp_sum inside); evaluate on every rank.
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
            " phi0=", phi0, ' tau=', alpha_acc, ' sum|dn|=', sum_dn, ' Ne=', ne_now
       CALL rdmft_inner_log_record('occ', inner, 0, etot_now, phi0, &
            alpha_acc, sum_dn, ne_now, kkt_resid=kkt_post)
       !
       IF (sum_dn < 0.0_DP .OR. cviol > MAX(10.0_DP * rdmft_occ_tol, 1.0e-10_DP)) THEN
          rdmft_n = n_save
          CALL rdmft_set_wg_from_n(1)
          CALL rdmft_refresh_x_occupation()
          CALL rdmft_invalidate_xi_cache()
          CALL rdmft_total_energy(etot_now)
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A)') &
               TRIM(tag)//' block: rejecting invalid step (sum|dn|=', sum_dn, &
               ', Ne violation=', cviol, '); keeping current n.'
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
    DEALLOCATE(grad_n, grad_dir, gamma, n_save, elk_scale, occ_flat, grad_flat)
    !
  END SUBROUTINE rdmft_bgd_occ_block
  !
END MODULE rdmft_bgd
