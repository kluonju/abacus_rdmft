!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_spg
  !------------------------------------------------------------------
  !! SPG2 occupation inner loop (Euclidean L2 proximal projection ``P_w``,
  !! BMR spectral steplength, optional ELK preconditioner).
  !
  USE kinds, ONLY : DP
  USE mp,        ONLY : mp_sum
  USE mp_pools,  ONLY : inter_pool_comm
  USE rdmft_linesearch, ONLY : rdmft_bb_state
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_spg_occ_block
  PUBLIC :: rdmft_spg_spectral_alpha
  !! Exposed for the joint SPG product-manifold driver
  !! (:subroutine:`rdmft_run_joint_spg`), which reuses the BMR spectral
  !! steplength for its occupation chord.
  !
  REAL(DP), ALLOCATABLE, SAVE :: spg_n0(:,:), spg_grad_dir(:,:), spg_trial(:,:)
  REAL(DP), ALLOCATABLE, SAVE :: spg_grad(:,:), spg_grad_trial(:,:)
  REAL(DP), ALLOCATABLE, SAVE :: spg_occ_flat(:), spg_grad_flat(:)
  REAL(DP), ALLOCATABLE, SAVE :: spg_wk(:)
  REAL(DP), ALLOCATABLE, SAVE :: spg_prev_dir(:,:), spg_prev_grad(:,:)
  REAL(DP), ALLOCATABLE, SAVE :: spg_dir(:,:)
  !! Birgin-Martinez SPG search direction
  !!   d_k = P(n_0 - alpha_BB * grad_dir) - n_0
  !! computed ONCE per inner iteration; the line search then walks
  !! the **straight** segment ``n_lambda = n_0 + lambda * d_k`` for
  !! ``lambda in [0, 1]`` (no extra projections needed because both
  !! endpoints are feasible and the feasible set is convex).
  INTEGER, SAVE :: spg_nbnd = 0, spg_nks = 0
  LOGICAL, SAVE :: spg_use_wolfe_grad = .FALSE.
  REAL(DP), SAVE :: spg_phi0 = 0.0_DP
  !! Chord slope \(\phi'(0)=\langle g(n_0), d_k\rangle\) computed once
  !! per inner iteration in :subroutine:`rdmft_spg_occ_block`.  The
  !! Armijo evaluator returns it directly instead of re-reducing the
  !! same inner product (one ``mp_sum`` per trial) -- the returned
  !! slope is ignored by the Armijo driver anyway.
  !
CONTAINS
  !
  SUBROUTINE rdmft_spg_ls_eval(alpha, f, g, ierr)
    !! Line-search evaluator for the **standard projected gradient
    !! method** (Birgin-Martinez SPG2): walks the straight segment
    !!
    !!   n(lambda) = n_0 + lambda * d_k,
    !!   d_k       = P(n_0 - alpha_BB * grad_dir) - n_0,
    !!
    !! with ``lambda = alpha`` here.  ``d_k`` is the SPG descent
    !! direction: it is built ONCE per outer inner iteration (in
    !! :subroutine:`rdmft_spg_occ_block`) and then re-used for every
    !! line-search trial.  Because the feasible set
    !! \(\{0 \le n_{ik} \le 1, \sum w_k n_{ik} = N_e\}\) is convex
    !! and ``n_0`` and ``n_0 + d_k = P(n_0 - alpha_BB * grad_dir)`` both
    !! lie in it, every \(\lambda \in [0, 1]\) gives a feasible
    !! point WITHOUT any further projection -- the line is straight
    !! in occupation space.  For \(\lambda > 1\) the trial may
    !! escape the feasible set and we project, falling back to the
    !! same projection-aware secant slope the old curved-line code
    !! used.
    !!
    !! Compared with the previous
    !!   n(alpha) = P(n_0 - alpha * grad_dir)
    !! ("PG line search", curved line through the proximal arc),
    !! the SPG straight line gives:
    !!   * an exact, finite slope ``phi'(0) = <grad, d_k>`` (no FD
    !!     probe needed);
    !!   * a well-defined directional derivative
    !!     ``phi'(lambda) = <grad(n_lambda), d_k>`` for every
    !!     \(\lambda \in [0, 1]\), even when the projection clipped
    !!     bands when building ``d_k`` -- the line search itself is
    !!     projection-free in this regime;
    !!   * standard Strong-Wolfe / Armijo guarantees.
    USE rdmft_module, ONLY : rdmft_n
    USE rdmft_energy, ONLY : rdmft_grad_n, rdmft_total_energy
    USE rdmft_occupation, ONLY : rdmft_proximal_project_occ, &
                                  rdmft_constraint_violation_occ
    REAL(DP), INTENT(IN)  :: alpha
    REAL(DP), INTENT(OUT) :: f, g
    INTEGER, INTENT(OUT) :: ierr
    REAL(DP) :: cviol, phi_prime, inv_alpha
    INTEGER :: ib, ik
    LOGICAL :: extrapolated
    !
    ierr = 0
    extrapolated = (alpha > 1.0_DP)
    !
    ! Straight-line step n(lambda) = n_0 + lambda * d_k.
    DO ik = 1, spg_nks
       DO ib = 1, spg_nbnd
          spg_trial(ib, ik) = spg_n0(ib, ik) + alpha * spg_dir(ib, ik)
       ENDDO
    ENDDO
    !
    ! For lambda in [0, 1] the trial is automatically feasible (convex
    ! combination of two feasible points).  For lambda > 1 we may have
    ! left the feasible set; project to enforce 0 <= n <= 1 and
    ! sum w_k n = N_e.
    IF (extrapolated) THEN
       CALL rdmft_proximal_project_occ(spg_trial, spg_wk, spg_nbnd, spg_nks)
       cviol = rdmft_constraint_violation_occ(spg_trial, spg_wk, spg_nbnd, spg_nks)
       IF (cviol > 1.0e-6_DP) THEN
          ierr = 1
          f = HUGE(1.0_DP)
          g = 0.0_DP
          RETURN
       ENDIF
    ENDIF
    !
    IF (spg_use_wolfe_grad) THEN
       rdmft_n = spg_trial
       CALL rdmft_grad_n(spg_grad_trial)
       IF (.NOT. extrapolated) THEN
          ! Straight-line regime: phi'(lambda) = <g(n_lambda), d_k>.
          phi_prime = 0.0_DP
          DO ik = 1, spg_nks
             DO ib = 1, spg_nbnd
                phi_prime = phi_prime + spg_grad_trial(ib, ik) * spg_dir(ib, ik)
             ENDDO
          ENDDO
          CALL mp_sum(phi_prime, inter_pool_comm)
       ELSE IF (alpha > 0.0_DP) THEN
          ! Re-projected past lambda = 1: secant slope along the arc.
          inv_alpha = 1.0_DP / alpha
          phi_prime = 0.0_DP
          DO ik = 1, spg_nks
             DO ib = 1, spg_nbnd
                phi_prime = phi_prime + spg_grad_trial(ib, ik) &
                     * (spg_trial(ib, ik) - spg_n0(ib, ik)) * inv_alpha
             ENDDO
          ENDDO
          CALL mp_sum(phi_prime, inter_pool_comm)
       ELSE
          phi_prime = 0.0_DP
          DO ik = 1, spg_nks
             DO ib = 1, spg_nbnd
                phi_prime = phi_prime + spg_grad_trial(ib, ik) * spg_dir(ib, ik)
             ENDDO
          ENDDO
          CALL mp_sum(phi_prime, inter_pool_comm)
       ENDIF
    ELSE
       ! No trial gradient: the slope is only consumed by the Wolfe
       ! drivers (the Armijo driver ignores it), so return the cached
       ! phi'(0) = <g(n_0), d_k> computed once per inner iteration
       ! instead of redoing the O(nbnd*nks) sum + inter-pool mp_sum
       ! on every backtracking trial.
       phi_prime = spg_phi0
    ENDIF
    rdmft_n = spg_trial
    CALL rdmft_total_energy(f)
    g = phi_prime
  END SUBROUTINE rdmft_spg_ls_eval
  !
  SUBROUTINE rdmft_spg_ls_post_reject()
    !! Restore the SPG line origin after a rejected trial.
    USE rdmft_module, ONLY : rdmft_n, rdmft_invalidate_xi_cache
    USE rdmft_energy, ONLY : rdmft_set_wg_from_n
    !
    rdmft_n = spg_n0
    CALL rdmft_set_wg_from_n(1)
    CALL rdmft_invalidate_xi_cache()
    !
  END SUBROUTINE rdmft_spg_ls_post_reject
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_spg_spectral_alpha(bb, x, g, alpha_min, alpha_max, &
                                              alpha_init)
    !---------------------------------------------------------------
    !! BMR SPG2 spectral steplength (Birgin--Mart\'inez--Raydan 1999/2000,
    !! Algorithm~2.2 Step~3).  Given the last accepted pair
    !! ``(bb%x_prev, bb%g_prev)`` and the current ``(x, g)``,
    !!
    !!   \alpha_{k+1} = \min\{\alpha_{\max},
    !!       \max\{\alpha_{\min}, \langle s,s\rangle / \langle s,y\rangle\}\}
    !!
    !! with ``s = x - x_{\rm prev}``, ``y = g - g_{\rm prev}``, and the
    !! Euclidean inner product ``\langle a,b\rangle = \sum_{ik} a_{ik} b_{ik}``.
    !! Here ``g`` is the physical gradient ``grad_n = \partial E/\partial n``.
    !! When ``\langle s,y\rangle \le 0`` fall back to ``alpha_init``.
    TYPE(rdmft_bb_state), INTENT(IN) :: bb
    REAL(DP), INTENT(IN) :: x(:), g(:), alpha_min, alpha_max, alpha_init
    REAL(DP) :: sy, ss, alpha, si, yi
    INTEGER :: i, n
    !
    rdmft_spg_spectral_alpha = MIN(alpha_max, MAX(alpha_min, alpha_init))
    IF (.NOT. bb%have_prev) RETURN
    IF (.NOT. ALLOCATED(bb%x_prev) .OR. .NOT. ALLOCATED(bb%g_prev)) RETURN
    n = SIZE(x)
    IF (n < 1) RETURN
    IF (n /= SIZE(g) .OR. n /= SIZE(bb%x_prev) .OR. n /= SIZE(bb%g_prev)) RETURN
    sy = 0.0_DP
    ss = 0.0_DP
    DO i = 1, n
       si = x(i) - bb%x_prev(i)
       yi = g(i) - bb%g_prev(i)
       sy = sy + si * yi
       ss = ss + si * si
    ENDDO
    CALL mp_sum(sy, inter_pool_comm)
    CALL mp_sum(ss, inter_pool_comm)
    IF (sy <= 0.0_DP) THEN
       rdmft_spg_spectral_alpha = MIN(alpha_max, MAX(alpha_min, alpha_init))
       RETURN
    ENDIF
    alpha = ss / sy
    rdmft_spg_spectral_alpha = MIN(alpha_max, MAX(alpha_min, alpha))
  END FUNCTION rdmft_spg_spectral_alpha
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_spg_occ_block(etot_now, occ_converged, max_inner, log_tag)
    USE io_global, ONLY : stdout, ionode
    USE wvfct, ONLY : nbnd
    USE klist, ONLY : nks, wk
    USE lsda_mod, ONLY : isk, lsda
    USE rdmft_module
    USE rdmft_energy, ONLY : rdmft_grad_n, rdmft_total_energy, &
                              rdmft_compute_elk_occ_scale, rdmft_build_occ_grad_dir, &
                              rdmft_set_wg_from_n
    USE rdmft_occupation, ONLY : rdmft_proximal_project_occ, &
                                  rdmft_constraint_violation_occ, &
                                  rdmft_pg_map_grad_inf, rdmft_pg_kkt_residual, &
                                  rdmft_weighted_sum
    USE rdmft_linesearch, ONLY : rdmft_gll_state, rdmft_gll_init, rdmft_gll_push, &
                                  rdmft_gll_ref, &
                                  rdmft_zhang_hager_state, rdmft_zhang_hager_init, &
                                  rdmft_zhang_hager_update, rdmft_zhang_hager_ref, &
                                  rdmft_bb_state, rdmft_bb_init, rdmft_bb_reset, &
                                  rdmft_bb_record, &
                                  rdmft_quad_ls_history, rdmft_quad_ls_reset, &
                                  rdmft_quad_ls_update, rdmft_ls_result, &
                                  rdmft_strong_wolfe_ls, &
                                  rdmft_ls_is_wolfe, rdmft_ls_is_nm, rdmft_ls_is_gll, &
                                  rdmft_armijo_ls, rdmft_resolve_occ_ls_type
    USE rdmft_lbfgs, ONLY : rdmft_lbfgs_state, rdmft_lbfgs_init, rdmft_lbfgs_reset, &
                            rdmft_lbfgs_direction, rdmft_lbfgs_update, rdmft_lbfgs_finalize
    USE rdmft_inner_log, ONLY : rdmft_inner_log_record
    USE rdmft_occ_eps_log, ONLY : rdmft_occ_eps_log_record
    USE mp,        ONLY : mp_sum
    USE mp_pools,  ONLY : inter_pool_comm
    !
    REAL(DP), INTENT(INOUT) :: etot_now
    LOGICAL, INTENT(OUT)   :: occ_converged
    INTEGER, INTENT(IN)    :: max_inner
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: log_tag
    !
    REAL(DP), ALLOCATABLE :: grad_n(:,:), grad_dir(:,:), n_save(:,:), elk_scale(:,:)
    REAL(DP), ALLOCATABLE :: occ_flat(:), grad_flat(:), grad_n_flat(:), dir_flat(:)
    REAL(DP), ALLOCATABLE :: step_flat(:), grad_diff_flat(:)
    REAL(DP) :: g1_inf, g1_post, kkt_resid, kkt_post, alpha0, alpha_spectral, alpha_acc, phi0, phi_trial, E0
    REAL(DP) :: sum_dn, cviol, ref_energy, dE, ne_now, alpha_init_spg
    REAL(DP) :: dn2_e, lemma_rhs_e
    REAL(DP) :: beta_pr, denom_pr, dd_check
    REAL(DP), PARAMETER :: spg2_g1_step = 1.0_DP
    !! BMR SPG2 stopping map: \(g_1 = P_w(n - t\,\partial E/\partial n) - n\), \(t=1\).
    INTEGER :: inner, ndim, ik, ib
    LOGICAL :: ls_ok, step_applied, use_wolfe, use_nm, use_gll, recompute_elk
    LOGICAL :: use_cg, use_lbfgs, grad_ready
    CHARACTER(LEN=24) :: ls_eff
    TYPE(rdmft_gll_state) :: gll_occ
    TYPE(rdmft_zhang_hager_state) :: zh_occ
    TYPE(rdmft_bb_state) :: bb
    TYPE(rdmft_quad_ls_history) :: qhist
    TYPE(rdmft_ls_result) :: lsres
    TYPE(rdmft_lbfgs_state) :: occ_lbfgs
    CHARACTER(LEN=32) :: tag
    CHARACTER(LEN=9) :: slbl
    !
    occ_converged = .TRUE.
    IF (max_inner <= 0) RETURN
    !
    tag = 'occ'
    IF (PRESENT(log_tag)) tag = TRIM(log_tag)
    !
    ALLOCATE(grad_n(nbnd, nks), grad_dir(nbnd, nks), &
         n_save(nbnd, nks), elk_scale(nbnd, nks))
    ndim = nbnd * nks
    ALLOCATE(occ_flat(ndim), grad_flat(ndim), grad_n_flat(ndim), dir_flat(ndim))
    ALLOCATE(step_flat(ndim), grad_diff_flat(ndim))
    !
    spg_nbnd = nbnd
    spg_nks = nks
    IF (.NOT. ALLOCATED(spg_n0)) ALLOCATE(spg_n0(nbnd, nks), spg_grad_dir(nbnd, nks), &
         spg_trial(nbnd, nks), spg_grad(nbnd, nks), spg_grad_trial(nbnd, nks))
    IF (.NOT. ALLOCATED(spg_dir)) ALLOCATE(spg_dir(nbnd, nks))
    IF (.NOT. ALLOCATED(spg_wk)) ALLOCATE(spg_wk(nks))
    spg_wk = wk
    !
    CALL rdmft_bb_init(bb, ndim, rdmft_spg_alpha_min, rdmft_spg_alpha_max)
    CALL rdmft_quad_ls_reset(qhist)
    use_gll = rdmft_ls_is_gll(rdmft_occ_ls_type)
    IF (use_gll) THEN
       CALL rdmft_gll_init(gll_occ, etot_now, memory=rdmft_spg_ls_memory)
    ENDIF
    !
    recompute_elk = .TRUE.
    occ_converged = .FALSE.
    grad_ready = .FALSE.
    !
    use_lbfgs = (TRIM(rdmft_occ_optimizer) == 'lbfgs') &
           .OR. (TRIM(rdmft_occ_optimizer) == 'bfgs')  &
           .OR. (TRIM(rdmft_occ_optimizer) == 'l-bfgs')
    use_cg = (TRIM(rdmft_occ_optimizer) == 'cg')
    IF (TRIM(rdmft_constraint) == 'projected_gradient') THEN
       SELECT CASE (TRIM(rdmft_occ_optimizer))
       CASE ('spg2', 'sd', '')
          IF (ionode) WRITE(stdout, '(7X,A)') &
               'SPG2 (canonical): BMR spectral steplength + straight chord LS (lambda=1)'
       CASE ('cg')
          IF (ionode) WRITE(stdout, '(7X,A)') &
               'SPG2 + Polak-Ribiere CG direction'
       CASE ('lbfgs', 'bfgs', 'l-bfgs')
          IF (ionode) WRITE(stdout, '(7X,A,I0,A)') &
               'SPG2 + L-BFGS direction (memory=', rdmft_lbfgs_memory, ')'
       CASE DEFAULT
          IF (ionode) WRITE(stdout, '(7X,A,A,A)') &
               'SPG2: unknown rdmft_occ_optimizer ''', &
               TRIM(rdmft_occ_optimizer), &
               ''' -- falling back to canonical SPG2 (spectral steepest descent).'
       END SELECT
    ENDIF
    !
    IF (use_cg .AND. .NOT. ALLOCATED(spg_prev_dir)) ALLOCATE(spg_prev_dir(nbnd, nks))
    IF (use_cg .AND. .NOT. ALLOCATED(spg_prev_grad)) ALLOCATE(spg_prev_grad(nbnd, nks))
    IF (use_lbfgs) CALL rdmft_lbfgs_init(occ_lbfgs, ndim, rdmft_lbfgs_memory)
    !
    CALL rdmft_resolve_occ_ls_type(rdmft_occ_ls_type, ls_eff)
    use_wolfe = rdmft_ls_is_wolfe(ls_eff)
    use_nm = rdmft_ls_is_nm(ls_eff)
    IF (use_wolfe .AND. use_nm) &
         CALL rdmft_zhang_hager_init(zh_occ, etot_now, rdmft_zhang_hager_eta)
    spg_use_wolfe_grad = use_wolfe
    IF (rdmft_verbose >= 1 .AND. ionode) THEN
       IF (use_gll .AND. .NOT. use_wolfe) THEN
          WRITE(stdout, '(7X,A,I0,A,I0,A)') &
               TRIM(tag)//' block (SPG2 / GLL nonmonotone Armijo, M=', &
               rdmft_spg_ls_memory, ', max ', max_inner, ' inner iters)'
       ELSE
          WRITE(stdout, '(7X,A,A,A,I0,A)') &
               TRIM(tag)//' block (SPG / ', TRIM(ls_eff), ', max ', max_inner, ' inner iters)'
       ENDIF
    ENDIF
    !
    DO inner = 1, max_inner
       ! Reuse the post-step gradient from the previous inner iteration
       ! when nothing (occupations, orbitals, density) has changed since
       ! it was computed: the accepted-step bookkeeping below evaluates
       ! rdmft_grad_n at the new n for the BB pair / KKT report, and the
       ! state is untouched between that call and this loop head.  This
       ! removes one full occupation-gradient evaluation (one-body
       ! h_psi sweep + per-channel exchange build) per inner iteration.
       ! The reuse decision is uniform across MPI ranks (step_applied
       ! and the loop control flow are collective-consistent), so all
       ! ranks skip the same collectives together.
       IF (grad_ready) THEN
          grad_ready = .FALSE.
          ! g1 / KKT were evaluated post-step from the same (n, grad):
          ! reuse them too (rdmft_pg_map_grad_inf runs a collective
          ! proximal bisection that need not be repeated).
          g1_inf = g1_post
          kkt_resid = kkt_post
       ELSE
          CALL rdmft_grad_n(grad_n)
          !
          ! SPG2 stopping map and KKT residual use physical ``grad_n`` only.
          CALL rdmft_pg_map_grad_inf(rdmft_n, grad_n, wk, nbnd, nks, &
               spg2_g1_step, g1_inf)
          CALL rdmft_pg_kkt_residual(rdmft_n, grad_n, wk, nbnd, nks, resid=kkt_resid)
       ENDIF
       spg_grad = grad_n
       ! NOTE: rdmft_weighted_sum performs an mp_sum over inter_pool_comm
       ! and is therefore COLLECTIVE.  It MUST be evaluated on every rank
       ! (the IF condition below is uniform across ranks because
       ! rdmft_verbose is a broadcast control variable); restricting the
       ! call to ionode would desynchronise the inter-pool reductions
       ! between k-point pools and corrupt every subsequent reduction.
       IF (rdmft_verbose >= 2) THEN
          ne_now = rdmft_weighted_sum(rdmft_n, wk, nbnd, nks)
          IF (ionode) WRITE(stdout, &
               '(9X,A,I3,A,1PE10.2,A,1PE10.2,0P,A,F12.6)') &
               '['//TRIM(tag)//' inner ', inner, '] KKT_resid=', kkt_resid, &
               ' ||g_1||_inf=', g1_inf, ' Ne=', ne_now
       ENDIF
       !
       ! Near an SPG2 stationary point the projected-gradient map vanishes.
       ! Do not build a chord or run a line search (avoids alpha ~ 1/||g_1||
       ! and spurious occupation changes when bands sit on the box).
       IF (rdmft_occ_grad_tol > 0.0_DP .AND. g1_inf <= rdmft_occ_grad_tol) THEN
          occ_converged = .TRUE.
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, &
               '(7X,A,1PE10.2,A,1PE10.2,A)') &
               TRIM(tag)//' block converged: ||g_1||_inf=', g1_inf, &
               ' <= ', rdmft_occ_grad_tol, ' (rdmft_occ_grad_tol)'
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, &
               '(7X,A,I4,A,F18.10,A,1PE10.2,A,1PE10.2,A,1PE10.2,A,1PE10.2)') &
               TRIM(tag)//' inner ', inner, ' E=', etot_now, ' KKT=', kkt_resid, &
               ' ||g_1||_inf=', g1_inf, ' lambda=', 0.0_DP, ' sum|dn|=', 0.0_DP
          ne_now = rdmft_weighted_sum(rdmft_n, wk, nbnd, nks)
          CALL rdmft_inner_log_record('occ', inner, 0, etot_now, g1_inf, &
               0.0_DP, 0.0_DP, ne_now, kkt_resid=kkt_resid)
          CALL rdmft_occ_eps_log_record(inner, grad_n)
          EXIT
       ENDIF
       !
       IF (rdmft_occ_precond) THEN
          IF (recompute_elk) THEN
             CALL rdmft_compute_elk_occ_scale(elk_scale)
             recompute_elk = .FALSE.
          ENDIF
       ENDIF
       !
       ! ---- Search direction for SPG proximal step ----
       !
       ! ``rdmft_grad_n`` returns the physical gradient
       ! ``g_{ik}=\partial E/\partial n_{ik}`` (includes ``w_k``).  SPG2
       ! uses the Euclidean L2 proximal map ``P_w`` (Section 5.4) and
       ! BMR Euclidean spectral steplength.  Chord slopes
       ! ``\phi'(0)=\sum_{ik} g_{ik} d_{ik}``.
       CALL rdmft_build_occ_grad_dir(grad_n, elk_scale, grad_dir)
       grad_flat = RESHAPE(grad_dir, [ndim])
       IF (use_lbfgs) THEN
          CALL rdmft_lbfgs_direction(occ_lbfgs, grad_flat, dir_flat)
          ! ``dir_flat = -H grad_dir`` (minimisation sense).
          grad_dir = RESHAPE(-dir_flat, [nbnd, nks])
          dd_check = 0.0_DP
          DO ik = 1, nks
             DO ib = 1, nbnd
                dd_check = dd_check + grad_dir(ib, ik) * grad_n(ib, ik)
             ENDDO
          ENDDO
          CALL mp_sum(dd_check, inter_pool_comm)
          IF (dd_check <= 0.0_DP) THEN
             CALL rdmft_lbfgs_reset(occ_lbfgs)
             CALL rdmft_build_occ_grad_dir(grad_n, elk_scale, grad_dir)
          ENDIF
       ELSE IF (use_cg) THEN
          IF (inner > 1) THEN
          beta_pr = 0.0_DP
          denom_pr = 0.0_DP
          DO ik = 1, nks
             DO ib = 1, nbnd
                beta_pr = beta_pr + grad_dir(ib,ik) * (grad_dir(ib,ik) - spg_prev_grad(ib,ik))
                denom_pr = denom_pr + spg_prev_grad(ib,ik) * spg_prev_grad(ib,ik)
             ENDDO
          ENDDO
          ! Occupation arrays are k-pool distributed: reduce the CG inner
          ! products so beta is identical on every pool (otherwise the
          ! search direction -- and hence the collective line-search trial
          ! sequence -- diverges between pools).
          CALL mp_sum(beta_pr, inter_pool_comm)
          CALL mp_sum(denom_pr, inter_pool_comm)
          IF (denom_pr > 1.0e-30_DP) THEN
             beta_pr = MAX(0.0_DP, beta_pr / denom_pr)  ! Hager safeguard
          ELSE
             beta_pr = 0.0_DP
          ENDIF
          ! grad_dir = g - beta * grad_dir_prev
          DO ik = 1, nks
             DO ib = 1, nbnd
                grad_dir(ib,ik) = grad_dir(ib,ik) - beta_pr * spg_prev_dir(ib,ik)
             ENDDO
          ENDDO
          ! Descent check: grad_dir^T grad_n > 0.
          dd_check = 0.0_DP
          DO ik = 1, nks
             DO ib = 1, nbnd
                dd_check = dd_check + grad_dir(ib,ik) * grad_n(ib,ik)
             ENDDO
          ENDDO
          CALL mp_sum(dd_check, inter_pool_comm)
          IF (dd_check <= 0.0_DP) THEN
             CALL rdmft_build_occ_grad_dir(grad_n, elk_scale, grad_dir)
          ENDIF
          ENDIF
          spg_prev_dir = grad_dir
          spg_prev_grad = RESHAPE(grad_flat, [nbnd, nks])
       ENDIF
       spg_grad_dir = grad_dir
       !
       occ_flat = RESHAPE(rdmft_n, [ndim])
       grad_flat = RESHAPE(grad_dir, [ndim])
       grad_n_flat = RESHAPE(grad_n, [ndim])
       !
       ! BMR SPG2 spectral steplength alpha_k (Algorithm 2.2 Step 3).
       ! First inner step: alpha_0 = 1 / ||g_1(x_0)||_inf (BMR experiments).
       ! Later steps: safeguarded inverse Rayleigh quotient from accepted
       ! (n, g) pairs.  Uses raw grad_n, not grad_dir.
       IF (g1_inf > 1.0e-30_DP) THEN
          alpha_init_spg = 1.0_DP / g1_inf
       ELSE
          alpha_init_spg = rdmft_spg_alpha_max
       ENDIF
       alpha_init_spg = MIN(rdmft_spg_alpha_max, &
            MAX(rdmft_spg_alpha_min, alpha_init_spg))
       alpha_spectral = rdmft_spg_spectral_alpha(bb, occ_flat, grad_n_flat, &
            rdmft_spg_alpha_min, rdmft_spg_alpha_max, alpha_init_spg)
       IF (alpha_spectral <= 0.0_DP) alpha_spectral = alpha_init_spg
       IF (rdmft_verbose >= 2 .AND. ionode) WRITE(stdout, &
            '(11X,A,1PE10.2,A,1PE10.2,A,1PE10.2,A,1PE10.2,A)') &
            'SPG2 spectral alpha_k=', alpha_spectral, ' ||g_1||_inf=', g1_inf, &
            ' clamp [', rdmft_spg_alpha_min, ',', rdmft_spg_alpha_max, ']'
       !
       n_save = rdmft_n
       spg_n0 = n_save
       spg_grad_dir = grad_dir
       spg_grad = grad_n
       E0 = etot_now
       !
       ! ---- Standard projected-gradient (SPG2 / Birgin-Martinez) -----
       !
       ! The earlier implementation did a curved-line search:
       !
       !   n(alpha) = P(n_0 - alpha * grad_dir)
       !
       ! Each Wolfe / Armijo trial re-projected onto the bound /
       ! sum-of-occupations constraint, so the line was not straight in
       ! occupation space and the directional derivative
       ! ``phi'(0) = -<grad, grad_dir>`` overstated the actual slope of the
       ! projected arc by orders of magnitude on near-stationary
       ! metallic systems (NiO 2x2x2 / 1x1x1 inner #10 traces).
       !
       ! Standard SPG2 builds
       !
       !   d_k = P_w( n_0 - alpha_k * g_k ) - n_0,
       !   g_{ik} = \partial E/\partial n_{ik}  (physical).
       !
       ! once per inner iteration and walks the straight chord
       ! n(lambda) = n_0 + lambda d_k.  CG / L-BFGS substitute
       ! ``grad_dir`` for ``g`` in the proximal preimage.
       DO ik = 1, nks
          DO ib = 1, nbnd
             IF (use_cg .OR. use_lbfgs) THEN
                spg_trial(ib, ik) = spg_n0(ib, ik) &
                     - alpha_spectral * spg_grad_dir(ib, ik)
             ELSE
                spg_trial(ib, ik) = spg_n0(ib, ik) &
                     - alpha_spectral * grad_n(ib, ik)
             ENDIF
          ENDDO
       ENDDO
       CALL rdmft_proximal_project_occ(spg_trial, wk, nbnd, nks)
       DO ik = 1, nks
          DO ib = 1, nbnd
             spg_dir(ib, ik) = spg_trial(ib, ik) - spg_n0(ib, ik)
          ENDDO
       ENDDO
       !
       ! phi'(0) = <g, d_k>.
       phi0 = 0.0_DP
       DO ik = 1, nks
          DO ib = 1, nbnd
             phi0 = phi0 + grad_n(ib, ik) * spg_dir(ib, ik)
          ENDDO
       ENDDO
       CALL mp_sum(phi0, inter_pool_comm)
       spg_phi0 = phi0
       !
       ! BMR Lemma 2.1: with Euclidean ``P_w``, ``d = P_w(n - alpha u) - n``,
       ! ``g^T d <= -||d||^2 / alpha`` (``alpha`` = spectral steplength).
       ! ``dn2_e`` (and its ``mp_sum``) feed ONLY the verbose>=2 sanity
       ! print below: gate the whole O(nbnd*nks) reduction + collective
       ! behind the same check as its consumers instead of paying for
       ! an extra inter-pool ``mp_sum`` on every inner iteration
       ! regardless of verbosity (``rdmft_verbose`` is a broadcast
       ! control variable, so this is uniform across ranks).
       IF (rdmft_verbose >= 2) THEN
          dn2_e = 0.0_DP
          DO ik = 1, nks
             DO ib = 1, nbnd
                dn2_e = dn2_e + spg_dir(ib, ik) * spg_dir(ib, ik)
             ENDDO
          ENDDO
          CALL mp_sum(dn2_e, inter_pool_comm)
          IF (alpha_spectral > 0.0_DP) THEN
             lemma_rhs_e = -dn2_e / alpha_spectral
          ELSE
             lemma_rhs_e = 0.0_DP
          ENDIF
          IF (ionode) WRITE(stdout, &
               '(11X,A,1PE12.4,A,1PE12.4)') &
               'SPG chord lemma (BMR): <g,d>=', phi0, '  -||d||^2/alpha=', lemma_rhs_e
          IF (ionode .AND. phi0 > lemma_rhs_e + 1.0e-10_DP) &
               WRITE(stdout, '(11X,A)') &
               'WARNING: Euclidean projection lemma violated (P_w bug?)'
       ENDIF
       !
       ! SPG chord uphill (``phi0 >= 0``): the proximal projection
       ! clipped the pre-image so ``d_k`` is not a descent direction.
       IF (phi0 >= 0.0_DP) THEN
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A)') &
               TRIM(tag)//' block: SPG chord uphill (phi0=', phi0, ').'
          ls_ok = .FALSE.
          alpha_acc = 0.0_DP
          step_applied = .FALSE.
          GOTO 1001
       ENDIF
       !
       ls_ok = .FALSE.
       alpha_acc = 0.0_DP
       step_applied = .FALSE.
       !
       ! BMR SPG2 Step 2.1: first chord trial at lambda = 1; backtrack
       ! along the same straight direction if the nonmonotone test fails.
       alpha0 = 1.0_DP
       IF (use_wolfe) THEN
          IF (use_nm) THEN
             ref_energy = rdmft_zhang_hager_ref(zh_occ)
             CALL rdmft_strong_wolfe_ls(rdmft_spg_ls_eval, E0, phi0, alpha0, &
                  rdmft_spg_gamma, rdmft_line_search_c2, &
                  rdmft_line_search_max_iter, rdmft_line_search_max_zoom, lsres, &
                  f_ref=ref_energy)
          ELSE IF (INDEX(ls_eff, 'weak') > 0 .OR. ls_eff == 'wolfe') THEN
             CALL rdmft_strong_wolfe_ls(rdmft_spg_ls_eval, E0, phi0, alpha0, &
                  rdmft_spg_gamma, rdmft_line_search_c2, &
                  rdmft_line_search_max_iter, rdmft_line_search_max_zoom, lsres, &
                  weak_curv=.TRUE.)
          ELSE
             CALL rdmft_strong_wolfe_ls(rdmft_spg_ls_eval, E0, phi0, alpha0, &
                  rdmft_spg_gamma, rdmft_line_search_c2, &
                  rdmft_line_search_max_iter, rdmft_line_search_max_zoom, lsres)
          ENDIF
          IF (lsres%success .AND. lsres%step > 0.0_DP) THEN
             alpha_acc = lsres%step
             etot_now = lsres%f_new
             rdmft_n = spg_trial
             ls_ok = .TRUE.
          ELSE IF (rdmft_verbose >= 1 .AND. ionode) THEN
             WRITE(stdout, '(7X,A)') &
                  TRIM(tag)//' block: strong Wolfe failed; Armijo fallback.'
          ENDIF
       ENDIF
       IF (.NOT. ls_ok) THEN
          ! Armijo backtracking on the straight SPG chord (GLL nonmonotone
          ! when rdmft_occ_ls_type = gll/nm_armijo; monotone when auto/armijo).
          alpha0 = 1.0_DP
          IF (rdmft_verbose >= 2 .AND. ionode) WRITE(stdout, &
               '(11X,A,1PE10.2)') 'Armijo LS seed lambda0=', alpha0
          ref_energy = E0
          IF (use_gll) THEN
             ref_energy = rdmft_gll_ref(gll_occ)
          ELSE IF (use_nm) THEN
             ref_energy = rdmft_zhang_hager_ref(zh_occ)
          ENDIF
          CALL rdmft_armijo_ls(rdmft_spg_ls_eval, E0, phi0, alpha0, &
               rdmft_spg_gamma, rdmft_line_search_rho, &
               rdmft_line_search_max_iter, lsres, rdmft_spg_ls_post_reject, &
               f_ref=ref_energy, alpha_floor=rdmft_spg_alpha_min, &
               allow_stall=.TRUE., use_polynomial=rdmft_line_search_polynomial, &
               log_tag=tag, sigma_lo=rdmft_spg_sigma1, sigma_hi=rdmft_spg_sigma2, &
               bmr_halve_outside=.TRUE.)
          IF (lsres%success) THEN
             alpha_acc = lsres%step
             etot_now = lsres%f_new
             ls_ok = .TRUE.
          ENDIF
       ENDIF
       !
       ! Jumped-to label for the chord-uphill branch (``phi0 >= 0``).
       1001 CONTINUE
       !
       IF (ls_ok) THEN
          IF (alpha_acc > 0.0_DP) THEN
             IF (.NOT. use_wolfe) rdmft_n = spg_trial
             step_applied = .TRUE.
          ELSE
             ! Stall accept at the numerical noise floor (zero step).
             rdmft_n = n_save
             CALL rdmft_set_wg_from_n(1)
             CALL rdmft_invalidate_xi_cache()
             step_applied = .FALSE.
          ENDIF
          IF (step_applied) THEN
             ! ELK occupation preconditioner ``elk_scale`` is a diagonal
             ! conditioning heuristic obtained from O(nbnd*nks) full
             ! ``rdmft_grad_n`` probes at n=0.5.  Within this occupation
             ! block the orbitals are frozen, so it is kept from the first
             ! inner iteration and NOT rebuilt after every accepted step
             ! (recompute_elk left .FALSE.), saving nbnd*nks gradient
             ! evaluations per inner iteration.  It is refreshed once at the
             ! next block entry (recompute_elk = .TRUE.) after the orbital
             ! block moves the wavefunctions.  Off by default
             ! (rdmft_occ_precond = .FALSE.).
             !
             ! verbose >= 2: per-band/k occupation changes (dn) after accepted
             ! step.  Gather across pools (collective) so spin-down rows print.
             IF (rdmft_verbose >= 2) THEN
                IF (ionode) WRITE(stdout, '(9X,A,I3,A,1PE10.2,A)') &
                     '['//TRIM(tag)//' inner ', inner, &
                     '] step accepted (alpha=', alpha_acc, '). occ change (ik, spin, ib, n_new, dn):'
                CALL rdmft_spg_dump_occ('', inner, '', spg_trial, n_save, .TRUE.)
             ENDIF
             IF (use_gll) THEN
                CALL rdmft_gll_push(gll_occ, etot_now)
             ELSE IF (use_nm) THEN
                CALL rdmft_zhang_hager_update(zh_occ, etot_now)
             ENDIF
             CALL rdmft_quad_ls_update(qhist, alpha_acc, E0, etot_now, phi0)
          ENDIF
       ELSE
          rdmft_n = n_save
          CALL rdmft_total_energy(etot_now)
          IF (use_lbfgs) CALL rdmft_lbfgs_reset(occ_lbfgs)
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A)') &
               TRIM(tag)//' block: line search failed; keeping current n.'
          EXIT
       ENDIF
       !
       IF (ls_ok .AND. .NOT. step_applied) THEN
          kkt_post = kkt_resid
          g1_post = g1_inf
          sum_dn = 0.0_DP
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A)') &
               TRIM(tag)//' block: occ line search at numerical floor; no step.'
          EXIT
       ENDIF
       !
       sum_dn = SUM(ABS(rdmft_n - n_save))
       CALL mp_sum(sum_dn, inter_pool_comm)
       IF (step_applied) THEN
          CALL rdmft_grad_n(grad_n)
          ! grad_n is now the gradient at the accepted iterate; the next
          ! inner iteration reuses it instead of recomputing (see the
          ! grad_ready check at the loop head).
          grad_ready = .TRUE.
          occ_flat = RESHAPE(rdmft_n, [ndim])
          grad_n_flat = RESHAPE(grad_n, [ndim])
          CALL rdmft_bb_record(bb, occ_flat, grad_n_flat)
          IF (use_lbfgs) THEN
             step_flat = RESHAPE(rdmft_n - n_save, [ndim])
             CALL rdmft_build_occ_grad_dir(grad_n, elk_scale, grad_dir)
             grad_diff_flat = RESHAPE(grad_dir, [ndim]) - grad_flat
             CALL rdmft_lbfgs_update(occ_lbfgs, step_flat, grad_diff_flat)
          ENDIF
          CALL rdmft_pg_kkt_residual(rdmft_n, grad_n, wk, nbnd, nks, resid=kkt_post)
          CALL rdmft_pg_map_grad_inf(rdmft_n, grad_n, wk, nbnd, nks, &
               spg2_g1_step, g1_post)
       ENDIF
       IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, &
            '(7X,A,I4,A,F18.10,A,1PE10.2,A,1PE10.2,A,1PE10.2,A,1PE10.2)') &
            TRIM(tag)//' inner ', inner, ' E=', etot_now, ' KKT=', kkt_post, &
            ' ||g_1||_inf=', g1_post, ' lambda=', alpha_acc, ' sum|dn|=', sum_dn
       IF (step_applied) THEN
          ne_now = rdmft_weighted_sum(rdmft_n, wk, nbnd, nks)
          CALL rdmft_inner_log_record('occ', inner, 0, etot_now, g1_post, &
               alpha_acc, sum_dn, ne_now, kkt_resid=kkt_post)
          CALL rdmft_occ_eps_log_record(inner, grad_n)
       ENDIF
       ! SPG2 convergence is governed solely by the projected-gradient
       ! (BMR stopping) map ``||g_1||_inf <= rdmft_occ_grad_tol`` or by
       ! exhausting ``max_inner`` inner iterations.  The check is done here
       ! against the freshly computed post-step ``g1_post`` (and again at
       ! the loop head, which reuses it), so a step accepted on the final
       ! inner iteration still reports convergence.
       !
       ! The occupation-change criterion ``sum|dn| < rdmft_occ_tol`` was
       ! removed: an accepted step can be arbitrarily small (spectral
       ! steplength clamp, projection clipping bands onto the box) far from
       ! an SPG2 stationary point, so ``sum|dn|`` declaring convergence
       ! reported the block as converged while ``||g_1||`` was still large.
       IF (step_applied .AND. rdmft_occ_grad_tol > 0.0_DP .AND. &
            g1_post <= rdmft_occ_grad_tol) THEN
          occ_converged = .TRUE.
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A)') &
               TRIM(tag)//' block converged: ||g_1||_inf=', g1_post, &
               ' <= ', rdmft_occ_grad_tol, ' (rdmft_occ_grad_tol)'
          EXIT
       ENDIF
    ENDDO
    !
    IF (use_lbfgs) CALL rdmft_lbfgs_finalize(occ_lbfgs)
    DEALLOCATE(grad_n, grad_dir, n_save, elk_scale, occ_flat, grad_flat, grad_n_flat, dir_flat)
    DEALLOCATE(step_flat, grad_diff_flat)
    !
  END SUBROUTINE rdmft_spg_occ_block
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_spg_dump_occ(header, inner, title, occ_local, ref_local, show_dn)
    !---------------------------------------------------------------
    !! Print the per-(ik, spin, ib) occupation table for the verbose
    !! SPG diagnostics.  k-points are distributed across MPI pools, so
    !! with nspin=2 + npool>1 the spin-up and spin-down rows live on
    !! different pools.  We therefore gather the per-pool occupation
    !! arrays into a global (nbnd, nkstot) table with \texttt{poolcollect}
    !! (collective over inter_pool_comm -- every rank must call it) and
    !! only ionode writes the gathered table.  This is what makes BOTH
    !! spin channels appear in the log for a spin-polarised run.
    USE io_global, ONLY : stdout, ionode
    USE wvfct,     ONLY : nbnd
    USE klist,     ONLY : nks, nkstot
    USE lsda_mod,  ONLY : lsda, isk
    !
    CHARACTER(LEN=*), INTENT(IN) :: header, title
    INTEGER, INTENT(IN) :: inner
    REAL(DP), INTENT(IN) :: occ_local(nbnd, nks), ref_local(nbnd, nks)
    LOGICAL, INTENT(IN)  :: show_dn
    !
    REAL(DP), ALLOCATABLE :: occ_g(:,:), ref_g(:,:), isk_loc(:,:), isk_g(:,:)
    INTEGER :: ik, ib, ispin
    CHARACTER(LEN=9) :: slbl
    !
    ALLOCATE(occ_g(nbnd, nkstot), isk_loc(1, nks), isk_g(1, nkstot))
    DO ik = 1, nks
       isk_loc(1, ik) = REAL(isk(ik), DP)
    ENDDO
    CALL poolcollect(nbnd, nks, occ_local, nkstot, occ_g)
    CALL poolcollect(1, nks, isk_loc, nkstot, isk_g)
    IF (show_dn) THEN
       ALLOCATE(ref_g(nbnd, nkstot))
       CALL poolcollect(nbnd, nks, ref_local, nkstot, ref_g)
    ENDIF
    !
    IF (ionode) THEN
       IF (LEN_TRIM(header) > 0) WRITE(stdout, '(9X,A,I3,A)') TRIM(header), inner, TRIM(title)
       DO ik = 1, nkstot
          ispin = NINT(isk_g(1, ik))
          slbl = MERGE('spin up  ', 'spin down', .NOT. (lsda .AND. ispin == 2))
          IF (.NOT. lsda) slbl = 'spinless '
          DO ib = 1, nbnd
             IF (show_dn) THEN
                ! ``n_old``, ``n_new`` and ``dn = n_new - n_old`` on the
                ! same line: enough information to spot every accepted
                ! occupation update at a glance without redundant
                ! per-iteration "n=" dumps.
                WRITE(stdout, '(11X,A,I3,2X,A9,A,I3,A,F12.7,A,F12.7,A,1PE11.3)') &
                     'ik=', ik, slbl, '  ib=', ib, &
                     '  n_old=', ref_g(ib, ik), &
                     '  n_new=', occ_g(ib, ik), &
                     '  dn=', occ_g(ib, ik) - ref_g(ib, ik)
             ELSE
                ! Standalone ``n=`` dump is intentionally a no-op: the
                ! caller is expected to use show_dn=.TRUE. so each line
                ! carries n_old + n_new + dn together.  Kept as a
                ! compile-time stub so the show_dn signature does not
                ! need to change throughout the SPG driver.
                CONTINUE
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    IF (show_dn) DEALLOCATE(ref_g)
    DEALLOCATE(occ_g, isk_loc, isk_g)
    !
  END SUBROUTINE rdmft_spg_dump_occ
  !
END MODULE rdmft_spg
