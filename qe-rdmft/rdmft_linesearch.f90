!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_linesearch
  !------------------------------------------------------------------
  !! Line-search utilities for the RDMFT solver (Zhang-Hager, Wolfe,
  !! Barzilai-Borwein initial step).  Ported from ABACUS
  !! \texttt{rdmft\_optimizer.h}.
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE rdmft_module, ONLY : rdmft_verbose
  USE mp,        ONLY : mp_sum
  USE mp_pools,  ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_zhang_hager_state, rdmft_zhang_hager_init, &
            rdmft_zhang_hager_update, rdmft_zhang_hager_ref, &
            rdmft_zhang_hager_reset
  PUBLIC :: rdmft_gll_state, rdmft_gll_init, rdmft_gll_push, &
            rdmft_gll_ref, rdmft_gll_reset
  PUBLIC :: rdmft_ls_result, rdmft_cubic_min, rdmft_zoom_interp
  PUBLIC :: rdmft_bb_state, rdmft_bb_init, rdmft_bb_reset, rdmft_bb_suggest, &
            rdmft_bb_record
  PUBLIC :: rdmft_quad_ls_history, rdmft_quad_ls_reset, rdmft_quad_ls_update, &
            rdmft_occ_ls_alpha0
  PUBLIC :: rdmft_default_ls_alpha_init
  PUBLIC :: rdmft_ls_is_wolfe, rdmft_ls_is_nm, rdmft_ls_is_gll, &
            rdmft_resolve_occ_ls_type, rdmft_resolve_orb_ls_type
  PUBLIC :: rdmft_wolfe_zoom, rdmft_strong_wolfe_ls
  PUBLIC :: rdmft_armijo_next_alpha, rdmft_armijo_ls, rdmft_ls_stall_accept
  !
  REAL(DP), PARAMETER :: rdmft_default_ls_alpha_init = 1.0_DP
  REAL(DP), PARAMETER :: rdmft_ls_energy_rtol = 1.0e-11_DP
  !
  TYPE :: rdmft_zhang_hager_state
     REAL(DP) :: eta = 0.85_DP
     REAL(DP) :: Q   = 1.0_DP
     REAL(DP) :: C   = 0.0_DP
     LOGICAL  :: initialised = .FALSE.
  END TYPE rdmft_zhang_hager_state
  !
  TYPE :: rdmft_gll_state
     !! Grippo--Lampariello--Lucidi nonmonotone reference for BMR SPG2:
     !! \(f_{\max} = \max\{f(x_{k-j})\mid 0\le j\le \min\{k,M-1\}\}\).
     INTEGER :: memory = 10
     INTEGER :: count = 0
     INTEGER :: head = 0
     REAL(DP), ALLOCATABLE :: f_hist(:)
     LOGICAL :: initialised = .FALSE.
  END TYPE rdmft_gll_state
  !
  TYPE :: rdmft_ls_result
     REAL(DP) :: step = 0.0_DP
     REAL(DP) :: f_new = 0.0_DP
     REAL(DP) :: alpha_init = 0.0_DP
     INTEGER  :: n_feval = 0
     LOGICAL  :: success = .FALSE.
  END TYPE rdmft_ls_result
  !
  TYPE :: rdmft_bb_state
     LOGICAL :: have_prev = .FALSE.
     REAL(DP), ALLOCATABLE :: x_prev(:), g_prev(:)
     REAL(DP) :: alpha_min = 1.0e-8_DP
     REAL(DP) :: alpha_max = 1.0e2_DP
     !! Standard SPG ``alpha`` clamp: the BB ratio s^T s / s^T y is
     !! clamped to \([\alpha_{\min}, \alpha_{\max}]\) (Birgin-Martinez
     !! 2000, equation (10)).  No additional "trust region" cap on
     !! ``alpha_prev``: the original (now-removed) 5x-previous-step
     !! heuristic was a non-standard kludge that the straight-line
     !! SPG2 line search (PR #48) made unnecessary.
  END TYPE rdmft_bb_state
  !
  TYPE :: rdmft_quad_ls_history
     LOGICAL :: valid = .FALSE.
     REAL(DP) :: alpha = 0.0_DP
     REAL(DP) :: f0 = 0.0_DP
     REAL(DP) :: f1 = 0.0_DP
     REAL(DP) :: dd0 = 0.0_DP
  END TYPE rdmft_quad_ls_history
  !
  !! Line-search eval contract (Nocedal & Wright): \(\phi(\alpha)=f(x+\alpha p)\);
  !! `g` must be \(\phi'(\alpha)=\nabla f(x+\alpha p)^\top p\), not
  !! \(\nabla f\cdot(x_{\mathrm{trial}}-x_0)\).
  ABSTRACT INTERFACE
     SUBROUTINE rdmft_ls_eval_if(alpha, f, g, ierr)
        IMPORT :: DP
        REAL(DP), INTENT(IN)  :: alpha
        REAL(DP), INTENT(OUT) :: f, g
        INTEGER, INTENT(OUT) :: ierr
     END SUBROUTINE rdmft_ls_eval_if
     SUBROUTINE rdmft_ls_post_reject_if()
     END SUBROUTINE rdmft_ls_post_reject_if
  END INTERFACE
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  LOGICAL FUNCTION rdmft_ls_is_wolfe(ls_type)
    CHARACTER(LEN=*), INTENT(IN) :: ls_type
    CHARACTER(LEN=24) :: t
    t = ADJUSTL(ls_type)
    rdmft_ls_is_wolfe = (INDEX(t, 'wolfe') > 0) .OR. (INDEX(t, 'sw') > 0) &
         .OR. (INDEX(t, 'strong') > 0) .OR. (INDEX(t, 'nm') > 0)
  END FUNCTION rdmft_ls_is_wolfe
  !
  LOGICAL FUNCTION rdmft_ls_is_nm(ls_type)
    CHARACTER(LEN=*), INTENT(IN) :: ls_type
    CHARACTER(LEN=24) :: t
    t = ADJUSTL(ls_type)
    rdmft_ls_is_nm = (INDEX(t, 'nm') > 0) .OR. (INDEX(t, 'nonmonotone') > 0)
  END FUNCTION rdmft_ls_is_nm
  !
  LOGICAL FUNCTION rdmft_ls_is_gll(ls_type)
    !! True when the occupation driver should use the BMR/GLL
    !! nonmonotone reference (``gll``, ``nm_armijo``).
    CHARACTER(LEN=*), INTENT(IN) :: ls_type
    CHARACTER(LEN=24) :: t
    t = ADJUSTL(ls_type)
    rdmft_ls_is_gll = (INDEX(t, 'gll') > 0) .OR. (TRIM(t) == 'nm_armijo')
  END FUNCTION rdmft_ls_is_gll
  !
  SUBROUTINE rdmft_resolve_occ_ls_type(ls_type_in, ls_eff)
    !! Resolve the occupation line-search type.  ``auto`` -> monotone
    !! Armijo backtracking (BMR SPG2 default in
    !! :subroutine:`rdmft_spg_occ_block`).  Set ``gll``/``nm_armijo`` for
    !! the GLL nonmonotone reference.  Other choices resolve to
    !! ``armijo`` (monotone backtracking),
    !! ``nm_sw`` (non-monotone strong Wolfe / Zhang-Hager reference),
    !! or ``weak`` (weak-Wolfe variant).
    CHARACTER(LEN=*), INTENT(IN)  :: ls_type_in
    CHARACTER(LEN=24), INTENT(OUT) :: ls_eff
    ls_eff = ADJUSTL(ls_type_in)
    IF (TRIM(ls_eff) == 'auto') ls_eff = 'armijo'
  END SUBROUTINE rdmft_resolve_occ_ls_type
  !
  SUBROUTINE rdmft_resolve_orb_ls_type(ls_type_in, ls_eff, use_wolfe, use_armijo)
    CHARACTER(LEN=*), INTENT(IN)  :: ls_type_in
    CHARACTER(LEN=24), INTENT(OUT) :: ls_eff
    LOGICAL, INTENT(OUT) :: use_wolfe, use_armijo
    ls_eff = ADJUSTL(ls_type_in)
    IF (TRIM(ls_eff) == 'auto') ls_eff = 'strong_wolfe'
    use_wolfe = rdmft_ls_is_wolfe(ls_eff)
    use_armijo = .NOT. use_wolfe
  END SUBROUTINE rdmft_resolve_orb_ls_type
  SUBROUTINE rdmft_armijo_next_alpha(alpha_cur, f0, phi0, f_trial, &
                                     rho_geom, alpha_next, &
                                     sigma_lo, sigma_hi, bmr_halve_outside)
    !! Polynomial (quadratic-fit) refinement of an Armijo backtracking
    !! trial step.  Optional ``sigma_lo`` / ``sigma_hi`` (BMR
    !! \(\sigma_1, \sigma_2\)) replace the default ``[0.1, 0.5]`` band.
    !! When ``bmr_halve_outside`` is true and the unconstrained quadratic
    !! minimiser lies outside that band, the next trial is ``alpha/2``.
    REAL(DP), INTENT(IN)  :: alpha_cur, f0, phi0, f_trial, rho_geom
    REAL(DP), INTENT(OUT) :: alpha_next
    REAL(DP), INTENT(IN), OPTIONAL :: sigma_lo, sigma_hi
    LOGICAL, INTENT(IN), OPTIONAL :: bmr_halve_outside
    REAL(DP) :: numer, denom, alpha_q, scale
    REAL(DP) :: rho_min, rho_max
    LOGICAL :: bmr_halve
    !
    rho_min = 0.1_DP
    rho_max = 0.5_DP
    IF (PRESENT(sigma_lo)) rho_min = sigma_lo
    IF (PRESENT(sigma_hi)) rho_max = sigma_hi
    bmr_halve = .FALSE.
    IF (PRESENT(bmr_halve_outside)) bmr_halve = bmr_halve_outside
    !
    IF (phi0 >= 0.0_DP) THEN
       alpha_next = alpha_cur * rho_geom
       RETURN
    ENDIF
    numer = -phi0 * alpha_cur * alpha_cur
    denom = 2.0_DP * (f_trial - f0 - phi0 * alpha_cur)
    scale = ABS(phi0 * alpha_cur) + ABS(f_trial - f0)
    IF (denom > 1.0e-30_DP .AND. ABS(denom) > 1.0e-14_DP * (1.0_DP + scale)) THEN
       alpha_q = numer / denom
       IF (bmr_halve .AND. &
            (alpha_q < rho_min * alpha_cur .OR. alpha_q > rho_max * alpha_cur)) THEN
          alpha_next = 0.5_DP * alpha_cur
       ELSE
          alpha_next = MIN(MAX(alpha_q, rho_min * alpha_cur), rho_max * alpha_cur)
       ENDIF
    ELSE
       IF (bmr_halve) THEN
          alpha_next = 0.5_DP * alpha_cur
       ELSE
          alpha_next = alpha_cur * rho_geom
       ENDIF
    ENDIF
  END SUBROUTINE rdmft_armijo_next_alpha
  !
  SUBROUTINE rdmft_armijo_ls(eval, f0, phi0, alpha_init, c1, rho, max_iter, &
                             result, post_reject, f_ref, alpha_floor, &
                             allow_stall, use_polynomial, log_tag, &
                             sigma_lo, sigma_hi, bmr_halve_outside)
    !! Monotone Armijo backtracking with optional quadratic step
    !! refinement (:func:`rdmft_armijo_next_alpha`).  Shared by the SPG,
    !! GD, and EBI occupation inner blocks.
    PROCEDURE(rdmft_ls_eval_if) :: eval
    PROCEDURE(rdmft_ls_post_reject_if), OPTIONAL :: post_reject
    REAL(DP), INTENT(IN) :: f0, phi0, alpha_init, c1, rho
    INTEGER, INTENT(IN) :: max_iter
    TYPE(rdmft_ls_result), INTENT(OUT) :: result
    REAL(DP), INTENT(IN), OPTIONAL :: f_ref, alpha_floor
    LOGICAL, INTENT(IN), OPTIONAL :: allow_stall, use_polynomial
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: log_tag
    REAL(DP), INTENT(IN), OPTIONAL :: sigma_lo, sigma_hi
    LOGICAL, INTENT(IN), OPTIONAL :: bmr_halve_outside
    REAL(DP) :: alpha, f_trial, g_trial, f_ref_loc, armijo_rhs, alpha_floor_loc
    INTEGER :: ls, ierr
    LOGICAL :: poly, stall, have_post_reject
    !
    result%alpha_init = alpha_init
    result%step = 0.0_DP
    result%f_new = f0
    result%n_feval = 0
    result%success = .FALSE.
    f_ref_loc = f0
    IF (PRESENT(f_ref)) f_ref_loc = f_ref
    alpha_floor_loc = 0.0_DP
    IF (PRESENT(alpha_floor)) alpha_floor_loc = alpha_floor
    poly = .TRUE.
    IF (PRESENT(use_polynomial)) poly = use_polynomial
    stall = .FALSE.
    IF (PRESENT(allow_stall)) stall = allow_stall
    have_post_reject = PRESENT(post_reject)
    !
    IF (.NOT. (phi0 < 0.0_DP)) RETURN
    IF (alpha_init <= 0.0_DP) RETURN
    !
    alpha = alpha_init
    DO ls = 1, max_iter
       IF (alpha_floor_loc > 0.0_DP .AND. alpha < alpha_floor_loc) EXIT
       CALL eval(alpha, f_trial, g_trial, ierr)
       result%n_feval = result%n_feval + 1
       IF (ierr /= 0) THEN
          IF (have_post_reject) CALL post_reject()
          alpha = alpha * rho
          CYCLE
       ENDIF
       IF (rdmft_verbose >= 2 .AND. ionode) THEN
          IF (PRESENT(log_tag) .AND. LEN_TRIM(log_tag) > 0) THEN
             WRITE(stdout, '(11X,A,I3,A,1PE10.2,0P,A,F18.10,A,F18.10,A,A)') &
                  'Armijo trial ', ls, ': alpha=', alpha, &
                  '  E_trial=', f_trial, '  E_ref=', f_ref_loc, '  ', TRIM(log_tag)
          ELSE
             WRITE(stdout, '(11X,A,I3,A,1PE10.2,0P,A,F18.10,A,F18.10)') &
                  'Armijo trial ', ls, ': alpha=', alpha, &
                  '  E_trial=', f_trial, '  E_ref=', f_ref_loc
          ENDIF
          FLUSH(stdout)
       ENDIF
       armijo_rhs = f_ref_loc + c1 * alpha * phi0
       IF (f_trial <= armijo_rhs) THEN
          result%step = alpha
          result%f_new = f_trial
          result%success = .TRUE.
          RETURN
       ENDIF
       IF (stall .AND. rdmft_ls_stall_accept(f_trial, f_ref_loc)) THEN
          IF (f_trial > f_ref_loc) THEN
             IF (have_post_reject) CALL post_reject()
             result%step = 0.0_DP
             result%f_new = f_ref_loc
          ELSE
             result%step = alpha
             result%f_new = f_trial
          ENDIF
          result%success = .TRUE.
          RETURN
       ENDIF
       IF (have_post_reject) CALL post_reject()
       IF (poly) THEN
          IF (PRESENT(sigma_lo) .OR. PRESENT(sigma_hi) .OR. PRESENT(bmr_halve_outside)) THEN
             ! Quadratic fit uses f(x_k)=f0; acceptance uses f_ref_loc.
             CALL rdmft_armijo_next_alpha(alpha, f0, phi0, f_trial, rho, alpha, &
                  sigma_lo=sigma_lo, sigma_hi=sigma_hi, &
                  bmr_halve_outside=bmr_halve_outside)
          ELSE
             CALL rdmft_armijo_next_alpha(alpha, f0, phi0, f_trial, rho, alpha)
          ENDIF
       ELSE
          alpha = alpha * rho
       ENDIF
    ENDDO
    IF (have_post_reject) CALL post_reject()
  END SUBROUTINE rdmft_armijo_ls
  !
  LOGICAL FUNCTION rdmft_ls_stall_accept(f_trial, f_ref)
    !! True when a line-search trial energy matches ``f_ref`` within
    !! working precision (stationary point / numerical noise floor).
    REAL(DP), INTENT(IN) :: f_trial, f_ref
    rdmft_ls_stall_accept = ABS(f_trial - f_ref) &
         <= rdmft_ls_energy_rtol * (1.0_DP + ABS(f_ref))
  END FUNCTION rdmft_ls_stall_accept
  !
  SUBROUTINE rdmft_ls_log_wolfe_trial(phase, iter, alpha, f, g, tag)
    CHARACTER(LEN=*), INTENT(IN) :: phase
    INTEGER, INTENT(IN) :: iter
    REAL(DP), INTENT(IN) :: alpha, f, g
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: tag
    IF (rdmft_verbose < 2 .OR. .NOT. ionode) RETURN
    IF (PRESENT(tag) .AND. LEN_TRIM(tag) > 0) THEN
       WRITE(stdout, '(11X,A,A,A,I0,A,1PE10.2,0P,A,F18.10,A,1PE10.2,A,A)') &
            'Wolfe ', TRIM(phase), ' ', iter, ': alpha=', alpha, &
            '  E=', f, "  phi'=", g, '  ', TRIM(tag)
    ELSE
       WRITE(stdout, '(11X,A,A,A,I0,A,1PE10.2,0P,A,F18.10,A,1PE10.2)') &
            'Wolfe ', TRIM(phase), ' ', iter, ': alpha=', alpha, &
            '  E=', f, "  phi'=", g
    ENDIF
  END SUBROUTINE rdmft_ls_log_wolfe_trial
  !
  SUBROUTINE rdmft_ls_log_wolfe_summary(alpha_init, step, n_feval, f0, f1, phi0, &
       c1, c2, success, weak_curv, has_f_ref, f_ref)
    REAL(DP), INTENT(IN) :: alpha_init, step, f0, f1, phi0, c1, c2
    INTEGER, INTENT(IN) :: n_feval
    LOGICAL, INTENT(IN) :: success, weak_curv, has_f_ref
    REAL(DP), INTENT(IN) :: f_ref
    CHARACTER(LEN=6) :: stat
    IF (rdmft_verbose < 1 .OR. .NOT. ionode) RETURN
    stat = '  fail'
    IF (success) stat = '  ok  '
    ! NB: the `1P` scale factor introduced by `1PE10.2` (alpha_init, step)
    ! persists for all subsequent F/E descriptors until an explicit `0P`
    ! reset.  Without that reset, the `F18.10` descriptors for f0 and f1
    ! print the values multiplied by 10 -- which is what produced the
    ! mysterious "factor of 10" between the Wolfe summary energies and
    ! the outer RDMFT energy.  Insert `0P` before the F18.10 fields.
    IF (has_f_ref) THEN
       IF (weak_curv) THEN
          WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A,I0,A,0PF18.10,A,F18.10,A,1PE10.2,0P,A,F8.4,A,F8.4,A,F18.10,A,A,A)') &
               'Wolfe LS: alpha_init=', alpha_init, ' step=', step, ' n_feval=', n_feval, &
               ' f0=', f0, ' f1=', f1, " phi0=", phi0, ' c1=', c1, ' c2=', c2, &
               ' f_ref=', f_ref, ' weak', stat
       ELSE
          WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A,I0,A,0PF18.10,A,F18.10,A,1PE10.2,0P,A,F8.4,A,F8.4,A,F18.10,A,A)') &
               'Wolfe LS: alpha_init=', alpha_init, ' step=', step, ' n_feval=', n_feval, &
               ' f0=', f0, ' f1=', f1, " phi0=", phi0, ' c1=', c1, ' c2=', c2, &
               ' f_ref=', f_ref, stat
       ENDIF
    ELSE
       IF (weak_curv) THEN
          WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A,I0,A,0PF18.10,A,F18.10,A,1PE10.2,0P,A,F8.4,A,F8.4,A,A)') &
               'Wolfe LS: alpha_init=', alpha_init, ' step=', step, ' n_feval=', n_feval, &
               ' f0=', f0, ' f1=', f1, " phi0=", phi0, ' c1=', c1, ' c2=', c2, &
               ' weak', stat
       ELSE
          WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A,I0,A,0PF18.10,A,F18.10,A,1PE10.2,0P,A,F8.4,A,F8.4,A,A)') &
               'Wolfe LS: alpha_init=', alpha_init, ' step=', step, ' n_feval=', n_feval, &
               ' f0=', f0, ' f1=', f1, " phi0=", phi0, ' c1=', c1, ' c2=', c2, stat
       ENDIF
    ENDIF
  END SUBROUTINE rdmft_ls_log_wolfe_summary
  !
  SUBROUTINE rdmft_wolfe_zoom(eval, alpha_lo, f_lo, g_lo, alpha_hi, f_hi, g_hi, &
                              f0, g0, f_ref, c1, c2, weak_curv, max_zoom, &
                              alpha_out, f_out, g_out, nfeval, success)
    PROCEDURE(rdmft_ls_eval_if) :: eval
    REAL(DP), INTENT(IN)  :: alpha_lo, f_lo, g_lo, alpha_hi, f_hi, g_hi
    REAL(DP), INTENT(IN)  :: f0, g0, f_ref, c1, c2
    LOGICAL, INTENT(IN)  :: weak_curv
    INTEGER, INTENT(IN)  :: max_zoom
    REAL(DP), INTENT(OUT) :: alpha_out, f_out, g_out
    INTEGER, INTENT(INOUT) :: nfeval
    LOGICAL, INTENT(OUT) :: success
    REAL(DP) :: al_lo, fl_lo, gl_lo, al_hi, fh_hi, gh_hi
    REAL(DP) :: alpha_j, f_j, g_j, alpha_j_prev, width
    REAL(DP) :: best_alpha, best_f, best_g
    REAL(DP), PARAMETER :: stag_rel = 1.0e-10_DP
    INTEGER :: j, ierr
    LOGICAL :: have_best
    !
    success = .FALSE.
    al_lo = alpha_lo
    fl_lo = f_lo
    gl_lo = g_lo
    al_hi = alpha_hi
    fh_hi = f_hi
    gh_hi = g_hi
    alpha_out = al_lo
    f_out = fl_lo
    g_out = gl_lo
    alpha_j_prev = HUGE(1.0_DP)
    ! Track the best Armijo-satisfying trial seen during zoom, so we can
    ! always return a usable step even when the bracket collapses on a
    ! plateau / kinked region and the textbook lo / hi endpoints become
    ! degenerate.
    have_best = .FALSE.
    best_alpha = 0.0_DP
    best_f = HUGE(1.0_DP)
    best_g = 0.0_DP
    IF (al_lo > 0.0_DP .AND. fl_lo <= f_ref + c1 * al_lo * g0) THEN
       have_best = .TRUE.
       best_alpha = al_lo
       best_f = fl_lo
       best_g = gl_lo
    ENDIF
    ! 0P needed before F18.10 to undo the 1P scale set by 1PE10.2 (would
    ! otherwise print f_lo / f_hi multiplied by 10).
    IF (rdmft_verbose >= 2 .AND. ionode) WRITE(stdout, &
         '(11X,A,1PE10.2,A,1PE10.2,A,0PF18.10,A,F18.10)') &
         'Wolfe zoom start: alpha_lo=', al_lo, ' alpha_hi=', al_hi, &
         ' f_lo=', fl_lo, ' f_hi=', fh_hi
    DO j = 1, max_zoom
       ! Safeguarded higher-order interpolation: try cubic Hermite first
       ! (uses (al_lo, fl_lo, gl_lo) + (al_hi, fh_hi, gh_hi)); if the
       ! cubic interpolant is degenerate or its predicted minimum falls
       ! outside the safeguard band, automatically fall back to a 3-pt
       ! quadratic (function values at both endpoints + gradient at one)
       ! and then to a 2-pt quadratic (linear-gradient model).  Bisects
       ! only when all three interpolants reject.  See Nocedal & Wright
       ! Algorithm 3.6 and More-Thuente 1994 for the strategy.
       alpha_j = rdmft_zoom_interp(al_lo, fl_lo, gl_lo, al_hi, fh_hi, gh_hi)
       !
       ! Stagnation guard: when the bracket interval has collapsed (e.g.
       ! the function is on a plateau because a proximal projection has
       ! saturated and (f,g) no longer change with alpha), the cubic
       ! interpolant produces the same trial alpha forever and Wolfe's
       ! curvature condition cannot be satisfied.  In that case, accept
       ! the current lo endpoint as a weak-Wolfe / Armijo step rather
       ! than burning the full max_zoom budget and returning failure.
       width = ABS(al_hi - al_lo)
       IF (width <= stag_rel * MAX(1.0_DP, ABS(al_lo), ABS(al_hi)) .OR. &
           ABS(alpha_j - alpha_j_prev) <= &
              stag_rel * MAX(1.0_DP, ABS(alpha_j), ABS(alpha_j_prev)) .OR. &
           (j >= 8 .AND. have_best)) THEN
          IF (have_best) THEN
             alpha_out = best_alpha
             f_out = best_f
             g_out = best_g
             success = .TRUE.
             CALL rdmft_ls_log_wolfe_trial('zoom', j, best_alpha, best_f, best_g, &
                  'stagnated_accept_best')
          ELSE
             CALL rdmft_ls_log_wolfe_trial('zoom', j, al_lo, fl_lo, gl_lo, &
                  'stagnated_no_armijo')
          ENDIF
          RETURN
       ENDIF
       alpha_j_prev = alpha_j
       CALL eval(alpha_j, f_j, g_j, ierr)
       nfeval = nfeval + 1
       IF (ierr /= 0) THEN
          CALL rdmft_ls_log_wolfe_trial('zoom', j, alpha_j, f_j, g_j, 'eval_err')
          EXIT
       ENDIF
       ! Record alpha_j as a candidate fall-back step if it satisfies
       ! Armijo and improves on our current best.
       IF (alpha_j > 0.0_DP .AND. &
           f_j <= f_ref + c1 * alpha_j * g0 .AND. &
           f_j < best_f) THEN
          have_best = .TRUE.
          best_alpha = alpha_j
          best_f = f_j
          best_g = g_j
       ENDIF
       IF (f_j > f_ref + c1 * alpha_j * g0 .OR. f_j >= fl_lo) THEN
          CALL rdmft_ls_log_wolfe_trial('zoom', j, alpha_j, f_j, g_j, 'shrink_hi')
          al_hi = alpha_j
          fh_hi = f_j
          gh_hi = g_j
       ELSE
          IF (weak_curv) THEN
             IF (g_j >= c2 * g0) THEN
                CALL rdmft_ls_log_wolfe_trial('zoom', j, alpha_j, f_j, g_j, 'curv_ok')
                alpha_out = alpha_j
                f_out = f_j
                g_out = g_j
                success = .TRUE.
                RETURN
             ENDIF
          ELSE
             IF (ABS(g_j) <= c2 * ABS(g0)) THEN
                CALL rdmft_ls_log_wolfe_trial('zoom', j, alpha_j, f_j, g_j, 'curv_ok')
                alpha_out = alpha_j
                f_out = f_j
                g_out = g_j
                success = .TRUE.
                RETURN
             ENDIF
          ENDIF
          IF (g_j * (al_hi - al_lo) >= 0.0_DP) THEN
             CALL rdmft_ls_log_wolfe_trial('zoom', j, alpha_j, f_j, g_j, 'swap_hi')
             al_hi = al_lo
             fh_hi = fl_lo
             gh_hi = gl_lo
          ELSE
             CALL rdmft_ls_log_wolfe_trial('zoom', j, alpha_j, f_j, g_j, 'update_lo')
          ENDIF
          al_lo = alpha_j
          fl_lo = f_j
          gl_lo = g_j
       ENDIF
    ENDDO
    !
    ! Curvature could not be satisfied within max_zoom iterations.  This
    ! is the typical failure mode when the function being searched has a
    ! kink (e.g. a proximal projection that saturates) or a divergent
    ! gradient on one side of the optimum (e.g. the Muller / sqrt(n)
    ! gradient blowing up at clipped occupations): the cubic interpolant
    ! happily converges to the right alpha but |phi'(alpha)| stays large
    ! and Wolfe's curvature test never fires.  In that case, fall back
    ! to the best Armijo-satisfying alpha we recorded during zoom (this
    ! is exactly the step a pure Armijo line search would have accepted)
    ! so the alternating optimiser keeps making progress instead of
    ! stalling.
    IF (have_best) THEN
       alpha_out = best_alpha
       f_out = best_f
       g_out = best_g
       success = .TRUE.
       CALL rdmft_ls_log_wolfe_trial('zoom', max_zoom, best_alpha, best_f, best_g, &
            'maxzoom_armijo_best')
    ENDIF
  END SUBROUTINE rdmft_wolfe_zoom
  !
  SUBROUTINE rdmft_strong_wolfe_ls(eval, f0, g0, alpha_init, c1, c2, max_iter, max_zoom, &
                                   result, f_ref, weak_curv)
    PROCEDURE(rdmft_ls_eval_if) :: eval
    REAL(DP), INTENT(IN) :: f0, g0, alpha_init, c1, c2
    INTEGER, INTENT(IN) :: max_iter, max_zoom
    TYPE(rdmft_ls_result), INTENT(OUT) :: result
    REAL(DP), INTENT(IN), OPTIONAL :: f_ref
    LOGICAL, INTENT(IN), OPTIONAL :: weak_curv
    REAL(DP) :: alpha, alpha_prev, f_prev, g_prev, f_i, g_i, alpha_new
    REAL(DP) :: f_out, g_out, alpha_out, f_ref_loc
    INTEGER :: i, ierr
    LOGICAL :: ok, weak_loc, has_f_ref
    !
    f_ref_loc = f0
    has_f_ref = PRESENT(f_ref)
    IF (has_f_ref) f_ref_loc = f_ref
    weak_loc = .FALSE.
    IF (PRESENT(weak_curv)) weak_loc = weak_curv
    result%alpha_init = alpha_init
    result%n_feval = 0
    result%success = .FALSE.
    IF (.NOT. (g0 < 0.0_DP)) THEN
       result%step = 0.0_DP
       result%f_new = f0
       IF (rdmft_verbose >= 2 .AND. ionode) WRITE(stdout, &
            '(11X,A,1PE10.2,A,A)') &
            'Wolfe bracket: phi0=', g0, '  no descent (g0>=0)'
       CALL rdmft_ls_log_wolfe_summary(alpha_init, result%step, result%n_feval, &
            f0, result%f_new, g0, c1, c2, result%success, weak_loc, has_f_ref, f_ref_loc)
       RETURN
    ENDIF
    alpha = alpha_init
    alpha_prev = 0.0_DP
    f_prev = f0
    g_prev = g0
    DO i = 1, max_iter
       CALL eval(alpha, f_i, g_i, ierr)
       result%n_feval = result%n_feval + 1
       IF (ierr /= 0) THEN
          CALL rdmft_ls_log_wolfe_trial('bracket', i, alpha, f_i, g_i, 'eval_err')
          EXIT
       ENDIF
       IF (f_i > f_ref_loc + c1 * alpha * g0 .OR. (i > 1 .AND. f_i >= f_prev)) THEN
          CALL rdmft_ls_log_wolfe_trial('bracket', i, alpha, f_i, g_i, 'zoom')
          CALL rdmft_wolfe_zoom(eval, alpha_prev, f_prev, g_prev, alpha, f_i, g_i, &
               f0, g0, f_ref_loc, c1, c2, weak_loc, max_zoom, &
               alpha_out, f_out, g_out, result%n_feval, ok)
          IF (ok) THEN
             result%step = alpha_out
             result%f_new = f_out
             result%success = .TRUE.
             CALL rdmft_ls_log_wolfe_summary(alpha_init, result%step, result%n_feval, &
                  f0, result%f_new, g0, c1, c2, result%success, weak_loc, has_f_ref, f_ref_loc)
             RETURN
          ENDIF
          result%step = alpha_prev
          result%f_new = f_prev
          CALL rdmft_ls_log_wolfe_summary(alpha_init, result%step, result%n_feval, &
               f0, result%f_new, g0, c1, c2, result%success, weak_loc, has_f_ref, f_ref_loc)
          RETURN
       ENDIF
       IF (weak_loc) THEN
          IF (g_i >= c2 * g0) THEN
             CALL rdmft_ls_log_wolfe_trial('bracket', i, alpha, f_i, g_i, 'curv_ok')
             result%step = alpha
             result%f_new = f_i
             result%success = .TRUE.
             CALL rdmft_ls_log_wolfe_summary(alpha_init, result%step, result%n_feval, &
                  f0, result%f_new, g0, c1, c2, result%success, weak_loc, has_f_ref, f_ref_loc)
             RETURN
          ENDIF
       ELSE
          IF (ABS(g_i) <= c2 * ABS(g0)) THEN
             CALL rdmft_ls_log_wolfe_trial('bracket', i, alpha, f_i, g_i, 'curv_ok')
             result%step = alpha
             result%f_new = f_i
             result%success = .TRUE.
             CALL rdmft_ls_log_wolfe_summary(alpha_init, result%step, result%n_feval, &
                  f0, result%f_new, g0, c1, c2, result%success, weak_loc, has_f_ref, f_ref_loc)
             RETURN
          ENDIF
       ENDIF
       IF (g_i >= 0.0_DP) THEN
          CALL rdmft_ls_log_wolfe_trial('bracket', i, alpha, f_i, g_i, 'zoom_sign')
          CALL rdmft_wolfe_zoom(eval, alpha, f_i, g_i, alpha_prev, f_prev, g_prev, &
               f0, g0, f_ref_loc, c1, c2, weak_loc, max_zoom, &
               alpha_out, f_out, g_out, result%n_feval, ok)
          IF (ok) THEN
             result%step = alpha_out
             result%f_new = f_out
             result%success = .TRUE.
             CALL rdmft_ls_log_wolfe_summary(alpha_init, result%step, result%n_feval, &
                  f0, result%f_new, g0, c1, c2, result%success, weak_loc, has_f_ref, f_ref_loc)
             RETURN
          ENDIF
          result%step = alpha_prev
          result%f_new = f_prev
          CALL rdmft_ls_log_wolfe_summary(alpha_init, result%step, result%n_feval, &
               f0, result%f_new, g0, c1, c2, result%success, weak_loc, has_f_ref, f_ref_loc)
          RETURN
       ENDIF
       CALL rdmft_ls_log_wolfe_trial('bracket', i, alpha, f_i, g_i, 'double')
       alpha_new = 2.0_DP * alpha
       IF (alpha_new <= alpha) THEN
          result%step = alpha
          result%f_new = f_i
          CALL rdmft_ls_log_wolfe_summary(alpha_init, result%step, result%n_feval, &
               f0, result%f_new, g0, c1, c2, result%success, weak_loc, has_f_ref, f_ref_loc)
          RETURN
       ENDIF
       alpha_prev = alpha
       f_prev = f_i
       g_prev = g_i
       alpha = alpha_new
    ENDDO
    result%step = alpha_prev
    result%f_new = f_prev
    CALL rdmft_ls_log_wolfe_summary(alpha_init, result%step, result%n_feval, &
         f0, result%f_new, g0, c1, c2, result%success, weak_loc, has_f_ref, f_ref_loc)
  END SUBROUTINE rdmft_strong_wolfe_ls
  !
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_zhang_hager_init(state, etot0, eta)
    TYPE(rdmft_zhang_hager_state), INTENT(OUT) :: state
    REAL(DP), INTENT(IN) :: etot0
    REAL(DP), INTENT(IN), OPTIONAL :: eta
    state%Q = 1.0_DP
    state%C = etot0
    IF (PRESENT(eta)) state%eta = MAX(0.0_DP, MIN(1.0_DP, eta))
    state%initialised = .TRUE.
  END SUBROUTINE rdmft_zhang_hager_init
  !
  SUBROUTINE rdmft_zhang_hager_update(state, etot_new)
    TYPE(rdmft_zhang_hager_state), INTENT(INOUT) :: state
    REAL(DP), INTENT(IN) :: etot_new
    REAL(DP) :: Qnew
    IF (.NOT. state%initialised) THEN
       CALL rdmft_zhang_hager_init(state, etot_new)
       RETURN
    ENDIF
    Qnew    = state%eta * state%Q + 1.0_DP
    state%C = (state%eta * state%Q * state%C + etot_new) / Qnew
    state%Q = Qnew
  END SUBROUTINE rdmft_zhang_hager_update
  !
  REAL(DP) FUNCTION rdmft_zhang_hager_ref(state)
    TYPE(rdmft_zhang_hager_state), INTENT(IN) :: state
    rdmft_zhang_hager_ref = state%C
  END FUNCTION rdmft_zhang_hager_ref
  !
  SUBROUTINE rdmft_zhang_hager_reset(state)
    TYPE(rdmft_zhang_hager_state), INTENT(INOUT) :: state
    state%Q = 1.0_DP
    state%C = 0.0_DP
    state%initialised = .FALSE.
  END SUBROUTINE rdmft_zhang_hager_reset
  !
  SUBROUTINE rdmft_gll_init(state, etot0, memory)
    TYPE(rdmft_gll_state), INTENT(OUT) :: state
    REAL(DP), INTENT(IN) :: etot0
    INTEGER, INTENT(IN), OPTIONAL :: memory
    INTEGER :: mem
    !
    mem = 10
    IF (PRESENT(memory)) mem = MAX(1, memory)
    state%memory = mem
    IF (ALLOCATED(state%f_hist)) DEALLOCATE(state%f_hist)
    ALLOCATE(state%f_hist(mem))
    state%f_hist(1) = etot0
    state%count = 1
    state%head = 1
    state%initialised = .TRUE.
  END SUBROUTINE rdmft_gll_init
  !
  SUBROUTINE rdmft_gll_push(state, etot_new)
    TYPE(rdmft_gll_state), INTENT(INOUT) :: state
    REAL(DP), INTENT(IN) :: etot_new
    INTEGER :: next
    !
    IF (.NOT. state%initialised) THEN
       CALL rdmft_gll_init(state, etot_new)
       RETURN
    ENDIF
    IF (.NOT. ALLOCATED(state%f_hist)) THEN
       CALL rdmft_gll_init(state, etot_new, memory=state%memory)
       RETURN
    ENDIF
    IF (state%count < state%memory) THEN
       state%count = state%count + 1
       state%head = state%count
       state%f_hist(state%head) = etot_new
    ELSE
       next = MOD(state%head, state%memory) + 1
       state%head = next
       state%f_hist(next) = etot_new
    ENDIF
  END SUBROUTINE rdmft_gll_push
  !
  REAL(DP) FUNCTION rdmft_gll_ref(state)
    TYPE(rdmft_gll_state), INTENT(IN) :: state
    !
    rdmft_gll_ref = 0.0_DP
    IF (.NOT. state%initialised) RETURN
    IF (.NOT. ALLOCATED(state%f_hist)) RETURN
    IF (state%count < 1) RETURN
    rdmft_gll_ref = MAXVAL(state%f_hist(1:state%count))
  END FUNCTION rdmft_gll_ref
  !
  SUBROUTINE rdmft_gll_reset(state)
    TYPE(rdmft_gll_state), INTENT(INOUT) :: state
    IF (ALLOCATED(state%f_hist)) DEALLOCATE(state%f_hist)
    state%count = 0
    state%head = 0
    state%initialised = .FALSE.
  END SUBROUTINE rdmft_gll_reset
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_cubic_min(a, fa, ga, b, fb, gb)
    !---------------------------------------------------------------
    !! Safeguarded cubic-Hermite minimiser for the interval
    !! \([\min(a,b), \max(a,b)]\) given function values and
    !! derivatives at both endpoints.  Falls back to the bracket
    !! midpoint when the cubic is degenerate (negative discriminant,
    !! zero leading coefficient, ...).
    !!
    !! Kept for backward compatibility -- the zoom phase of the
    !! strong-Wolfe line search now calls :function:`rdmft_zoom_interp`,
    !! which adds quadratic fall-backs when the cubic step is invalid.
    REAL(DP), INTENT(IN) :: a, fa, ga, b, fb, gb
    REAL(DP) :: d1, d2, d2_sq, alpha_star, lo, hi, margin, h, sgn, denom
    IF (ABS(b - a) <= 1.0e-30_DP * (ABS(a) + ABS(b) + 1.0_DP)) THEN
       rdmft_cubic_min = 0.5_DP * (a + b)
       RETURN
    ENDIF
    h = b - a
    sgn = SIGN(1.0_DP, h)
    d1 = ga + gb - 3.0_DP * (fb - fa) / h
    d2_sq = d1 * d1 - ga * gb
    IF (d2_sq < 0.0_DP) THEN
       rdmft_cubic_min = 0.5_DP * (a + b)
       RETURN
    ENDIF
    ! sign-correct d2 (Nocedal & Wright Eq 3.59); without this the
    ! formula returns the local *maximum* of the cubic when b < a,
    ! which the zoom loop never wants.
    d2 = sgn * SQRT(d2_sq)
    denom = gb - ga + 2.0_DP * d2
    IF (ABS(denom) <= 1.0e-30_DP * (ABS(gb) + ABS(ga) + ABS(d2) + 1.0_DP)) THEN
       rdmft_cubic_min = 0.5_DP * (a + b)
       RETURN
    ENDIF
    alpha_star = b - h * (gb + d2 - d1) / denom
    lo = MIN(a, b)
    hi = MAX(a, b)
    margin = 0.1_DP * (hi - lo)
    rdmft_cubic_min = MIN(MAX(alpha_star, lo + margin), hi - margin)
  END FUNCTION rdmft_cubic_min
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_zoom_interp(a, fa, ga, b, fb, gb)
    !---------------------------------------------------------------
    !! Safeguarded zoom interpolation with automatic cubic /
    !! quadratic / bisection fall-back, following Nocedal \& Wright
    !! Algorithm 3.6 and the safeguard strategy used by More-Thuente
    !! 1994 (\textit{ACM Trans.\ Math.\ Software} \textbf{20}, 286).
    !!
    !! Given the bracket endpoints ``(a, fa, ga)`` and ``(b, fb, gb)``
    !! the routine tries three candidate interpolants in order and
    !! returns the first one whose predicted minimiser lies strictly
    !! inside the safeguard band ``[lo + margin, hi - margin]``
    !! with ``margin = 0.1 * (hi - lo)``:
    !!
    !!   1. **Cubic Hermite** -- the unique cubic through the four
    !!      pieces of information ``{fa, ga, fb, gb}``.  Used when
    !!      sufficient curvature is present (positive discriminant
    !!      and non-zero denominator) and the minimum is well inside
    !!      the bracket.  This is the most accurate choice when the
    !!      objective is locally smooth.
    !!   2. **Quadratic-3pt** -- the unique quadratic through
    !!      ``(a, fa, ga)`` and ``(b, fb)``; minimum at
    !!      ``a - ga / (2 c)`` with ``c = (fb - fa - ga (b-a)) / (b-a)^2``.
    !!      Used when the cubic's curvature information is unreliable
    !!      (e.g.\ ``ga`` and ``gb`` nearly equal -- the cubic
    !!      becomes ill-conditioned).  Requires positive curvature
    !!      (``c > 0``).
    !!   3. **Quadratic-2pt** -- linear-gradient model, minimum at
    !!      ``a - ga * (b-a) / (gb - ga)``.  Used when neither the
    !!      cubic nor the 3-pt quadratic gives a usable step;
    !!      requires ``gb > ga`` (positive curvature).
    !!
    !! When all three reject, the routine bisects.  All candidates
    !! are clamped to the safeguard band on output so an aggressive
    !! step never collapses the bracket to either endpoint.
    !!
    !! With ``rdmft_verbose >= 3`` the ionode reports which branch
    !! produced the accepted step (``cubic`` / ``quad3`` / ``quad2``
    !! / ``bisect``); useful when tuning convergence of the
    !! occupation / orbital line searches.
    REAL(DP), INTENT(IN) :: a, fa, ga, b, fb, gb
    REAL(DP) :: h, lo, hi, margin, alpha
    REAL(DP) :: d1, d2_sq, d2, sgn, denom_c, alpha_cub
    REAL(DP) :: c_q3, alpha_q3, alpha_q2, dg
    REAL(DP), PARAMETER :: tiny_h    = 1.0e-30_DP
    REAL(DP), PARAMETER :: tiny_pos  = 1.0e-30_DP
    LOGICAL :: ok_cub, ok_q3, ok_q2
    CHARACTER(LEN=6) :: branch
    !
    h = b - a
    IF (ABS(h) <= tiny_h * (ABS(a) + ABS(b) + 1.0_DP)) THEN
       rdmft_zoom_interp = 0.5_DP * (a + b)
       IF (rdmft_verbose >= 3 .AND. ionode) WRITE(stdout, '(13X,A)') &
            'zoom_interp: degenerate bracket; midpoint.'
       RETURN
    ENDIF
    lo = MIN(a, b)
    hi = MAX(a, b)
    margin = 0.1_DP * (hi - lo)
    sgn = SIGN(1.0_DP, h)
    !
    ! Candidate 1: cubic Hermite interpolant.
    d1 = ga + gb - 3.0_DP * (fb - fa) / h
    d2_sq = d1 * d1 - ga * gb
    ok_cub = .FALSE.
    alpha_cub = 0.5_DP * (a + b)
    IF (d2_sq >= 0.0_DP) THEN
       d2 = sgn * SQRT(d2_sq)
       denom_c = gb - ga + 2.0_DP * d2
       IF (ABS(denom_c) > tiny_pos * (ABS(gb) + ABS(ga) + ABS(d2) + 1.0_DP)) THEN
          alpha_cub = b - h * (gb + d2 - d1) / denom_c
          ok_cub = (alpha_cub > lo + margin .AND. alpha_cub < hi - margin)
       ENDIF
    ENDIF
    !
    ! Candidate 2: quadratic through (a, fa, ga) and (b, fb).
    !   q(x) = fa + ga*(x-a) + c_q3*(x-a)^2;   c_q3 = (fb - fa - ga*h)/h^2.
    !   Minimum at  x = a - ga / (2 c_q3),  valid only when c_q3 > 0.
    c_q3 = (fb - fa - ga * h) / (h * h)
    ok_q3 = .FALSE.
    alpha_q3 = 0.5_DP * (a + b)
    IF (c_q3 > tiny_pos * (ABS(ga / MAX(ABS(h), tiny_h)) + 1.0_DP)) THEN
       alpha_q3 = a - 0.5_DP * ga / c_q3
       ok_q3 = (alpha_q3 > lo + margin .AND. alpha_q3 < hi - margin)
    ENDIF
    !
    ! Candidate 3: quadratic from linear gradient (gb - ga)/h.
    !   g(x) = ga + (gb - ga)*(x-a)/h;  zero at x = a - ga*h/(gb-ga).
    !   Valid only when gb > ga (positive secant curvature).
    dg = gb - ga
    ok_q2 = .FALSE.
    alpha_q2 = 0.5_DP * (a + b)
    IF (dg > tiny_pos * (ABS(ga) + ABS(gb) + 1.0_DP)) THEN
       alpha_q2 = a - ga * h / dg
       ok_q2 = (alpha_q2 > lo + margin .AND. alpha_q2 < hi - margin)
    ENDIF
    !
    ! Selection: cubic first (highest order), then quad3 (function
    ! values across the bracket), then quad2 (linear gradient), then
    ! bisection.  We keep the same safeguard band for all three so a
    ! "good cubic that lands on the margin" is not preferred over a
    ! "well-centered quadratic".
    IF (ok_cub) THEN
       alpha = alpha_cub
       branch = 'cubic '
    ELSE IF (ok_q3) THEN
       alpha = alpha_q3
       branch = 'quad3 '
    ELSE IF (ok_q2) THEN
       alpha = alpha_q2
       branch = 'quad2 '
    ELSE
       alpha = 0.5_DP * (a + b)
       branch = 'bisect'
    ENDIF
    rdmft_zoom_interp = MIN(MAX(alpha, lo + margin), hi - margin)
    !
    IF (rdmft_verbose >= 3 .AND. ionode) WRITE(stdout, &
         '(13X,A,A,A,1PE12.4,A,1PE12.4,A,1PE12.4)') &
         'zoom_interp ', TRIM(branch), ': alpha=', rdmft_zoom_interp, &
         '  bracket=[', lo, ', ', hi
  END FUNCTION rdmft_zoom_interp
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_bb_init(bb, ndim, alpha_min, alpha_max)
    TYPE(rdmft_bb_state), INTENT(OUT) :: bb
    INTEGER, INTENT(IN) :: ndim
    REAL(DP), INTENT(IN), OPTIONAL :: alpha_min, alpha_max
    CALL rdmft_bb_reset(bb)
    ALLOCATE(bb%x_prev(ndim), bb%g_prev(ndim))
    IF (PRESENT(alpha_min)) bb%alpha_min = alpha_min
    IF (PRESENT(alpha_max)) bb%alpha_max = alpha_max
  END SUBROUTINE rdmft_bb_init
  !
  SUBROUTINE rdmft_bb_reset(bb)
    TYPE(rdmft_bb_state), INTENT(INOUT) :: bb
    bb%have_prev = .FALSE.
    IF (ALLOCATED(bb%x_prev)) bb%x_prev = 0.0_DP
    IF (ALLOCATED(bb%g_prev)) bb%g_prev = 0.0_DP
  END SUBROUTINE rdmft_bb_reset
  !
  SUBROUTINE rdmft_bb_record(bb, x, g)
    TYPE(rdmft_bb_state), INTENT(INOUT) :: bb
    REAL(DP), INTENT(IN) :: x(:), g(:)
    IF (.NOT. ALLOCATED(bb%x_prev) .OR. .NOT. ALLOCATED(bb%g_prev)) RETURN
    IF (SIZE(x) /= SIZE(bb%x_prev) .OR. SIZE(g) /= SIZE(bb%g_prev)) RETURN
    bb%x_prev = x
    bb%g_prev = g
    bb%have_prev = .TRUE.
  END SUBROUTINE rdmft_bb_record
  !
  REAL(DP) FUNCTION rdmft_bb_suggest(bb, x, g, fallback)
    !---------------------------------------------------------------
    !! Barzilai--Borwein **BB1** steplength in the **Euclidean** metric
    !! (line-search initial trial only):
    !!
    !!   s = x - x_{\mathrm{prev}},\quad
    !!   y = g - g_{\mathrm{prev}},\quad
    !!   \alpha = (s^{\mathsf T}s) / (s^{\mathsf T}y).
    !!
    !! Used by ``rdmft_occ_ls_alpha0`` for Armijo/Wolfe backtracking along
    !! a fixed search direction.
    TYPE(rdmft_bb_state), INTENT(IN) :: bb
    REAL(DP), INTENT(IN) :: x(:), g(:), fallback
    REAL(DP), ALLOCATABLE :: xprev(:), gprev(:)
    REAL(DP) :: sTy, sTs, alpha, si
    INTEGER :: i, n
    !
    rdmft_bb_suggest = fallback
    IF (.NOT. bb%have_prev) RETURN
    IF (.NOT. ALLOCATED(bb%x_prev) .OR. .NOT. ALLOCATED(bb%g_prev)) RETURN
    xprev = bb%x_prev
    gprev = bb%g_prev
    n = SIZE(x)
    IF (n < 1) RETURN
    IF (n /= SIZE(g) .OR. n /= SIZE(xprev) .OR. n /= SIZE(gprev)) RETURN
    sTs = 0.0_DP
    sTy = 0.0_DP
    DO i = 1, n
       si = x(i) - xprev(i)
       sTs = sTs + si * (g(i) - gprev(i))
       sTy = sTy + si * si
    ENDDO
    ! The Barzilai-Borwein step is built from the occupation vector,
    ! which is distributed across k-point pools; reduce the two inner
    ! products so every pool proposes the SAME initial step (otherwise
    ! the first line-search trial differs per pool and desynchronises
    ! the collective trial-energy evaluations).
    CALL mp_sum(sTs, inter_pool_comm)
    CALL mp_sum(sTy, inter_pool_comm)
    IF (ABS(sTs) > 1.0e-30_DP) THEN
       alpha = sTy / sTs
       IF (alpha > 0.0_DP) THEN
          rdmft_bb_suggest = MIN(bb%alpha_max, MAX(bb%alpha_min, alpha))
       ENDIF
    ENDIF
  END FUNCTION rdmft_bb_suggest
  !
  SUBROUTINE rdmft_quad_ls_reset(qhist)
    TYPE(rdmft_quad_ls_history), INTENT(OUT) :: qhist
    qhist%valid = .FALSE.
  END SUBROUTINE rdmft_quad_ls_reset
  !
  SUBROUTINE rdmft_quad_ls_update(qhist, alpha, f0, f1, dd0)
    TYPE(rdmft_quad_ls_history), INTENT(OUT) :: qhist
    REAL(DP), INTENT(IN) :: alpha, f0, f1, dd0
    qhist%valid = .TRUE.
    qhist%alpha = alpha
    qhist%f0 = f0
    qhist%f1 = f1
    qhist%dd0 = dd0
  END SUBROUTINE rdmft_quad_ls_update
  !
  REAL(DP) FUNCTION rdmft_occ_ls_alpha0(init_mode, occ_flat, grad_flat, fallback, &
                                         bb, qhist, bb_alpha_min, bb_alpha_max)
    !---------------------------------------------------------------
    !! Initial line-search trial ``alpha_0`` for occupation / packed
    !! vectors.  The ``barzilai_borwein`` branch uses **Euclidean BB1**
    !! via ``rdmft_bb_suggest`` (Euclidean BB1; independent of BZ weights).
    CHARACTER(LEN=*), INTENT(IN) :: init_mode
    REAL(DP), INTENT(IN) :: occ_flat(:), grad_flat(:), fallback
    TYPE(rdmft_bb_state), INTENT(IN) :: bb
    TYPE(rdmft_quad_ls_history), INTENT(IN) :: qhist
    REAL(DP), INTENT(IN), OPTIONAL :: bb_alpha_min, bb_alpha_max
    REAL(DP) :: denom, a
    !
    rdmft_occ_ls_alpha0 = fallback
    IF (TRIM(init_mode) == 'fixed') THEN
       rdmft_occ_ls_alpha0 = fallback
       RETURN
    ENDIF
    IF (TRIM(init_mode) == 'barzilai_borwein' .OR. TRIM(init_mode) == 'bb') THEN
       IF (PRESENT(bb_alpha_min) .AND. PRESENT(bb_alpha_max)) THEN
          rdmft_occ_ls_alpha0 = rdmft_bb_suggest(bb, occ_flat, grad_flat, fallback)
          rdmft_occ_ls_alpha0 = MIN(bb_alpha_max, MAX(bb_alpha_min, rdmft_occ_ls_alpha0))
       ELSE
          rdmft_occ_ls_alpha0 = rdmft_bb_suggest(bb, occ_flat, grad_flat, fallback)
       ENDIF
       RETURN
    ENDIF
    IF (TRIM(init_mode) == 'quadratic') THEN
       IF (qhist%valid .AND. qhist%alpha > 0.0_DP .AND. qhist%dd0 < 0.0_DP) THEN
          denom = 2.0_DP * (qhist%f1 - qhist%f0 - qhist%dd0 * qhist%alpha)
          IF (ABS(denom) > 1.0e-16_DP * (1.0_DP + ABS(qhist%dd0 * qhist%alpha))) THEN
             a = -qhist%dd0 * qhist%alpha * qhist%alpha / denom
             IF (a > 0.0_DP) THEN
                IF (PRESENT(bb_alpha_min) .AND. PRESENT(bb_alpha_max)) THEN
                   rdmft_occ_ls_alpha0 = MIN(bb_alpha_max, MAX(bb_alpha_min, a))
                ELSE
                   rdmft_occ_ls_alpha0 = a
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    ENDIF
  END FUNCTION rdmft_occ_ls_alpha0
  !
END MODULE rdmft_linesearch
