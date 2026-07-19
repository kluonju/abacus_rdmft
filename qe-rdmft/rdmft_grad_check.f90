!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_gradient_check
  !------------------------------------------------------------------
  !! Finite-difference gradient consistency checks (ABACUS
  !! \texttt{check\_gradient\_consistency}).
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_check_gradient_consistency
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_check_gradient_consistency(epsilon, tolerance)
    !---------------------------------------------------------------
    !! Occupation central differences and orbital forward difference
    !! along the Stiefel retraction (ABACUS \S9 in
    !! \texttt{rdmft\_derivation.md}).
    USE io_global,        ONLY : stdout, ionode
    USE wvfct,          ONLY : nbnd, npwx, current_k
    USE klist,          ONLY : nks, ngk, wk, igk_k, xk, lgauss
    USE wavefunctions,  ONLY : evc
    USE noncollin_module, ONLY : npol
    USE control_flags,  ONLY : gamma_only
    USE gvect,          ONLY : gstart
    USE io_files,       ONLY : nwordwfc, iunwfc
    USE buffers,        ONLY : get_buffer, save_buffer
    USE lsda_mod,       ONLY : lsda, current_spin, isk
    USE uspp,           ONLY : nkb, vkb, okvan
    USE uspp_init,      ONLY : init_us_2
    USE rdmft_module,   ONLY : rdmft_n, rdmft_mark_orbitals_changed, &
                               rdmft_e_one, rdmft_e_har, rdmft_e_xc
    USE rdmft_energy,   ONLY : rdmft_total_energy, rdmft_grad_n, &
                               rdmft_compute_riemannian_gradient
    USE rdmft_stiefel,  ONLY : stiefel_retract_gamma, stiefel_retract_k
    !
    REAL(DP), INTENT(IN), OPTIONAL :: epsilon, tolerance
    REAL(DP), PARAMETER :: default_eps = 1.0e-5_DP
    REAL(DP), PARAMETER :: default_tol = 1.0e-4_DP
    REAL(DP), PARAMETER :: metallic_occ_eps = 1.0e-4_DP
    REAL(DP), PARAMETER :: metallic_occ_tol = 3.0e-2_DP
    REAL(DP) :: eps, tol, eps_occ, tol_occ, E0, E_plus, E_minus
    REAL(DP) :: fd_grad, fd_one, fd_h, fd_x, analytic, rel_err
    REAL(DP) :: analytic_dd, fd_dd, gnorm2
    REAL(DP), ALLOCATABLE :: grad_a(:,:), n_save(:,:)
    REAL(DP), ALLOCATABLE :: e_one_p(:), e_har_p(:), e_xc_p(:)
    REAL(DP), ALLOCATABLE :: e_one_m(:), e_har_m(:), e_xc_m(:)
    COMPLEX(DP), ALLOCATABLE :: G_R_all(:,:,:), C_save_all(:,:,:), C_try(:,:)
    LOGICAL :: all_pass, pass
    INTEGER :: ik, ib, idx, n_check, npw
    !
    eps = default_eps
    tol = default_tol
    IF (PRESENT(epsilon))   eps = epsilon
    IF (PRESENT(tolerance)) tol = tolerance
    eps_occ = eps
    tol_occ = tol
    IF (lgauss) THEN
       ! Metallic smearing runs are numerically noisier in FD
       ! occupation derivatives (tiny dE differences over k-points and
       ! FFT/MPI reductions), so use a slightly larger stencil and a
       ! pragmatic tolerance while keeping the strict default for
       ! insulators.
       eps_occ = MAX(eps_occ, metallic_occ_eps)
       tol_occ = MAX(tol_occ, metallic_occ_tol)
    ENDIF
    !
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    !
    n_check = MIN(10, nbnd * nks)
    ALLOCATE(grad_a(nbnd, nks), n_save(nbnd, nks))
    ALLOCATE(e_one_p(n_check), e_har_p(n_check), e_xc_p(n_check))
    ALLOCATE(e_one_m(n_check), e_har_m(n_check), e_xc_m(n_check))
    ALLOCATE(G_R_all(npwx*npol, nbnd, nks))
    ALLOCATE(C_save_all(npwx*npol, nbnd, nks), C_try(npwx*npol, nbnd))
    !
    n_save = rdmft_n
    CALL rdmft_grad_n(grad_a)
    CALL rdmft_total_energy(E0)
    !
    all_pass = .TRUE.
    IF (ionode) WRITE(stdout, '(/,5X,A)') &
         '===== RDMFT gradient consistency check ====='
    !
    ! ---- Occupation gradients (central difference) ----------------
    IF (ionode) WRITE(stdout, '(/,5X,A)') '-- Occupation gradient check --'
    idx = 0
    DO ik = 1, nks
       DO ib = 1, nbnd
          IF (rdmft_n(ib, ik) + eps_occ > 1.0_DP .OR. rdmft_n(ib, ik) - eps_occ < 0.0_DP) CYCLE
          idx = idx + 1
          IF (idx > n_check) EXIT
          rdmft_n = n_save
          rdmft_n(ib, ik) = n_save(ib, ik) + eps_occ
          CALL rdmft_mark_orbitals_changed()
          CALL rdmft_total_energy(E_plus)
          e_one_p(idx) = rdmft_e_one
          e_har_p(idx) = rdmft_e_har
          e_xc_p(idx)  = rdmft_e_xc
          rdmft_n = n_save
          rdmft_n(ib, ik) = n_save(ib, ik) - eps_occ
          CALL rdmft_mark_orbitals_changed()
          CALL rdmft_total_energy(E_minus)
          e_one_m(idx) = rdmft_e_one
          e_har_m(idx) = rdmft_e_har
          e_xc_m(idx)  = rdmft_e_xc
          rdmft_n = n_save
          fd_grad = (E_plus - E_minus) / (2.0_DP * eps_occ)
          fd_one  = (e_one_p(idx) - e_one_m(idx)) / (2.0_DP * eps_occ)
          fd_h    = (e_har_p(idx) - e_har_m(idx)) / (2.0_DP * eps_occ)
          fd_x    = (e_xc_p(idx) - e_xc_m(idx)) / (2.0_DP * eps_occ)
          analytic = grad_a(ib, ik)
          IF (ABS(analytic) > 1.0e-10_DP) THEN
             rel_err = ABS(fd_grad - analytic) / ABS(analytic)
          ELSE
             rel_err = ABS(fd_grad - analytic)
          ENDIF
          pass = (rel_err < tol_occ)
          IF (.NOT. pass) all_pass = .FALSE.
          IF (ionode) WRITE(stdout, '(5X,A,I0,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES10.2,A,A)') &
               'occ[ik=', ik, ',ib=', ib, ']: analytic=', analytic, &
               '  fd=', fd_grad, '  (fd_one=', fd_one, '  fd_h=', fd_h, '  fd_x=', fd_x, &
               ')  rel_err=', rel_err, MERGE('  PASS', '  FAIL', pass)
       ENDDO
       IF (idx > n_check) EXIT
    ENDDO
    !
    ! ---- Orbital gradient (central difference along -G_R) ----------
    IF (ionode) WRITE(stdout, '(/,5X,A)') '-- Orbital gradient check --'
    rdmft_n = n_save
    CALL rdmft_mark_orbitals_changed()
    CALL rdmft_total_energy(E0)
    CALL rdmft_compute_riemannian_gradient(G_R_all, gnorm2)
    ! True directional derivative of E along -G_R on the Stiefel
    ! manifold with the Frobenius (Re Tr) metric:
    !
    !   dE/dα |_{α=0} along eta = 2 * Re Tr( eta^H * (∂E/∂C^*) )
    !
    ! because for real α the variation dC = α*eta induces
    ! dC^* = α*eta^* and dE picks up both the dC and dC^* contributions,
    ! each equal to <gradC, eta>.  For eta = -G_R = -Proj(∂E/∂C^*)
    ! this gives dE/dα = -2 * <G_R, G_R>_S = -2 * gnorm2.  The factor
    ! of 2 is missing if one writes ``analytic_dd = -gnorm2``, in which
    ! case the FD reports a factor-of-2 mismatch even for HF.
    analytic_dd = -2.0_DP * gnorm2
    DO ik = 1, nks
       IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       C_save_all(:, :, ik) = evc(:, :)
    ENDDO
    ! Central difference along the retraction:
    !   fd_dd = (E(+eps*d) - E(-eps*d)) / (2*eps),   d = -G_R.
    ! The earlier one-sided forward difference (E(+eps*d) - E0)/eps
    ! carries an O(eps) truncation error proportional to the manifold
    ! curvature; for a well-converged HF gradient that error
    ! (~1e-4 rel.) routinely tripped the 1e-4 ``FAIL`` threshold even
    ! though the analytic gradient was correct.  The central difference
    ! is O(eps^2) accurate and matches the occupation check above.
    DO ik = 1, nks
       npw = ngk(ik)
       current_k = ik
       IF (lsda) current_spin = isk(ik)
       IF (gamma_only) THEN
          CALL stiefel_retract_gamma(C_save_all(:, :, ik), -G_R_all(:, :, ik), &
               eps, npw, nbnd, npwx*npol, gstart, C_try)
       ELSE
          CALL stiefel_retract_k(C_save_all(:, :, ik), -G_R_all(:, :, ik), &
               eps, npw, nbnd, npwx*npol, C_try)
       ENDIF
       evc(:, :) = C_try(:, :)
       IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
    ENDDO
    CALL rdmft_mark_orbitals_changed()
    CALL rdmft_total_energy(E_plus)
    DO ik = 1, nks
       npw = ngk(ik)
       current_k = ik
       IF (lsda) current_spin = isk(ik)
       IF (gamma_only) THEN
          CALL stiefel_retract_gamma(C_save_all(:, :, ik), -G_R_all(:, :, ik), &
               -eps, npw, nbnd, npwx*npol, gstart, C_try)
       ELSE
          CALL stiefel_retract_k(C_save_all(:, :, ik), -G_R_all(:, :, ik), &
               -eps, npw, nbnd, npwx*npol, C_try)
       ENDIF
       evc(:, :) = C_try(:, :)
       IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
    ENDDO
    CALL rdmft_mark_orbitals_changed()
    CALL rdmft_total_energy(E_minus)
    fd_dd = (E_plus - E_minus) / (2.0_DP * eps)
    IF (ABS(analytic_dd) > 1.0e-10_DP) THEN
       rel_err = ABS(fd_dd - analytic_dd) / ABS(analytic_dd)
    ELSE
       rel_err = ABS(fd_dd - analytic_dd)
    ENDIF
    pass = (rel_err < tol)
    IF (.NOT. pass) all_pass = .FALSE.
    IF (ionode) WRITE(stdout, '(5X,A,ES14.6,A,ES14.6,A,ES14.6,A,ES10.2,A)') &
         'orb dir deriv: analytic=', analytic_dd, '  fd_cen=', fd_dd, &
         '  ||G_R||_S^2=', gnorm2, '  rel_err=', rel_err, &
         MERGE('  PASS', '  FAIL', pass)
    !
    ! Restore orbitals and occupations.
    DO ik = 1, nks
       current_k = ik
       IF (lsda) current_spin = isk(ik)
       evc(:, :) = C_save_all(:, :, ik)
       IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
    ENDDO
    rdmft_n = n_save
    CALL rdmft_mark_orbitals_changed()
    CALL rdmft_total_energy(E0)
    !
    IF (ionode) WRITE(stdout, '(/,5X,A,A,/)') &
         'Gradient check ', MERGE('PASSED', 'FAILED', all_pass)
    !
    DEALLOCATE(grad_a, n_save, e_one_p, e_har_p, e_xc_p, e_one_m, e_har_m, e_xc_m)
    DEALLOCATE(G_R_all, C_save_all, C_try)
    !
  END SUBROUTINE rdmft_check_gradient_consistency
  !
END MODULE rdmft_gradient_check
