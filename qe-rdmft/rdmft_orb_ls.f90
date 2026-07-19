!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_orb_ls
  !------------------------------------------------------------------
  !! Line-search trial evaluation for alternating orbital blocks
  !! (product Stiefel joint step and per-k block_k step).
  !
  USE kinds, ONLY : DP
  USE mp,        ONLY : mp_sum
  USE mp_pools,  ONLY : inter_pool_comm
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_orb_ls_setup_joint, rdmft_orb_ls_setup_block_k, rdmft_orb_ls_eval
  !
  INTEGER, PARAMETER :: ORB_LS_JOINT = 1, ORB_LS_BLOCK_K = 2
  !
  INTEGER, SAVE :: orb_ls_mode = 0
  INTEGER, SAVE :: orb_nks = 0, orb_nbnd = 0, orb_npwx = 0, orb_npol = 0
  INTEGER, SAVE :: orb_ik = 0, orb_npw = 0, orb_nch = 0
  LOGICAL, SAVE :: orb_gamma = .FALSE.
  INTEGER, SAVE :: orb_gstart = 1
  REAL(DP), SAVE :: orb_wk_ik = 0.0_DP
  !
  COMPLEX(DP), POINTER, SAVE :: orb_C_save_all(:,:,:), orb_dir_all(:,:,:)
  COMPLEX(DP), POINTER, SAVE :: orb_C_save(:,:), orb_dir(:,:)
  COMPLEX(DP), ALLOCATABLE, SAVE :: orb_G_R(:,:,:), orb_C_try(:,:)
  COMPLEX(DP), ALLOCATABLE, SAVE :: orb_hpsi(:,:), orb_gradC(:,:), orb_eta(:,:)
  COMPLEX(DP), POINTER, SAVE :: orb_vxpsi(:,:,:)
  REAL(DP), POINTER, SAVE :: orb_vx_diag(:,:)
  !
CONTAINS
  !
  SUBROUTINE rdmft_orb_ls_setup_joint(C_save_all, dir_all, nks, npwx, nbnd, npol, &
                                      gamma_only, gstart, nch)
    COMPLEX(DP), TARGET, INTENT(IN) :: C_save_all(:,:,:), dir_all(:,:,:)
    INTEGER, INTENT(IN) :: nks, npwx, nbnd, npol, gstart, nch
    LOGICAL, INTENT(IN) :: gamma_only
    !
    orb_ls_mode = ORB_LS_JOINT
    orb_C_save_all => C_save_all
    orb_dir_all => dir_all
    orb_nks = nks
    orb_npwx = npwx
    orb_nbnd = nbnd
    orb_npol = npol
    orb_gamma = gamma_only
    orb_gstart = gstart
    orb_nch = nch
    IF (.NOT. ALLOCATED(orb_G_R)) &
         ALLOCATE(orb_G_R(npwx*npol, nbnd, nks), orb_C_try(npwx*npol, nbnd))
    IF (SIZE(orb_G_R, 3) /= nks) THEN
       DEALLOCATE(orb_G_R, orb_C_try)
       ALLOCATE(orb_G_R(npwx*npol, nbnd, nks), orb_C_try(npwx*npol, nbnd))
    ENDIF
  END SUBROUTINE rdmft_orb_ls_setup_joint
  !
  SUBROUTINE rdmft_orb_ls_setup_block_k(ik, C_save, dir, npw, npwx, nbnd, npol, &
                                      gamma_only, gstart, nch, wk_ik, vxpsi, vx_diag)
    INTEGER, INTENT(IN) :: ik, npw, npwx, nbnd, npol, gstart, nch
    REAL(DP), INTENT(IN) :: wk_ik
    LOGICAL, INTENT(IN) :: gamma_only
    COMPLEX(DP), TARGET, INTENT(IN) :: C_save(:,:), dir(:,:), vxpsi(:,:,:)
    REAL(DP), TARGET, INTENT(IN) :: vx_diag(:,:)
    !
    orb_ls_mode = ORB_LS_BLOCK_K
    orb_ik = ik
    orb_npw = npw
    orb_C_save => C_save
    orb_dir => dir
    orb_vxpsi => vxpsi
    orb_vx_diag => vx_diag
    orb_npwx = npwx
    orb_nbnd = nbnd
    orb_npol = npol
    orb_gamma = gamma_only
    orb_gstart = gstart
    orb_nch = nch
    orb_wk_ik = wk_ik
    IF (.NOT. ALLOCATED(orb_hpsi)) &
         ALLOCATE(orb_hpsi(npwx*npol, nbnd), orb_gradC(npwx*npol, nbnd), &
                 orb_eta(npwx*npol, nbnd), orb_C_try(npwx*npol, nbnd))
    IF (SIZE(orb_hpsi, 1) /= npwx*npol .OR. SIZE(orb_hpsi, 2) /= nbnd) THEN
       DEALLOCATE(orb_hpsi, orb_gradC, orb_eta, orb_C_try)
       ALLOCATE(orb_hpsi(npwx*npol, nbnd), orb_gradC(npwx*npol, nbnd), &
                orb_eta(npwx*npol, nbnd), orb_C_try(npwx*npol, nbnd))
    ENDIF
  END SUBROUTINE rdmft_orb_ls_setup_block_k
  !
  SUBROUTINE rdmft_orb_ls_eval(alpha, f, g, ierr)
    !! Wolfe callback: \(\phi'(\alpha)=\langle G_R(\alpha), d\rangle_{\mathrm S}\)
    !! for retraction along tangent direction \(d\) (textbook \(\nabla E^\top p\)).
    USE wvfct,             ONLY : current_k
    USE klist,             ONLY : nks, ngk, igk_k, xk, wk
    USE wavefunctions,     ONLY : evc
    USE io_files,          ONLY : nwordwfc, iunwfc
    USE buffers,           ONLY : get_buffer, save_buffer
    USE lsda_mod,          ONLY : lsda, current_spin, isk
    USE uspp,              ONLY : nkb, vkb
    USE uspp_init,         ONLY : init_us_2
    USE rdmft_module,      ONLY : rdmft_n, rdmft_mark_orbitals_changed
    USE rdmft_stiefel
    USE rdmft_energy,      ONLY : rdmft_total_energy, rdmft_compute_riemannian_gradient, &
                                  rdmft_apply_h_one_psi, rdmft_compute_xc_channel
    USE rdmft_xc,          ONLY : rdmft_xc_channel
    !
    REAL(DP), INTENT(IN)  :: alpha
    REAL(DP), INTENT(OUT) :: f, g
    INTEGER, INTENT(OUT)  :: ierr
    !
    REAL(DP) :: coef, w, dw, e_t, gnorm2
    INTEGER :: ik, ib, it, npw
    !
    ierr = 0
    IF (orb_ls_mode == ORB_LS_JOINT) THEN
       DO ik = 1, orb_nks
          npw = ngk(ik)
          current_k = ik
          IF (lsda) current_spin = isk(ik)
          IF (orb_gamma) THEN
             CALL stiefel_retract_gamma(orb_C_save_all(:, :, ik), orb_dir_all(:, :, ik), &
                  alpha, npw, orb_nbnd, orb_npwx*orb_npol, orb_gstart, orb_C_try)
          ELSE
             CALL stiefel_retract_k(orb_C_save_all(:, :, ik), orb_dir_all(:, :, ik), &
                  alpha, npw, orb_nbnd, orb_npwx*orb_npol, orb_C_try)
          ENDIF
          evc(:, :) = orb_C_try(:, :)
          IF (orb_nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
       ENDDO
       CALL rdmft_mark_orbitals_changed()
       CALL rdmft_total_energy(f)
       CALL rdmft_compute_riemannian_gradient(orb_G_R, gnorm2)
       g = 0.0_DP
       DO ik = 1, orb_nks
          npw = ngk(ik)
          IF (orb_gamma) THEN
             g = g + stiefel_inner_product_gamma(orb_G_R(:, :, ik), orb_dir_all(:, :, ik), &
                  npw, orb_nbnd, orb_npwx*orb_npol, orb_gstart)
          ELSE
             g = g + stiefel_inner_product_k(orb_G_R(:, :, ik), orb_dir_all(:, :, ik), &
                  npw, orb_nbnd, orb_npwx*orb_npol)
          ENDIF
       ENDDO
       ! The k-loop above only covers this pool's k-points; reduce the
       ! product-manifold directional derivative over the k-point pools
       ! so the Wolfe curvature test is consistent on every pool (the
       ! trial energy f is already global via rdmft_total_energy).
       CALL mp_sum(g, inter_pool_comm)
       ! Note: do NOT set ierr=1 when g >= 0.  Strong-Wolfe bracketing
       ! needs a valid (f, g) here so it can enter the zoom phase; the
       ! eval contract only signals ierr/=0 for genuine evaluation
       ! failures, not for a positive directional derivative.
       RETURN
    ENDIF
    !
    IF (orb_ls_mode == ORB_LS_BLOCK_K) THEN
       ik = orb_ik
       npw = orb_npw
       current_k = ik
       IF (lsda) current_spin = isk(ik)
       IF (orb_gamma) THEN
          CALL stiefel_retract_gamma(orb_C_save, orb_dir, alpha, npw, orb_nbnd, &
               orb_npwx*orb_npol, orb_gstart, orb_C_try)
       ELSE
          CALL stiefel_retract_k(orb_C_save, orb_dir, alpha, npw, orb_nbnd, &
               orb_npwx*orb_npol, orb_C_try)
       ENDIF
       evc(:, :) = orb_C_try(:, :)
       IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
       CALL rdmft_mark_orbitals_changed()
       CALL rdmft_total_energy(f)
       IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       CALL g2_kin(ik)
       IF (nkb > 0) CALL init_us_2(npw, igk_k(1, ik), xk(:, ik), vkb)
       CALL rdmft_apply_h_one_psi(npw, orb_nbnd, evc, orb_hpsi)
       orb_gradC = (0.0_DP, 0.0_DP)
       DO ib = 1, orb_nbnd
          orb_gradC(:, ib) = orb_wk_ik * rdmft_n(ib, ik) * orb_hpsi(:, ib)
       ENDDO
       DO it = 1, orb_nch
          ! Single-k aceinit_k build (only_k=ik): the BLOCK_K line
          ! search only changes ``evc(:, :, ik)``, so the all-k ACE
          ! projector ``xi`` is otherwise reusable.  Rebuilding ``xi``
          ! only at this k is ``O(N_k)`` instead of the ``O(N_k^2)``
          ! full-aceinit rebuild that would otherwise happen here on
          ! every Wolfe / Armijo backtrack.
          CALL rdmft_compute_xc_channel(it, e_t, orb_vx_diag, orb_vxpsi, only_k=ik)
          DO ib = 1, orb_nbnd
             CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
             orb_gradC(:, ib) = orb_gradC(:, ib) + coef * orb_wk_ik * w * orb_vxpsi(:, ib, ik)
          ENDDO
       ENDDO
       current_k = ik
       IF (lsda) current_spin = isk(ik)
       evc(:, :) = orb_C_try(:, :)
       IF (orb_gamma) THEN
          CALL stiefel_project_tangent_gamma(evc, orb_gradC, npw, orb_nbnd, &
               orb_npwx*orb_npol, orb_gstart, orb_eta)
          g = stiefel_inner_product_gamma(orb_eta, orb_dir, npw, orb_nbnd, &
               orb_npwx*orb_npol, orb_gstart)
       ELSE
          CALL stiefel_project_tangent_k(evc, orb_gradC, npw, orb_nbnd, &
               orb_npwx*orb_npol, orb_eta)
          g = stiefel_inner_product_k(orb_eta, orb_dir, npw, orb_nbnd, orb_npwx*orb_npol)
       ENDIF
       ! Do NOT set ierr=1 when g >= 0 — see comment in the JOINT branch.
       RETURN
    ENDIF
    !
    ierr = 1
    f = HUGE(1.0_DP)
    g = 0.0_DP
  END SUBROUTINE rdmft_orb_ls_eval
  !
END MODULE rdmft_orb_ls
