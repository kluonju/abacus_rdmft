!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_krefine_mod
  !------------------------------------------------------------------
  !! Denser-k RDMFT restart: interpolate coarse saved \((n, C^k)\)
  !! onto the finer Monkhorst--Pack mesh from the current ``pw.x``
  !! ``K_POINTS`` input before the RDMFT optimiser runs.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_krefine_for_restart, rdmft_krefine_requested, &
            rdmft_resolve_krefine_from_source, rdmft_prepare_krefine_skip_scf
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  LOGICAL FUNCTION rdmft_krefine_requested()
    !---------------------------------------------------------------
    USE rdmft_module, ONLY : rdmft_do_krefine
    rdmft_krefine_requested = rdmft_do_krefine
  END FUNCTION rdmft_krefine_requested
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_prepare_krefine_skip_scf()
    !---------------------------------------------------------------
    !! Called from \texttt{electrons} before the KS loop.  When a
    !! denser-k RDMFT restart is pending, force \texttt{niter = 0} so
    !! PWscf does not diagonalize on the fine mesh; occupations and
    !! orbitals are bootstrapped from the coarse save inside
    !! \texttt{rdmft\_run}.
    USE io_global,     ONLY : stdout, ionode
    USE control_flags, ONLY : niter
    USE rdmft_module,  ONLY : do_rdmft, rdmft_do_krefine
    !
    IF (.NOT. do_rdmft) RETURN
    CALL rdmft_resolve_krefine_from_source()
    IF (.NOT. rdmft_do_krefine) RETURN
    niter = 0
    IF (ionode) WRITE(stdout, '(5X,A)') &
         'RDMFT k-refine: skipping KS SCF; bootstrap (n, C^k) from coarse save.'
  END SUBROUTINE rdmft_prepare_krefine_skip_scf
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_resolve_krefine_from_source()
    !---------------------------------------------------------------
    !! Decide whether to run denser-k interpolation before the RDMFT
    !! optimiser.  When \texttt{rdmft\_source\_prefix/outdir} point at a
    !! coarse save and the current automatic Monkhorst--Pack mesh is a
    !! denser Monkhorst--Pack mesh (e.g.\ 2$\times$2$\times$2
    !! $\to$ 3$\times$3$\times$3), set \texttt{rdmft\_do\_krefine}; when
    !! the meshes match, fall back to the normal KS-seeded bootstrap.
    USE io_global,     ONLY : stdout, ionode
    USE io_files,      ONLY : prefix, tmp_dir
    USE klist,         ONLY : nks
    USE start_k,       ONLY : nk1, nk2, nk3, k1, k2, k3
    USE rdmft_module,  ONLY : rdmft_source_prefix, rdmft_source_outdir, &
                               rdmft_do_krefine, rdmft_krefine_resolved
    USE rdmft_io,      ONLY : rdmft_read_source_kmesh
    !
    CHARACTER(LEN=256) :: src_prefix, src_outdir
    INTEGER :: nk1_c, nk2_c, nk3_c, k1_c, k2_c, k3_c
    INTEGER :: nk1_f, nk2_f, nk3_f, k1_f, k2_f, k3_f
    INTEGER :: nkstot_c, ierr
    REAL(DP), ALLOCATABLE :: xk_coarse_g(:,:)
    LOGICAL :: same_mp
    !
    IF (rdmft_krefine_resolved) RETURN
    rdmft_krefine_resolved = .TRUE.
    rdmft_do_krefine = .FALSE.
    !
    IF (LEN_TRIM(rdmft_source_prefix) == 0 .OR. &
        LEN_TRIM(rdmft_source_outdir) == 0) RETURN
    !
    src_prefix = TRIM(rdmft_source_prefix)
    src_outdir = TRIM(rdmft_source_outdir)
    IF (rdmft_krefine_same_save_path(src_outdir, src_prefix, tmp_dir, prefix)) &
         CALL errore('rdmft_resolve_krefine_from_source', &
              'fine run must use a different prefix/outdir than the coarse source save', 2)
    !
    nk1_f = nk1; nk2_f = nk2; nk3_f = nk3
    k1_f = k1; k2_f = k2; k3_f = k3
    IF (nk1_f <= 0 .OR. nk2_f <= 0 .OR. nk3_f <= 0) &
         CALL errore('rdmft_resolve_krefine_from_source', &
              'k-refine requires automatic Monkhorst-Pack K_POINTS on the fine run', 3)
    !
    CALL rdmft_read_source_kmesh(src_prefix, src_outdir, &
         nk1_c, nk2_c, nk3_c, k1_c, k2_c, k3_c, nkstot_c, xk_coarse_g, ierr)
    IF (ierr /= 0) &
         CALL errore('rdmft_resolve_krefine_from_source', &
              'could not read coarse k-mesh from source save', ierr)
    IF (ALLOCATED(xk_coarse_g)) DEALLOCATE(xk_coarse_g)
    !
    same_mp = (nk1_f == nk1_c .AND. nk2_f == nk2_c .AND. nk3_f == nk3_c .AND. &
               k1_f == k1_c .AND. k2_f == k2_c .AND. k3_f == k3_c)
    IF (same_mp) THEN
       IF (ionode) WRITE(stdout, '(5X,A)') &
            'RDMFT: source save k-mesh matches current K_POINTS; skipping k-refine.'
       RETURN
    ENDIF
    !
    IF (.NOT. rdmft_krefine_is_denser_mp(nk1_c, nk2_c, nk3_c, &
         nk1_f, nk2_f, nk3_f)) THEN
       IF (ionode) THEN
          WRITE(stdout, '(5X,A)') &
               'RDMFT k-refine: fine Monkhorst-Pack grid is not denser than the coarse save.'
          WRITE(stdout, '(5X,A,3(I0,1X),A,3(I0,1X))') &
               '  coarse MP grid  ', nk1_c, nk2_c, nk3_c, '  offsets ', k1_c, k2_c, k3_c
          WRITE(stdout, '(5X,A,3(I0,1X),A,3(I0,1X))') &
               '  fine   MP grid  ', nk1_f, nk2_f, nk3_f, '  offsets ', k1_f, k2_f, k3_f
       ENDIF
       CALL errore('rdmft_resolve_krefine_from_source', &
            'fine k mesh must be denser than the coarse mesh (each nk >= coarse, at least one >)', 4)
    ENDIF
    IF (k1_f /= k1_c .OR. k2_f /= k2_c .OR. k3_f /= k3_c) THEN
       IF (ionode) THEN
          WRITE(stdout, '(5X,A)') &
               'RDMFT k-refine: Monkhorst-Pack offsets differ between coarse and fine runs.'
          WRITE(stdout, '(5X,A,3(I0,1X),A,3(I0,1X))') &
               '  coarse MP grid  ', nk1_c, nk2_c, nk3_c, '  offsets ', k1_c, k2_c, k3_c
          WRITE(stdout, '(5X,A,3(I0,1X),A,3(I0,1X))') &
               '  fine   MP grid  ', nk1_f, nk2_f, nk3_f, '  offsets ', k1_f, k2_f, k3_f
       ENDIF
       CALL errore('rdmft_resolve_krefine_from_source', &
            'fine and coarse MP offsets must match', 5)
    ENDIF
    !
    rdmft_do_krefine = .TRUE.
    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') 'RDMFT: denser-k restart from coarse source save'
       WRITE(stdout, '(5X,A,A,A)') '  source  ', TRIM(src_outdir), TRIM(src_prefix)
       WRITE(stdout, '(5X,A,3(I0,1X),A,3(I0,1X))') &
            '  coarse MP grid  ', nk1_c, nk2_c, nk3_c, '  offsets ', k1_c, k2_c, k3_c
       WRITE(stdout, '(5X,A,3(I0,1X),A,3(I0,1X))') &
            '  fine   MP grid  ', nk1_f, nk2_f, nk3_f, '  offsets ', k1_f, k2_f, k3_f
       WRITE(stdout, '(5X,A,I0,A,I0)') &
            '  coarse k-list nks = ', nkstot_c, ',  fine k-list nks = ', nks
       IF (MOD(nks, nkstot_c) /= 0) WRITE(stdout, '(5X,A)') &
            '  (irreducible k counts need not divide evenly; interpolation uses the full MP grid)'
    ENDIF
    !
  END SUBROUTINE rdmft_resolve_krefine_from_source
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_krefine_for_restart()
    !---------------------------------------------------------------
    USE io_global,        ONLY : stdout, ionode
    USE io_files,         ONLY : prefix, tmp_dir, iunwfc, nwordwfc
    USE control_flags,    ONLY : gamma_only
    USE buffers,          ONLY : save_buffer
    USE cell_base,        ONLY : at
    USE wvfct,            ONLY : nbnd, npwx, current_k, g2kin, btype
    USE wavefunctions,    ONLY : evc
    USE klist,            ONLY : nks, xk, nelec, &
                                  nelup, neldw, two_fermi_energies, ngk, igk_k
    USE lsda_mod,         ONLY : nspin, lsda, isk, current_spin
    USE noncollin_module, ONLY : npol, noncolin, domag
    USE gvect,            ONLY : mill, gstart, g, ngm
    USE gvecw,            ONLY : gcutw
    USE start_k,          ONLY : nk1, nk2, nk3, k1, k2, k3
    USE rdmft_module,     ONLY : rdmft_n, rdmft_allocate, &
                                  rdmft_source_prefix, rdmft_source_outdir, &
                                  rdmft_krefine_done, rdmft_spin_factor, &
                                  rdmft_n_target, rdmft_fix_magnetization, &
                                  rdmft_n_target_up, rdmft_n_target_down
    USE rdmft_io,         ONLY : rdmft_load_state_raw, rdmft_save_filename_from, &
                                  rdmft_read_coarse_collected_wfc, &
                                  rdmft_read_source_kmesh
    USE rdmft_energy,     ONLY : rdmft_set_wg_from_n, rdmft_clear_band_diagnostics
    USE rdmft_stiefel,    ONLY : stiefel_orthonormalize_k, stiefel_orthonormalize_gamma
    !
    CHARACTER(LEN=256) :: src_prefix, src_outdir, src_fname, src_wfc_dir
    INTEGER :: nk1_c, nk2_c, nk3_c, k1_c, k2_c, k3_c
    INTEGER :: nk1_f, nk2_f, nk3_f, k1_f, k2_f, k3_f
    INTEGER :: nbnd_r, nkstot_c, nspin_r, nk_spatial, ik_s
    INTEGER :: nkr_c, ik, ib, ig, ispin, ierr, ic, ig_g, ig_cc, npw, mp_n
    INTEGER :: j1, j2, j3, j1p, j2p, j3p, ik_c(8), c1(8), c2(8), c3(8)
    INTEGER, ALLOCATABLE :: ngk_c(:), igk_c(:,:)
    REAL(DP), ALLOCATABLE :: n_coarse(:,:), wk_coarse(:)
    INTEGER, ALLOCATABLE :: isk_coarse(:)
    REAL(DP), ALLOCATABLE :: xkg_c(:,:), xk_cryst(:,:), xk_coarse_g(:,:), &
                             xk_spatial(:,:)
    INTEGER, ALLOCATABLE :: mp_to_ik_c(:,:)
    REAL(DP) :: w8(8), val8(8)
    REAL(DP) :: mill_f(3)
    LOGICAL :: loaded
    COMPLEX(DP), ALLOCATABLE :: wfc_stencil(:,:,:), wfc_coarse(:,:,:)
    REAL(DP), ALLOCATABLE :: gk_tmp(:)
    INTEGER :: npwx_c
    !
    IF (.NOT. rdmft_krefine_requested()) RETURN
    !
    src_prefix = TRIM(rdmft_source_prefix)
    src_outdir = TRIM(rdmft_source_outdir)
    IF (rdmft_krefine_same_save_path(src_outdir, src_prefix, tmp_dir, prefix)) &
         CALL errore('rdmft_krefine_for_restart', &
              'fine run must use a different prefix/outdir than the coarse source save', 2)
    !
    nk1_f = nk1; nk2_f = nk2; nk3_f = nk3
    k1_f = k1; k2_f = k2; k3_f = k3
    IF (nk1_f <= 0 .OR. nk2_f <= 0 .OR. nk3_f <= 0) &
         CALL errore('rdmft_krefine_for_restart', &
              'fine mesh must use automatic Monkhorst-Pack K_POINTS', 1)
    !
    CALL rdmft_read_source_kmesh(src_prefix, src_outdir, &
         nk1_c, nk2_c, nk3_c, k1_c, k2_c, k3_c, nk_spatial, xk_spatial, ierr)
    IF (ierr /= 0) &
         CALL errore('rdmft_krefine_for_restart', &
              'could not read coarse k-mesh from source save', ierr)
    !
    IF (.NOT. rdmft_krefine_is_denser_mp(nk1_c, nk2_c, nk3_c, &
         nk1_f, nk2_f, nk3_f)) &
         CALL errore('rdmft_krefine_for_restart', &
              'fine k mesh must be denser than the coarse mesh (each nk >= coarse, at least one >)', 2)
    IF (k1_f /= k1_c .OR. k2_f /= k2_c .OR. k3_f /= k3_c) &
         CALL errore('rdmft_krefine_for_restart', &
              'fine and coarse MP offsets must match', 3)
    !
    src_fname = rdmft_save_filename_from(src_prefix, src_outdir)
    CALL rdmft_load_state_raw(src_fname, nbnd, nspin, nbnd_r, nkstot_c, &
         nspin_r, n_coarse, wk_coarse, isk_coarse, loaded)
    IF (.NOT. loaded) &
         CALL errore('rdmft_krefine_for_restart', &
              'could not load coarse RDMFT occupations', 5)
    !
    ! LSDA RDMFT saves use nkstot = 2 * n_spatial (spin-up block then
    ! spin-down block) while the xml lists only the spatial k-points.
    IF (nspin_r == 2 .AND. nkstot_c == 2 * nk_spatial) THEN
       ALLOCATE(xk_coarse_g(3, nkstot_c))
       DO ik = 1, nkstot_c
          ik_s = MOD(ik - 1, nk_spatial) + 1
          xk_coarse_g(:, ik) = xk_spatial(:, ik_s)
       ENDDO
       DEALLOCATE(xk_spatial)
    ELSE IF (nkstot_c == nk_spatial) THEN
       xk_coarse_g = xk_spatial
       DEALLOCATE(xk_spatial)
    ELSE
       CALL errore('rdmft_krefine_for_restart', &
            'coarse RDMFT save k-count does not match source xml k-mesh', 6)
    ENDIF
    !
    src_wfc_dir = TRIM(src_outdir) // TRIM(src_prefix) // '.save/'
    !
    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)') 'RDMFT denser-k restart refinement'
       WRITE(stdout, '(5X,A,A,A)') '  coarse source save  ', TRIM(src_outdir), &
            TRIM(src_prefix)
       WRITE(stdout, '(5X,A,A,A)') '  fine output save    ', TRIM(tmp_dir), &
            TRIM(prefix)
       WRITE(stdout, '(5X,A,3(I0,1X),A,3(I0,1X))') &
            '  coarse MP grid  ', nk1_c, nk2_c, nk3_c, '  offsets ', k1_c, k2_c, k3_c
       WRITE(stdout, '(5X,A,3(I0,1X),A,3(I0,1X))') &
            '  fine   MP grid  ', nk1_f, nk2_f, nk3_f, '  offsets ', k1_f, k2_f, k3_f
       WRITE(stdout, '(5X,A,I0,A,I0,A,I0,A)') &
            '  coarse k-list nks = ', nkstot_c, ',  fine k-list nks = ', &
            nks, ',  fine full MP = ', nk1_f * nk2_f * nk3_f
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
    ENDIF
    !
    nkr_c = nk1_c * nk2_c * nk3_c
    ALLOCATE(xkg_c(3, nkr_c), xk_cryst(3, nkstot_c))
    CALL rdmft_krefine_build_mp_grid(nk1_c, nk2_c, nk3_c, k1_c, k2_c, k3_c, &
         xkg_c, nkr_c)
    xk_cryst = xk_coarse_g(:, 1:nkstot_c)
    CALL cryst_to_cart(nkstot_c, xk_cryst, at, -1)
    ALLOCATE(mp_to_ik_c(nkr_c, nspin))
    CALL rdmft_krefine_build_mp_to_ik(nkr_c, nk1_c, nk2_c, nk3_c, &
         k1_c, k2_c, k3_c, xkg_c, xk_cryst, isk_coarse, mp_to_ik_c)
    DO ik = 1, nkr_c
       DO ispin = 1, nspin
          IF (mp_to_ik_c(ik, ispin) == 0) &
               CALL errore('rdmft_krefine_for_restart', &
                    'coarse MP point could not be mapped to an irreducible k', ik)
       ENDDO
    ENDDO
    !
    ! ``npwx`` is the pool-local maximum over the fine k-list; coarse
    ! stencil corners can include k-points (e.g.\ $\Gamma$ on another
    ! pool) with a larger sphere.  Size the coarse G tables from the
    ! full coarse k-mesh instead of reusing pool-local ``npwx``.
    BLOCK
       INTEGER :: ngk_try, npwx_needed
       REAL(DP), ALLOCATABLE :: gk_work(:)
       INTEGER, ALLOCATABLE :: igk_work(:)
       npwx_needed = npwx
       ALLOCATE(gk_work(ngm), igk_work(ngm))
       DO ik = 1, nkstot_c
          CALL rdmft_krefine_gk_sort(xk_coarse_g(1, ik), ngm, g, gcutw, &
               ngk_try, igk_work, gk_work)
          npwx_needed = MAX(npwx_needed, ngk_try)
       ENDDO
       DEALLOCATE(gk_work, igk_work)
       npwx_c = npwx_needed
    END BLOCK
    ALLOCATE(ngk_c(nkstot_c), igk_c(npwx_c, nkstot_c), gk_tmp(npwx_c))
    DO ik = 1, nkstot_c
       CALL rdmft_krefine_gk_sort(xk_coarse_g(1, ik), ngm, g, gcutw, &
            ngk_c(ik), igk_c(:, ik), gk_tmp)
    ENDDO
    DEALLOCATE(gk_tmp)
    !
    ALLOCATE(wfc_coarse(npwx_c*npol, nbnd, nkstot_c))
    DO ik = 1, nkstot_c
       CALL rdmft_read_coarse_collected_wfc(src_wfc_dir, ik, isk_coarse(ik), &
            nkstot_c, ngk_c(ik), igk_c(:, ik), npwx_c, wfc_coarse(:,:,ik), ierr)
       IF (ierr /= 0) &
            CALL errore('rdmft_krefine_for_restart', &
                 'could not read coarse natural orbital from save dir', ik)
    ENDDO
    !
    CALL rdmft_allocate(nbnd, nks)
    rdmft_spin_factor = 2.0_DP
    IF (noncolin) rdmft_spin_factor = 1.0_DP
    rdmft_n_target = nelec
    rdmft_fix_magnetization = (lsda .AND. two_fermi_energies)
    IF (rdmft_fix_magnetization) THEN
       rdmft_n_target_up   = nelup
       rdmft_n_target_down = neldw
    ELSE
       rdmft_n_target_up   = 0.0_DP
       rdmft_n_target_down = 0.0_DP
    ENDIF
    !
    ALLOCATE(wfc_stencil(npwx_c*npol, nbnd, 8))
    !
    DO ik = 1, nks
       current_k = ik
       IF (lsda) current_spin = isk(ik)
       npw = ngk(ik)
       ispin = 1
       IF (lsda) ispin = isk(ik)
       !
       CALL rdmft_krefine_mp_stencil(xk(:, ik), nk1_c, nk2_c, nk3_c, &
            k1_c, k2_c, k3_c, j1, j2, j3, j1p, j2p, j3p, w8)
       CALL rdmft_krefine_stencil_ik(j1, j2, j3, j1p, j2p, j3p, nk1_c, nk2_c, &
            nk3_c, ispin, mp_to_ik_c, ik_c, ierr)
       IF (ierr /= 0) THEN
          c1 = (/ j1, j1p, j1, j1p, j1, j1p, j1, j1p /)
          c2 = (/ j2, j2, j2, j2, j2p, j2p, j2p, j2p /)
          c3 = (/ j3, j3, j3, j3, j3, j3, j3p, j3p /)
          mp_n = (c3(ierr) - 1) + (c2(ierr) - 1) * nk3_c &
               + (c1(ierr) - 1) * nk2_c * nk3_c + 1
          IF (ionode) WRITE(stdout, '(5X,A,I0,A,I0,A,3(I0,1X),A,I0,A,I0,A,I0)') &
               'k-refine stencil miss: fine ik=', ik, ' corner=', ierr, &
               ' mp=(', c1(ierr), c2(ierr), c3(ierr), ') mp_n=', mp_n, &
               ' ispin=', ispin, ' map=', mp_to_ik_c(mp_n, ispin)
          CALL errore('rdmft_krefine_for_restart', &
               'could not map fine k to coarse MP stencil', ik)
       ENDIF
       !
       DO ib = 1, nbnd
          DO ic = 1, 8
             val8(ic) = n_coarse(ib, ik_c(ic))
          ENDDO
          rdmft_n(ib, ik) = rdmft_krefine_interp8(w8, val8)
       ENDDO
       !
       DO ic = 1, 8
          wfc_stencil(:,:,ic) = wfc_coarse(:,:,ik_c(ic))
       ENDDO
       !
       evc = (0.0_DP, 0.0_DP)
       DO ig = 1, npw
          ig_g = igk_k(ig, ik)
          mill_f(1) = mill(1, ig_g)
          mill_f(2) = mill(2, ig_g)
          mill_f(3) = mill(3, ig_g)
          DO ib = 1, nbnd
             evc(ig, ib) = (0.0_DP, 0.0_DP)
             DO ic = 1, 8
                ig_cc = 0
                DO ig_g = 1, ngk_c(ik_c(ic))
                   IF (mill(1, igk_c(ig_g, ik_c(ic))) == mill_f(1) .AND. &
                       mill(2, igk_c(ig_g, ik_c(ic))) == mill_f(2) .AND. &
                       mill(3, igk_c(ig_g, ik_c(ic))) == mill_f(3)) THEN
                      ig_cc = ig_g
                      EXIT
                   ENDIF
                ENDDO
                IF (ig_cc > 0) &
                     evc(ig, ib) = evc(ig, ib) + w8(ic) * wfc_stencil(ig_cc, ib, ic)
             ENDDO
          ENDDO
       ENDDO
       !
       IF (gamma_only) THEN
          CALL stiefel_orthonormalize_gamma(evc, npw, nbnd, npwx*npol, gstart)
       ELSE
          CALL stiefel_orthonormalize_k(evc, npw, nbnd, npwx*npol)
       ENDIF
       CALL save_buffer(evc, nwordwfc, iunwfc, ik)
    ENDDO
    !
    DEALLOCATE(wfc_stencil, wfc_coarse, ngk_c, igk_c, xkg_c, xk_cryst, mp_to_ik_c)
    DEALLOCATE(xk_coarse_g)
    DEALLOCATE(n_coarse, wk_coarse, isk_coarse)
    !
    CALL rdmft_set_wg_from_n(1)
    CALL rdmft_clear_band_diagnostics()
    !
    rdmft_krefine_done = .TRUE.
    !
    IF (ionode) WRITE(stdout, '(/,5X,A,I0,A)') &
         'RDMFT k-refine: interpolated onto fine nks = ', nks, &
         ' (optional orbital polish runs after EXX init)'
    !
  END SUBROUTINE rdmft_krefine_for_restart
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_krefine_build_mp_grid(nk1i, nk2i, nk3i, k1i, k2i, k3i, &
       xkg, nkr)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN) :: nk1i, nk2i, nk3i, k1i, k2i, k3i, nkr
    REAL(DP), INTENT(OUT) :: xkg(3, nkr)
    INTEGER :: i, j, k, n
    !
    n = 0
    DO i = 1, nk1i
       DO j = 1, nk2i
          DO k = 1, nk3i
             n = n + 1
             xkg(1, n) = DBLE(i - 1) / nk1i + DBLE(k1i) / (2.0_DP * nk1i)
             xkg(2, n) = DBLE(j - 1) / nk2i + DBLE(k2i) / (2.0_DP * nk2i)
             xkg(3, n) = DBLE(k - 1) / nk3i + DBLE(k3i) / (2.0_DP * nk3i)
          ENDDO
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_krefine_build_mp_grid
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_krefine_mp_stencil(xk_cart, nk1i, nk2i, nk3i, k1i, k2i, k3i, &
       j1, j2, j3, j1p, j2p, j3p, w8)
    !---------------------------------------------------------------
    USE cell_base, ONLY : at
    REAL(DP), INTENT(IN)  :: xk_cart(3)
    INTEGER, INTENT(IN)  :: nk1i, nk2i, nk3i, k1i, k2i, k3i
    INTEGER, INTENT(OUT) :: j1, j2, j3, j1p, j2p, j3p
    REAL(DP), INTENT(OUT) :: w8(8)
    REAL(DP) :: xc(3), xx, yy, zz, fx, fy, fz
    INTEGER :: i1, i2, i3, ii
    !
    xc = xk_cart
    CALL cryst_to_cart(1, xc, at, -1)
    DO ii = 1, 3
       xc(ii) = xc(ii) - NINT(xc(ii))
    ENDDO
    !
    ! Same convention as ``kpoint_grid`` (xx = xk_cryst * nk - k/2).
    xx = xc(1) * DBLE(nk1i) - 0.5_DP * DBLE(k1i)
    yy = xc(2) * DBLE(nk2i) - 0.5_DP * DBLE(k2i)
    zz = xc(3) * DBLE(nk3i) - 0.5_DP * DBLE(k3i)
    !
    xx = rdmft_krefine_wrap_grid_coord(xx, nk1i)
    yy = rdmft_krefine_wrap_grid_coord(yy, nk2i)
    zz = rdmft_krefine_wrap_grid_coord(zz, nk3i)
    !
    i1 = INT(xx)
    IF (i1 >= nk1i) i1 = nk1i - 1
    IF (i1 < 0) i1 = 0
    i2 = INT(yy)
    IF (i2 >= nk2i) i2 = nk2i - 1
    IF (i2 < 0) i2 = 0
    i3 = INT(zz)
    IF (i3 >= nk3i) i3 = nk3i - 1
    IF (i3 < 0) i3 = 0
    !
    j1 = i1 + 1
    j2 = i2 + 1
    j3 = i3 + 1
    fx = xx - DBLE(i1)
    fy = yy - DBLE(i2)
    fz = zz - DBLE(i3)
    IF (fx >= 1.0_DP - 1.0e-12_DP) fx = 0.0_DP
    IF (fy >= 1.0_DP - 1.0e-12_DP) fy = 0.0_DP
    IF (fz >= 1.0_DP - 1.0e-12_DP) fz = 0.0_DP
    !
    j1p = j1 + 1; IF (j1p > nk1i) j1p = 1
    j2p = j2 + 1; IF (j2p > nk2i) j2p = 1
    j3p = j3 + 1; IF (j3p > nk3i) j3p = 1
    !
    w8(1) = (1.0_DP - fx) * (1.0_DP - fy) * (1.0_DP - fz)
    w8(2) = fx * (1.0_DP - fy) * (1.0_DP - fz)
    w8(3) = (1.0_DP - fx) * fy * (1.0_DP - fz)
    w8(4) = fx * fy * (1.0_DP - fz)
    w8(5) = (1.0_DP - fx) * (1.0_DP - fy) * fz
    w8(6) = fx * (1.0_DP - fy) * fz
    w8(7) = (1.0_DP - fx) * fy * fz
    w8(8) = fx * fy * fz
    !
  END SUBROUTINE rdmft_krefine_mp_stencil
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_krefine_wrap_grid_coord(t, nki)
    !---------------------------------------------------------------
    REAL(DP), INTENT(IN) :: t
    INTEGER, INTENT(IN)  :: nki
    REAL(DP) :: u
    !
    u = t - DBLE(nki) * INT(t / DBLE(nki))
    IF (u < 0.0_DP) u = u + DBLE(nki)
    IF (u >= DBLE(nki)) u = u - DBLE(nki)
    rdmft_krefine_wrap_grid_coord = u
  END FUNCTION rdmft_krefine_wrap_grid_coord
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_krefine_stencil_ik(j1, j2, j3, j1p, j2p, j3p, &
       nk1i, nk2i, nk3i, ispin, mp_to_ik, ik_c, ierr)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN) :: j1, j2, j3, j1p, j2p, j3p
    INTEGER, INTENT(IN) :: nk1i, nk2i, nk3i, ispin
    INTEGER, INTENT(IN) :: mp_to_ik(:,:)
    INTEGER, INTENT(OUT) :: ik_c(8)
    INTEGER, INTENT(OUT) :: ierr
    !
    INTEGER :: corners(3, 8), ic, mp_n
    !
    corners = reshape((/ &
         j1,  j2,  j3, &
         j1p, j2,  j3, &
         j1,  j2p, j3, &
         j1p, j2p, j3, &
         j1,  j2,  j3p, &
         j1p, j2,  j3p, &
         j1,  j2p, j3p, &
         j1p, j2p, j3p /), (/ 3, 8 /))
    !
    DO ic = 1, 8
       mp_n = (corners(3, ic) - 1) + (corners(2, ic) - 1) * nk3i &
            + (corners(1, ic) - 1) * nk2i * nk3i + 1
       ik_c(ic) = mp_to_ik(mp_n, ispin)
       IF (ik_c(ic) == 0) THEN
          ierr = ic
          RETURN
       ENDIF
    ENDDO
    ierr = 0
    !
  END SUBROUTINE rdmft_krefine_stencil_ik
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_krefine_build_mp_to_ik(nkr, nk1i, nk2i, nk3i, k1i, k2i, k3i, &
       xkg_c, xk_cryst, isk_c, mp_to_ik)
    !---------------------------------------------------------------
    !! Map each saved coarse irreducible k to a full-grid MP index.
    !! Primary route: invert the ``kpoint_grid`` indexing on folded
    !! crystal coordinates.  Fallback: symmetry search as in
    !! ``rdmft_brzint_setup_mp``.
    USE symm_base, ONLY : nsym, s, t_rev, time_reversal
    USE noncollin_module, ONLY : colin_mag
    USE lsda_mod, ONLY : nspin
    INTEGER, INTENT(IN) :: nkr, nk1i, nk2i, nk3i, k1i, k2i, k3i
    REAL(DP), INTENT(IN) :: xkg_c(3, nkr), xk_cryst(3, *)
    INTEGER, INTENT(IN) :: isk_c(:)
    INTEGER, INTENT(OUT) :: mp_to_ik(nkr, nspin)
    INTEGER :: mp_n, ispin, ik, ns, i, ii, j, k
    REAL(DP) :: xkr(3), deltap(3), deltam(3)
    LOGICAL :: on_grid
    REAL(DP), PARAMETER :: eps_loc = 1.0e-5_DP
    !
    mp_to_ik = 0
    !
    ! Direct grid-index map for each saved irreducible k.
    DO ik = 1, SIZE(isk_c)
       ispin = isk_c(ik)
       xkr = xk_cryst(:, ik)
       DO ii = 1, 3
          xkr(ii) = xkr(ii) - NINT(xkr(ii))
       ENDDO
       CALL rdmft_krefine_mp_index(xkr, nk1i, nk2i, nk3i, k1i, k2i, k3i, &
            on_grid, i, j, k)
       IF (.NOT. on_grid) CYCLE
       mp_n = (k - 1) + (j - 1) * nk3i + (i - 1) * nk2i * nk3i + 1
       IF (mp_n >= 1 .AND. mp_n <= nkr) mp_to_ik(mp_n, ispin) = ik
    ENDDO
    !
    ! Symmetry fallback for any still-unmapped (mp_n, ispin) pairs.
    DO mp_n = 1, nkr
       DO ispin = 1, nspin
          IF (mp_to_ik(mp_n, ispin) /= 0) CYCLE
          DO ik = 1, SIZE(isk_c)
             IF (isk_c(ik) /= ispin) CYCLE
             DO ns = 1, nsym
                DO i = 1, 3
                   xkr(i) = s(i, 1, ns) * xk_cryst(1, ik) + &
                            s(i, 2, ns) * xk_cryst(2, ik) + &
                            s(i, 3, ns) * xk_cryst(3, ik)
                ENDDO
                IF (t_rev(ns) == 1 .AND. colin_mag < 2) xkr = -xkr
                deltap(1) = xkr(1) - xkg_c(1, mp_n) - &
                     NINT(xkr(1) - xkg_c(1, mp_n))
                deltap(2) = xkr(2) - xkg_c(2, mp_n) - &
                     NINT(xkr(2) - xkg_c(2, mp_n))
                deltap(3) = xkr(3) - xkg_c(3, mp_n) - &
                     NINT(xkr(3) - xkg_c(3, mp_n))
                IF (DSQRT(deltap(1)**2 + deltap(2)**2 + deltap(3)**2) < eps_loc) THEN
                   mp_to_ik(mp_n, ispin) = ik
                   GOTO 10
                ENDIF
                IF (time_reversal) THEN
                   deltam(1) = xkr(1) + xkg_c(1, mp_n) - &
                        NINT(xkr(1) + xkg_c(1, mp_n))
                   deltam(2) = xkr(2) + xkg_c(2, mp_n) - &
                        NINT(xkr(2) + xkg_c(2, mp_n))
                   deltam(3) = xkr(3) + xkg_c(3, mp_n) - &
                        NINT(xkr(3) + xkg_c(3, mp_n))
                   IF (DSQRT(deltam(1)**2 + deltam(2)**2 + deltam(3)**2) < eps_loc) THEN
                      mp_to_ik(mp_n, ispin) = ik
                      GOTO 10
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
10       CONTINUE
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_krefine_build_mp_to_ik
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_krefine_mp_index(xkr, nk1i, nk2i, nk3i, k1i, k2i, k3i, &
       on_grid, i, j, k)
    !---------------------------------------------------------------
    REAL(DP), INTENT(IN)  :: xkr(3)
    INTEGER, INTENT(IN)   :: nk1i, nk2i, nk3i, k1i, k2i, k3i
    LOGICAL, INTENT(OUT)  :: on_grid
    INTEGER, INTENT(OUT)  :: i, j, k
    REAL(DP) :: xx, yy, zz
    REAL(DP), PARAMETER :: eps_loc = 1.0e-5_DP
    !
    xx = xkr(1) * DBLE(nk1i) - 0.5_DP * DBLE(k1i)
    yy = xkr(2) * DBLE(nk2i) - 0.5_DP * DBLE(k2i)
    zz = xkr(3) * DBLE(nk3i) - 0.5_DP * DBLE(k3i)
    on_grid = ABS(xx - NINT(xx)) <= eps_loc .AND. &
         ABS(yy - NINT(yy)) <= eps_loc .AND. &
         ABS(zz - NINT(zz)) <= eps_loc
    IF (.NOT. on_grid) RETURN
    i = MOD(NINT(xkr(1) * DBLE(nk1i) - 0.5_DP * DBLE(k1i) + 2 * DBLE(nk1i)), nk1i) + 1
    j = MOD(NINT(xkr(2) * DBLE(nk2i) - 0.5_DP * DBLE(k2i) + 2 * DBLE(nk2i)), nk2i) + 1
    k = MOD(NINT(xkr(3) * DBLE(nk3i) - 0.5_DP * DBLE(k3i) + 2 * DBLE(nk3i)), nk3i) + 1
    !
  END SUBROUTINE rdmft_krefine_mp_index
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_krefine_interp8(w8, val8)
    !---------------------------------------------------------------
    REAL(DP), INTENT(IN) :: w8(8), val8(8)
    INTEGER :: ic
    rdmft_krefine_interp8 = 0.0_DP
    DO ic = 1, 8
       rdmft_krefine_interp8 = rdmft_krefine_interp8 + w8(ic) * val8(ic)
    ENDDO
  END FUNCTION rdmft_krefine_interp8
  !
  !-----------------------------------------------------------------
  LOGICAL FUNCTION rdmft_krefine_same_save_path(outdir_a, prefix_a, &
       outdir_b, prefix_b)
    !---------------------------------------------------------------
    !! Return \texttt{.TRUE.} when two \texttt{(outdir, prefix)} pairs
    !! refer to the same RDMFT/PWscf save tree.
    CHARACTER(LEN=*), INTENT(IN) :: outdir_a, prefix_a, outdir_b, prefix_b
    CHARACTER(LEN=512) :: dir_a, dir_b
    !
    dir_a = rdmft_krefine_norm_outdir(outdir_a)
    dir_b = rdmft_krefine_norm_outdir(outdir_b)
    rdmft_krefine_same_save_path = (TRIM(dir_a) == TRIM(dir_b) .AND. &
         TRIM(prefix_a) == TRIM(prefix_b))
  END FUNCTION rdmft_krefine_same_save_path
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_krefine_gk_sort(k, ngm, g, ecut, ngk, igk, gk)
    !---------------------------------------------------------------
    !! Like QE ``gk_sort`` but bounds checks against ``SIZE(gk)``
    !! instead of the pool-local module ``npwx``.
    USE kinds,     ONLY : DP
    USE constants, ONLY : eps8
    REAL(DP), INTENT(IN)  :: k(3)
    INTEGER,  INTENT(IN)  :: ngm
    REAL(DP), INTENT(IN)  :: g(3, ngm)
    REAL(DP), INTENT(IN)  :: ecut
    INTEGER,  INTENT(OUT) :: ngk
    INTEGER,  INTENT(OUT) :: igk(:)
    REAL(DP), INTENT(OUT) :: gk(:)
    INTEGER :: ng, nk, npwx_loc
    REAL(DP) :: q, q2x
    !
    npwx_loc = SIZE(gk)
    IF (SIZE(igk) < npwx_loc) &
         CALL errore('rdmft_krefine_gk_sort', 'igk smaller than gk', 1)
    q2x = (SQRT(SUM(k(:)**2)) + SQRT(ecut))**2
    ngk = 0
    igk(:) = 0
    gk(:) = 0.0_DP
    DO ng = 1, ngm
       q = SUM((k(:) + g(:, ng))**2)
       IF (q <= eps8) q = 0.0_DP
       IF (q <= ecut) THEN
          ngk = ngk + 1
          IF (ngk > npwx_loc) &
               CALL errore('rdmft_krefine_gk_sort', 'gk workspace too small', 2)
          gk(ngk) = q
          igk(ngk) = ng
       ELSE
          IF (SUM(g(:, ng)**2) > (q2x + eps8)) EXIT
       ENDIF
    ENDDO
    IF (k(1)**2 + k(2)**2 + k(3)**2 > eps8) THEN
       CALL hpsort_eps(ngk, gk, igk, eps8)
       DO nk = 1, ngk
          gk(nk) = SUM((k(:) + g(:, igk(nk)))**2)
       ENDDO
    ENDIF
  END SUBROUTINE rdmft_krefine_gk_sort
  !
  !-----------------------------------------------------------------
  LOGICAL FUNCTION rdmft_krefine_is_denser_mp(nk1_c, nk2_c, nk3_c, &
       nk1_f, nk2_f, nk3_f)
    !---------------------------------------------------------------
    !! \texttt{.TRUE.} when the fine Monkhorst--Pack grid is strictly
    !! denser than the coarse one (each \texttt{nk\_f >= nk\_c}, at
    !! least one strict inequality).  Trilinear interpolation on the
    !! coarse full grid does not require integer commensurability
    !! (2$\times$2$\times$2 $\to$ 3$\times$3$\times$3 is allowed).
    INTEGER, INTENT(IN) :: nk1_c, nk2_c, nk3_c
    INTEGER, INTENT(IN) :: nk1_f, nk2_f, nk3_f
    !
    IF (nk1_f < nk1_c .OR. nk2_f < nk2_c .OR. nk3_f < nk3_c) THEN
       rdmft_krefine_is_denser_mp = .FALSE.
       RETURN
    ENDIF
    rdmft_krefine_is_denser_mp = (nk1_f > nk1_c .OR. &
         nk2_f > nk2_c .OR. nk3_f > nk3_c)
  END FUNCTION rdmft_krefine_is_denser_mp
  !
  !-----------------------------------------------------------------
  FUNCTION rdmft_krefine_norm_outdir(outdir_in) RESULT(outdir_out)
    !---------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(IN) :: outdir_in
    CHARACTER(LEN=512) :: outdir_out
    INTEGER :: n
    !
    outdir_out = TRIM(outdir_in)
    n = LEN_TRIM(outdir_out)
    IF (n > 0 .AND. outdir_out(n:n) /= '/') outdir_out = TRIM(outdir_out) // '/'
  END FUNCTION rdmft_krefine_norm_outdir
  !
END MODULE rdmft_krefine_mod
