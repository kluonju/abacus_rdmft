!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_pp_bootstrap
  !------------------------------------------------------------------
  !! Shared RDMFT post-processing bootstrap for ``pp.x`` and
  !! ``rdmft_dos.x``: load natural orbitals and occupations from a
  !! converged ``pw.x`` RDMFT run, initialise EXX/ACE, and optionally
  !! rebuild the RDMFT particle density.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_pp_is_rdmft_plot, rdmft_pp_bootstrap_state, &
            rdmft_pp_setup, rdmft_pp_teardown
  !
  LOGICAL, SAVE :: rdmft_pp_ready = .FALSE.
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  LOGICAL FUNCTION rdmft_pp_is_rdmft_plot(plot_num)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN) :: plot_num
    rdmft_pp_is_rdmft_plot = (plot_num == 126 .OR. plot_num == 127 &
         .OR. plot_num == 128)
  END FUNCTION rdmft_pp_is_rdmft_plot
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_pp_bootstrap_state(caller)
    !---------------------------------------------------------------
    USE io_global,        ONLY : stdout, ionode
    USE io_files,         ONLY : tmp_dir, nwordwfc, iunwfc, restart_dir, wfc_dir
    USE control_flags,    ONLY : io_level, gamma_only
    USE symme,            ONLY : sym_rho_init
    USE ions_base,        ONLY : nat, tau, ityp
    USE symm_base,        ONLY : set_sym, nsym
    USE buffers,          ONLY : open_buffer, save_buffer
    USE pw_restart_new,   ONLY : read_collected_wfc
    USE wavefunctions,    ONLY : evc
    USE wvfct,            ONLY : nbnd, npwx, btype
    USE klist,            ONLY : nks, nkstot
    USE noncollin_module, ONLY : npol, nspin_mag
    USE uspp,             ONLY : nkb
    USE becmod,           ONLY : becp, allocate_bec_type, is_allocated_bec_type
    USE gvect,            ONLY : ecutrho
    USE gvecw,            ONLY : ecutwfc
    USE xc_lib,           ONLY : xclib_set_exx_fraction, xclib_set_dft_IDs, &
                                  xclib_set_auxiliary_flags, exx_is_active, stop_exx, &
                                  xclib_dft_is
    USE exx,              ONLY : use_ace, ecutfock
    USE exx_base,         ONLY : exxalfa, exx_grid_initialized, exx_bgrp_type, EXX_BGRP_BANDS, &
                                  exxdiv_treatment, x_gamma_extrapolation
    USE rdmft_module,     ONLY : do_rdmft, rdmft_active, rdmft_functional, &
                                  rdmft_restart_file, rdmft_exxbuff_stale
    USE rdmft_xc,         ONLY : rdmft_xc_set_type
    USE rdmft_io,         ONLY : rdmft_load_state
    USE rdmft_energy,     ONLY : rdmft_initial_n_from_ks, rdmft_exxinit_once, &
                                  rdmft_clear_band_diagnostics, &
                                  rdmft_ensure_exxbuff_tsm_capacity
    CHARACTER(LEN=*), INTENT(IN) :: caller
    LOGICAL :: needwf, exst, loaded, dft_set_ok, repaired_tabxx
    INTEGER :: ik
    REAL(DP) :: ecutfock_old, ecutfock_pw
    REAL(DP), ALLOCATABLE :: m_loc_tmp(:,:)
    !
    do_rdmft = .TRUE.
    needwf = .TRUE.
    wfc_dir = tmp_dir
    nwordwfc = nbnd * npwx * npol
    IF (io_level /= 0) io_level = 1
    CALL open_buffer(iunwfc, 'wfc', nwordwfc, io_level, exst)
    IF (needwf) THEN
       WRITE(stdout, '(5X,A,A)') 'Reading collected natural orbitals from ', &
            TRIM(restart_dir())
       DO ik = 1, nks
          CALL read_collected_wfc(restart_dir(), ik, evc)
          CALL save_buffer(evc, nwordwfc, iunwfc, ik)
       ENDDO
    ELSE
       CALL errore(caller, 'collected natural orbitals not found in restart dir', 1)
    ENDIF
    !
    IF (nkb > 0 .AND. .NOT. is_allocated_bec_type(becp)) &
         CALL allocate_bec_type(nkb, nbnd, becp)
    CALL sym_rho_init(gamma_only)
    ! read_file_new does not call set_sym; brzint / EMD k-mesh unfolding
    ! needs the full crystal symmetry group (nsym > 1 for symmetry-
    ! reduced meshes such as bcc 2x2x2 with 3 IBZ k-points).
    ALLOCATE(m_loc_tmp(3, nat))
    m_loc_tmp = 0.0_DP
    CALL set_sym(nat, tau, ityp, nspin_mag, m_loc_tmp)
    DEALLOCATE(m_loc_tmp)
    IF (ionode) WRITE(stdout, '(5X,A,I0)') &
         'RDMFT post-proc: crystal symmetries nsym = ', nsym
    IF (.NOT. ALLOCATED(btype)) THEN
       ALLOCATE(btype(nbnd, nkstot))
       btype = 1
    ENDIF
    !
    ! Fock cutoff: match ``pw.x`` RDMFT on the same save.  Must run
    ! while the loaded SCF DFT is still active (before the HF switch
    ! below).  For a non-hybrid SCF starter (e.g. PBE), ``iosys`` sets
    ! ``ecutfock = MIN(ecutrho, 4*ecutwfc)``; the xml ``<hybrid>``
    ! block may still carry ``ecutfock`` in Hartree without the Ry
    ! conversion that ``read_xml_file`` only applies for hybrid DFT.
    ecutfock_pw = MIN(ecutrho, 4.0_DP * ecutwfc)
    ecutfock_old = ecutfock
    IF (.NOT. xclib_dft_is('hybrid')) THEN
       ecutfock = ecutfock_pw
    ELSE IF (ecutfock <= 0.0_DP) THEN
       ecutfock = ecutfock_pw
    ENDIF
    IF (ionode .AND. ABS(ecutfock - ecutfock_old) > 1.0e-8_DP) THEN
       WRITE(stdout, '(5X,A,F8.2,A,F8.2,A)') &
            'RDMFT post-proc: ecutfock adjusted ', ecutfock_old, &
            ' -> ', ecutfock, ' Ry (pw.x RDMFT match)'
    ENDIF
    CALL xclib_set_exx_fraction(1.0_DP)
    exxalfa = 1.0_DP
    dft_set_ok = xclib_set_dft_IDs(5, 0, 0, 0, 0, 0)
    CALL xclib_set_auxiliary_flags(.FALSE.)
    use_ace = .TRUE.
    exx_bgrp_type = EXX_BGRP_BANDS
    ! PBE-started RDMFT saves carry no hybrid block, so exxdiv_treatment
    ! is blank when read_file_new loads the xml.  setup_exx needs a valid
    ! singularity treatment (same default as ppacf.f90 for Fock plots).
    IF (LEN_TRIM(exxdiv_treatment) == 0) THEN
       exxdiv_treatment = 'gygi-baldereschi'
       x_gamma_extrapolation = .FALSE.
       IF (ionode) WRITE(stdout, '(5X,A)') &
            'RDMFT post-proc: exxdiv_treatment unset; using gygi-baldereschi'
    ENDIF
    IF (.NOT. exx_grid_initialized) CALL setup_exx()
    IF (exx_is_active()) CALL stop_exx()
    CALL rdmft_clear_band_diagnostics()
    rdmft_active = .TRUE.
    !
    CALL rdmft_initial_n_from_ks()
    CALL rdmft_load_state(loaded)
    IF (.NOT. loaded) &
         CALL errore(caller, 'could not load occupations from the RDMFT save file', 1)
    ! After ``rdmft_load_state`` (v2 files restore ``power_alpha``).
    CALL rdmft_xc_set_type(rdmft_functional)
    ! ``rdmft_exxinit_once`` sizes exxbuff from the widest channel
    ! support (``rdmft_set_wg_for_exxinit`` inside), not ``wk*n`` alone.
    CALL rdmft_exxinit_once()
    ! ``exx_fft_create`` aliases ``tabxx => tabp`` when
    ! ``ecutfock == ecutrho``, but ``dfftt`` is sized from
    ! ``max(ecutfock, gkcut)`` and is not always layout-compatible
    ! with ``dfftp`` even when the Fock and charge cutoffs match.
    ! The ``read_file`` post-proc path then underestimates US Fock
    ! (``E_xc`` too small).  Build a dedicated EXX-grid ``tabxx`` and
    ! refresh ``exxbuff`` when that alias was used.
    repaired_tabxx = rdmft_pp_repair_exx_tabxx()
    IF (repaired_tabxx) THEN
       rdmft_exxbuff_stale = .TRUE.
       CALL rdmft_exxinit_once()
    ENDIF
    CALL rdmft_ensure_exxbuff_tsm_capacity()
  END SUBROUTINE rdmft_pp_bootstrap_state
  !
  !-----------------------------------------------------------------
  LOGICAL FUNCTION rdmft_pp_repair_exx_tabxx()
    !---------------------------------------------------------------
    !! When ``exx_fft_create`` set ``tabxx => tabp``, replace it with
    !! augmentation tables built on the EXX FFT grid ``dfftt``.
    !! Returns .TRUE. if ``exxbuff`` must be rebuilt.
    USE io_global,     ONLY : stdout, ionode
    USE control_flags, ONLY : tqr
    USE uspp,          ONLY : okvan
    USE realus,        ONLY : qpointlist, tabxx, tabp
    USE exx,           ONLY : dfftt, exx_fft_initialized
    !
    rdmft_pp_repair_exx_tabxx = .FALSE.
    IF (.NOT. tqr .OR. .NOT. okvan .OR. .NOT. exx_fft_initialized) RETURN
    IF (.NOT. ASSOCIATED(tabxx, tabp)) RETURN
    IF (ionode) WRITE(stdout, '(5X,A)') &
         'RDMFT post-proc: rebuilding EXX-grid US augmentation (tabxx)'
    tabxx => NULL()
    CALL qpointlist(dfftt, tabxx)
    rdmft_pp_repair_exx_tabxx = .TRUE.
  END FUNCTION rdmft_pp_repair_exx_tabxx
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_pp_setup(caller, build_density, update_save_flag)
    !---------------------------------------------------------------
    USE io_global,     ONLY : stdout, ionode
    USE io_rho_xml,    ONLY : write_scf
    USE control_flags, ONLY : io_level
    USE scf,           ONLY : rho
    USE lsda_mod,      ONLY : nspin, lsda
    USE noncollin_module, ONLY : noncolin
    USE fft_base,      ONLY : dfftp
    USE klist,         ONLY : nelec
    USE rdmft_density, ONLY : rdmft_build_particle_density, rdmft_print_density_integrals
    CHARACTER(LEN=*), INTENT(IN) :: caller
    LOGICAL, INTENT(IN) :: build_density, update_save_flag
    !
    IF (rdmft_pp_ready) RETURN
    IF (noncolin .OR. (lsda .AND. nspin == 2)) &
         CALL errore(caller, 'RDMFT pp.x density plots not implemented for spin yet', 1)
    !
    CALL rdmft_pp_bootstrap_state(caller)
    IF (build_density) THEN
       IF (ionode) WRITE(stdout, '(/,5X,A)') 'RDMFT pp.x: rebuilding particle density'
       CALL rdmft_build_particle_density()
       IF (update_save_flag .AND. io_level > -2) CALL write_scf(rho, nspin)
       CALL rdmft_print_density_integrals()
    ENDIF
    rdmft_pp_ready = .TRUE.
  END SUBROUTINE rdmft_pp_setup
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_pp_teardown()
    !---------------------------------------------------------------
    USE rdmft_module, ONLY : rdmft_active
    IF (.NOT. rdmft_pp_ready) RETURN
    rdmft_active = .FALSE.
    rdmft_pp_ready = .FALSE.
  END SUBROUTINE rdmft_pp_teardown
  !
END MODULE rdmft_pp_bootstrap
