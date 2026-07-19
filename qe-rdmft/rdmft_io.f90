!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_io
  !------------------------------------------------------------------
  !! Checkpoint / restart support for the RDMFT solver.
  !!
  !! Two complementary pieces of data are persisted across runs:
  !!
  !!   * the natural occupations ``rdmft_n(nbnd, nks)`` -- gathered
  !!     across k-point pools and written by the I/O node to a small
  !!     formatted file (default
  !!     ``{outdir}/{prefix}.rdmft.save``);
  !!
  !!   * the natural orbitals -- the QE wavefunction buffer
  !!     (``iunwfc``) already holds the current ``evc`` for every
  !!     k-point owned by this pool, and is automatically flushed
  !!     to disk by QE when the run exits.  To make the buffer
  !!     durable for any intermediate save we also call
  !!     ``punch('all')`` so the collected wavefunctions land in
  !!     the standard ``{outdir}/{prefix}.save/`` directory.  On
  !!     restart, the user requests
  !!     ``startingwfc = 'file'`` in the ``&electrons`` namelist
  !!     (or relies on the QE ``calculation = 'nscf'`` / 'bands'
  !!     restart logic) so that PW's standard ``wfcinit`` re-fills
  !!     the wavefunction buffer with the saved natural orbitals
  !!     before ``rdmft_run`` is entered.
  !!
  !! The save file is a tiny formatted text file with a header
  !! identifying the layout (``nbnd``, ``nkstot``, ``nspin``) so
  !! that a mismatch is detected and reported up front.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_save_state, rdmft_load_state, rdmft_load_state_raw, &
            rdmft_save_filename, rdmft_save_filename_from, &
            rdmft_read_coarse_collected_wfc, rdmft_read_source_kmesh
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  FUNCTION rdmft_save_filename() RESULT(fname)
    !---------------------------------------------------------------
    !! Build the absolute path to the RDMFT save file.  Honors the
    !! optional \texttt{rdmft\_restart\_file} keyword: if non-blank
    !! and contains a path separator the literal value is returned,
    !! otherwise the filename is placed under ``outdir``.
    USE io_files,     ONLY : prefix, tmp_dir
    USE rdmft_module, ONLY : rdmft_restart_file
    CHARACTER(LEN=512) :: fname
    !
    IF (LEN_TRIM(rdmft_restart_file) > 0) THEN
       IF (INDEX(rdmft_restart_file, '/') > 0) THEN
          fname = TRIM(rdmft_restart_file)
       ELSE
          fname = TRIM(tmp_dir) // TRIM(rdmft_restart_file)
       ENDIF
    ELSE
       fname = TRIM(tmp_dir) // TRIM(prefix) // '.rdmft.save'
    ENDIF
    !
  END FUNCTION rdmft_save_filename
  !
  !-----------------------------------------------------------------
  FUNCTION rdmft_save_filename_from(prefix_in, outdir_in) RESULT(fname)
    !---------------------------------------------------------------
    !! Build the path to ``{outdir}/{prefix}.rdmft.save`` (ignores
    !! \texttt{rdmft\_restart\_file} override; used when loading a
    !! coarse save during k-grid refinement in ``pw.x`` restart).
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in, outdir_in
    CHARACTER(LEN=512) :: fname
    !
    fname = TRIM(outdir_in) // TRIM(prefix_in) // '.rdmft.save'
    !
  END FUNCTION rdmft_save_filename_from
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_load_state_raw(fname, nbnd_expect, nspin_expect, &
       nbnd_r, nkstot_r, nspin_r, n_g, wk_g, isk_g, loaded)
    !---------------------------------------------------------------
    !! Read natural occupations from \texttt{fname} without requiring
    !! \texttt{nkstot\_r == nkstot} of the current run.  On success
    !! returns the global tables \texttt{n\_g(nbnd\_r,nkstot\_r)},
    !! \texttt{wk\_g(nkstot\_r)}, \texttt{isk\_g(nkstot\_r)} on every
    !! MPI rank.  Used by the k-mesh refinement restart path.
    USE io_global,    ONLY : ionode, ionode_id, stdout
    USE mp,           ONLY : mp_bcast
    USE mp_images,    ONLY : intra_image_comm
    !
    CHARACTER(LEN=*), INTENT(IN)  :: fname
    INTEGER, INTENT(IN)           :: nbnd_expect, nspin_expect
    INTEGER, INTENT(OUT)          :: nbnd_r, nkstot_r, nspin_r
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: n_g(:,:)
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: wk_g(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: isk_g(:)
    LOGICAL, INTENT(OUT)          :: loaded
    !
    INTEGER :: iu, ios, ib_r, ik_r, isk_r, outer_r, nread, ntot
    REAL(DP) :: wk_r, n_r, etot_r
    LOGICAL :: noncolin_r, ok_local
    CHARACTER(LEN=512) :: line
    !
    loaded = .FALSE.
    nbnd_r = 0
    nkstot_r = 0
    nspin_r = 0
    ok_local = .FALSE.
    !
    IF (ionode) THEN
       OPEN(NEWUNIT=iu, FILE=TRIM(fname), STATUS='OLD', &
            FORM='FORMATTED', ACTION='READ', IOSTAT=ios)
       IF (ios /= 0) THEN
          WRITE(stdout, '(/,5X,A,A)') &
               '** ERROR: RDMFT save file not found: ', TRIM(fname)
       ELSE
          READ(iu, '(A)', IOSTAT=ios) line
          IF (ios == 0) READ(iu, *, IOSTAT=ios) nbnd_r, nkstot_r, nspin_r
          IF (ios == 0) READ(iu, *, IOSTAT=ios) noncolin_r
          IF (ios == 0) READ(iu, *, IOSTAT=ios) outer_r
          IF (ios == 0) READ(iu, *, IOSTAT=ios) etot_r
          IF (ios == 0) READ(iu, '(A)', IOSTAT=ios) line
          IF (ios == 0 .AND. nbnd_r == nbnd_expect .AND. &
              nspin_r == nspin_expect .AND. nkstot_r > 0) THEN
             ntot = nbnd_r * nkstot_r
             ALLOCATE(n_g(nbnd_r, nkstot_r))
             n_g = 0.0_DP
             nread = 0
             DO
                READ(iu, *, IOSTAT=ios) ik_r, ib_r, isk_r, wk_r, n_r
                IF (ios /= 0) EXIT
                IF (ik_r >= 1 .AND. ik_r <= nkstot_r .AND. &
                    ib_r >= 1 .AND. ib_r <= nbnd_r) THEN
                   n_g(ib_r, ik_r) = n_r
                   nread = nread + 1
                ENDIF
             ENDDO
             IF (nread == ntot) THEN
                ALLOCATE(wk_g(nkstot_r), isk_g(nkstot_r))
                wk_g = 0.0_DP
                isk_g = 1
                REWIND(iu)
                READ(iu, '(A)', IOSTAT=ios) line
                READ(iu, *, IOSTAT=ios) nbnd_r, nkstot_r, nspin_r
                READ(iu, *, IOSTAT=ios) noncolin_r
                READ(iu, *, IOSTAT=ios) outer_r
                READ(iu, *, IOSTAT=ios) etot_r
                READ(iu, '(A)', IOSTAT=ios) line
                DO
                   READ(iu, *, IOSTAT=ios) ik_r, ib_r, isk_r, wk_r, n_r
                   IF (ios /= 0) EXIT
                   IF (ik_r >= 1 .AND. ik_r <= nkstot_r) THEN
                      wk_g(ik_r) = wk_r
                      isk_g(ik_r) = isk_r
                   ENDIF
                ENDDO
                ok_local = .TRUE.
                WRITE(stdout, '(/,5X,A,A)') &
                     'RDMFT k-refine: loaded coarse occupations from ', TRIM(fname)
                WRITE(stdout, '(5X,A,I0,A,I0,A,I0,A)') &
                     '  coarse (nbnd, nkstot, nspin) = (', nbnd_r, ',', &
                     nkstot_r, ',', nspin_r, ')'
             ELSE
                WRITE(stdout, '(/,5X,A)') &
                     '** ERROR: RDMFT save file truncated.'
                DEALLOCATE(n_g)
             ENDIF
          ELSE IF (ios == 0) THEN
             WRITE(stdout, '(/,5X,A)') &
                  '** ERROR: RDMFT save file header mismatch:'
             WRITE(stdout, '(5X,A,I0,A,I0,A,I0,A)') &
                  '   file has (nbnd, nkstot, nspin) = (', nbnd_r, ',', &
                  nkstot_r, ',', nspin_r, ')'
             WRITE(stdout, '(5X,A,I0,A,I0,A)') &
                  '   expect (nbnd, nspin) = (', nbnd_expect, ',', &
                  nspin_expect, ')'
          ENDIF
          CLOSE(iu)
       ENDIF
    ENDIF
    !
    CALL mp_bcast(ok_local, ionode_id, intra_image_comm)
    IF (.NOT. ok_local) THEN
       IF (ALLOCATED(n_g)) DEALLOCATE(n_g)
       IF (ALLOCATED(wk_g)) DEALLOCATE(wk_g)
       IF (ALLOCATED(isk_g)) DEALLOCATE(isk_g)
       RETURN
    ENDIF
    !
    CALL mp_bcast(nbnd_r, ionode_id, intra_image_comm)
    CALL mp_bcast(nkstot_r, ionode_id, intra_image_comm)
    CALL mp_bcast(nspin_r, ionode_id, intra_image_comm)
    IF (.NOT. ionode) THEN
       ALLOCATE(n_g(nbnd_r, nkstot_r))
       ALLOCATE(wk_g(nkstot_r), isk_g(nkstot_r))
    ENDIF
    CALL mp_bcast(n_g, ionode_id, intra_image_comm)
    CALL mp_bcast(wk_g, ionode_id, intra_image_comm)
    CALL mp_bcast(isk_g, ionode_id, intra_image_comm)
    loaded = .TRUE.
    !
  END SUBROUTINE rdmft_load_state_raw
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_save_state(outer_iter, etot_now)
    !---------------------------------------------------------------
    !! Save the natural occupations to disk and flush the
    !! wavefunction buffer (so the natural orbitals survive an
    !! abort or crash between this save and the next).
    !!
    !! Called from the alternating / joint outer loop every
    !! ``rdmft_save_every`` cycles, and once unconditionally at the
    !! end of \texttt{rdmft\_run}.
    USE io_global,    ONLY : ionode, stdout
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks, nkstot, wk, xk
    USE lsda_mod,     ONLY : nspin, isk
    USE noncollin_module, ONLY : noncolin
    USE control_flags, ONLY : io_level
    USE buffers,      ONLY : save_buffer
    USE wavefunctions, ONLY : evc
    USE io_files,     ONLY : iunwfc, nwordwfc, prefix
    USE rdmft_module, ONLY : rdmft_n, rdmft_power_alpha
    !
    INTEGER, INTENT(IN) :: outer_iter
    REAL(DP), INTENT(IN) :: etot_now
    !
    REAL(DP), ALLOCATABLE :: n_g(:,:), wk_g(:,:), wk_loc(:,:)
    INTEGER, ALLOCATABLE :: isk_g(:)
    REAL(DP), ALLOCATABLE :: isk_loc(:,:)
    REAL(DP), ALLOCATABLE :: isk_g_real(:,:)
    INTEGER :: ik, ib, iu, ios
    CHARACTER(LEN=512) :: fname
    !
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    !
    ! Collect occupations and per-k metadata across k-point pools so
    ! that ionode has the full global table.
    ALLOCATE(n_g(nbnd, nkstot))
    ALLOCATE(wk_loc(1, nks), wk_g(1, nkstot))
    ALLOCATE(isk_loc(1, nks), isk_g_real(1, nkstot), isk_g(nkstot))
    DO ik = 1, nks
       wk_loc(1, ik)  = wk(ik)
       isk_loc(1, ik) = REAL(isk(ik), DP)
    ENDDO
    CALL poolcollect(nbnd, nks, rdmft_n, nkstot, n_g)
    CALL poolcollect(1, nks, wk_loc, nkstot, wk_g)
    CALL poolcollect(1, nks, isk_loc, nkstot, isk_g_real)
    DO ik = 1, nkstot
       isk_g(ik) = NINT(isk_g_real(1, ik))
    ENDDO
    !
    ! Make the wavefunction buffer durable on disk.  When pw.x was
    ! launched with io_level <= 0 the buffer is held in memory and
    ! lost on a crash; calling save_buffer with the buffer's regular
    ! unit forces a flush.  For nks > 1 the orbital block of the
    ! solver already calls save_buffer at every accepted retraction,
    ! so this is a defensive no-op for those cases.  For nks == 1
    ! save_buffer keeps the buffer (and the on-disk file when present)
    ! in sync; combined with the QE iunwfc that QE keeps open
    ! throughout the run, the wavefunction state at the time of the
    ! save is preserved.
    IF (nks == 1 .AND. io_level >= 1) CALL save_buffer(evc, nwordwfc, iunwfc, 1)
    !
    fname = rdmft_save_filename()
    IF (ionode) THEN
       OPEN(NEWUNIT=iu, FILE=TRIM(fname), STATUS='REPLACE', &
            FORM='FORMATTED', ACTION='WRITE', IOSTAT=ios)
       IF (ios /= 0) THEN
          WRITE(stdout, '(/,5X,A,A)') &
               '** WARNING: RDMFT could not open restart file for writing: ', TRIM(fname)
       ELSE
          ! v1 layout, list-directed (free-format) so the reader is
          ! tolerant of width changes and of compiler-specific
          ! whitespace.  Header lines start with "#"; the version
          ! line is the first token of line 1.
          WRITE(iu, '(A)') '# RDMFT_RESTART v2'
          WRITE(iu, *)    nbnd, nkstot, nspin
          WRITE(iu, *)    noncolin
          WRITE(iu, *)    outer_iter
          WRITE(iu, *)    etot_now
          WRITE(iu, *)    rdmft_power_alpha
          WRITE(iu, '(A)') '# rows: ik ib isk wk n_ik'
          DO ik = 1, nkstot
             DO ib = 1, nbnd
                WRITE(iu, '(2I8, I4, ES20.12, ES24.16)') &
                     ik, ib, isk_g(ik), wk_g(1, ik), n_g(ib, ik)
             ENDDO
          ENDDO
          CLOSE(iu)
          WRITE(stdout, '(5X,A,A,A,I0,A)') &
               'RDMFT restart: wrote occupations to ', TRIM(fname), &
               ' (after outer iter ', outer_iter, ')'
       ENDIF
    ENDIF
    !
    ! Persist the natural orbitals in collected format.  punch('all')
    ! is collective so it must be called by every MPI rank; the
    ! routine handles ionode-only file I/O internally.  When the QE
    ! disk_io option is 'none' or 'minimal' the call falls through and
    ! only the in-memory wfc buffer is preserved -- which is still
    ! enough to resume within the same pw.x invocation but not across
    ! restarts; warn the user in that case.
    IF (io_level >= 0) THEN
       CALL punch('config')
    ELSE IF (ionode) THEN
       WRITE(stdout, '(5X,A)') &
            '** WARNING: disk_io low/none -- natural orbitals are NOT'
       WRITE(stdout, '(5X,A)') &
            '   persisted to disk.  Set disk_io = ''medium'' or ''high'' to'
       WRITE(stdout, '(5X,A)') &
            '   enable a true RDMFT restart between pw.x invocations.'
    ENDIF
    !
    DEALLOCATE(n_g, wk_loc, wk_g, isk_loc, isk_g_real, isk_g)
    !
  END SUBROUTINE rdmft_save_state
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_load_state(loaded)
    !---------------------------------------------------------------
    !! Read the natural occupations from the RDMFT restart file and
    !! scatter them across the k-point pools so each pool fills the
    !! ``rdmft_n(:, 1:nks)`` columns it owns.  Returns
    !! ``loaded = .FALSE.`` (and leaves ``rdmft_n`` untouched) when
    !! the file is missing, unreadable, or inconsistent with the
    !! current run's (nbnd, nkstot, nspin).  On success the natural
    !! orbitals are expected to already be in the wavefunction buffer
    !! (loaded by PWscf's standard ``startingwfc = 'file'`` /
    !! ``read_collected_wfc`` path).
    USE io_global,    ONLY : ionode, ionode_id, stdout
    USE mp,           ONLY : mp_bcast
    USE mp_images,    ONLY : intra_image_comm
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks, nkstot
    USE lsda_mod,     ONLY : nspin
    USE rdmft_module, ONLY : rdmft_n, rdmft_allocate, rdmft_power_alpha
    !
    LOGICAL, INTENT(OUT) :: loaded
    !
    REAL(DP), ALLOCATABLE :: n_g(:,:)
    INTEGER :: iu, ios, ib_r, ik_r, isk_r
    INTEGER :: nbnd_r, nkstot_r, nspin_r, outer_r, nread, ntot
    REAL(DP) :: wk_r, n_r, etot_r, alpha_r, alpha_before
    LOGICAL :: noncolin_r, ok_local, file_v2
    CHARACTER(LEN=512) :: fname, line
    !
    loaded = .FALSE.
    fname = rdmft_save_filename()
    !
    IF (.NOT. ALLOCATED(rdmft_n)) CALL rdmft_allocate(nbnd, nks)
    !
    ALLOCATE(n_g(nbnd, nkstot))
    n_g = 0.0_DP
    ok_local = .FALSE.
    !
    IF (ionode) THEN
       OPEN(NEWUNIT=iu, FILE=TRIM(fname), STATUS='OLD', &
            FORM='FORMATTED', ACTION='READ', IOSTAT=ios)
       IF (ios /= 0) THEN
          WRITE(stdout, '(/,5X,A,A)') &
               '** WARNING: RDMFT restart file not found: ', TRIM(fname)
          WRITE(stdout, '(5X,A)') &
               '   Falling back to KS-derived initial occupations.'
       ELSE
          ! v1 layout written by rdmft_save_state.  Use list-directed
          ! reads for the numeric header so we are tolerant of small
          ! format changes (whitespace, integer widths, etc.).  We
          ! validate (nbnd, nkstot, nspin) before accepting the file.
          READ(iu, '(A)', IOSTAT=ios) line
          file_v2 = (ios == 0 .AND. INDEX(line, 'v2') > 0)
          IF (ios == 0) READ(iu, *, IOSTAT=ios) nbnd_r, nkstot_r, nspin_r
          IF (ios == 0) READ(iu, *, IOSTAT=ios) noncolin_r
          IF (ios == 0) READ(iu, *, IOSTAT=ios) outer_r
          IF (ios == 0) READ(iu, *, IOSTAT=ios) etot_r
          IF (ios == 0 .AND. file_v2) READ(iu, *, IOSTAT=ios) alpha_r
          IF (ios == 0) READ(iu, '(A)', IOSTAT=ios) line ! column header comment
          !
          IF (ios == 0 .AND. nbnd_r == nbnd .AND. nkstot_r == nkstot .AND. &
              nspin_r == nspin) THEN
             ntot = nbnd * nkstot
             nread = 0
             DO
                READ(iu, *, IOSTAT=ios) ik_r, ib_r, isk_r, wk_r, n_r
                IF (ios /= 0) EXIT
                IF (ik_r >= 1 .AND. ik_r <= nkstot .AND. &
                    ib_r >= 1 .AND. ib_r <= nbnd) THEN
                   n_g(ib_r, ik_r) = n_r
                   nread = nread + 1
                ENDIF
             ENDDO
             IF (nread == ntot) THEN
                ok_local = .TRUE.
                IF (file_v2) THEN
                   alpha_before = rdmft_power_alpha
                   rdmft_power_alpha = alpha_r
                   IF (ABS(alpha_before - alpha_r) > 1.0e-8_DP) THEN
                      WRITE(stdout, '(5X,A,F8.5,A,F8.5,A)') &
                           'RDMFT restart: power_alpha restored ', alpha_r, &
                           ' (was ', alpha_before, ' from input)'
                   ENDIF
                ENDIF
                WRITE(stdout, '(/,5X,A,A)') &
                     'RDMFT restart: loaded occupations from ', TRIM(fname)
                WRITE(stdout, '(5X,A,I0,A,ES16.8,A)') &
                     '  previous outer iter = ', outer_r, &
                     ',  previous etot = ', etot_r, ' Ry'
                IF (.NOT. file_v2) THEN
                   WRITE(stdout, '(5X,A)') &
                        'NOTE: v1 save has no power_alpha; match dos.in to pw rdmft.in'
                ENDIF
             ELSE
                WRITE(stdout, '(/,5X,A)') &
                     '** WARNING: RDMFT restart file truncated; ignoring.'
                WRITE(stdout, '(5X,A,I0,A,I0,A)') &
                     '   read ', nread, ' rows, expected ', ntot, '.'
             ENDIF
          ELSE IF (ios == 0) THEN
             WRITE(stdout, '(/,5X,A)') &
                  '** WARNING: RDMFT restart file shape mismatch:'
             WRITE(stdout, '(5X,A,I0,A,I0,A,I0,A)') &
                  '   file has (nbnd, nkstot, nspin) = (', nbnd_r, ',', &
                  nkstot_r, ',', nspin_r, ')'
             WRITE(stdout, '(5X,A,I0,A,I0,A,I0,A)') &
                  '   run  has (nbnd, nkstot, nspin) = (', nbnd, ',', &
                  nkstot, ',', nspin, ')'
             WRITE(stdout, '(5X,A)') &
                  '   Falling back to KS-derived initial occupations.'
          ELSE
             WRITE(stdout, '(/,5X,A,A)') &
                  '** WARNING: malformed RDMFT restart file: ', TRIM(fname)
             WRITE(stdout, '(5X,A)') &
                  '   Falling back to KS-derived initial occupations.'
          ENDIF
          CLOSE(iu)
       ENDIF
    ENDIF
    !
    CALL mp_bcast(ok_local, ionode_id, intra_image_comm)
    IF (ok_local) THEN
       CALL mp_bcast(n_g, ionode_id, intra_image_comm)
       ! Scatter the (nbnd, nkstot) table to per-pool (nbnd, nks) slices.
       CALL poolscatter(nbnd, nkstot, n_g, nks, rdmft_n)
       loaded = .TRUE.
    ENDIF
    !
    DEALLOCATE(n_g)
    !
  END SUBROUTINE rdmft_load_state
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_read_coarse_collected_wfc(dirname, ik_g, ispin_k, &
       nkstot_g, ngk_k, igk_k, npwx_k, arr, ierr)
    !---------------------------------------------------------------
    !! Read one coarse collected natural orbital by global irreducible
    !! k index, using explicit G tables (for denser-k interpolation
    !! after the fine k mesh is already loaded).
    USE control_flags,    ONLY : gamma_only
    USE io_files,         ONLY : iunpun
    USE io_base,          ONLY : read_wfc
    USE lsda_mod,         ONLY : nspin
    USE gvect,            ONLY : ig_l2g
    USE mp,               ONLY : mp_sum, mp_max
    USE mp_bands,         ONLY : root_bgrp, intra_bgrp_comm
    USE pw_restart_new,   ONLY : gk_l2gmap_kdip
    USE wvfct,            ONLY : nbnd
    USE noncollin_module, ONLY : npol
    CHARACTER(LEN=*), INTENT(IN)  :: dirname
    INTEGER, INTENT(IN)           :: ik_g, ispin_k, nkstot_g, ngk_k, npwx_k
    INTEGER, INTENT(IN)           :: igk_k(:)
    COMPLEX(DP), INTENT(OUT)      :: arr(:,:)
    INTEGER, INTENT(OUT)          :: ierr
    CHARACTER(LEN=2), DIMENSION(2) :: updw = (/ 'up', 'dw' /)
    CHARACTER(LEN=320) :: filename
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    INTEGER :: ik_file, ispin, ig, npw_g, ngk_glob, nbnd_
    INTEGER, ALLOCATABLE :: igk_l2g(:), igk_l2g_kdip(:), mill_k(:,:)
    REAL(DP) :: xk_(3), b1(3), b2(3), b3(3), scalef
    INTEGER :: npol_
    !
    IF (nspin == 2) THEN
       ik_file = MOD(ik_g - 1, nkstot_g / 2) + 1
       ispin = ispin_k
       filename = TRIM(dirname) // 'wfc' // updw(ispin) // &
            TRIM(int_to_char(ik_file))
    ELSE
       ik_file = ik_g
       ispin = 1
       filename = TRIM(dirname) // 'wfc' // TRIM(int_to_char(ik_file))
    ENDIF
    !
    ALLOCATE(igk_l2g(npwx_k), igk_l2g_kdip(npwx_k), mill_k(3, npwx_k))
    igk_l2g = 0
    DO ig = 1, ngk_k
       igk_l2g(ig) = ig_l2g(igk_k(ig))
    ENDDO
    npw_g = MAXVAL(igk_l2g(1:ngk_k))
    CALL mp_max(npw_g, intra_bgrp_comm)
    ngk_glob = ngk_k
    CALL mp_sum(ngk_glob, intra_bgrp_comm)
    igk_l2g_kdip = 0
    CALL gk_l2gmap_kdip(npw_g, ngk_glob, ngk_k, igk_l2g, igk_l2g_kdip)
    !
    arr = (0.0_DP, 0.0_DP)
    CALL read_wfc(iunpun, filename, root_bgrp, intra_bgrp_comm, &
         ik_file, xk_, ispin, npol_, arr, npw_g, gamma_only, nbnd_, &
         igk_l2g_kdip(:), ngk_k, b1, b2, b3, mill_k, scalef, ierr)
    DEALLOCATE(igk_l2g, igk_l2g_kdip, mill_k)
    !
  END SUBROUTINE rdmft_read_coarse_collected_wfc
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_read_source_kmesh(src_prefix, src_outdir, &
       nk1_c, nk2_c, nk3_c, k1_c, k2_c, k3_c, nkstot_c, xk_coarse_g, ierr)
    !---------------------------------------------------------------
    !! Read the Monkhorst--Pack grid and **spatial** irreducible k
    !! coordinates from a coarse PWscf save ``data-file-schema.xml``.
    !! ``nkstot_c`` is the number of spatial k-points (``ndim``); for
    !! LSDA the RDMFT save may list ``2 * nkstot_c`` spin-resolved
    !! entries, expanded later in k-refinement.
    USE io_global,        ONLY : ionode, ionode_id
    USE mp,               ONLY : mp_bcast
    USE mp_images,        ONLY : intra_image_comm
    USE qexsd_module,     ONLY : qexsd_readschema
    USE qes_types_module, ONLY : output_type, parallel_info_type, &
                                  general_info_type, input_type
    USE qes_bcast_module, ONLY : qes_bcast
    USE qes_reset_module, ONLY : qes_reset
    USE cell_base,        ONLY : at
    CHARACTER(LEN=*), INTENT(IN)  :: src_prefix, src_outdir
    INTEGER, INTENT(OUT)          :: nk1_c, nk2_c, nk3_c
    INTEGER, INTENT(OUT)          :: k1_c, k2_c, k3_c, nkstot_c, ierr
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: xk_coarse_g(:,:)
    TYPE(output_type)        :: output_obj
    TYPE(parallel_info_type) :: parinfo_obj
    TYPE(general_info_type ) :: geninfo_obj
    TYPE(input_type)         :: input_obj
    CHARACTER(LEN=512) :: fname
    INTEGER :: ios, ik, ndim
    REAL(DP), ALLOCATABLE :: xk_cryst(:,:)
    !
    ierr = 0
    nk1_c = 0; nk2_c = 0; nk3_c = 0
    k1_c = 0; k2_c = 0; k3_c = 0
    nkstot_c = 0
    IF (ionode) THEN
       fname = TRIM(src_outdir) // TRIM(src_prefix) // &
            '.save/data-file-schema.xml'
       CALL qexsd_readschema(fname, ios, output_obj, parinfo_obj, &
            geninfo_obj, input_obj)
       IF (ios > 0) THEN
          ierr = ios
       ELSE
          IF (.NOT. output_obj%band_structure%starting_k_points% &
               monkhorst_pack_ispresent) THEN
             ierr = 1
          ELSE
             nk1_c = output_obj%band_structure%starting_k_points% &
                  monkhorst_pack%nk1
             nk2_c = output_obj%band_structure%starting_k_points% &
                  monkhorst_pack%nk2
             nk3_c = output_obj%band_structure%starting_k_points% &
                  monkhorst_pack%nk3
             k1_c = output_obj%band_structure%starting_k_points% &
                  monkhorst_pack%k1
             k2_c = output_obj%band_structure%starting_k_points% &
                  monkhorst_pack%k2
             k3_c = output_obj%band_structure%starting_k_points% &
                  monkhorst_pack%k3
             ndim = output_obj%band_structure%ndim_ks_energies
             nkstot_c = ndim
             IF (ndim <= 0) THEN
                ierr = 2
             ELSE
                ALLOCATE(xk_coarse_g(3, ndim), xk_cryst(3, ndim))
                DO ik = 1, ndim
                   xk_cryst(:, ik) = output_obj%band_structure% &
                        ks_energies(ik)%k_point%k_point(:)
                ENDDO
                CALL cryst_to_cart(ndim, xk_cryst, at, +1)
                xk_coarse_g = xk_cryst
                DEALLOCATE(xk_cryst)
             ENDIF
          ENDIF
          CALL qes_reset(output_obj)
       ENDIF
    ENDIF
    !
    CALL mp_bcast(ierr, ionode_id, intra_image_comm)
    IF (ierr /= 0) RETURN
    CALL mp_bcast(nk1_c, ionode_id, intra_image_comm)
    CALL mp_bcast(nk2_c, ionode_id, intra_image_comm)
    CALL mp_bcast(nk3_c, ionode_id, intra_image_comm)
    CALL mp_bcast(k1_c, ionode_id, intra_image_comm)
    CALL mp_bcast(k2_c, ionode_id, intra_image_comm)
    CALL mp_bcast(k3_c, ionode_id, intra_image_comm)
    CALL mp_bcast(nkstot_c, ionode_id, intra_image_comm)
    IF (.NOT. ALLOCATED(xk_coarse_g)) ALLOCATE(xk_coarse_g(3, nkstot_c))
    CALL mp_bcast(xk_coarse_g, ionode_id, intra_image_comm)
    !
  END SUBROUTINE rdmft_read_source_kmesh
  !
END MODULE rdmft_io
