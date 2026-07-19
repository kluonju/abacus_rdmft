!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_dos
  !------------------------------------------------------------------
  !! RDMFT density of states and per-atom magnetization.
  !!
  !! Uses natural occupations \texttt{rdmft\_n} and the transition-state
  !! model (TSM) band energies -- the occupation gradient evaluated at
  !! half filling of each natural orbital,
  !! \((\partial E/\partial n_{ik})|_{n=1/2}\)
  !! (Slater/Janak transition state; the quantity ELK's \texttt{rdmeval}
  !! probes) -- for Gaussian-broadened total/partial DOS, plus
  !! Mulliken-style spin moments per atom site and per species.
  !!
    !! The spectral function follows Sharma \textit{et al.}, PRL
    !! \textbf{110}, 116403 (2013), Eq.~(7) (DER/TSM):
    !! \(\mathrm{DOS}=2\pi\sum_\lambda[n_\lambda\delta(\omega-\varepsilon^-_\lambda)
    !! +(1-n_\lambda)\delta(\omega+\varepsilon^+_\lambda)]\), with
    !! \(\varepsilon^\pm=(\partial E/\partial n)|_{n=1/2}\) (Ry) and
    !! \(\omega\) measured relative to the RDMFT chemical potential \(\mu\).
    !! Gaussian mode multiplies each branch by \(w_k\); \texttt{brzint} uses
    !! ELK \texttt{occmax}.  \texttt{rdmft\_dos\_spectral='elk'} (default)
    !! follows ELK \texttt{dos.f90} (single \(\delta\) per state);
    !! \texttt{'sharma'} uses PRL Eq.~(7) two-branch sum.
    !! \texttt{rdmft\_dos\_occ\_weighted=.true.} maps to ELK
    !! \texttt{dosocc=.true.} in \texttt{'elk'} mode.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_compute_dos_and_mag, rdmft_compute_dos, &
            rdmft_compute_magnetization, rdmft_dos_collect_and_print_band_tables, &
            rdmft_print_tsm_energies, rdmft_print_koopman_energies, &
            rdmft_print_band_diagonal_table, rdmft_report_chemical_potential
  !
  TYPE wfc_label
     INTEGER :: na, n, l, m, ind
     REAL(DP) :: jj
     CHARACTER(LEN=2) :: els
  END TYPE wfc_label
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_fill_nlmchi(natomwfc, lmax_wfc, nlmchi)
    !---------------------------------------------------------------
    !! Build atomic-wfc labels (collinear and noncollinear without SO).
    USE ions_base,          ONLY : ityp, nat
    USE upf_ions,           ONLY : n_atom_wfc
    USE upf_utils,          ONLY : l_to_spdf
    USE uspp_param,         ONLY : upf
    USE noncollin_module,   ONLY : noncolin, lspinorb
    !
    INTEGER, INTENT(OUT) :: natomwfc, lmax_wfc
    TYPE(wfc_label), ALLOCATABLE, INTENT(OUT) :: nlmchi(:)
    !
    INTEGER :: nwfc, na, nt, n, l, m
    CHARACTER(LEN=2) :: label
    INTEGER :: nn(0:3)
    !
    IF (noncolin .AND. lspinorb) &
         CALL errore('rdmft_fill_nlmchi', 'noncollinear PDOS with lspinorb not yet supported', 1)
    !
    natomwfc = n_atom_wfc(nat, ityp, noncolin)
    ALLOCATE(nlmchi(natomwfc))
    lmax_wfc = 0
    nwfc = 0
    DO na = 1, nat
       nt = ityp(na)
       nn = [1, 2, 3, 4]
       DO n = 1, upf(nt)%nwfc
          IF (upf(nt)%oc(n) >= 0.0_DP) THEN
             label = upf(nt)%els(n)
             l = upf(nt)%lchi(n)
             IF (label == 'Xn') THEN
                WRITE(label, '(I1,A1)') nn(l), l_to_spdf(l)
                nn(l) = nn(l) + 1
             END IF
             lmax_wfc = MAX(lmax_wfc, l)
             DO m = 1, 2 * l + 1
                nwfc = nwfc + 1
                nlmchi(nwfc)%na = na
                nlmchi(nwfc)%n = n
                nlmchi(nwfc)%l = l
                nlmchi(nwfc)%m = m
                nlmchi(nwfc)%ind = m
                nlmchi(nwfc)%jj = 0.0_DP
                nlmchi(nwfc)%els = label
             END DO
             IF (noncolin) THEN
                DO m = 1, 2 * l + 1
                   nwfc = nwfc + 1
                   nlmchi(nwfc)%na = na
                   nlmchi(nwfc)%n = n
                   nlmchi(nwfc)%l = l
                   nlmchi(nwfc)%m = m
                   nlmchi(nwfc)%ind = m + 2 * l + 1
                   nlmchi(nwfc)%jj = 0.0_DP
                   nlmchi(nwfc)%els = label
                END DO
             END IF
          END IF
       END DO
    END DO
    IF (lmax_wfc > 3) CALL errore('rdmft_fill_nlmchi', 'l > 3 not implemented', 1)
    IF (nwfc /= natomwfc) CALL errore('rdmft_fill_nlmchi', 'wrong # of atomic wfcs', 1)
    !
  END SUBROUTINE rdmft_fill_nlmchi
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ortho_atomic_wfc(npw, natomwfc, wfcatom, swfcatom, use_s, ierr, lda_in, n_inner)
    !---------------------------------------------------------------
    !! Orthonormalize atomic wfcs: O^{-1/2} |phi> with
    !! O_ij = <phi_i|S|phi_j> when \texttt{use_s}, else <phi_i|phi_j>.
    USE wvfct,    ONLY : npwx
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp,       ONLY : mp_sum
    !
    INTEGER, INTENT(IN)    :: npw, natomwfc
    COMPLEX(DP), INTENT(INOUT) :: wfcatom(:, :)
    COMPLEX(DP), INTENT(IN)    :: swfcatom(:, :)
    LOGICAL, INTENT(IN)    :: use_s
    INTEGER, INTENT(OUT)   :: ierr
    INTEGER, INTENT(IN), OPTIONAL :: lda_in, n_inner
    !
    COMPLEX(DP), ALLOCATABLE :: overlap(:,:), work(:,:), wfc_out(:,:)
    REAL(DP), ALLOCATABLE :: e(:)
    INTEGER :: i, j, k, lda, nin
    !
    lda = npwx
    nin = npw
    IF (PRESENT(lda_in)) lda = lda_in
    IF (PRESENT(n_inner)) nin = n_inner
    !
    ierr = 0
    ALLOCATE(overlap(natomwfc, natomwfc), work(natomwfc, natomwfc), e(natomwfc))
    ALLOCATE(wfc_out(lda, natomwfc))
    !
    overlap = (0.0_DP, 0.0_DP)
    IF (use_s) THEN
       CALL ZGEMM('C', 'N', natomwfc, natomwfc, nin, (1.0_DP, 0.0_DP), &
                  wfcatom, lda, swfcatom, lda, (0.0_DP, 0.0_DP), overlap, natomwfc)
    ELSE
       CALL ZGEMM('C', 'N', natomwfc, natomwfc, nin, (1.0_DP, 0.0_DP), &
                  wfcatom, lda, wfcatom, lda, (0.0_DP, 0.0_DP), overlap, natomwfc)
    END IF
    CALL mp_sum(overlap, intra_bgrp_comm)
    !
    CALL cdiagh(natomwfc, overlap, natomwfc, e, work)
    DO i = 1, natomwfc
       IF (ABS(e(i)) < 1.0e-10_DP) THEN
          ierr = 1
          DEALLOCATE(overlap, work, e, wfc_out)
          RETURN
       END IF
       e(i) = 1.0_DP / DSQRT(e(i))
    END DO
    overlap = (0.0_DP, 0.0_DP)
    DO i = 1, natomwfc
       DO j = 1, natomwfc
          DO k = 1, natomwfc
             overlap(i, j) = overlap(i, j) + e(k) * work(i, k) * DCONJG(work(j, k))
          END DO
       END DO
    END DO
    wfc_out = (0.0_DP, 0.0_DP)
    CALL ZGEMM('N', 'N', nin, natomwfc, natomwfc, (1.0_DP, 0.0_DP), &
               wfcatom, lda, overlap, natomwfc, (0.0_DP, 0.0_DP), wfc_out, lda)
    wfcatom(1:nin, 1:natomwfc) = wfc_out(1:nin, 1:natomwfc)
    DEALLOCATE(overlap, work, e, wfc_out)
    !
  END SUBROUTINE rdmft_ortho_atomic_wfc
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_atomic_proj(proj)
    !---------------------------------------------------------------
    !! |<phi_atom|psi_ib>|^2 for each k owned by the current pool.
    USE wvfct,            ONLY : nbnd, npwx, current_k
    USE wavefunctions,    ONLY : evc
    USE klist,            ONLY : nks, ngk, xk, igk_k
    USE control_flags,    ONLY : gamma_only
    USE noncollin_module, ONLY : noncolin, npol
    USE lsda_mod,         ONLY : nspin, isk, current_spin
    USE uspp,             ONLY : nkb, vkb
    USE uspp_init,        ONLY : init_us_2
    USE becmod,           ONLY : becp, allocate_bec_type, deallocate_bec_type, &
                                  is_allocated_bec_type, calbec
    USE io_files,         ONLY : nwordwfc, iunwfc
    USE buffers,          ONLY : get_buffer
    !
    REAL(DP), INTENT(OUT) :: proj(:, :, :)
    !
    TYPE(wfc_label), ALLOCATABLE :: nlmchi(:)
    COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:), swfcatom(:,:), proj_c(:,:)
    COMPLEX(DP), ALLOCATABLE :: wfcatom3(:,:,:)
    REAL(DP), ALLOCATABLE :: rproj(:,:)
    INTEGER :: ik, ib, npw, nwfc, natomwfc, lmax_wfc, ierr, npw_, lda
    !
    CALL rdmft_fill_nlmchi(natomwfc, lmax_wfc, nlmchi)
    IF (noncolin) THEN
       lda = npwx * npol
       ALLOCATE(wfcatom(lda, natomwfc), swfcatom(lda, natomwfc))
       ALLOCATE(wfcatom3(npwx, npol, natomwfc))
    ELSE
       lda = npwx
       ALLOCATE(wfcatom(npwx, natomwfc), swfcatom(npwx, natomwfc))
    ENDIF
    !
    DO ik = 1, nks
       npw = ngk(ik)
       current_k = ik
       IF (nspin == 2) current_spin = isk(ik)
       CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       wfcatom = (0.0_DP, 0.0_DP)
       IF (noncolin) THEN
          CALL atomic_wfc_nc_proj(ik, wfcatom3)
          wfcatom = RESHAPE(wfcatom3, (/ lda, natomwfc /))
          npw_ = npol * npwx
       ELSE
          CALL atomic_wfc(ik, wfcatom)
          npw_ = npwx
       END IF
       IF (nkb > 0) THEN
          IF (.NOT. is_allocated_bec_type(becp)) CALL allocate_bec_type(nkb, natomwfc, becp)
          CALL init_us_2(npw, igk_k(1, ik), xk(:, ik), vkb)
          CALL calbec(npw, vkb, wfcatom, becp)
          CALL s_psi(npwx, npw, natomwfc, wfcatom, swfcatom)
          IF (noncolin) THEN
             CALL rdmft_ortho_atomic_wfc(npw, natomwfc, wfcatom, swfcatom, .TRUE., ierr, &
                  lda_in=lda, n_inner=npw_)
          ELSE
             CALL rdmft_ortho_atomic_wfc(npw, natomwfc, wfcatom, swfcatom, .TRUE., ierr)
          END IF
       ELSE
          swfcatom = wfcatom
          IF (noncolin) THEN
             CALL rdmft_ortho_atomic_wfc(npw, natomwfc, wfcatom, swfcatom, .FALSE., ierr, &
                  lda_in=lda, n_inner=npw_)
          ELSE
             CALL rdmft_ortho_atomic_wfc(npw, natomwfc, wfcatom, swfcatom, .FALSE., ierr)
          END IF
       END IF
       IF (ierr /= 0) CALL errore('rdmft_compute_atomic_proj', &
            'atomic overlap matrix is singular', ik)
       IF (gamma_only .AND. .NOT. noncolin) THEN
          ALLOCATE(rproj(natomwfc, nbnd))
          CALL calbec(npw, wfcatom, evc, rproj)
          DO ib = 1, nbnd
             DO nwfc = 1, natomwfc
                proj(nwfc, ib, ik) = rproj(nwfc, ib) ** 2
             END DO
          END DO
          DEALLOCATE(rproj)
       ELSE
          ALLOCATE(proj_c(natomwfc, nbnd))
          CALL calbec(npw_, wfcatom, evc, proj_c)
          DO ib = 1, nbnd
             DO nwfc = 1, natomwfc
                proj(nwfc, ib, ik) = ABS(proj_c(nwfc, ib)) ** 2
             END DO
          END DO
          DEALLOCATE(proj_c)
       END IF
    END DO
    !
    IF (nkb > 0 .AND. is_allocated_bec_type(becp)) CALL deallocate_bec_type(becp)
    DEALLOCATE(wfcatom, swfcatom, nlmchi)
    IF (noncolin) DEALLOCATE(wfcatom3)
    !
  END SUBROUTINE rdmft_compute_atomic_proj
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_dos_and_mag()
    !---------------------------------------------------------------
    !! Legacy wrapper calling \texttt{rdmft\_compute\_dos} and
    !! \texttt{rdmft\_compute\_magnetization}.  Not invoked from
    !! \texttt{pw.x}; \texttt{rdmft\_dos.x} calls the two routines
    !! directly.
    USE rdmft_module, ONLY : rdmft_n, rdmft_compute_dos_flag => rdmft_compute_dos
    !
    IF (.NOT. rdmft_compute_dos_flag) RETURN
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    !
    CALL rdmft_compute_dos()
    CALL rdmft_compute_magnetization()
    !
  END SUBROUTINE rdmft_compute_dos_and_mag
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_dos_collect_and_print_band_tables()
    !---------------------------------------------------------------
    !! Collect per-band RDMFT diagnostics (one-body + XC diagonals and
    !! optional TSM energies) and print the TSM, Koopmans-Fock, and full
    !! band-diagonal tables.  Used by ``rdmft_dos.x`` only; ``pw.x``
    !! does not call this routine.
    USE rdmft_module, ONLY : rdmft_n
    USE rdmft_energy, ONLY : rdmft_collect_band_diagnostics
    !
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    CALL rdmft_collect_band_diagnostics(.TRUE.)
    CALL rdmft_print_tsm_energies()
    CALL rdmft_print_koopman_energies()
    CALL rdmft_print_band_diagonal_table()
    !
  END SUBROUTINE rdmft_dos_collect_and_print_band_tables
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_print_tsm_energies()
    !---------------------------------------------------------------
    !! Print the transition-state-model (TSM) band energies
    !! \(\mathrm{dedn}^{\mathrm{TSM}}_{ik}=(\partial E/\partial n_{ik})|
    !! _{n_{ik}=1/2}\) that drive the RDMFT DOS.  Every MPI rank
    !! must enter (the underlying probe is collective); only ionode
    !! writes to stdout.
    USE io_global,    ONLY : stdout, ionode
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks, nkstot, wk
    USE lsda_mod,     ONLY : lsda, isk
    USE rdmft_module, ONLY : rdmft_n
    USE rdmft_energy, ONLY : rdmft_diag_tsm, rdmft_diag_computed
    !
    REAL(DP) :: sum_n
    CHARACTER(LEN=9) :: slbl
    INTEGER  :: ik, ib, ispin
    REAL(DP), ALLOCATABLE :: n_g(:,:), tsm_g(:,:)
    REAL(DP), ALLOCATABLE :: wk_loc(:,:), wk_g(:,:), isk_loc(:,:), isk_g(:,:)
    !
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    IF (.NOT. rdmft_diag_computed) RETURN
    IF (.NOT. ALLOCATED(rdmft_diag_tsm)) &
         CALL errore('rdmft_print_tsm_energies', &
         'transition-state energies not available; call rdmft_collect_band_diagnostics(.TRUE.) first', 1)
    !
    ALLOCATE(n_g(nbnd, nkstot), tsm_g(nbnd, nkstot))
    ALLOCATE(wk_loc(1, nks), wk_g(1, nkstot), isk_loc(1, nks), isk_g(1, nkstot))
    DO ik = 1, nks
       wk_loc(1, ik)  = wk(ik)
       isk_loc(1, ik) = REAL(isk(ik), DP)
    ENDDO
    CALL poolcollect(nbnd, nks, rdmft_n,        nkstot, n_g)
    CALL poolcollect(nbnd, nks, rdmft_diag_tsm, nkstot, tsm_g)
    CALL poolcollect(1, nks, wk_loc,  nkstot, wk_g)
    CALL poolcollect(1, nks, isk_loc, nkstot, isk_g)
    !
    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') &
            'Transition-state-model (TSM) band energies (Ry):'
       WRITE(stdout, '(5X,A)') &
            '  ik    ib   spin      n_ik   dedn_tsm(n=1/2) (Ry)  wk*n_ik'
       DO ik = 1, nkstot
          ispin = NINT(isk_g(1, ik))
          IF (lsda) THEN
             slbl = MERGE('spin up  ', 'spin down', ispin == 1)
          ELSE
             slbl = 'spinless '
          ENDIF
          sum_n = 0.0_DP
          DO ib = 1, nbnd
             WRITE(stdout, '(5X,2I6,2X,A9,1X,F10.7,3X,F14.7,3X,F10.7)') &
                  ik, ib, slbl, n_g(ib, ik), tsm_g(ib, ik), &
                  wk_g(1, ik) * n_g(ib, ik)
             sum_n = sum_n + wk_g(1, ik) * n_g(ib, ik)
          ENDDO
          WRITE(stdout, '(5X,A,I0,A,F12.8)') &
               'sum_i wk*n_(i,k=', ik, ') = ', sum_n
       ENDDO
    ENDIF
    !
    DEALLOCATE(n_g, tsm_g, wk_loc, wk_g, isk_loc, isk_g)
    !
  END SUBROUTINE rdmft_print_tsm_energies
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_print_koopman_energies()
    !---------------------------------------------------------------
    !! Print Koopmans-Fock band energies at the converged natural
    !! occupations:
    !! \(\varepsilon^{\mathrm{F}}_{ik}=h_{ik}+v_{x,ik}\).
    !! Every MPI rank must enter; only ionode writes to stdout.
    USE io_global,    ONLY : stdout, ionode
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks, nkstot, wk
    USE lsda_mod,     ONLY : lsda, isk
    USE rdmft_module, ONLY : rdmft_n
    USE rdmft_energy, ONLY : rdmft_diag_h, rdmft_diag_vx, rdmft_diag_computed
    !
    REAL(DP) :: sum_n
    CHARACTER(LEN=9) :: slbl
    INTEGER  :: ik, ib, ispin
    REAL(DP), ALLOCATABLE :: n_g(:,:), h_g(:,:), vx_g(:,:), koop_g(:,:)
    REAL(DP), ALLOCATABLE :: wk_loc(:,:), wk_g(:,:), isk_loc(:,:), isk_g(:,:)
    !
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    IF (.NOT. rdmft_diag_computed) RETURN
    IF (.NOT. ALLOCATED(rdmft_diag_h) .OR. .NOT. ALLOCATED(rdmft_diag_vx)) &
         CALL errore('rdmft_print_koopman_energies', &
         'Koopmans-Fock diagonals not available; call rdmft_collect_band_diagnostics first', 1)
    !
    ALLOCATE(n_g(nbnd, nkstot), h_g(nbnd, nkstot), vx_g(nbnd, nkstot), &
             koop_g(nbnd, nkstot))
    ALLOCATE(wk_loc(1, nks), wk_g(1, nkstot), isk_loc(1, nks), isk_g(1, nkstot))
    DO ik = 1, nks
       wk_loc(1, ik)  = wk(ik)
       isk_loc(1, ik) = REAL(isk(ik), DP)
    ENDDO
    CALL poolcollect(nbnd, nks, rdmft_n,        nkstot, n_g)
    CALL poolcollect(nbnd, nks, rdmft_diag_h,   nkstot, h_g)
    CALL poolcollect(nbnd, nks, rdmft_diag_vx,  nkstot, vx_g)
    CALL poolcollect(1, nks, wk_loc,  nkstot, wk_g)
    CALL poolcollect(1, nks, isk_loc, nkstot, isk_g)
    koop_g = h_g + vx_g
    !
    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') &
            'Koopmans-Fock band energies (Ry):'
       WRITE(stdout, '(5X,A)') &
            '  ik    ib   spin      n_ik   eps^F=h+vx (Ry)  wk*n_ik'
       DO ik = 1, nkstot
          ispin = NINT(isk_g(1, ik))
          IF (lsda) THEN
             slbl = MERGE('spin up  ', 'spin down', ispin == 1)
          ELSE
             slbl = 'spinless '
          ENDIF
          sum_n = 0.0_DP
          DO ib = 1, nbnd
             WRITE(stdout, '(5X,2I6,2X,A9,1X,F10.7,3X,F14.7,3X,F10.7)') &
                  ik, ib, slbl, n_g(ib, ik), koop_g(ib, ik), &
                  wk_g(1, ik) * n_g(ib, ik)
             sum_n = sum_n + wk_g(1, ik) * n_g(ib, ik)
          ENDDO
          WRITE(stdout, '(5X,A,I0,A,F12.8)') &
               'sum_i wk*n_(i,k=', ik, ') = ', sum_n
       ENDDO
    ENDIF
    !
    DEALLOCATE(n_g, h_g, vx_g, koop_g, wk_loc, wk_g, isk_loc, isk_g)
    !
  END SUBROUTINE rdmft_print_koopman_energies
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_print_band_diagonal_table()
    !---------------------------------------------------------------
    !! Print h_diag, vx_diag, and Koopmans eps^F from cached
    !! diagnostics in :mod:`rdmft_energy`.  Does **not** compute them;
    !! call :subroutine:`rdmft_collect_band_diagnostics` first (via
    !! :subroutine:`rdmft_dos_collect_and_print_band_tables` in the
    !! ``rdmft_dos.x`` driver).  Every MPI rank must enter; only
    !! ionode writes.
    USE io_global,    ONLY : stdout, ionode
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks, nkstot, wk
    USE lsda_mod,     ONLY : lsda, isk
    USE rdmft_module, ONLY : rdmft_n
    USE rdmft_energy, ONLY : rdmft_diag_h, rdmft_diag_vx, rdmft_diag_computed
    !
    REAL(DP) :: sum_n
    CHARACTER(LEN=9) :: slbl
    INTEGER  :: ik, ib, ispin
    REAL(DP), ALLOCATABLE :: n_g(:,:), h_g(:,:), vx_g(:,:)
    REAL(DP), ALLOCATABLE :: wk_loc(:,:), wk_g(:,:), isk_loc(:,:), isk_g(:,:)
    !
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    IF (.NOT. rdmft_diag_computed) RETURN
    IF (.NOT. ALLOCATED(rdmft_diag_h) .OR. .NOT. ALLOCATED(rdmft_diag_vx)) RETURN
    !
    ALLOCATE(n_g(nbnd, nkstot), h_g(nbnd, nkstot), vx_g(nbnd, nkstot))
    ALLOCATE(wk_loc(1, nks), wk_g(1, nkstot), isk_loc(1, nks), isk_g(1, nkstot))
    DO ik = 1, nks
       wk_loc(1, ik)  = wk(ik)
       isk_loc(1, ik) = REAL(isk(ik), DP)
    ENDDO
    CALL poolcollect(nbnd, nks, rdmft_n,       nkstot, n_g)
    CALL poolcollect(nbnd, nks, rdmft_diag_h,  nkstot, h_g)
    CALL poolcollect(nbnd, nks, rdmft_diag_vx, nkstot, vx_g)
    CALL poolcollect(1, nks, wk_loc,  nkstot, wk_g)
    CALL poolcollect(1, nks, isk_loc, nkstot, isk_g)
    !
    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') &
            'Band diagonals at converged (n, C):'
       WRITE(stdout, '(5X,A)') &
            '  ik    ib   spin      n_ik       h_diag (Ry)   vx_diag (Ry)'// &
            '  eps^F=h+vx (Ry)  wk*n_ik'
       DO ik = 1, nkstot
          ispin = NINT(isk_g(1, ik))
          IF (lsda) THEN
             slbl = MERGE('spin up  ', 'spin down', ispin == 1)
          ELSE
             slbl = 'spinless '
          ENDIF
          sum_n = 0.0_DP
          DO ib = 1, nbnd
             WRITE(stdout, '(5X,2I6,2X,A9,1X,F10.7,3X,F12.7,3X,F12.7,3X,F14.7,3X,F10.7)') &
                  ik, ib, slbl, n_g(ib, ik), h_g(ib, ik), &
                  vx_g(ib, ik), h_g(ib, ik) + vx_g(ib, ik), &
                  wk_g(1, ik) * n_g(ib, ik)
             sum_n = sum_n + wk_g(1, ik) * n_g(ib, ik)
          ENDDO
          WRITE(stdout, '(5X,A,I0,A,F12.8)') &
               'sum_i wk*n_(i,k=', ik, ') = ', sum_n
       ENDDO
    ENDIF
    !
    DEALLOCATE(n_g, h_g, vx_g, wk_loc, wk_g, isk_loc, isk_g)
    !
  END SUBROUTINE rdmft_print_band_diagonal_table
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_dos_occmax()
    !---------------------------------------------------------------
    !! ELK \texttt{occmax}: 2 for closed-shell, 1 for lsda.
    USE lsda_mod, ONLY : nspin
    IF (nspin == 2) THEN
       rdmft_dos_occmax = 1.0_DP
    ELSE
       rdmft_dos_occmax = 2.0_DP
    ENDIF
  END FUNCTION rdmft_dos_occmax
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_r3cross(x, y, z)
    !---------------------------------------------------------------
    REAL(DP), INTENT(IN)  :: x(3), y(3)
    REAL(DP), INTENT(OUT) :: z(3)
    z(1) = x(2) * y(3) - x(3) * y(2)
    z(2) = x(3) * y(1) - x(1) * y(3)
    z(3) = x(1) * y(2) - x(2) * y(1)
  END SUBROUTINE rdmft_r3cross
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_axangsu2(v, th, su2)
    !---------------------------------------------------------------
    !! SU(2) spin rotation matrix (ELK \texttt{axangsu2}).
    REAL(DP), INTENT(IN)  :: v(3), th
    COMPLEX(DP), INTENT(OUT) :: su2(2, 2)
    REAL(DP) :: x, y, zz, cs, sn, t1
    x = v(1)
    y = v(2)
    zz = v(3)
    t1 = DSQRT(x*x + y*y + zz*zz)
    IF (t1 < 1.0e-8_DP) THEN
       su2 = RESHAPE((/ (1.0_DP, 0.0_DP), (0.0_DP, 0.0_DP), &
            (0.0_DP, 0.0_DP), (1.0_DP, 0.0_DP) /), (/ 2, 2 /))
       RETURN
    END IF
    t1 = 1.0_DP / t1
    x = x * t1
    y = y * t1
    zz = zz * t1
    cs = DCOS(0.5_DP * th)
    sn = DSIN(0.5_DP * th)
    su2(1, 1) = CMPLX(cs, -zz * sn, KIND=DP)
    su2(2, 1) = CMPLX(y * sn, -x * sn, KIND=DP)
    su2(1, 2) = CMPLX(-y * sn, -x * sn, KIND=DP)
    su2(2, 2) = CMPLX(cs, zz * sn, KIND=DP)
  END SUBROUTINE rdmft_axangsu2
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_sqasu2(sqaxis, su2, tsqaz)
    !---------------------------------------------------------------
    !! SU(2) operator rotating +z to \texttt{sqaxis} (ELK \texttt{sqasu2}).
    REAL(DP), INTENT(IN)  :: sqaxis(3)
    COMPLEX(DP), INTENT(OUT) :: su2(2, 2)
    LOGICAL, INTENT(OUT) :: tsqaz
    REAL(DP) :: v1(3), v2(3), v3(3), th, t1
    v1 = sqaxis
    t1 = DSQRT(v1(1)**2 + v1(2)**2 + v1(3)**2)
    IF (t1 <= 1.0e-8_DP) CALL errore('rdmft_sqasu2', &
         'spin-quantisation axis (sqaxis) has zero length', 1)
    v1 = v1 / t1
    IF (ABS(v1(3) - 1.0_DP) < 1.0e-8_DP) THEN
       tsqaz = .TRUE.
       su2 = RESHAPE((/ (1.0_DP, 0.0_DP), (0.0_DP, 0.0_DP), &
            (0.0_DP, 0.0_DP), (1.0_DP, 0.0_DP) /), (/ 2, 2 /))
    ELSE
       tsqaz = .FALSE.
       v2 = (/ 0.0_DP, 0.0_DP, 1.0_DP /)
       CALL rdmft_r3cross(v1, v2, v3)
       th = -DACOS(v1(3))
       CALL rdmft_axangsu2(v3, th, su2)
    END IF
  END SUBROUTINE rdmft_sqasu2
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_spin_weights(sc, nspinor)
    !---------------------------------------------------------------
    !! Spin-density-matrix diagonal weights (ELK \texttt{gensdmat} / sc).
    USE wvfct,            ONLY : nbnd, npwx, current_k
    USE wavefunctions,    ONLY : evc
    USE klist,            ONLY : nks, ngk
    USE noncollin_module, ONLY : npol
    USE io_files,         ONLY : nwordwfc, iunwfc
    USE buffers,          ONLY : get_buffer
    USE rdmft_module,     ONLY : rdmft_dos_sqaxis
    !
    REAL(DP), INTENT(OUT) :: sc(:, :, :)
    INTEGER, INTENT(IN)  :: nspinor
    !
    COMPLEX(DP) :: sd(2, 2), su2(2, 2), b(2, 2)
    COMPLEX(DP) :: cup, cdn
    REAL(DP) :: sqax(3)
    LOGICAL :: tsqaz
    INTEGER :: ik, ib, ig, ispn
    !
    sc = 0.0_DP
    sqax = rdmft_dos_sqaxis
    CALL rdmft_sqasu2(sqax, su2, tsqaz)
    !
    DO ik = 1, nks
       current_k = ik
       CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       DO ib = 1, nbnd
          sd = (0.0_DP, 0.0_DP)
          DO ig = 1, ngk(ik)
             cup = evc(ig, ib)
             cdn = evc(ig + npwx, ib)
             sd(1, 1) = sd(1, 1) + DCONJG(cup) * cup
             sd(2, 2) = sd(2, 2) + DCONJG(cdn) * cdn
             sd(1, 2) = sd(1, 2) + DCONJG(cup) * cdn
          ENDDO
          sd(2, 1) = DCONJG(sd(1, 2))
          IF (.NOT. tsqaz) THEN
             b = MATMUL(su2, sd)
             sd = MATMUL(b, CONJG(TRANSPOSE(su2)))
          END IF
          sc(ib, ik, 1) = REAL(sd(1, 1), DP)
          sc(ib, ik, 2) = REAL(sd(2, 2), DP)
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_compute_spin_weights
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_dos_state_weight(ib, ik, ispin, occ_g, wk_g, sc, &
                                           has_sc, use_brzint, use_sharma, &
                                           occupied)
    !---------------------------------------------------------------
    !! Spectral weight for ELK (\texttt{sc}$\times$\texttt{occmax}[/\texttt{n}])
    !! or Sharma (\(n\)/(1-\(n\))) branch sums.
    USE rdmft_module, ONLY : rdmft_dos_occ_weighted
    !
    INTEGER, INTENT(IN) :: ib, ik, ispin
    REAL(DP), INTENT(IN) :: occ_g(:, :), wk_g(:), sc(:, :, :)
    LOGICAL, INTENT(IN) :: has_sc, use_brzint, use_sharma, occupied
    !
    REAL(DP) :: occmax, sc_val, occ_fac, wt
    !
    occmax = rdmft_dos_occmax()
    sc_val = 1.0_DP
    IF (has_sc) sc_val = sc(ib, ik, ispin)
    IF (use_sharma) THEN
       IF (occupied) THEN
          occ_fac = occ_g(ib, ik)
       ELSE
          occ_fac = 1.0_DP - occ_g(ib, ik)
       END IF
       wt = sc_val * occmax * occ_fac
    ELSE
       IF (rdmft_dos_occ_weighted) THEN
          wt = sc_val * occmax * occ_g(ib, ik)
       ELSE
          wt = sc_val * occmax
       END IF
    END IF
    IF (use_brzint) THEN
       rdmft_dos_state_weight = wt
    ELSE
       rdmft_dos_state_weight = wk_g(ik) * wt
    END IF
  END FUNCTION rdmft_dos_state_weight
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_dos_efermi(eps, wk_g, isk_g, degauss, efermi)
    !---------------------------------------------------------------
    !! Fermi energy (Ry) from TSM band energies, ELK \texttt{occupy}
    !! bisection analogue using QE \texttt{efermig}.
    USE klist, ONLY : nelec, ngauss
    !
    REAL(DP), INTENT(IN)  :: eps(:, :), wk_g(:), degauss
    INTEGER, INTENT(IN)   :: isk_g(:)
    REAL(DP), INTENT(OUT) :: efermi
    !
    REAL(DP), ALLOCATABLE :: et(:,:)
    REAL(DP), EXTERNAL :: efermig
    INTEGER :: ik, ib
    !
    ALLOCATE(et(SIZE(eps, 1), SIZE(wk_g)))
    DO ik = 1, SIZE(wk_g)
       DO ib = 1, SIZE(eps, 1)
          et(ib, ik) = eps(ib, ik)
       ENDDO
    ENDDO
    efermi = efermig(et, SIZE(eps, 1), SIZE(wk_g), nelec, wk_g, degauss, &
         ngauss, 0, isk_g)
    DEALLOCATE(et)
  END SUBROUTINE rdmft_compute_dos_efermi
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_dedn_band_averages(grad_n, occ_g, wk_g, occ_band, dedn_band)
    !---------------------------------------------------------------
    !! \(w_k\)-averaged occupation \(n_i\) and per-band ELK-style
    !! ``dedn`` (\(\partial E/\partial n\), no extra \(w_k\) factor).
    USE rdmft_energy, ONLY : rdmft_dedn_from_grad_n
    REAL(DP), INTENT(IN)  :: grad_n(:,:), occ_g(:,:), wk_g(:)
    REAL(DP), INTENT(OUT) :: occ_band(:), dedn_band(:)
    !
    INTEGER :: ik, ib, nkstot
    REAL(DP) :: sum_w
    REAL(DP), ALLOCATABLE :: dedn_g(:,:)
    !
    nkstot = SIZE(wk_g)
    sum_w = 0.0_DP
    DO ik = 1, nkstot
       sum_w = sum_w + wk_g(ik)
    ENDDO
    IF (sum_w <= 0.0_DP) &
         CALL errore('rdmft_dedn_band_averages', 'zero total k-point weight', 1)
    ALLOCATE(dedn_g(SIZE(grad_n, 1), nkstot))
    CALL rdmft_dedn_from_grad_n(grad_n, wk_g, SIZE(grad_n, 1), nkstot, dedn_g)
    occ_band = 0.0_DP
    dedn_band = 0.0_DP
    DO ib = 1, SIZE(grad_n, 1)
       DO ik = 1, nkstot
          occ_band(ib) = occ_band(ib) + wk_g(ik) * occ_g(ib, ik)
          dedn_band(ib) = dedn_band(ib) + wk_g(ik) * dedn_g(ib, ik)
       ENDDO
       occ_band(ib) = occ_band(ib) / sum_w
       dedn_band(ib) = dedn_band(ib) / sum_w
    ENDDO
    DEALLOCATE(dedn_g)
    !
  END SUBROUTINE rdmft_dedn_band_averages
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_efermi_from_dedn_table(occ_band, dedn_band, ib_homo, &
       efermi, print_msg)
    !---------------------------------------------------------------
    !! RDMFT chemical potential from \texttt{RDM\_DEDN.OUT}: Janak
    !! energy \texttt{dedn} on the last band index with nonzero
    !! occupation \(n\) (lowest occupation among occupied bands).
    USE io_global, ONLY : stdout, ionode
    !
    REAL(DP), INTENT(IN)  :: occ_band(:), dedn_band(:)
    INTEGER, INTENT(OUT) :: ib_homo
    REAL(DP), INTENT(OUT) :: efermi
    LOGICAL, INTENT(IN), OPTIONAL :: print_msg
    !
    INTEGER :: ib
    REAL(DP), PARAMETER :: occ_thresh = 1.0e-6_DP
    LOGICAL :: do_print
    !
    do_print = .TRUE.
    IF (PRESENT(print_msg)) do_print = print_msg
    ib_homo = 0
    DO ib = 1, SIZE(occ_band)
       IF (occ_band(ib) > occ_thresh) ib_homo = ib
    ENDDO
    IF (ib_homo <= 0) &
         CALL errore('rdmft_efermi_from_dedn_table', &
              'no occupied band found for chemical potential', 1)
    efermi = dedn_band(ib_homo)
    !
    IF (do_print .AND. ionode) WRITE(stdout, '(5X,A,I0,A,F14.8,A)') &
         'mu = dedn(ib=', ib_homo, ') = ', efermi, &
         ' Ry  (RDM_DEDN.OUT)'
  END SUBROUTINE rdmft_efermi_from_dedn_table
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_write_efermi_file(fprefix, efermi)
    !---------------------------------------------------------------
    !! Write a single-line ELK \texttt{EFERMI.OUT} chemical potential (Ry).
    USE io_global, ONLY : ionode
    !
    CHARACTER(LEN=*), INTENT(IN) :: fprefix
    REAL(DP), INTENT(IN) :: efermi
    !
    INTEGER, EXTERNAL :: find_free_unit
    INTEGER :: iunit
    CHARACTER(LEN=256) :: fname
    !
    IF (.NOT. ionode) RETURN
    iunit = find_free_unit()
    fname = TRIM(fprefix) // '.efermi'
    OPEN(UNIT=iunit, FILE=TRIM(fname), STATUS='UNKNOWN', FORM='FORMATTED')
    WRITE(iunit, '(ES24.16)') efermi
    CLOSE(iunit)
  END SUBROUTINE rdmft_write_efermi_file
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_dos_mu(grad_n, occ_g, wk_g, mu)
    !---------------------------------------------------------------
    !! RDMFT chemical potential \(\mu\) (Ry) from the converged
    !! occupation gradient, following ELK \texttt{rdmvaryn} / Sharma
    !! \textit{et al.}, Phys.\ Rev.\ Lett.\ \textbf{110}, 116403
    !! (2013): find \(\mu\) such that \(\sum_{ik} w_k\,\gamma_{ik}=0\)
    !! with \(\gamma_{ik}=(\partial E/\partial n_{ik}-\mu)\,n_{ik}\) if
    !! \(\partial E/\partial n_{ik}<\mu\) and
    !! \((\partial E/\partial n_{ik}-\mu)\,(1-n_{ik})\) otherwise.
    REAL(DP), INTENT(IN)  :: grad_n(:,:), occ_g(:,:), wk_g(:)
    REAL(DP), INTENT(OUT) :: mu
    !
    REAL(DP) :: dkapa, gsp, gs, dgs, sm, t1, eps_mu
    INTEGER :: ik, ib, it
    !
    mu = 0.0_DP
    dkapa = 0.1_DP
    gsp = 0.0_DP
    eps_mu = 1.0e-8_DP
    DO it = 1, 200
       gs = 0.0_DP
       sm = 0.0_DP
       DO ik = 1, SIZE(wk_g)
          DO ib = 1, SIZE(grad_n, 1)
             t1 = grad_n(ib, ik) - mu
             IF (t1 > 0.0_DP) THEN
                gs = gs + wk_g(ik) * t1 * (1.0_DP - occ_g(ib, ik))
             ELSE
                gs = gs + wk_g(ik) * t1 * occ_g(ib, ik)
             ENDIF
             sm = sm + wk_g(ik) * t1**2
          ENDDO
       ENDDO
       sm = MAX(SQRT(sm), 1.0_DP)
       IF (ABS(gs) / sm < eps_mu) RETURN
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
       mu = mu + dkapa
    ENDDO
    CALL errore('rdmft_compute_dos_mu', 'could not determine RDMFT chemical potential', 1)
  END SUBROUTINE rdmft_compute_dos_mu
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_dos_state_omega(ib, ik, use_sharma, occupied, &
                                           eps_plus, eps_minus, wk_g, e_ref)
    !---------------------------------------------------------------
    !! Energy argument \(\omega\) (Ry) relative to \texttt{e\_ref}.
    INTEGER, INTENT(IN) :: ib, ik
    LOGICAL, INTENT(IN) :: use_sharma, occupied
    REAL(DP), INTENT(IN) :: eps_plus(:,:), eps_minus(:,:), wk_g(:), e_ref
    !
    REAL(DP) :: eps_abs
    !
    IF (wk_g(ik) < 1.0e-30_DP) THEN
       rdmft_dos_state_omega = 0.0_DP
       RETURN
    END IF
    IF (use_sharma) THEN
       IF (occupied) THEN
          eps_abs = eps_minus(ib, ik)
          rdmft_dos_state_omega = eps_abs - e_ref
       ELSE
          eps_abs = eps_plus(ib, ik)
          rdmft_dos_state_omega = -(eps_abs - e_ref)
       END IF
    ELSE
       eps_abs = eps_minus(ib, ik)
       rdmft_dos_state_omega = eps_abs - e_ref
    END IF
  END FUNCTION rdmft_dos_state_omega
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_dos_scan_omega_range(eps_arr, wk_g, nbnd, nkstot, use_sharma, &
       occ_weighted, e_ref, Elw, Eup)
    !---------------------------------------------------------------
    !! Scan \(\omega=\varepsilon-\mu\) extrema for a band-energy table.
    REAL(DP), INTENT(IN) :: eps_arr(:,:), wk_g(:), e_ref
    INTEGER, INTENT(IN) :: nbnd, nkstot
    LOGICAL, INTENT(IN) :: use_sharma, occ_weighted
    REAL(DP), INTENT(OUT) :: Elw, Eup
    !
    INTEGER :: ik, ib, ib_brz, n_branch
    LOGICAL :: occupied
    REAL(DP) :: omega_ik
    !
    n_branch = 1
    IF (use_sharma .AND. .NOT. occ_weighted) n_branch = 2
    Elw = 1.0e6_DP
    Eup = -1.0e6_DP
    DO ik = 1, nkstot
       DO ib = 1, nbnd
          DO ib_brz = 1, n_branch
             occupied = (ib_brz == 1)
             omega_ik = rdmft_dos_state_omega(ib, ik, use_sharma, occupied, &
                  eps_arr, eps_arr, wk_g, e_ref)
             Elw = MIN(Elw, omega_ik)
             Eup = MAX(Eup, omega_ik)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE rdmft_dos_scan_omega_range
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_write_band_energies_file(fprefix, eps_tsm, eps_koop, &
       occ_g, wk_g, isk_g)
    !---------------------------------------------------------------
    !! Write per-state TSM and Koopmans-Fock band energies
    !! (\texttt{<prefix>.rdmft.bands}); sole band-energy table for
    !! \texttt{rdmft\_dos.x}.
    USE io_global, ONLY : ionode
    USE wvfct,     ONLY : nbnd
    USE klist,     ONLY : nkstot
    USE lsda_mod,  ONLY : lsda
    !
    REAL(DP), INTENT(IN) :: eps_tsm(:, :), eps_koop(:, :), occ_g(:, :), wk_g(:)
    INTEGER, INTENT(IN) :: isk_g(:)
    CHARACTER(LEN=*), INTENT(IN) :: fprefix
    !
    INTEGER, EXTERNAL :: find_free_unit
    INTEGER :: iunit, ik, ib, ispin
    CHARACTER(LEN=9) :: slbl
    CHARACTER(LEN=256) :: fname
    REAL(DP) :: occ
    !
    IF (.NOT. ionode) RETURN
    iunit = find_free_unit()
    fname = TRIM(fprefix) // '.bands'
    OPEN(UNIT=iunit, FILE=TRIM(fname), STATUS='UNKNOWN', FORM='FORMATTED')
    WRITE(iunit, '(A)') '# ik  ib  spin  occ  eps_tsm(Ry)  eps_koopman(Ry)  wk'
    WRITE(iunit, '(A)') '# occ = natural occupation n in [0,1]; eps_tsm at n=1/2; eps_koopman = h+vx at converged n'
    DO ik = 1, nkstot
       ispin = isk_g(ik)
       IF (lsda) THEN
          IF (ispin == 1) THEN
             slbl = 'up  '
          ELSE
             slbl = 'down'
          END IF
       ELSE
          slbl = '--'
       END IF
       DO ib = 1, nbnd
          occ = occ_g(ib, ik)
          WRITE(iunit, '(2I6,2X,A4,1X,F10.7,2X,F14.7,2X,F14.7,2X,F10.7)') &
               ik, ib, slbl, occ, eps_tsm(ib, ik), &
               eps_koop(ib, ik), wk_g(ik)
       ENDDO
    ENDDO
    CLOSE(iunit)
  END SUBROUTINE rdmft_write_band_energies_file
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_dos()
    !---------------------------------------------------------------
    !! Build and write the RDMFT total and partial densities of states
    !! from Sharma \textit{et al.}, PRL \textbf{110}, 116403 (2013) Eq.~(7)
    !! with DER/TSM energies
    !! \(\varepsilon^\pm_{ik}=(\partial E/\partial n_{ik})|_{n_{ik}=1/2}\) (Ry),
    !! \(\omega=\varepsilon^--\mu\) or \(-(\varepsilon^+-\mu)\) per PRL Eq.~(7).
    !! Writes:
    !!
    !!   * \texttt{<prefix>.rdmft.dos}        -- total DOS from TSM and
    !!     Koopmans-Fock band energies (one file, both curves);
    !!   * \texttt{<prefix>.rdmft.bands}     -- TSM + Koopman energies
    !!     and natural occupation \(n\) per state;
    !!   * \texttt{<prefix>.rdmft.pdos\_at<S>\_<N>} -- l-resolved
    !!     partial DOS per atom site from TSM energies only (only
    !!     when the pseudopotentials carry atomic wave-functions and
    !!     \texttt{rdmft\_dos\_pdos} is enabled).
    USE io_global,        ONLY : stdout, ionode
    USE io_files,         ONLY : prefix
    USE wvfct,            ONLY : nbnd
    USE klist,            ONLY : nks, nkstot, wk, xk
    USE lsda_mod,         ONLY : nspin, isk, lsda
    USE noncollin_module, ONLY : noncolin
    USE ions_base,        ONLY : nat, ityp, atm
    USE start_k,          ONLY : nk1, nk2, nk3
    USE mp_pools,         ONLY : me_pool, root_pool, my_pool_id
    USE rdmft_module,     ONLY : rdmft_n, rdmft_compute_dos_flag => rdmft_compute_dos, &
                                  rdmft_dos_spectral, rdmft_dos_eref, &
                                  rdmft_dos_emin, rdmft_dos_emax, rdmft_dos_deltae, &
                                  rdmft_dos_degauss, rdmft_dos_pdos, rdmft_dos_prefix, &
                                  rdmft_dos_occ_weighted, rdmft_dos_integration, &
                                  rdmft_dos_nwplot, rdmft_dos_nswplot, rdmft_dos_ngrkf, &
                                  rdmft_dos_msum, rdmft_dos_ssum, &
                                  rdmft_exxbuff_stale, rdmft_xi_cache_valid
    USE rdmft_energy,     ONLY : rdmft_collect_band_diagnostics, rdmft_grad_n, &
                                  rdmft_grad_n_from_cached_diag, &
                                  rdmft_diag_tsm, rdmft_diag_computed, &
                                  rdmft_diag_h, rdmft_diag_vx, &
                                  rdmft_write_dedn
    USE rdmft_xc,         ONLY : rdmft_xc_n_channels
    USE rdmft_brzint_mod, ONLY : rdmft_brzint, rdmft_brzint_setup_mp
    !
    REAL(DP), PARAMETER :: E_UNSET = 1.0e6_DP
    REAL(DP), EXTERNAL :: w0gauss
    INTEGER, EXTERNAL :: find_free_unit
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    REAL(DP), ALLOCATABLE :: eps(:,:), eps_plus(:,:), eps_minus(:,:), grad_mu(:,:)
    REAL(DP), ALLOCATABLE :: occ_g(:,:), koop_g(:,:), proj(:,:,:), sc(:,:,:)
    REAL(DP), ALLOCATABLE :: tdos(:,:), tdos_tsm(:,:), tdos_koop(:,:), pdos(:,:,:,:), pdos_lm(:,:,:,:)
    REAL(DP), ALLOCATABLE :: pdos_sum(:,:), idos(:,:), wk_g(:)
    REAL(DP), ALLOCATABLE :: e_brz(:,:), f_brz(:,:), g_brz(:), wint(:)
    REAL(DP), ALLOCATABLE :: xk_g(:,:)
    REAL(DP), ALLOCATABLE :: grad_loc(:,:), eps_loc(:,:), occ_loc(:,:)
    REAL(DP), ALLOCATABLE :: occ_band(:), dedn_band(:)
    REAL(DP), ALLOCATABLE :: wk_loc(:,:), isk_loc(:,:), wk_col(:,:), isk_col(:,:)
    REAL(DP), ALLOCATABLE :: proj_pack(:,:), proj_col(:,:)
    INTEGER,  ALLOCATABLE :: isk_g(:), equiv(:), ivkik(:,:,:)
    TYPE(wfc_label), ALLOCATABLE :: nlmchi(:)
    !
    REAL(DP) :: Emin, Emax, DeltaE, degauss, Elw, Eup, E, delta, weight, occmax, omega_ik
    REAL(DP) :: e_ref_model, e_ref_tsm, e_ref_koop, mu_rdm
    REAL(DP) :: Elw_tsm, Eup_tsm, Elw_koop, Eup_koop
    REAL(DP) :: pdos_w
    CHARACTER(LEN=256) :: fprefix, fname
    CHARACTER(LEN=64) :: dos_label
    INTEGER :: natomwfc, lmax_wfc, lmmax, ne, ie, ie_mid, ik, ib, nwfc, na, nt, l, lm
    INTEGER :: is, ispn, nsd, iunit, ie_delta, nkr, ierr, ngridk(3), nsk(3)
    INTEGER :: nw_brz, ikr, ie_b, ld_brz, n_brz, ib_brz, n_branch, im
    LOGICAL :: do_pdos, use_brzint, has_sc, do_msum, occupied, use_sharma, &
               use_efermig, do_pdos_run
    INTEGER :: ib_homo
    !
    IF (.NOT. rdmft_compute_dos_flag) RETURN
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    !
    ! Janak \texttt{dedn} table for \texttt{RDM\_DEDN.OUT} and chemical potential.
    ALLOCATE(grad_mu(nbnd, nkstot), grad_loc(nbnd, nks))
    ! ``rdmft_collect_band_diagnostics`` (called by the driver before
    ! the print routines) caches one-body / XC diagonals; reuse them
    ! instead of a second all-k ACE rebuild.
    IF (rdmft_diag_computed .AND. ALLOCATED(rdmft_diag_h) &
        .AND. ALLOCATED(rdmft_diag_vx) .AND. rdmft_xc_n_channels() == 1) THEN
       CALL rdmft_grad_n_from_cached_diag(grad_loc)
    ELSE
       CALL rdmft_grad_n(grad_loc)
    ENDIF
    CALL poolcollect(nbnd, nks, grad_loc, nkstot, grad_mu)
    DEALLOCATE(grad_loc)
    ! When the driver already collected TSM energies, reuse the cache.
    IF (.NOT. rdmft_diag_computed .OR. .NOT. ALLOCATED(rdmft_diag_tsm)) THEN
       CALL rdmft_collect_band_diagnostics(.TRUE.)
    ENDIF
    IF (.NOT. rdmft_diag_computed) RETURN
    !
    ALLOCATE(eps(nbnd, nkstot), eps_plus(nbnd, nkstot), eps_minus(nbnd, nkstot), &
             occ_g(nbnd, nkstot), eps_loc(nbnd, nks), occ_loc(nbnd, nks))
    IF (.NOT. ALLOCATED(rdmft_diag_tsm)) CALL errore('rdmft_compute_dos', &
         'transition-state energies not available', 1)
    eps_loc = rdmft_diag_tsm(:, 1:nks)
    occ_loc = rdmft_n(:, 1:nks)
    CALL poolcollect(nbnd, nks, eps_loc, nkstot, eps)
    CALL poolcollect(nbnd, nks, occ_loc, nkstot, occ_g)
    DEALLOCATE(eps_loc, occ_loc)
    eps_plus = eps
    eps_minus = eps
    !
    ALLOCATE(wk_g(nkstot), isk_g(nkstot))
    ALLOCATE(wk_loc(1, nks), isk_loc(1, nks), wk_col(1, nkstot), isk_col(1, nkstot))
    DO ik = 1, nks
       wk_loc(1, ik)  = wk(ik)
       isk_loc(1, ik) = REAL(isk(ik), DP)
    ENDDO
    CALL poolcollect(1, nks, wk_loc, nkstot, wk_col)
    CALL poolcollect(1, nks, isk_loc, nkstot, isk_col)
    wk_g = wk_col(1, :)
    DO ik = 1, nkstot
       isk_g(ik) = NINT(isk_col(1, ik))
    ENDDO
    DEALLOCATE(wk_loc, isk_loc, wk_col, isk_col)
    !
    ! Koopmans-Fock energies at converged occupations (cheap; reuse cached h/vx).
    ALLOCATE(koop_g(nbnd, nkstot))
    IF (ALLOCATED(rdmft_diag_h) .AND. ALLOCATED(rdmft_diag_vx)) THEN
       BLOCK
          REAL(DP), ALLOCATABLE :: h_loc(:,:), vx_loc(:,:), h_g(:,:), vx_g(:,:)
          ALLOCATE(h_loc(nbnd, nks), vx_loc(nbnd, nks))
          ALLOCATE(h_g(nbnd, nkstot), vx_g(nbnd, nkstot))
          h_loc  = rdmft_diag_h(:, 1:nks)
          vx_loc = rdmft_diag_vx(:, 1:nks)
          CALL poolcollect(nbnd, nks, h_loc,  nkstot, h_g)
          CALL poolcollect(nbnd, nks, vx_loc, nkstot, vx_g)
          DEALLOCATE(h_loc, vx_loc)
          koop_g = h_g + vx_g
          DEALLOCATE(h_g, vx_g)
       END BLOCK
    ELSE
       koop_g = 0.0_DP
    ENDIF
    !
    use_sharma = (INDEX(rdmft_dos_spectral, 'sharma') > 0)
    use_efermig = (INDEX(rdmft_dos_eref, 'efermig') > 0)
    use_brzint = (INDEX(rdmft_dos_integration, 'brzint') > 0)
    IF (.NOT. use_sharma) THEN
       IF (.NOT. use_brzint) THEN
          IF (ionode) WRITE(stdout, '(5X,A)') &
               'RDMFT DOS: elk spectral mode requires brzint; switching.'
          use_brzint = .TRUE.
       END IF
    END IF
    has_sc = noncolin
    occmax = rdmft_dos_occmax()
    do_msum = rdmft_dos_msum
    IF (has_sc) THEN
       ALLOCATE(sc(nbnd, nkstot, 2))
       CALL rdmft_compute_spin_weights(sc, 2)
       BLOCK
          REAL(DP), ALLOCATABLE :: sc_loc(:,:), sc_col(:,:)
          INTEGER :: is_c
          ALLOCATE(sc_loc(nbnd, nks), sc_col(nbnd, nkstot))
          DO is_c = 1, 2
             sc_loc = sc(:, 1:nks, is_c)
             CALL poolcollect(nbnd, nks, sc_loc, nkstot, sc_col)
             sc(:, :, is_c) = sc_col
          ENDDO
          DEALLOCATE(sc_loc, sc_col)
       END BLOCK
    ELSE
       ALLOCATE(sc(1, 1, 1))
    ENDIF
    !
    CALL rdmft_fill_nlmchi(natomwfc, lmax_wfc, nlmchi)
    do_pdos = (natomwfc > 0) .AND. rdmft_dos_pdos
    IF (natomwfc > 0) THEN
       ALLOCATE(proj(natomwfc, nbnd, nkstot))
       proj = 0.0_DP
       CALL rdmft_compute_atomic_proj(proj)
       ALLOCATE(proj_pack(natomwfc * nbnd, nks), proj_col(natomwfc * nbnd, nkstot))
       DO ik = 1, nks
          proj_pack(:, ik) = RESHAPE(proj(1:natomwfc, 1:nbnd, ik), (/ natomwfc * nbnd /))
       ENDDO
       CALL poolcollect(natomwfc * nbnd, nks, proj_pack, nkstot, proj_col)
       DO ik = 1, nkstot
          proj(1:natomwfc, 1:nbnd, ik) = RESHAPE(proj_col(:, ik), (/ natomwfc, nbnd /))
       ENDDO
       DEALLOCATE(proj_pack, proj_col)
    ELSE
       ALLOCATE(proj(1, 1, 1))
       proj = 0.0_DP
       IF (ionode) THEN
          WRITE(stdout, '(/,5X,A)') 'RDMFT DOS: pseudopotentials have no atomic'
          WRITE(stdout, '(5X,A)')    'wave-functions (PP_CHI blocks); skipping'
          WRITE(stdout, '(5X,A)')    'partial DOS.  Total DOS will still be written.'
       END IF
    END IF
    !
    ! Full k-list for brzint MP unfolding (must run on every pool before the
    ! root_pool early return; poolcollect uses inter_pool_comm).
    ALLOCATE(xk_g(3, nkstot))
    CALL poolcollect(3, nks, xk(:, 1:nks), nkstot, xk_g)
    !
    ! Pool-collect above uses inter_pool_comm on every rank; only the first
    ! k-pool root assembles and writes the DOS tables.
    IF (me_pool /= root_pool .OR. my_pool_id /= 0) THEN
       DEALLOCATE(eps, eps_plus, eps_minus, grad_mu, occ_g, koop_g, proj, nlmchi, wk_g, isk_g, sc, xk_g)
       RETURN
    END IF
    !
    ALLOCATE(occ_band(nbnd), dedn_band(nbnd))
    CALL rdmft_dedn_band_averages(grad_mu, occ_g, wk_g, occ_band, dedn_band)
    CALL rdmft_efermi_from_dedn_table(occ_band, dedn_band, &
         ib_homo, mu_rdm, print_msg=.FALSE.)
    !
    IF (TRIM(rdmft_dos_prefix) == ' ') THEN
       fprefix = TRIM(prefix) // '.rdmft'
    ELSE
       fprefix = TRIM(rdmft_dos_prefix)
    END IF
    !
    IF (ionode) THEN
       CALL rdmft_write_dedn(fprefix, grad_mu, occ_g, wk_g, xk_g)
       CALL rdmft_write_efermi_file(fprefix, mu_rdm)
       CALL rdmft_write_band_energies_file(fprefix, eps, koop_g, occ_g, wk_g, isk_g)
    ENDIF
    DEALLOCATE(grad_mu, occ_band, dedn_band)
    !
    ! Unified energy mesh and references for TSM + Koopman total DOS.
    IF (use_efermig) THEN
       degauss = rdmft_dos_degauss
       IF (degauss <= 0.0_DP) degauss = 1.0e-3_DP
       CALL rdmft_compute_dos_efermi(eps, wk_g, isk_g, degauss, e_ref_tsm)
       CALL rdmft_compute_dos_efermi(koop_g, wk_g, isk_g, degauss, e_ref_koop)
    ELSE
       e_ref_tsm  = mu_rdm
       e_ref_koop = mu_rdm
    ENDIF
    CALL rdmft_dos_scan_omega_range(eps, wk_g, nbnd, nkstot, use_sharma, &
         rdmft_dos_occ_weighted, e_ref_tsm, Elw_tsm, Eup_tsm)
    CALL rdmft_dos_scan_omega_range(koop_g, wk_g, nbnd, nkstot, use_sharma, &
         rdmft_dos_occ_weighted, e_ref_koop, Elw_koop, Eup_koop)
    Elw = MIN(Elw_tsm, Elw_koop)
    Eup = MAX(Eup_tsm, Eup_koop)
    degauss = rdmft_dos_degauss
    IF (.NOT. use_brzint .AND. degauss > 0.0_DP) THEN
       Elw = Elw - 3.0_DP * degauss
       Eup = Eup + 3.0_DP * degauss
    END IF
    IF (ABS(rdmft_dos_emin) < E_UNSET - 1.0_DP) THEN
       Emin = rdmft_dos_emin
    ELSE
       Emin = Elw
    ENDIF
    IF (ABS(rdmft_dos_emax) < E_UNSET - 1.0_DP) THEN
       Emax = rdmft_dos_emax
    ELSE
       Emax = Eup
    ENDIF
    IF (rdmft_dos_nwplot >= 2) THEN
       ne = rdmft_dos_nwplot - 1
       IF (ne < 1) CALL errore('rdmft_compute_dos', 'rdmft_dos_nwplot < 2', &
            rdmft_dos_nwplot)
       DeltaE = (Emax - Emin) / DBLE(ne)
    ELSE
       DeltaE = rdmft_dos_deltae
       ne = NINT((Emax - Emin) / DeltaE + 0.500001_DP)
    ENDIF
    IF (use_brzint) THEN
       ie_delta = ne
    ELSEIF (degauss > 0.0_DP .AND. DeltaE > 0.0_DP) THEN
       ie_delta = MAX(1, NINT(5.0_DP * degauss / DeltaE) + 1)
    ELSE
       ie_delta = ne
    END IF
    nsd = 1
    IF (.NOT. rdmft_dos_ssum .AND. (nspin == 2 .OR. noncolin)) nsd = 2
    ALLOCATE(tdos_tsm(0:ne, 2), tdos_koop(0:ne, 2))
    tdos_tsm  = 0.0_DP
    tdos_koop = 0.0_DP
    lmmax = (lmax_wfc + 1)**2 - 1
    IF (do_pdos) THEN
       IF (do_msum) THEN
          ALLOCATE(pdos(0:ne, nat, 0:lmax_wfc, 2))
          pdos = 0.0_DP
       ELSE
          ALLOCATE(pdos_lm(0:ne, nat, 0:lmmax, 2))
          pdos_lm = 0.0_DP
       ENDIF
       ALLOCATE(pdos_sum(0:ne, 2))
       pdos_sum = 0.0_DP
    END IF
    !
    DO im = 1, 2
       IF (im == 1) THEN
          eps_plus = eps
          eps_minus = eps
          e_ref_model = e_ref_tsm
          dos_label = 'TSM'
       ELSE
          eps_plus = koop_g
          eps_minus = koop_g
          e_ref_model = e_ref_koop
          dos_label = 'Koopmans-Fock'
       ENDIF
       do_pdos_run = do_pdos .AND. (im == 1)
       !
       n_branch = 1
       IF (use_sharma .AND. .NOT. rdmft_dos_occ_weighted) n_branch = 2
       ALLOCATE(tdos(0:ne, 2))
       tdos = 0.0_DP
       !
       IF (use_brzint) THEN
       nw_brz = ne + 1
       ld_brz = nbnd
       IF (use_sharma .AND. .NOT. rdmft_dos_occ_weighted) ld_brz = 2 * nbnd
       ALLOCATE(wint(2))
       wint(1) = Emin
       wint(2) = Emax
       IF (ionode) THEN
          WRITE(stdout, '(/,5X,A)') 'RDMFT DOS: brzint BZ integration (' // &
               TRIM(dos_label) // '; may take a while)...'
          WRITE(stdout, '(5X,A,I0,A,I0,A,I0,A,I0)') &
               '  MP grid ', nk1, ' x ', nk2, ' x ', nk3, ',  ngrkf = ', rdmft_dos_ngrkf
       ENDIF
       DO ispn = 1, MERGE(2, 1, nsd == 2)
          IF (ionode .AND. nsd == 2) WRITE(stdout, '(5X,A,I0)') &
               '  brzint spin channel ', ispn
          CALL rdmft_brzint_setup_mp(ngridk, nsk, ivkik, nkr, equiv, ispn, isk_g, &
               xk_g, ierr)
          IF (ierr /= 0) CALL errore('rdmft_compute_dos', &
               'brzint MP unfolding failed (ierr=' // TRIM(int_to_char(ierr)) // &
               '): use automatic K_POINTS on a uniform Monkhorst-Pack mesh and '// &
               'matching vkloff / k1,k2,k3 offsets', ierr)
          ALLOCATE(e_brz(ld_brz, nkr), f_brz(ld_brz, nkr), g_brz(nw_brz))
          e_brz = 0.0_DP
          f_brz = 0.0_DP
          DO ikr = 1, nkr
             ik = equiv(ikr)
             IF (ik <= 0) CYCLE
             DO ib = 1, nbnd
                e_brz(ib, ikr) = rdmft_dos_state_omega(ib, ik, use_sharma, .TRUE., &
                     eps_plus, eps_minus, wk_g, e_ref_model)
                f_brz(ib, ikr) = rdmft_dos_state_weight(ib, ik, ispn, occ_g, wk_g, &
                     sc, has_sc, .TRUE., use_sharma, .TRUE.)
                IF (use_sharma .AND. .NOT. rdmft_dos_occ_weighted) THEN
                   e_brz(nbnd + ib, ikr) = rdmft_dos_state_omega(ib, ik, use_sharma, &
                        .FALSE., eps_plus, eps_minus, wk_g, e_ref_model)
                   f_brz(nbnd + ib, ikr) = rdmft_dos_state_weight(ib, ik, ispn, occ_g, &
                        wk_g, sc, has_sc, .TRUE., use_sharma, .FALSE.)
                END IF
             ENDDO
          ENDDO
          n_brz = ld_brz
          CALL rdmft_brzint(rdmft_dos_nswplot, ngridk, nsk, ivkik, nw_brz, wint, &
               n_brz, ld_brz, e_brz, f_brz, g_brz)
          DO ie_b = 0, ne
             IF (rdmft_dos_ssum) THEN
                tdos(ie_b, 1) = tdos(ie_b, 1) + g_brz(ie_b + 1)
             ELSE
                tdos(ie_b, ispn) = g_brz(ie_b + 1)
             ENDIF
          ENDDO
          IF (do_pdos_run) THEN
             DO nwfc = 1, natomwfc
                IF (ionode .AND. MOD(nwfc, MAX(natomwfc / 10, 1)) == 1) &
                     WRITE(stdout, '(5X,A,I0,A,I0)') &
                     '  brzint PDOS projector ', nwfc, ' / ', natomwfc
                na = nlmchi(nwfc)%na
                l = nlmchi(nwfc)%l
                lm = l**2 + nlmchi(nwfc)%m
                is = ispn
                IF (noncolin) THEN
                   IF (nlmchi(nwfc)%ind <= 2 * l + 1) THEN
                      is = 1
                   ELSE
                      is = 2
                   ENDIF
                   IF (is /= ispn) CYCLE
                ENDIF
                DO ikr = 1, nkr
                   ik = equiv(ikr)
                   IF (ik <= 0) CYCLE
                   DO ib = 1, nbnd
                      e_brz(ib, ikr) = rdmft_dos_state_omega(ib, ik, use_sharma, .TRUE., &
                           eps_plus, eps_minus, wk_g, e_ref_model)
                      f_brz(ib, ikr) = rdmft_dos_state_weight(ib, ik, ispn, occ_g, wk_g, &
                           sc, has_sc, .TRUE., use_sharma, .TRUE.) * proj(nwfc, ib, ik)
                      IF (use_sharma .AND. .NOT. rdmft_dos_occ_weighted) THEN
                         e_brz(nbnd + ib, ikr) = rdmft_dos_state_omega(ib, ik, use_sharma, &
                              .FALSE., eps_plus, eps_minus, wk_g, e_ref_model)
                         f_brz(nbnd + ib, ikr) = rdmft_dos_state_weight(ib, ik, ispn, &
                              occ_g, wk_g, sc, has_sc, .TRUE., use_sharma, .FALSE.) * &
                              proj(nwfc, ib, ik)
                      END IF
                   ENDDO
                ENDDO
                CALL rdmft_brzint(rdmft_dos_nswplot, ngridk, nsk, ivkik, nw_brz, wint, &
                     n_brz, ld_brz, e_brz, f_brz, g_brz)
                DO ie_b = 0, ne
                   IF (do_msum) THEN
                      pdos(ie_b, na, l, ispn) = pdos(ie_b, na, l, ispn) + g_brz(ie_b + 1)
                   ELSE
                      pdos_lm(ie_b, na, lm, ispn) = pdos_lm(ie_b, na, lm, ispn) + g_brz(ie_b + 1)
                   ENDIF
                   IF (rdmft_dos_ssum) THEN
                      pdos_sum(ie_b, 1) = pdos_sum(ie_b, 1) + g_brz(ie_b + 1)
                   ELSE
                      pdos_sum(ie_b, ispn) = pdos_sum(ie_b, ispn) + g_brz(ie_b + 1)
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          DEALLOCATE(e_brz, f_brz, g_brz, ivkik, equiv)
       ENDDO
       DEALLOCATE(wint)
    ELSE
       DO ik = 1, nkstot
          DO ib = 1, nbnd
             DO ib_brz = 1, n_branch
                occupied = (ib_brz == 1)
                IF (use_sharma .AND. rdmft_dos_occ_weighted .AND. .NOT. occupied) CYCLE
                omega_ik = rdmft_dos_state_omega(ib, ik, use_sharma, occupied, &
                     eps_plus, eps_minus, wk_g, e_ref_model)
                IF (noncolin) THEN
                   DO ispn = 1, 2
                      weight = rdmft_dos_state_weight(ib, ik, ispn, occ_g, wk_g, &
                           sc, has_sc, .FALSE., use_sharma, occupied)
                      IF (weight <= 0.0_DP) CYCLE
                      ie_mid = NINT((omega_ik - Emin) / DeltaE)
                      DO ie = MAX(ie_mid - ie_delta, 0), MIN(ie_mid + ie_delta, ne)
                         E = Emin + DeltaE * ie
                         delta = w0gauss((E - omega_ik) / degauss, 0) / degauss
                         IF (rdmft_dos_ssum) THEN
                            tdos(ie, 1) = tdos(ie, 1) + weight * delta
                         ELSE
                            tdos(ie, ispn) = tdos(ie, ispn) + weight * delta
                         ENDIF
                         IF (do_pdos_run) THEN
                            DO nwfc = 1, natomwfc
                               na = nlmchi(nwfc)%na
                               l = nlmchi(nwfc)%l
                               lm = l**2 + nlmchi(nwfc)%m
                               is = ispn
                               IF (nlmchi(nwfc)%ind > 2 * l + 1) is = 2
                               IF (nlmchi(nwfc)%ind <= 2 * l + 1) is = 1
                               IF (is /= ispn) CYCLE
                               pdos_w = weight * proj(nwfc, ib, ik)
                               IF (do_msum) THEN
                                  pdos(ie, na, l, ispn) = pdos(ie, na, l, ispn) + pdos_w * delta
                               ELSE
                                  pdos_lm(ie, na, lm, ispn) = pdos_lm(ie, na, lm, ispn) + &
                                       pdos_w * delta
                               ENDIF
                               IF (rdmft_dos_ssum) THEN
                                  pdos_sum(ie, 1) = pdos_sum(ie, 1) + pdos_w * delta
                               ELSE
                                  pdos_sum(ie, ispn) = pdos_sum(ie, ispn) + pdos_w * delta
                               ENDIF
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ELSE
                   IF (nspin == 2) THEN
                      is = isk_g(ik)
                   ELSE
                      is = 1
                   ENDIF
                   weight = rdmft_dos_state_weight(ib, ik, is, occ_g, wk_g, sc, &
                        has_sc, .FALSE., use_sharma, occupied)
                   IF (weight <= 0.0_DP) CYCLE
                   ie_mid = NINT((omega_ik - Emin) / DeltaE)
                   DO ie = MAX(ie_mid - ie_delta, 0), MIN(ie_mid + ie_delta, ne)
                      E = Emin + DeltaE * ie
                      delta = w0gauss((E - omega_ik) / degauss, 0) / degauss
                      IF (rdmft_dos_ssum) THEN
                         tdos(ie, 1) = tdos(ie, 1) + weight * delta
                      ELSE
                         tdos(ie, is) = tdos(ie, is) + weight * delta
                      ENDIF
                      IF (do_pdos_run) THEN
                         DO nwfc = 1, natomwfc
                            na = nlmchi(nwfc)%na
                            l = nlmchi(nwfc)%l
                            lm = l**2 + nlmchi(nwfc)%m
                            pdos_w = weight * proj(nwfc, ib, ik)
                            IF (do_msum) THEN
                               pdos(ie, na, l, is) = pdos(ie, na, l, is) + pdos_w * delta
                            ELSE
                               pdos_lm(ie, na, lm, is) = pdos_lm(ie, na, lm, is) + &
                                    pdos_w * delta
                            ENDIF
                            IF (rdmft_dos_ssum) THEN
                               pdos_sum(ie, 1) = pdos_sum(ie, 1) + pdos_w * delta
                            ELSE
                               pdos_sum(ie, is) = pdos_sum(ie, is) + pdos_w * delta
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !
    IF (do_pdos_run) THEN
       ALLOCATE(idos(0:ne, 2))
       DO ispn = 1, MERGE(2, 1, rdmft_dos_ssum)
          idos(:, ispn) = tdos(:, ispn) - pdos_sum(:, ispn)
       ENDDO
    END IF
    !
    IF (TRIM(rdmft_dos_prefix) == ' ') THEN
       fprefix = TRIM(prefix) // '.rdmft'
    ELSE
       fprefix = TRIM(rdmft_dos_prefix)
    END IF
    !
    IF (ionode) THEN
       IF (im == 1) THEN
          WRITE(stdout, '(/,5X,A)') 'RDMFT density of states (TSM + Koopmans-Fock):'
          WRITE(stdout, '(5X,A,A)') '  spectral model = ', TRIM(rdmft_dos_spectral)
          WRITE(stdout, '(5X,A,A)') '  integration    = ', TRIM(rdmft_dos_integration)
          IF (use_brzint) THEN
             WRITE(stdout, '(5X,A,I0)') '  ngrkf (brzint) = ', rdmft_dos_ngrkf
          ELSE
             WRITE(stdout, '(5X,A,F8.4,A)') '  degauss (Ry)   = ', rdmft_dos_degauss, ' Ry'
          ENDIF
          IF (use_efermig) THEN
             WRITE(stdout, '(5X,A,F12.6,A)') '  mu (RDM_DEDN)  = ', mu_rdm, ' Ry'
             WRITE(stdout, '(5X,A,F12.6,A)') '  TSM reference  = ', e_ref_tsm, ' Ry'
             WRITE(stdout, '(5X,A,F12.6,A)') '  Koop reference = ', e_ref_koop, ' Ry'
          ELSE IF (use_sharma) THEN
             WRITE(stdout, '(5X,A,F12.6,A,I0,A)') '  mu (RDM_DEDN)  = ', mu_rdm, &
                  ' Ry  (dedn at ib=', ib_homo, ')'
             WRITE(stdout, '(5X,A)') &
                  '  Sharma DOS: omega relative to mu (PRL 110, 116403 Eq. 7)'
          ELSE
             WRITE(stdout, '(5X,A,F12.6,A,I0,A)') '  mu (RDM_DEDN)  = ', mu_rdm, &
                  ' Ry  (dedn at ib=', ib_homo, ')'
          END IF
          WRITE(stdout, '(5X,A,A)') '  ELK files      = ', &
               TRIM(fprefix) // '.RDM_DEDN.OUT, ' // TRIM(fprefix) // '.efermi, ' // &
               TRIM(fprefix) // '.bands, ' // TRIM(fprefix) // '.dos'
          IF (use_sharma) THEN
             IF (rdmft_dos_occ_weighted) THEN
                WRITE(stdout, '(5X,A)') &
                     '  spectral sum   = occupied branch only: n*delta(omega-eps-)'
             ELSE
                WRITE(stdout, '(5X,A)') &
                     '  spectral sum   = n*delta(omega-eps-) + (1-n)*delta(omega+eps+)'
             END IF
          ELSE
             IF (rdmft_dos_occ_weighted) THEN
                WRITE(stdout, '(5X,A)') &
                     '  spectral sum   = ELK dosocc: sc*occsv per state'
             ELSE
                WRITE(stdout, '(5X,A)') &
                     '  spectral sum   = ELK default: sc*occmax per state'
             END IF
          END IF
          WRITE(stdout, '(5X,A,F8.4,A,F8.4,A)') '  window (Ry)    = ', Emin, ' : ', Emax
       END IF
       WRITE(stdout, '(5X,A,A)') '  finished pass: ', TRIM(dos_label)
       !
       IF (do_pdos_run) THEN
          DO na = 1, nat
             nt = ityp(na)
             WRITE(fname, '(A,A,I0.3,A,I0.3)') TRIM(fprefix), '.pdos_at', nt, '_', na
             iunit = find_free_unit()
             OPEN(UNIT=iunit, FILE=TRIM(fname), STATUS='UNKNOWN', FORM='FORMATTED')
             IF (do_msum) THEN
                DO ispn = 1, MERGE(2, 1, nsd == 2)
                   DO l = 0, lmax_wfc
                      WRITE(iunit, '(A,I0,A,I0,A,I0,A,A,A)') '# l=', l, ' spin=', ispn, &
                           ' atom=', na, ' (', TRIM(atm(nt)), ')'
                      DO ie = 0, ne
                         E = Emin + DeltaE * ie
                         IF (ispn == 1) THEN
                            WRITE(iunit, '(2ES16.8)') E, pdos(ie, na, l, ispn)
                         ELSE
                            WRITE(iunit, '(2ES16.8)') E, -pdos(ie, na, l, ispn)
                         ENDIF
                      ENDDO
                      WRITE(iunit, '(A)') ''
                   ENDDO
                ENDDO
             ELSE
                DO ispn = 1, MERGE(2, 1, nsd == 2)
                   DO lm = 0, lmmax
                      WRITE(iunit, '(A,I0,A,I0,A,I0,A,A,A)') '# lm=', lm, ' spin=', ispn, &
                           ' atom=', na, ' (', TRIM(atm(nt)), ')'
                      DO ie = 0, ne
                         E = Emin + DeltaE * ie
                         IF (ispn == 1) THEN
                            WRITE(iunit, '(2ES16.8)') E, pdos_lm(ie, na, lm, ispn)
                         ELSE
                            WRITE(iunit, '(2ES16.8)') E, -pdos_lm(ie, na, lm, ispn)
                         ENDIF
                      ENDDO
                      WRITE(iunit, '(A)') ''
                   ENDDO
                ENDDO
             ENDIF
             CLOSE(iunit)
          ENDDO
          iunit = find_free_unit()
          fname = TRIM(fprefix) // '.idos'
          OPEN(UNIT=iunit, FILE=TRIM(fname), STATUS='UNKNOWN', FORM='FORMATTED')
          WRITE(iunit, '(A)') '# E(Ry)  IDOS  (states/Ry/cell, interstitial)'
          DO ispn = 1, MERGE(2, 1, nsd == 2)
             DO ie = 0, ne
                E = Emin + DeltaE * ie
                IF (ispn == 1) THEN
                   WRITE(iunit, '(2ES16.8)') E, idos(ie, ispn)
                ELSE
                   WRITE(iunit, '(2ES16.8)') E, -idos(ie, ispn)
                ENDIF
             ENDDO
             WRITE(iunit, '(A)') ''
          ENDDO
          CLOSE(iunit)
          WRITE(stdout, '(5X,A)') '  partial DOS written to *.pdos_at<S>_<N> files'
          WRITE(stdout, '(5X,A,A)') '  interstitial DOS written to ', TRIM(fprefix) // '.idos'
       END IF
    END IF
    !
       IF (im == 1) THEN
          tdos_tsm = tdos
       ELSE
          tdos_koop = tdos
       ENDIF
       DEALLOCATE(tdos)
       IF (do_pdos_run) DEALLOCATE(idos)
    ENDDO
    !
    IF (ionode) THEN
       iunit = find_free_unit()
       fname = TRIM(fprefix) // '.dos'
       OPEN(UNIT=iunit, FILE=TRIM(fname), STATUS='UNKNOWN', FORM='FORMATTED')
       IF (lsda .AND. nsd == 2) THEN
          WRITE(iunit, '(A)') &
               '# E(Ry)  DOS_tsm_up  DOS_tsm_down  DOS_koop_up  DOS_koop_down'
          WRITE(iunit, '(A)') '# spin-down columns negative; TSM + Koopmans-Fock on one grid'
       ELSE
          WRITE(iunit, '(A)') '# E(Ry)  DOS_tsm  DOS_koop'
       ENDIF
       DO ie = 0, ne
          E = Emin + DeltaE * ie
          IF (nsd == 2) THEN
             WRITE(iunit, '(5ES16.8)') E, tdos_tsm(ie, 1), -tdos_tsm(ie, 2), &
                  tdos_koop(ie, 1), -tdos_koop(ie, 2)
          ELSE
             WRITE(iunit, '(3ES16.8)') E, tdos_tsm(ie, 1), tdos_koop(ie, 1)
          ENDIF
       ENDDO
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') '  total DOS written to ', TRIM(fname)
    END IF
    !
    IF (do_pdos) THEN
       IF (do_msum) DEALLOCATE(pdos)
       IF (.NOT. do_msum) DEALLOCATE(pdos_lm)
       DEALLOCATE(pdos_sum)
    ENDIF
    DEALLOCATE(tdos_tsm, tdos_koop)
    !
    DEALLOCATE(eps, eps_plus, eps_minus, occ_g, koop_g, proj, nlmchi, wk_g, isk_g, sc, xk_g)
    !
  END SUBROUTINE rdmft_compute_dos
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_magnetization()
    !---------------------------------------------------------------
    !! Compute and write the per-atom and per-species Mulliken
    !! magnetic moments \(M_z = Q_\uparrow - Q_\downarrow\) from
    !! the RDMFT natural occupations projected onto the atomic
    !! wave-function basis.  Writes \texttt{<prefix>.rdmft.mag} and
    !! prints the table on \texttt{stdout}.
    !!
    !! Inactive for non-spin-polarized (\texttt{lsda = .false.}) runs
    !! and for pseudopotentials with no atomic wave-functions
    !! (\texttt{n\_atom\_wfc} returns zero).  Also inactive in the
    !! noncollinear regime (the sigma-resolved spinor PDOS / Mz is
    !! deferred to the standalone \texttt{projwfc.x} post-processor;
    !! see the README notes).  For the lsda + no-PP_CHI corner case
    !! the routine still prints the global \(M_z\) computed directly
    !! from the natural occupations and spin indices.
    USE io_global,        ONLY : stdout, ionode
    USE io_files,         ONLY : prefix
    USE wvfct,            ONLY : nbnd
    USE klist,            ONLY : nks, nkstot, wk
    USE lsda_mod,         ONLY : isk, lsda
    USE noncollin_module, ONLY : noncolin
    USE ions_base,        ONLY : nat, ityp, atm, nsp
    USE mp_pools,         ONLY : me_pool, root_pool
    USE rdmft_module,     ONLY : rdmft_n, rdmft_compute_dos_flag => rdmft_compute_dos, &
                                  rdmft_dos_prefix
    USE rdmft_energy,     ONLY : rdmft_collect_band_diagnostics, rdmft_diag_computed
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    REAL(DP), ALLOCATABLE :: occ_g(:,:), proj(:,:,:)
    REAL(DP), ALLOCATABLE :: q_up(:), q_dn(:), mz_site(:), mz_species(:)
    REAL(DP), ALLOCATABLE :: wk_g(:)
    INTEGER,  ALLOCATABLE :: isk_g(:)
    TYPE(wfc_label), ALLOCATABLE :: nlmchi(:)
    !
    REAL(DP) :: occ_w, mz_global
    CHARACTER(LEN=256) :: fprefix, fname
    INTEGER :: natomwfc, lmax_wfc, ik, ib, nwfc, na, nt, iunit
    LOGICAL :: do_mag
    !
    IF (.NOT. rdmft_compute_dos_flag) RETURN
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    !
    ! Mz is only meaningful for collinear spin-polarised runs.  Skip
    ! silently for nspin = 1 and noncolin; the latter also requires
    ! a spinor-aware projector.
    IF (noncolin .OR. .NOT. lsda) RETURN
    !
    CALL rdmft_collect_band_diagnostics(.TRUE.)
    IF (.NOT. rdmft_diag_computed) RETURN
    !
    ALLOCATE(occ_g(nbnd, nkstot))
    occ_g(:, 1:nks) = rdmft_n(:, 1:nks)
    CALL poolrecover(occ_g, nbnd, nkstot, nks)
    !
    ALLOCATE(wk_g(nkstot), isk_g(nkstot))
    wk_g = 0.0_DP
    isk_g = 1
    wk_g(1:nks)  = wk(1:nks)
    isk_g(1:nks) = isk(1:nks)
    CALL poolrecover(wk_g, 1, nkstot, nks)
    CALL ipoolrecover(isk_g, 1, nkstot, nks)
    !
    CALL rdmft_fill_nlmchi(natomwfc, lmax_wfc, nlmchi)
    do_mag = (natomwfc > 0)
    IF (do_mag) THEN
       ALLOCATE(proj(natomwfc, nbnd, nkstot))
       proj = 0.0_DP
       CALL rdmft_compute_atomic_proj(proj)
       CALL poolrecover(proj, natomwfc * nbnd, nkstot, nks)
    ELSE
       ALLOCATE(proj(1, 1, 1))
       proj = 0.0_DP
       IF (ionode) THEN
          WRITE(stdout, '(/,5X,A)') 'RDMFT Mz: pseudopotentials have no atomic'
          WRITE(stdout, '(5X,A)')    'wave-functions (PP_CHI blocks); skipping'
          WRITE(stdout, '(5X,A)')    'per-atom Mulliken magnetization.'
       END IF
    END IF
    !
    IF (me_pool /= root_pool) THEN
       DEALLOCATE(occ_g, proj, nlmchi, wk_g, isk_g)
       RETURN
    END IF
    !
    ALLOCATE(q_up(nat), q_dn(nat), mz_site(nat), mz_species(nsp))
    q_up = 0.0_DP
    q_dn = 0.0_DP
    IF (do_mag) THEN
       DO ik = 1, nkstot
          DO ib = 1, nbnd
             occ_w = wk_g(ik) * occ_g(ib, ik)
             DO nwfc = 1, natomwfc
                na = nlmchi(nwfc)%na
                IF (isk_g(ik) == 1) THEN
                   q_up(na) = q_up(na) + occ_w * proj(nwfc, ib, ik)
                ELSE
                   q_dn(na) = q_dn(na) + occ_w * proj(nwfc, ib, ik)
                END IF
             END DO
          END DO
       END DO
    END IF
    mz_site = q_up - q_dn
    mz_species = 0.0_DP
    DO na = 1, nat
       mz_species(ityp(na)) = mz_species(ityp(na)) + mz_site(na)
    END DO
    !
    IF (TRIM(rdmft_dos_prefix) == ' ') THEN
       fprefix = TRIM(prefix) // '.rdmft'
    ELSE
       fprefix = TRIM(rdmft_dos_prefix)
    END IF
    !
    IF (ionode) THEN
       IF (do_mag) THEN
          iunit = find_free_unit()
          fname = TRIM(fprefix) // '.mag'
          OPEN(UNIT=iunit, FILE=TRIM(fname), STATUS='UNKNOWN', FORM='FORMATTED')
          WRITE(iunit, '(A)') '# atom  species  Mz(electrons)  Q_up  Q_dn'
          DO na = 1, nat
             WRITE(iunit, '(I6,2X,I4,2X,3F14.8)') na, ityp(na), mz_site(na), q_up(na), q_dn(na)
          END DO
          WRITE(iunit, '(A)') '# species totals (Mz):'
          DO nt = 1, nsp
             WRITE(iunit, '(A,A,A,2X,F14.8)') '# ', TRIM(atm(nt)), ' sum Mz = ', mz_species(nt)
          END DO
          CLOSE(iunit)
          !
          WRITE(stdout, '(/,5X,A)') 'Per-atom magnetization (Mulliken, electron units):'
          WRITE(stdout, '(5X,A)') '  atom  species      Mz        Q_up       Q_dn'
          DO na = 1, nat
             WRITE(stdout, '(5X,I4,2X,A4,3F12.6)') na, atm(ityp(na)), &
                  mz_site(na), q_up(na), q_dn(na)
          END DO
          WRITE(stdout, '(5X,A)') 'Species-summed Mz:'
          DO nt = 1, nsp
             WRITE(stdout, '(5X,A,A,A,F12.6)') '  ', TRIM(atm(nt)), ' : ', mz_species(nt)
          END DO
          mz_global = rdmft_dos_global_mz(occ_g, wk_g, isk_g, nbnd, nkstot)
          WRITE(stdout, '(5X,A,F12.6)') 'Global Mz (from occ_g) = ', mz_global
          WRITE(stdout, '(5X,A,F12.6)') 'Sum of site Mz         = ', SUM(mz_site)
          WRITE(stdout, '(5X,A,A)') '  magnetization table written to ', TRIM(fprefix) // '.mag'
       ELSE
          mz_global = rdmft_dos_global_mz(occ_g, wk_g, isk_g, nbnd, nkstot)
          WRITE(stdout, '(/,5X,A,F12.6)') 'RDMFT Global Mz (from occ_g) = ', mz_global
       END IF
    END IF
    !
    DEALLOCATE(occ_g, proj, nlmchi, q_up, q_dn, mz_site, mz_species)
    DEALLOCATE(wk_g, isk_g)
    !
  END SUBROUTINE rdmft_compute_magnetization
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_dos_global_mz(occ, wk_g, isk_g, nb, nk)
    !---------------------------------------------------------------
    !! Global magnetization \(\sum_{k,\sigma} \pm w_k \sum_i n_{ik}\)
    !! from the already pool-GATHERED (1:nkstot) occupation / weight /
    !! spin arrays held on ionode.  Computed as a plain local sum --
    !! it must NOT call any pool reduction, because this routine runs
    !! inside an ``IF (ionode)`` block; an ``mp_sum(inter_pool_comm)``
    !! there would be entered by ionode only and deadlock every other
    !! k-point pool (the cause of the RDMFT-DOS hang under npool > 1).
    REAL(DP), INTENT(IN) :: occ(nb, nk), wk_g(nk)
    INTEGER,  INTENT(IN) :: isk_g(nk), nb, nk
    INTEGER  :: ib, ik
    REAL(DP) :: s
    s = 0.0_DP
    DO ik = 1, nk
       DO ib = 1, nb
          IF (isk_g(ik) == 1) THEN
             s = s + wk_g(ik) * occ(ib, ik)
          ELSE
             s = s - wk_g(ik) * occ(ib, ik)
          END IF
       END DO
    END DO
    rdmft_dos_global_mz = s
  END FUNCTION rdmft_dos_global_mz
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_report_chemical_potential()
    !---------------------------------------------------------------
    !! Print RDMFT chemical potential \(\mu=\texttt{dedn}(i_\mu)\) from
    !! the converged \texttt{RDM\_DEDN.OUT} band table.
    USE io_global,    ONLY : stdout, ionode
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks, nkstot, wk
    USE rdmft_module, ONLY : rdmft_n
    USE rdmft_energy, ONLY : rdmft_grad_n, rdmft_grad_n_from_cached_diag, &
                              rdmft_collect_band_diagnostics, rdmft_diag_h, &
                              rdmft_diag_vx, rdmft_diag_computed
    USE rdmft_xc,     ONLY : rdmft_xc_n_channels
    !
    REAL(DP), ALLOCATABLE :: grad_loc(:,:), grad_g(:,:), occ_g(:,:), wk_g(:)
    REAL(DP), ALLOCATABLE :: wk_loc(:,:), wk_col(:,:)
    REAL(DP), ALLOCATABLE :: occ_band(:), dedn_band(:)
    REAL(DP) :: mu
    INTEGER :: ib_mu, ik
    !
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    !
    CALL rdmft_collect_band_diagnostics(.FALSE.)
    !
    ALLOCATE(grad_loc(nbnd, nks), grad_g(nbnd, nkstot), occ_g(nbnd, nkstot))
    ALLOCATE(wk_g(nkstot), wk_loc(1, nks), wk_col(1, nkstot))
    IF (rdmft_diag_computed .AND. ALLOCATED(rdmft_diag_h) &
        .AND. ALLOCATED(rdmft_diag_vx) .AND. rdmft_xc_n_channels() == 1) THEN
       CALL rdmft_grad_n_from_cached_diag(grad_loc)
    ELSE
       CALL rdmft_grad_n(grad_loc)
    ENDIF
    CALL poolcollect(nbnd, nks, grad_loc, nkstot, grad_g)
    CALL poolcollect(nbnd, nks, rdmft_n, nkstot, occ_g)
    DO ik = 1, nks
       wk_loc(1, ik) = wk(ik)
    ENDDO
    CALL poolcollect(1, nks, wk_loc, nkstot, wk_col)
    wk_g = wk_col(1, :)
    DEALLOCATE(grad_loc, wk_loc, wk_col)
    !
    ALLOCATE(occ_band(nbnd), dedn_band(nbnd))
    CALL rdmft_dedn_band_averages(grad_g, occ_g, wk_g, occ_band, dedn_band)
    CALL rdmft_efermi_from_dedn_table(occ_band, dedn_band, ib_mu, mu)
    !
    DEALLOCATE(grad_g, occ_g, wk_g, occ_band, dedn_band)
    !
  END SUBROUTINE rdmft_report_chemical_potential
  !
END MODULE rdmft_dos
