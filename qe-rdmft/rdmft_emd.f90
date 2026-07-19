!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_emd
  !------------------------------------------------------------------
  !! RDMFT electron momentum density (EMD), Compton profiles, and
  !! directional autocorrelation functions.
  !!
  !! Port of ELK tasks 170--171 (\texttt{writeemd.f90},
  !! \texttt{emdplot1d.f90}, \texttt{rfhkintp.f90}) adapted to the
  !! PW plane-wave representation of natural orbitals.  Adds
  !! physical-\(q\) Compton output (\texttt{*.emd\_compton}) and
  !! \(B(r)\) via cosine transform (\texttt{*.emd\_af}).
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_compute_emd, rdmft_emd_plot1d, rdmft_emd_compton, rdmft_emd_autocorr
  !
  REAL(DP), PARAMETER :: EMD_EPS = 1.0e-5_DP
  !
  ! Cached momentum density for plotting (filled by rdmft_compute_emd).
  INTEGER,  SAVE :: emd_nspin_ch = 1
  REAL(DP), ALLOCATABLE, SAVE :: emd_pcryst(:,:,:), emd_val(:,:,:)
  REAL(DP), ALLOCATABLE, SAVE :: emd_val_band(:,:,:,:)
  INTEGER,  SAVE :: emd_ib_start = 1
  INTEGER,  SAVE :: emd_ib_end = 0
  LOGICAL,  SAVE :: emd_band_cache_valid = .FALSE.
  INTEGER,  SAVE :: emd_ngridk(3) = 0
  INTEGER,  ALLOCATABLE, SAVE :: emd_ivkik(:,:,:), emd_equiv(:)
  REAL(DP), ALLOCATABLE, SAVE :: emd_xkg(:,:)
  INTEGER,  SAVE :: emd_nkr = 0
  REAL(DP), SAVE :: emd_pmax = 0.0_DP
  REAL(DP), SAVE :: emd_pmax_pw = 0.0_DP
  REAL(DP), SAVE :: emd_ecutwfc = 0.0_DP
  REAL(DP), SAVE :: emd_kf_total = 0.0_DP
  REAL(DP), SAVE :: emd_kf_spin(2) = 0.0_DP
  REAL(DP), SAVE :: emd_ne_expected = 0.0_DP
  LOGICAL,  SAVE :: emd_kf_total_valid = .FALSE.
  LOGICAL,  SAVE :: emd_kf_spin_valid(2) = .FALSE.
  !
  ! Cached directional Compton profile J(q) for compton + autocorr reuse.
  INTEGER,  SAVE :: emd_cp_nq = 0
  INTEGER,  SAVE :: emd_cp_nch = 0
  REAL(DP), SAVE :: emd_cp_qmin = 0.0_DP
  REAL(DP), SAVE :: emd_cp_qmax = 0.0_DP
  REAL(DP), SAVE :: emd_cp_dq = 0.0_DP
  REAL(DP), SAVE :: emd_cp_qscale = 0.0_DP
  REAL(DP), SAVE :: emd_cp_dir(3) = 0.0_DP
  REAL(DP), ALLOCATABLE, SAVE :: emd_cp_q(:), emd_cp_j(:,:)
  LOGICAL,  SAVE :: emd_cp_valid = .FALSE.
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE r3cross(a, b, c)
    REAL(DP), INTENT(IN)  :: a(3), b(3)
    REAL(DP), INTENT(OUT) :: c(3)
    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  END SUBROUTINE r3cross
  !
  !-----------------------------------------------------------------
  SUBROUTINE axangsu2(axis, angle, su2)
    REAL(DP), INTENT(IN)  :: axis(3)
    REAL(DP), INTENT(IN)  :: angle
    COMPLEX(DP), INTENT(OUT) :: su2(2, 2)
    REAL(DP) :: x, y, zz, cs, sn, t1
    x = axis(1)
    y = axis(2)
    zz = axis(3)
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
    cs = DCOS(0.5_DP * angle)
    sn = DSIN(0.5_DP * angle)
    su2(1, 1) = CMPLX(cs, -zz * sn, KIND=DP)
    su2(2, 1) = CMPLX(y * sn, -x * sn, KIND=DP)
    su2(1, 2) = CMPLX(-y * sn, -x * sn, KIND=DP)
    su2(2, 2) = CMPLX(cs, zz * sn, KIND=DP)
  END SUBROUTINE axangsu2
  !
  !-----------------------------------------------------------------
  SUBROUTINE sqasu2_emd(sqaxis, su2, tsqaz)
    REAL(DP), INTENT(IN)  :: sqaxis(3)
    COMPLEX(DP), INTENT(OUT) :: su2(2, 2)
    LOGICAL, INTENT(OUT) :: tsqaz
    REAL(DP) :: v1(3), v2(3), v3(3), th, t1
    v1 = sqaxis
    t1 = DSQRT(v1(1)**2 + v1(2)**2 + v1(3)**2)
    IF (t1 <= 1.0e-8_DP) CALL errore('sqasu2_emd', &
         'spin-quantisation axis (sqaxis) has zero length', 1)
    v1 = v1 / t1
    IF (ABS(v1(3) - 1.0_DP) < 1.0e-8_DP) THEN
       tsqaz = .TRUE.
       su2 = RESHAPE((/ (1.0_DP, 0.0_DP), (0.0_DP, 0.0_DP), &
            (0.0_DP, 0.0_DP), (1.0_DP, 0.0_DP) /), (/ 2, 2 /))
    ELSE
       tsqaz = .FALSE.
       v2 = (/ 0.0_DP, 0.0_DP, 1.0_DP /)
       CALL r3cross(v1, v2, v3)
       th = -DACOS(v1(3))
       CALL axangsu2(v3, th, su2)
    END IF
  END SUBROUTINE sqasu2_emd
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION emd_band_weight(ik, ib, ispn, weight, su2, tsqaz, &
       evc_up, evc_dn, ig)
    USE rdmft_module, ONLY : rdmft_emd_spin_sum
    INTEGER, INTENT(IN) :: ik, ib, ispn, ig
    REAL(DP), INTENT(IN) :: weight
    COMPLEX(DP), INTENT(IN) :: su2(2, 2)
    LOGICAL, INTENT(IN) :: tsqaz
    COMPLEX(DP), INTENT(IN) :: evc_up, evc_dn
    REAL(DP) :: nup, ndn, n1, n2
    nup = ABS(evc_up)**2
    ndn = ABS(evc_dn)**2
    IF (rdmft_emd_spin_sum) THEN
       emd_band_weight = weight * (nup + ndn)
    ELSE IF (tsqaz) THEN
       IF (ispn == 1) THEN
          emd_band_weight = weight * nup
       ELSE
          emd_band_weight = weight * ndn
       END IF
    ELSE
       n1 = REAL(su2(1, 1) * evc_up + su2(1, 2) * evc_dn, DP)**2 + &
            AIMAG(su2(1, 1) * evc_up + su2(1, 2) * evc_dn)**2
       n2 = REAL(su2(2, 1) * evc_up + su2(2, 2) * evc_dn, DP)**2 + &
            AIMAG(su2(2, 1) * evc_up + su2(2, 2) * evc_dn)**2
       IF (ispn == 1) THEN
          emd_band_weight = weight * n1
       ELSE
          emd_band_weight = weight * n2
       END IF
    END IF
  END FUNCTION emd_band_weight
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION momentum_cryst_diff(a, b)
    REAL(DP), INTENT(IN) :: a(3), b(3)
    REAL(DP) :: d(3)
    d = a - b
    momentum_cryst_diff = DSQRT(d(1)**2 + d(2)**2 + d(3)**2)
  END FUNCTION momentum_cryst_diff
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION periodic_cryst_diff(a, b)
    REAL(DP), INTENT(IN) :: a(3), b(3)
    REAL(DP) :: d(3)
    INTEGER :: i
    d = a - b
    DO i = 1, 3
       d(i) = d(i) - NINT(d(i))
    END DO
    periodic_cryst_diff = DSQRT(d(1)**2 + d(2)**2 + d(3)**2)
  END FUNCTION periodic_cryst_diff
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_free_compton_cache()
    IF (ALLOCATED(emd_cp_q)) DEALLOCATE(emd_cp_q)
    IF (ALLOCATED(emd_cp_j)) DEALLOCATE(emd_cp_j)
    emd_cp_nq = 0
    emd_cp_nch = 0
    emd_cp_valid = .FALSE.
  END SUBROUTINE emd_free_compton_cache
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_free_cache()
    IF (ALLOCATED(emd_pcryst)) DEALLOCATE(emd_pcryst)
    IF (ALLOCATED(emd_val)) DEALLOCATE(emd_val)
    IF (ALLOCATED(emd_val_band)) DEALLOCATE(emd_val_band)
    IF (ALLOCATED(emd_ivkik)) DEALLOCATE(emd_ivkik)
    IF (ALLOCATED(emd_equiv)) DEALLOCATE(emd_equiv)
    IF (ALLOCATED(emd_xkg)) DEALLOCATE(emd_xkg)
    CALL emd_free_compton_cache()
    emd_ngridk = 0
    emd_nkr = 0
    emd_band_cache_valid = .FALSE.
    emd_kf_total = 0.0_DP
    emd_kf_spin = 0.0_DP
    emd_kf_total_valid = .FALSE.
    emd_kf_spin_valid = .FALSE.
  END SUBROUTINE emd_free_cache
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_compute_kf(ib_start, ib_end, ne_expected)
    !---------------------------------------------------------------
    !! Selected-band equivalent Fermi momentum from integrated density.
    !
    USE klist,      ONLY : nks, wk
    USE lsda_mod,   ONLY : lsda, nspin, isk
    USE cell_base,  ONLY : omega
    USE constants,  ONLY : pi
    USE rdmft_module, ONLY : rdmft_n
    !
    INTEGER, INTENT(IN) :: ib_start, ib_end
    REAL(DP), INTENT(OUT) :: ne_expected
    !
    INTEGER :: ik, ispn
    REAL(DP) :: ne_spin(2), rho_tot, rho_spin, kf_eq
    !
    ne_expected = SUM(wk(1:nks) * SUM(rdmft_n(ib_start:ib_end, 1:nks), 1))
    emd_ne_expected = ne_expected
    emd_kf_total_valid = .FALSE.
    emd_kf_spin_valid = .FALSE.
    emd_kf_total = 0.0_DP
    emd_kf_spin = 0.0_DP
    IF (omega <= 1.0e-16_DP) RETURN
    rho_tot = ne_expected / omega
    IF (rho_tot > 1.0e-16_DP) THEN
       emd_kf_total = (3.0_DP * pi**2 * rho_tot)**(1.0_DP / 3.0_DP)
       emd_kf_total_valid = .TRUE.
    END IF
    IF (lsda .AND. nspin == 2) THEN
       ne_spin = 0.0_DP
       DO ik = 1, nks
          IF (isk(ik) >= 1 .AND. isk(ik) <= 2) &
               ne_spin(isk(ik)) = ne_spin(isk(ik)) + &
               wk(ik) * SUM(rdmft_n(ib_start:ib_end, ik))
       END DO
       DO ispn = 1, 2
          rho_spin = ne_spin(ispn) / omega
          IF (rho_spin > 1.0e-16_DP) THEN
             emd_kf_spin(ispn) = (6.0_DP * pi**2 * rho_spin)**(1.0_DP / 3.0_DP)
             emd_kf_spin_valid(ispn) = .TRUE.
          END IF
       END DO
       IF (ALL(emd_kf_spin_valid)) THEN
          kf_eq = ((emd_kf_spin(1)**3 + emd_kf_spin(2)**3) / 2.0_DP)**(1.0_DP / 3.0_DP)
          emd_kf_total = kf_eq
          emd_kf_total_valid = .TRUE.
       END IF
    END IF
  END SUBROUTINE emd_compute_kf
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_emd_setup_mp(ierr)
    USE io_global,        ONLY : stdout, ionode
    USE klist,            ONLY : nks, nkstot, xk
    USE lsda_mod,         ONLY : isk
    USE mp_pools,         ONLY : npool
    USE rdmft_brzint_mod, ONLY : rdmft_brzint_setup_mp, rdmft_brzint_vkoff
    INTEGER, INTENT(OUT) :: ierr
    REAL(DP), ALLOCATABLE :: xk_g(:,:), isk_loc(:,:), isk_col(:,:)
    INTEGER, ALLOCATABLE :: isk_g(:)
    REAL(DP) :: vkoff(3)
    INTEGER :: ik, ispin, nsk(3), i, j, k, n
    !
    ierr = 0
    CALL emd_free_cache()
    IF (npool > 1) THEN
       IF (ionode) WRITE(stdout, '(5X,A)') &
            'RDMFT EMD: k-pool parallelization is not supported'
       ierr = 3
       RETURN
    ENDIF
    ALLOCATE(xk_g(3, nkstot), isk_loc(1, nks), isk_col(1, nkstot), isk_g(nkstot))
    CALL poolcollect(3, nks, xk, nkstot, xk_g)
    DO ik = 1, nks
       isk_loc(1, ik) = REAL(isk(ik), DP)
    END DO
    CALL poolcollect(1, nks, isk_loc, nkstot, isk_col)
    DO ik = 1, nkstot
       isk_g(ik) = NINT(isk_col(1, ik))
    END DO
    ispin = 1
    CALL rdmft_brzint_setup_mp(emd_ngridk, nsk, emd_ivkik, emd_nkr, emd_equiv, &
         ispin, isk_g, xk_g, ierr)
    IF (ierr /= 0) THEN
       IF (ionode) WRITE(stdout, '(5X,A,I0)') &
            'RDMFT EMD: could not unfold k-mesh onto MP grid, ierr=', ierr
       DEALLOCATE(xk_g, isk_loc, isk_col, isk_g)
       CALL emd_free_cache()
       RETURN
    ENDIF
    CALL rdmft_brzint_vkoff(vkoff)
    ALLOCATE(emd_xkg(3, emd_nkr))
    DO i = 1, emd_ngridk(1)
       DO j = 1, emd_ngridk(2)
          DO k = 1, emd_ngridk(3)
             n = (k - 1) + (j - 1) * emd_ngridk(3) + (i - 1) * emd_ngridk(2) * emd_ngridk(3) + 1
             emd_xkg(1, n) = DBLE(i - 1) / emd_ngridk(1) + vkoff(1)
             emd_xkg(2, n) = DBLE(j - 1) / emd_ngridk(2) + vkoff(2)
             emd_xkg(3, n) = DBLE(k - 1) / emd_ngridk(3) + vkoff(3)
          END DO
       END DO
    END DO
    DEALLOCATE(xk_g, isk_loc, isk_col, isk_g)
  END SUBROUTINE rdmft_emd_setup_mp
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_emd()
    !---------------------------------------------------------------
    USE io_global,        ONLY : stdout, ionode
    USE io_files,         ONLY : prefix, nwordwfc, iunwfc
    USE wvfct,            ONLY : nbnd, npwx, current_k
    USE klist,            ONLY : nks, nkstot, ngk, igk_k, xk, wk
    USE wavefunctions,    ONLY : evc
    USE lsda_mod,         ONLY : nspin, lsda, isk
    USE noncollin_module, ONLY : noncolin
    USE gvect,            ONLY : g
    USE gvecw,            ONLY : ecutwfc, gcutw
    USE cell_base,        ONLY : tpiba, at, bg
    USE buffers,          ONLY : get_buffer
    USE rdmft_module,        ONLY : rdmft_n, rdmft_emd_prefix, &
                                 rdmft_emd_write_elk_fmt, rdmft_emd_spin_sum, &
                                 rdmft_emd_band_min, rdmft_emd_band_max, &
                                 rdmft_emd_write_band_resolved, rdmft_emd_do_compton, &
                                 rdmft_emd_sqaxis, &
                                 rdmft_emd_pmax_kf_factor
    !
    CHARACTER(LEN=256) :: fprefix, fname
    INTEGER :: ik, ig, ib, ispn, igl, iunit, iunit_band, ios, ierr, nch
    INTEGER :: ib_start, ib_end
    REAL(DP) :: g2, pmag, weight, sm, ne_check, ne_expected, contrib, gpc_cart(3)
    REAL(DP) :: pcryst(3)
    COMPLEX(DP) :: su2(2, 2)
    LOGICAL :: tsqaz, band_unit_open, cache_band
    REAL(DP), ALLOCATABLE :: emd_loc(:,:), band_ch(:,:)
    !
    emd_ecutwfc = ecutwfc
    emd_pmax_pw = DSQRT(gcutw) * tpiba
    IF (rdmft_emd_spin_sum) THEN
       nch = 1
    ELSE IF (noncolin .OR. (lsda .AND. nspin == 2)) THEN
       nch = 2
    ELSE
       nch = 1
    ENDIF
    emd_nspin_ch = nch
    ib_start = MAX(rdmft_emd_band_min, 1)
    IF (rdmft_emd_band_max <= 0) THEN
       ib_end = nbnd
    ELSE
       ib_end = MIN(rdmft_emd_band_max, nbnd)
    ENDIF
    IF (ib_start > nbnd .OR. ib_start > ib_end) CALL errore('rdmft_compute_emd', &
         'invalid EMD band window: require 1 <= band_min <= band_max <= nbnd', 1)
    CALL emd_compute_kf(ib_start, ib_end, ne_expected)
    IF (emd_kf_total_valid .AND. rdmft_emd_pmax_kf_factor > 0.0_DP) THEN
       emd_pmax = rdmft_emd_pmax_kf_factor * emd_kf_total
    ELSE
       emd_pmax = emd_pmax_pw
    END IF
    IF (ionode) THEN
       WRITE(stdout, '(5X,A,F8.4,A)') 'RDMFT EMD: ecutwfc = ', emd_ecutwfc, ' Ry'
       WRITE(stdout, '(5X,A,F8.4,A)') 'RDMFT EMD: |p|_max (PW sphere) = ', &
            emd_pmax_pw, ' Bohr^-1'
       IF (emd_kf_total_valid .AND. rdmft_emd_pmax_kf_factor > 0.0_DP) THEN
          WRITE(stdout, '(5X,A,F8.4,A,F8.4,A,F8.4,A)') &
               'RDMFT EMD: |p|_max = ', emd_pmax, ' Bohr^-1  (', &
               rdmft_emd_pmax_kf_factor, ' * k_F)'
       ELSE
          WRITE(stdout, '(5X,A,F8.4,A)') 'RDMFT EMD: |p|_max = ', emd_pmax, &
               ' Bohr^-1  (from ecutwfc)'
       END IF
       IF (emd_kf_total_valid) WRITE(stdout, '(5X,A,F12.6,A)') &
            'RDMFT EMD: k_F (selected-band equivalent) = ', emd_kf_total, ' Bohr^-1'
       IF (lsda .AND. nspin == 2) THEN
          IF (emd_kf_spin_valid(1)) WRITE(stdout, '(5X,A,F12.6,A)') &
               'RDMFT EMD: k_F(up) = ', emd_kf_spin(1), ' Bohr^-1'
          IF (emd_kf_spin_valid(2)) WRITE(stdout, '(5X,A,F12.6,A)') &
               'RDMFT EMD: k_F(dn) = ', emd_kf_spin(2), ' Bohr^-1'
       END IF
       IF (emd_pmax > emd_pmax_pw + 1.0e-8_DP) WRITE(stdout, '(5X,A)') &
            '  WARNING: |p|_max exceeds ecutwfc sphere; high-|p| tails may be truncated'
    ENDIF
    IF (ionode) THEN
       WRITE(stdout, '(5X,A,I0,A,I0,A,I0,A)') &
            'RDMFT EMD: band window = ', ib_start, ':', ib_end, ' of ', nbnd, ' bands'
       IF (rdmft_emd_write_band_resolved) &
            WRITE(stdout, '(5X,A)') 'RDMFT EMD: writing band-resolved EMD'
       IF (rdmft_emd_do_compton) &
            WRITE(stdout, '(5X,A)') 'RDMFT EMD: caching band-resolved EMD for Compton'
    ENDIF
    cache_band = rdmft_emd_write_band_resolved .OR. rdmft_emd_do_compton
    emd_band_cache_valid = .FALSE.
    IF (TRIM(rdmft_emd_prefix) == ' ') THEN
       fprefix = TRIM(prefix) // '.rdmft'
    ELSE
       fprefix = TRIM(rdmft_emd_prefix)
    ENDIF
    !
    CALL sqasu2_emd(rdmft_emd_sqaxis, su2, tsqaz)
    CALL rdmft_emd_setup_mp(ierr)
    emd_ne_expected = ne_expected
    emd_ib_start = ib_start
    emd_ib_end = ib_end
    IF (ierr == 3) CALL errore('rdmft_compute_emd', &
         'k-pool parallelization not supported; run rdmft_emd.x with -nk 1', 3)
    IF (ierr /= 0) CALL errore('rdmft_compute_emd', &
         'could not unfold symmetry-reduced k-mesh onto the full MP grid (try nosym=.true. in pw.x)', ierr)
    !
    ALLOCATE(emd_pcryst(3, npwx, nkstot), emd_val(npwx, nkstot, nch))
    emd_pcryst = 0.0_DP
    emd_val = 0.0_DP
    IF (cache_band) THEN
       ALLOCATE(emd_val_band(npwx, nkstot, nch, ib_end - ib_start + 1))
       emd_val_band = 0.0_DP
    END IF
    ALLOCATE(emd_loc(npwx, nch))
    !
    ALLOCATE(band_ch(nbnd, nch))
    band_unit_open = .FALSE.
    IF (rdmft_emd_write_band_resolved .AND. ionode) THEN
       fname = TRIM(fprefix) // '.emd_band'
       OPEN(NEWUNIT=iunit_band, FILE=TRIM(fname), STATUS='REPLACE', FORM='FORMATTED', &
            IOSTAT=ios)
       IF (ios /= 0) CALL errore('rdmft_compute_emd', 'opening '//TRIM(fname), ABS(ios))
       band_unit_open = .TRUE.
       WRITE(iunit_band, '(A)') '# RDMFT band-resolved electron momentum density n_ib(p)'
       WRITE(iunit_band, '(A,I0,A,I0,A,I0)') '# band window = ', ib_start, ':', ib_end, &
            ' of ', nbnd
       WRITE(iunit_band, '(A)') '# n_ib(p) excludes k-weights w_k; integrate as sum_k w_k sum_G n_ib(p)'
       WRITE(iunit_band, '(A)') '# ik ib ig p1 p2 p3 (crystal reciprocal) |p| (Bohr^-1) n_ib(p) [n_up n_dn]'
    ENDIF
    !
    DO ik = 1, nks
       current_k = ik
       CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       emd_loc = 0.0_DP
       DO ig = 1, ngk(ik)
          band_ch = 0.0_DP
          igl = igk_k(ig, ik)
          gpc_cart(1) = xk(1, ik) + g(1, igl)
          gpc_cart(2) = xk(2, ik) + g(2, igl)
          gpc_cart(3) = xk(3, ik) + g(3, igl)
          g2 = gpc_cart(1)**2 + gpc_cart(2)**2 + gpc_cart(3)**2
          pmag = DSQRT(g2) * tpiba
          pcryst = gpc_cart
          CALL cryst_to_cart(1, pcryst, at, -1)
          emd_pcryst(:, ig, ik) = pcryst
          DO ispn = 1, nch
             sm = 0.0_DP
             DO ib = ib_start, ib_end
                weight = rdmft_n(ib, ik)
                IF (weight < 1.0e-30_DP) CYCLE
                contrib = 0.0_DP
                IF (noncolin) THEN
                   contrib = emd_band_weight(ik, ib, ispn, weight, su2, tsqaz, &
                        evc(ig, ib), evc(ig + npwx, ib), ig)
                ELSE IF (lsda .AND. nspin == 2) THEN
                   IF (.NOT. rdmft_emd_spin_sum) THEN
                      IF (isk(ik) /= ispn) CYCLE
                   ENDIF
                   contrib = weight * ABS(evc(ig, ib))**2
                ELSE
                   contrib = weight * ABS(evc(ig, ib))**2
                END IF
                sm = sm + contrib
                band_ch(ib, ispn) = contrib
             END DO
             IF (cache_band) THEN
                DO ib = ib_start, ib_end
                   emd_val_band(ig, ik, ispn, ib - ib_start + 1) = band_ch(ib, ispn)
                END DO
             END IF
             emd_loc(ig, ispn) = sm
          END DO
          IF (band_unit_open) THEN
             DO ib = ib_start, ib_end
                IF (ALL(band_ch(ib, :) < 1.0e-30_DP)) CYCLE
                IF (nch == 1) THEN
                   WRITE(iunit_band, '(3I6,4F14.8,ES18.10)') ik, ib, ig, pcryst, pmag, &
                        band_ch(ib, 1)
                ELSE
                   WRITE(iunit_band, '(3I6,4F14.8,2ES18.10)') ik, ib, ig, pcryst, pmag, &
                        band_ch(ib, 1), band_ch(ib, 2)
                ENDIF
             ENDDO
          ENDIF
       END DO
       emd_val(:, ik, :) = emd_loc
       IF (ionode) WRITE(stdout, '(5X,A,I0,A,I0,A)') &
            'RDMFT EMD: k-point ', ik, ' / ', nks
    END DO
    IF (band_unit_open) THEN
       CLOSE(iunit_band)
       IF (ionode) WRITE(stdout, '(5X,A,A)') 'RDMFT band-resolved EMD written to ', &
            TRIM(TRIM(fprefix)//'.emd_band')
    ENDIF
    IF (cache_band) emd_band_cache_valid = .TRUE.
    DEALLOCATE(band_ch)
    DEALLOCATE(emd_loc)
    !
    ! Electron-count sanity check: sum_k w_k sum_G n(p) = N_e (wk not in n(p)).
    ne_check = 0.0_DP
    DO ik = 1, nks
       DO ig = 1, ngk(ik)
          ne_check = ne_check + wk(ik) * emd_val(ig, ik, 1)
          IF (.NOT. rdmft_emd_spin_sum .AND. nch >= 2) &
               ne_check = ne_check + wk(ik) * emd_val(ig, ik, 2)
       END DO
    END DO
    IF (ionode) THEN
       WRITE(stdout, '(5X,A,F12.6)') 'RDMFT EMD: sum_k w_k sum_G n(p) = ', ne_check
       WRITE(stdout, '(5X,A,F12.6)') 'RDMFT EMD: selected-band occupation sum = ', ne_expected
       IF (ABS(ne_check - ne_expected) > 0.05_DP) &
            WRITE(stdout, '(5X,A)') '  WARNING: EMD integral deviates from N_e'
    ENDIF
    !
    IF (ionode) THEN
       fname = TRIM(fprefix) // '.emd'
       OPEN(NEWUNIT=iunit, FILE=TRIM(fname), STATUS='REPLACE', FORM='FORMATTED', &
            IOSTAT=ios)
       IF (ios /= 0) CALL errore('rdmft_compute_emd', 'opening '//TRIM(fname), ABS(ios))
       WRITE(iunit, '(A)') '# RDMFT electron momentum density n(p)'
       WRITE(iunit, '(A,F12.6,A,F12.6,A,F12.6,A)') '# ecutwfc (Ry) = ', emd_ecutwfc, &
            '  |p|_max_pw (Bohr^-1) = ', emd_pmax_pw, '  (from ecutwfc)'
       IF (emd_kf_total_valid .AND. rdmft_emd_pmax_kf_factor > 0.0_DP) THEN
          WRITE(iunit, '(A,F12.6,A,F12.6,A,F12.6,A)') '# |p|_max (Bohr^-1) = ', emd_pmax, &
               '  (', rdmft_emd_pmax_kf_factor, ' * k_F)'
       ELSE
          WRITE(iunit, '(A,F12.6,A)') '# |p|_max (Bohr^-1) = ', emd_pmax, '  (from ecutwfc)'
       END IF
       IF (emd_kf_total_valid) WRITE(iunit, '(A,F12.6,A)') &
            '# k_F (selected-band equivalent, Bohr^-1) = ', emd_kf_total, &
            '  [free-electron mapping]'
       IF (emd_kf_spin_valid(1)) WRITE(iunit, '(A,F12.6)') '# k_F_up (Bohr^-1) = ', emd_kf_spin(1)
       IF (emd_kf_spin_valid(2)) WRITE(iunit, '(A,F12.6)') '# k_F_dn (Bohr^-1) = ', emd_kf_spin(2)
       WRITE(iunit, '(A,I0)') '# spin channels = ', nch
       WRITE(iunit, '(A,I0,A,I0,A,I0)') '# band window = ', ib_start, ':', ib_end, &
            ' of ', nbnd
       WRITE(iunit, '(A)') '# n(p) excludes k-weights w_k; integrate as sum_k w_k sum_G n(p)'
       WRITE(iunit, '(A)') '# ik ig p1 p2 p3 (crystal reciprocal) |p| (Bohr^-1) n(p) [n_up n_dn]'
       DO ik = 1, nkstot
          DO ig = 1, ngk(ik)
             IF (emd_val(ig, ik, 1) < 1.0e-30_DP .AND. &
                  ALL(emd_val(ig, ik, :) < 1.0e-30_DP)) CYCLE
             pcryst = emd_pcryst(:, ig, ik)
             gpc_cart = pcryst
             CALL cryst_to_cart(1, gpc_cart, bg, +1)
             pmag = DSQRT(gpc_cart(1)**2 + gpc_cart(2)**2 + gpc_cart(3)**2) * tpiba
             IF (nch == 1) THEN
                WRITE(iunit, '(2I6,4F14.8,ES18.10)') ik, ig, pcryst, pmag, &
                     emd_val(ig, ik, 1)
             ELSE
                WRITE(iunit, '(2I6,4F14.8,2ES18.10)') ik, ig, pcryst, pmag, &
                     emd_val(ig, ik, 1), emd_val(ig, ik, 2)
             END IF
          END DO
       END DO
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'RDMFT EMD written to ', TRIM(fname)
    ENDIF
    !
    IF (rdmft_emd_write_elk_fmt .AND. ionode) THEN
       CALL rdmft_write_emd_elk(fprefix)
    ENDIF
  END SUBROUTINE rdmft_compute_emd
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_write_emd_elk(fprefix)
    USE io_global, ONLY : ionode
    USE klist,     ONLY : nkstot, ngk, xk
    USE cell_base, ONLY : at
    CHARACTER(LEN=*), INTENT(IN) :: fprefix
    INTEGER :: ik, ig, iunit, ios, nhk
    REAL(DP) :: xk_cryst(3)
    CHARACTER(LEN=256) :: fname
    IF (.NOT. ionode) RETURN
    fname = 'EMD.OUT'
    OPEN(NEWUNIT=iunit, FILE=TRIM(fname), STATUS='REPLACE', FORM='UNFORMATTED', &
         ACCESS='STREAM', IOSTAT=ios)
    IF (ios /= 0) CALL errore('rdmft_write_emd_elk', 'opening EMD.OUT', ABS(ios))
    DO ik = 1, nkstot
       nhk = ngk(ik)
       xk_cryst = xk(:, ik)
       CALL cryst_to_cart(1, xk_cryst, at, -1)
       WRITE(iunit) xk_cryst(1), xk_cryst(2), xk_cryst(3), nhk
       DO ig = 1, nhk
          WRITE(iunit) emd_val(ig, ik, 1)
       END DO
    END DO
    CLOSE(iunit)
    WRITE(6, '(5X,A,A)') 'ELK-format EMD written to ', TRIM(fname)
  END SUBROUTINE rdmft_write_emd_elk
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION emd_lookup_at_ik(ik, target, ispn)
    USE klist, ONLY : ngk
    INTEGER, INTENT(IN) :: ik, ispn
    REAL(DP), INTENT(IN) :: target(3)
    INTEGER :: ig
    emd_lookup_at_ik = 0.0_DP
    IF (ik < 1 .OR. ik > SIZE(emd_val, 2)) RETURN
    DO ig = 1, ngk(ik)
       IF (momentum_cryst_diff(emd_pcryst(:, ig, ik), target) < EMD_EPS) THEN
          IF (ispn >= 1 .AND. ispn <= emd_nspin_ch) THEN
             emd_lookup_at_ik = emd_val(ig, ik, ispn)
          ELSE
             emd_lookup_at_ik = emd_val(ig, ik, 1)
          END IF
          RETURN
       END IF
    END DO
  END FUNCTION emd_lookup_at_ik
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION emd_lookup_at_ik_band(ik, target, ispn, ib)
    USE klist, ONLY : ngk
    INTEGER, INTENT(IN) :: ik, ispn, ib
    REAL(DP), INTENT(IN) :: target(3)
    INTEGER :: ig, ib_loc
    emd_lookup_at_ik_band = 0.0_DP
    IF (.NOT. emd_band_cache_valid) RETURN
    IF (ib < emd_ib_start .OR. ib > emd_ib_end) RETURN
    IF (ik < 1 .OR. ik > SIZE(emd_val_band, 2)) RETURN
    ib_loc = ib - emd_ib_start + 1
    DO ig = 1, ngk(ik)
       IF (momentum_cryst_diff(emd_pcryst(:, ig, ik), target) < EMD_EPS) THEN
          IF (ispn >= 1 .AND. ispn <= emd_nspin_ch) THEN
             emd_lookup_at_ik_band = emd_val_band(ig, ik, ispn, ib_loc)
          ELSE
             emd_lookup_at_ik_band = emd_val_band(ig, ik, 1, ib_loc)
          END IF
          RETURN
       END IF
    END DO
  END FUNCTION emd_lookup_at_ik_band
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_emd_lookup_band(ik, target, ispn, ib)
    USE lsda_mod,  ONLY : nspin, lsda
    USE rdmft_module, ONLY : rdmft_emd_spin_sum
    INTEGER, INTENT(IN) :: ik, ispn, ib
    REAL(DP), INTENT(IN) :: target(3)
    INTEGER :: ik2
    rdmft_emd_lookup_band = emd_lookup_at_ik_band(ik, target, ispn, ib)
    IF (rdmft_emd_spin_sum .AND. lsda .AND. nspin == 2) THEN
       ik2 = emd_spin_partner(ik)
       IF (ik2 > 0) &
            rdmft_emd_lookup_band = rdmft_emd_lookup_band + &
                 emd_lookup_at_ik_band(ik2, target, ispn, ib)
    END IF
  END FUNCTION rdmft_emd_lookup_band
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_emd_rfhkintp_band(vhpl, ispn, ib)
    USE klist,     ONLY : ngk
    USE cell_base, ONLY : tpiba, bg
    USE symm_base, ONLY : nsym, s
    REAL(DP), INTENT(IN) :: vhpl(3)
    INTEGER, INTENT(IN)  :: ispn, ib
    INTEGER :: ivh0(3), ivk0(3), ivku(3), ivkb(3, 0:1, 0:1, 0:1)
    INTEGER :: i, j, k, mp_n, ik, ns
    REAL(DP) :: vpl(3), v1(3), v2(3), fb(0:1, 0:1, 0:1)
    REAL(DP) :: f00, f01, f10, f11, f0, f1, t1, t2, v0(3), xkr(3)
    REAL(DP) :: pcart(3), pmax2
    !
    rdmft_emd_rfhkintp_band = 0.0_DP
    IF (emd_nkr <= 0) RETURN
    pmax2 = emd_pmax**2
    ivh0 = FLOOR(vhpl)
    vpl = vhpl - DBLE(ivh0)
    ivk0 = FLOOR(vpl * DBLE(emd_ngridk))
    DO i = 0, 1
       DO j = 0, 1
          DO k = 0, 1
             ivkb(1, i, j, k) = ivk0(1) + i
             ivkb(2, i, j, k) = ivk0(2) + j
             ivkb(3, i, j, k) = ivk0(3) + k
             ivkb(1, i, j, k) = MOD(ivkb(1, i, j, k), emd_ngridk(1))
             ivkb(2, i, j, k) = MOD(ivkb(2, i, j, k), emd_ngridk(2))
             ivkb(3, i, j, k) = MOD(ivkb(3, i, j, k), emd_ngridk(3))
          END DO
       END DO
    END DO
    DO i = 0, 1
       DO j = 0, 1
          DO k = 0, 1
             fb(i, j, k) = 0.0_DP
             mp_n = emd_ivkik(ivkb(1, i, j, k), ivkb(2, i, j, k), ivkb(3, i, j, k))
             IF (mp_n < 1 .OR. mp_n > emd_nkr) CYCLE
             ik = emd_equiv(mp_n)
             IF (ik <= 0) CYCLE
             ivku = (/ ivk0(1) + i, ivk0(2) + j, ivk0(3) + k /)
             v1 = DBLE(ivh0) + emd_xkg(:, mp_n)
             v1 = v1 + DBLE(ivku - ivkb(:, i, j, k)) / DBLE(emd_ngridk)
             IF (i == 0 .AND. j == 0 .AND. k == 0) v0 = v1
             pcart = v1
             CALL cryst_to_cart(1, pcart, bg, +1)
             pcart = pcart * tpiba
             IF (pcart(1)**2 + pcart(2)**2 + pcart(3)**2 > pmax2) CYCLE
             fb(i, j, k) = rdmft_emd_lookup_band(ik, v1, ispn, ib)
             IF (fb(i, j, k) < 1.0e-30_DP) THEN
                DO ns = 1, nsym
                   xkr(1) = s(1, 1, ns) * v1(1) + s(1, 2, ns) * v1(2) + s(1, 3, ns) * v1(3)
                   xkr(2) = s(2, 1, ns) * v1(1) + s(2, 2, ns) * v1(2) + s(2, 3, ns) * v1(3)
                   xkr(3) = s(3, 1, ns) * v1(1) + s(3, 2, ns) * v1(2) + s(3, 3, ns) * v1(3)
                   fb(i, j, k) = rdmft_emd_lookup_band(ik, xkr, ispn, ib)
                   IF (fb(i, j, k) > 1.0e-30_DP) EXIT
                END DO
             END IF
          END DO
       END DO
    END DO
    t2 = (vhpl(1) - v0(1)) * DBLE(emd_ngridk(1))
    t1 = 1.0_DP - t2
    f00 = fb(0, 0, 0) * t1 + fb(1, 0, 0) * t2
    f01 = fb(0, 0, 1) * t1 + fb(1, 0, 1) * t2
    f10 = fb(0, 1, 0) * t1 + fb(1, 1, 0) * t2
    f11 = fb(0, 1, 1) * t1 + fb(1, 1, 1) * t2
    t2 = (vhpl(2) - v0(2)) * DBLE(emd_ngridk(2))
    t1 = 1.0_DP - t2
    f0 = f00 * t1 + f10 * t2
    f1 = f01 * t1 + f11 * t2
    t2 = (vhpl(3) - v0(3)) * DBLE(emd_ngridk(3))
    t1 = 1.0_DP - t2
    rdmft_emd_rfhkintp_band = f0 * t1 + f1 * t2
  END FUNCTION rdmft_emd_rfhkintp_band
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION emd_compton_at_band(vpl, ispn, ib, vl2, vl3, n, x, wx)
    INTEGER, INTENT(IN) :: ispn, ib, n
    REAL(DP), INTENT(IN) :: vpl(3), vl2(3), vl3(3), x(n), wx(n)
    REAL(DP) :: f1(n), f2(n), vpt(3)
    INTEGER :: i, j
    DO i = 1, n
       DO j = 1, n
          vpt = vpl + x(i) * vl2 + x(j) * vl3
          f1(j) = rdmft_emd_rfhkintp_band(vpt, ispn, ib)
       END DO
       f2(i) = DOT_PRODUCT(wx, f1)
    END DO
    emd_compton_at_band = DOT_PRODUCT(wx, f2)
  END FUNCTION emd_compton_at_band
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION emd_ne_band(ib)
    USE klist, ONLY : nks, wk
    USE rdmft_module, ONLY : rdmft_n
    INTEGER, INTENT(IN) :: ib
    INTEGER :: ik
    emd_ne_band = 0.0_DP
    DO ik = 1, nks
       emd_ne_band = emd_ne_band + wk(ik) * rdmft_n(ib, ik)
    END DO
  END FUNCTION emd_ne_band
  !
  !-----------------------------------------------------------------
  INTEGER FUNCTION emd_spin_partner(ik)
    USE klist,    ONLY : nkstot, xk
    USE lsda_mod, ONLY : isk
    USE cell_base, ONLY : at
    INTEGER, INTENT(IN) :: ik
    INTEGER :: ik2
    REAL(DP) :: xk1(3), xk2(3)
    emd_spin_partner = 0
    xk1 = xk(:, ik)
    CALL cryst_to_cart(1, xk1, at, -1)
    DO ik2 = 1, nkstot
       IF (ik2 == ik) CYCLE
       IF (isk(ik2) == isk(ik)) CYCLE
       xk2 = xk(:, ik2)
       CALL cryst_to_cart(1, xk2, at, -1)
       IF (periodic_cryst_diff(xk2, xk1) < EMD_EPS) THEN
          emd_spin_partner = ik2
          RETURN
       END IF
    END DO
  END FUNCTION emd_spin_partner
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_emd_lookup(ik, target, ispn)
    USE lsda_mod,  ONLY : nspin, lsda
    USE rdmft_module, ONLY : rdmft_emd_spin_sum
    INTEGER, INTENT(IN) :: ik, ispn
    REAL(DP), INTENT(IN) :: target(3)
    INTEGER :: ik2
    rdmft_emd_lookup = emd_lookup_at_ik(ik, target, ispn)
    IF (rdmft_emd_spin_sum .AND. lsda .AND. nspin == 2) THEN
       ik2 = emd_spin_partner(ik)
       IF (ik2 > 0) &
            rdmft_emd_lookup = rdmft_emd_lookup + emd_lookup_at_ik(ik2, target, ispn)
    END IF
  END FUNCTION rdmft_emd_lookup
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_emd_rfhkintp(vhpl, ispn)
    USE klist,     ONLY : ngk
    USE cell_base, ONLY : tpiba, bg
    USE symm_base, ONLY : nsym, s
    REAL(DP), INTENT(IN) :: vhpl(3)
    INTEGER, INTENT(IN)  :: ispn
    INTEGER :: ivh0(3), ivk0(3), ivku(3), ivkb(3, 0:1, 0:1, 0:1)
    INTEGER :: i, j, k, mp_n, ik, ns
    REAL(DP) :: vpl(3), v1(3), v2(3), fb(0:1, 0:1, 0:1)
    REAL(DP) :: f00, f01, f10, f11, f0, f1, t1, t2, v0(3), xkr(3)
    REAL(DP) :: pcart(3), pmax2
    !
    rdmft_emd_rfhkintp = 0.0_DP
    IF (emd_nkr <= 0) RETURN
    pmax2 = emd_pmax**2
    ivh0 = FLOOR(vhpl)
    vpl = vhpl - DBLE(ivh0)
    ivk0 = FLOOR(vpl * DBLE(emd_ngridk))
    DO i = 0, 1
       DO j = 0, 1
          DO k = 0, 1
             ivkb(1, i, j, k) = ivk0(1) + i
             ivkb(2, i, j, k) = ivk0(2) + j
             ivkb(3, i, j, k) = ivk0(3) + k
             ivkb(1, i, j, k) = MOD(ivkb(1, i, j, k), emd_ngridk(1))
             ivkb(2, i, j, k) = MOD(ivkb(2, i, j, k), emd_ngridk(2))
             ivkb(3, i, j, k) = MOD(ivkb(3, i, j, k), emd_ngridk(3))
          END DO
       END DO
    END DO
    DO i = 0, 1
       DO j = 0, 1
          DO k = 0, 1
             fb(i, j, k) = 0.0_DP
             mp_n = emd_ivkik(ivkb(1, i, j, k), ivkb(2, i, j, k), ivkb(3, i, j, k))
             IF (mp_n < 1 .OR. mp_n > emd_nkr) CYCLE
             ik = emd_equiv(mp_n)
             IF (ik <= 0) CYCLE
             ivku = (/ ivk0(1) + i, ivk0(2) + j, ivk0(3) + k /)
             v1 = DBLE(ivh0) + emd_xkg(:, mp_n)
             v1 = v1 + DBLE(ivku - ivkb(:, i, j, k)) / DBLE(emd_ngridk)
             IF (i == 0 .AND. j == 0 .AND. k == 0) v0 = v1
             pcart = v1
             CALL cryst_to_cart(1, pcart, bg, +1)
             pcart = pcart * tpiba
             IF (pcart(1)**2 + pcart(2)**2 + pcart(3)**2 > pmax2) CYCLE
             fb(i, j, k) = rdmft_emd_lookup(ik, v1, ispn)
             IF (fb(i, j, k) < 1.0e-30_DP) THEN
                DO ns = 1, nsym
                   xkr(1) = s(1, 1, ns) * v1(1) + s(1, 2, ns) * v1(2) + s(1, 3, ns) * v1(3)
                   xkr(2) = s(2, 1, ns) * v1(1) + s(2, 2, ns) * v1(2) + s(2, 3, ns) * v1(3)
                   xkr(3) = s(3, 1, ns) * v1(1) + s(3, 2, ns) * v1(2) + s(3, 3, ns) * v1(3)
                   fb(i, j, k) = rdmft_emd_lookup(ik, xkr, ispn)
                   IF (fb(i, j, k) > 1.0e-30_DP) EXIT
                END DO
             END IF
          END DO
       END DO
    END DO
    t2 = (vhpl(1) - v0(1)) * DBLE(emd_ngridk(1))
    t1 = 1.0_DP - t2
    f00 = fb(0, 0, 0) * t1 + fb(1, 0, 0) * t2
    f01 = fb(0, 0, 1) * t1 + fb(1, 0, 1) * t2
    f10 = fb(0, 1, 0) * t1 + fb(1, 1, 0) * t2
    f11 = fb(0, 1, 1) * t1 + fb(1, 1, 1) * t2
    t2 = (vhpl(2) - v0(2)) * DBLE(emd_ngridk(2))
    t1 = 1.0_DP - t2
    f0 = f00 * t1 + f10 * t2
    f1 = f01 * t1 + f11 * t2
    t2 = (vhpl(3) - v0(3)) * DBLE(emd_ngridk(3))
    t1 = 1.0_DP - t2
    rdmft_emd_rfhkintp = f0 * t1 + f1 * t2
  END FUNCTION rdmft_emd_rfhkintp
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_wsplint(n, x, w)
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: x(n)
    REAL(DP), INTENT(OUT) :: w(n)
    INTEGER :: i
    REAL(DP) :: f(9)
    DO i = 1, MIN(n, 9)
       f = 0.0_DP
       f(i) = 1.0_DP
       w(i) = emd_splint(n, x, f)
    END DO
    IF (n <= 9) RETURN
    DO i = 1, 4
       f = 0.0_DP
       f(i) = 1.0_DP
       w(i) = emd_splint(9, x, f)
    END DO
    f = 0.0_DP
    f(5) = 1.0_DP
    DO i = 5, n - 4
       w(i) = emd_splint(9, x(i - 4), f)
    END DO
    DO i = 1, 4
       f = 0.0_DP
       f(i + 5) = 1.0_DP
       w(n - 4 + i) = emd_splint(9, x(n - 8), f)
    END DO
  END SUBROUTINE emd_wsplint
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION emd_splint(n, x, f)
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: x(n), f(n)
    INTEGER :: i
    REAL(DP) :: x0, x1, x2, x3, y0, y1, y2, y3
    REAL(DP) :: t0, t1, t2, t3, t4, t5, t6, t7
    IF (n <= 4) THEN
       emd_splint = emd_polynm(-1, n, x, f, x(n))
       RETURN
    END IF
    x0 = x(1)
    x1 = x(2) - x0; x2 = x(3) - x0; x3 = x(4) - x0
    t4 = x1 - x2; t5 = x1 - x3; t6 = x2 - x3
    y0 = f(1)
    y1 = f(2) - y0; y2 = f(3) - y0; y3 = f(4) - y0
    t1 = x1 * x2 * y3; t2 = x2 * x3 * y1; t3 = x1 * x3
    t0 = 0.5_DP / (t3 * t4 * t5 * t6)
    t3 = t3 * y2
    t7 = t1 * t4 + t2 * t6 - t3 * t5
    emd_splint = x2 * (y0 + t0 * (x1 * (x1 * y1 + x2 * y2 + x3 * y3) + &
         x2 * (0.5_DP * t7 * x2 - (2.0_DP / 3.0_DP) * (y1 + y2 + y3))))
    DO i = 3, n - 3
       x0 = x(i)
       x1 = x(i - 1) - x0; x2 = x(i + 1) - x0; x3 = x(i + 2) - x0
       t4 = x1 - x2; t5 = x1 - x3; t6 = x2 - x3; t3 = x1 * x3
       y0 = f(i)
       y1 = f(i - 1) - y0; y2 = f(i + 1) - y0; y3 = f(i + 2) - y0
       t1 = x1 * x2 * y3; t2 = x2 * x3 * y1
       t0 = 0.5_DP / (t3 * t4 * t5 * t6)
       t3 = t3 * y2
       t7 = t1 * t4 + t2 * t6 - t3 * t5
       emd_splint = emd_splint + x2 * (y0 + t0 * (x1 * (x1 * y1 + x2 * y2 + x3 * y3) + &
            x2 * (0.5_DP * t7 * x2 - (2.0_DP / 3.0_DP) * (y1 + y2 + y3))))
    END DO
    x0 = x(n - 2)
    x1 = x(n - 3) - x0; x2 = x(n - 1) - x0; x3 = x(n) - x0
    t4 = x1 - x2; t5 = x1 - x3; t6 = x2 - x3
    y0 = f(n - 2)
    y1 = f(n - 3) - y0; y2 = f(n - 1) - y0; y3 = f(n) - y0
    t1 = x1 * x2; t2 = x2 * x3 * y1; t3 = x1 * x3 * y2
    t0 = 0.5_DP / (t1 * t4 * t5 * t6)
    t1 = t1 * y3
    t7 = t1 * t4 + t2 * t6 - t3 * t5
    emd_splint = emd_splint + x3 * (y0 + t0 * (x1 * (x1 * y1 + x2 * y2 + x3 * y3) + &
         x3 * (0.5_DP * t7 * x3 - (2.0_DP / 3.0_DP) * (y1 + y2 + y3))))
  END FUNCTION emd_splint
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION emd_polynm(m, n, x, y, t)
    INTEGER, INTENT(IN) :: m, n
    REAL(DP), INTENT(IN) :: x(n), y(n), t
    REAL(DP) :: c(0:64), z
    INTEGER :: i, j
    IF (n < 1) THEN
       emd_polynm = 0.0_DP
       RETURN
    END IF
    c(0:n - 1) = y(1:n)
    DO j = 1, n - 1
       DO i = 0, n - 1 - j
          c(i) = (c(i + 1) - c(i)) / (x(i + j + 1) - x(i + 1))
       END DO
    END DO
    z = t - x(1)
    emd_polynm = c(0)
    DO i = 1, n - 1
       emd_polynm = emd_polynm + c(i) * z
       z = z * (t - x(i + 1))
    END DO
    IF (m >= 0) THEN
       z = x(1)
       emd_polynm = 0.0_DP
       DO i = 1, n
          emd_polynm = emd_polynm + c(i - 1) * (t - z)**i / DBLE(i)
          z = x(i)
       END DO
    END IF
  END FUNCTION emd_polynm
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_file_prefix(fprefix)
    USE io_files, ONLY : prefix
    USE rdmft_module, ONLY : rdmft_emd_prefix
    CHARACTER(LEN=256), INTENT(OUT) :: fprefix
    IF (TRIM(rdmft_emd_prefix) == ' ') THEN
       fprefix = TRIM(prefix) // '.rdmft'
    ELSE
       fprefix = TRIM(rdmft_emd_prefix)
    END IF
  END SUBROUTINE emd_file_prefix
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_write_pmax_header(iunit)
    USE rdmft_module, ONLY : rdmft_emd_pmax_kf_factor
    INTEGER, INTENT(IN) :: iunit
    WRITE(iunit, '(A,F12.6,A,F12.6,A,F12.6,A)') '# ecutwfc (Ry) = ', emd_ecutwfc, &
         '  |p|_max_pw (Bohr^-1) = ', emd_pmax_pw, '  (from ecutwfc)'
    IF (emd_kf_total_valid .AND. rdmft_emd_pmax_kf_factor > 0.0_DP) THEN
       WRITE(iunit, '(A,F12.6,A,F12.6,A,F12.6,A)') '# |p|_max (Bohr^-1) = ', emd_pmax, &
            '  (', rdmft_emd_pmax_kf_factor, ' * k_F)'
    ELSE
       WRITE(iunit, '(A,F12.6,A)') '# |p|_max (Bohr^-1) = ', emd_pmax, '  (from ecutwfc)'
    END IF
    IF (emd_kf_total_valid) WRITE(iunit, '(A,F12.6,A)') &
         '# k_F (selected-band equivalent, Bohr^-1) = ', emd_kf_total, &
         '  [free-electron mapping]'
    IF (emd_kf_spin_valid(1)) WRITE(iunit, '(A,F12.6)') '# k_F_up (Bohr^-1) = ', emd_kf_spin(1)
    IF (emd_kf_spin_valid(2)) WRITE(iunit, '(A,F12.6)') '# k_F_dn (Bohr^-1) = ', emd_kf_spin(2)
  END SUBROUTINE emd_write_pmax_header
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_setup_direction(vl_dir, vl2, vl3, q_scale, r_scale)
    !---------------------------------------------------------------
    !! Unit direction in crystal reciprocal coords and perpendicular
    !! integration frame; q_scale (Bohr^-1) and r_scale (Bohr) per
    !! unit fractional step along the same Miller indices.
    !
    USE cell_base, ONLY : at, bg, tpiba, alat
    USE rdmft_module, ONLY : rdmft_emd_line_start, rdmft_emd_line_end
    REAL(DP), INTENT(OUT) :: vl_dir(3), vl2(3), vl3(3), q_scale, r_scale
    REAL(DP) :: vc1(3), vc2(3), vc3(3), vdiff(3), pcart(3), rcart(3), t1
    INTEGER :: i
    !
    vdiff = rdmft_emd_line_end - rdmft_emd_line_start
    t1 = DSQRT(vdiff(1)**2 + vdiff(2)**2 + vdiff(3)**2)
    IF (t1 < 1.0e-12_DP) CALL errore('emd_setup_direction', 'zero-length direction', 1)
    vl_dir = vdiff / t1
    emd_cp_dir = vl_dir
    !
    pcart = vl_dir
    CALL cryst_to_cart(1, pcart, bg, +1)
    pcart = pcart * tpiba
    q_scale = DSQRT(pcart(1)**2 + pcart(2)**2 + pcart(3)**2)
    IF (q_scale < 1.0e-12_DP) CALL errore('emd_setup_direction', 'zero momentum scale', 1)
    vc1 = pcart / q_scale
    !
    rcart = vl_dir
    CALL cryst_to_cart(1, rcart, at, +1)
    r_scale = alat * DSQRT(rcart(1)**2 + rcart(2)**2 + rcart(3)**2)
    IF (r_scale < 1.0e-12_DP) CALL errore('emd_setup_direction', 'zero real-space scale', 1)
    !
    i = 1
    IF (ABS(vc1(2)) < ABS(vc1(i))) i = 2
    IF (ABS(vc1(3)) < ABS(vc1(i))) i = 3
    vc2 = 0.0_DP
    vc2(i) = 1.0_DP
    t1 = DOT_PRODUCT(vc1, vc2)
    vc2 = vc2 - t1 * vc1
    t1 = DSQRT(vc2(1)**2 + vc2(2)**2 + vc2(3)**2)
    IF (t1 < 1.0e-12_DP) CALL errore('emd_setup_direction', 'degenerate perpendicular frame', 1)
    vc2 = vc2 / t1
    CALL r3cross(vc1, vc2, vc3)
    vl2(1) = vc2(1) * bg(1, 1) + vc2(2) * bg(1, 2) + vc2(3) * bg(1, 3)
    vl2(2) = vc2(1) * bg(2, 1) + vc2(2) * bg(2, 2) + vc2(3) * bg(2, 3)
    vl2(3) = vc2(1) * bg(3, 1) + vc2(2) * bg(3, 2) + vc2(3) * bg(3, 3)
    vl3(1) = vc3(1) * bg(1, 1) + vc3(2) * bg(1, 2) + vc3(3) * bg(1, 3)
    vl3(2) = vc3(1) * bg(2, 1) + vc3(2) * bg(2, 2) + vc3(3) * bg(2, 3)
    vl3(3) = vc3(1) * bg(3, 1) + vc3(2) * bg(3, 2) + vc3(3) * bg(3, 3)
    vl2 = vl2 / tpiba
    vl3 = vl3 / tpiba
  END SUBROUTINE emd_setup_direction
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_setup_perp_grid(pmax, n, x, wx)
    USE cell_base, ONLY : at, tpiba
    USE constants, ONLY : pi
    INTEGER, INTENT(OUT) :: n
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: x(:), wx(:)
    REAL(DP), INTENT(IN) :: pmax
    INTEGER :: i
    n = 2 * MAX(NINT(pmax * DSQRT(at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2) / pi), 8) + 1
    ALLOCATE(x(n), wx(n))
    DO i = 1, n
       x(i) = (2.0_DP * DBLE(i - 1) / DBLE(n - 1) - 1.0_DP) * pmax / tpiba
    END DO
    CALL emd_wsplint(n, x, wx)
  END SUBROUTINE emd_setup_perp_grid
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION emd_compton_at(vpl, ispn, vl2, vl3, n, x, wx)
    INTEGER, INTENT(IN) :: ispn, n
    REAL(DP), INTENT(IN) :: vpl(3), vl2(3), vl3(3), x(n), wx(n)
    REAL(DP) :: f1(n), f2(n), vpt(3)
    INTEGER :: i, j
    DO i = 1, n
       DO j = 1, n
          vpt = vpl + x(i) * vl2 + x(j) * vl3
          f1(j) = rdmft_emd_rfhkintp(vpt, ispn)
       END DO
       f2(i) = DOT_PRODUCT(wx, f1)
    END DO
    emd_compton_at = DOT_PRODUCT(wx, f2)
  END FUNCTION emd_compton_at
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_build_compton_profile(nq, qmin, qmax, nch, force_rebuild)
    !---------------------------------------------------------------
    !! Build J(q) on a uniform q-grid (Bohr^-1) and cache the result.
    !
    USE rdmft_module, ONLY : rdmft_emd_cp_nq, rdmft_emd_cp_qmax, rdmft_emd_cp_qmin
    INTEGER, INTENT(IN) :: nq, nch
    REAL(DP), INTENT(IN) :: qmin, qmax
    LOGICAL, INTENT(IN) :: force_rebuild
    REAL(DP) :: vl_dir(3), vl2(3), vl3(3), q_scale, r_scale, pmax
    REAL(DP), ALLOCATABLE :: x(:), wx(:)
    REAL(DP) :: vpl(3)
    INTEGER :: iq, ispn, n, nq_use
    REAL(DP) :: dq
    LOGICAL :: need_build
    !
    IF (.NOT. ALLOCATED(emd_val)) &
         CALL errore('emd_build_compton_profile', 'call rdmft_compute_emd first', 1)
    nq_use = MAX(nq, 64)
    IF (MOD(nq_use, 2) /= 0) nq_use = nq_use + 1
    need_build = force_rebuild .OR. .NOT. emd_cp_valid .OR. emd_cp_nq /= nq_use &
         .OR. emd_cp_nch /= nch .OR. ABS(emd_cp_qmin - qmin) > 1.0e-10_DP &
         .OR. ABS(emd_cp_qmax - qmax) > 1.0e-10_DP
    IF (.NOT. need_build) RETURN
    !
    CALL emd_free_compton_cache()
    pmax = emd_pmax
    CALL emd_setup_direction(vl_dir, vl2, vl3, q_scale, r_scale)
    CALL emd_setup_perp_grid(pmax, n, x, wx)
    emd_cp_nq = nq_use
    emd_cp_nch = nch
    emd_cp_qmin = qmin
    emd_cp_qmax = qmax
    emd_cp_qscale = q_scale
    dq = (qmax - qmin) / DBLE(nq_use - 1)
    emd_cp_dq = dq
    ALLOCATE(emd_cp_q(nq_use), emd_cp_j(nq_use, nch))
    DO iq = 1, nq_use
       emd_cp_q(iq) = qmin + DBLE(iq - 1) * dq
    END DO
    DO ispn = 1, nch
       DO iq = 1, nq_use
          vpl = (emd_cp_q(iq) / q_scale) * vl_dir
          emd_cp_j(iq, ispn) = emd_compton_at(vpl, ispn, vl2, vl3, n, x, wx)
       END DO
    END DO
    emd_cp_valid = .TRUE.
    DEALLOCATE(x, wx)
  END SUBROUTINE emd_build_compton_profile
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_resolve_compton_qrange(qmin, qmax)
    !---------------------------------------------------------------
    !! Resolve Compton $q$-grid bounds: default $q \in [0, q_\mathrm{max}]$.
    !
    USE rdmft_module, ONLY : rdmft_emd_cp_qmax, rdmft_emd_cp_qmin
    REAL(DP), INTENT(OUT) :: qmin, qmax
    qmax = rdmft_emd_cp_qmax
    IF (qmax <= 0.0_DP) qmax = emd_pmax
    qmin = MAX(rdmft_emd_cp_qmin, 0.0_DP)
    IF (qmax <= qmin) CALL errore('emd_resolve_compton_qrange', &
         'require cp_qmax > cp_qmin (with cp_qmin >= 0)', 1)
  END SUBROUTINE emd_resolve_compton_qrange
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION emd_compton_sym_factor()
    !! Factor 2 when $J(q)$ is sampled on $q \ge 0$ only (even profile).
    IF (emd_cp_qmin >= -1.0e-10_DP) THEN
       emd_compton_sym_factor = 2.0_DP
    ELSE
       emd_compton_sym_factor = 1.0_DP
    END IF
  END FUNCTION emd_compton_sym_factor
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION emd_integrate_j(ispn)
    INTEGER, INTENT(IN) :: ispn
    INTEGER :: iq
    IF (.NOT. emd_cp_valid) THEN
       emd_integrate_j = 0.0_DP
       RETURN
    END IF
    emd_integrate_j = 0.5_DP * emd_cp_dq * (emd_cp_j(1, ispn) + emd_cp_j(emd_cp_nq, ispn))
    DO iq = 2, emd_cp_nq - 1
       emd_integrate_j = emd_integrate_j + emd_cp_dq * emd_cp_j(iq, ispn)
    END DO
  END FUNCTION emd_integrate_j
  !
  !-----------------------------------------------------------------
  SUBROUTINE emd_dct_autocorr(j_in, nq, dq, npt, r_start, r_end, b_out)
    !---------------------------------------------------------------
    !! B(r) = sum_k J(q_k) cos(q_k r) dq on a uniform r grid.
    !! When the $q$-grid starts at 0, include the $q<0$ mirror ($J$ even).
    !
    INTEGER, INTENT(IN) :: nq, npt
    REAL(DP), INTENT(IN) :: j_in(nq), dq, r_start, r_end
    REAL(DP), INTENT(OUT) :: b_out(npt)
    INTEGER :: ir, iq
    REAL(DP) :: r, arg, symfac
    symfac = emd_compton_sym_factor()
    IF (npt < 1) RETURN
    DO ir = 1, npt
       IF (npt == 1) THEN
          r = r_start
       ELSE
          r = r_start + DBLE(ir - 1) / DBLE(npt - 1) * (r_end - r_start)
       END IF
       b_out(ir) = 0.0_DP
       DO iq = 1, nq
          arg = emd_cp_q(iq) * r
          b_out(ir) = b_out(ir) + j_in(iq) * DCOS(arg)
       END DO
       b_out(ir) = b_out(ir) * dq * symfac
    END DO
  END SUBROUTINE emd_dct_autocorr
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_emd_compton()
    !---------------------------------------------------------------
    USE io_global, ONLY : stdout, ionode
    USE rdmft_module, ONLY : rdmft_emd_spin_sum, rdmft_emd_write_elk_fmt, &
         rdmft_emd_line_start, rdmft_emd_line_end, rdmft_emd_cp_nq, &
         rdmft_emd_do_compton
    CHARACTER(LEN=256) :: fprefix, fname
    INTEGER :: iq, ispn, nch, iunit, ios, iunit_elk, ib
    INTEGER :: iunit_band, n_perp, ib_loc
    REAL(DP) :: qmax, qmin, j_int, j_band, symfac
    REAL(DP) :: vl_dir(3), vl2(3), vl3(3), q_scale, r_scale, pmax
    REAL(DP), ALLOCATABLE :: x(:), wx(:)
    REAL(DP) :: vpl(3)
    !
    IF (.NOT. ALLOCATED(emd_val)) &
         CALL errore('rdmft_emd_compton', 'call rdmft_compute_emd first', 1)
    IF (rdmft_emd_do_compton .AND. .NOT. emd_band_cache_valid) &
         CALL errore('rdmft_emd_compton', &
              'band-resolved Compton requires band cache from rdmft_compute_emd', 1)
    nch = 1
    IF (.NOT. rdmft_emd_spin_sum) nch = emd_nspin_ch
    CALL emd_resolve_compton_qrange(qmin, qmax)
    !
    CALL emd_build_compton_profile(rdmft_emd_cp_nq, qmin, qmax, nch, .TRUE.)
    symfac = emd_compton_sym_factor()
    CALL emd_file_prefix(fprefix)
    !
    IF (ionode) THEN
       fname = TRIM(fprefix) // '.emd_compton'
       OPEN(NEWUNIT=iunit, FILE=TRIM(fname), STATUS='REPLACE', FORM='FORMATTED', &
            IOSTAT=ios)
       IF (ios /= 0) CALL errore('rdmft_emd_compton', 'opening '//TRIM(fname), ABS(ios))
       WRITE(iunit, '(A)') '# directional Compton profile J(q)'
       WRITE(iunit, '(A,3F12.6)') '# line_start (crystal reciprocal) = ', rdmft_emd_line_start
       WRITE(iunit, '(A,3F12.6)') '# line_end   (crystal reciprocal) = ', rdmft_emd_line_end
       WRITE(iunit, '(A,I0,A,I0,A)') '# band window = ', emd_ib_start, ':', emd_ib_end
       WRITE(iunit, '(A)') '# q (Bohr^-1)  J(q)  [spin]  (q >= 0)'
       CALL emd_write_pmax_header(iunit)
       DO ispn = 1, nch
          j_int = emd_integrate_j(ispn) * symfac
          WRITE(iunit, '(A,I0,A,F12.6,A,F12.6,A)') '# spin ', ispn, &
               '  integral(J) dq = ', j_int, '  (x2 mirror for q>=0)  N_e = ', emd_ne_expected
       END DO
       DO ispn = 1, nch
          DO iq = 1, emd_cp_nq
             IF (nch == 1) THEN
                WRITE(iunit, '(2ES18.10)') emd_cp_q(iq), emd_cp_j(iq, ispn)
             ELSE
                WRITE(iunit, '(3ES18.10)') emd_cp_q(iq), emd_cp_j(iq, ispn), DBLE(ispn)
             END IF
             IF (MOD(iq, MAX(emd_cp_nq / 10, 1)) == 0 .OR. iq == emd_cp_nq) &
                  WRITE(stdout, '(5X,A,I0,A,I0,A,I0)') &
                  'RDMFT Compton: spin ', ispn, ' q-point ', iq, ' / ', emd_cp_nq
          END DO
       END DO
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'RDMFT Compton profile written to ', TRIM(fname)
    END IF
    !
    IF (rdmft_emd_do_compton .AND. ionode) THEN
       IF (emd_ib_end < emd_ib_start) CALL errore('rdmft_emd_compton', &
            'empty band window for band-resolved Compton', 1)
       pmax = emd_pmax
       CALL emd_setup_direction(vl_dir, vl2, vl3, q_scale, r_scale)
       CALL emd_setup_perp_grid(pmax, n_perp, x, wx)
       fname = TRIM(fprefix) // '.emd_compton_band'
       OPEN(NEWUNIT=iunit_band, FILE=TRIM(fname), STATUS='REPLACE', FORM='FORMATTED', &
            IOSTAT=ios)
       IF (ios /= 0) CALL errore('rdmft_emd_compton', 'opening '//TRIM(fname), ABS(ios))
       WRITE(iunit_band, '(A)') '# band-resolved directional Compton profile J_ib(q)'
       WRITE(iunit_band, '(A,3F12.6)') '# line_start (crystal reciprocal) = ', rdmft_emd_line_start
       WRITE(iunit_band, '(A,3F12.6)') '# line_end   (crystal reciprocal) = ', rdmft_emd_line_end
       WRITE(iunit_band, '(A,I0,A,I0,A)') '# band window = ', emd_ib_start, ':', emd_ib_end
       WRITE(iunit_band, '(A)') '# ib  q (Bohr^-1)  J_ib(q)  [spin]  (q >= 0)'
       CALL emd_write_pmax_header(iunit_band)
       DO ib = emd_ib_start, emd_ib_end
          ib_loc = ib - emd_ib_start + 1
          WRITE(iunit_band, '(A,I0,A,F12.6)') '# band ', ib, '  n_ib = ', emd_ne_band(ib)
          DO ispn = 1, nch
             DO iq = 1, emd_cp_nq
                vpl = (emd_cp_q(iq) / q_scale) * vl_dir
                j_band = emd_compton_at_band(vpl, ispn, ib, vl2, vl3, n_perp, x, wx)
                IF (nch == 1) THEN
                   WRITE(iunit_band, '(I6,2ES18.10)') ib, emd_cp_q(iq), j_band
                ELSE
                   WRITE(iunit_band, '(I6,3ES18.10)') ib, emd_cp_q(iq), j_band, DBLE(ispn)
                END IF
             END DO
             IF (MOD(ib_loc, MAX((emd_ib_end - emd_ib_start + 1) / 5, 1)) == 0 .OR. &
                  ib == emd_ib_end) &
                  WRITE(stdout, '(5X,A,I0,A,I0,A,I0)') &
                  'RDMFT Compton (band): band ', ib, ' spin ', ispn, ' / ', nch
          END DO
       END DO
       CLOSE(iunit_band)
       DEALLOCATE(x, wx)
       WRITE(stdout, '(5X,A,A)') 'RDMFT band-resolved Compton written to ', TRIM(fname)
    END IF
    !
    IF (rdmft_emd_write_elk_fmt .AND. ionode .AND. nch == 1) THEN
       OPEN(NEWUNIT=iunit_elk, FILE='EMDCOMPTON.OUT', STATUS='REPLACE', FORM='FORMATTED', &
            IOSTAT=ios)
       IF (ios == 0) THEN
          DO iq = 1, emd_cp_nq
             WRITE(iunit_elk, '(2ES18.10)') emd_cp_q(iq), emd_cp_j(iq, 1)
          END DO
          CLOSE(iunit_elk)
          WRITE(stdout, '(5X,A)') 'ELK-format Compton profile written to EMDCOMPTON.OUT'
       END IF
    END IF
  END SUBROUTINE rdmft_emd_compton
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_emd_autocorr()
    !---------------------------------------------------------------
    USE io_global, ONLY : stdout, ionode
    USE rdmft_module, ONLY : rdmft_emd_spin_sum, rdmft_emd_line_start, &
         rdmft_emd_line_end, rdmft_emd_cp_nq, rdmft_emd_af_npt, &
         rdmft_emd_af_r_start, rdmft_emd_af_r_end, rdmft_emd_do_compton
    CHARACTER(LEN=256) :: fprefix, fname
    INTEGER :: ir, ispn, nch, npt, iunit, ios
    REAL(DP) :: qmax, qmin, r_val
    REAL(DP), ALLOCATABLE :: b_prof(:)
    !
    IF (.NOT. ALLOCATED(emd_val)) &
         CALL errore('rdmft_emd_autocorr', 'call rdmft_compute_emd first', 1)
    nch = 1
    IF (.NOT. rdmft_emd_spin_sum) nch = emd_nspin_ch
    CALL emd_resolve_compton_qrange(qmin, qmax)
    npt = MAX(rdmft_emd_af_npt, 2)
    !
    CALL emd_build_compton_profile(rdmft_emd_cp_nq, qmin, qmax, nch, &
         .NOT. rdmft_emd_do_compton)
    CALL emd_file_prefix(fprefix)
    ALLOCATE(b_prof(npt))
    !
    IF (ionode) THEN
       fname = TRIM(fprefix) // '.emd_af'
       OPEN(NEWUNIT=iunit, FILE=TRIM(fname), STATUS='REPLACE', FORM='FORMATTED', &
            IOSTAT=ios)
       IF (ios /= 0) CALL errore('rdmft_emd_autocorr', 'opening '//TRIM(fname), ABS(ios))
       WRITE(iunit, '(A)') '# directional autocorrelation function B(r)'
       WRITE(iunit, '(A,3F12.6)') '# line_start (crystal reciprocal) = ', rdmft_emd_line_start
       WRITE(iunit, '(A,3F12.6)') '# line_end   (crystal reciprocal) = ', rdmft_emd_line_end
       WRITE(iunit, '(A)') '# r (Bohr)  B(r)  [spin]'
       WRITE(iunit, '(A,F12.6,A,F12.6,A,F12.6)') '# q range (Bohr^-1) = ', &
            emd_cp_qmin, ' : ', emd_cp_qmax, '  dq = ', emd_cp_dq
       CALL emd_write_pmax_header(iunit)
    END IF
    !
    DO ispn = 1, nch
       CALL emd_dct_autocorr(emd_cp_j(:, ispn), emd_cp_nq, emd_cp_dq, npt, &
            rdmft_emd_af_r_start, rdmft_emd_af_r_end, b_prof)
       IF (ionode) THEN
          DO ir = 1, npt
             IF (npt == 1) THEN
                r_val = rdmft_emd_af_r_start
             ELSE
                r_val = rdmft_emd_af_r_start + DBLE(ir - 1) / DBLE(npt - 1) &
                     * (rdmft_emd_af_r_end - rdmft_emd_af_r_start)
             END IF
             IF (nch == 1) THEN
                WRITE(iunit, '(2ES18.10)') r_val, b_prof(ir)
             ELSE
                WRITE(iunit, '(3ES18.10)') r_val, b_prof(ir), DBLE(ispn)
             END IF
             IF (MOD(ir, MAX(npt / 10, 1)) == 0 .OR. ir == npt) &
                  WRITE(stdout, '(5X,A,I0,A,I0,A,I0)') &
                  'RDMFT AF: spin ', ispn, ' r-point ', ir, ' / ', npt
          END DO
       END IF
    END DO
    IF (ionode) THEN
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'RDMFT autocorrelation written to ', TRIM(fname)
       WRITE(stdout, '(5X,A,F12.6,A,F12.6)') 'RDMFT AF: B(0) = ', b_prof(1), &
            '  N_e (expected) = ', emd_ne_expected
    END IF
    DEALLOCATE(b_prof)
  END SUBROUTINE rdmft_emd_autocorr
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_emd_plot1d()
    !---------------------------------------------------------------
    !! Legacy Compton profile along a parametric line (ELK emdplot1d).
    !
    USE io_global,   ONLY : stdout, ionode
    USE cell_base,   ONLY : at, bg, tpiba
    USE constants,   ONLY : pi
    USE rdmft_module,   ONLY : rdmft_emd_line_npt, rdmft_emd_line_start, &
                            rdmft_emd_line_end, rdmft_emd_spin_sum, rdmft_emd_write_elk_fmt
    CHARACTER(LEN=256) :: fprefix, fname
    INTEGER :: ip, i, n, npt, iunit, iunit_elk, ios, ispn, nch
    REAL(DP) :: vl_dir(3), vl2(3), vl3(3), q_scale, r_scale, dist, pmax, t1
    REAL(DP) :: vpl(3)
    REAL(DP), ALLOCATABLE :: x(:), wx(:), line_dist(:), profile(:)
    !
    IF (.NOT. ALLOCATED(emd_val)) &
         CALL errore('rdmft_emd_plot1d', 'call rdmft_compute_emd first', 1)
    pmax = emd_pmax
    npt = MAX(rdmft_emd_line_npt, 2)
    nch = 1
    IF (.NOT. rdmft_emd_spin_sum) nch = emd_nspin_ch
    CALL emd_setup_direction(vl_dir, vl2, vl3, q_scale, r_scale)
    CALL emd_setup_perp_grid(pmax, n, x, wx)
    ALLOCATE(line_dist(npt), profile(npt))
    DO ip = 1, npt
       line_dist(ip) = DBLE(ip - 1) / DBLE(npt - 1)
    END DO
    CALL emd_file_prefix(fprefix)
    IF (ionode) THEN
       fname = TRIM(fprefix) // '.emd1d'
       OPEN(NEWUNIT=iunit, FILE=TRIM(fname), STATUS='REPLACE', FORM='FORMATTED', &
            IOSTAT=ios)
       IF (ios /= 0) CALL errore('rdmft_emd_plot1d', 'opening '//TRIM(fname), ABS(ios))
       WRITE(iunit, '(A)') '# twice-integrated EMD (Compton profile, legacy parametric line)'
       CALL emd_write_pmax_header(iunit)
    ENDIF
    !
    DO ispn = 1, nch
       DO ip = 1, npt
          vpl = rdmft_emd_line_start + line_dist(ip) * &
               (rdmft_emd_line_end - rdmft_emd_line_start)
          t1 = emd_compton_at(vpl, ispn, vl2, vl3, n, x, wx)
          dist = line_dist(ip) * DSQRT(SUM((rdmft_emd_line_end - rdmft_emd_line_start)**2))
          profile(ip) = t1
          IF (ionode) THEN
             IF (nch == 1) THEN
                WRITE(iunit, '(2ES18.10)') dist, t1
             ELSE
                WRITE(iunit, '(3ES18.10)') dist, t1, DBLE(ispn)
             END IF
             IF (MOD(ip, MAX(npt / 10, 1)) == 0 .OR. ip == npt) &
                  WRITE(stdout, '(5X,A,I0,A,I0,A,I0)') &
                  'RDMFT EMD 1D: spin ', ispn, ' point ', ip, ' / ', npt
          END IF
       END DO
    END DO
    IF (ionode) CLOSE(iunit)
    IF (rdmft_emd_write_elk_fmt .AND. ionode .AND. rdmft_emd_spin_sum) THEN
       OPEN(NEWUNIT=iunit_elk, FILE='EMD1D.OUT', STATUS='REPLACE', FORM='FORMATTED', &
            IOSTAT=ios)
       IF (ios == 0) THEN
          DO ip = 1, npt
             dist = line_dist(ip) * DSQRT(SUM((rdmft_emd_line_end - rdmft_emd_line_start)**2))
             WRITE(iunit_elk, '(2ES18.10)') dist, profile(ip)
          END DO
          CLOSE(iunit_elk)
          WRITE(stdout, '(5X,A)') 'ELK-format Compton profile written to EMD1D.OUT'
       END IF
    END IF
    IF (ionode) &
         WRITE(stdout, '(5X,A,A)') 'RDMFT Compton profile written to ', TRIM(TRIM(fprefix)//'.emd1d')
    DEALLOCATE(x, wx, line_dist, profile)
  END SUBROUTINE rdmft_emd_plot1d
  !
END MODULE rdmft_emd
