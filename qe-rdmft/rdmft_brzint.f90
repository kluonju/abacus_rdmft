!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_brzint_mod
  !------------------------------------------------------------------
  !! ELK-style Brillouin-zone integration for the RDMFT DOS.
  !!
  !! Ports the histogram integration of ELK's \texttt{brzint.f90} and
  !! the optional \texttt{fsmooth.f90} post-smoothing, plus a
  !! Monkhorst-Pack grid adapter that unfolds symmetry-reduced QE
  !! k-meshes onto the full uniform grid (same idea as
  !! \texttt{tetra\_init} in \texttt{PW/src/tetra.f90}).
  !
  USE kinds,              ONLY : DP
  USE start_k,            ONLY : nk1, nk2, nk3, k1, k2, k3
  USE klist,              ONLY : nkstot
  USE cell_base,          ONLY : at
  USE symm_base,          ONLY : nsym, s, time_reversal, t_rev
  USE noncollin_module,   ONLY : colin_mag
  USE lsda_mod,           ONLY : nspin
  USE rdmft_module,       ONLY : rdmft_dos_ngrkf
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_fsmooth, rdmft_brzint, rdmft_brzint_setup_mp, &
            rdmft_brzint_vkoff
  !
  REAL(DP), SAVE :: brzint_vkoff(3) = 0.0_DP
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_brzint_vkoff(vkoff)
    !---------------------------------------------------------------
    !! Crystal fractional Monkhorst--Pack shift used by the last
    !! ``rdmft_brzint_setup_mp`` call (``vkoff`` = ``k/(2*n)`` for
    !! integer QE offsets, or inferred from the saved IBZ k-list).
    REAL(DP), INTENT(OUT) :: vkoff(3)
    vkoff = brzint_vkoff
  END SUBROUTINE rdmft_brzint_vkoff
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_brzint_infer_vkoff(xk_cryst, ispin, isk_g, vkoff)
    !---------------------------------------------------------------
    !! Infer a uniform MP shift from the symmetry-reduced k-points.
    !! PWscf xml often stores ``k1=k2=k3=0`` even when the run used a
    !! fractional vkoff (e.g.\ ``0.25 0.5 0.625`` in ``K_POINTS``).
    USE io_global, ONLY : stdout, ionode
    REAL(DP), INTENT(IN)    :: xk_cryst(3, nkstot)
    INTEGER,  INTENT(IN)    :: ispin, isk_g(:)
    REAL(DP), INTENT(INOUT) :: vkoff(3)
    REAL(DP) :: x, vk_try, vk_sum(3), vk_seed(3)
    INTEGER :: ik, idim, nv, ncount(3)
    !
    vk_seed = vkoff
    vk_sum = 0.0_DP
    ncount = 0
    DO ik = 1, nkstot
       IF (nspin == 2 .AND. isk_g(ik) /= ispin) CYCLE
       DO idim = 1, 3
          IF (idim == 1) nv = nk1
          IF (idim == 2) nv = nk2
          IF (idim == 3) nv = nk3
          x = xk_cryst(idim, ik) - NINT(xk_cryst(idim, ik))
          vk_try = x - DBLE(NINT(x * DBLE(nv))) / DBLE(nv)
          vk_try = vk_try - NINT(vk_try)
          vk_sum(idim) = vk_sum(idim) + vk_try
          ncount(idim) = ncount(idim) + 1
       ENDDO
    ENDDO
    DO idim = 1, 3
       IF (ncount(idim) == 0) CYCLE
       vk_try = vk_sum(idim) / DBLE(ncount(idim))
       vk_try = vk_try - NINT(vk_try)
       IF (ABS(vk_seed(idim)) <= 1.0e-8_DP .OR. &
            ABS(vk_try - vk_seed(idim)) > 1.0e-4_DP) vkoff(idim) = vk_try
    ENDDO
    IF (ionode .AND. ANY(ABS(vkoff) > 1.0e-6_DP)) &
         WRITE(stdout, '(5X,A,3F8.4)') &
         'RDMFT brzint: MP vkoff (crystal) = ', vkoff(1), vkoff(2), vkoff(3)
  END SUBROUTINE rdmft_brzint_infer_vkoff
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_fsmooth(m, n, f)
    !---------------------------------------------------------------
    !! $m$ successive 3-point running averages (ELK \texttt{fsmooth}).
    INTEGER, INTENT(IN)    :: m, n
    REAL(DP), INTENT(INOUT) :: f(n)
    INTEGER :: i, j
    REAL(DP) :: f1, f2, f3
    DO i = 1, m
       f1 = f(1)
       f2 = f(2)
       DO j = 2, n - 1
          f3 = f(j + 1)
          f(j) = (f1 + f2 + f3) / 3.0_DP
          f1 = f2
          f2 = f3
       ENDDO
    ENDDO
  END SUBROUTINE rdmft_fsmooth
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_brzint(nsm, ngridk, nsk, ivkik, nw, wint, n, ld, e, f, g)
    !---------------------------------------------------------------
    !! Port of ELK \texttt{brzint}: trilinear k-interpolation plus
    !! nearest-bin histogram integration.
    INTEGER, INTENT(IN)  :: nsm, ngridk(3), nsk(3)
    INTEGER, INTENT(IN)  :: ivkik(0:ngridk(1)-1, 0:ngridk(2)-1, 0:ngridk(3)-1)
    INTEGER, INTENT(IN)  :: nw, n, ld
    REAL(DP), INTENT(IN) :: wint(2), e(ld, *), f(ld, *)
    REAL(DP), INTENT(OUT) :: g(nw)
    INTEGER :: nk, i1, i2, i3, j1, j2, j3, k1, k2, k3, i, iw
    INTEGER :: i000, i001, i010, i011, i100, i101, i110, i111
    REAL(DP) :: wd, dw, dwi, w1, t1, t2, norm
    REAL(DP) :: f0(n), f1(n), e0(n), e1(n)
    REAL(DP) :: f00(n), f01(n), f10(n), f11(n)
    REAL(DP) :: e00(n), e01(n), e10(n), e11(n)
    !
    IF (ngridk(1) < 1 .OR. ngridk(2) < 1 .OR. ngridk(3) < 1) &
         CALL errore('rdmft_brzint', 'ngridk < 1', 1)
    IF (nsk(1) < 1 .OR. nsk(2) < 1 .OR. nsk(3) < 1) &
         CALL errore('rdmft_brzint', 'nsk < 1', 1)
    !
    nk = ngridk(1) * ngridk(2) * ngridk(3)
    wd = wint(2) - wint(1)
    dw = wd / DBLE(nw)
    dwi = 1.0_DP / dw
    w1 = wint(1) * dwi
    g = 0.0_DP
    !
    DO j1 = 0, ngridk(1) - 1
       DO j2 = 0, ngridk(2) - 1
          DO j3 = 0, ngridk(3) - 1
             k1 = MOD(j1 + 1, ngridk(1))
             k2 = MOD(j2 + 1, ngridk(2))
             k3 = MOD(j3 + 1, ngridk(3))
             i000 = ivkik(j1, j2, j3)
             i001 = ivkik(j1, j2, k3)
             i010 = ivkik(j1, k2, j3)
             i011 = ivkik(j1, k2, k3)
             i100 = ivkik(k1, j2, j3)
             i101 = ivkik(k1, j2, k3)
             i110 = ivkik(k1, k2, j3)
             i111 = ivkik(k1, k2, k3)
             DO i1 = 0, nsk(1) - 1
                t2 = DBLE(i1) / DBLE(nsk(1))
                t1 = 1.0_DP - t2
                f00(1:n) = f(1:n, i000) * t1 + f(1:n, i100) * t2
                f01(1:n) = f(1:n, i001) * t1 + f(1:n, i101) * t2
                f10(1:n) = f(1:n, i010) * t1 + f(1:n, i110) * t2
                f11(1:n) = f(1:n, i011) * t1 + f(1:n, i111) * t2
                t1 = t1 * dwi
                t2 = t2 * dwi
                e00(1:n) = e(1:n, i000) * t1 + e(1:n, i100) * t2 - w1
                e01(1:n) = e(1:n, i001) * t1 + e(1:n, i101) * t2 - w1
                e10(1:n) = e(1:n, i010) * t1 + e(1:n, i110) * t2 - w1
                e11(1:n) = e(1:n, i011) * t1 + e(1:n, i111) * t2 - w1
                DO i2 = 0, nsk(2) - 1
                   t2 = DBLE(i2) / DBLE(nsk(2))
                   t1 = 1.0_DP - t2
                   f0(1:n) = f00(1:n) * t1 + f10(1:n) * t2
                   f1(1:n) = f01(1:n) * t1 + f11(1:n) * t2
                   e0(1:n) = e00(1:n) * t1 + e10(1:n) * t2
                   e1(1:n) = e01(1:n) * t1 + e11(1:n) * t2
                   DO i3 = 0, nsk(3) - 1
                      t2 = DBLE(i3) / DBLE(nsk(3))
                      t1 = 1.0_DP - t2
                      DO i = 1, n
                         iw = NINT(e0(i) * t1 + e1(i) * t2) + 1
                         IF (iw >= 1 .AND. iw <= nw) &
                              g(iw) = g(iw) + f0(i) * t1 + f1(i) * t2
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    norm = dw * DBLE(nk) * DBLE(nsk(1) * nsk(2) * nsk(3))
    g = g / norm
    IF (nsm > 0) CALL rdmft_fsmooth(nsm, nw, g)
    !
  END SUBROUTINE rdmft_brzint
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_brzint_mp_index(xkr, on_grid, i, j, k)
    !---------------------------------------------------------------
    !! Locate a crystal-coordinate k-point on the saved Monkhorst-Pack
    !! grid (same convention as ``kpoint_grid`` / ``rdmft_krefine_mp_index``).
    REAL(DP), INTENT(IN)  :: xkr(3)
    LOGICAL, INTENT(OUT)  :: on_grid
    INTEGER, INTENT(OUT)  :: i, j, k
    REAL(DP) :: xx, yy, zz
    REAL(DP), PARAMETER :: eps_loc = 1.0e-5_DP
    !
    xx = xkr(1) * DBLE(nk1) - brzint_vkoff(1) * DBLE(nk1)
    yy = xkr(2) * DBLE(nk2) - brzint_vkoff(2) * DBLE(nk2)
    zz = xkr(3) * DBLE(nk3) - brzint_vkoff(3) * DBLE(nk3)
    on_grid = ABS(xx - NINT(xx)) <= eps_loc .AND. &
         ABS(yy - NINT(yy)) <= eps_loc .AND. &
         ABS(zz - NINT(zz)) <= eps_loc
    IF (.NOT. on_grid) RETURN
    i = MOD(NINT(xkr(1) * DBLE(nk1) - brzint_vkoff(1) * DBLE(nk1) + 2 * DBLE(nk1)), nk1) + 1
    j = MOD(NINT(xkr(2) * DBLE(nk2) - brzint_vkoff(2) * DBLE(nk2) + 2 * DBLE(nk2)), nk2) + 1
    k = MOD(NINT(xkr(3) * DBLE(nk3) - brzint_vkoff(3) * DBLE(nk3) + 2 * DBLE(nk3)), nk3) + 1
    !
  END SUBROUTINE rdmft_brzint_mp_index
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_brzint_setup_mp(ngridk, nsk, ivkik, nkr, equiv, ispin, isk_g, &
       xk_full, ierr)
    !---------------------------------------------------------------
    !! Build the uniform Monkhorst-Pack grid and map each grid point
    !! to an index in the caller's \texttt{nkstot} arrays.
    !!
    !! \texttt{equiv(n)} is the matching k-index in the spectral data
    !! (0 if no match for the requested spin channel).
    !!
    !! \texttt{xk_full} must hold the full \texttt{nkstot} k-list
    !! (pool-collected by the caller; do not read \texttt{klist\%xk}
    !! here after k-pool division).
    INTEGER, INTENT(OUT) :: ngridk(3), nsk(3)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ivkik(:,:,:), equiv(:)
    INTEGER, INTENT(OUT) :: nkr, ierr
    INTEGER, INTENT(IN)  :: ispin, isk_g(:)
    REAL(DP), INTENT(IN) :: xk_full(3, nkstot)
    !
    REAL(DP), PARAMETER :: eps = 1.0e-5_DP
    REAL(DP), ALLOCATABLE :: xkg(:,:), xk_cryst(:,:)
    REAL(DP) :: xkr(3), deltap(3), deltam(3)
    INTEGER :: i, j, k, n, nk, ns, ik, mp_n, ii
    LOGICAL :: on_grid
    !
    ierr = 0
    nkr = 0
    ngridk = 0
    nsk = 1
    !
    IF (nk1 <= 0 .OR. nk2 <= 0 .OR. nk3 <= 0) THEN
       ! Gamma-only ({gamma} input): treat as a 1x1x1 Monkhorst-Pack mesh.
       IF (nkstot == 1) THEN
          ngridk = 1
          nkr = 1
          nsk = 1
          ALLOCATE(ivkik(0:0, 0:0, 0:0), equiv(1))
          ivkik(0, 0, 0) = 1
          equiv(1) = 1
          IF (nspin == 2 .AND. isk_g(1) /= ispin) equiv(1) = 0
          IF (equiv(1) == 0) ierr = 2
          RETURN
       ENDIF
       ierr = 1
       RETURN
    ENDIF
    !
    ngridk(1) = nk1
    ngridk(2) = nk2
    ngridk(3) = nk3
    nkr = nk1 * nk2 * nk3
    nsk(1) = MAX(rdmft_dos_ngrkf / nk1, 1)
    nsk(2) = MAX(rdmft_dos_ngrkf / nk2, 1)
    nsk(3) = MAX(rdmft_dos_ngrkf / nk3, 1)
    !
    ALLOCATE(ivkik(0:nk1-1, 0:nk2-1, 0:nk3-1))
    ALLOCATE(equiv(nkr), xkg(3, nkr), xk_cryst(3, nkstot))
    equiv = 0
    !
    xk_cryst = xk_full
    ! Match EXX (``exx_grid_init``): ``klist%xk`` / pool-collected k are in
    ! Cartesian units; Monkhorst-Pack indexing and symmetry ops use crystal
    ! fractional coordinates.
    CALL cryst_to_cart(nkstot, xk_cryst, at, -1)
    brzint_vkoff(1) = DBLE(k1) / (2.0_DP * nk1)
    brzint_vkoff(2) = DBLE(k2) / (2.0_DP * nk2)
    brzint_vkoff(3) = DBLE(k3) / (2.0_DP * nk3)
    CALL rdmft_brzint_infer_vkoff(xk_cryst, ispin, isk_g, brzint_vkoff)
    !
    DO i = 1, nk1
       DO j = 1, nk2
          DO k = 1, nk3
             n = (k - 1) + (j - 1) * nk3 + (i - 1) * nk2 * nk3 + 1
             xkg(1, n) = DBLE(i - 1) / nk1 + brzint_vkoff(1)
             xkg(2, n) = DBLE(j - 1) / nk2 + brzint_vkoff(2)
             xkg(3, n) = DBLE(k - 1) / nk3 + brzint_vkoff(3)
             ivkik(i - 1, j - 1, k - 1) = n
          ENDDO
       ENDDO
    ENDDO
    !
    ! Primary route: map each saved irreducible k to its MP grid cell via
    ! the ``kpoint_grid`` indexing (``rdmft_krefine_build_mp_to_ik``).
    ! A symmetry-only search with first-match assignment can miss later
    ! IBZ points that fold onto the same star (common for 3x3x3 meshes).
    DO ik = 1, nkstot
       IF (nspin == 2 .AND. isk_g(ik) /= ispin) CYCLE
       xkr = xk_cryst(:, ik)
       DO i = 1, 3
          xkr(i) = xkr(i) - NINT(xkr(i))
       ENDDO
       CALL rdmft_brzint_mp_index(xkr, on_grid, i, j, k)
       IF (.NOT. on_grid) CYCLE
       mp_n = (k - 1) + (j - 1) * nk3 + (i - 1) * nk2 * nk3 + 1
       IF (mp_n >= 1 .AND. mp_n <= nkr) equiv(mp_n) = ik
    ENDDO
    !
    ! Symmetry fallback for any still-unmapped full-grid points.
    DO nk = 1, nkr
       IF (equiv(nk) /= 0) CYCLE
       DO ik = 1, nkstot
          IF (nspin == 2 .AND. isk_g(ik) /= ispin) CYCLE
          DO ns = 1, nsym
             DO i = 1, 3
                xkr(i) = s(i, 1, ns) * xk_cryst(1, ik) + &
                         s(i, 2, ns) * xk_cryst(2, ik) + &
                         s(i, 3, ns) * xk_cryst(3, ik)
             ENDDO
             IF (t_rev(ns) == 1 .AND. colin_mag < 2) xkr = -xkr
             DO ii = 1, 3
                xkr(ii) = xkr(ii) - NINT(xkr(ii))
             ENDDO
             DO i = 1, 3
                deltap(i) = xkr(i) - xkg(i, nk) - NINT(xkr(i) - xkg(i, nk))
                deltam(i) = xkr(i) + xkg(i, nk) - NINT(xkr(i) + xkg(i, nk))
             ENDDO
             IF (DSQRT(deltap(1)**2 + deltap(2)**2 + deltap(3)**2) < eps) THEN
                equiv(nk) = ik
                GOTO 10
             ENDIF
             IF (time_reversal) THEN
                IF (DSQRT(deltam(1)**2 + deltam(2)**2 + deltam(3)**2) < eps) THEN
                   equiv(nk) = ik
                   GOTO 10
                ENDIF
             ENDIF
          ENDDO
       ENDDO
10    CONTINUE
    ENDDO
    !
    DO nk = 1, nkr
       IF (equiv(nk) == 0) ierr = 2
    ENDDO
    !
    DEALLOCATE(xkg, xk_cryst)
    !
  END SUBROUTINE rdmft_brzint_setup_mp
  !
END MODULE rdmft_brzint_mod
