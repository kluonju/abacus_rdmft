!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_stiefel
  !------------------------------------------------------------------
  !! Plane-wave Stiefel manifold operations for the RDMFT solver.
  !!
  !! For PWscf, the natural orbitals are stored as
  !! \texttt{evc(npwx*npol, nbnd)} with the orthonormality constraint
  !!
  !!   evc^H * evc = I
  !!
  !! at the (norm-conserving) plane-wave level.  This is the standard
  !! Stiefel manifold ``St(nbnd, npw; I)``; the AO overlap that the
  !! ABACUS LCAO implementation has to carry around degenerates to
  !! ``S = I`` here, so the projection and retraction are particularly
  !! simple:
  !!
  !!   project_tangent(C, G) = G - C * sym(C^H G)
  !!   retract(C, eta, alpha) = QR factor of (C + alpha*eta)
  !!
  !! For the gamma-only path the wavefunctions are real-coefficient
  !! and we use the gamma-trick conventions of QE.
  !!
  !! Mirror of \texttt{rdmft\_stiefel.h} from the ABACUS reference.
  !
  USE kinds,    ONLY : DP
  USE mp,       ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: stiefel_project_tangent_k, stiefel_retract_k, &
            stiefel_inner_product_k, stiefel_orthonormalize_k, &
            stiefel_project_tangent_gamma, stiefel_retract_gamma, &
            stiefel_inner_product_gamma, stiefel_orthonormalize_gamma
  !
CONTAINS
  !
  !=================================================================
  !  Generic-k (complex) routines
  !=================================================================
  !
  SUBROUTINE stiefel_project_tangent_k(C, G, npw, nbnd, lda, eta)
    !---------------------------------------------------------------
    !! Tangent projection at \(C\):
    !!
    !!   eta = G - C * sym(C^H G)  with sym(A) = (A + A^H)/2.
    !!
    !! Returns a tangent vector satisfying ``C^H eta + eta^H C = 0``.
    INTEGER, INTENT(IN) :: npw, nbnd, lda
    COMPLEX(DP), INTENT(IN)  :: C(lda, nbnd), G(lda, nbnd)
    COMPLEX(DP), INTENT(OUT) :: eta(lda, nbnd)
    !
    COMPLEX(DP), ALLOCATABLE :: A(:,:), Asym(:,:)
    INTEGER :: i, j
    !
    ALLOCATE(A(nbnd, nbnd), Asym(nbnd, nbnd))
    ! A = C^H * G  (npw x nbnd  ->  nbnd x nbnd).
    ! ZGEMM only touches the local G-vector slice, so plane-wave
    ! parallel runs need an explicit reduction across the band
    ! group to recover the full matrix.
    CALL ZGEMM('C', 'N', nbnd, nbnd, npw, &
               (1.0_DP, 0.0_DP), C, lda, G, lda, &
               (0.0_DP, 0.0_DP), A, nbnd)
    CALL mp_sum(A, intra_bgrp_comm)
    ! Symmetrise: Asym = 0.5*(A + A^H).
    DO j = 1, nbnd
       DO i = 1, nbnd
          Asym(i, j) = 0.5_DP * (A(i, j) + CONJG(A(j, i)))
       ENDDO
    ENDDO
    eta(:, :) = G(:, :)
    CALL ZGEMM('N', 'N', npw, nbnd, nbnd, &
               (-1.0_DP, 0.0_DP), C, lda, Asym, nbnd, &
               (1.0_DP, 0.0_DP), eta, lda)
    DEALLOCATE(A, Asym)
    !
  END SUBROUTINE stiefel_project_tangent_k
  !
  !-----------------------------------------------------------------
  SUBROUTINE stiefel_orthonormalize_k(Y, npw, nbnd, lda)
    !---------------------------------------------------------------
    !! Cholesky / QR-style polar reorthonormalisation:
    !!   Y^H Y = L L^H,  Y_new = Y * L^{-H}.
    !! In place.
    INTEGER, INTENT(IN) :: npw, nbnd, lda
    COMPLEX(DP), INTENT(INOUT) :: Y(lda, nbnd)
    !
    COMPLEX(DP), ALLOCATABLE :: M(:,:), tau_vec(:), work(:)
    INTEGER :: info, lwork, i, j
    !
    ALLOCATE(M(nbnd, nbnd))
    ! M = Y^H Y (partial over the local G-vector slice)
    CALL ZGEMM('C', 'N', nbnd, nbnd, npw, &
               (1.0_DP, 0.0_DP), Y, lda, Y, lda, &
               (0.0_DP, 0.0_DP), M, nbnd)
    CALL mp_sum(M, intra_bgrp_comm)
    ! Hermitian symmetrise to fight roundoff.
    DO j = 1, nbnd
       DO i = 1, j - 1
          M(i, j) = 0.5_DP * (M(i, j) + CONJG(M(j, i)))
          M(j, i) = CONJG(M(i, j))
       ENDDO
       M(j, j) = CMPLX(REAL(M(j, j), KIND=DP), 0.0_DP, KIND=DP)
    ENDDO
    ! Cholesky M = L L^H (lower).
    CALL ZPOTRF('L', nbnd, M, nbnd, info)
    IF (info /= 0) THEN
       ! Add a tiny ridge and retry (defensive against roundoff).
       DO j = 1, nbnd
          M(j, j) = M(j, j) + CMPLX(1.0e-10_DP, 0.0_DP, KIND=DP)
       ENDDO
       CALL ZPOTRF('L', nbnd, M, nbnd, info)
       IF (info /= 0) CALL errore('stiefel_orthonormalize_k', &
            'ZPOTRF failed (singular Y^H Y).', info)
    ENDIF
    ! Solve Y * L^{-H}, i.e. find X with X * L^H = Y.
    ! ZTRSM: X*op(A) = alpha*B with X overwriting B.
    CALL ZTRSM('R', 'L', 'C', 'N', npw, nbnd, &
               (1.0_DP, 0.0_DP), M, nbnd, Y, lda)
    DEALLOCATE(M)
    !
  END SUBROUTINE stiefel_orthonormalize_k
  !
  !-----------------------------------------------------------------
  SUBROUTINE stiefel_retract_k(C, eta, alpha, npw, nbnd, lda, C_new)
    !---------------------------------------------------------------
    !! \(R_C(\alpha\eta) = (C + \alpha\eta) L^{-H}\) with
    !! \((C + \alpha\eta)^H (C + \alpha\eta) = L L^H\).
    INTEGER, INTENT(IN) :: npw, nbnd, lda
    COMPLEX(DP), INTENT(IN)  :: C(lda, nbnd), eta(lda, nbnd)
    REAL(DP),    INTENT(IN)  :: alpha
    COMPLEX(DP), INTENT(OUT) :: C_new(lda, nbnd)
    !
    INTEGER :: i, j
    !
    DO j = 1, nbnd
       DO i = 1, npw
          C_new(i, j) = C(i, j) + CMPLX(alpha, 0.0_DP, KIND=DP) * eta(i, j)
       ENDDO
       ! pad zeros above npw if lda > npw
       DO i = npw + 1, lda
          C_new(i, j) = (0.0_DP, 0.0_DP)
       ENDDO
    ENDDO
    CALL stiefel_orthonormalize_k(C_new, npw, nbnd, lda)
    !
  END SUBROUTINE stiefel_retract_k
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION stiefel_inner_product_k(eta1, eta2, npw, nbnd, lda)
    !---------------------------------------------------------------
    !! Frobenius (Euclidean) inner product Re Tr(eta1^H eta2).
    INTEGER, INTENT(IN) :: npw, nbnd, lda
    COMPLEX(DP), INTENT(IN) :: eta1(lda, nbnd), eta2(lda, nbnd)
    !
    INTEGER :: i, j
    REAL(DP) :: s
    COMPLEX(DP) :: c
    !
    s = 0.0_DP
    DO j = 1, nbnd
       DO i = 1, npw
          c = CONJG(eta1(i, j)) * eta2(i, j)
          s = s + REAL(c, KIND=DP)
       ENDDO
    ENDDO
    ! Plane-wave parallel reduction: each rank only owns its slice of
    ! the G-vector axis.
    CALL mp_sum(s, intra_bgrp_comm)
    stiefel_inner_product_k = s
    !
  END FUNCTION stiefel_inner_product_k
  !
  !=================================================================
  !  Gamma-only (real) routines
  !=================================================================
  !  In QE's gamma-only convention, evc(:, ibnd) holds half of the
  !  reciprocal-space coefficients (G and -G are conjugates).  The
  !  inner product is
  !
  !     <psi_i | psi_j> = 2 Re sum_{G != 0} c_i(G)^* c_j(G)
  !                       + c_i(G=0)^* c_j(G=0)
  !
  !  -- the so-called ``gamma trick''.  We expose helper routines that
  !  accept the gamma-trick array directly (passed in as
  !  COMPLEX(DP)(npw, nbnd)) and use ``gstart`` from QE's gvect module
  !  to handle the G=0 component without double-counting.  For the
  !  alternating loop we always go through the complex path above with
  !  the appropriate inner-product subroutines selected by the solver.
  !=================================================================
  !
  SUBROUTINE stiefel_project_tangent_gamma(C, G, npw, nbnd, lda, gstart, eta)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN) :: npw, nbnd, lda, gstart
    COMPLEX(DP), INTENT(IN)  :: C(lda, nbnd), G(lda, nbnd)
    COMPLEX(DP), INTENT(OUT) :: eta(lda, nbnd)
    !
    REAL(DP), ALLOCATABLE :: A(:,:), Asym(:,:)
    INTEGER :: i, j
    !
    ALLOCATE(A(nbnd, nbnd), Asym(nbnd, nbnd))
    ! A = 2 Re(C^H G), with G=0 component counted only once.
    ! The DGEMM here only touches the local G-vector slice; sum the
    ! partial result across the band group before symmetrising.
    CALL DGEMM('T', 'N', nbnd, nbnd, 2*npw, &
               2.0_DP, C, 2*lda, G, 2*lda, 0.0_DP, A, nbnd)
    IF (gstart == 2) THEN
       DO j = 1, nbnd
          DO i = 1, nbnd
             A(i, j) = A(i, j) - REAL(C(1, i), KIND=DP) * REAL(G(1, j), KIND=DP)
          ENDDO
       ENDDO
    ENDIF
    CALL mp_sum(A, intra_bgrp_comm)
    ! Symmetrise (real symmetric).
    DO j = 1, nbnd
       DO i = 1, nbnd
          Asym(i, j) = 0.5_DP * (A(i, j) + A(j, i))
       ENDDO
    ENDDO
    eta(:, :) = G(:, :)
    ! eta = G - C * Asym
    DO j = 1, nbnd
       DO i = 1, nbnd
          IF (Asym(i, j) /= 0.0_DP) &
             CALL ZAXPY(npw, CMPLX(-Asym(i, j), 0.0_DP, KIND=DP), &
                        C(1, i), 1, eta(1, j), 1)
       ENDDO
    ENDDO
    DEALLOCATE(A, Asym)
    !
  END SUBROUTINE stiefel_project_tangent_gamma
  !
  !-----------------------------------------------------------------
  SUBROUTINE stiefel_orthonormalize_gamma(Y, npw, nbnd, lda, gstart)
    !---------------------------------------------------------------
    !! Real (gamma-trick) Cholesky orthonormalisation in place.
    INTEGER, INTENT(IN) :: npw, nbnd, lda, gstart
    COMPLEX(DP), INTENT(INOUT) :: Y(lda, nbnd)
    !
    REAL(DP), ALLOCATABLE :: M(:,:)
    INTEGER :: info, i, j
    !
    ALLOCATE(M(nbnd, nbnd))
    ! Partial Gram matrix on the local G-vector slice.
    CALL DGEMM('T', 'N', nbnd, nbnd, 2*npw, &
               2.0_DP, Y, 2*lda, Y, 2*lda, 0.0_DP, M, nbnd)
    IF (gstart == 2) THEN
       DO j = 1, nbnd
          DO i = 1, nbnd
             M(i, j) = M(i, j) - REAL(Y(1, i), KIND=DP) * REAL(Y(1, j), KIND=DP)
          ENDDO
       ENDDO
    ENDIF
    ! Sum the partials over the band group before Cholesky.
    CALL mp_sum(M, intra_bgrp_comm)
    DO j = 1, nbnd
       DO i = 1, j - 1
          M(i, j) = 0.5_DP * (M(i, j) + M(j, i))
          M(j, i) = M(i, j)
       ENDDO
    ENDDO
    CALL DPOTRF('L', nbnd, M, nbnd, info)
    IF (info /= 0) THEN
       DO j = 1, nbnd
          M(j, j) = M(j, j) + 1.0e-10_DP
       ENDDO
       CALL DPOTRF('L', nbnd, M, nbnd, info)
       IF (info /= 0) CALL errore('stiefel_orthonormalize_gamma', &
            'DPOTRF failed.', info)
    ENDIF
    ! Y = Y * L^{-T}: solve Y * L^T = Y_new -> use ZTRSM treating Y as 2*npw x nbnd real.
    CALL DTRSM('R', 'L', 'T', 'N', 2*npw, nbnd, &
               1.0_DP, M, nbnd, Y, 2*lda)
    DEALLOCATE(M)
    !
  END SUBROUTINE stiefel_orthonormalize_gamma
  !
  !-----------------------------------------------------------------
  SUBROUTINE stiefel_retract_gamma(C, eta, alpha, npw, nbnd, lda, gstart, C_new)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN) :: npw, nbnd, lda, gstart
    COMPLEX(DP), INTENT(IN)  :: C(lda, nbnd), eta(lda, nbnd)
    REAL(DP),    INTENT(IN)  :: alpha
    COMPLEX(DP), INTENT(OUT) :: C_new(lda, nbnd)
    !
    INTEGER :: i, j
    !
    DO j = 1, nbnd
       DO i = 1, npw
          C_new(i, j) = C(i, j) + CMPLX(alpha, 0.0_DP, KIND=DP) * eta(i, j)
       ENDDO
       DO i = npw + 1, lda
          C_new(i, j) = (0.0_DP, 0.0_DP)
       ENDDO
    ENDDO
    CALL stiefel_orthonormalize_gamma(C_new, npw, nbnd, lda, gstart)
    !
  END SUBROUTINE stiefel_retract_gamma
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION stiefel_inner_product_gamma(eta1, eta2, npw, nbnd, lda, gstart)
    !---------------------------------------------------------------
    !! Real inner product with the gamma trick.
    INTEGER, INTENT(IN) :: npw, nbnd, lda, gstart
    COMPLEX(DP), INTENT(IN) :: eta1(lda, nbnd), eta2(lda, nbnd)
    !
    INTEGER :: i, j
    REAL(DP) :: s
    !
    s = 0.0_DP
    DO j = 1, nbnd
       DO i = 1, npw
          s = s + 2.0_DP * REAL(CONJG(eta1(i, j)) * eta2(i, j), KIND=DP)
       ENDDO
    ENDDO
    IF (gstart == 2) THEN
       DO j = 1, nbnd
          s = s - REAL(eta1(1, j), KIND=DP) * REAL(eta2(1, j), KIND=DP)
       ENDDO
    ENDIF
    ! Plane-wave parallel reduction over the band group.
    CALL mp_sum(s, intra_bgrp_comm)
    stiefel_inner_product_gamma = s
    !
  END FUNCTION stiefel_inner_product_gamma
  !
END MODULE rdmft_stiefel
