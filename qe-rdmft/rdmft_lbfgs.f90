!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_lbfgs
  !------------------------------------------------------------------
  !! Generic limited-memory BFGS in the standard Euclidean metric.
  !!
  !! Used by the SPG occupation block (\texttt{rdmft\_occ\_optimizer =
  !! 'lbfgs'}), the joint solver, and the orbital block (via
  !! re-projection on the Stiefel tangent space).  The two-loop
  !! recursion is the textbook Nocedal \& Wright Algorithm 7.4.
  !!
  !! API:
  !!   CALL rdmft_lbfgs_init(state, n, memory)
  !!   CALL rdmft_lbfgs_reset(state)
  !!   CALL rdmft_lbfgs_direction(state, grad, dir)        ! dir = -H grad
  !!   CALL rdmft_lbfgs_update(state, step_vec, grad_diff)
  !!   CALL rdmft_lbfgs_finalize(state)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_lbfgs_state, rdmft_lbfgs_init, rdmft_lbfgs_reset, &
            rdmft_lbfgs_direction, rdmft_lbfgs_update, rdmft_lbfgs_finalize
  !
  TYPE :: rdmft_lbfgs_state
     INTEGER :: ndim = 0
     INTEGER :: m_max = 10
     INTEGER :: n_stored = 0
     REAL(DP) :: gamma_init = 1.0_DP
     ! Ring buffer of (s, y, rho).  Slot 1..n_stored holds the most
     ! recent n_stored pairs in chronological order (oldest first).
     REAL(DP), ALLOCATABLE :: s(:,:)
     REAL(DP), ALLOCATABLE :: y(:,:)
     REAL(DP), ALLOCATABLE :: rho_buf(:)
  END TYPE rdmft_lbfgs_state
  !
CONTAINS
  !
  SUBROUTINE rdmft_lbfgs_init(state, ndim, memory)
    TYPE(rdmft_lbfgs_state), INTENT(INOUT) :: state
    INTEGER, INTENT(IN) :: ndim, memory
    state%ndim   = ndim
    state%m_max  = MAX(1, memory)
    state%n_stored = 0
    state%gamma_init = 1.0_DP
    IF (ALLOCATED(state%s)) DEALLOCATE(state%s)
    IF (ALLOCATED(state%y)) DEALLOCATE(state%y)
    IF (ALLOCATED(state%rho_buf)) DEALLOCATE(state%rho_buf)
    ALLOCATE(state%s(ndim, state%m_max))
    ALLOCATE(state%y(ndim, state%m_max))
    ALLOCATE(state%rho_buf(state%m_max))
    state%s = 0.0_DP
    state%y = 0.0_DP
    state%rho_buf = 0.0_DP
  END SUBROUTINE rdmft_lbfgs_init
  !
  SUBROUTINE rdmft_lbfgs_reset(state)
    TYPE(rdmft_lbfgs_state), INTENT(INOUT) :: state
    state%n_stored = 0
    state%gamma_init = 1.0_DP
    IF (ALLOCATED(state%s)) state%s = 0.0_DP
    IF (ALLOCATED(state%y)) state%y = 0.0_DP
    IF (ALLOCATED(state%rho_buf)) state%rho_buf = 0.0_DP
  END SUBROUTINE rdmft_lbfgs_reset
  !
  SUBROUTINE rdmft_lbfgs_finalize(state)
    TYPE(rdmft_lbfgs_state), INTENT(INOUT) :: state
    state%n_stored = 0
    state%ndim = 0
    IF (ALLOCATED(state%s)) DEALLOCATE(state%s)
    IF (ALLOCATED(state%y)) DEALLOCATE(state%y)
    IF (ALLOCATED(state%rho_buf)) DEALLOCATE(state%rho_buf)
  END SUBROUTINE rdmft_lbfgs_finalize
  !
  SUBROUTINE rdmft_lbfgs_direction(state, grad, dir)
    !! Two-loop recursion: returns dir = -H grad with H built from the
    !! stored pairs.  When the buffer is empty falls back to dir = -grad.
    TYPE(rdmft_lbfgs_state), INTENT(IN) :: state
    REAL(DP), INTENT(IN)  :: grad(:)
    REAL(DP), INTENT(OUT) :: dir(:)
    REAL(DP), ALLOCATABLE :: q(:), alpha_vec(:)
    REAL(DP) :: alpha_k, beta_k
    INTEGER :: k
    !
    IF (state%n_stored == 0) THEN
       dir = -grad
       RETURN
    ENDIF
    ALLOCATE(q(state%ndim), alpha_vec(state%n_stored))
    q = grad
    DO k = state%n_stored, 1, -1
       alpha_k = state%rho_buf(k) * SUM(state%s(:, k) * q)
       alpha_vec(k) = alpha_k
       q = q - alpha_k * state%y(:, k)
    ENDDO
    q = state%gamma_init * q
    DO k = 1, state%n_stored
       beta_k = state%rho_buf(k) * SUM(state%y(:, k) * q)
       q = q + (alpha_vec(k) - beta_k) * state%s(:, k)
    ENDDO
    dir = -q
    DEALLOCATE(q, alpha_vec)
  END SUBROUTINE rdmft_lbfgs_direction
  !
  SUBROUTINE rdmft_lbfgs_update(state, step_vec, grad_diff)
    !! Append a new (s, y) pair.  Skips when s.y < 1e-12 (curvature
    !! safeguard).  Updates the diagonal scaling gamma_init from the
    !! most recent pair.
    TYPE(rdmft_lbfgs_state), INTENT(INOUT) :: state
    REAL(DP), INTENT(IN) :: step_vec(:), grad_diff(:)
    REAL(DP) :: sy, yy
    INTEGER :: k
    sy = SUM(step_vec * grad_diff)
    yy = SUM(grad_diff * grad_diff)
    IF (sy <= 1.0e-12_DP) RETURN
    IF (state%n_stored == state%m_max) THEN
       ! drop the oldest pair (slot 1) by shifting left
       DO k = 1, state%m_max - 1
          state%s(:, k) = state%s(:, k+1)
          state%y(:, k) = state%y(:, k+1)
          state%rho_buf(k) = state%rho_buf(k+1)
       ENDDO
       state%s(:, state%m_max) = step_vec
       state%y(:, state%m_max) = grad_diff
       state%rho_buf(state%m_max) = 1.0_DP / sy
    ELSE
       state%n_stored = state%n_stored + 1
       state%s(:, state%n_stored) = step_vec
       state%y(:, state%n_stored) = grad_diff
       state%rho_buf(state%n_stored) = 1.0_DP / sy
    ENDIF
    IF (yy > 1.0e-30_DP) state%gamma_init = sy / yy
  END SUBROUTINE rdmft_lbfgs_update
  !
END MODULE rdmft_lbfgs
