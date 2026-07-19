!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_inner_log
  !------------------------------------------------------------------
  !! Per-inner-iteration energy trace for the alternating RDMFT
  !! occupation and orbital blocks.
  !!
  !! Writes ``{outdir}/{prefix}.rdmft.inner_energy`` (text, ionode
  !! only) with one row per accepted inner occupation step and per
  !! orbital inner iteration (joint or block-k).
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_inner_log_open, rdmft_inner_log_close, &
            rdmft_inner_log_set_outer, rdmft_inner_log_record
  !
  INTEGER, SAVE :: inner_log_iu = 0
  LOGICAL, SAVE :: inner_log_open = .FALSE.
  INTEGER, SAVE :: inner_log_outer = 0
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_inner_log_open()
    !---------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE io_files,  ONLY : prefix, tmp_dir
    !
    CHARACTER(LEN=512) :: fname
    INTEGER :: ios
    !
    IF (.NOT. ionode) RETURN
    IF (inner_log_open) RETURN
    fname = TRIM(tmp_dir) // TRIM(prefix) // '.rdmft.inner_energy'
    OPEN(NEWUNIT=inner_log_iu, FILE=TRIM(fname), STATUS='REPLACE', &
         FORM='FORMATTED', ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) THEN
       WRITE(stdout, '(5X,A,A)') &
            'RDMFT: could not open inner energy log ', TRIM(fname)
       RETURN
    ENDIF
    inner_log_open = .TRUE.
    WRITE(inner_log_iu, '(A)') '# RDMFT inner occupation / orbital energies'
    WRITE(inner_log_iu, '(A)') &
         '# outer  block  inner  ik  E(Ry)  metric1  step  sum_dn  Ne  KKT'
    WRITE(inner_log_iu, '(A)') &
         '# block = occ | orb; ik = 0 (joint orb) or k-index (block_k orb)'
    WRITE(inner_log_iu, '(A)') &
         '# occ bgd: metric1 = phi0,        step = tau'
    WRITE(inner_log_iu, '(A)') &
         '# occ spg: metric1 = ||g_1||_inf, step = alpha'
    WRITE(inner_log_iu, '(A)') &
         '# orb:     metric1 = ||G_R|| or ||G_R^k||, step = alpha'
    !
  END SUBROUTINE rdmft_inner_log_open
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_inner_log_close()
    !---------------------------------------------------------------
    USE io_global, ONLY : ionode
    !
    IF (.NOT. ionode .OR. .NOT. inner_log_open) RETURN
    CLOSE(inner_log_iu)
    inner_log_open = .FALSE.
    inner_log_iu = 0
    !
  END SUBROUTINE rdmft_inner_log_close
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_inner_log_set_outer(outer)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN) :: outer
    inner_log_outer = outer
  END SUBROUTINE rdmft_inner_log_set_outer
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_inner_log_record(block, inner, ik, etot, metric1, &
                                  step, sum_dn, ne_val, kkt_resid)
    !---------------------------------------------------------------
  !! Append one inner-iteration row.  ``metric1`` / ``step`` mirror the
  !! stdout diagnostics (phi0 or gradient norm, line-search step).
  !! ``kkt_resid`` is the box+Ne KKT residual from ``rdmft_pg_kkt_residual``
  !! (occupation block only; zero for orbital rows).
    USE io_global, ONLY : ionode
    !
    CHARACTER(LEN=*), INTENT(IN) :: block
    INTEGER, INTENT(IN) :: inner, ik
    REAL(DP), INTENT(IN) :: etot, metric1, step, sum_dn, ne_val
    REAL(DP), INTENT(IN), OPTIONAL :: kkt_resid
    REAL(DP) :: kkt_loc
    !
    IF (.NOT. ionode .OR. .NOT. inner_log_open) RETURN
    kkt_loc = 0.0_DP
    IF (PRESENT(kkt_resid)) kkt_loc = kkt_resid
    WRITE(inner_log_iu, '(I6,2X,A3,2X,I6,2X,I6,2X,F24.16,2X,5ES18.10)') &
         inner_log_outer, block, inner, ik, etot, metric1, step, sum_dn, ne_val, kkt_loc
    !
  END SUBROUTINE rdmft_inner_log_record
  !
END MODULE rdmft_inner_log
