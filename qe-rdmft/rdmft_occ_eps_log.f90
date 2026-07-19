!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_occ_eps_log
  !------------------------------------------------------------------
  !! Per-state occupation derivatives during RDMFT occupation
  !! optimisation.
  !!
  !! Writes ``{tmpdir}/{prefix}.rdmft.occ_eps`` (text, ionode
  !! only) with one row per (ik, ib) at each accepted inner
  !! occupation step and a final converged snapshot (``inner = 0``).
  !! Each snapshot is appended as a new block (file is created once
  !! at the start of ``rdmft_run``; subsequent records use
  !! ``POSITION='APPEND'``).
  !!
  !! ``dedn`` matches ELK ``RDM_DEDN.OUT`` column 3:
  !! \(\partial E / \partial n_{ik}\) without the BZ weight \(w_k\)
  !! Logged column is physical \(\partial E/\partial n_{ik}\) from
  !! :func:`rdmft_grad_n`.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_occ_eps_log_open, rdmft_occ_eps_log_close, &
            rdmft_occ_eps_log_set_outer, rdmft_occ_eps_log_record
  !
  CHARACTER(LEN=512), SAVE :: occ_eps_fname = ' '
  LOGICAL, SAVE :: occ_eps_active = .FALSE.
  INTEGER, SAVE :: occ_eps_outer = 0
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_occ_eps_log_open()
    !---------------------------------------------------------------
    !! Create the log file and write the header once.  Collective:
    !! every rank must call; ``occ_eps_active`` is broadcast from
    !! ionode so all ranks know whether logging succeeded.
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE io_files,  ONLY : prefix, tmp_dir
    USE mp,        ONLY : mp_bcast
    USE mp_images, ONLY : intra_image_comm
    !
    INTEGER :: ios, iu
    LOGICAL :: active_loc
    !
    active_loc = .FALSE.
    IF (ionode) THEN
       occ_eps_fname = TRIM(tmp_dir) // TRIM(prefix) // '.rdmft.occ_eps'
       OPEN(NEWUNIT=iu, FILE=TRIM(occ_eps_fname), STATUS='REPLACE', &
            FORM='FORMATTED', ACTION='WRITE', IOSTAT=ios)
       IF (ios /= 0) THEN
          WRITE(stdout, '(5X,A,A,A,I0)') &
               'RDMFT: could not open occupation grad log ', TRIM(occ_eps_fname), &
               ' (ios=', ios, ')'
       ELSE
          WRITE(iu, '(A)') '# RDMFT occupation dedn per inner step (ELK RDM_DEDN.OUT convention)'
          WRITE(iu, '(A)') &
               '# outer  inner  ik  ib  spin  n  wk  dedn(Ry)'
          WRITE(iu, '(A)') &
               '# dedn = (dE/dn)_ik = grad_n/wk (ELK RDM_DEDN.OUT column 3)'
          WRITE(iu, '(A)') &
               '# inner = -1 initial (pre-optimisation); inner = 0 final (post-sort)'
          CLOSE(iu)
          active_loc = .TRUE.
          WRITE(stdout, '(5X,A,A)') &
               'RDMFT: occupation gradient log -> ', TRIM(occ_eps_fname)
       ENDIF
    ENDIF
    CALL mp_bcast(active_loc, ionode_id, intra_image_comm)
    occ_eps_active = active_loc
    !
  END SUBROUTINE rdmft_occ_eps_log_open
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_occ_eps_log_close()
    !---------------------------------------------------------------
    occ_eps_active = .FALSE.
    occ_eps_fname = ' '
    !
  END SUBROUTINE rdmft_occ_eps_log_close
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_occ_eps_log_set_outer(outer)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN) :: outer
    occ_eps_outer = outer
  END SUBROUTINE rdmft_occ_eps_log_set_outer
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_occ_eps_log_record(inner, grad_local)
    !---------------------------------------------------------------
  !! Gather per-pool ``grad_local`` and ``rdmft_n`` and append one
  !! snapshot block.  Every MPI rank must call (``poolcollect`` is
  !! collective); only ionode appends to disk.
    USE io_global, ONLY : ionode, stdout
    USE wvfct,     ONLY : nbnd
    USE klist,     ONLY : nks, nkstot, wk
    USE lsda_mod,  ONLY : lsda, isk
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : world_comm
    USE rdmft_module, ONLY : rdmft_n
    USE rdmft_energy, ONLY : rdmft_dedn_from_grad_n
    !
    INTEGER, INTENT(IN) :: inner
    REAL(DP), INTENT(IN) :: grad_local(nbnd, nks)
    !
    REAL(DP), ALLOCATABLE :: n_g(:,:), grad_g(:,:), dedn_g(:,:)
    REAL(DP), ALLOCATABLE :: wk_loc(:,:), wk_g(:,:), isk_loc(:,:), isk_g(:,:)
    INTEGER :: ik, ib, ispin, iu, ios
    CHARACTER(LEN=9) :: slbl
    !
    IF (.NOT. occ_eps_active) RETURN
    !
    ! Cross-pool synchronisation before any ``inter_pool_comm`` allreduce.
    ! Without this barrier, the per-image-pool ``exxinit_std`` broadcasts
    ! inside ``rdmft_grad_n`` can leave one image-pool finishing earlier
    ! than another; the early image-pool then enters the
    ! ``poolcollect`` allreduce and blocks on its ``inter_pool_comm``
    ! partner that is still inside ``exxinit_std``.  The barrier on
    ! ``world_comm`` is the only comm shared by every rank regardless
    ! of the ``-nk`` / image decomposition, so it forces the call
    ! frame to be uniform before any cross-pool collective.
    CALL mp_barrier(world_comm)
    !
    ALLOCATE(n_g(nbnd, nkstot), grad_g(nbnd, nkstot), dedn_g(nbnd, nkstot))
    ALLOCATE(wk_loc(1, nks), wk_g(1, nkstot), isk_loc(1, nks), isk_g(1, nkstot))
    DO ik = 1, nks
       wk_loc(1, ik)  = wk(ik)
       isk_loc(1, ik) = REAL(isk(ik), DP)
    ENDDO
    CALL poolcollect(nbnd, nks, rdmft_n,    nkstot, n_g)
    CALL poolcollect(nbnd, nks, grad_local,  nkstot, grad_g)
    CALL poolcollect(1, nks, wk_loc,  nkstot, wk_g)
    CALL poolcollect(1, nks, isk_loc, nkstot, isk_g)
    !
    BLOCK
       REAL(DP) :: wk_vec(nkstot)
       wk_vec = wk_g(1, :)
       CALL rdmft_dedn_from_grad_n(grad_g, wk_vec, nbnd, nkstot, dedn_g)
    END BLOCK
    !
    IF (ionode) THEN
       OPEN(NEWUNIT=iu, FILE=TRIM(occ_eps_fname), STATUS='OLD', &
            POSITION='APPEND', FORM='FORMATTED', ACTION='WRITE', IOSTAT=ios)
       IF (ios /= 0) THEN
          WRITE(stdout, '(5X,A,A,A,I0)') &
               'RDMFT: occ_eps append failed for ', TRIM(occ_eps_fname), &
               ' (ios=', ios, ')'
       ELSE
          WRITE(iu, '(A,I0,A,I0,A)') '# --- outer=', occ_eps_outer, &
               ' inner=', inner, ' ---'
          DO ik = 1, nkstot
             ispin = NINT(isk_g(1, ik))
             slbl = MERGE('spin up  ', 'spin down', .NOT. (lsda .AND. ispin == 2))
             IF (.NOT. lsda) slbl = 'spinless '
             DO ib = 1, nbnd
                WRITE(iu, '(I6,2X,I6,2X,I6,2X,I6,2X,A9,2X,ES18.10,2X,ES18.10,2X,ES18.10)') &
                     occ_eps_outer, inner, ik, ib, slbl, &
                     n_g(ib, ik), wk_g(1, ik), dedn_g(ib, ik)
             ENDDO
          ENDDO
          CLOSE(iu)
       ENDIF
    ENDIF
    !
    DEALLOCATE(n_g, grad_g, dedn_g, wk_loc, wk_g, isk_loc, isk_g)
    !
  END SUBROUTINE rdmft_occ_eps_log_record
  !
END MODULE rdmft_occ_eps_log
