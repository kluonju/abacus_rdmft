!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_input
  !------------------------------------------------------------------
  !! Read the optional ``&RDMFT`` namelist that controls the RDMFT
  !! solver.  The namelist is appended at the end of the standard
  !! ``pw.x`` input (after the SYSTEM / CONTROL / ELECTRONS blocks)
  !! and is only read here -- not in
  !! \texttt{Modules/read\_namelists.f90} -- to keep the patch local
  !! to the new \texttt{PW/src/rdmft/} folder.
  !!
  !! Example::
  !!
  !!     &control
  !!       calculation = 'scf'
  !!       prefix      = 'h2'
  !!       outdir      = './tmp'
  !!     /
  !!     ... &system, &electrons, ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS ...
  !!     &rdmft
  !!       do_rdmft           = .true.
  !!       rdmft_functional   = 'muller'
  !!       rdmft_constraint   = 'projected_gradient'
  !!       rdmft_outer_maxiter = 30
  !!       rdmft_grad_check    = .true.
  !!     /
  !!
  !! Keyword names mirror the ABACUS ``rdmft_*`` INPUT_PARAMETERS
  !! conventions one-to-one (see PW/src/rdmft/README.md).
  !
  USE kinds, ONLY : DP
  USE rdmft_module
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_read_input
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_read_input()
    !---------------------------------------------------------------
    !! Read the &rdmft namelist from the pw.x input file.
    !!
    !! QE reads its input into a private file unit (``qestdin``) that
    !! is already closed by the time ``iosys()`` calls this routine.
    !! To find the optional ``&rdmft`` block we reopen, by name, the
    !! very file QE connected to ``qestdin`` (returned by
    !! ``get_input_file_name``):
    !!
    !!   * ``pw.x -i input.in``  -> the named file (left on disk);
    !!   * ``pw.x < input.in``   -> the ``input_tmp.in`` scratch copy
    !!     of standard input, which ``read_input_file`` now keeps alive
    !!     for us and which we delete again here.
    !!
    !! Reopening by name (rather than the old ``REWIND(5)`` trick)
    !! is essential: unit 5 is a non-seekable pipe under
    !! ``mpirun ... pw.x < input.in``, so ``REWIND(5)`` aborts with a
    !! fatal "Illegal seek" runtime error -- which broke every MPI
    !! (``-nk`` / ``-npool``) run started from stdin.
    !!
    !! If no ``&rdmft`` block is found the routine returns silently
    !! with \texttt{do\_rdmft = .FALSE.} leaving all defaults intact.
    !
    USE io_global,              ONLY : ionode, ionode_id, stdout
    USE mp,                     ONLY : mp_bcast
    USE mp_images,              ONLY : intra_image_comm
    USE open_close_input_file,  ONLY : get_input_file_name
    !
    NAMELIST /rdmft/ do_rdmft, rdmft_functional, rdmft_power_alpha, &
                      rdmft_reg_eps, rdmft_solver_strategy, rdmft_block_order, &
                      rdmft_constraint, rdmft_occ_optimizer, &
                      rdmft_orb_optimizer, rdmft_joint_optimizer, &
                      rdmft_lbfgs_memory, rdmft_bgd_tau, rdmft_bgd_backtrack, &
                      rdmft_outer_maxiter, &
                      rdmft_occ_maxiter, rdmft_orb_maxiter, &
                      rdmft_energy_tol, rdmft_orb_grad_tol, &
                      rdmft_occ_grad_tol, rdmft_occ_tol, &
                      rdmft_line_search, rdmft_zhang_hager_eta, &
                      rdmft_line_search_c1, rdmft_line_search_c2, &
                      rdmft_line_search_rho, rdmft_line_search_max_iter, &
                      rdmft_line_search_polynomial, &
                      rdmft_line_search_max_zoom, &
                      rdmft_occ_ls_type, rdmft_orb_ls_type, rdmft_occ_ls_init_step, &
                      rdmft_orb_ls_init_step, &
                      rdmft_bb_alpha_min, rdmft_bb_alpha_max, &
                      rdmft_orb_ls_stepsize, rdmft_occ_ls_stepsize, &
                      rdmft_orb_no_ls, rdmft_orb_strategy, &
                      rdmft_occ_init_mode, rdmft_occ_init_perturb, &
                      rdmft_occ_init_nbands_top, rdmft_verbose, &
                      rdmft_grad_check, rdmft_temp, &
                      rdmft_occ_precond, rdmft_occ_precond_shift, &
                      rdmft_orb_precond, rdmft_orb_precond_shift, &
                      rdmft_restart, rdmft_save_every, rdmft_restart_file, &
                      rdmft_source_prefix, &
                      rdmft_source_outdir, rdmft_krefine_orb_maxiter, &
                      rdmft_krefine_orb_tol
    !
    INTEGER            :: ios, iu
    LOGICAL            :: is_tmp
    CHARACTER(LEN=256) :: input_fname
    !
    ios = -1
    IF (ionode) THEN
       !
       ! Reopen, on a private unit, the file that QE connected to
       ! qestdin while reading the standard namelists, and scan it for
       ! the optional &rdmft block.  This file is either the
       ! command-line "-i input.in" file or the "input_tmp.in" scratch
       ! copy of standard input (kept alive for us by read_input_file).
       !
       ! Reopening by name is deliberate: the previous implementation
       ! used REWIND(5) on standard input, which raises a fatal
       ! "Illegal seek" runtime error whenever stdin is a non-seekable
       ! pipe -- exactly what happens under "mpirun ... pw.x < input.in"
       ! (and hence with any -nk / -npool run driven from stdin).
       !
       input_fname = get_input_file_name()
       is_tmp = ( TRIM(input_fname) == 'input_tmp.in' )
       iu = 98
       IF (TRIM(input_fname) /= ' ') THEN
          OPEN(UNIT=iu, FILE=TRIM(input_fname), FORM='FORMATTED', &
               STATUS='OLD', IOSTAT=ios)
          IF (ios == 0) THEN
             READ(iu, rdmft, IOSTAT=ios)
             ! Delete the stdin scratch copy now that we are done with
             ! it (restoring QE's usual cleanup); never touch a real,
             ! user-provided input file.
             IF (is_tmp) THEN
                CLOSE(iu, STATUS='delete')
             ELSE
                CLOSE(iu, STATUS='keep')
             END IF
          ELSE IF (is_tmp) THEN
             ! Could not reopen the scratch file: remove it if present
             ! so it does not linger.
             OPEN(UNIT=iu, FILE=TRIM(input_fname), FORM='FORMATTED', &
                  STATUS='OLD', IOSTAT=ios)
             IF (ios == 0) CLOSE(iu, STATUS='delete')
             ios = 1
          END IF
       END IF
       IF (ios /= 0) do_rdmft = .FALSE.
    ENDIF
    !
    ! Broadcast all RDMFT control variables.
    CALL mp_bcast(do_rdmft,                  ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_functional,          ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_power_alpha,         ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_reg_eps,             ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_solver_strategy,     ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_block_order,         ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_constraint,          ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_optimizer,       ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_orb_optimizer,       ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_joint_optimizer,     ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_lbfgs_memory,        ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_bgd_tau,              ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_bgd_backtrack,        ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_outer_maxiter,       ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_maxiter,         ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_orb_maxiter,         ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_energy_tol,          ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_orb_grad_tol,        ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_grad_tol,        ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_tol,             ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_line_search,         ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_zhang_hager_eta,     ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_line_search_c1,      ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_line_search_c2,      ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_line_search_rho,     ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_line_search_polynomial, ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_line_search_max_iter, ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_line_search_max_zoom, ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_ls_type,         ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_orb_ls_type,         ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_ls_init_step,   ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_orb_ls_init_step,   ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_bb_alpha_min,       ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_bb_alpha_max,       ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_orb_ls_stepsize,     ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_ls_stepsize,     ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_orb_no_ls,           ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_orb_strategy,        ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_init_mode,       ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_init_perturb,    ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_init_nbands_top, ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_verbose,             ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_grad_check,          ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_temp,                ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_precond,         ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_occ_precond_shift,   ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_orb_precond,         ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_orb_precond_shift,   ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_restart,           ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_save_every,        ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_restart_file,      ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_source_prefix,     ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_source_outdir,       ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_krefine_orb_maxiter, ionode_id, intra_image_comm)
    CALL mp_bcast(rdmft_krefine_orb_tol,     ionode_id, intra_image_comm)
    !
    IF (do_rdmft .AND. ionode) THEN
       WRITE(stdout, '(/,5X,A)') '&rdmft namelist found: post-SCF RDMFT optimisation enabled.'
       WRITE(stdout, '(5X,A,A)') '   rdmft_functional      = ', TRIM(rdmft_functional)
       WRITE(stdout, '(5X,A,A)') '   rdmft_solver_strategy = ', TRIM(rdmft_solver_strategy)
       WRITE(stdout, '(5X,A,A)') '   rdmft_constraint      = ', TRIM(rdmft_constraint)
       IF (TRIM(rdmft_solver_strategy) == 'joint') THEN
          WRITE(stdout, '(5X,A,A)') '   rdmft_joint_optimizer = ', TRIM(rdmft_joint_optimizer)
       ELSE
          WRITE(stdout, '(5X,A,A)') '   rdmft_block_order     = ', TRIM(rdmft_block_order)
          WRITE(stdout, '(5X,A,A)') '   rdmft_occ_ls_type     = ', TRIM(rdmft_occ_ls_type)
          WRITE(stdout, '(5X,A,A)') '   rdmft_orb_ls_type     = ', TRIM(rdmft_orb_ls_type)
          WRITE(stdout, '(5X,A,A)') '   rdmft_orb_optimizer   = ', TRIM(rdmft_orb_optimizer)
       ENDIF
       WRITE(stdout, '(5X,A,L1)') '   rdmft_occ_precond     = ', rdmft_occ_precond
       WRITE(stdout, '(5X,A)') '   occ XC gradient       = closed-form c_t wk w_t''(n) D_ii'
       IF (rdmft_restart) &
            WRITE(stdout, '(5X,A,L1)') '   rdmft_restart         = ', rdmft_restart
       IF (LEN_TRIM(rdmft_source_prefix) > 0) &
            WRITE(stdout, '(5X,A,A)') '   rdmft_source_prefix   = ', &
            TRIM(rdmft_source_prefix)
       IF (LEN_TRIM(rdmft_source_outdir) > 0) &
            WRITE(stdout, '(5X,A,A)') '   rdmft_source_outdir   = ', &
            TRIM(rdmft_source_outdir)
       IF (LEN_TRIM(rdmft_source_prefix) > 0 .AND. &
           LEN_TRIM(rdmft_source_outdir) > 0) &
            WRITE(stdout, '(5X,A)') &
                 '   (k-refine auto-enabled when K_POINTS is a finer mesh)'
       IF (rdmft_save_every /= 1) &
            WRITE(stdout, '(5X,A,I0)') '   rdmft_save_every      = ', rdmft_save_every
    ENDIF
    !
  END SUBROUTINE rdmft_read_input
  !
END MODULE rdmft_input
