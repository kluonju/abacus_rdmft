!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_solver
  !------------------------------------------------------------------
  !! Top-level entry point for the alternating RDMFT optimisation.
  !!
  !! Implements the user-requested algorithmic combination:
  !!
  !! * **Outer loop**: alternate between an occupation step and an
  !!   orbital step (order set by \texttt{rdmft\_block\_order}).
  !! * **Occupation block**: ABACUS-style SPG using \(\partial E/\partial n_{ik}\)
  !!   from \texttt{rdmft\_grad\_n}, uniform-shift proximal projection,
  !!   optional ELK preconditioner, and KKT-residual stopping (monotone /
  !!   Wolfe line search).
  !! * **Orbital block**: Stiefel-manifold steepest-descent or
  !!   conjugate-gradient on the per-k-point natural orbital
  !!   coefficient matrix \(C^k = \mathrm{evc}(:,:,k)\), using QE's
  !!   plane-wave basis.  The retraction is the Cholesky-based polar
  !!   reorthonormalisation in \texttt{rdmft\_stiefel}; the line search
  !!   is monotone Armijo.
  !!
  !! On entry, the converged KS-SCF orbitals and occupations are taken
  !! as the initial guess.  On exit, the natural occupations and
  !! orbitals (plus the RDMFT total energy) are stored in
  !! \texttt{rdmft\_module}; QE's \texttt{etot} and \texttt{wg} are
  !! restored to their KS-SCF values so any post-RDMFT ``forces`` /
  !! ``stress`` driver still uses the proper KS quantities.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_run, rdmft_krefine_orbital_polish
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_run()
    !---------------------------------------------------------------
    !! Top-level RDMFT entry.  Dispatches between the alternating and
    !! the joint (product-manifold) drivers based on
    !! \texttt{rdmft\_solver\_strategy}.
    USE io_global,        ONLY : stdout, ionode
    USE wvfct,            ONLY : nbnd, wg
    USE klist,            ONLY : nks, nkstot, wk, xk
    USE io_files,         ONLY : prefix, tmp_dir
    USE uspp,             ONLY : nkb
    USE becmod,           ONLY : becp, allocate_bec_type, &
                                  deallocate_bec_type, is_allocated_bec_type
    USE xc_lib,           ONLY : xclib_get_exx_fraction, xclib_set_exx_fraction, &
                                  xclib_get_id, xclib_set_dft_IDs, xclib_dft_is, &
                                  exx_is_active, xclib_set_auxiliary_flags
    USE exx,              ONLY : use_ace
    USE exx_base,         ONLY : exxalfa, exx_grid_initialized, &
                                  exx_bgrp_type, EXX_BGRP_BANDS
    USE rdmft_module
    USE rdmft_xc,         ONLY : rdmft_xc_set_type, rdmft_xc_alpha, &
                                  RDMFT_XC_HF
    USE rdmft_energy,     ONLY : rdmft_initial_n_from_ks, &
                                  rdmft_total_energy, rdmft_exxinit_once, &
                                  rdmft_set_wg_for_exxinit, &
                                  rdmft_clear_band_diagnostics, rdmft_grad_n, &
                                  rdmft_write_dedn
    USE rdmft_gradient_check, ONLY : rdmft_check_gradient_consistency
    USE rdmft_io,         ONLY : rdmft_save_state, rdmft_load_state
    USE rdmft_inner_log,    ONLY : rdmft_inner_log_open, rdmft_inner_log_close
    USE rdmft_occ_eps_log,  ONLY : rdmft_occ_eps_log_open, rdmft_occ_eps_log_close, &
                                  rdmft_occ_eps_log_record
    USE rdmft_krefine_mod, ONLY : rdmft_krefine_for_restart, &
                                   rdmft_resolve_krefine_from_source, &
                                   rdmft_krefine_requested
    USE lsda_mod,         ONLY : lsda
    USE rdmft_occupation, ONLY : rdmft_weighted_sum
    !
    REAL(DP) :: exxalfa_save
    INTEGER  :: iexch_save, icorr_save, igcx_save, igcc_save, &
                imeta_save, imetac_save
    LOGICAL  :: dft_set_ok
    !
    REAL(DP), ALLOCATABLE :: wg_save(:,:)
    REAL(DP) :: etot_now
    !
    IF (.NOT. do_rdmft) RETURN
    !
    WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
    IF (TRIM(rdmft_solver_strategy) == 'joint') THEN
       WRITE(stdout, '(5X,A)')   'RDMFT joint (product-manifold) optimisation'
    ELSE
       WRITE(stdout, '(5X,A)')   'RDMFT alternating optimisation'
    ENDIF
    WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
    WRITE(stdout, '(5X,A,A)')        'rdmft_functional      = ', TRIM(rdmft_functional)
    CALL rdmft_xc_set_type(rdmft_functional)
    WRITE(stdout, '(5X,A,F8.4)')     'rdmft_power_alpha     = ', rdmft_xc_alpha()
    WRITE(stdout, '(5X,A,A)')        'rdmft_solver_strategy = ', TRIM(rdmft_solver_strategy)
    WRITE(stdout, '(5X,A,A)')        'rdmft_constraint      = ', TRIM(rdmft_constraint)
    IF (TRIM(rdmft_solver_strategy) == 'joint') THEN
       WRITE(stdout, '(5X,A,A)')     'rdmft_joint_optimizer = ', TRIM(rdmft_joint_optimizer)
       IF (TRIM(rdmft_joint_optimizer) == 'spg') &
            WRITE(stdout, '(5X,A)') &
            '  (SPG2 occupation chord x Stiefel CG orbital step, one joint line search)'
    ELSE
       WRITE(stdout, '(5X,A,A)')     'rdmft_block_order     = ', TRIM(rdmft_block_order)
       WRITE(stdout, '(5X,A,A)')     'rdmft_occ_optimizer   = ', TRIM(rdmft_occ_optimizer)
       WRITE(stdout, '(5X,A,A)')     'rdmft_orb_optimizer   = ', TRIM(rdmft_orb_optimizer)
    ENDIF
    WRITE(stdout, '(5X,A,I6)')       'rdmft_outer_maxiter   = ', rdmft_outer_maxiter
    WRITE(stdout, '(5X,A,I6)')       'rdmft_occ_maxiter     = ', rdmft_occ_maxiter
    WRITE(stdout, '(5X,A,I6)')       'rdmft_orb_maxiter     = ', rdmft_orb_maxiter
    WRITE(stdout, '(5X,A,1PE10.2)')  'rdmft_energy_tol      = ', rdmft_energy_tol
    WRITE(stdout, '(5X,A,1PE10.2)')  'rdmft_occ_tol         = ', rdmft_occ_tol
    WRITE(stdout, '(5X,A,A)')        'rdmft_line_search     = ', TRIM(rdmft_line_search)
    IF (TRIM(rdmft_line_search) == 'zhang_hager' .OR. &
        TRIM(rdmft_line_search) == 'zh'          .OR. &
        TRIM(rdmft_line_search) == 'nonmonotone') &
       WRITE(stdout, '(5X,A,F8.4)')     'rdmft_zhang_hager_eta = ', rdmft_zhang_hager_eta
    WRITE(stdout, '(5X,A,A)')        'rdmft_occ_ls_type     = ', TRIM(rdmft_occ_ls_type)
    WRITE(stdout, '(5X,A,A)')        'rdmft_orb_ls_type     = ', TRIM(rdmft_orb_ls_type)
    WRITE(stdout, '(5X,A,A)')        'rdmft_occ_ls_init_step = ', TRIM(rdmft_occ_ls_init_step)
    WRITE(stdout, '(5X,A,A)')        'rdmft_orb_ls_init_step = ', TRIM(rdmft_orb_ls_init_step)
    WRITE(stdout, '(5X,A,L1)') 'occ ELK precond.      = ', rdmft_occ_precond
    WRITE(stdout, '(5X,A)') 'occ XC gradient      = closed-form  c_t wk w_t''(n) D_ii'
    WRITE(stdout, '(5X,A,L1,A,F8.4)') 'orb level-shift prec. = ', rdmft_orb_precond, &
                                       '  shift = ', rdmft_orb_precond_shift
    !
    ! ``wg`` is allocated (nbnd, nkstot) but only its first nks columns
    ! belong to this pool; slice explicitly so the save/restore is
    ! conformable under k-point pools (npool > 1).
    ALLOCATE(wg_save(nbnd, nks))
    wg_save = wg(1:nbnd, 1:nks)
    !
    ! Take ownership of ``wg`` for the entire RDMFT solve: this stops
    ! ``sum_band`` from re-running ``weights()`` (smearing, fixed_occ,
    ! tetrahedra) every time the RDMFT energy evaluator refreshes the
    ! density.  Without this, the natural-occupation weights we set in
    ! :subroutine:`rdmft_set_wg_from_n` are silently overwritten by the
    ! standard KS Fermi-Dirac fill, so ``rho`` (and therefore ``ehart``
    ! and the one-body diagonal) become independent of ``rdmft_n`` and
    ! the optimiser cannot make any progress on the occupation block.
    rdmft_active = .TRUE.
    CALL rdmft_clear_band_diagnostics()
    !
    ! Allocate the (nkb, nbnd) becp workspace ONCE for the whole
    ! RDMFT run.  All subsequent h_psi calls reuse it; some libc
    ! allocators trip on heap corruption when becp is alloc /
    ! deallocated inside a tight inner loop.
    IF (nkb > 0 .AND. .NOT. is_allocated_bec_type(becp)) &
       CALL allocate_bec_type(nkb, nbnd, becp)
    !
    ! For RDMFT functionals that mix HF-like channels (all eight
    ! supported XC types do), force the global EXX fraction to 1.0
    ! during the RDMFT solve so that vexxace_* returns the full
    ! (un-scaled) Fock operator.  When the SCF was run with a hybrid
    ! functional (e.g. PBE0 with exxalfa = 0.25) we need to override
    ! this; otherwise the channel weights would be off by a factor
    ! exxalfa.  Save the original fraction and restore it on exit.
    !
    ! Note that exxinit only refreshes exxalfa from xclib on the
    ! FIRST call (when EXX is not yet active); afterwards we have
    ! to set the cached exxalfa in exx_base directly.
    exxalfa_save = xclib_get_exx_fraction()
    CALL xclib_set_exx_fraction(1.0_DP)
    exxalfa = 1.0_DP
    !
    ! All RDMFT XC kernels supported here build the XC energy purely
    ! from Fock-like channels; the semilocal V_xc that the SCF used
    ! (e.g. PBE for the PBE0 starter) MUST be removed from the KS
    ! Hamiltonian, otherwise the orbital gradient picks up a spurious
    ! V_xc^semilocal piece and the line search refuses to descend.
    ! Switch the global XC IDs to (5,0,0,0,0,0) = HF for the duration
    ! of the RDMFT run; this makes v_of_rho return V_xc = 0 while
    ! keeping ehart, vtxc, etxc consistent.  Restore on exit.
    iexch_save  = xclib_get_id('LDA','EXCH')
    icorr_save  = xclib_get_id('LDA','CORR')
    igcx_save   = xclib_get_id('GGA','EXCH')
    igcc_save   = xclib_get_id('GGA','CORR')
    imeta_save  = xclib_get_id('MGGA','EXCH')
    imetac_save = xclib_get_id('MGGA','CORR')
    dft_set_ok = xclib_set_dft_IDs(5, 0, 0, 0, 0, 0)
    ! Refresh the auxiliary flags (ishybrid, dft_is_gradient, ...)
    ! so that subsequent xc_lib calls see the new IDs.
    CALL xclib_set_auxiliary_flags(.FALSE.)
    !
    ! When the SCF was run with a non-hybrid functional (e.g. pure
    ! PBE), QE never ran setup_exx and the EXX module is dormant.
    ! Set up the EXX grid + symmetry machinery here so the
    ! per-channel ACE machinery used by RDMFT can compute V_x; this
    ! is essentially a manual replay of the work that QE's
    ! hybrid-loop driver does at the first EXX iteration.
    !
    ! Use `exx_grid_initialized` (NOT `exx_is_active()`) as the
    ! guard.  Two scenarios:
    !
    !   * pure PBE SCF: setup_exx was never called; both
    !     exx_grid_initialized and exx_is_active() are FALSE -- we
    !     have to run setup_exx here.
    !   * hybrid SCF (e.g. input_dft = 'hf' or 'pbe0'): setup.f90
    !     called setup_exx() at startup so exx_grid_initialized is
    !     already TRUE, but exx_is_active() can still be FALSE if
    !     the SCF exited before the first EXX iteration (e.g. with
    !     electron_maxstep small enough that the non-EXX first
    !     iteration did not converge).  In that case the old
    !     `.NOT. exx_is_active()` guard would invoke setup_exx a
    !     second time and exx_grid_init aborts with "grid already
    !     initialized".
    !
    ! NOTE: do *not* call start_exx() here.  exxinit's first-time
    ! branch is what computes exxdiv = exx_divergence() (which the
    ! Gygi-Baldereschi correction needs); that branch is only
    ! entered when ``.NOT. exx_is_active()``.  Let the very first
    ! exxinit() call (driven by ``rdmft_exxinit_once`` /
    ! ``rdmft_compute_xc_channel`` below) do the activation.
    use_ace = .TRUE.
    ! Force the standard band-parallel EXX scheme.  Without this,
    ! Fortran zero-init or input defaults can leave exx_bgrp_type on
    ! band_pairs and route vexx through exx_bp.
    exx_bgrp_type = EXX_BGRP_BANDS
    IF (ionode) WRITE(stdout, '(5X,A)') &
         'RDMFT: EXX band parallelism forced to bands (EXX_BGRP_BANDS)'
    IF (.NOT. exx_grid_initialized) THEN
       CALL setup_exx()
    ENDIF
    !
    ! Seed the natural occupations from the KS converged weights, OR
    ! restore them from a previous run when ``rdmft_restart = .true.``.
    ! The natural orbitals are expected to already be in QE's ``evc``
    ! buffer at this point: with ``rdmft_restart = .true.`` the user
    ! must arrange for the standard PWscf restart path (e.g.
    ! ``startingwfc = 'file'`` in &electrons) so that
    ! :subroutine:`wfcinit` re-fills the buffer from the saved
    ! collected wavefunctions before electrons -> rdmft_run is
    ! entered.  Either way, ``rdmft_initial_n_from_ks`` first lays
    ! out target electron counts / spin constraints; ``rdmft_load_state``
    ! then overwrites ``rdmft_n`` with the previously saved
    ! occupations when the save file is valid.
    !
    ! IMPORTANT: state-loading MUST run BEFORE ``rdmft_exxinit_once``
    ! below.  ``exxinit`` sizes ``exxbuff`` from ``x_nbnd_occ`` (= the
    ! upper-band bound of the current ``wg``).  If we let ``exxinit``
    ! run with the SCF-flavoured ``wg`` only (just the few KS-occupied
    ! bands) and then -- inside the alternating loop -- the
    ! channel-1 RDMFT ``wg`` widens (because ``rdmft_load_state``
    ! restored a previously-converged set of natural occupations
    ! spread over many more bands), the next time ``exxinit`` runs
    ! it segfaults with an out-of-bounds write into ``exxbuff(1:nrxxs,
    ! ibnd, ikq)`` (``ibnd`` past the original ``ibnd_buff_end``).
    ! Loading the state and installing ``wg = wk * n`` here makes the
    ! first ``exxinit`` see the full RDMFT band support and allocate
    ! ``exxbuff`` to fit.  When ``rdmft_restart = .false.`` (the
    ! default) the saved-state branch is skipped, ``rdmft_n`` is just
    ! the KS occupations, and the ``wg`` pattern is identical to what
    ! the SCF left behind -- so this reorder is a no-op for non-restart
    ! runs.
    rdmft_krefine_done = .FALSE.
    CALL rdmft_resolve_krefine_from_source()
    IF (rdmft_krefine_requested()) THEN
       CALL rdmft_krefine_for_restart()
    ELSE
       CALL rdmft_initial_n_from_ks()
       IF (rdmft_restart) THEN
          BLOCK
            LOGICAL :: loaded
            CALL rdmft_load_state(loaded)
            IF (.NOT. loaded .AND. ionode) THEN
               WRITE(stdout, '(5X,A)') &
                    '** RDMFT: ``rdmft_restart = .true.`` but no usable'
               WRITE(stdout, '(5X,A)') &
                    '   restart file found; continuing with KS-seeded'
               WRITE(stdout, '(5X,A)') &
                    '   natural occupations.'
            ENDIF
          END BLOCK
       ENDIF
    ENDIF
    !
    ! Size ``exxbuff`` for the channel-weighted EXX support (``wk *
    ! w_t(n)``), not ``wk * n`` alone: ``pow_g`` etc. can be nonzero
    ! below ``rdmft_reg_eps`` when ``n=0``.
    CALL rdmft_set_wg_for_exxinit()
    !
    ! Heavy EXX setup (symmetry rotations, FFT grid, exxbuff) is done
    ! ONCE here; the per-channel updates inside the inner loops only
    ! refresh the lightweight x_occupation array.  Calling exxinit
    ! every iteration causes glibc free() heap corruption after a few
    ! hundred trips on multi-k complex pseudopotential paths.
    CALL rdmft_exxinit_once()
    !
    IF (rdmft_krefine_done .AND. rdmft_krefine_orb_maxiter > 0) &
         CALL rdmft_krefine_orbital_polish(rdmft_krefine_orb_maxiter, &
              rdmft_krefine_orb_tol, etot_now)
    !
    IF (ionode) THEN
       IF (rdmft_fix_magnetization .AND. lsda) THEN
          WRITE(stdout, '(5X,A)') 'occupation constraints: fixed nelup / neldw (tot_magnetization)'
          WRITE(stdout, '(5X,A,F12.6,A,F12.6)') '   N_up target = ', rdmft_n_target_up, &
               '   N_down target = ', rdmft_n_target_down
       ELSE
          WRITE(stdout, '(5X,A)') 'occupation constraints: total Ne only'
       ENDIF
    ENDIF
    !
    CALL rdmft_total_energy(etot_now)
    CALL rdmft_inner_log_open()
    CALL rdmft_occ_eps_log_open()
    BLOCK
       REAL(DP), ALLOCATABLE :: grad_init(:,:)
       ALLOCATE(grad_init(nbnd, nks))
       CALL rdmft_grad_n(grad_init)
       CALL rdmft_occ_eps_log_record(-1, grad_init)
       DEALLOCATE(grad_init)
    END BLOCK
    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A,F18.10,A)') '  Initial RDMFT total energy = ', etot_now, ' Ry'
       WRITE(stdout, '(5X,A,F18.10)')     '  E_one (T+V_loc+V_NL)       = ', rdmft_e_one
       WRITE(stdout, '(5X,A,F18.10)')     '  E_Hartree                  = ', rdmft_e_har
       WRITE(stdout, '(5X,A,F18.10)')     '  E_xc_RDMFT                 = ', rdmft_e_xc
       WRITE(stdout, '(5X,A,F18.10)')     '  E_Ewald                    = ', rdmft_e_const
       IF (rdmft_e_entropy /= 0.0_DP) &
            WRITE(stdout, '(5X,A,F8.2,A,F18.10,A)') &
            '  E_entropy (HF, T=', rdmft_temp, ' K) = ', rdmft_e_entropy, ' Ry'
    ENDIF
    !
    IF (rdmft_grad_check) CALL rdmft_check_gradient_consistency()
    !
    BLOCK
       REAL(DP) :: t_rdmft0
       IF (rdmft_verbose >= 1) t_rdmft0 = rdmft_wall_time()
       IF (TRIM(rdmft_solver_strategy) == 'joint') THEN
          CALL rdmft_run_joint(etot_now)
       ELSE
          CALL rdmft_run_alternating(etot_now)
       ENDIF
       IF (rdmft_verbose >= 1) &
            CALL rdmft_report_wall_time('RDMFT optimisation total', t_rdmft0)
    END BLOCK
    !
    ! Sort the natural orbitals so the printed occupations match the
    ! conventional descending order (most-occupied first).  Without
    ! this step, the orbital block's Stiefel rotation produces a
    ! natural-orbital basis whose 1-RDM eigenvalues n_ik are NOT
    ! monotonic in the band index ib -- mathematically correct, but
    ! a common source of confusion ("I expected n = [1,1,1,0] but
    ! got [0,0,1,1]").  The sort is a unitary reshuffling on a
    ! single k-point, so the total energy is invariant.
    CALL rdmft_sort_natocc()
    !
    BLOCK
       REAL(DP), ALLOCATABLE :: grad_final(:,:), grad_g(:,:), occ_g(:,:), wk_g(:)
       REAL(DP), ALLOCATABLE :: xk_g(:,:), wk_loc(:,:), wk_col(:,:)
       CHARACTER(LEN=256) :: fprefix
       INTEGER :: ik
       !
       ALLOCATE(grad_final(nbnd, nks), grad_g(nbnd, nkstot), occ_g(nbnd, nkstot))
       ALLOCATE(wk_g(nkstot), wk_loc(1, nks), wk_col(1, nkstot), xk_g(3, nkstot))
       CALL rdmft_grad_n(grad_final)
       CALL rdmft_occ_eps_log_record(0, grad_final)
       CALL poolcollect(nbnd, nks, grad_final, nkstot, grad_g)
       CALL poolcollect(nbnd, nks, rdmft_n, nkstot, occ_g)
       DO ik = 1, nks
          wk_loc(1, ik) = wk(ik)
       ENDDO
       CALL poolcollect(1, nks, wk_loc, nkstot, wk_col)
       wk_g = wk_col(1, :)
       CALL poolcollect(3, nks, xk(:, 1:nks), nkstot, xk_g)
       fprefix = TRIM(tmp_dir) // TRIM(prefix)
       CALL rdmft_write_dedn(fprefix, grad_g, occ_g, wk_g, xk_g)
       DEALLOCATE(grad_final, grad_g, occ_g, wk_g, wk_loc, wk_col, xk_g)
    END BLOCK
    !
    ! Final report.
    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A,F18.10,A)') '!RDMFT_ETOTAL  = ', rdmft_etot, ' Ry'
       WRITE(stdout, '(5X,A,F18.10,A)')   ' E_one+E_H+E_NL = ', rdmft_e_one, ' Ry (T+V_loc+V_NL)'
       WRITE(stdout, '(5X,A,F18.10,A)')   ' E_Hartree      = ', rdmft_e_har, ' Ry'
       WRITE(stdout, '(5X,A,F18.10,A)')   ' E_xc_RDMFT     = ', rdmft_e_xc , ' Ry'
       WRITE(stdout, '(5X,A,F18.10,A)')   ' E_Ewald        = ', rdmft_e_const, ' Ry'
    ENDIF
    !
    ! Optional lightweight occupation table (no EXX / TSM work).
    IF (rdmft_verbose >= 1) CALL rdmft_print_natural_occupations()
    !
    ! Final RDMFT state checkpoint: dump occupations and natural
    ! orbitals to disk regardless of ``rdmft_save_every`` so the
    ! converged solution is always recoverable.
    CALL rdmft_inner_log_close()
    CALL rdmft_occ_eps_log_close()
    CALL rdmft_save_state(-1, rdmft_etot)
    !
    ! Release the becp workspace allocated at the top.
    IF (nkb > 0 .AND. is_allocated_bec_type(becp)) CALL deallocate_bec_type(becp)
    !
    ! Restore the original EXX fraction in case downstream code
    ! (forces, stress, post-processing) relies on it.
    CALL xclib_set_exx_fraction(exxalfa_save)
    exxalfa = exxalfa_save
    !
    ! Restore the original XC functional IDs.
    dft_set_ok = xclib_set_dft_IDs(iexch_save, icorr_save, igcx_save, &
                                    igcc_save, imeta_save, imetac_save)
    !
    ! Restore the KS occupation weights and release the RDMFT-owned
    ! ``wg`` guard so subsequent SCF / post-processing calls run the
    ! normal ``weights()`` path inside ``sum_band`` again.
    wg(1:nbnd, 1:nks) = wg_save
    DEALLOCATE(wg_save)
    rdmft_active = .FALSE.
    !
  END SUBROUTINE rdmft_run
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_alternating_occ_block(etot_now, occ_converged)
    !---------------------------------------------------------------
    !! One occupation (SPG) inner block of an alternating outer cycle.
    USE io_global,        ONLY : stdout, ionode
    USE rdmft_module
    USE rdmft_spg,        ONLY : rdmft_spg_occ_block
    USE rdmft_bgd,         ONLY : rdmft_bgd_occ_block
    USE rdmft_ebi,        ONLY : rdmft_ebi_occ_block
    !
    REAL(DP), INTENT(INOUT) :: etot_now
    LOGICAL, INTENT(OUT)    :: occ_converged
    !
    occ_converged = (rdmft_occ_maxiter <= 0)
    IF (rdmft_occ_maxiter > 0) THEN
       IF (TRIM(rdmft_occ_optimizer) == 'bgd' .OR. TRIM(rdmft_occ_optimizer) == 'gd') THEN
          ! ELK-style box-aware occupation gradient descent (rdmvaryn) + line search.
          ! ``gd`` is a deprecated alias for ``bgd``.
          CALL rdmft_bgd_occ_block(etot_now, occ_converged, rdmft_occ_maxiter, 'occ')
       ELSE IF (TRIM(rdmft_occ_optimizer) == 'ebi') THEN
          ! Explicit-by-implicit (EBI) erf parameterisation + gradient descent.
          CALL rdmft_ebi_occ_block(etot_now, occ_converged, rdmft_occ_maxiter, 'occ')
       ELSE
          CALL rdmft_spg_occ_block(etot_now, occ_converged, rdmft_occ_maxiter, 'occ')
       ENDIF
    ENDIF
  END SUBROUTINE rdmft_alternating_occ_block
  !
  SUBROUTINE rdmft_alternating_orb_block(etot_now, orb_converged)
    !---------------------------------------------------------------
    !! One orbital (Stiefel) inner block of an alternating outer cycle.
    USE io_global, ONLY : stdout, ionode
    USE rdmft_module
    !
    REAL(DP), INTENT(INOUT) :: etot_now
    LOGICAL, INTENT(OUT)    :: orb_converged
    !
    IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,A,A,A,A,I0,A)') &
         'ORB block (', TRIM(rdmft_orb_optimizer), &
         ' / Stiefel CG, strategy=', TRIM(rdmft_orb_strategy), &
         ', max ', rdmft_orb_maxiter, ' inner iters)'
    orb_converged = (rdmft_orb_maxiter <= 0)
    IF (rdmft_orb_maxiter > 0) THEN
       CALL rdmft_orbital_step(etot_now, orb_converged)
    ENDIF
  END SUBROUTINE rdmft_alternating_orb_block
  !
  SUBROUTINE rdmft_run_alternating(etot_now)
    !---------------------------------------------------------------
    !! Alternating RDMFT loop with the projected-gradient inner block
    !! on occupations and the Stiefel orbital block.
    !!
    !! Outer convergence (``outer > 1``): stop when **either**
    !! ``|dE| < rdmft_energy_tol`` **or** both inner blocks report
    !! convergence (OCC: ``||g_1||_inf < rdmft_occ_grad_tol``; ORB:
    !! ``||G_R|| < rdmft_orb_grad_tol``).
    USE io_global,        ONLY : stdout, ionode
    USE wvfct,            ONLY : nbnd
    USE klist,            ONLY : nks, wk
    USE lsda_mod,         ONLY : lsda, nspin
    USE rdmft_module
    USE rdmft_energy,     ONLY : rdmft_total_energy
    USE rdmft_spg,        ONLY : rdmft_spg_occ_block
    USE rdmft_occupation, ONLY : rdmft_weighted_sum, rdmft_weighted_sum_ispin
    USE rdmft_io,         ONLY : rdmft_save_state
    USE rdmft_inner_log,   ONLY : rdmft_inner_log_set_outer
    USE rdmft_occ_eps_log, ONLY : rdmft_occ_eps_log_set_outer
    USE mp,               ONLY : mp_sum
    USE mp_pools,         ONLY : inter_pool_comm
    !
    REAL(DP), INTENT(INOUT) :: etot_now
    !
    REAL(DP), ALLOCATABLE :: n_outer_save(:,:)
    REAL(DP) :: etot_prev, dE, occ_outer_dn_sum, ne_now, ne_up_now, ne_down_now, mz_now
    REAL(DP) :: t_outer0, t_occ0, t_orb0
    INTEGER :: outer
    LOGICAL :: occ_converged, orb_converged, orb_first, do_spin_report
    LOGICAL :: loop_converged
    !
    loop_converged = .FALSE.
    orb_first = (INDEX(rdmft_block_order, 'orb_occ') > 0)
    IF (.NOT. orb_first .AND. INDEX(rdmft_block_order, 'occ_orb') == 0) &
         CALL errore('rdmft_run_alternating', &
              'rdmft_block_order must be ''occ_orb'' or ''orb_occ''', 1)
    !
    ALLOCATE(n_outer_save(nbnd, nks))
    etot_prev = etot_now
    !
    DO outer = 1, rdmft_outer_maxiter
       CALL rdmft_inner_log_set_outer(outer)
       CALL rdmft_occ_eps_log_set_outer(outer)
       !
       IF (rdmft_verbose >= 1 .AND. ionode) THEN
          IF (orb_first) THEN
             WRITE(stdout, '(/,5X,A,I4,A,A,A,A)') &
                  '--- RDMFT alternating outer ', outer, ' (', &
                  TRIM(rdmft_orb_optimizer), ' orb / SPG occ) ---'
          ELSE
             WRITE(stdout, '(/,5X,A,I4,A,A,A)') &
                  '--- RDMFT alternating outer ', outer, ' (SPG occ / ', &
                  TRIM(rdmft_orb_optimizer), ' orb) ---'
          ENDIF
       ENDIF
       IF (rdmft_verbose >= 1) t_outer0 = rdmft_wall_time()
       !
       n_outer_save = rdmft_n
       !
       IF (orb_first) THEN
          IF (rdmft_verbose >= 1) t_orb0 = rdmft_wall_time()
          CALL rdmft_alternating_orb_block(etot_now, orb_converged)
          IF (rdmft_verbose >= 1) CALL rdmft_report_wall_time('ORB block', t_orb0)
          IF (rdmft_verbose >= 1) t_occ0 = rdmft_wall_time()
          CALL rdmft_alternating_occ_block(etot_now, occ_converged)
          IF (rdmft_verbose >= 1) CALL rdmft_report_wall_time('OCC block', t_occ0)
       ELSE
          IF (rdmft_verbose >= 1) t_occ0 = rdmft_wall_time()
          CALL rdmft_alternating_occ_block(etot_now, occ_converged)
          IF (rdmft_verbose >= 1) CALL rdmft_report_wall_time('OCC block', t_occ0)
          IF (rdmft_verbose >= 1) t_orb0 = rdmft_wall_time()
          CALL rdmft_alternating_orb_block(etot_now, orb_converged)
          IF (rdmft_verbose >= 1) CALL rdmft_report_wall_time('ORB block', t_orb0)
       ENDIF
       !
       occ_outer_dn_sum = SUM(ABS(rdmft_n - n_outer_save))
       CALL mp_sum(occ_outer_dn_sum, inter_pool_comm)
       !
       dE = ABS(etot_now - etot_prev)
       ! rdmft_weighted_sum / rdmft_magnetization_electrons reduce over
       ! inter_pool_comm and are therefore COLLECTIVE: evaluate them on
       ! every rank (the guarding conditions are uniform across ranks)
       ! BEFORE restricting the WRITE to ionode.  Calling them only on
       ! ionode desynchronises the inter-pool reductions between k-point
       ! pools, which corrupts every later energy / gradient reduction
       ! and ultimately aborts in ZPOTRF ("Cholesky failed in invchol").
       ne_now = rdmft_weighted_sum(rdmft_n, wk, nbnd, nks)
       ne_up_now = 0.0_DP
       ne_down_now = 0.0_DP
       mz_now = 0.0_DP
       do_spin_report = (nspin > 1 .AND. lsda)
       IF (do_spin_report) THEN
          ne_up_now = rdmft_weighted_sum_ispin(rdmft_n, wk, nbnd, nks, 1)
          ne_down_now = rdmft_weighted_sum_ispin(rdmft_n, wk, nbnd, nks, 2)
          mz_now = ne_up_now - ne_down_now
       ENDIF
       IF (rdmft_verbose >= 1 .AND. ionode) THEN
          IF (do_spin_report) THEN
             WRITE(stdout, '(5X,A,I4,A,F18.10,A,1PE10.2,0P,A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,1PE10.2,0P,A,L1,A,L1)') &
                  'RDMFT outer ', outer, '  E = ', etot_now, &
                  '  |dE| = ', dE, '  Ne = ', ne_now, &
                  '  N_up = ', ne_up_now, '  N_down = ', ne_down_now, &
                  '  Mz = ', mz_now, &
                  '  sum|dn|_out = ', occ_outer_dn_sum, &
                  '  occ_conv = ', occ_converged, '  orb_conv = ', orb_converged
          ELSE
             WRITE(stdout, '(5X,A,I4,A,F18.10,A,1PE10.2,0P,A,F12.6,A,1PE10.2,0P,A,L1,A,L1)') &
                  'RDMFT outer ', outer, '  E = ', etot_now, &
                  '  |dE| = ', dE, '  Ne = ', ne_now, &
                  '  sum|dn|_out = ', occ_outer_dn_sum, &
                  '  occ_conv = ', occ_converged, '  orb_conv = ', orb_converged
          ENDIF
       ENDIF
       IF (rdmft_verbose >= 1) &
            CALL rdmft_report_outer_wall_time(outer, t_outer0)
       IF (dE == 0.0_DP) THEN
          IF (ionode) WRITE(stdout, '(5X,A)') &
               'RDMFT alternating optimisation converged (|dE| = 0).'
          loop_converged = .TRUE.
          EXIT
       ELSE IF (outer > 1 .AND. dE < rdmft_energy_tol) THEN
          IF (ionode) WRITE(stdout, '(5X,A)') &
               'RDMFT alternating optimisation converged (|dE| below tolerance).'
          loop_converged = .TRUE.
          EXIT
       ELSE IF (occ_converged .AND. orb_converged) THEN
          IF (ionode) WRITE(stdout, '(5X,A)') &
               'RDMFT alternating optimisation converged (occ and orb blocks).'
          loop_converged = .TRUE.
          EXIT
       ENDIF
       ! Periodic checkpoint: save the current natural occupations
       ! (and persist the orbital buffer) every ``rdmft_save_every``
       ! outer cycles, so a crash / wall-time abort can be resumed
       ! later without losing the macro-iteration history.
       IF (rdmft_save_every > 0) THEN
          IF (MOD(outer, rdmft_save_every) == 0) &
               CALL rdmft_save_state(outer, etot_now)
       ENDIF
       etot_prev = etot_now
    ENDDO
    !
    ! Refresh the global energy decomposition so the final
    ! '!RDMFT_ETOTAL' line reflects the post-loop (n, evc) state and
    ! not whatever (possibly rejected) trial value was last left in
    ! rdmft_etot by the inner Armijo line search.
    !
    ! NOTE: with \texttt{rdmft\_block\_order='occ\_orb'} the loop ends
    ! on the ORB block; with \texttt{'orb\_occ'} it ends on OCC.  No
    ! extra post-loop occupation polish is performed -- bump
    ! \texttt{rdmft\_outer\_maxiter} for another outer cycle instead.
    !
    ! Warn loudly when the outer loop ran out of iterations WITHOUT
    ! reaching any convergence criterion.  In that case the reported
    ! total energy is NOT a converged stationary point of the
    ! functional: it still depends on the initial occupations / orbitals
    ! and on the (occ,orb) block ordering.  This is the usual cause of
    ! "different initial occupations give very different energies"
    ! reports -- e.g. running a single outer cycle
    ! (``rdmft_outer_maxiter = 1``) only does one partial alternating
    ! sweep.  Increase ``rdmft_outer_maxiter`` until ``|dE|`` drops
    ! below ``rdmft_energy_tol``.
    IF (.NOT. loop_converged .AND. rdmft_outer_maxiter >= 1) THEN
       IF (ionode) THEN
          WRITE(stdout, '(/,5X,A)') &
               'WARNING: RDMFT alternating optimisation did NOT converge.'
          WRITE(stdout, '(5X,A,I0,A)') &
               '         The outer loop hit rdmft_outer_maxiter = ', &
               rdmft_outer_maxiter, ' without |dE| < rdmft_energy_tol.'
          WRITE(stdout, '(5X,A,1PE10.2,A,1PE10.2,A)') &
               '         Last |dE| = ', dE, ' Ry (tol = ', &
               rdmft_energy_tol, ' Ry).'
          WRITE(stdout, '(5X,A)') &
               '         The reported total energy is NOT a converged'
          WRITE(stdout, '(5X,A)') &
               '         minimum and may still depend on the initial'
          WRITE(stdout, '(5X,A)') &
               '         occupations.  Increase rdmft_outer_maxiter.'
       ENDIF
    ENDIF
    CALL rdmft_total_energy(etot_now)
    !
    DEALLOCATE(n_outer_save)
    !
  END SUBROUTINE rdmft_run_alternating
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_run_joint(etot_now)
    !---------------------------------------------------------------
    !! Joint (product-manifold) RDMFT optimisation.
    !!
    !! The natural occupations are mapped onto unconstrained parameters
    !! ``p`` via the cosine-squared parameterisation
    !!   ``n_ik = cos^2(p_ik),   dn/dp = -sin(2 p)``
    !! so the box constraint ``0 <= n <= 1`` is automatically
    !! satisfied.  The augmented-Lagrangian penalty on
    !! ``sum_k w_k sum_i n_ik = N_e`` keeps the equality constraint
    !! satisfied (we keep \texttt{rdmft\_constraint = 'projected\_gradient'}
    !! semantics for the alternating path; for joint we always run the
    !! ALM penalty so the unified line search has a single merit
    !! function).
    !!
    !! At each outer iteration the solver:
    !!   1. assembles the full energy + gradients ``(dE/dn, dE/dC)``;
    !!   2. chain-rules ``dE/dn -> dE/dp`` through the parameterisation;
    !!   3. projects ``dE/dC`` onto the Stiefel tangent space at C;
    !!   4. packs the gradient into one vector, builds an SD or
    !!      Polak-Ribiere CG direction, re-projects the orbital block
    !!      onto the tangent space at C, and falls back to packed SD
    !!      if the joint directional derivative is non-negative;
    !!   5. runs **one** Armijo line search along the packed direction
    !!      (linear update for ``p``, Stiefel retraction for each
    !!      ``C^k``);
    !!   6. updates the augmented-Lagrangian multiplier.
    !
    USE io_global,        ONLY : stdout, ionode
    USE wvfct,            ONLY : nbnd, npwx, current_k, wg
    USE klist,            ONLY : nks, ngk, igk_k, xk, wk
    USE wavefunctions,    ONLY : evc
    USE noncollin_module, ONLY : npol, noncolin
    USE control_flags,    ONLY : gamma_only
    USE gvect,            ONLY : gstart
    USE io_files,         ONLY : nwordwfc, iunwfc
    USE buffers,          ONLY : get_buffer, save_buffer
    USE lsda_mod,         ONLY : lsda, current_spin, isk, nspin
    USE uspp,             ONLY : nkb, vkb, okvan
    USE uspp_init,        ONLY : init_us_2
    USE rdmft_module
    USE rdmft_lbfgs
    USE rdmft_linesearch, ONLY : rdmft_armijo_next_alpha, rdmft_occ_ls_alpha0, &
                                  rdmft_bb_init, rdmft_bb_record, rdmft_quad_ls_reset, &
                                  rdmft_bb_state, rdmft_quad_ls_history
    USE rdmft_xc,         ONLY : rdmft_xc_n_channels, rdmft_xc_channel
    USE rdmft_stiefel
    USE rdmft_energy,     ONLY : rdmft_apply_h_one_psi, &
                                  rdmft_compute_xc_channel, &
                                  rdmft_total_energy, &
                                  rdmft_grad_n
    USE rdmft_occupation, ONLY : rdmft_proximal_project_occ, rdmft_weighted_sum, &
                                  rdmft_weighted_sum_ispin
    USE rdmft_io,         ONLY : rdmft_save_state
    USE rdmft_occ_eps_log, ONLY : rdmft_occ_eps_log_set_outer, rdmft_occ_eps_log_record
    USE mp,               ONLY : mp_sum
    USE mp_pools,         ONLY : inter_pool_comm
    !
    REAL(DP), INTENT(INOUT) :: etot_now
    !
    REAL(DP), PARAMETER :: PI = 3.14159265358979323846_DP
    REAL(DP), PARAMETER :: LOGIT_EPS = 1.0e-12_DP
    !
    REAL(DP), ALLOCATABLE :: params(:,:), grad_p(:,:), dir_p(:,:)
    REAL(DP), ALLOCATABLE :: grad_n(:,:), prev_grad_p(:,:), prev_dir_p(:,:)
    COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), gradC(:,:), eta(:,:)
    COMPLEX(DP), ALLOCATABLE :: vxpsi(:,:,:), vxpsi_total(:,:,:), C_save(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: dir_C(:,:,:), prev_grad_C(:,:,:)
    REAL(DP), ALLOCATABLE :: vx_diag(:,:)
    REAL(DP), ALLOCATABLE :: pack_grad(:), pack_dir(:), pack_step(:), pack_diff(:)
    REAL(DP), ALLOCATABLE :: pack_iterate(:)
    REAL(DP) :: degspin, jac, n_clip, gnorm_p, gnorm_C, slope, alpha, etot_trial
    REAL(DP) :: coef, w, dw, e_t, etot_prev, dE, beta, num, den, occ_dn
    REAL(DP) :: ne_now, mz_now, t_outer0
    REAL(DP) :: ne_up_now, ne_down_now
    INTEGER :: outer, ik, ib, ls, it, nch, npw
    INTEGER :: ndim_pack, n_occ_params, orb_flat_per_k, ndim_orb_total
    LOGICAL :: line_ok, use_cg, use_lbfgs, do_spin_report
    LOGICAL :: loop_converged
    !
    TYPE(rdmft_lbfgs_state) :: joint_lbfgs
    TYPE(rdmft_bb_state)        :: joint_bb
    TYPE(rdmft_quad_ls_history) :: joint_qhist
    !
    ! Dispatch: the SPG product-manifold mode optimises occupations in
    ! native n-space (Euclidean proximal projection P_w + BMR spectral
    ! steplength + straight feasible chord) rather than through the
    ! cosine-squared reparameterisation used below, so it lives in its
    ! own driver.
    IF (TRIM(rdmft_joint_optimizer) == 'spg') THEN
       CALL rdmft_run_joint_spg(etot_now)
       RETURN
    ENDIF
    !
    loop_converged = .FALSE.
    use_cg    = (TRIM(rdmft_joint_optimizer) == 'cg')
    use_lbfgs = (TRIM(rdmft_joint_optimizer) == 'lbfgs') &
           .OR. (TRIM(rdmft_joint_optimizer) == 'bfgs')  &
           .OR. (TRIM(rdmft_joint_optimizer) == 'l-bfgs')
    degspin = 2.0_DP
    IF (noncolin) degspin = 1.0_DP
    nch = rdmft_xc_n_channels()
    !
    ALLOCATE(params(nbnd, nks), grad_p(nbnd, nks), dir_p(nbnd, nks))
    ALLOCATE(grad_n(nbnd, nks))
    ALLOCATE(prev_grad_p(nbnd, nks), prev_dir_p(nbnd, nks))
    ALLOCATE(hpsi(npwx*npol, nbnd), gradC(npwx*npol, nbnd))
    ALLOCATE(eta(npwx*npol, nbnd))
    ALLOCATE(vxpsi(npwx*npol, nbnd, nks), vx_diag(nbnd, nks))
    ! Channel-summed Vx|psi> over all k, recomputed once per outer
    ! iteration before the per-k Stiefel projection loop.  See the
    ! comment in step 3 below for the cost rationale.
    ALLOCATE(vxpsi_total(npwx*npol, nbnd, nks))
    ALLOCATE(C_save(npwx*npol, nbnd, nks), dir_C(npwx*npol, nbnd, nks))
    ALLOCATE(prev_grad_C(npwx*npol, nbnd, nks))
    ! Sizes for the packed (params, evc) view used by CG / L-BFGS.
    n_occ_params   = nbnd * nks
    orb_flat_per_k = 2 * npwx * npol * nbnd
    ndim_orb_total = orb_flat_per_k * nks
    ndim_pack      = n_occ_params + ndim_orb_total
    ALLOCATE(pack_grad(ndim_pack), pack_dir(ndim_pack), pack_step(ndim_pack), pack_diff(ndim_pack))
    ALLOCATE(pack_iterate(ndim_pack))
    IF (use_lbfgs) CALL rdmft_lbfgs_init(joint_lbfgs, ndim_pack, rdmft_lbfgs_memory)
    CALL rdmft_bb_init(joint_bb, ndim_pack, rdmft_bb_alpha_min, rdmft_bb_alpha_max)
    CALL rdmft_quad_ls_reset(joint_qhist)
    !
    pack_step = 0.0_DP
    pack_diff = 0.0_DP
    prev_grad_p = 0.0_DP
    prev_dir_p  = 0.0_DP
    prev_grad_C = (0.0_DP, 0.0_DP)
    !
    ! Initial parameters from the natural occupations: p = arccos(sqrt(n)).
    DO ik = 1, nks
       DO ib = 1, nbnd
          n_clip = MIN(MAX(LOGIT_EPS, rdmft_n(ib, ik)), 1.0_DP - LOGIT_EPS)
          params(ib, ik) = ACOS(SQRT(n_clip))
       ENDDO
    ENDDO
    !
    prev_grad_p = 0.0_DP
    prev_dir_p  = 0.0_DP
    etot_prev = etot_now
    !
    DO outer = 1, rdmft_outer_maxiter
       !
       CALL rdmft_occ_eps_log_set_outer(outer)
       IF (rdmft_verbose >= 1) t_outer0 = rdmft_wall_time()
       !
       ! 1. Refresh occupations from the parameters and rebuild density.
       DO ik = 1, nks
          DO ib = 1, nbnd
             rdmft_n(ib, ik) = COS(params(ib, ik))**2
          ENDDO
       ENDDO
       ! Project so the equality constraint is exactly satisfied (the
       ! cosine^2 mapping handles 0 <= n <= 1 alone, but the dual
       ! shift fixes the linear equality without touching the box.)
       CALL rdmft_proximal_project_occ(rdmft_n, wk, nbnd, nks)
       ! Mirror back to params consistent with the projected n.
       DO ik = 1, nks
          DO ib = 1, nbnd
             n_clip = MIN(MAX(LOGIT_EPS, rdmft_n(ib, ik)), 1.0_DP - LOGIT_EPS)
             params(ib, ik) = ACOS(SQRT(n_clip))
          ENDDO
       ENDDO
       !
       ! 2. Compute energy + dE/dn at the current point.
       CALL rdmft_total_energy(etot_now)
       CALL rdmft_grad_n(grad_n)
       CALL rdmft_occ_eps_log_record(outer, grad_n)
       !
       ! Chain rule: dE/dp = dE/dn * dn/dp = dE/dn * (-sin(2p)).
       DO ik = 1, nks
          DO ib = 1, nbnd
             jac = -SIN(2.0_DP * params(ib, ik))
             grad_p(ib, ik) = grad_n(ib, ik) * jac
          ENDDO
       ENDDO
       !
       ! 3. Build the orbital ambient gradient and project to Stiefel
       !    per k-point.  Save the gradients into dir_C for now (dir_C
       !    will become the search direction below).
       !
       !    Cost note: ``rdmft_compute_xc_channel`` does an all-k
       !    aceinit (``O(N_k^2)``) and was previously called inside the
       !    per-k loop, blowing the EXX cost up to ``O(N_k^3)``.  We now
       !    pre-compute the channel-summed ``Vx|psi>`` for every k once
       !    (``vxpsi_total``) and reuse it inside the per-k loop.
       vxpsi_total = (0.0_DP, 0.0_DP)
       DO it = 1, nch
          CALL rdmft_compute_xc_channel(it, e_t, vx_diag, vxpsi)
          ! Distinct (ib, ik) slices: thread-safe channel accumulation.
!$omp parallel do collapse(2) default(shared) private(ik, ib, coef, w, dw)
          DO ik = 1, nks
             DO ib = 1, nbnd
                CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
                vxpsi_total(:, ib, ik) = vxpsi_total(:, ib, ik) &
                                       + coef * w * vxpsi(:, ib, ik)
             ENDDO
          ENDDO
!$omp end parallel do
       ENDDO
       !
       DO ik = 1, nks
          npw = ngk(ik)
          current_k = ik
          IF (lsda) current_spin = isk(ik)
          IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
          CALL g2_kin(ik)
          IF (okvan .OR. nkb > 0) CALL init_us_2(npw, igk_k(1, ik), xk(:, ik), vkb)
          !
          CALL rdmft_apply_h_one_psi(npw, nbnd, evc, hpsi)
          DO ib = 1, nbnd
             gradC(:, ib) = wk(ik) * (rdmft_n(ib, ik) * hpsi(:, ib) &
                                      + vxpsi_total(:, ib, ik))
          ENDDO
          !
          ! Project to the Stiefel tangent space.
          IF (gamma_only) THEN
             CALL stiefel_project_tangent_gamma(evc, gradC, npw, nbnd, npwx*npol, gstart, eta)
          ELSE
             CALL stiefel_project_tangent_k(evc, gradC, npw, nbnd, npwx*npol, eta)
          ENDIF
          ! Save Riemannian gradient as dir_C and current evc as C_save.
          dir_C(:, :, ik) = eta(:, :)
          C_save(:, :, ik) = evc(:, :)
       ENDDO
       !
       ! 4. Packed search direction.  We support three optimizers on
       ! the packed (p, R-grad) vector:
       !   * sd     -> dir = -grad
       !   * cg     -> Polak-Ribiere with previous packed gradient
       !               (transport-by-reprojection on the orbital block)
       !   * lbfgs  -> Euclidean L-BFGS on the packed vector; the
       !               orbital block of the returned direction is
       !               re-projected onto the current Stiefel tangent.
       IF (use_lbfgs .OR. use_cg) THEN
          ! Pack the gradient: (grad_p, R-grad blocks).  dir_C is
          ! holding the projected Riemannian gradient at this point.
          CALL pack_joint(grad_p, dir_C, n_occ_params, orb_flat_per_k, nks, &
                          npwx*npol, nbnd, pack_grad)
          IF (use_lbfgs) THEN
             ! Update history with the (step, grad-diff) pair from the
             ! previous accepted iteration.  The very first iteration
             ! (outer == 1) has pack_step = pack_diff = 0 and the
             ! L-BFGS update is skipped by the curvature safeguard.
             IF (outer > 1) THEN
                pack_diff = pack_grad - pack_diff  ! grad_new - grad_old
                CALL rdmft_lbfgs_update(joint_lbfgs, pack_step, pack_diff)
             ENDIF
             CALL rdmft_lbfgs_direction(joint_lbfgs, pack_grad, pack_dir)
             ! Save current pack_grad for the next outer iteration's
             ! grad-difference computation.
             pack_diff = pack_grad
          ELSE
             ! Polak-Ribiere on the packed Euclidean vector.
             IF (outer > 1) THEN
                CALL pack_joint(prev_grad_p, prev_grad_C, n_occ_params, &
                                orb_flat_per_k, nks, npwx*npol, nbnd, pack_step)
                den = SUM(pack_step * pack_step)
                num = SUM(pack_grad * (pack_grad - pack_step))
                ! The packed vector is k-pool distributed; reduce the CG
                ! inner products over the pools so beta is identical on
                ! every pool (a per-pool direction would desynchronise the
                ! shared line search).
                CALL mp_sum(den, inter_pool_comm)
                CALL mp_sum(num, inter_pool_comm)
                IF (den > 1.0e-30_DP) THEN
                   beta = MAX(0.0_DP, num / den)
                ELSE
                   beta = 0.0_DP
                ENDIF
                pack_dir = -pack_grad + beta * (-pack_step)
             ELSE
                pack_dir = -pack_grad
             ENDIF
          ENDIF
          ! Unpack: dir_p (linear), dir_C (re-project to tangent at C).
          CALL unpack_joint(pack_dir, n_occ_params, orb_flat_per_k, nks, &
                            npwx*npol, nbnd, dir_p, dir_C)
          DO ik = 1, nks
             npw = ngk(ik)
             current_k = ik
             IF (lsda) current_spin = isk(ik)
             IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
             eta = dir_C(:, :, ik)
             IF (gamma_only) THEN
                CALL stiefel_project_tangent_gamma(evc, eta, npw, nbnd, &
                     npwx*npol, gstart, dir_C(:, :, ik))
             ELSE
                CALL stiefel_project_tangent_k(evc, eta, npw, nbnd, &
                     npwx*npol, dir_C(:, :, ik))
             ENDIF
          ENDDO
          ! Save grad for next-iter CG/L-BFGS history and convert eta
          ! (currently the Riemannian gradient) into prev_grad_C.
          ! At this point dir_C is the search direction; we need the
          ! Riemannian gradient saved separately.  We re-fetch by
          ! reading back from pack_grad.
          CALL unpack_joint(pack_grad, n_occ_params, orb_flat_per_k, nks, &
                            npwx*npol, nbnd, prev_grad_p, prev_grad_C)
       ELSE
          ! Steepest descent: dir_C(:, :, ik) currently holds the
          ! Riemannian gradient G_R^k from step 3.  Save it into
          ! prev_grad_C BEFORE flipping the sign so the joint
          ! directional-derivative computation below has a uniform
          ! representation of G_R regardless of optimiser branch.
          prev_grad_C = dir_C
          dir_p = -grad_p
          dir_C = -dir_C
          prev_grad_p = grad_p
       ENDIF
       !
       ! Joint directional derivative dd = <grad_p, dir_p>
       !                              + sum_k <G_R^k, dir_C^k>_S.
       !
       ! We use prev_grad_C as the representation of G_R^k (set above
       ! in every branch).  The previous implementation evaluated
       ! ``<-dir_C, dir_C>_S = -|dir_C|_S^2`` instead, which only
       ! coincides with ``<G_R, dir_C>_S`` in the SD branch where
       ! ``dir_C = -G_R``.  For CG / L-BFGS, ``dir_C`` is an
       ! arbitrary search direction (no longer parallel to G_R), so
       ! the old expression returned the wrong slope: the descent
       ! check ``slope >= 0`` then never triggered (-|d|^2 is
       ! trivially negative) and the Armijo line search was driven
       ! by a slope of the wrong magnitude / sign.
       slope = SUM(grad_p * dir_p)
       DO ik = 1, nks
          IF (gamma_only) THEN
             slope = slope + stiefel_inner_product_gamma(prev_grad_C(:, :, ik), &
                                                         dir_C(:, :, ik), &
                                                         ngk(ik), nbnd, npwx*npol, gstart)
          ELSE
             slope = slope + stiefel_inner_product_k(prev_grad_C(:, :, ik), &
                                                     dir_C(:, :, ik), &
                                                     ngk(ik), nbnd, npwx*npol)
          ENDIF
       ENDDO
       ! grad_p / dir_p and the per-k inner products only cover this
       ! pool's k-points; reduce the joint slope over the k-point pools.
       CALL mp_sum(slope, inter_pool_comm)
       IF (slope >= 0.0_DP) THEN
          ! Direction is not a descent on the joint manifold; fall
          ! back to packed steepest descent.  Note that dir_C must
          ! be reset to ``-G_R`` (= -prev_grad_C) -- the previous
          ! implementation only updated dir_p and left dir_C at the
          ! non-descent CG / L-BFGS direction, so the line search
          ! retracted along a direction whose true slope was
          ! positive and the energy increased.
          dir_p = -grad_p
          DO ik = 1, nks
             dir_C(:, :, ik) = -prev_grad_C(:, :, ik)
          ENDDO
          slope = -SUM(grad_p * grad_p)
          DO ik = 1, nks
             IF (gamma_only) THEN
                slope = slope - stiefel_inner_product_gamma(prev_grad_C(:, :, ik), &
                                                            prev_grad_C(:, :, ik), &
                                                            ngk(ik), nbnd, npwx*npol, gstart)
             ELSE
                slope = slope - stiefel_inner_product_k(prev_grad_C(:, :, ik), &
                                                        prev_grad_C(:, :, ik), &
                                                        ngk(ik), nbnd, npwx*npol)
             ENDIF
          ENDDO
          CALL mp_sum(slope, inter_pool_comm)
       ENDIF
       prev_grad_p = grad_p
       prev_dir_p  = dir_p
       !
       ! 5. Single Armijo line search along the packed direction.
       !    Initial steplength from ``rdmft_orb_ls_init_step`` (default
       !    ``bb``): Barzilai-Borwein on the packed (p, C) iterate,
       !    falling back to ``rdmft_orb_ls_stepsize`` on the first step.
       CALL pack_joint(params, C_save, n_occ_params, orb_flat_per_k, nks, &
                       npwx*npol, nbnd, pack_iterate)
       alpha = rdmft_occ_ls_alpha0(rdmft_orb_ls_init_step, pack_iterate, pack_grad, &
            rdmft_orb_ls_stepsize, joint_bb, joint_qhist, &
            rdmft_bb_alpha_min, rdmft_bb_alpha_max)
       IF (alpha <= 0.0_DP) alpha = rdmft_orb_ls_stepsize
       CALL rdmft_bb_record(joint_bb, pack_iterate, pack_grad)
       line_ok = .FALSE.
       DO ls = 1, rdmft_line_search_max_iter
          ! Update parameters linearly.
          DO ik = 1, nks
             DO ib = 1, nbnd
                rdmft_n(ib, ik) = COS(params(ib, ik) + alpha * dir_p(ib, ik))**2
             ENDDO
          ENDDO
          CALL rdmft_proximal_project_occ(rdmft_n, wk, nbnd, nks)
          ! Update orbitals via Stiefel retraction at each k-point.
          DO ik = 1, nks
             npw = ngk(ik)
             current_k = ik
             IF (lsda) current_spin = isk(ik)
             IF (gamma_only) THEN
                CALL stiefel_retract_gamma(C_save(:, :, ik), dir_C(:, :, ik), &
                                           alpha, npw, nbnd, npwx*npol, gstart, evc)
             ELSE
                CALL stiefel_retract_k(C_save(:, :, ik), dir_C(:, :, ik), &
                                        alpha, npw, nbnd, npwx*npol, evc)
             ENDIF
             IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
          ENDDO
          CALL rdmft_mark_orbitals_changed()
          CALL rdmft_total_energy(etot_trial)
          IF (etot_trial <= etot_now + rdmft_line_search_c1 * alpha * slope) THEN
             line_ok = .TRUE.
             etot_now = etot_trial
             EXIT
          ENDIF
          IF (rdmft_line_search_polynomial) THEN
             CALL rdmft_armijo_next_alpha(alpha, etot_now, slope, etot_trial, &
                  rdmft_line_search_rho, alpha)
          ELSE
             alpha = alpha * rdmft_line_search_rho
          ENDIF
       ENDDO
       !
       IF (.NOT. line_ok) THEN
          ! Fallback: revert orbitals; take a tiny SD step on params only.
          DO ik = 1, nks
             evc(:, :) = C_save(:, :, ik)
             IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
          ENDDO
          CALL rdmft_mark_orbitals_changed()
          DO ik = 1, nks
             DO ib = 1, nbnd
                rdmft_n(ib, ik) = COS(params(ib, ik) - 1.0e-3_DP * grad_p(ib, ik))**2
             ENDDO
          ENDDO
          CALL rdmft_proximal_project_occ(rdmft_n, wk, nbnd, nks)
          CALL rdmft_total_energy(etot_now)
          IF (use_lbfgs) CALL rdmft_lbfgs_reset(joint_lbfgs)
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(5X,A,I4)') &
               'RDMFT joint: line search failed; tiny SD fallback at outer ', outer
       ELSE
          ! Update params from the accepted step.
          params = params + alpha * dir_p
          ! Save the accepted packed step so the next outer iteration
          ! can build the L-BFGS curvature pair.
          IF (use_lbfgs) pack_step = alpha * pack_dir
       ENDIF
       !
       occ_dn = SUM(ABS(rdmft_n - COS(params - alpha * dir_p)**2))
       CALL mp_sum(occ_dn, inter_pool_comm)
       dE = ABS(etot_now - etot_prev)
       ! Collective inter-pool reductions must run on every rank; see the
       ! matching comment in rdmft_run_alternating.
       ne_now = rdmft_weighted_sum(rdmft_n, wk, nbnd, nks)
       ne_up_now = 0.0_DP
       ne_down_now = 0.0_DP
       mz_now = 0.0_DP
       do_spin_report = (nspin > 1 .AND. lsda)
       IF (do_spin_report) THEN
          ne_up_now = rdmft_weighted_sum_ispin(rdmft_n, wk, nbnd, nks, 1)
          ne_down_now = rdmft_weighted_sum_ispin(rdmft_n, wk, nbnd, nks, 2)
          mz_now = ne_up_now - ne_down_now
       ENDIF
       IF (rdmft_verbose >= 1 .AND. ionode) THEN
          IF (do_spin_report) THEN
             WRITE(stdout, '(5X,A,I4,A,F18.10,A,1PE10.2,0P,A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,1PE10.2)') &
                  'RDMFT joint outer ', outer, '  E = ', etot_now, &
                  '  |dE| = ', dE, '  Ne = ', ne_now, &
                  '  N_up = ', ne_up_now, '  N_down = ', ne_down_now, &
                  '  Mz = ', mz_now, &
                  '  alpha = ', alpha
          ELSE
             WRITE(stdout, '(5X,A,I4,A,F18.10,A,1PE10.2,0P,A,F12.6,A,1PE10.2)') &
                  'RDMFT joint outer ', outer, '  E = ', etot_now, &
                  '  |dE| = ', dE, '  Ne = ', ne_now, &
                  '  alpha = ', alpha
          ENDIF
       ENDIF
       IF (rdmft_verbose >= 1) &
            CALL rdmft_report_outer_wall_time(outer, t_outer0)
       IF (dE == 0.0_DP) THEN
          IF (ionode) WRITE(stdout, '(5X,A)') &
               'RDMFT joint optimisation converged (|dE| = 0).'
          loop_converged = .TRUE.
          EXIT
       ELSE IF (outer > 1 .AND. dE < rdmft_energy_tol) THEN
          IF (ionode) WRITE(stdout, '(5X,A)') &
               'RDMFT joint optimisation converged.'
          loop_converged = .TRUE.
          EXIT
       ENDIF
       ! Periodic restart checkpoint (see rdmft_run_alternating).
       IF (rdmft_save_every > 0) THEN
          IF (MOD(outer, rdmft_save_every) == 0) &
               CALL rdmft_save_state(outer, etot_now)
       ENDIF
       etot_prev = etot_now
    ENDDO
    !
    ! Warn when the joint outer loop ran out of iterations without
    ! reaching ``|dE| < rdmft_energy_tol`` (see the analogous note in
    ! rdmft_run_alternating): the reported total energy is then not a
    ! converged minimum and can still depend on the initial guess.
    IF (.NOT. loop_converged .AND. rdmft_outer_maxiter >= 1) THEN
       IF (ionode) THEN
          WRITE(stdout, '(/,5X,A)') &
               'WARNING: RDMFT joint optimisation did NOT converge.'
          WRITE(stdout, '(5X,A,I0,A)') &
               '         The outer loop hit rdmft_outer_maxiter = ', &
               rdmft_outer_maxiter, ' without |dE| < rdmft_energy_tol.'
          WRITE(stdout, '(5X,A,1PE10.2,A,1PE10.2,A)') &
               '         Last |dE| = ', dE, ' Ry (tol = ', &
               rdmft_energy_tol, ' Ry).'
          WRITE(stdout, '(5X,A)') &
               '         The reported total energy is NOT a converged'
          WRITE(stdout, '(5X,A)') &
               '         minimum and may still depend on the initial'
          WRITE(stdout, '(5X,A)') &
               '         occupations.  Increase rdmft_outer_maxiter.'
       ENDIF
    ENDIF
    CALL rdmft_total_energy(etot_now)
    !
    DEALLOCATE(params, grad_p, dir_p, grad_n, prev_grad_p, prev_dir_p)
    DEALLOCATE(hpsi, gradC, eta, vxpsi, vxpsi_total, vx_diag, C_save, dir_C, prev_grad_C)
    DEALLOCATE(pack_grad, pack_dir, pack_step, pack_diff, pack_iterate)
    IF (use_lbfgs) CALL rdmft_lbfgs_finalize(joint_lbfgs)
    !
  END SUBROUTINE rdmft_run_joint
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_run_joint_spg(etot_now)
    !---------------------------------------------------------------
    !! Joint (product-manifold) RDMFT optimisation, SPG variant
    !! (``rdmft_solver_strategy = 'joint'`` with
    !! ``rdmft_joint_optimizer = 'spg'``).
    !!
    !! Combines the canonical SPG2 projected-gradient occupation step
    !! with the Stiefel CG orbital step and drives both through a single
    !! Armijo line search on the product manifold
    !!
    !!   \(\mathcal{F} \times \prod_{\mathbf{k}} \mathrm{St}(N_b, N_{\rm pw}; I)\),
    !!
    !! where \(\mathcal{F} = \{0 \le n \le 1,\ \sum_{i\mathbf{k}} w_{\mathbf{k}}
    !! n_{i\mathbf{k}} = N_e\}\) is the occupation polytope (Section 5).
    !!
    !! Per outer iteration at \((\mathbf{n}_0, C_0)\):
    !!   1. project \(\mathbf{n}_0\) onto \(\mathcal{F}\) with \(P_w\); assemble
    !!      \(E\), the occupation gradient \(\mathbf{g}_n = \partial E/\partial n\)
    !!      and the Stiefel Riemannian gradient \(G_C^{\mathbf{k}}\) per k;
    !!   2. occupation direction (SPG2 straight chord, native n-space):
    !!      \(\mathbf{d}_n = P_w(\mathbf{n}_0 - \alpha_n \mathbf{g}_n) - \mathbf{n}_0\)
    !!      with \(\alpha_n\) the BMR spectral steplength
    !!      (:func:`rdmft_spg_spectral_alpha`);
    !!   3. orbital direction (Stiefel Polak--Ribiere CG with
    !!      transport-by-reprojection and steepest-descent fallback);
    !!   4. joint slope \(\phi'(0) = \langle \mathbf{g}_n, \mathbf{d}_n\rangle
    !!      + \sum_{\mathbf{k}} \langle G_C^{\mathbf{k}}, \mathbf{d}_C^{\mathbf{k}}
    !!      \rangle_S\) (both terms are descent contributions);
    !!   5. **one** Armijo line search coupling the feasible chord
    !!      \(\mathbf{n}(\lambda) = \mathbf{n}_0 + \lambda \mathbf{d}_n\)
    !!      (projection-free for \(\lambda \in [0,1]\) by convexity) with the
    !!      Stiefel retraction \(C^{\mathbf{k}}(\lambda) =
    !!      R_{C_0^{\mathbf{k}}}(\lambda \mathbf{d}_C^{\mathbf{k}})\), first
    !!      trial \(\lambda = 1\) (BMR SPG2 convention).
    !!
    !! Convergence uses the unified stationarity map
    !! \(\max(\|g_1\|_\infty, \|G_C\|) \le\) tol (occupation SPG2 stopping map
    !! + Stiefel gradient norm) or \(|dE| < \) ``rdmft_energy_tol``.
    !
    USE io_global,        ONLY : stdout, ionode
    USE wvfct,            ONLY : nbnd, npwx, current_k
    USE klist,            ONLY : nks, ngk, igk_k, xk, wk
    USE wavefunctions,    ONLY : evc
    USE noncollin_module, ONLY : npol
    USE control_flags,    ONLY : gamma_only
    USE gvect,            ONLY : gstart
    USE io_files,         ONLY : nwordwfc, iunwfc
    USE buffers,          ONLY : get_buffer, save_buffer
    USE lsda_mod,         ONLY : lsda, current_spin, isk, nspin
    USE uspp,             ONLY : nkb, vkb, okvan
    USE uspp_init,        ONLY : init_us_2
    USE rdmft_module
    USE rdmft_linesearch, ONLY : rdmft_bb_state, rdmft_bb_init, rdmft_bb_reset, &
                                  rdmft_bb_record, rdmft_armijo_next_alpha
    USE rdmft_xc,         ONLY : rdmft_xc_n_channels, rdmft_xc_channel
    USE rdmft_stiefel
    USE rdmft_spg,        ONLY : rdmft_spg_spectral_alpha
    USE rdmft_energy,     ONLY : rdmft_apply_h_one_psi, rdmft_compute_xc_channel, &
                                  rdmft_total_energy, rdmft_grad_n
    USE rdmft_occupation, ONLY : rdmft_proximal_project_occ, rdmft_pg_map_grad_inf, &
                                  rdmft_weighted_sum, rdmft_weighted_sum_ispin
    USE rdmft_io,         ONLY : rdmft_save_state
    USE rdmft_occ_eps_log, ONLY : rdmft_occ_eps_log_set_outer, rdmft_occ_eps_log_record
    USE mp,               ONLY : mp_sum
    USE mp_pools,         ONLY : inter_pool_comm
    !
    REAL(DP), INTENT(INOUT) :: etot_now
    !
    REAL(DP), ALLOCATABLE :: grad_n(:,:), d_n(:,:), n_save(:,:)
    REAL(DP), ALLOCATABLE :: occ_flat(:), grad_flat(:)
    COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), gradC(:,:), eta(:,:)
    COMPLEX(DP), ALLOCATABLE :: vxpsi(:,:,:), vxpsi_total(:,:,:)
    REAL(DP), ALLOCATABLE :: vx_diag(:,:)
    COMPLEX(DP), ALLOCATABLE :: C_save(:,:,:), dir_C(:,:,:), G_C(:,:,:), prev_G_C(:,:,:), prev_dir_C(:,:,:)
    REAL(DP) :: coef, w, dw, e_t, etot_prev, dE, etot_trial
    REAL(DP) :: alpha_n, alpha_init, g1_inf, gnorm_c2, gnorm_c, slope_n, slope_c
    REAL(DP) :: phi0, alpha, beta, num, den, occ_dn
    REAL(DP) :: ne_now, ne_up_now, ne_down_now, mz_now, t_outer0
    INTEGER :: outer, ik, ib, ls, it, nch, npw, ndim
    LOGICAL :: line_ok, do_spin_report, loop_converged, occ_conv, orb_conv
    TYPE(rdmft_bb_state) :: occ_bb
    !
    loop_converged = .FALSE.
    nch = rdmft_xc_n_channels()
    ndim = nbnd * nks
    !
    ALLOCATE(grad_n(nbnd, nks), d_n(nbnd, nks), n_save(nbnd, nks))
    ALLOCATE(occ_flat(ndim), grad_flat(ndim))
    ALLOCATE(hpsi(npwx*npol, nbnd), gradC(npwx*npol, nbnd), eta(npwx*npol, nbnd))
    ALLOCATE(vxpsi(npwx*npol, nbnd, nks), vxpsi_total(npwx*npol, nbnd, nks))
    ALLOCATE(vx_diag(nbnd, nks))
    ALLOCATE(C_save(npwx*npol, nbnd, nks), dir_C(npwx*npol, nbnd, nks), G_C(npwx*npol, nbnd, nks))
    ALLOCATE(prev_G_C(npwx*npol, nbnd, nks), prev_dir_C(npwx*npol, nbnd, nks))
    !
    CALL rdmft_bb_init(occ_bb, ndim, rdmft_spg_alpha_min, rdmft_spg_alpha_max)
    prev_G_C   = (0.0_DP, 0.0_DP)
    prev_dir_C = (0.0_DP, 0.0_DP)
    etot_prev  = etot_now
    !
    IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(5X,A)') &
         'RDMFT joint (SPG2 x Stiefel CG) product-manifold optimisation'
    !
    DO outer = 1, rdmft_outer_maxiter
       !
       CALL rdmft_occ_eps_log_set_outer(outer)
       IF (rdmft_verbose >= 1) t_outer0 = rdmft_wall_time()
       !
       ! 1. Restore feasibility and assemble energy + gradients.
       CALL rdmft_proximal_project_occ(rdmft_n, wk, nbnd, nks)
       CALL rdmft_total_energy(etot_now)
       CALL rdmft_grad_n(grad_n)
       CALL rdmft_occ_eps_log_record(outer, grad_n)
       n_save = rdmft_n
       !
       ! Orbital ambient gradient -> Stiefel Riemannian gradient per k.
       ! Hoist the channel-summed Vx|psi> for all k once (see the cost
       ! note in rdmft_run_joint).
       vxpsi_total = (0.0_DP, 0.0_DP)
       DO it = 1, nch
          CALL rdmft_compute_xc_channel(it, e_t, vx_diag, vxpsi)
!$omp parallel do collapse(2) default(shared) private(ik, ib, coef, w, dw)
          DO ik = 1, nks
             DO ib = 1, nbnd
                CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
                vxpsi_total(:, ib, ik) = vxpsi_total(:, ib, ik) + coef * w * vxpsi(:, ib, ik)
             ENDDO
          ENDDO
!$omp end parallel do
       ENDDO
       !
       gnorm_c2 = 0.0_DP
       DO ik = 1, nks
          npw = ngk(ik)
          current_k = ik
          IF (lsda) current_spin = isk(ik)
          IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
          CALL g2_kin(ik)
          IF (okvan .OR. nkb > 0) CALL init_us_2(npw, igk_k(1, ik), xk(:, ik), vkb)
          CALL rdmft_apply_h_one_psi(npw, nbnd, evc, hpsi)
          DO ib = 1, nbnd
             gradC(:, ib) = wk(ik) * (rdmft_n(ib, ik) * hpsi(:, ib) + vxpsi_total(:, ib, ik))
          ENDDO
          IF (gamma_only) THEN
             CALL stiefel_project_tangent_gamma(evc, gradC, npw, nbnd, npwx*npol, gstart, eta)
          ELSE
             CALL stiefel_project_tangent_k(evc, gradC, npw, nbnd, npwx*npol, eta)
          ENDIF
          G_C(:, :, ik) = eta(:, :)
          C_save(:, :, ik) = evc(:, :)
          IF (gamma_only) THEN
             gnorm_c2 = gnorm_c2 + stiefel_inner_product_gamma(eta, eta, npw, nbnd, npwx*npol, gstart)
          ELSE
             gnorm_c2 = gnorm_c2 + stiefel_inner_product_k(eta, eta, npw, nbnd, npwx*npol)
          ENDIF
       ENDDO
       ! Per-k inner products only cover this pool's k-points.
       CALL mp_sum(gnorm_c2, inter_pool_comm)
       gnorm_c = SQRT(MAX(0.0_DP, gnorm_c2))
       !
       ! SPG2 occupation stopping map ||g_1||_inf (step t = 1).
       CALL rdmft_pg_map_grad_inf(rdmft_n, grad_n, wk, nbnd, nks, 1.0_DP, g1_inf)
       !
       ! Unified stationarity check (skip on the very first iteration so
       ! at least one step is taken).
       occ_conv = (rdmft_occ_grad_tol <= 0.0_DP) .OR. (g1_inf <= rdmft_occ_grad_tol)
       orb_conv = (rdmft_orb_grad_tol <= 0.0_DP) .OR. (gnorm_c <= rdmft_orb_grad_tol)
       IF (outer > 1 .AND. occ_conv .AND. orb_conv) THEN
          IF (ionode) WRITE(stdout, '(5X,A,1PE10.2,A,1PE10.2,A)') &
               'RDMFT joint (SPG) converged: ||g_1||_inf=', g1_inf, &
               ' ||G_C||=', gnorm_c, ' (stationarity).'
          loop_converged = .TRUE.
          EXIT
       ENDIF
       !
       ! 2. Occupation direction: SPG2 straight feasible chord.
       occ_flat  = RESHAPE(rdmft_n, [ndim])
       grad_flat = RESHAPE(grad_n,  [ndim])
       IF (g1_inf > 1.0e-30_DP) THEN
          alpha_init = 1.0_DP / g1_inf
       ELSE
          alpha_init = rdmft_spg_alpha_max
       ENDIF
       alpha_init = MIN(rdmft_spg_alpha_max, MAX(rdmft_spg_alpha_min, alpha_init))
       alpha_n = rdmft_spg_spectral_alpha(occ_bb, occ_flat, grad_flat, &
            rdmft_spg_alpha_min, rdmft_spg_alpha_max, alpha_init)
       IF (alpha_n <= 0.0_DP) alpha_n = alpha_init
       ! Record the current (n, g) pair for the next spectral steplength.
       CALL rdmft_bb_record(occ_bb, occ_flat, grad_flat)
       ! d_n = P_w(n_0 - alpha_n * g_n) - n_0.
       DO ik = 1, nks
          DO ib = 1, nbnd
             d_n(ib, ik) = n_save(ib, ik) - alpha_n * grad_n(ib, ik)
          ENDDO
       ENDDO
       CALL rdmft_proximal_project_occ(d_n, wk, nbnd, nks)
       d_n = d_n - n_save
       slope_n = SUM(grad_n * d_n)
       CALL mp_sum(slope_n, inter_pool_comm)
       !
       ! 3. Orbital direction: Stiefel Polak-Ribiere CG (transport by
       !    reprojection onto the current tangent), SD on the first
       !    iteration or when the CG direction is not a descent.
       IF (outer > 1) THEN
          ! Reproject the stored previous gradient/direction onto the
          ! tangent space at the current C (vector transport), then form
          ! Polak-Ribiere beta with the transported quantities.
          num = 0.0_DP
          den = 0.0_DP
          DO ik = 1, nks
             npw = ngk(ik)
             IF (gamma_only) THEN
                CALL stiefel_project_tangent_gamma(C_save(:,:,ik), prev_G_C(:,:,ik), &
                     npw, nbnd, npwx*npol, gstart, eta)
             ELSE
                CALL stiefel_project_tangent_k(C_save(:,:,ik), prev_G_C(:,:,ik), &
                     npw, nbnd, npwx*npol, eta)
             ENDIF
             prev_G_C(:, :, ik) = eta(:, :)
             IF (gamma_only) THEN
                CALL stiefel_project_tangent_gamma(C_save(:,:,ik), prev_dir_C(:,:,ik), &
                     npw, nbnd, npwx*npol, gstart, eta)
             ELSE
                CALL stiefel_project_tangent_k(C_save(:,:,ik), prev_dir_C(:,:,ik), &
                     npw, nbnd, npwx*npol, eta)
             ENDIF
             prev_dir_C(:, :, ik) = eta(:, :)
             IF (gamma_only) THEN
                num = num + stiefel_inner_product_gamma(G_C(:,:,ik), &
                     G_C(:,:,ik) - prev_G_C(:,:,ik), npw, nbnd, npwx*npol, gstart)
                den = den + stiefel_inner_product_gamma(prev_G_C(:,:,ik), &
                     prev_G_C(:,:,ik), npw, nbnd, npwx*npol, gstart)
             ELSE
                num = num + stiefel_inner_product_k(G_C(:,:,ik), &
                     G_C(:,:,ik) - prev_G_C(:,:,ik), npw, nbnd, npwx*npol)
                den = den + stiefel_inner_product_k(prev_G_C(:,:,ik), &
                     prev_G_C(:,:,ik), npw, nbnd, npwx*npol)
             ENDIF
          ENDDO
          CALL mp_sum(num, inter_pool_comm)
          CALL mp_sum(den, inter_pool_comm)
          IF (den > 1.0e-30_DP) THEN
             beta = MAX(0.0_DP, num / den)   ! Polak-Ribiere+ (Hager safeguard)
          ELSE
             beta = 0.0_DP
          ENDIF
          DO ik = 1, nks
             dir_C(:, :, ik) = -G_C(:, :, ik) + beta * prev_dir_C(:, :, ik)
          ENDDO
       ELSE
          DO ik = 1, nks
             dir_C(:, :, ik) = -G_C(:, :, ik)
          ENDDO
       ENDIF
       ! Re-project the CG direction onto the tangent space at C and
       ! compute the orbital slope <G_C, dir_C>_S.
       slope_c = 0.0_DP
       DO ik = 1, nks
          npw = ngk(ik)
          IF (gamma_only) THEN
             CALL stiefel_project_tangent_gamma(C_save(:,:,ik), dir_C(:,:,ik), &
                  npw, nbnd, npwx*npol, gstart, eta)
          ELSE
             CALL stiefel_project_tangent_k(C_save(:,:,ik), dir_C(:,:,ik), &
                  npw, nbnd, npwx*npol, eta)
          ENDIF
          dir_C(:, :, ik) = eta(:, :)
          IF (gamma_only) THEN
             slope_c = slope_c + stiefel_inner_product_gamma(G_C(:,:,ik), &
                  dir_C(:,:,ik), npw, nbnd, npwx*npol, gstart)
          ELSE
             slope_c = slope_c + stiefel_inner_product_k(G_C(:,:,ik), &
                  dir_C(:,:,ik), npw, nbnd, npwx*npol)
          ENDIF
       ENDDO
       CALL mp_sum(slope_c, inter_pool_comm)
       ! Descent safeguard: fall back to steepest descent on the orbital
       ! block if the CG direction is not downhill.
       IF (slope_c >= 0.0_DP) THEN
          DO ik = 1, nks
             dir_C(:, :, ik) = -G_C(:, :, ik)
          ENDDO
          slope_c = -gnorm_c2
       ENDIF
       ! Save the (gradient, direction) for the next CG iteration.
       prev_G_C   = G_C
       prev_dir_C = dir_C
       !
       phi0 = slope_n + slope_c
       !
       ! 4. Single Armijo line search on the product manifold, lambda0 = 1.
       alpha = 1.0_DP
       line_ok = .FALSE.
       DO ls = 1, rdmft_line_search_max_iter
          ! Occupation chord (feasible for lambda in [0,1]).
          DO ik = 1, nks
             DO ib = 1, nbnd
                rdmft_n(ib, ik) = n_save(ib, ik) + alpha * d_n(ib, ik)
             ENDDO
          ENDDO
          ! Orbital Stiefel retraction at every k.
          DO ik = 1, nks
             npw = ngk(ik)
             current_k = ik
             IF (lsda) current_spin = isk(ik)
             IF (gamma_only) THEN
                CALL stiefel_retract_gamma(C_save(:,:,ik), dir_C(:,:,ik), &
                     alpha, npw, nbnd, npwx*npol, gstart, evc)
             ELSE
                CALL stiefel_retract_k(C_save(:,:,ik), dir_C(:,:,ik), &
                     alpha, npw, nbnd, npwx*npol, evc)
             ENDIF
             IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
          ENDDO
          CALL rdmft_mark_orbitals_changed()
          CALL rdmft_total_energy(etot_trial)
          IF (etot_trial <= etot_now + rdmft_spg_gamma * alpha * phi0) THEN
             line_ok = .TRUE.
             etot_now = etot_trial
             EXIT
          ENDIF
          IF (rdmft_line_search_polynomial) THEN
             CALL rdmft_armijo_next_alpha(alpha, etot_now, phi0, etot_trial, &
                  rdmft_line_search_rho, alpha)
          ELSE
             alpha = alpha * rdmft_line_search_rho
          ENDIF
       ENDDO
       !
       IF (.NOT. line_ok) THEN
          ! Revert orbitals; take a tiny projected SD step on occupations
          ! only, and restart the orbital CG history.
          DO ik = 1, nks
             evc(:, :) = C_save(:, :, ik)
             IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
          ENDDO
          CALL rdmft_mark_orbitals_changed()
          DO ik = 1, nks
             DO ib = 1, nbnd
                rdmft_n(ib, ik) = n_save(ib, ik) - 1.0e-3_DP * grad_n(ib, ik)
             ENDDO
          ENDDO
          CALL rdmft_proximal_project_occ(rdmft_n, wk, nbnd, nks)
          CALL rdmft_total_energy(etot_now)
          CALL rdmft_bb_reset(occ_bb)
          prev_dir_C = (0.0_DP, 0.0_DP)
          prev_G_C   = (0.0_DP, 0.0_DP)
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(5X,A,I4)') &
               'RDMFT joint (SPG): line search failed; tiny SD fallback at outer ', outer
       ENDIF
       !
       occ_dn = SUM(ABS(rdmft_n - n_save))
       CALL mp_sum(occ_dn, inter_pool_comm)
       dE = ABS(etot_now - etot_prev)
       ne_now = rdmft_weighted_sum(rdmft_n, wk, nbnd, nks)
       ne_up_now = 0.0_DP
       ne_down_now = 0.0_DP
       mz_now = 0.0_DP
       do_spin_report = (nspin > 1 .AND. lsda)
       IF (do_spin_report) THEN
          ne_up_now = rdmft_weighted_sum_ispin(rdmft_n, wk, nbnd, nks, 1)
          ne_down_now = rdmft_weighted_sum_ispin(rdmft_n, wk, nbnd, nks, 2)
          mz_now = ne_up_now - ne_down_now
       ENDIF
       IF (rdmft_verbose >= 1 .AND. ionode) THEN
          IF (do_spin_report) THEN
             WRITE(stdout, '(5X,A,I4,A,F18.10,A,1PE10.2,0P,A,F12.6,A,F12.6,A,F12.6,A,F12.6,A,1PE10.2)') &
                  'RDMFT joint(SPG) outer ', outer, '  E = ', etot_now, &
                  '  |dE| = ', dE, '  Ne = ', ne_now, &
                  '  N_up = ', ne_up_now, '  N_down = ', ne_down_now, &
                  '  Mz = ', mz_now, '  alpha = ', alpha
          ELSE
             WRITE(stdout, '(5X,A,I4,A,F18.10,A,1PE10.2,0P,A,F12.6,A,1PE10.2,A,1PE10.2,A,1PE10.2)') &
                  'RDMFT joint(SPG) outer ', outer, '  E = ', etot_now, &
                  '  |dE| = ', dE, '  Ne = ', ne_now, &
                  '  ||g_1||=', g1_inf, '  ||G_C||=', gnorm_c, '  alpha = ', alpha
          ENDIF
       ENDIF
       IF (rdmft_verbose >= 1) &
            CALL rdmft_report_outer_wall_time(outer, t_outer0)
       !
       IF (dE == 0.0_DP) THEN
          IF (ionode) WRITE(stdout, '(5X,A)') &
               'RDMFT joint (SPG) optimisation converged (|dE| = 0).'
          loop_converged = .TRUE.
          EXIT
       ELSE IF (outer > 1 .AND. dE < rdmft_energy_tol) THEN
          IF (ionode) WRITE(stdout, '(5X,A)') &
               'RDMFT joint (SPG) optimisation converged.'
          loop_converged = .TRUE.
          EXIT
       ENDIF
       IF (rdmft_save_every > 0) THEN
          IF (MOD(outer, rdmft_save_every) == 0) &
               CALL rdmft_save_state(outer, etot_now)
       ENDIF
       etot_prev = etot_now
    ENDDO
    !
    IF (.NOT. loop_converged .AND. rdmft_outer_maxiter >= 1) THEN
       IF (ionode) THEN
          WRITE(stdout, '(/,5X,A)') &
               'WARNING: RDMFT joint (SPG) optimisation did NOT converge.'
          WRITE(stdout, '(5X,A,I0,A)') &
               '         The outer loop hit rdmft_outer_maxiter = ', &
               rdmft_outer_maxiter, ' without |dE| < rdmft_energy_tol.'
       ENDIF
    ENDIF
    CALL rdmft_total_energy(etot_now)
    !
    DEALLOCATE(grad_n, d_n, n_save, occ_flat, grad_flat)
    DEALLOCATE(hpsi, gradC, eta, vxpsi, vxpsi_total, vx_diag)
    DEALLOCATE(C_save, dir_C, G_C, prev_G_C, prev_dir_C)
    !
  END SUBROUTINE rdmft_run_joint_spg
  !
  !-----------------------------------------------------------------
  SUBROUTINE pack_joint(p, C, n_p, n_orb_per_k, nks, lda, nb, vec)
    !---------------------------------------------------------------
    !! Pack ``(p, C(:,:,1..nks))`` into a single real vector.
    INTEGER, INTENT(IN) :: n_p, n_orb_per_k, nks, lda, nb
    REAL(DP), INTENT(IN) :: p(n_p)
    COMPLEX(DP), INTENT(IN) :: C(lda, nb, nks)
    REAL(DP), INTENT(OUT) :: vec(n_p + n_orb_per_k * nks)
    INTEGER :: i, j, k, idx
    DO i = 1, n_p
       vec(i) = p(i)
    ENDDO
    idx = n_p
    DO k = 1, nks
       DO j = 1, nb
          DO i = 1, lda
             idx = idx + 1
             vec(idx) = REAL(C(i, j, k), KIND=DP)
             idx = idx + 1
             vec(idx) = AIMAG(C(i, j, k))
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE pack_joint
  !
  !-----------------------------------------------------------------
  SUBROUTINE unpack_joint(vec, n_p, n_orb_per_k, nks, lda, nb, p, C)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN) :: n_p, n_orb_per_k, nks, lda, nb
    REAL(DP), INTENT(IN) :: vec(n_p + n_orb_per_k * nks)
    REAL(DP), INTENT(OUT) :: p(n_p)
    COMPLEX(DP), INTENT(OUT) :: C(lda, nb, nks)
    INTEGER :: i, j, k, idx
    DO i = 1, n_p
       p(i) = vec(i)
    ENDDO
    idx = n_p
    DO k = 1, nks
       DO j = 1, nb
          DO i = 1, lda
             idx = idx + 1
             C(i, j, k) = CMPLX(vec(idx), 0.0_DP, KIND=DP)
             idx = idx + 1
             C(i, j, k) = C(i, j, k) + CMPLX(0.0_DP, vec(idx), KIND=DP)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE unpack_joint
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_orbital_step(etot_io, orb_converged)
    !---------------------------------------------------------------
    !! Dispatch to either the joint multi-k Stiefel descent (default)
    !! or the per-k block-coordinate variant based on the
    !! ``rdmft_orb_strategy`` keyword.
    USE rdmft_module, ONLY : rdmft_orb_strategy
    REAL(DP), INTENT(INOUT) :: etot_io
    LOGICAL, INTENT(OUT), OPTIONAL :: orb_converged
    LOGICAL :: conv
    conv = .FALSE.
    SELECT CASE (TRIM(rdmft_orb_strategy))
    CASE ('block_k', 'blockk', 'k')
       CALL rdmft_orbital_step_block_k(etot_io, conv)
    CASE DEFAULT
       CALL rdmft_orbital_step_joint(etot_io, conv)
    END SELECT
    IF (PRESENT(orb_converged)) orb_converged = conv
  END SUBROUTINE rdmft_orbital_step
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_orbital_step_joint(etot_io, orb_converged)
    !---------------------------------------------------------------
    !! Joint multi-k Stiefel descent on the natural orbital
    !! coefficients ``\{C^k\}_{k=1..nks}``.
    !!
    !! Per inner iteration, the routine:
    !!   1. computes the Riemannian gradient ``G_R^k`` at every
    !!      k-point (one h_psi + one ace channel per k);
    !!   2. builds the search direction ``d^k`` (SD / CG / L-BFGS) on
    !!      the **product** Stiefel manifold
    !!      \(\prod_k \mathrm{St}(\text{nbnd}, \text{npwx}; I)\);
    !!   3. runs **one Armijo line search per inner iteration** that
    !!      retracts ``C^k \to R_{C^k}(\alpha\, d^k)`` simultaneously
    !!      at every k-point and asks ``rdmft_total_energy`` for the
    !!      single global trial energy;
    !!   4. on acceptance, updates the per-k buffers; on failure,
    !!      restores all k-points and resets the optimiser state.
    !!
    !! This is the proper product-manifold formulation: the energy
    !! is a function of all ``C^k`` together (through the Hartree
    !! density and the cross-k Fock exchange) and only a global
    !! Armijo can give a monotone descent.  The previous per-k local
    !! line search did not see the cross-k coupling and could
    !! degrade the global energy even when each per-k step
    !! "succeeded".
    !
    USE io_global,         ONLY : stdout, ionode
    USE wvfct,             ONLY : nbnd, npwx, wg, current_k
    USE klist,             ONLY : nks, ngk, igk_k, xk, wk
    USE wavefunctions,     ONLY : evc
    USE noncollin_module,  ONLY : npol, noncolin
    USE control_flags,     ONLY : gamma_only
    USE gvect,             ONLY : gstart
    USE io_files,          ONLY : nwordwfc, iunwfc
    USE buffers,           ONLY : get_buffer, save_buffer
    USE lsda_mod,          ONLY : lsda, current_spin, isk
    USE uspp,              ONLY : nkb, vkb, okvan
    USE uspp_init,         ONLY : init_us_2
    USE rdmft_module
    USE rdmft_lbfgs
    USE rdmft_linesearch
    USE rdmft_orb_ls
    USE rdmft_stiefel
    USE rdmft_energy
    USE rdmft_xc,          ONLY : rdmft_xc_n_channels, rdmft_xc_channel
    USE rdmft_inner_log,   ONLY : rdmft_inner_log_record
    USE mp,                ONLY : mp_sum
    USE mp_pools,          ONLY : inter_pool_comm
    !
    REAL(DP), INTENT(INOUT) :: etot_io
    LOGICAL, INTENT(OUT), OPTIONAL :: orb_converged
    !
    ! Per-k arrays kept across the inner loop on the **product**
    ! Stiefel manifold.  Memory: 4 * npwx*npol*nbnd*nks complex
    ! (gradient + direction + save + try) + 1 channel buffer per k
    ! for vxpsi.  For Si 2x2x2 (nks=3, nbnd=8, npwx~350) this is
    ! ~270 KB per array, which fits comfortably in any modern
    ! machine.
    COMPLEX(DP), ALLOCATABLE :: G_R_all(:,:,:)     ! Riemannian gradient per k
    COMPLEX(DP), ALLOCATABLE :: G_R_pc_all(:,:,:)  ! preconditioned R-grad per k
    COMPLEX(DP), ALLOCATABLE :: prev_G_R_pc_all(:,:,:)! previous prec. Riem. grad
    COMPLEX(DP), ALLOCATABLE :: prev_G_R_all(:,:,:)! previous Riem. grad (raw)
    COMPLEX(DP), ALLOCATABLE :: dir_all(:,:,:)     ! search direction per k
    COMPLEX(DP), ALLOCATABLE :: prev_dir_all(:,:,:)
    COMPLEX(DP), ALLOCATABLE :: C_save_all(:,:,:)  ! evc snapshot per k
    COMPLEX(DP), ALLOCATABLE :: C_try(:,:)         ! retraction scratch
    COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), vxpsi(:,:,:), vxpsi_total(:,:,:), gradC(:,:), eta(:,:)
    COMPLEX(DP), ALLOCATABLE :: eta_pc(:,:)
    REAL(DP),    ALLOCATABLE :: vx_diag(:,:)
    REAL(DP),    ALLOCATABLE :: h_diag_all(:,:)
    REAL(DP),    ALLOCATABLE :: h_diag_k(:)
    !
    ! L-BFGS Re/Im packed view of all (npwx*npol, nbnd, nks) blocks.
    REAL(DP), ALLOCATABLE :: grad_flat(:), dir_flat(:), step_flat(:), grad_diff_flat(:)
    !
    REAL(DP) :: alpha, gnorm_total, slope, etot_trial, etot_start
    REAL(DP) :: gnorm_prev, beta, num, den, coef, w, dw, e_t
    REAL(DP) :: ref_energy
    INTEGER :: ik, npw, ib, inner, ls, it, nch, ndim_orb_total
    LOGICAL :: line_ok, use_cg_orb, use_lbfgs_orb, use_zh
    LOGICAL :: use_orb_wolfe, use_orb_armijo
    CHARACTER(LEN=24) :: orb_ls_eff
    TYPE(rdmft_ls_result) :: lsres
    !
    TYPE(rdmft_lbfgs_state)       :: orb_lbfgs
    TYPE(rdmft_zhang_hager_state) :: zh_orb
    !
    ! ---- BB / quadratic alpha0 state for the orbital line search ----
    ! Same machinery the SPG occupation block uses (see
    ! :var:`rdmft_orb_ls_init_step`).  The BB pair (s, y) lives in the
    ! Euclidean packed (Re/Im, per-k flattened) view of the orbital
    ! iterate -- standard practice for a Stiefel optimiser coupled with
    ! Stiefel retraction.  ``orb_iterate_flat`` holds the current
    ! packed evc; the BB recorder snapshots (orb_iterate_flat,
    ! grad_flat) BEFORE each accepted step so iter k+1 can compute
    !    s = orb_iterate_flat_{k+1} - orb_iterate_flat_k
    !    y = grad_flat_{k+1}      - grad_flat_k
    ! and propose alpha_{k+1} = (s^T s) / (s^T y).
    TYPE(rdmft_bb_state)          :: orb_bb
    TYPE(rdmft_quad_ls_history)   :: orb_qhist
    REAL(DP), ALLOCATABLE :: orb_iterate_flat(:)
    CALL rdmft_resolve_orb_ls_type(rdmft_orb_ls_type, orb_ls_eff, use_orb_wolfe, use_orb_armijo)
    use_cg_orb    = (TRIM(rdmft_orb_optimizer) == 'cg')
    use_lbfgs_orb = (TRIM(rdmft_orb_optimizer) == 'lbfgs') &
              .OR. (TRIM(rdmft_orb_optimizer) == 'bfgs')  &
              .OR. (TRIM(rdmft_orb_optimizer) == 'l-bfgs')
    !
    nch = rdmft_xc_n_channels()
    !
    ALLOCATE(G_R_all     (npwx*npol, nbnd, nks))
    ALLOCATE(G_R_pc_all  (npwx*npol, nbnd, nks))
    ALLOCATE(prev_G_R_all(npwx*npol, nbnd, nks))
    ALLOCATE(prev_G_R_pc_all(npwx*npol, nbnd, nks))
    ALLOCATE(dir_all     (npwx*npol, nbnd, nks))
    ALLOCATE(prev_dir_all(npwx*npol, nbnd, nks))
    ALLOCATE(C_save_all  (npwx*npol, nbnd, nks))
    ALLOCATE(C_try       (npwx*npol, nbnd))
    ALLOCATE(hpsi (npwx*npol, nbnd))
    ALLOCATE(gradC(npwx*npol, nbnd))
    ALLOCATE(eta  (npwx*npol, nbnd))
    ALLOCATE(eta_pc(npwx*npol, nbnd))
    ALLOCATE(vxpsi(npwx*npol, nbnd, nks))
    ! Channel-summed Vx|psi> over all k, recomputed once before the
    ! per-k Stiefel projection in the L-BFGS gradient build below.
    ALLOCATE(vxpsi_total(npwx*npol, nbnd, nks))
    ALLOCATE(vx_diag(nbnd, nks))
    ALLOCATE(h_diag_all(nbnd, nks), h_diag_k(nbnd))
    ndim_orb_total = 2 * npwx * npol * nbnd * nks
    ALLOCATE(grad_flat     (ndim_orb_total))
    ALLOCATE(dir_flat      (ndim_orb_total))
    ALLOCATE(step_flat     (ndim_orb_total))
    ALLOCATE(grad_diff_flat(ndim_orb_total))
    !
    use_zh = (TRIM(rdmft_line_search) == 'zhang_hager') &
        .OR. (TRIM(rdmft_line_search) == 'zh') &
        .OR. (TRIM(rdmft_line_search) == 'nonmonotone')
    !
    prev_G_R_all    = (0.0_DP, 0.0_DP)
    prev_G_R_pc_all = (0.0_DP, 0.0_DP)
    prev_dir_all    = (0.0_DP, 0.0_DP)
    gnorm_prev = 0.0_DP
    IF (use_lbfgs_orb) CALL rdmft_lbfgs_init(orb_lbfgs, ndim_orb_total, rdmft_lbfgs_memory)
    CALL rdmft_zhang_hager_init(zh_orb, etot_io, rdmft_zhang_hager_eta)
    !
    ! Initialise the orbital BB / quadratic-fit alpha0 state.  The
    ! packed iterate / gradient live in the same Re/Im flattening
    ! used by ``pack_psi_real`` for L-BFGS, so the existing
    ! ``ndim_orb_total`` is the right size.
    CALL rdmft_bb_init(orb_bb, ndim_orb_total, rdmft_bb_alpha_min, rdmft_bb_alpha_max)
    CALL rdmft_quad_ls_reset(orb_qhist)
    ALLOCATE(orb_iterate_flat(ndim_orb_total))
    !
    ! Pre-compute the per-band h_diag_k used by the orbital
    ! level-shift preconditioner.  Cheap (one h_psi sweep) compared
    ! to the inner-loop ACE rebuilds, and held constant during one
    ! orbital block (eps_band is updated at the start of the next
    ! outer cycle when h_diag changes).
    IF (rdmft_orb_precond) CALL rdmft_band_energies(h_diag_all)
    !
    IF (PRESENT(orb_converged)) orb_converged = .FALSE.
    !
    ! ---- Joint multi-k inner loop -----------------------------
    DO inner = 1, rdmft_orb_maxiter
       !
       ! 1) Riemannian gradient at every k; save evc for line search.
       DO ik = 1, nks
          current_k = ik
          IF (lsda) current_spin = isk(ik)
          IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
          C_save_all(:, :, ik) = evc(:, :)
       ENDDO
       CALL rdmft_compute_riemannian_gradient(G_R_all, gnorm_total)
       gnorm_total = SQRT(MAX(0.0_DP, gnorm_total))
       DO ik = 1, nks
          npw = ngk(ik)
          current_k = ik
          IF (lsda) current_spin = isk(ik)
          IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
          evc(:, :) = C_save_all(:, :, ik)
          eta(:, :) = G_R_all(:, :, ik)
          IF (rdmft_orb_precond) THEN
             h_diag_k(:) = h_diag_all(:, ik)
             CALL rdmft_apply_orb_precond(eta, eta_pc, npw, nbnd, &
                                          npwx*npol, h_diag_k, &
                                          rdmft_orb_precond_shift)
             IF (gamma_only) THEN
                CALL stiefel_project_tangent_gamma(evc, eta_pc, npw, nbnd, &
                     npwx*npol, gstart, eta)
             ELSE
                CALL stiefel_project_tangent_k(evc, eta_pc, npw, nbnd, &
                     npwx*npol, eta)
             ENDIF
             G_R_pc_all(:, :, ik) = eta(:, :)
          ELSE
             G_R_pc_all(:, :, ik) = G_R_all(:, :, ik)
          ENDIF
       ENDDO
       IF (gnorm_total < rdmft_orb_grad_tol) THEN
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,1PE10.2,A,1PE10.2,A)') &
               'ORB block converged: ||G_R|| = ', gnorm_total, &
               ' < ', rdmft_orb_grad_tol, ' (rdmft_orb_grad_tol)'
          IF (PRESENT(orb_converged)) orb_converged = .TRUE.
          EXIT
       ENDIF
       !
       ! 2) Build the search direction on the product manifold.
       !    SD: dir = -P*G_R.
       !    CG: Polak-Ribière on the packed (per-k) **preconditioned**
       !        gradient with transport-by-reprojection.
       !    L-BFGS: Euclidean two-loop on the Re/Im flattened
       !        preconditioned product vector, re-projected per-k onto
       !        the tangent.
       IF (use_lbfgs_orb) THEN
          DO ik = 1, nks
             CALL pack_psi_real(G_R_pc_all(:, :, ik), npwx*npol, nbnd, &
                  grad_flat( (ik-1)*2*npwx*npol*nbnd + 1 :  ik*2*npwx*npol*nbnd ))
          ENDDO
          CALL rdmft_lbfgs_direction(orb_lbfgs, grad_flat, dir_flat)
          DO ik = 1, nks
             CALL unpack_psi_real( &
                  dir_flat( (ik-1)*2*npwx*npol*nbnd + 1 :  ik*2*npwx*npol*nbnd ), &
                  npwx*npol, nbnd, gradC)
             ! Re-project onto the tangent at the current C^k.
             current_k = ik
             IF (lsda) current_spin = isk(ik)
             IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
             npw = ngk(ik)
             IF (gamma_only) THEN
                CALL stiefel_project_tangent_gamma(evc, gradC, npw, nbnd, &
                     npwx*npol, gstart, dir_all(:, :, ik))
             ELSE
                CALL stiefel_project_tangent_k(evc, gradC, npw, nbnd, &
                     npwx*npol, dir_all(:, :, ik))
             ENDIF
          ENDDO
       ELSE IF (use_cg_orb .AND. inner > 1 .AND. gnorm_prev > 0.0_DP) THEN
          ! Preconditioned Polak-Ribière on the product manifold:
          !   beta = max(0,
          !       sum_k <G_pc^k - G_pc^k_prev, G_pc^k>_S
          !       / sum_k <G^k_prev, G_pc^k_prev>_S)
          ! Numerator and denominator both live in the same
          ! (preconditioned) metric so beta has the right scale.
          num = 0.0_DP
          DO ik = 1, nks
             npw = ngk(ik)
             ! Use the workspace `eta` (already allocated) to hold the
             ! preconditioned gradient difference for this k-point.
             eta(:, :) = G_R_pc_all(:, :, ik) - prev_G_R_pc_all(:, :, ik)
             IF (gamma_only) THEN
                num = num + stiefel_inner_product_gamma( &
                     eta, G_R_pc_all(:, :, ik), npw, nbnd, npwx*npol, gstart)
             ELSE
                num = num + stiefel_inner_product_k( &
                     eta, G_R_pc_all(:, :, ik), npw, nbnd, npwx*npol)
             ENDIF
          ENDDO
          ! Preconditioned CG denominator: <G_prev, G_pc_prev>_S.
          den = 0.0_DP
          DO ik = 1, nks
             npw = ngk(ik)
             IF (gamma_only) THEN
                den = den + stiefel_inner_product_gamma( &
                     prev_G_R_all(:, :, ik), prev_G_R_pc_all(:, :, ik), &
                     npw, nbnd, npwx*npol, gstart)
             ELSE
                den = den + stiefel_inner_product_k( &
                     prev_G_R_all(:, :, ik), prev_G_R_pc_all(:, :, ik), &
                     npw, nbnd, npwx*npol)
             ENDIF
          ENDDO
          ! num / den sum the per-k Stiefel inner products over this
          ! pool's k-points only; reduce over the k-point pools so beta
          ! is identical on every pool (a per-pool beta would build a
          ! different search direction per pool and desynchronise the
          ! collective line search).
          CALL mp_sum(num, inter_pool_comm)
          CALL mp_sum(den, inter_pool_comm)
          IF (den > 1.0e-30_DP) THEN
             beta = MAX(0.0_DP, num / den)
          ELSE
             beta = 0.0_DP
          ENDIF
          DO ik = 1, nks
             dir_all(:, :, ik) = -G_R_pc_all(:, :, ik) + beta * prev_dir_all(:, :, ik)
          ENDDO
       ELSE
          dir_all = -G_R_pc_all
       ENDIF
       !
       ! Riemannian directional derivative on the product manifold:
       !   slope = sum_k <G_R^k, dir^k>_S
       ! (note: the slope uses the **raw** Riemannian gradient
       ! G_R_all, not the preconditioned one; the Armijo test below
       ! is a condition on the directional derivative of E, which
       ! is independent of any preconditioner.)
       slope = 0.0_DP
       DO ik = 1, nks
          npw = ngk(ik)
          IF (gamma_only) THEN
             slope = slope + stiefel_inner_product_gamma( &
                  G_R_all(:, :, ik), dir_all(:, :, ik), npw, nbnd, npwx*npol, gstart)
          ELSE
             slope = slope + stiefel_inner_product_k( &
                  G_R_all(:, :, ik), dir_all(:, :, ik), npw, nbnd, npwx*npol)
          ENDIF
       ENDDO
       ! Sum the per-k product-manifold slope over the k-point pools.
       CALL mp_sum(slope, inter_pool_comm)
       IF (slope >= 0.0_DP) THEN
          ! Direction is not a descent on the product manifold;
          ! fall back to packed (preconditioned) steepest descent.
          dir_all = -G_R_pc_all
          slope = 0.0_DP
          DO ik = 1, nks
             npw = ngk(ik)
             IF (gamma_only) THEN
                slope = slope - stiefel_inner_product_gamma( &
                     G_R_all(:, :, ik), G_R_pc_all(:, :, ik), &
                     npw, nbnd, npwx*npol, gstart)
             ELSE
                slope = slope - stiefel_inner_product_k( &
                     G_R_all(:, :, ik), G_R_pc_all(:, :, ik), &
                     npw, nbnd, npwx*npol)
             ENDIF
          ENDDO
          CALL mp_sum(slope, inter_pool_comm)
          IF (use_lbfgs_orb) CALL rdmft_lbfgs_reset(orb_lbfgs)
       ENDIF
       prev_G_R_all    = G_R_all
       prev_G_R_pc_all = G_R_pc_all
       prev_dir_all    = dir_all
       gnorm_prev      = gnorm_total
       !
       ! 3) Line search along the product direction (strong Wolfe by
       !    default, or Armijo / Zhang-Hager when rdmft_orb_ls_type
       !    requests Armijo).
       etot_start = etot_io
       !
       ! Pack the current orbital iterate (= C_save_all, the pre-step
       ! evc snapshot) and the raw Riemannian gradient (G_R_all) into
       ! their Re/Im, per-k flattened views.  ``rdmft_occ_ls_alpha0``
       ! is shared with the SPG occupation block: with
       ! ``rdmft_orb_ls_init_step = 'barzilai_borwein'`` (the default)
       ! it consumes the previous (x_{k-1}, g_{k-1}) snapshot in
       ! ``orb_bb`` and proposes the BB step
       !
       !     alpha_0 = (s^T s) / (s^T y),
       !     s = x_k - x_{k-1},   y = g_k - g_{k-1}
       !
       ! clamped to [rdmft_bb_alpha_min, rdmft_bb_alpha_max].  On the
       ! very first inner step ``orb_bb`` has no history, so the
       ! routine falls back to ``rdmft_orb_ls_stepsize``; the same
       ! holds with ``rdmft_orb_ls_init_step = 'fixed'``, which
       ! reproduces the legacy behaviour of always starting at
       ! ``rdmft_orb_ls_stepsize``.
       !
       ! BB on a Stiefel manifold is computed in the **Euclidean**
       ! ambient packed view, exactly the same way ABACUS's
       ! ``joint_bb`` (see ``abacus_source/source_lcao/module_rdmft/
       ! rdmft_solver.cpp``) handles the joint product-manifold
       ! step, and the same way the L-BFGS update in this routine
       ! already operates on (s, y) flattened across all k-points.
       DO ik = 1, nks
          CALL pack_psi_real(C_save_all(:, :, ik), npwx*npol, nbnd, &
               orb_iterate_flat( (ik-1)*2*npwx*npol*nbnd + 1 :  ik*2*npwx*npol*nbnd ))
          CALL pack_psi_real(G_R_all(:, :, ik), npwx*npol, nbnd, &
               grad_flat( (ik-1)*2*npwx*npol*nbnd + 1 :  ik*2*npwx*npol*nbnd ))
       ENDDO
       alpha = rdmft_occ_ls_alpha0(rdmft_orb_ls_init_step, orb_iterate_flat, &
            grad_flat, rdmft_orb_ls_stepsize, orb_bb, orb_qhist, &
            rdmft_bb_alpha_min, rdmft_bb_alpha_max)
       IF (alpha <= 0.0_DP) alpha = rdmft_orb_ls_stepsize
       !
       ! Snapshot (x_k, g_k) into ``orb_bb`` BEFORE the step so that
       ! iter k+1 can form s = x_{k+1} - x_k, y = g_{k+1} - g_k (same
       ! BB-bookkeeping order as in the SPG occupation block; see the
       ! comment there for why the post-step record was a bug).
       CALL rdmft_bb_record(orb_bb, orb_iterate_flat, grad_flat)
       !
       line_ok = .FALSE.
       IF (.NOT. rdmft_orb_no_ls .AND. use_orb_wolfe) THEN
          CALL rdmft_orb_ls_setup_joint(C_save_all, dir_all, nks, npwx, nbnd, npol, &
               gamma_only, gstart, nch)
          CALL rdmft_strong_wolfe_ls(rdmft_orb_ls_eval, etot_start, slope, alpha, &
               rdmft_line_search_c1, rdmft_line_search_c2, &
               rdmft_line_search_max_iter, rdmft_line_search_max_zoom, lsres)
          IF (lsres%success .AND. lsres%step > 0.0_DP) THEN
             alpha = lsres%step
             etot_io = lsres%f_new
             line_ok = .TRUE.
          ELSE IF (rdmft_verbose >= 1 .AND. ionode) THEN
             WRITE(stdout, '(7X,A)') &
                  'ORB block: strong Wolfe failed; Armijo fallback.'
          ENDIF
       ENDIF
       IF (.NOT. rdmft_orb_no_ls .AND. .NOT. line_ok) THEN
          IF (use_zh) THEN
             ref_energy = MAX(rdmft_zhang_hager_ref(zh_orb), etot_io)
          ELSE
             ref_energy = etot_start
          ENDIF
          DO ls = 1, rdmft_line_search_max_iter
             DO ik = 1, nks
                npw = ngk(ik)
                current_k = ik
                IF (lsda) current_spin = isk(ik)
                IF (gamma_only) THEN
                   CALL stiefel_retract_gamma(C_save_all(:, :, ik), dir_all(:, :, ik), &
                        alpha, npw, nbnd, npwx*npol, gstart, C_try)
                ELSE
                   CALL stiefel_retract_k(C_save_all(:, :, ik), dir_all(:, :, ik), &
                        alpha, npw, nbnd, npwx*npol, C_try)
                ENDIF
                evc(:, :) = C_try(:, :)
                IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
             ENDDO
             CALL rdmft_mark_orbitals_changed()
             CALL rdmft_total_energy(etot_trial)
             IF (etot_trial <= ref_energy + rdmft_line_search_c1 * alpha * slope) THEN
                line_ok = .TRUE.
                etot_io = etot_trial
                EXIT
             ENDIF
             IF (rdmft_line_search_polynomial) THEN
                CALL rdmft_armijo_next_alpha(alpha, ref_energy, slope, etot_trial, &
                     rdmft_line_search_rho, alpha)
             ELSE
                alpha = alpha * rdmft_line_search_rho
             ENDIF
          ENDDO
          IF (use_zh .AND. line_ok) CALL rdmft_zhang_hager_update(zh_orb, etot_io)
       ENDIF
       !
       IF (.NOT. line_ok) THEN
          ! Restore all k-points to the start of this inner iteration.
          DO ik = 1, nks
             current_k = ik
             IF (lsda) current_spin = isk(ik)
             evc(:, :) = C_save_all(:, :, ik)
             IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
          ENDDO
          CALL rdmft_mark_orbitals_changed()
          ! Refresh etot_io at the restored point.
          CALL rdmft_total_energy(etot_io)
          IF (use_lbfgs_orb) CALL rdmft_lbfgs_reset(orb_lbfgs)
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,I0,A,1PE10.2)') &
               'ORB block stopped: line search failed at inner ', inner, &
               '; ||G_R|| = ', gnorm_total
          EXIT
       ENDIF
       !
       ! 4) Update L-BFGS history with the accepted product step.
       IF (use_lbfgs_orb) THEN
          ! Build (s, y) on the packed product vector.  s = alpha *
          ! dir; y = G_R(new) - G_R(old).  We need the gradient at
          ! the new point.  Recompute it (this is the same work as
          ! step 1 of the next inner iteration; for inner < maxiter
          ! it's cheap to do here so the curvature pair is current).
          DO ik = 1, nks
             CALL pack_psi_real(G_R_all(:, :, ik), npwx*npol, nbnd, &
                  grad_diff_flat( (ik-1)*2*npwx*npol*nbnd + 1 :  ik*2*npwx*npol*nbnd ))
          ENDDO
          ! The new gradient will be computed at the start of the
          ! next iteration (step 1).  We store the current gradient
          ! here and use it as "old" next time round; the solver then
          ! computes (G_new - G_old) and calls rdmft_lbfgs_update.
          ! That bookkeeping is handled by `prev_G_R_all` already, so
          ! all we need to remember here is the actual step.
          DO ik = 1, nks
             dir_all(:, :, ik) = alpha * dir_all(:, :, ik)
             CALL pack_psi_real(dir_all(:, :, ik), npwx*npol, nbnd, &
                  step_flat( (ik-1)*2*npwx*npol*nbnd + 1 :  ik*2*npwx*npol*nbnd ))
             ! Restore dir_all so prev_dir_all (saved earlier) is
             ! consistent with subsequent CG iterations.
             dir_all(:, :, ik) = dir_all(:, :, ik) / MAX(alpha, 1.0e-30_DP)
          ENDDO
          ! Compute the new gradient exactly so we can build y = G_new - G_old.
          !
          ! Cost note: ``rdmft_compute_xc_channel`` does an all-k aceinit
          ! (``O(N_k^2)``); calling it inside the per-k loop below would
          ! make the gradient build cost ``O(N_k^3)``.  Pre-compute the
          ! channel-summed ``Vx|psi>`` for every k once, then reuse it.
          vxpsi_total = (0.0_DP, 0.0_DP)
          DO it = 1, nch
             CALL rdmft_compute_xc_channel(it, e_t, vx_diag, vxpsi)
             ! Distinct (ib, ik) slices: thread-safe channel accumulation.
!$omp parallel do collapse(2) default(shared) private(ik, ib, coef, w, dw)
             DO ik = 1, nks
                DO ib = 1, nbnd
                   CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
                   vxpsi_total(:, ib, ik) = vxpsi_total(:, ib, ik) &
                                          + coef * w * vxpsi(:, ib, ik)
                ENDDO
             ENDDO
!$omp end parallel do
          ENDDO
          DO ik = 1, nks
             npw = ngk(ik)
             current_k = ik
             IF (lsda) current_spin = isk(ik)
             IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
             CALL g2_kin(ik)
             IF (nkb > 0) CALL init_us_2(npw, igk_k(1, ik), xk(:, ik), vkb)
             CALL rdmft_apply_h_one_psi(npw, nbnd, evc, hpsi)
             DO ib = 1, nbnd
                gradC(:, ib) = wk(ik) * (rdmft_n(ib, ik) * hpsi(:, ib) &
                                         + vxpsi_total(:, ib, ik))
             ENDDO
             IF (gamma_only) THEN
                CALL stiefel_project_tangent_gamma(evc, gradC, npw, nbnd, &
                     npwx*npol, gstart, eta)
             ELSE
                CALL stiefel_project_tangent_k(evc, gradC, npw, nbnd, &
                     npwx*npol, eta)
             ENDIF
             CALL pack_psi_real(eta, npwx*npol, nbnd, &
                  grad_flat( (ik-1)*2*npwx*npol*nbnd + 1 :  ik*2*npwx*npol*nbnd ))
          ENDDO
          ! y = grad_new - grad_old; we have grad_old in grad_diff_flat
          ! and grad_new in grad_flat now.
          grad_diff_flat = grad_flat - grad_diff_flat
          CALL rdmft_lbfgs_update(orb_lbfgs, step_flat, grad_diff_flat)
       ENDIF
       !
       IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,I4,A,F18.10,A,1PE10.2,0P,A,1PE10.2)') &
            'RDMFT orb inner ', inner, ' end: global E = ', etot_io, &
            '  ||G_R|| = ', gnorm_total, '  alpha = ', alpha
       CALL rdmft_inner_log_record('orb', inner, 0, etot_io, gnorm_total, &
            alpha, 0.0_DP, 0.0_DP)
    ENDDO
    !
    IF (use_lbfgs_orb) CALL rdmft_lbfgs_finalize(orb_lbfgs)
    !
    DEALLOCATE(G_R_all, G_R_pc_all, prev_G_R_all, prev_G_R_pc_all)
    DEALLOCATE(dir_all, prev_dir_all, C_save_all, C_try)
    DEALLOCATE(hpsi, gradC, eta, eta_pc, vxpsi, vxpsi_total, vx_diag)
    DEALLOCATE(h_diag_all, h_diag_k)
    DEALLOCATE(grad_flat, dir_flat, step_flat, grad_diff_flat)
    IF (ALLOCATED(orb_iterate_flat)) DEALLOCATE(orb_iterate_flat)
    !
  END SUBROUTINE rdmft_orbital_step_joint
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_orbital_step_block_k(etot_io, orb_converged)
    !---------------------------------------------------------------
    !! Per-k block-coordinate Stiefel descent.
    !!
    !! Strategy (Gauss-Seidel sweep over k-points):
    !!
    !!   DO inner = 1, rdmft_orb_maxiter
    !!     DO ik = 1, nks
    !!       1. load evc[ik] from the wfc buffer; snapshot C_save = evc.
    !!       2. compute the Riemannian gradient G_R^k = proj(grad_C^k)
    !!          AT THIS k ONLY, holding all other k's fixed (their
    !!          evc remains in the buffer and contributes to rho /
    !!          gamma_xc unchanged).
    !!       3. build the search direction with this k's own
    !!          optimiser state (per-k SD / CG / L-BFGS history,
    !!          persisted across outer iterations).
    !!       4. line search on the GLOBAL energy
    !!          (rdmft_total_energy summed over all k); only
    !!          evc[ik] changes between trials -- the other k-points
    !!          stay fixed in their buffers.
    !!       5. on acceptance, save evc[ik] to its buffer and update
    !!          the per-k L-BFGS curvature pair.  On rejection,
    !!          restore evc[ik] to C_save and reset this k's
    !!          optimiser state.
    !!     ENDDO
    !!   ENDDO
    !!
    !! Why this matters for parallelism: in QE's k-point pool
    !! parallelisation each MPI pool owns a subset of the k-points.
    !! Block_k naturally maps to this layout: every pool optimises
    !! its own k-points independently, and only the global energy
    !! evaluation at each line-search trial requires a cross-pool
    !! ``mp_sum`` (already performed inside ``rdmft_total_energy``
    !! through ``sum_band`` / ``v_of_rho`` / the energy reductions).
    !!
    !! For nks = 1 (gamma-only) block_k reduces to single-k Stiefel
    !! descent and is essentially the same as ``joint`` (the
    !! product manifold collapses to a single Stiefel factor).
    !
    USE io_global,         ONLY : stdout, ionode
    USE wvfct,             ONLY : nbnd, npwx, current_k
    USE klist,             ONLY : nks, ngk, igk_k, xk, wk
    USE wavefunctions,     ONLY : evc
    USE noncollin_module,  ONLY : npol
    USE control_flags,     ONLY : gamma_only
    USE gvect,             ONLY : gstart
    USE io_files,          ONLY : nwordwfc, iunwfc
    USE buffers,           ONLY : get_buffer, save_buffer
    USE lsda_mod,          ONLY : lsda, current_spin, isk
    USE uspp,              ONLY : nkb, vkb
    USE uspp_init,         ONLY : init_us_2
    USE rdmft_module
    USE rdmft_lbfgs
    USE rdmft_linesearch
    USE rdmft_orb_ls
    USE rdmft_stiefel
    USE rdmft_energy
    USE rdmft_xc,          ONLY : rdmft_xc_n_channels, rdmft_xc_channel
    USE rdmft_inner_log,   ONLY : rdmft_inner_log_record
    !
    REAL(DP), INTENT(INOUT) :: etot_io
    LOGICAL, INTENT(OUT), OPTIONAL :: orb_converged
    !
    ! Per-k state (persisted across outer iterations).
    COMPLEX(DP), ALLOCATABLE, SAVE :: prev_G_R_pk(:,:,:)
    COMPLEX(DP), ALLOCATABLE, SAVE :: prev_dir_pk(:,:,:)
    REAL(DP),    ALLOCATABLE, SAVE :: gnorm_prev_pk(:)
    LOGICAL,     ALLOCATABLE, SAVE :: have_prev_pk(:)
    !
    ! Per-k BB state for the orbital line-search alpha0 (one slot per
    ! k-point because the BLOCK_K sweep updates one k at a time, so
    ! the BB increment s = x_{k+1} - x_k naturally lives on a single k).
    TYPE(rdmft_bb_state), ALLOCATABLE, SAVE :: orb_bb_pk(:)
    LOGICAL,                              SAVE :: bb_pk_inited = .FALSE.
    !
    COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), gradC(:,:), eta(:,:)
    COMPLEX(DP), ALLOCATABLE :: C_save(:,:), C_try(:,:), dir(:,:)
    COMPLEX(DP), ALLOCATABLE :: vxpsi(:,:,:)
    REAL(DP),    ALLOCATABLE :: vx_diag(:,:)
    REAL(DP),    ALLOCATABLE :: grad_flat(:), dir_flat(:), step_flat(:), grad_diff_flat(:)
    REAL(DP),    ALLOCATABLE :: orb_iterate_flat(:)
    TYPE(rdmft_quad_ls_history) :: orb_qhist  ! unused for now (only BB / fixed)
    !
    ! Per-k L-BFGS state arrays (one bag per k; only when
    ! rdmft_orb_optimizer = lbfgs).
    TYPE(rdmft_lbfgs_state), ALLOCATABLE, SAVE :: orb_lbfgs_pk(:)
    LOGICAL, SAVE :: lbfgs_pk_inited = .FALSE.
    !
    REAL(DP) :: alpha, gnorm, slope, etot_trial, etot_start
    REAL(DP) :: coef, w, dw, e_t, beta, num, den
    INTEGER :: ik, npw, ib, inner, ls, it, nch, ndim_orb_per_k
    LOGICAL :: line_ok, use_cg_orb, use_lbfgs_orb, any_active_k
    LOGICAL :: use_orb_wolfe, use_orb_armijo
    CHARACTER(LEN=24) :: orb_ls_eff
    TYPE(rdmft_ls_result) :: lsres
    !
    CALL rdmft_resolve_orb_ls_type(rdmft_orb_ls_type, orb_ls_eff, use_orb_wolfe, use_orb_armijo)
    use_cg_orb    = (TRIM(rdmft_orb_optimizer) == 'cg')
    use_lbfgs_orb = (TRIM(rdmft_orb_optimizer) == 'lbfgs') &
              .OR. (TRIM(rdmft_orb_optimizer) == 'bfgs')  &
              .OR. (TRIM(rdmft_orb_optimizer) == 'l-bfgs')
    !
    nch = rdmft_xc_n_channels()
    ndim_orb_per_k = 2 * npwx * npol * nbnd
    !
    ! ---- Allocate / re-init the per-k state if needed ----------
    IF (.NOT. ALLOCATED(prev_G_R_pk) .OR. SIZE(prev_G_R_pk, 3) /= nks) THEN
       IF (ALLOCATED(prev_G_R_pk)) DEALLOCATE(prev_G_R_pk)
       IF (ALLOCATED(prev_dir_pk)) DEALLOCATE(prev_dir_pk)
       IF (ALLOCATED(gnorm_prev_pk)) DEALLOCATE(gnorm_prev_pk)
       IF (ALLOCATED(have_prev_pk)) DEALLOCATE(have_prev_pk)
       ALLOCATE(prev_G_R_pk(npwx*npol, nbnd, nks))
       ALLOCATE(prev_dir_pk(npwx*npol, nbnd, nks))
       ALLOCATE(gnorm_prev_pk(nks))
       ALLOCATE(have_prev_pk(nks))
       prev_G_R_pk   = (0.0_DP, 0.0_DP)
       prev_dir_pk   = (0.0_DP, 0.0_DP)
       gnorm_prev_pk = 0.0_DP
       have_prev_pk  = .FALSE.
    ENDIF
    IF (use_lbfgs_orb) THEN
       IF (.NOT. lbfgs_pk_inited .OR. .NOT. ALLOCATED(orb_lbfgs_pk) &
            .OR. SIZE(orb_lbfgs_pk) /= nks) THEN
          IF (ALLOCATED(orb_lbfgs_pk)) THEN
             DO ik = 1, SIZE(orb_lbfgs_pk)
                CALL rdmft_lbfgs_finalize(orb_lbfgs_pk(ik))
             ENDDO
             DEALLOCATE(orb_lbfgs_pk)
          ENDIF
          ALLOCATE(orb_lbfgs_pk(nks))
          DO ik = 1, nks
             CALL rdmft_lbfgs_init(orb_lbfgs_pk(ik), ndim_orb_per_k, rdmft_lbfgs_memory)
          ENDDO
          lbfgs_pk_inited = .TRUE.
       ENDIF
    ENDIF
    !
    ! Per-k BB state for ``rdmft_orb_ls_init_step = 'barzilai_borwein'``.
    ! Always allocated (cheap), independent of the orb optimiser, so
    ! BB also seeds alpha0 for SD / CG / fixed paths.
    IF (.NOT. bb_pk_inited .OR. .NOT. ALLOCATED(orb_bb_pk) &
         .OR. SIZE(orb_bb_pk) /= nks) THEN
       IF (ALLOCATED(orb_bb_pk)) DEALLOCATE(orb_bb_pk)
       ALLOCATE(orb_bb_pk(nks))
       DO ik = 1, nks
          CALL rdmft_bb_init(orb_bb_pk(ik), ndim_orb_per_k, &
               rdmft_bb_alpha_min, rdmft_bb_alpha_max)
       ENDDO
       bb_pk_inited = .TRUE.
    ENDIF
    CALL rdmft_quad_ls_reset(orb_qhist)
    ALLOCATE(orb_iterate_flat(ndim_orb_per_k))
    !
    ALLOCATE(hpsi(npwx*npol, nbnd))
    ALLOCATE(gradC(npwx*npol, nbnd))
    ALLOCATE(eta  (npwx*npol, nbnd))
    ALLOCATE(C_save(npwx*npol, nbnd))
    ALLOCATE(C_try (npwx*npol, nbnd))
    ALLOCATE(dir   (npwx*npol, nbnd))
    ALLOCATE(vxpsi (npwx*npol, nbnd, nks))
    ALLOCATE(vx_diag(nbnd, nks))
    ALLOCATE(grad_flat     (ndim_orb_per_k))
    ALLOCATE(dir_flat      (ndim_orb_per_k))
    ALLOCATE(step_flat     (ndim_orb_per_k))
    ALLOCATE(grad_diff_flat(ndim_orb_per_k))
    !
    IF (PRESENT(orb_converged)) orb_converged = .FALSE.
    !
    ! ---- Block-coordinate sweep -------------------------------
    DO inner = 1, rdmft_orb_maxiter
       any_active_k = .FALSE.
       DO ik = 1, nks
          npw = ngk(ik)
          current_k = ik
          IF (lsda) current_spin = isk(ik)
          IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
          CALL g2_kin(ik)
          IF (nkb > 0) CALL init_us_2(npw, igk_k(1, ik), xk(:, ik), vkb)
          C_save(:, :) = evc(:, :)
          !
          ! 1) Riemannian gradient at this k (one h_psi + one ace
          !    channel per k).
          CALL rdmft_apply_h_one_psi(npw, nbnd, evc, hpsi)
          gradC = (0.0_DP, 0.0_DP)
          DO ib = 1, nbnd
             gradC(:, ib) = wk(ik) * rdmft_n(ib, ik) * hpsi(:, ib)
          ENDDO
          ! Single-k aceinit_k build (only_k=ik): the BLOCK_K sweep
          ! changes only evc(:,:,ik) per inner step, so the all-k
          ! ACE projector ``xi`` is otherwise reusable.  The
          ! ``only_k`` fast path costs ``O(N_k)`` instead of the
          ! ``O(N_k^2)`` full-aceinit rebuild that would otherwise
          ! happen here.
          DO it = 1, nch
             CALL rdmft_compute_xc_channel(it, e_t, vx_diag, vxpsi, only_k=ik)
             DO ib = 1, nbnd
                CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
                gradC(:, ib) = gradC(:, ib) + coef * wk(ik) * w * vxpsi(:, ib, ik)
             ENDDO
          ENDDO
          ! Restore evc / current_k after the per-channel single-k build.
          current_k = ik
          IF (lsda) current_spin = isk(ik)
          evc(:, :) = C_save(:, :)
          IF (gamma_only) THEN
             CALL stiefel_project_tangent_gamma(evc, gradC, npw, nbnd, &
                  npwx*npol, gstart, eta)
             gnorm = SQRT(MAX(0.0_DP, &
                  stiefel_inner_product_gamma(eta, eta, npw, nbnd, npwx*npol, gstart)))
          ELSE
             CALL stiefel_project_tangent_k(evc, gradC, npw, nbnd, &
                  npwx*npol, eta)
             gnorm = SQRT(MAX(0.0_DP, &
                  stiefel_inner_product_k(eta, eta, npw, nbnd, npwx*npol)))
          ENDIF
          IF (gnorm < rdmft_orb_grad_tol) CYCLE  ! this k already at stationary
          any_active_k = .TRUE.
          !
          ! 2) Search direction with per-k optimiser state.
          IF (use_lbfgs_orb) THEN
             CALL pack_psi_real(eta, npwx*npol, nbnd, grad_flat)
             CALL rdmft_lbfgs_direction(orb_lbfgs_pk(ik), grad_flat, dir_flat)
             CALL unpack_psi_real(dir_flat, npwx*npol, nbnd, gradC)
             IF (gamma_only) THEN
                CALL stiefel_project_tangent_gamma(evc, gradC, npw, nbnd, &
                     npwx*npol, gstart, dir)
             ELSE
                CALL stiefel_project_tangent_k(evc, gradC, npw, nbnd, &
                     npwx*npol, dir)
             ENDIF
             ! Slope on the tangent at this k.
             IF (gamma_only) THEN
                slope = stiefel_inner_product_gamma(eta, dir, npw, nbnd, npwx*npol, gstart)
             ELSE
                slope = stiefel_inner_product_k(eta, dir, npw, nbnd, npwx*npol)
             ENDIF
             IF (slope >= 0.0_DP) THEN
                dir = -eta
                slope = -gnorm * gnorm
                CALL rdmft_lbfgs_reset(orb_lbfgs_pk(ik))
             ENDIF
          ELSE IF (use_cg_orb .AND. have_prev_pk(ik) .AND. gnorm_prev_pk(ik) > 0.0_DP) THEN
             ! Use the Stiefel inner product so the gamma-trick factors
             ! and the plane-wave parallel mp_sum match the ones in
             ! gnorm_prev_pk; see the comment in
             ! rdmft_orbital_step_joint for the failure mode of the
             ! naive `SUM(Re*Re + Im*Im)` Euclidean version.
             dir(:, :) = eta(:, :) - prev_G_R_pk(:, :, ik)
             IF (gamma_only) THEN
                num = stiefel_inner_product_gamma( &
                     dir, eta, npw, nbnd, npwx*npol, gstart)
             ELSE
                num = stiefel_inner_product_k( &
                     dir, eta, npw, nbnd, npwx*npol)
             ENDIF
             den = gnorm_prev_pk(ik) * gnorm_prev_pk(ik)
             IF (den > 1.0e-30_DP) THEN
                beta = MAX(0.0_DP, num / den)
             ELSE
                beta = 0.0_DP
             ENDIF
             dir = -eta + beta * prev_dir_pk(:, :, ik)
             IF (gamma_only) THEN
                slope = stiefel_inner_product_gamma(eta, dir, npw, nbnd, npwx*npol, gstart)
             ELSE
                slope = stiefel_inner_product_k(eta, dir, npw, nbnd, npwx*npol)
             ENDIF
             IF (slope >= 0.0_DP) THEN
                dir = -eta
                slope = -gnorm * gnorm
             ENDIF
          ELSE
             dir = -eta
             slope = -gnorm * gnorm
          ENDIF
          !
          ! Save current per-k Riemannian gradient for the next CG /
          ! L-BFGS curvature pair.
          prev_G_R_pk (:, :, ik) = eta
          prev_dir_pk(:, :, ik) = dir
          gnorm_prev_pk(ik)     = gnorm
          have_prev_pk (ik)     = .TRUE.
          !
          ! 3) GLOBAL line search: only evc[ik] changes per trial, but
          !    the energy is the global one.
          etot_start = etot_io
          !
          ! Pack this k's iterate (C_save) and Riemannian gradient
          ! (eta) for BB / quadratic-fit alpha0 estimation.  Each k
          ! has its own (orb_bb_pk(ik)) state because the BLOCK_K
          ! sweep updates one k at a time -- the BB increment
          ! s = x_{k+1}^{(ik)} - x_k^{(ik)} naturally lives on this
          ! single k.
          CALL pack_psi_real(C_save, npwx*npol, nbnd, orb_iterate_flat)
          CALL pack_psi_real(eta,    npwx*npol, nbnd, grad_flat)
          alpha = rdmft_occ_ls_alpha0(rdmft_orb_ls_init_step, &
               orb_iterate_flat, grad_flat, rdmft_orb_ls_stepsize, &
               orb_bb_pk(ik), orb_qhist, &
               rdmft_bb_alpha_min, rdmft_bb_alpha_max)
          IF (alpha <= 0.0_DP) alpha = rdmft_orb_ls_stepsize
          ! Snapshot (x_k, g_k) into the per-k BB state BEFORE the
          ! step so iter k+1 (next BLOCK_K visit to this k) can form
          ! s = x_{k+1} - x_k, y = g_{k+1} - g_k.
          CALL rdmft_bb_record(orb_bb_pk(ik), orb_iterate_flat, grad_flat)
          line_ok = .FALSE.
          IF (use_orb_wolfe) THEN
             CALL rdmft_orb_ls_setup_block_k(ik, C_save, dir, npw, npwx, nbnd, npol, &
                  gamma_only, gstart, nch, wk(ik), vxpsi, vx_diag)
             CALL rdmft_strong_wolfe_ls(rdmft_orb_ls_eval, etot_start, slope, alpha, &
                  rdmft_line_search_c1, rdmft_line_search_c2, &
                  rdmft_line_search_max_iter, rdmft_line_search_max_zoom, lsres)
             IF (lsres%success .AND. lsres%step > 0.0_DP) THEN
                alpha = lsres%step
                etot_io = lsres%f_new
                line_ok = .TRUE.
             ELSE IF (rdmft_verbose >= 1 .AND. ionode) THEN
                WRITE(stdout, '(7X,A,I0,A)') &
                     'ORB block (ik=', ik, '): strong Wolfe failed; Armijo fallback.'
             ENDIF
          ENDIF
          IF (.NOT. line_ok) THEN
             DO ls = 1, rdmft_line_search_max_iter
                IF (gamma_only) THEN
                   CALL stiefel_retract_gamma(C_save, dir, alpha, npw, nbnd, &
                        npwx*npol, gstart, C_try)
                ELSE
                   CALL stiefel_retract_k(C_save, dir, alpha, npw, nbnd, &
                        npwx*npol, C_try)
                ENDIF
                evc(:, :) = C_try(:, :)
                IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
                CALL rdmft_mark_orbitals_changed()
                CALL rdmft_total_energy(etot_trial)
                IF (etot_trial <= etot_start + rdmft_line_search_c1 * alpha * slope) THEN
                   line_ok = .TRUE.
                   etot_io = etot_trial
                   EXIT
                ENDIF
                IF (rdmft_line_search_polynomial) THEN
                   CALL rdmft_armijo_next_alpha(alpha, etot_start, slope, etot_trial, &
                        rdmft_line_search_rho, alpha)
                ELSE
                   alpha = alpha * rdmft_line_search_rho
                ENDIF
             ENDDO
          ENDIF
          !
          IF (.NOT. line_ok) THEN
             ! Restore this k and reset its optimiser state; move on
             ! to the next k.
             evc(:, :) = C_save(:, :)
             IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
             CALL rdmft_mark_orbitals_changed()
             CALL rdmft_total_energy(etot_io)
             IF (use_lbfgs_orb) CALL rdmft_lbfgs_reset(orb_lbfgs_pk(ik))
             have_prev_pk(ik) = .FALSE.
             IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,I0,A,I0,A,1PE10.2)') &
                  'RDMFT orb inner ', inner, ' (ik=', ik, &
                  ') line search failed; ||G_R^k|| = ', gnorm
             CYCLE
          ENDIF
          !
          ! 4) Update the per-k L-BFGS curvature pair with the
          !    accepted step and the gradient at the new evc[ik].
          IF (use_lbfgs_orb) THEN
             CALL pack_psi_real(C_try - C_save, npwx*npol, nbnd, step_flat)
             CALL rdmft_apply_h_one_psi(npw, nbnd, evc, hpsi)
             gradC = (0.0_DP, 0.0_DP)
             DO ib = 1, nbnd
                gradC(:, ib) = wk(ik) * rdmft_n(ib, ik) * hpsi(:, ib)
             ENDDO
             ! Single-k aceinit_k build (only_k=ik): see the ``only_k``
             ! comment at the gradient build above.  ``O(N_k)`` instead
             ! of ``O(N_k^2)`` for the per-k post-step gradient.
             DO it = 1, nch
                CALL rdmft_compute_xc_channel(it, e_t, vx_diag, vxpsi, only_k=ik)
                DO ib = 1, nbnd
                   CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
                   gradC(:, ib) = gradC(:, ib) + coef * wk(ik) * w * vxpsi(:, ib, ik)
                ENDDO
             ENDDO
             ! Restore evc / current_k after the per-channel single-k build.
             current_k = ik
             IF (lsda) current_spin = isk(ik)
             evc(:, :) = C_try(:, :)
             IF (gamma_only) THEN
                CALL stiefel_project_tangent_gamma(evc, gradC, npw, nbnd, &
                     npwx*npol, gstart, eta)
             ELSE
                CALL stiefel_project_tangent_k(evc, gradC, npw, nbnd, &
                     npwx*npol, eta)
             ENDIF
             CALL pack_psi_real(eta, npwx*npol, nbnd, grad_diff_flat)
             grad_diff_flat = grad_diff_flat - grad_flat
             CALL rdmft_lbfgs_update(orb_lbfgs_pk(ik), step_flat, grad_diff_flat)
          ENDIF
          !
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A,I4,A,I3,A,F18.10,A,1PE10.2,0P,A,1PE10.2)') &
               'RDMFT orb inner ', inner, ' (ik=', ik, ') end: global E = ', etot_io, &
               '  ||G_R^k|| = ', gnorm, '  alpha = ', alpha
          CALL rdmft_inner_log_record('orb', inner, ik, etot_io, gnorm, &
               alpha, 0.0_DP, 0.0_DP)
       ENDDO  ! ik
       IF (.NOT. any_active_k) THEN
          IF (rdmft_verbose >= 1 .AND. ionode) WRITE(stdout, '(7X,A)') &
               'ORB block converged: all k-points stationary (rdmft_orb_grad_tol)'
          IF (PRESENT(orb_converged)) orb_converged = .TRUE.
          EXIT
       ENDIF
    ENDDO  ! inner
    !
    DEALLOCATE(hpsi, gradC, eta, C_save, C_try, dir, vxpsi, vx_diag)
    DEALLOCATE(grad_flat, dir_flat, step_flat, grad_diff_flat)
    IF (ALLOCATED(orb_iterate_flat)) DEALLOCATE(orb_iterate_flat)
    !
  END SUBROUTINE rdmft_orbital_step_block_k
  !
  !-----------------------------------------------------------------
  SUBROUTINE pack_psi_real(psi, lda, nb, vec)
    !---------------------------------------------------------------
    !! Re/Im interleaved packing of a complex (lda, nb) array into
    !! a real vector of length 2*lda*nb.  Used by the L-BFGS view of
    !! the orbital coefficients.
    INTEGER, INTENT(IN) :: lda, nb
    COMPLEX(DP), INTENT(IN) :: psi(lda, nb)
    REAL(DP), INTENT(OUT) :: vec(2 * lda * nb)
    INTEGER :: i, j, p
    p = 0
    DO j = 1, nb
       DO i = 1, lda
          p = p + 1
          vec(p) = REAL(psi(i, j), KIND=DP)
          p = p + 1
          vec(p) = AIMAG(psi(i, j))
       ENDDO
    ENDDO
  END SUBROUTINE pack_psi_real
  !
  !-----------------------------------------------------------------
  SUBROUTINE unpack_psi_real(vec, lda, nb, psi)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN) :: lda, nb
    REAL(DP), INTENT(IN) :: vec(2 * lda * nb)
    COMPLEX(DP), INTENT(OUT) :: psi(lda, nb)
    INTEGER :: i, j, p
    p = 0
    DO j = 1, nb
       DO i = 1, lda
          p = p + 1
          psi(i, j) = CMPLX(vec(p), 0.0_DP, KIND=DP)
          p = p + 1
          psi(i, j) = psi(i, j) + CMPLX(0.0_DP, vec(p), KIND=DP)
       ENDDO
    ENDDO
  END SUBROUTINE unpack_psi_real
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_sort_natocc()
    !---------------------------------------------------------------
    !! Per-k-point permutation of the natural orbitals (and their
    !! occupations) so that ``rdmft_n(1:nbnd, ik)`` is monotonically
    !! non-increasing in the band index ``ib``.  The Stiefel orbital
    !! block produces orbitals that are 1-RDM eigenvectors but NOT
    !! eigenstates of any one-body Hamiltonian, so the natural band
    !! index has no a-priori energy ordering.  Sorting by occupation
    !! puts the occupations in the "expected" descending order
    !! (most-occupied band first).  The total energy is invariant
    !! under this purely-unitary reshuffling.
    USE wvfct,         ONLY : nbnd, npwx
    USE klist,         ONLY : nks
    USE noncollin_module, ONLY : npol
    USE wavefunctions, ONLY : evc
    USE io_files,      ONLY : nwordwfc, iunwfc
    USE buffers,       ONLY : get_buffer, save_buffer
    USE mp,            ONLY : mp_sum
    USE mp_pools,      ONLY : inter_pool_comm
    USE rdmft_module,  ONLY : rdmft_n, rdmft_mark_orbitals_changed
    !
    INTEGER, ALLOCATABLE :: perm(:)
    COMPLEX(DP), ALLOCATABLE :: evc_tmp(:,:)
    REAL(DP), ALLOCATABLE :: n_tmp(:)
    INTEGER :: ik, ib, jb, kmax
    REAL(DP) :: nmax
    LOGICAL :: any_swap
    INTEGER :: any_swap_i
    !
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    any_swap = .FALSE.
    ALLOCATE(perm(nbnd), evc_tmp(npwx*npol, nbnd), n_tmp(nbnd))
    DO ik = 1, nks
       IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       ! Selection sort -- nbnd is small (10s), no need for std::sort
       ! plumbing; deterministic and stable for ties.
       DO ib = 1, nbnd
          perm(ib) = ib
       ENDDO
       DO ib = 1, nbnd - 1
          kmax = ib
          nmax = rdmft_n(perm(ib), ik)
          DO jb = ib + 1, nbnd
             IF (rdmft_n(perm(jb), ik) > nmax) THEN
                kmax = jb
                nmax = rdmft_n(perm(jb), ik)
             ENDIF
          ENDDO
          IF (kmax /= ib) THEN
             ! Swap perm(ib) <-> perm(kmax).
             jb       = perm(ib)
             perm(ib) = perm(kmax)
             perm(kmax) = jb
             any_swap = .TRUE.
          ENDIF
       ENDDO
       ! Apply permutation to (n, evc) at this k-point.
       DO ib = 1, nbnd
          n_tmp(ib)       = rdmft_n(perm(ib), ik)
          evc_tmp(:, ib)  = evc(:, perm(ib))
       ENDDO
       rdmft_n(:, ik) = n_tmp(:)
       evc(:, :)      = evc_tmp(:, :)
       IF (nks > 1) CALL save_buffer(evc, nwordwfc, iunwfc, ik)
    ENDDO
    DEALLOCATE(perm, evc_tmp, n_tmp)
    !
    ! The local ``any_swap`` flag is the OR over this pool's k-points
    ! only.  Two image-pools carry disjoint k-point slices, so each
    ! pool can independently conclude ``any_swap = .FALSE.`` (its slice
    ! is already sorted) while the OTHER pool sees ``.TRUE.`` (its
    ! slice was unsorted).  If we let the flag stay divergent the
    ! downstream post-loop ``rdmft_grad_n`` takes a different branch
    ! on each pool: the pool that set ``rdmft_exxbuff_stale = .TRUE.``
    ! enters ``rdmft_xc_reinit_exxbuff`` -> ``exxinit`` ->
    ! ``poolcollect`` (line 323 of ``exx.f90``), the other skips it
    ! entirely.  The two pools then rendezvous at DIFFERENT
    ! ``poolcollect`` calls and deadlock against each other on
    ! ``inter_pool_comm``.  Reduce the flag across pools (an OR)
    ! before deciding whether to mark the EXX buffer stale so both
    ! pools take the same branch in ``rdmft_xc_compute_exchange_channels``.
    any_swap_i = 0
    IF (any_swap) any_swap_i = 1
    CALL mp_sum(any_swap_i, inter_pool_comm)
    any_swap = (any_swap_i > 0)
    !
    ! Permuting evc invalidates the global EXX real-space buffer
    ! `exxbuff` (which holds an ordered copy of evc indexed by
    ! ib).  Any subsequent vexxace call would use the OLD ordering,
    ! returning meaningless <psi_sorted | V_x_unsorted | psi_sorted>
    ! values.  Mark exxbuff stale so the next ``rdmft_compute_xc_channel``
    ! rebuilds it against the sorted orbitals.
    IF (any_swap) CALL rdmft_mark_orbitals_changed()
    !
  END SUBROUTINE rdmft_sort_natocc
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_print_natural_occupations()
    !---------------------------------------------------------------
    !! Print the final natural occupations (no band diagonals).
    !! Every MPI rank must enter (``poolcollect`` is collective); only
    !! ionode writes to stdout.
    !
    USE io_global,    ONLY : stdout, ionode
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks, nkstot, wk
    USE lsda_mod,     ONLY : lsda, isk
    USE mp,           ONLY : mp_barrier
    USE mp_world,     ONLY : world_comm
    USE rdmft_module, ONLY : rdmft_n
    !
    REAL(DP) :: sum_n
    CHARACTER(LEN=9) :: slbl
    INTEGER  :: ik, ib, ispin
    REAL(DP), ALLOCATABLE :: n_g(:,:)
    REAL(DP), ALLOCATABLE :: wk_loc(:,:), wk_g(:,:), isk_loc(:,:), isk_g(:,:)
    !
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    !
    ! Cross-pool synchronisation (see ``rdmft_occ_eps_log_record`` for
    ! the rationale).  Same one-barrier-on-``world_comm`` pattern: any
    ! rank entering this routine must wait for the slowest image-pool
    ! to finish its ``exxinit_std`` broadcast loop.
    CALL mp_barrier(world_comm)
    !
    ALLOCATE(n_g(nbnd, nkstot))
    ALLOCATE(wk_loc(1, nks), wk_g(1, nkstot), isk_loc(1, nks), isk_g(1, nkstot))
    DO ik = 1, nks
       wk_loc(1, ik)  = wk(ik)
       isk_loc(1, ik) = REAL(isk(ik), DP)
    ENDDO
    CALL poolcollect(nbnd, nks, rdmft_n, nkstot, n_g)
    CALL poolcollect(1, nks, wk_loc,  nkstot, wk_g)
    CALL poolcollect(1, nks, isk_loc, nkstot, isk_g)
    !
    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') &
            'Final natural occupations (sorted by n_ik, descending):'
       WRITE(stdout, '(5X,A)') '  ik    ib   spin      n_ik       wk*n_ik'
       DO ik = 1, nkstot
          ispin = NINT(isk_g(1, ik))
          IF (lsda) THEN
             slbl = MERGE('spin up  ', 'spin down', ispin == 1)
          ELSE
             slbl = 'spinless '
          ENDIF
          sum_n = 0.0_DP
          DO ib = 1, nbnd
             WRITE(stdout, '(5X,2I6,2X,A9,1X,F10.7,3X,F10.7)') &
                  ik, ib, slbl, n_g(ib, ik), wk_g(1, ik) * n_g(ib, ik)
             sum_n = sum_n + wk_g(1, ik) * n_g(ib, ik)
          ENDDO
          WRITE(stdout, '(5X,A,I0,A,F12.8)') &
               'sum_i wk*n_(i,k=', ik, ') = ', sum_n
       ENDDO
    ENDIF
    !
    DEALLOCATE(n_g, wk_loc, wk_g, isk_loc, isk_g)
    !
  END SUBROUTINE rdmft_print_natural_occupations
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_krefine_orbital_polish(orb_maxiter, orb_tol, etot_out)
    !---------------------------------------------------------------
    !! Fixed-occupation orbital minimisation on the current k-mesh.
    !! Used after k-interpolation in ``pw.x`` restart to polish
    !! natural orbitals while keeping ``rdmft\_n`` frozen.
    USE io_global, ONLY : stdout, ionode
    USE rdmft_module, ONLY : rdmft_orb_maxiter, rdmft_occ_maxiter, &
                              rdmft_orb_strategy, rdmft_orb_grad_tol, &
                              rdmft_e_one, rdmft_e_har, rdmft_e_xc, rdmft_e_const
    USE rdmft_energy, ONLY : rdmft_total_energy
    !
    INTEGER, INTENT(IN)  :: orb_maxiter
    REAL(DP), INTENT(IN) :: orb_tol
    REAL(DP), INTENT(OUT) :: etot_out
    !
    INTEGER :: save_occ_max, save_orb_max, it
    CHARACTER(LEN=16) :: save_orb_strategy
    REAL(DP) :: etot_now
    LOGICAL :: orb_conv
    !
    IF (orb_maxiter <= 0) THEN
       CALL rdmft_total_energy(etot_out)
       RETURN
    ENDIF
    !
    save_occ_max = rdmft_occ_maxiter
    save_orb_max = rdmft_orb_maxiter
    save_orb_strategy = rdmft_orb_strategy
    rdmft_occ_maxiter = 0
    rdmft_orb_maxiter = 1
    rdmft_orb_strategy = 'block_k'
    !
    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A,I0,A,ES10.2,A)') &
            'RDMFT k-refine: fixed-n orbital polish (max ', orb_maxiter, &
            ' iters, tol ', orb_tol, ')'
       WRITE(stdout, '(5X,A)') &
            '  E below is the global RDMFT total (all k, incl. Ewald).'
       WRITE(stdout, '(5X,A)') &
            '  RDMFT orb inner lines report the same global E after each k line search.'
    ENDIF
    !
    CALL rdmft_total_energy(etot_now)
    IF (ionode) THEN
       WRITE(stdout, '(5X,A,F18.10,A)') &
            '  E before polish = ', etot_now, ' Ry'
       WRITE(stdout, '(5X,A,F18.10)') &
            '    E_one (T+V_loc+V_NL)       = ', rdmft_e_one
       WRITE(stdout, '(5X,A,F18.10)') &
            '    E_Hartree                  = ', rdmft_e_har
       WRITE(stdout, '(5X,A,F18.10)') &
            '    E_xc_RDMFT                 = ', rdmft_e_xc
       WRITE(stdout, '(5X,A,F18.10)') &
            '    E_Ewald                    = ', rdmft_e_const
    ENDIF
    !
    orb_conv = .FALSE.
    DO it = 1, orb_maxiter
       CALL rdmft_orbital_step(etot_now, orb_conv)
       IF (ionode) WRITE(stdout, '(5X,A,I0,A,F18.10,A)') &
            '  polish iter ', it, ': global E = ', etot_now, ' Ry'
       IF (orb_conv) THEN
          IF (ionode) WRITE(stdout, '(5X,A)') &
               '  orbital block converged (gradient below rdmft_orb_grad_tol)'
          EXIT
       ENDIF
    ENDDO
    !
    rdmft_occ_maxiter = save_occ_max
    rdmft_orb_maxiter = save_orb_max
    rdmft_orb_strategy = save_orb_strategy
    etot_out = etot_now
    !
  END SUBROUTINE rdmft_krefine_orbital_polish
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_wall_time()
    !---------------------------------------------------------------
    !! Wall-clock seconds from the QE ``cclock`` timer.
    INTERFACE
       FUNCTION f_wall() BIND(C, name="cclock") RESULT(t)
          USE ISO_C_BINDING
          REAL(KIND=c_double) :: t
       END FUNCTION f_wall
    END INTERFACE
    rdmft_wall_time = f_wall()
  END FUNCTION rdmft_wall_time
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_report_wall_time(label, t0)
    !---------------------------------------------------------------
    !! Print a block wall time when ``rdmft_verbose >= 1``.  The
    !! elapsed interval is the MPI-image maximum so the reported
    !! value matches the slowest rank.
    USE io_global,  ONLY : stdout, ionode
    USE mp,         ONLY : mp_max
    USE mp_images,  ONLY : intra_image_comm
    USE rdmft_module, ONLY : rdmft_verbose
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(DP), INTENT(IN)         :: t0
    REAL(DP) :: elapsed
    !
    IF (rdmft_verbose < 1) RETURN
    elapsed = rdmft_wall_time() - t0
    CALL mp_max(elapsed, intra_image_comm)
    IF (ionode) WRITE(stdout, '(7X,A,A,F10.2,A)') &
         TRIM(label), ' wall time = ', elapsed, ' s'
  END SUBROUTINE rdmft_report_wall_time
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_report_outer_wall_time(outer, t0)
    !---------------------------------------------------------------
    !! Macro-iteration wall time (``rdmft_verbose >= 1``).
    USE io_global,  ONLY : stdout, ionode
    USE mp,         ONLY : mp_max
    USE mp_images,  ONLY : intra_image_comm
    USE rdmft_module, ONLY : rdmft_verbose
    INTEGER, INTENT(IN)  :: outer
    REAL(DP), INTENT(IN) :: t0
    REAL(DP) :: elapsed
    !
    IF (rdmft_verbose < 1) RETURN
    elapsed = rdmft_wall_time() - t0
    CALL mp_max(elapsed, intra_image_comm)
    IF (ionode) WRITE(stdout, '(7X,A,I0,A,F10.2,A)') &
         'RDMFT outer ', outer, ' wall time = ', elapsed, ' s'
  END SUBROUTINE rdmft_report_outer_wall_time
  !
END MODULE rdmft_solver
