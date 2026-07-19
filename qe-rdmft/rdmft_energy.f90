!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_energy
  !------------------------------------------------------------------
  !! Energy and gradient evaluator for the PWscf RDMFT solver.
  !!
  !! The RDMFT total energy is
  !!
  !!   E[n, psi] = T + E_loc + E_NL + E_Hartree[rho_RDMFT]
  !!             + sum_t coef_t * E_x^HF[gamma_t]
  !!             + E_Ewald
  !!
  !! where the natural orbitals are stored in QE's \texttt{evc} (one
  !! k-point at a time) and the natural occupations are
  !! \texttt{rdmft\_n(i,k)} from \texttt{rdmft\_module}.
  !!
  !! Sign / factor conventions match QE's wg-based bookkeeping:
  !!
  !!   wg_RDMFT(i,k) = degspin * wk(k) * n_ik     (one-body, Hartree)
  !!   wg_t(i,k)     = degspin * wk(k) * w_t(n)   (channel t exchange)
  !!
  !! The exchange Hamiltonians \(V_x^t\) are built by reusing the
  !! existing PWscf EXX-via-ACE machinery: we temporarily overwrite
  !! \texttt{wg(:,:)} with the channel-specific weights, call
  !! \texttt{exxinit} + \texttt{aceinit} so the global \texttt{xi(:,:,:)}
  !! encodes the modified-DM Fock operator, then apply
  !! \texttt{vexxace\_*} to obtain \(V_x^t |\psi\rangle\) for the energy
  !! diagonals and the orbital gradient.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: rdmft_initial_n_from_ks, rdmft_update_density_and_pot
  PUBLIC :: rdmft_compute_one_body_diag, rdmft_compute_xc_channel
  PUBLIC :: rdmft_total_energy, rdmft_grad_n, rdmft_dedn_from_grad_n, &
            rdmft_write_dedn
  PUBLIC :: rdmft_compute_riemannian_gradient
  PUBLIC :: rdmft_apply_h_one_psi, rdmft_refresh_x_occupation
  PUBLIC :: rdmft_exxinit_once
  PUBLIC :: rdmft_set_wg_from_n, rdmft_set_wg_for_channel, rdmft_set_wg_for_exxinit
  PUBLIC :: rdmft_band_energies, rdmft_apply_occ_precond, &
            rdmft_compute_elk_occ_scale, rdmft_compute_tsm_probe_energies, &
            rdmft_build_occ_grad_dir, rdmft_apply_orb_precond, &
            rdmft_collect_band_diagnostics, rdmft_clear_band_diagnostics, &
            rdmft_grad_n_from_cached_diag, rdmft_ensure_exxbuff_tsm_capacity, &
            rdmft_diag_computed, rdmft_diag_h, rdmft_diag_vx, rdmft_diag_tsm
  !
  LOGICAL :: rdmft_diag_computed = .FALSE.
  REAL(DP), ALLOCATABLE :: rdmft_diag_h(:,:), rdmft_diag_vx(:,:), &
                            rdmft_diag_tsm(:,:)
  INTEGER :: rdmft_tsm_probe_ib = 0
  !! When \(>0\), :subroutine:`rdmft_compute_xc_channel` is building
  !! the ACE metric for a TSM occupation probe at band ``ib``; nearly
  !! empty bands at that k-point are dropped from ``wg`` so
  !! ``mexx`` stays well conditioned.
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_initial_n_from_ks()
    !---------------------------------------------------------------
    !! Seed the natural occupations from the converged KS weights.
    !!
    !!   n_ik = wg(i,k) / (degspin * wk(k))
    !!
  !! Optional occupation seed when \texttt{rdmft\_occ\_init\_mode} is
  !! \texttt{perturbed}, \texttt{binary}, or \texttt{uniform} (ABACUS
  !! \texttt{OccInitMode} mirror in \texttt{rdmft\_solver.cpp}):
  !!
  !! * \texttt{perturbed}: $\pm\delta$ on Fermi-window bands when
  !!   $K>0$ and $\delta>0$; legacy global push when $K\le 0$.
  !! * \texttt{binary}: set window bands to $\delta$ or $1-\delta$
  !!   (from KS $n\lessgtr 0.5$) when $K>0$ and $\delta>0$.
  !! * \texttt{uniform}: replace window bands by their mean when $K>0$.
  !! * \texttt{binary}/\texttt{uniform} with $K\le 0$ (or $\delta\le 0$
  !!   for binary) leave KS occupations unchanged.
    !
    USE klist, ONLY : nks, wk, nelec, nelup, neldw, two_fermi_energies
    USE wvfct, ONLY : nbnd, wg
    USE noncollin_module, ONLY : noncolin
    USE lsda_mod, ONLY : lsda
    USE rdmft_module
    USE rdmft_occupation, ONLY : rdmft_proximal_project_occ, &
                                  rdmft_weighted_sum
    !
    IMPLICIT NONE
    !
    INTEGER :: ib, ik, kfound
    REAL(DP) :: degspin, frac, sumw
    CHARACTER(LEN=16) :: occ_mode
    INTEGER :: i
    !
    degspin = 2.0_DP
    IF (noncolin) degspin = 1.0_DP
    rdmft_spin_factor = degspin
    !
    IF (.NOT. ALLOCATED(rdmft_n)) CALL rdmft_allocate(nbnd, nks)
    !
    ! NOTE on spin conventions: for nspin=1 (closed shell) QE bakes the
    ! spin doubling into wk so that wg(i,k) = wk(k) * f_per_spin with
    ! f_per_spin in [0,1].  The natural occupation we evolve is the
    ! same per-spin number.  For nspin=2 each spin channel has its own
    ! wk and we likewise evolve n in [0,1] per band per (k,spin) row.
    ! In all cases the equality target is the total electron count
    ! ``nelec`` because that is what ``sum_k wk * sum_i n_ik`` equals
    ! when wg = wk * n_ik (the convention used by the rest of QE).
    DO ik = 1, nks
       IF (wk(ik) > 0.0_DP) THEN
          DO ib = 1, nbnd
             rdmft_n(ib, ik) = wg(ib, ik) / wk(ik)
          ENDDO
       ELSE
          rdmft_n(:, ik) = 0.0_DP
       ENDIF
    ENDDO
    !
    rdmft_n_target = nelec
    rdmft_fix_magnetization = (lsda .AND. two_fermi_energies)
    IF (rdmft_fix_magnetization) THEN
       rdmft_n_target_up   = nelup
       rdmft_n_target_down = neldw
    ELSE
       rdmft_n_target_up   = 0.0_DP
       rdmft_n_target_down = 0.0_DP
    ENDIF
    !
    occ_mode = ADJUSTL(rdmft_occ_init_mode)
    DO i = 1, LEN_TRIM(occ_mode)
       IF (occ_mode(i:i) >= 'A' .AND. occ_mode(i:i) <= 'Z') &
            occ_mode(i:i) = ACHAR(IACHAR(occ_mode(i:i)) + 32)
    ENDDO
    !
    SELECT CASE (TRIM(occ_mode))
    CASE ('perturbed')
       IF (rdmft_occ_init_perturb > 0.0_DP) THEN
          IF (rdmft_occ_init_nbands_top > 0) THEN
             CALL rdmft_apply_perturbed_occ_fermi_window(nbnd, nks, &
                  rdmft_occ_init_nbands_top, rdmft_occ_init_perturb, rdmft_n)
          ELSE
             ! Legacy QE behaviour when ``rdmft_occ_init_nbands_top <= 0``:
             ! push **every** band off the {0,1} clip walls.  Prefer setting
             ! ``rdmft_occ_init_nbands_top > 0`` to match ABACUS production
             ! inputs (ABACUS LiH_HF uses K=4, delta=1e-3; QE LiH uses K=1).
             DO ik = 1, nks
                DO ib = 1, nbnd
                   IF (rdmft_n(ib, ik) >= 0.5_DP) THEN
                      rdmft_n(ib, ik) = MAX(rdmft_occ_init_perturb, &
                                             rdmft_n(ib, ik) - rdmft_occ_init_perturb)
                   ELSE
                      rdmft_n(ib, ik) = MIN(1.0_DP - rdmft_occ_init_perturb, &
                                             rdmft_n(ib, ik) + rdmft_occ_init_perturb)
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    CASE ('binary')
       IF (rdmft_occ_init_nbands_top > 0 .AND. rdmft_occ_init_perturb > 0.0_DP) THEN
          CALL rdmft_apply_binary_occ_fermi_window(nbnd, nks, &
               rdmft_occ_init_nbands_top, rdmft_occ_init_perturb, rdmft_n)
       ENDIF
    CASE ('uniform')
       IF (rdmft_occ_init_nbands_top > 0) THEN
          CALL rdmft_apply_uniform_occ_fermi_window(nbnd, nks, &
               rdmft_occ_init_nbands_top, rdmft_n)
       ENDIF
    CASE DEFAULT
       ! 'ks' or unknown: keep KS occupations (unknown modes silently
       ! fall back to ks; misspelled keywords are not fatal).
    END SELECT
    !
    ! Project onto the feasible set (clip + dual rescaling).
    CALL rdmft_proximal_project_occ(rdmft_n, wk, nbnd, nks)
    !
  END SUBROUTINE rdmft_initial_n_from_ks
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_apply_perturbed_occ_fermi_window(nbnd, nks, nbands_top_in, delta, nmat)
    !---------------------------------------------------------------
    !! ABACUS ``OccInitMode::Perturbed`` mirror: for each k (spin
    !! channel row) find an approximate Fermi boundary and apply
    !! ``±delta`` to the ``K`` bands immediately above and ``K`` bands
    !! immediately below, with ``K = min(nbands_top_in, nbnd)``.
    INTEGER, INTENT(IN) :: nbnd, nks, nbands_top_in
    REAL(DP), INTENT(IN) :: delta
    REAL(DP), INTENT(INOUT) :: nmat(nbnd, nks)
    !
    INTEGER :: ik, K, it_win, t, ib_fermi, ib_above, ib_below
    !
    DO ik = 1, nks
       K = MIN(nbands_top_in, nbnd)
       ib_fermi = rdmft_find_fermi_boundary_ib(nmat(:, ik), nbnd)
       DO it_win = 1, K
          t = it_win - 1
          ib_above = ib_fermi + 1 + t
          IF (ib_above >= 1 .AND. ib_above <= nbnd) THEN
             IF (nmat(ib_above, ik) < 0.5_DP) THEN
                nmat(ib_above, ik) = nmat(ib_above, ik) + delta
             ELSE
                nmat(ib_above, ik) = nmat(ib_above, ik) - delta
             ENDIF
          ENDIF
          ib_below = ib_fermi - t
          IF (ib_below >= 1 .AND. ib_below <= nbnd) THEN
             IF (nmat(ib_below, ik) < 0.5_DP) THEN
                nmat(ib_below, ik) = nmat(ib_below, ik) + delta
             ELSE
                nmat(ib_below, ik) = nmat(ib_below, ik) - delta
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_apply_perturbed_occ_fermi_window
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_apply_binary_occ_fermi_window(nbnd, nks, nbands_top_in, delta, nmat)
    !---------------------------------------------------------------
    !! ABACUS ``OccInitMode::Binary`` mirror: Fermi-window bands with
    !! KS $n<0.5$ are set to ``delta``, bands with $n\ge 0.5$ to
    !! ``1-delta``.
    INTEGER, INTENT(IN) :: nbnd, nks, nbands_top_in
    REAL(DP), INTENT(IN) :: delta
    REAL(DP), INTENT(INOUT) :: nmat(nbnd, nks)
    !
    INTEGER :: ik, K, it_win, t, ib_fermi, ib_above, ib_below
    REAL(DP) :: n_ks
    REAL(DP), ALLOCATABLE :: n_ks_col(:)
    !
    ALLOCATE(n_ks_col(nbnd))
    DO ik = 1, nks
       K = MIN(nbands_top_in, nbnd)
       n_ks_col = nmat(:, ik)
       ib_fermi = rdmft_find_fermi_boundary_ib(n_ks_col, nbnd)
       DO it_win = 1, K
          t = it_win - 1
          ib_above = ib_fermi + 1 + t
          IF (ib_above >= 1 .AND. ib_above <= nbnd) THEN
             n_ks = n_ks_col(ib_above)
             IF (n_ks < 0.5_DP) THEN
                nmat(ib_above, ik) = delta
             ELSE
                nmat(ib_above, ik) = 1.0_DP - delta
             ENDIF
          ENDIF
          ib_below = ib_fermi - t
          IF (ib_below >= 1 .AND. ib_below <= nbnd) THEN
             n_ks = n_ks_col(ib_below)
             IF (n_ks < 0.5_DP) THEN
                nmat(ib_below, ik) = delta
             ELSE
                nmat(ib_below, ik) = 1.0_DP - delta
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DEALLOCATE(n_ks_col)
    !
  END SUBROUTINE rdmft_apply_binary_occ_fermi_window
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_apply_uniform_occ_fermi_window(nbnd, nks, nbands_top_in, nmat)
    !---------------------------------------------------------------
    !! ABACUS ``OccInitMode::Uniform`` mirror: average the KS
    !! occupations over the Fermi window and assign the mean to each
    !! selected band.
    INTEGER, INTENT(IN) :: nbnd, nks, nbands_top_in
    REAL(DP), INTENT(INOUT) :: nmat(nbnd, nks)
    !
    INTEGER :: ik, K, it_win, t, ib_fermi, ib_above, ib_below, n_sel
    REAL(DP) :: n_top, n_uniform
    !
    DO ik = 1, nks
       K = MIN(nbands_top_in, nbnd)
       ib_fermi = rdmft_find_fermi_boundary_ib(nmat(:, ik), nbnd)
       n_top = 0.0_DP
       n_sel = 0
       DO it_win = 1, K
          t = it_win - 1
          ib_above = ib_fermi + 1 + t
          IF (ib_above >= 1 .AND. ib_above <= nbnd) THEN
             n_top = n_top + nmat(ib_above, ik)
             n_sel = n_sel + 1
          ENDIF
          ib_below = ib_fermi - t
          IF (ib_below >= 1 .AND. ib_below <= nbnd) THEN
             n_top = n_top + nmat(ib_below, ik)
             n_sel = n_sel + 1
          ENDIF
       ENDDO
       IF (n_sel > 0) THEN
          n_uniform = n_top / REAL(n_sel, DP)
          DO it_win = 1, K
             t = it_win - 1
             ib_above = ib_fermi + 1 + t
             IF (ib_above >= 1 .AND. ib_above <= nbnd) THEN
                nmat(ib_above, ik) = n_uniform
             ENDIF
             ib_below = ib_fermi - t
             IF (ib_below >= 1 .AND. ib_below <= nbnd) THEN
                nmat(ib_below, ik) = n_uniform
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    !
  END SUBROUTINE rdmft_apply_uniform_occ_fermi_window
  !
  !-----------------------------------------------------------------
  INTEGER FUNCTION rdmft_find_fermi_boundary_ib(n_col, nbnd)
    !---------------------------------------------------------------
    !! First band (1-based) with ``n < 0.5`` implies boundary ``ib-1``;
    !! otherwise the last band with ``n > 1e-8``.  Matches ABACUS
    !! ``find_fermi_boundary_index`` up to the 0-based / 1-based offset.
    INTEGER, INTENT(IN) :: nbnd
    REAL(DP), INTENT(IN) :: n_col(nbnd)
    INTEGER :: ib, last_occ
    !
    rdmft_find_fermi_boundary_ib = 0
    DO ib = 1, nbnd
       IF (n_col(ib) < 0.5_DP) THEN
          rdmft_find_fermi_boundary_ib = ib - 1
          RETURN
       ENDIF
    ENDDO
    last_occ = 0
    DO ib = 1, nbnd
       IF (n_col(ib) > 1.0e-8_DP) last_occ = ib
    ENDDO
    rdmft_find_fermi_boundary_ib = last_occ
    !
  END FUNCTION rdmft_find_fermi_boundary_ib
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_set_wg_from_n(weight_kind)
    !---------------------------------------------------------------
    !! Overwrite QE's \texttt{wg} so the rest of the SCF infrastructure
    !! sees the RDMFT weights.  ``weight_kind`` selects the weighting:
    !!
    !!   1 -> wg(i,k) = degspin * wk(k) * n_ik         (rho / Hartree)
    !!   2 -> wg(i,k) = degspin * wk(k) * g(n_ik)      (single-channel exchange)
    !
    USE klist,    ONLY : nks, wk
    USE wvfct,    ONLY : nbnd, wg
    USE noncollin_module, ONLY : noncolin
    USE rdmft_module
    USE rdmft_xc, ONLY : rdmft_g
    !
    INTEGER, INTENT(IN) :: weight_kind
    REAL(DP) :: degspin, w
    INTEGER  :: ib, ik
    !
    degspin = 2.0_DP
    IF (noncolin) degspin = 1.0_DP
    !
    DO ik = 1, nks
       DO ib = 1, nbnd
          SELECT CASE (weight_kind)
          CASE (1)
             w = wk(ik) * rdmft_n(ib, ik)
          CASE (2)
             w = wk(ik) * rdmft_g(rdmft_n(ib, ik))
          CASE DEFAULT
             w = wg(ib, ik)
          END SELECT
          wg(ib, ik) = w
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_set_wg_from_n
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_set_wg_for_channel(it)
    !---------------------------------------------------------------
    !! Thin wrapper around :func:`rdmft_xc_set_wg_for_channel`.
    USE rdmft_xc, ONLY : rdmft_xc_set_wg_for_channel
    INTEGER, INTENT(IN) :: it
    CALL rdmft_xc_set_wg_for_channel(it)
  END SUBROUTINE rdmft_set_wg_for_channel
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_set_wg_for_exxinit()
    !---------------------------------------------------------------
    !! Thin wrapper around :func:`rdmft_xc_set_exx_wg`.
    USE rdmft_xc, ONLY : rdmft_xc_set_exx_wg
    CALL rdmft_xc_set_exx_wg()
  END SUBROUTINE rdmft_set_wg_for_exxinit
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_reinit_exxbuff()
    !---------------------------------------------------------------
    !! Thin wrapper around :func:`rdmft_xc_reinit_exxbuff`.
    USE rdmft_xc, ONLY : rdmft_xc_reinit_exxbuff
    CALL rdmft_xc_reinit_exxbuff()
  END SUBROUTINE rdmft_reinit_exxbuff
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_update_density_and_pot(ehart_out, etxc_out, vtxc_out)
    !---------------------------------------------------------------
    !! With the current natural occupations and natural orbitals,
    !! rebuild the electron density via \texttt{sum\_band} and the
    !! local potential via \texttt{v\_of\_rho} / \texttt{set\_vrs}.
    !!
    !! Returns the Hartree, exchange-correlation and integrated XC
    !! potential energies of the RDMFT density (V_xc here is the KS
    !! semilocal V_xc; the RDMFT Fock-like piece is added separately).
    !
    USE wvfct,  ONLY : nbnd
    USE klist,  ONLY : nks
    USE scf,    ONLY : rho, rho_core, rhog_core, v, vltot, vrs, kedtau
    USE gvecs,  ONLY : doublegrid
    USE fft_base, ONLY : dfftp
    USE lsda_mod, ONLY : nspin
    USE ldaU,   ONLY : eth
    USE extfield, ONLY : etotefield
    !
    REAL(DP), INTENT(OUT) :: ehart_out, etxc_out, vtxc_out
    REAL(DP) :: charge
    !
    CALL rdmft_set_wg_from_n(1)
    CALL sum_band()
    CALL v_of_rho(rho, rho_core, rhog_core, &
                  ehart_out, etxc_out, vtxc_out, eth, etotefield, charge, v)
    CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
    !
  END SUBROUTINE rdmft_update_density_and_pot
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_apply_h_one_psi(npw, m, psi, hpsi)
    !---------------------------------------------------------------
    !! Apply the KS one-body+Hartree(+semilocal-V_xc) Hamiltonian to a
    !! batch of orbitals, with EXX **disabled** (so any RDMFT exchange
    !! contribution will be added separately by the caller).
    !!
    !! Uses QE's standard \texttt{h\_psi}.  Caller must already have
    !! set the proper k-point context (current_k, current_spin,
    !! init_us_2, etc.) so this routine just forwards to h_psi.
    USE wvfct, ONLY : npwx
    USE noncollin_module, ONLY : npol
    USE xc_lib, ONLY : exx_is_active, stop_exx, start_exx, xclib_dft_is
    !
    INTEGER, INTENT(IN) :: npw, m
    COMPLEX(DP), INTENT(IN)  :: psi(npwx*npol, m)
    COMPLEX(DP), INTENT(OUT) :: hpsi(npwx*npol, m)
    LOGICAL :: was_active
    !
    was_active = exx_is_active()
    IF (was_active) CALL stop_exx()
    CALL h_psi(npwx, npw, m, psi, hpsi)
    IF (was_active) CALL start_exx()
    !
  END SUBROUTINE rdmft_apply_h_one_psi
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_exxinit_once()
    !---------------------------------------------------------------
    !! Thin wrapper around :func:`rdmft_xc_exxinit_once`.
    USE rdmft_xc, ONLY : rdmft_xc_exxinit_once
    CALL rdmft_xc_exxinit_once()
  END SUBROUTINE rdmft_exxinit_once
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_refresh_x_occupation()
    !---------------------------------------------------------------
    !! Thin wrapper around :func:`rdmft_xc_refresh_x_occupation`.
    USE rdmft_xc, ONLY : rdmft_xc_refresh_x_occupation
    CALL rdmft_xc_refresh_x_occupation()
  END SUBROUTINE rdmft_refresh_x_occupation
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_one_body_diag(ik, npw, h_diag)
    !---------------------------------------------------------------
    !! Per-band \(\langle\psi_i|H_{\mathrm{one}}+V_H+V_{xc}^{KS}|\psi_i\rangle\)
    !! at k-point \texttt{ik}.  Excludes any EXX contribution.
    !!
    !! Sets the k-point context (current_k, current_spin, vkb,
    !! evc-from-buffer) before calling \texttt{h\_psi}, so this routine
    !! can be invoked stand-alone for any k-point in the pool.
    USE wvfct, ONLY : nbnd, npwx, current_k
    USE noncollin_module, ONLY : npol
    USE wavefunctions, ONLY : evc
    USE control_flags, ONLY : gamma_only
    USE gvect, ONLY : gstart
    USE klist, ONLY : igk_k, xk, nks
    USE lsda_mod, ONLY : lsda, current_spin, isk
    USE uspp, ONLY : okvan, vkb, nkb
    USE uspp_init, ONLY : init_us_2
    USE io_files, ONLY : nwordwfc, iunwfc
    USE buffers, ONLY : get_buffer
    USE mp,        ONLY : mp_sum
    USE mp_bands,  ONLY : intra_bgrp_comm
    !
    INTEGER, INTENT(IN)  :: ik, npw
    REAL(DP), INTENT(OUT) :: h_diag(nbnd)
    !
    COMPLEX(DP), ALLOCATABLE :: hpsi(:,:)
    INTEGER :: ib, ig
    !
    current_k = ik
    IF (lsda) current_spin = isk(ik)
    IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
    ! Update g2kin = (k+G)^2 for the current k-point (h_psi reads it
    ! through USE wvfct, ONLY : g2kin).  Without this, multi-k h_psi
    ! uses the previous k's kinetic operator and the per-band
    ! diagonals are wrong.
    CALL g2_kin(ik)
    IF (nkb > 0) CALL init_us_2(npw, igk_k(1, ik), xk(:, ik), vkb)
    !
    ALLOCATE(hpsi(npwx*npol, nbnd))
    hpsi = (0.0_DP, 0.0_DP)
    CALL rdmft_apply_h_one_psi(npw, nbnd, evc, hpsi)
    !
    ! Per-band inner products on this rank's G-vector slice.  Band
    ! parallelism (-nband > 1) is NOT supported: every band group
    ! computes the full diagonal redundantly.  See the README in
    ! PW/src/rdmft/ for the recommended parallel layout.
!$omp parallel do default(shared) private(ib)
    DO ib = 1, nbnd
       IF (gamma_only) THEN
          h_diag(ib) = 2.0_DP * REAL(SUM(CONJG(evc(1:npw, ib)) * hpsi(1:npw, ib)), KIND=DP)
          IF (gstart == 2) THEN
             h_diag(ib) = h_diag(ib) - REAL(evc(1, ib), KIND=DP) * REAL(hpsi(1, ib), KIND=DP)
          ENDIF
       ELSE
          h_diag(ib) = REAL(SUM(CONJG(evc(1:npw, ib)) * hpsi(1:npw, ib)), KIND=DP)
          IF (npol == 2) THEN
             h_diag(ib) = h_diag(ib) + REAL(SUM(CONJG(evc(npwx+1:npwx+npw, ib)) * &
                                                   hpsi(npwx+1:npwx+npw, ib)), KIND=DP)
          ENDIF
       ENDIF
    ENDDO
!$omp end parallel do
    ! Plane-wave parallelisation: each MPI rank in the band group has
    ! only its slice of (igk_k, evc, hpsi).  Sum the partial inner
    ! products across intra_bgrp_comm so every rank ends up with the
    ! complete per-band diagonal.
    CALL mp_sum(h_diag, intra_bgrp_comm)
    DEALLOCATE(hpsi)
    !
  END SUBROUTINE rdmft_compute_one_body_diag
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_xc_channel(it, e_t, vx_diag, vxpsi, only_k, ib_only)
    !---------------------------------------------------------------
    !! Thin wrapper around :func:`rdmft_xc_compute_exchange` (unified XC API).
    USE wvfct, ONLY : nbnd
    USE klist, ONLY : nks
    USE rdmft_xc, ONLY : rdmft_xc_compute_exchange
    !
    INTEGER, INTENT(IN)   :: it
    REAL(DP), INTENT(OUT) :: e_t
    REAL(DP), INTENT(OUT) :: vx_diag(nbnd, nks)
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: vxpsi(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: only_k, ib_only
    !
    IF (PRESENT(vxpsi)) THEN
       IF (PRESENT(only_k)) THEN
          IF (PRESENT(ib_only)) THEN
             CALL errore('rdmft_compute_xc_channel', &
                  'ib_only is incompatible with vxpsi output', 1)
          ELSE
             CALL rdmft_xc_compute_exchange(it, e_t, vx_diag, vxpsi=vxpsi, only_k=only_k)
          ENDIF
       ELSE
          CALL rdmft_xc_compute_exchange(it, e_t, vx_diag, vxpsi=vxpsi)
       ENDIF
    ELSEIF (PRESENT(only_k)) THEN
       IF (PRESENT(ib_only)) THEN
          CALL rdmft_xc_compute_exchange(it, e_t, vx_diag, only_k=only_k, ib_only=ib_only)
       ELSE
          CALL rdmft_xc_compute_exchange(it, e_t, vx_diag, only_k=only_k)
       ENDIF
    ELSE
       CALL rdmft_xc_compute_exchange(it, e_t, vx_diag)
    ENDIF
  END SUBROUTINE rdmft_compute_xc_channel
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_total_energy(etot_out)
    !---------------------------------------------------------------
    !! Evaluate the total RDMFT energy at the current (n, evc).
    USE rdmft_module
    USE wvfct,             ONLY : nbnd, wg
    USE klist,             ONLY : nks, ngk, wk
    USE control_flags,     ONLY : gamma_only
    USE constants,         ONLY : K_BOLTZMANN_RY
    USE ener,              ONLY : ewld, etxc, ehart, vtxc
    USE mp,                ONLY : mp_sum
    USE mp_pools,          ONLY : inter_pool_comm
    USE rdmft_xc,          ONLY : rdmft_binary_entropy_f, RDMFT_XC_HF, &
                                  rdmft_xc_total_energy
    !
    REAL(DP), INTENT(OUT) :: etot_out
    REAL(DP) :: e_one_total, e_xc_total, ehart_loc, etxc_loc, vtxc_loc
    REAL(DP) :: e_entropy, kbt
    REAL(DP), ALLOCATABLE :: h_diag_k(:)
    INTEGER  :: ik, ib, npw
    !
    !!
    !!   E = T + E_loc + E_NL + ehart + sum_t coef_t * E_x^HF[gamma_t]
    !!     + ewld
    !!
    !! Uses \texttt{rdmft\_update\_density\_and\_pot} (which calls
    !! \texttt{sum\_band} + \texttt{v\_of\_rho}) to refresh \texttt{ehart}
    !! and the local potential, then \texttt{rdmft\_compute\_one\_body\_diag}
    !! to read off the kinetic + local + NL piece per band, and finally
    !! \texttt{rdmft\_compute\_xc\_channel} for each XC channel.
    !
    ! NOTE: rdmft_compute_xc_channel already returns the full
    ! <psi|Vx|psi> diagonal (with the gamma-trick correction baked in),
    ! and rdmft_compute_one_body_diag does the same for h_diag, so the
    ! caller-side energy sum is just (1/2) sum wg <psi|Vx|psi> without
    ! an extra gamma-only doubling.
    !
    ALLOCATE(h_diag_k(nbnd))
    !
    ! Update density and KS Hartree+V_xc with current wg = wg_RDMFT.
    CALL rdmft_update_density_and_pot(ehart_loc, etxc_loc, vtxc_loc)
    !
    ! One-body+H+V_xc(KS) energy contribution.  We then **subtract** the
    ! KS V_xc and replace it by the RDMFT exchange below.  We also need
    ! to subtract the double-counted Hartree to get T+E_loc+E_NL alone.
    e_one_total = 0.0_DP
    DO ik = 1, nks
       npw = ngk(ik)
       CALL rdmft_compute_one_body_diag(ik, npw, h_diag_k)
       DO ib = 1, nbnd
          e_one_total = e_one_total + wg(ib, ik) * h_diag_k(ib)
       ENDDO
    ENDDO
    ! Inter-pool reduction: with k-point pools each pool owns a
    ! subset of k-points, so e_one_total is only the per-pool
    ! partial sum.  Combine across pools so every rank gets the
    ! total over all k-points.  The reduction is a no-op when
    ! npool = 1.
    CALL mp_sum(e_one_total, inter_pool_comm)
    !
    ! eband = T+E_loc+E_NL + 2*ehart + int rho V_xc.
    ! Subtract 2*ehart (will add back ehart) and int rho V_xc
    ! (replaced by the RDMFT Fock exchange below).
    ! ehart_loc and vtxc_loc are already global (computed by
    ! v_of_rho, which performs its own intra_bgrp_comm and
    ! inter_pool_comm reductions).
    e_one_total = e_one_total - 2.0_DP * ehart_loc - vtxc_loc
    !
    ! RDMFT exchange (delegated to rdmft_xc).
    CALL rdmft_xc_total_energy(e_xc_total)
    !
    e_entropy = 0.0_DP
    IF (rdmft_xc_id == RDMFT_XC_HF .AND. rdmft_temp > 0.0_DP) THEN
       kbt = K_BOLTZMANN_RY * rdmft_temp
       DO ik = 1, nks
          DO ib = 1, nbnd
             e_entropy = e_entropy + kbt * wk(ik) &
                         * rdmft_binary_entropy_f(rdmft_n(ib, ik))
          ENDDO
       ENDDO
       CALL mp_sum(e_entropy, inter_pool_comm)
    ENDIF
    !
    rdmft_e_one  = e_one_total
    rdmft_e_har  = ehart_loc
    rdmft_e_xc   = e_xc_total
    rdmft_e_entropy = e_entropy
    rdmft_e_const = ewld
    rdmft_etot   = e_one_total + ehart_loc + e_xc_total + e_entropy + ewld
    etot_out = rdmft_etot
    !
    DEALLOCATE(h_diag_k)
    !
  END SUBROUTINE rdmft_total_energy
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_grad_n(grad_n, only_k, ib_only)
    !---------------------------------------------------------------
    !! Compute the physical occupation gradient
    !! \(g_{ik}=\partial E/\partial n_{ik}\) (Section 9.4 of
    !! ``doc/rdmft_calculation.md``):
    !!
    !!   g_ik = w_k * h_diag(i,k)
    !!        + sum_t c_t * w_k * w_t'(n_ik) * <psi|Vx_t|psi>
    !!
    !! \(h_{\mathrm{diag}}\) here is the per-band diagonal of
    !! ``T + V_loc + V_NL + V_H[rho_RDMFT]`` (NOT including the KS
    !! V_xc, which would double-count when added on top of the RDMFT
    !! Fock-exchange channels).
    !!
    !! With the optional ``only_k`` argument, only the gradient column
    !! ``grad_n(:, only_k)`` is evaluated (the exchange is built at that
    !! k-point only).  This avoids an all-k ACE rebuild when only one
    !! k-point column of the occupation gradient is needed.
    !!
    !! With both ``only_k`` and ``ib_only``, only the single element
    !! ``grad_n(ib_only, only_k)`` is built and the underlying ``vexx``
    !! call uses ``m = 1`` -- the TSM-probe fast path, ``nbnd``-fold
    !! cheaper than building the full per-band Fock diagonal.
    !
    USE rdmft_module
    USE wvfct,            ONLY : nbnd
    USE klist,            ONLY : nks, ngk, wk
    USE noncollin_module, ONLY : noncolin
    USE constants,        ONLY : K_BOLTZMANN_RY
    USE rdmft_xc,         ONLY : rdmft_binary_entropy_dfdn, RDMFT_XC_HF, &
                                  rdmft_xc_add_occ_gradient
    USE ener,             ONLY : etxc, vtxc
    USE scf,              ONLY : rho, v
    USE mp,               ONLY : mp_barrier
    USE mp_world,         ONLY : world_comm
    !
    REAL(DP), INTENT(OUT) :: grad_n(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: only_k
    INTEGER, INTENT(IN), OPTIONAL :: ib_only
    !! When present (must combine with ``only_k``), only
    !! ``grad_n(ib_only, only_k)`` is built.  Other entries of
    !! ``grad_n`` are returned as zero.
    REAL(DP), ALLOCATABLE :: h_diag(:)
    REAL(DP) :: degspin, kbt
    INTEGER  :: ik, ib, npw, klo, khi, ib_lo, ib_hi
    !
    IF (PRESENT(ib_only) .AND. .NOT. PRESENT(only_k)) &
         CALL errore('rdmft_grad_n', 'ib_only requires only_k', 1)
    !
    degspin = 2.0_DP
    IF (noncolin) degspin = 1.0_DP
    !
    IF (PRESENT(only_k)) THEN
       klo = only_k
       khi = only_k
    ELSE
       klo = 1
       khi = nks
    ENDIF
    IF (PRESENT(ib_only)) THEN
       ib_lo = ib_only
       ib_hi = ib_only
    ELSE
       ib_lo = 1
       ib_hi = nbnd
    ENDIF
    !
    ALLOCATE(h_diag(nbnd))
    grad_n = 0.0_DP
    !
    DO ik = klo, khi
       npw = ngk(ik)
       ! ``rdmft_compute_one_body_diag`` returns the full per-band
       ! one-body diagonal at ``ik``.  Even in the single-band fast
       ! path we run it as-is: it has no ``vexx`` work, scales as
       ! ``O(nbnd * npw)``, and is vastly cheaper than the Fock kernel.
       CALL rdmft_compute_one_body_diag(ik, npw, h_diag)
       DO ib = ib_lo, ib_hi
          grad_n(ib, ik) = grad_n(ib, ik) + wk(ik) * h_diag(ib)
       ENDDO
    ENDDO
    !
    IF (PRESENT(only_k)) THEN
       IF (PRESENT(ib_only)) THEN
          CALL rdmft_xc_add_occ_gradient(grad_n, only_k=only_k, ib_only=ib_only)
       ELSE
          CALL rdmft_xc_add_occ_gradient(grad_n, only_k=only_k)
       ENDIF
    ELSE
       CALL rdmft_xc_add_occ_gradient(grad_n)
    ENDIF
    !
    IF (rdmft_xc_id == RDMFT_XC_HF .AND. rdmft_temp > 0.0_DP) THEN
       kbt = K_BOLTZMANN_RY * rdmft_temp
       DO ik = klo, khi
          DO ib = ib_lo, ib_hi
             grad_n(ib, ik) = grad_n(ib, ik) + kbt * wk(ik) &
                  * rdmft_binary_entropy_dfdn(rdmft_n(ib, ik))
          ENDDO
       ENDDO
    ENDIF
    !
    DEALLOCATE(h_diag)
    !
    ! The exchange-coupling rebuild in ``rdmft_xc_add_occ_gradient``
    ! calls ``exxinit_std`` which broadcasts ``exxbuff`` in a sequential
    ! ``DO ikq`` loop on ``intra_orthopool_comm``.  The number of
    ! ``ikq`` iterations and the per-broadcast payload are exactly
    ! the same on every rank of an image-pool (they depend on the
    ! k-point / q-point partition which is the same across pools),
    ! so the intra-pool broadcast time is well balanced.  However,
    ! the wall-time to reach this point differs across image-pools
    ! because some pools carry a heavier k-point slice: in those
    ! pools ``rdmft_compute_one_body_diag`` runs longer before
    ! ``rdmft_xc_add_occ_gradient`` is even entered, so the pool
    ! finishes the broadcast noticeably later than its peers.
    ! Without a global barrier here, a fast pool would race ahead to
    ! the next ``poolcollect`` (``inter_pool_comm`` allreduce) and
    ! deadlock against its still-busy ``inter_pool_comm`` partner.
    ! The barrier on ``world_comm`` is the cheapest point that
    ! guarantees uniform call frames across all pools regardless
    ! of the ``-nk`` / ``-ni`` image partition.
    CALL mp_barrier(world_comm)
    !
  END SUBROUTINE rdmft_grad_n
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_dedn_from_grad_n(grad_n, wk, nbnd, nks, dedn)
    !---------------------------------------------------------------
    !! Per-state ``dedn`` in the ELK ``RDM_DEDN.OUT`` convention.
    !!
    !! ``rdmft_grad_n`` returns ``grad_n = w_k (\partial E/\partial n)``.
    !! This routine returns ``dedn_{ik} = grad_n_{ik}/w_k`` so files and
    !! logs match ELK without manual rescaling.
    REAL(DP), INTENT(IN)  :: grad_n(:,:), wk(:)
    INTEGER,  INTENT(IN)  :: nbnd, nks
    REAL(DP), INTENT(OUT) :: dedn(:,:)
    !
    INTEGER :: ib, ik
    !
    dedn = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          IF (ABS(wk(ik)) > 0.0_DP) &
               dedn(ib, ik) = grad_n(ib, ik) / wk(ik)
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_dedn_from_grad_n
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_write_dedn(fprefix, grad_g, occ_g, wk_g, xk_g)
    !---------------------------------------------------------------
    !! Write ELK ``RDM_DEDN.OUT``: per k-point, state index,
    !! occupation, and ``dedn = dE/dn`` (third column).
    USE io_global, ONLY : ionode
    !
    CHARACTER(LEN=*), INTENT(IN) :: fprefix
    REAL(DP), INTENT(IN) :: grad_g(:,:), occ_g(:,:), wk_g(:), xk_g(:,:)
    !
    INTEGER, EXTERNAL :: find_free_unit
    INTEGER :: iunit, ik, ib, nb, nk
    CHARACTER(LEN=256) :: fname
    REAL(DP), ALLOCATABLE :: dedn_g(:,:)
    !
    IF (.NOT. ionode) RETURN
    nb = SIZE(grad_g, 1)
    nk = SIZE(wk_g)
    ALLOCATE(dedn_g(nb, nk))
    CALL rdmft_dedn_from_grad_n(grad_g, wk_g, nb, nk, dedn_g)
    iunit = find_free_unit()
    fname = TRIM(fprefix) // '.RDM_DEDN.OUT'
    OPEN(UNIT=iunit, FILE=TRIM(fname), STATUS='UNKNOWN', FORM='FORMATTED')
    WRITE(iunit, '(I6,A)') nk, ' : nkpt'
    WRITE(iunit, '(I6,A)') nb, ' : nbnd'
    DO ik = 1, nk
       WRITE(iunit, *)
       WRITE(iunit, '(I6,3G18.10,A)') ik, xk_g(1, ik), xk_g(2, ik), xk_g(3, ik), &
            ' : k-point, xk'
       WRITE(iunit, '(A)') '     (state, occupancy and derivative below)'
       DO ib = 1, nb
          WRITE(iunit, '(I6,2G18.10)') ib, occ_g(ib, ik), dedn_g(ib, ik)
       ENDDO
    ENDDO
    CLOSE(iunit)
    DEALLOCATE(dedn_g)
    !
  END SUBROUTINE rdmft_write_dedn
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_grad_n_from_cached_diag(grad_n)
    !---------------------------------------------------------------
    !! Reconstruct \(dE/dn_{ik}\) from cached one-body and XC
    !! diagonals (no ACE rebuild).  Valid for single-channel
    !! functionals (HF / Muller / Power / GU) after
    !! :subroutine:`rdmft_collect_band_diagnostics`.
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks, wk
    USE rdmft_module, ONLY : rdmft_n, rdmft_xc_id
    USE rdmft_xc,     ONLY : rdmft_xc_n_channels, rdmft_xc_add_occ_gradient_from_diag
    !
    REAL(DP), INTENT(OUT) :: grad_n(:,:)
    !
    INTEGER :: ik, ib, nch
    !
    IF (.NOT. rdmft_diag_computed .OR. .NOT. ALLOCATED(rdmft_diag_h) &
        .OR. .NOT. ALLOCATED(rdmft_diag_vx)) &
         CALL errore('rdmft_grad_n_from_cached_diag', &
              'band diagnostics not cached', 1)
    nch = rdmft_xc_n_channels()
    IF (nch /= 1) &
         CALL errore('rdmft_grad_n_from_cached_diag', &
              'only implemented for single-channel XC functionals', nch)
    !
    ! Reconstruct dE/dn from the cached one-body and Vx band diagonals
    ! (collected at the current n by rdmft_collect_band_diagnostics).
    ! Only used post-convergence (DOS / Koopmans-Fock), where the cached
    ! diagonals correspond to the current occupations, so this is exact.
    grad_n = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          grad_n(ib, ik) = wk(ik) * rdmft_diag_h(ib, ik)
       ENDDO
    ENDDO
    CALL rdmft_xc_add_occ_gradient_from_diag(grad_n, rdmft_diag_vx)
    !
  END SUBROUTINE rdmft_grad_n_from_cached_diag
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ensure_exxbuff_tsm_capacity()
    !---------------------------------------------------------------
    !! Thin wrapper around :func:`rdmft_xc_ensure_exxbuff_tsm_capacity`.
    USE rdmft_xc, ONLY : rdmft_xc_ensure_exxbuff_tsm_capacity
    CALL rdmft_xc_ensure_exxbuff_tsm_capacity()
  END SUBROUTINE rdmft_ensure_exxbuff_tsm_capacity
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_riemannian_gradient(G_R_all, gnorm2)
    !---------------------------------------------------------------
    !! Ambient gradient projected to the Stiefel tangent at each k,
    !! and the squared canonical norm \(\sum_k \langle G_R^k,G_R^k\rangle_S\).
    !!
    !! Cost: per call this routine should scale as ``O(N_k^2)`` (the
    !! cost of one all-k :subroutine:`exxinit` + :subroutine:`aceinit`
    !! pair, plus an ``O(N_k)`` :subroutine:`vexxace_*` sweep).  In an
    !! earlier version it scaled as ``O(N_k^3)`` instead, because the
    !! per-channel :subroutine:`rdmft_compute_xc_channel` call was
    !! placed INSIDE the per-k Stiefel projection loop -- so the full
    !! all-k ACE projector was rebuilt :math:`N_k` times per gradient,
    !! turning the ``O(N_k^2)`` aceinit into ``O(N_k^3)``.  The
    !! channel-summed ``Vx|psi>`` is now pre-computed once for all k
    !! before the per-k loop and reused inside it (pure memory
    !! gather), restoring the proper ``O(N_k^2)`` scaling.  See also
    !! the analogous hoist in :subroutine:`rdmft_run_joint`,
    !! :subroutine:`rdmft_orbital_step_joint` and the BLOCK_K
    !! gradient builders in :file:`rdmft_solver.f90`.
    USE wvfct,             ONLY : nbnd, npwx, current_k
    USE klist,             ONLY : nks, ngk, wk, igk_k, xk
    USE rdmft_module,      ONLY : rdmft_n
    USE wavefunctions,     ONLY : evc
    USE noncollin_module,  ONLY : npol
    USE control_flags,     ONLY : gamma_only
    USE gvect,             ONLY : gstart
    USE io_files,          ONLY : nwordwfc, iunwfc
    USE buffers,           ONLY : get_buffer
    USE lsda_mod,          ONLY : lsda, current_spin, isk
    USE uspp,              ONLY : nkb, vkb, okvan
    USE uspp_init,         ONLY : init_us_2
    USE rdmft_xc,          ONLY : rdmft_xc_n_channels, rdmft_xc_channel, &
                                  rdmft_xc_orb_exchange
    USE rdmft_stiefel,     ONLY : stiefel_project_tangent_gamma, &
                                  stiefel_project_tangent_k, &
                                  stiefel_inner_product_gamma, &
                                  stiefel_inner_product_k
    USE mp,                ONLY : mp_sum
    USE mp_pools,          ONLY : inter_pool_comm
    !
    COMPLEX(DP), INTENT(OUT) :: G_R_all(:,:,:)
    REAL(DP),    INTENT(OUT) :: gnorm2
    !
    COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), gradC(:,:), eta(:,:)
    COMPLEX(DP), ALLOCATABLE :: vxpsi_total(:,:,:), vxpsi_buf(:,:,:)
    REAL(DP), ALLOCATABLE :: vx_diag(:,:)
    REAL(DP) :: coef, w, dw, e_t
    INTEGER :: ik, ib, npw, it, nch
    !
    nch = rdmft_xc_n_channels()
    ALLOCATE(hpsi(npwx*npol, nbnd), gradC(npwx*npol, nbnd), eta(npwx*npol, nbnd))
    ALLOCATE(vxpsi_total(npwx*npol, nbnd, nks), vx_diag(nbnd, nks))
    !
    ! ---- Pre-compute the channel-summed XC contribution to the ----
    ! ---- ambient gradient for every k, ONCE.                   ----
    !
    ! For each XC channel ``it`` we want
    !
    !   contribution_t(:, ib, ik) = coef_t * w_t(n_ib_ik) * Vx[gamma_t]|psi_ib_ik>
    !
    ! and the gradient is the sum over channels.  The expensive part
    ! is the all-k ACE build inside ``rdmft_compute_xc_channel`` (one
    ! exxinit + one aceinit + one vexxace_* sweep, total ``O(N_k^2)``);
    ! it depends ONLY on ``it`` (via the channel weights wg_t = wk *
    ! w_t(n)), NOT on the per-k Stiefel projection logic that follows.
    ! Hoisting this loop out of the ``DO ik`` below saves a factor of
    ! :math:`N_k` of redundant aceinit work and is the dominant
    ! contribution to the (super-linear) EXX cost the user reported.
    !
    ! For ``nch = 1`` (HF / Muller / Power / GU -- the typical case)
    ! the loop runs once and ``vxpsi_total`` is just ``coef * w *
    ! vxpsi_buf``.  For ``nch > 1`` (CHF / CGA / GEO / HybOpt) we
    ! accumulate.  Memory cost is one extra (npwx*nbnd*nks) complex
    ! buffer in addition to ``vxpsi_total``; ``vxpsi_buf`` is reused
    ! across channels.
    ALLOCATE(vxpsi_buf(npwx*npol, nbnd, nks))
    vxpsi_total = (0.0_DP, 0.0_DP)
    IF (nch == 0) THEN
       CALL rdmft_xc_orb_exchange(e_t, vxpsi_total, vx_diag, vxpsi_buf)
    ELSE
    DO it = 1, nch
       CALL rdmft_compute_xc_channel(it, e_t, vx_diag, vxpsi_buf)
       ! Distinct (ib, ik) slices: thread-safe channel accumulation.
!$omp parallel do collapse(2) default(shared) private(ik, ib, coef, w, dw)
       DO ik = 1, nks
          DO ib = 1, nbnd
             CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
             vxpsi_total(:, ib, ik) = vxpsi_total(:, ib, ik) &
                                    + coef * w * vxpsi_buf(:, ib, ik)
          ENDDO
       ENDDO
!$omp end parallel do
    ENDDO
    ENDIF
    DEALLOCATE(vxpsi_buf)
    !
    gnorm2 = 0.0_DP
    DO ik = 1, nks
       npw = ngk(ik)
       current_k = ik
       IF (lsda) current_spin = isk(ik)
       IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       CALL g2_kin(ik)
       IF (okvan .OR. nkb > 0) CALL init_us_2(npw, igk_k(1, ik), xk(:, ik), vkb)
       CALL rdmft_apply_h_one_psi(npw, nbnd, evc, hpsi)
       ! Ambient gradient = wk * (n_ib * H_one |psi> + Vx_total |psi>),
       ! with Vx_total already channel-summed above.
       DO ib = 1, nbnd
          gradC(:, ib) = wk(ik) * (rdmft_n(ib, ik) * hpsi(:, ib) &
                                   + vxpsi_total(:, ib, ik))
       ENDDO
       IF (gamma_only) THEN
          CALL stiefel_project_tangent_gamma(evc, gradC, npw, nbnd, npwx*npol, gstart, eta)
          gnorm2 = gnorm2 + stiefel_inner_product_gamma(eta, eta, npw, nbnd, npwx*npol, gstart)
       ELSE
          CALL stiefel_project_tangent_k(evc, gradC, npw, nbnd, npwx*npol, eta)
          gnorm2 = gnorm2 + stiefel_inner_product_k(eta, eta, npw, nbnd, npwx*npol)
       ENDIF
       G_R_all(:, :, ik) = eta(:, :)
    ENDDO
    ! The per-k Stiefel inner products are reduced over the plane-wave
    ! band group inside stiefel_inner_product_*, but gnorm2 still
    ! accumulates only this pool's k-points.  Reduce over the k-point
    ! pools so every rank gets the global ||G_R||^2 (essential for a
    ! consistent orbital convergence test and line-search slope when
    ! npool > 1).
    CALL mp_sum(gnorm2, inter_pool_comm)
    DEALLOCATE(hpsi, gradC, eta, vxpsi_total, vx_diag)
  END SUBROUTINE rdmft_compute_riemannian_gradient
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_band_energies(h_diag_all)
    !---------------------------------------------------------------
    !! Per-band diagonal Hamiltonian \(h^{\mathrm{diag}}_{ik} =
    !! \langle\psi_{ik}|H_{\mathrm{one}}|\psi_{ik}\rangle\) for every
    !! (band, k-point) pair, returned in a single array.  Used by
    !! the occupation and orbital preconditioners.
    !
    USE wvfct, ONLY : nbnd
    USE klist, ONLY : nks, ngk
    !
    REAL(DP), INTENT(OUT) :: h_diag_all(nbnd, nks)
    REAL(DP), ALLOCATABLE :: h_diag(:)
    INTEGER :: ik, npw
    !
    ALLOCATE(h_diag(nbnd))
    DO ik = 1, nks
       npw = ngk(ik)
       CALL rdmft_compute_one_body_diag(ik, npw, h_diag)
       h_diag_all(:, ik) = h_diag(:)
    ENDDO
    DEALLOCATE(h_diag)
    !
  END SUBROUTINE rdmft_band_energies
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_apply_occ_precond(grad, h_diag_all, grad_pc)
    !---------------------------------------------------------------
    !! Diagonal preconditioner for the occupation gradient:
    !!
    !!   \(\tilde g_{ik} = g_{ik} / (|h^{\mathrm{diag}}_{ik}| + \tau)\)
    !!
    !! with \(\tau = \texttt{rdmft\_occ\_precond\_shift}\) (Ry).
    !! Rescales bands of very different energy to comparable
    !! gradient magnitudes so the projected line search converges
    !! in fewer backtracks.  No-op when
    !! \texttt{rdmft\_occ\_precond = .FALSE.}.
    !
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks
    USE rdmft_module, ONLY : rdmft_occ_precond, rdmft_occ_precond_shift
    !
    REAL(DP), INTENT(IN)  :: grad(nbnd, nks)
    REAL(DP), INTENT(IN)  :: h_diag_all(nbnd, nks)
    REAL(DP), INTENT(OUT) :: grad_pc(nbnd, nks)
    INTEGER  :: ib, ik
    REAL(DP) :: denom
    !
    IF (.NOT. rdmft_occ_precond) THEN
       grad_pc(:, :) = grad(:, :)
       RETURN
    ENDIF
    DO ik = 1, nks
       DO ib = 1, nbnd
          denom = ABS(h_diag_all(ib, ik)) + rdmft_occ_precond_shift
          grad_pc(ib, ik) = grad(ib, ik) / denom
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_apply_occ_precond
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_elk_occ_scale(elk_scale_sqrt)
    !---------------------------------------------------------------
    !! ELK-style diagonal occupation scaling (ABACUS mirror).
    !!
    !! When the fixed-orbital cache is active, uses one uniform probe
    !! at \(n_{ik}=0.5\) instead of \(\mathcal{O}(N_b N_k)\) full
    !! \texttt{rdmft\_grad\_n} calls.
    USE wvfct,            ONLY : nbnd
    USE klist,            ONLY : nks, wk
    USE rdmft_module, ONLY : rdmft_n
    !
    REAL(DP), INTENT(OUT) :: elk_scale_sqrt(nbnd, nks)
    !
    REAL(DP), ALLOCATABLE :: grad_probe(:,:), n_save(:,:)
    REAL(DP) :: mean_abs_grad, grad_floor, m_inv, abs_grad
    REAL(DP), PARAMETER :: k_inv_min = 0.1_DP, k_inv_max = 10.0_DP
    REAL(DP), PARAMETER :: n_probe = 0.5_DP
    INTEGER :: ik, ib, n_finite
    !
    ALLOCATE(grad_probe(nbnd, nks), n_save(nbnd, nks))
    n_save = rdmft_n
    mean_abs_grad = 0.0_DP
    n_finite = 0
    !
    DO ik = 1, nks
       DO ib = 1, nbnd
          rdmft_n(ib, ik) = n_probe
          CALL rdmft_grad_n(grad_probe)
          rdmft_n(ib, ik) = n_save(ib, ik)
          abs_grad = ABS(grad_probe(ib, ik))
          mean_abs_grad = mean_abs_grad + abs_grad
          n_finite = n_finite + 1
       ENDDO
    ENDDO
    !
    IF (n_finite > 0) THEN
       mean_abs_grad = mean_abs_grad / REAL(n_finite, DP)
    ELSE
       mean_abs_grad = 1.0_DP
    ENDIF
    mean_abs_grad = MAX(mean_abs_grad, 1.0e-8_DP)
    grad_floor = MAX(1.0e-8_DP, 1.0e-3_DP * mean_abs_grad)
    !
    DO ik = 1, nks
       DO ib = 1, nbnd
          abs_grad = ABS(grad_probe(ib, ik))
          m_inv = mean_abs_grad / MAX(abs_grad, grad_floor)
          m_inv = MIN(k_inv_max, MAX(k_inv_min, m_inv))
          elk_scale_sqrt(ib, ik) = SQRT(m_inv)
       ENDDO
    ENDDO
    !
    DEALLOCATE(grad_probe, n_save)
    !
  END SUBROUTINE rdmft_compute_elk_occ_scale
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_compute_tsm_probe_energies(tsm_grad_out)
    !---------------------------------------------------------------
    !! Transition-state-model (TSM) band energies
    !!   \((\partial E/\partial n_{ik})\big|_{n_{ik}=0.5}\)
    !! for every \((ib,ik)\) pair: temporarily set that natural
    !! orbital's occupation to 1/2 (Slater/Janak transition state) and
    !! evaluate the analytic occupation gradient.  Used by the
    !! ``rdmft_dos.x`` post-processor (same quantity as ELK
    !! \texttt{rdmeval}).
    !! Cost: \(\mathcal{O}(N_b N_k)\) calls to \texttt{rdmft\_grad\_n}.
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks, wk
    USE rdmft_module, ONLY : rdmft_n
    !
    REAL(DP), INTENT(OUT) :: tsm_grad_out(nbnd, nks)
    !
    REAL(DP), ALLOCATABLE :: grad_probe(:,:), n_save(:,:)
    ! QE stores n_ik in [0,1]; ELK stores occsv in [0,occmax] with
    ! occmax=2 (nspin=1) or 1 (lsda).  n_probe=0.5 is ELK occmax/2.
    REAL(DP), PARAMETER :: n_probe = 0.5_DP
    INTEGER :: ik, ib
    !
    ALLOCATE(grad_probe(nbnd, nks), n_save(nbnd, nks))
    n_save = rdmft_n
    tsm_grad_out = 0.0_DP
    !
    ! Size ``exxbuff`` once for the TSM occupation pattern (every band
    ! at \(n=1/2\)).  Without this, probing a nearly-empty band raises
    ! ``x_nbnd_occ`` above the prior bound and each probe triggers
    ! ``rdmft_reinit_exxbuff`` (very expensive).
    CALL rdmft_ensure_exxbuff_tsm_capacity()
    !
    ! The single-k probe path in ``rdmft_compute_xc_channel`` now
    ! refreshes the complete channel occupation table with a collective
    ! poolcollect before calling ``vexx``.  Keep every pool in this probe
    ! loop so the collective count/order remains identical across pools.
    ! NOTE: the (ib, ik) loop must run identically on every k-point pool
    ! -- rdmft_grad_n -> rdmft_compute_xc_channel refreshes x_occupation
    ! with a poolcollect (an inter-pool collective), so all pools must
    ! make the same number of probe calls.  Do NOT skip iterations on a
    ! per-pool condition (e.g. n_ik ~ 0): the pools split by spin would
    ! skip different bands, desynchronise the collective and abort with
    ! MPI_ERR_TRUNCATE.
    DO ik = 1, nks
       DO ib = 1, nbnd
          rdmft_n(ib, ik) = n_probe
          rdmft_tsm_probe_ib = ib
          ! Single-k single-band Fock build: pass ``ib_only = ib`` so
          ! ``vexx`` is called with ``m = 1`` instead of ``m = nbnd``.
          ! Each probe only ever reads ``grad_probe(ib, ik)`` so the
          ! other ``nbnd - 1`` band's worth of Fock work is wasted in
          ! the all-band path.  This is the dominant cost reduction
          ! for ``rdmft_dos.x`` -- on NiO 2x2x2 / nbnd = 32 the TSM
          ! probe loop drops from O(nbnd^2 * nks * cost(FFT)) to
          ! O(nbnd * nks * cost(FFT)), a roughly nbnd-fold speedup
          ! (~30x on NiO).
          CALL rdmft_grad_n(grad_probe, only_k=ik, ib_only=ib)
          rdmft_tsm_probe_ib = 0
          rdmft_n(ib, ik) = n_save(ib, ik)
          tsm_grad_out(ib, ik) = grad_probe(ib, ik)
       ENDDO
    ENDDO
    !
    rdmft_n = n_save
    DEALLOCATE(grad_probe, n_save)
    !
  END SUBROUTINE rdmft_compute_tsm_probe_energies
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_clear_band_diagnostics()
    !---------------------------------------------------------------
    !! Release cached per-band diagnostics from
    !! \texttt{rdmft\_collect\_band\_diagnostics}.
    rdmft_diag_computed = .FALSE.
    IF (ALLOCATED(rdmft_diag_h))   DEALLOCATE(rdmft_diag_h)
    IF (ALLOCATED(rdmft_diag_vx)) DEALLOCATE(rdmft_diag_vx)
    IF (ALLOCATED(rdmft_diag_tsm)) DEALLOCATE(rdmft_diag_tsm)
  END SUBROUTINE rdmft_clear_band_diagnostics
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_collect_band_diagnostics(need_tsm)
    !---------------------------------------------------------------
    !! Compute and cache per-band RDMFT diagnostics:
    !!
    !!   h_diag(ib,ik)   = one-body diagonal
    !!   vx_diag(ib,ik)  = weighted XC diagonal
    !!   tsm(ib,ik)      = optional transition-state energy (need_tsm)
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks
    USE rdmft_module, ONLY : rdmft_n
    USE rdmft_xc,     ONLY : rdmft_xc_n_channels, rdmft_xc_channel, &
                              rdmft_xc_diag_vx_weight
    !
    LOGICAL, INTENT(IN) :: need_tsm
    !
    REAL(DP), ALLOCATABLE :: h_diag(:,:), vx_ch(:,:), vx_tot(:,:), tsm(:,:)
    REAL(DP) :: coef, w, dw, e_t
    INTEGER :: ik, ib, it, nch
    !
    IF (.NOT. ALLOCATED(rdmft_n)) RETURN
    !
    ! Fast exit if the cached diagnostics already satisfy the request.
    ! If a cheap (need_tsm=.FALSE.) call computed h/vx earlier but a
    ! later call requests transition-state energies, compute ONLY the
    ! missing TSM part rather than redoing the one-body / vx diagonals.
    IF (rdmft_diag_computed) THEN
       IF (.NOT. need_tsm) RETURN
       IF (ALLOCATED(rdmft_diag_tsm)) RETURN
       ALLOCATE(tsm(nbnd, nks))
       CALL rdmft_compute_tsm_probe_energies(tsm)
       ALLOCATE(rdmft_diag_tsm(nbnd, nks))
       rdmft_diag_tsm = tsm
       DEALLOCATE(tsm)
       RETURN
    ENDIF
    !
    ALLOCATE(h_diag(nbnd, nks), vx_ch(nbnd, nks), vx_tot(nbnd, nks))
    CALL rdmft_band_energies(h_diag)
    vx_tot = 0.0_DP
    nch = rdmft_xc_n_channels()
    IF (nch == 0) THEN
       CALL rdmft_compute_xc_channel(1, e_t, vx_ch)
       vx_tot = vx_ch
    ELSE
    DO it = 1, nch
       CALL rdmft_compute_xc_channel(it, e_t, vx_ch)
       DO ik = 1, nks
          DO ib = 1, nbnd
             CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
             vx_tot(ib, ik) = vx_tot(ib, ik) + coef * vx_ch(ib, ik)
          ENDDO
       ENDDO
    ENDDO
    ENDIF
    IF (need_tsm) THEN
       ALLOCATE(tsm(nbnd, nks))
       CALL rdmft_compute_tsm_probe_energies(tsm)
    ENDIF
    !
    rdmft_diag_h   = h_diag
    rdmft_diag_vx  = vx_tot
    IF (need_tsm) THEN
       IF (ALLOCATED(rdmft_diag_tsm)) DEALLOCATE(rdmft_diag_tsm)
       ALLOCATE(rdmft_diag_tsm(nbnd, nks))
       rdmft_diag_tsm = tsm
       DEALLOCATE(tsm)
    ENDIF
    DEALLOCATE(h_diag, vx_ch, vx_tot)
    rdmft_diag_computed = .TRUE.
    !
  END SUBROUTINE rdmft_collect_band_diagnostics
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_build_occ_grad_dir(grad_n, elk_scale_sqrt, grad_dir)
    !---------------------------------------------------------------
    !! Occupation search direction from \(\partial E/\partial n_{ik}\).
    !! With \texttt{rdmft\_occ\_precond}, apply ELK diagonal scaling
    !! \(\sqrt{m^{-1}_{ik}}\) to \texttt{grad\_n}; otherwise copy as-is.
    USE wvfct, ONLY : nbnd
    USE klist, ONLY : nks
    USE rdmft_module, ONLY : rdmft_occ_precond
    !
    REAL(DP), INTENT(IN)  :: grad_n(nbnd, nks), elk_scale_sqrt(nbnd, nks)
    REAL(DP), INTENT(OUT) :: grad_dir(nbnd, nks)
    INTEGER :: ib, ik
    !
    IF (rdmft_occ_precond) THEN
       DO ik = 1, nks
          DO ib = 1, nbnd
             grad_dir(ib, ik) = elk_scale_sqrt(ib, ik) * grad_n(ib, ik)
          ENDDO
       ENDDO
    ELSE
       grad_dir(:, :) = grad_n(:, :)
    ENDIF
    !
  END SUBROUTINE rdmft_build_occ_grad_dir
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_apply_orb_precond(eta_in, eta_out, npw, nb, lda, &
                                      eps_band, shift_in)
    !---------------------------------------------------------------
    !! Level-shift / Teter-Payne-Allen rotation preconditioner for
    !! the Stiefel orbital gradient at the **current** k-point.
    !! Mirrors QE's ``g_psi`` band-by-band preconditioner: for each
    !! band ``ib`` and plane wave ``ig``, divide by
    !!
    !!   \(\mathrm{denm}_{ib,ig} = \tfrac12 (\mathbf k + \mathbf g)^2 + V_0
    !!                            - \epsilon_{ib} + \mathrm{shift}\),
    !!
    !! clamped so it stays strictly positive.  Without the
    !! preconditioner, high-G components of the gradient dominate
    !! and the Stiefel line search has to take very small steps.
    !!
    !! ``eta_in`` and ``eta_out`` may alias the same array.
    !! Caller must have called ``g2_kin(ik)`` first so
    !! :module:`wvfct::g2kin` holds the current k-point's
    !! \((k+g)^2\).
    !
    USE wvfct,            ONLY : g2kin
    USE noncollin_module, ONLY : npol
    USE control_flags,    ONLY : gamma_only
    USE rdmft_module,     ONLY : rdmft_orb_precond
    !
    INTEGER,     INTENT(IN)  :: npw, nb, lda
    COMPLEX(DP), INTENT(IN)  :: eta_in (lda, nb)
    COMPLEX(DP), INTENT(OUT) :: eta_out(lda, nb)
    REAL(DP),    INTENT(IN)  :: eps_band(nb)
    REAL(DP),    INTENT(IN)  :: shift_in
    !
    INTEGER  :: ib, ig
    REAL(DP) :: denm, denm_min
    !
    IF (.NOT. rdmft_orb_precond) THEN
       eta_out(:, :) = eta_in(:, :)
       RETURN
    ENDIF
    denm_min = MAX(1.0e-3_DP, 0.5_DP * shift_in)
    DO ib = 1, nb
       DO ig = 1, npw
          denm = 0.5_DP * g2kin(ig) + shift_in - eps_band(ib)
          IF (denm < denm_min) denm = denm_min
          eta_out(ig, ib) = eta_in(ig, ib) / denm
       ENDDO
       ! Pad above npw if lda > npw (Stiefel arrays are allocated
       ! at npwx*npol; the kinetic operator is undefined there).
       DO ig = npw + 1, lda
          eta_out(ig, ib) = (0.0_DP, 0.0_DP)
       ENDDO
    ENDDO
    !
  END SUBROUTINE rdmft_apply_orb_precond
  !
  !
END MODULE rdmft_energy
