!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_density
  !------------------------------------------------------------------
  !! RDMFT particle density and 1-RDM coherency diagnostics from the
  !! natural-orbital 1-RDM, for use by ``pp.x`` (plot_num 126--128).
  !!
  !!   \(\rho(\mathbf r)=\gamma(\mathbf r,\mathbf r)
  !!   = \sum_k w_k\sum_i n_{ik}|\psi_{ik}(\mathbf r)|^2\)
  !!
  !!   \(h(\mathbf r_0,\mathbf r)=\gamma(\mathbf r_0,\mathbf r)
  !!   -\rho(\mathbf r_0)\rho(\mathbf r)\)  (plot 127; NOT the XC hole)
  !!
  !!   \(h(\mathbf r,\mathbf r)=\rho(\mathbf r)-\rho(\mathbf r)^2\)
  !!   (plot 128; on-site 1-RDM diagnostic)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rdmft_build_particle_density, rdmft_fill_raux_rho, &
            rdmft_fill_raux_xhole_ref, rdmft_fill_raux_xhole_onsite, &
            rdmft_print_density_integrals
  !
  REAL(DP), PARAMETER :: XHOLE_ORIGIN_AUTO = -1.0_DP
  LOGICAL, SAVE :: xhole_fields_built = .FALSE.
  INTEGER, SAVE :: ref_i0 = 0, ref_j0 = 0, ref_k0 = 0
  REAL(DP), ALLOCATABLE, SAVE :: g_ref_cache(:), g_onsite_cache(:)
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_psi_at_ijk(i0, j0, k0, dfft, psic, psi_val)
    !---------------------------------------------------------------
    USE fft_types, ONLY : fft_type_descriptor, fft_index_to_3d
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: i0, j0, k0
    COMPLEX(DP), INTENT(IN) :: psic(:)
    COMPLEX(DP), INTENT(OUT) :: psi_val
    INTEGER :: ir, i, j, k
    LOGICAL :: off
    psi_val = (0.0_DP, 0.0_DP)
    DO ir = 1, dfft%nnr
       CALL fft_index_to_3d(ir, dfft, i, j, k, off)
       IF (.NOT. off .AND. i == i0 .AND. j == j0 .AND. k == k0) THEN
          psi_val = psic(ir)
          RETURN
       ENDIF
    ENDDO
  END SUBROUTINE rdmft_psi_at_ijk
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_spinor_at_ijk(i0, j0, k0, dfft, psic_up, psic_dn, psi_up, psi_dn)
    !---------------------------------------------------------------
    USE fft_types, ONLY : fft_type_descriptor
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: i0, j0, k0
    COMPLEX(DP), INTENT(IN) :: psic_up(:), psic_dn(:)
    COMPLEX(DP), INTENT(OUT) :: psi_up, psi_dn
    CALL rdmft_psi_at_ijk(i0, j0, k0, dfft, psic_up, psi_up)
    CALL rdmft_psi_at_ijk(i0, j0, k0, dfft, psic_dn, psi_dn)
  END SUBROUTINE rdmft_spinor_at_ijk
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_add_gamma_contribution(i0, j0, k0, dfft, weight, psic, gamma_r)
    !---------------------------------------------------------------
    USE mp, ONLY : mp_sum
    USE mp_images, ONLY : intra_image_comm
    USE fft_types, ONLY : fft_type_descriptor
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: i0, j0, k0
    REAL(DP), INTENT(IN) :: weight
    COMPLEX(DP), INTENT(IN) :: psic(:)
    REAL(DP), INTENT(INOUT) :: gamma_r(:)
    COMPLEX(DP) :: psi0, prod
    INTEGER :: ir
    CALL rdmft_psi_at_ijk(i0, j0, k0, dfft, psic, psi0)
    CALL mp_sum(psi0, intra_image_comm)
    IF (ABS(weight) < 1.0e-30_DP .OR. ABS(psi0) < 1.0e-30_DP) RETURN
    DO ir = 1, dfft%nnr
       prod = psi0 * CONJG(psic(ir))
       gamma_r(ir) = gamma_r(ir) + weight * REAL(prod, KIND=DP)
    ENDDO
  END SUBROUTINE rdmft_add_gamma_contribution
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_add_gamma_contribution_nc(i0, j0, k0, dfft, weight, &
                                            psic_up, psic_dn, gamma_r)
    !---------------------------------------------------------------
    USE mp, ONLY : mp_sum
    USE mp_images, ONLY : intra_image_comm
    USE fft_types, ONLY : fft_type_descriptor
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: i0, j0, k0
    REAL(DP), INTENT(IN) :: weight
    COMPLEX(DP), INTENT(IN) :: psic_up(:), psic_dn(:)
    REAL(DP), INTENT(INOUT) :: gamma_r(:)
    COMPLEX(DP) :: psi0_up, psi0_dn, prod
    INTEGER :: ir
    CALL rdmft_spinor_at_ijk(i0, j0, k0, dfft, psic_up, psic_dn, psi0_up, psi0_dn)
    CALL mp_sum(psi0_up, intra_image_comm)
    CALL mp_sum(psi0_dn, intra_image_comm)
    IF (ABS(weight) < 1.0e-30_DP) RETURN
    DO ir = 1, dfft%nnr
       prod = psi0_up * CONJG(psic_up(ir)) + psi0_dn * CONJG(psic_dn(ir))
       gamma_r(ir) = gamma_r(ir) + weight * REAL(prod, KIND=DP)
    ENDDO
  END SUBROUTINE rdmft_add_gamma_contribution_nc
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_gather_real_field(dfft, f_loc, f_g)
    !---------------------------------------------------------------
    USE scatter_mod, ONLY : gather_grid
    USE fft_types,   ONLY : fft_type_descriptor
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    REAL(DP), INTENT(IN)  :: f_loc(:)
    REAL(DP), INTENT(OUT) :: f_g(:)
#if defined(__MPI)
    CALL gather_grid(dfft, f_loc, f_g)
#else
    IF (SIZE(f_g) >= SIZE(f_loc)) f_g(1:SIZE(f_loc)) = f_loc
#endif
  END SUBROUTINE rdmft_gather_real_field
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_integrate_field(dfft, f_g)
    !---------------------------------------------------------------
    USE cell_base, ONLY : omega
    USE fft_types, ONLY : fft_type_descriptor
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    REAL(DP), INTENT(IN) :: f_g(:)
    REAL(DP) :: vol
    vol = REAL(dfft%nr1 * dfft%nr2 * dfft%nr3, KIND=DP)
    IF (vol < 1.0e-30_DP) THEN
       rdmft_integrate_field = 0.0_DP
    ELSE
       rdmft_integrate_field = SUM(f_g(1:dfft%nr1x * dfft%nr2x * dfft%nr3x)) &
            * omega / vol
    ENDIF
  END FUNCTION rdmft_integrate_field
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_select_reference_ijk(dfft, rho_loc, i0, j0, k0)
    !---------------------------------------------------------------
    USE io_global,  ONLY : ionode
    USE mp,         ONLY : mp_bcast
    USE mp_images,  ONLY : intra_image_comm, ionode_id
    USE fft_types,  ONLY : fft_type_descriptor
    USE rdmft_module, ONLY : rdmft_xhole_origin
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    REAL(DP), INTENT(IN) :: rho_loc(:)
    INTEGER, INTENT(OUT) :: i0, j0, k0
    REAL(DP), ALLOCATABLE :: rho_g(:)
    INTEGER :: idx, idx0, rem, nphys
    LOGICAL :: use_auto
    REAL(DP) :: frac(3)
    use_auto = (rdmft_xhole_origin(1) <= XHOLE_ORIGIN_AUTO + 1.0e-12_DP .AND. &
                rdmft_xhole_origin(2) <= XHOLE_ORIGIN_AUTO + 1.0e-12_DP .AND. &
                rdmft_xhole_origin(3) <= XHOLE_ORIGIN_AUTO + 1.0e-12_DP)
    i0 = 0
    j0 = 0
    k0 = 0
    IF (use_auto) THEN
       IF (ionode) THEN
          nphys = dfft%nr1 * dfft%nr2 * dfft%nr3
          ALLOCATE(rho_g(dfft%nr1x * dfft%nr2x * dfft%nr3x))
          CALL rdmft_gather_real_field(dfft, rho_loc, rho_g)
          idx = MAXLOC(rho_g(1:nphys), 1)
          idx0 = idx - 1
          i0 = idx0 / (dfft%nr2x * dfft%nr3x)
          rem = MOD(idx0, dfft%nr2x * dfft%nr3x)
          j0 = rem / dfft%nr3x
          k0 = MOD(rem, dfft%nr3x)
          DEALLOCATE(rho_g)
       ENDIF
    ELSE
       IF (ionode) THEN
          frac = rdmft_xhole_origin
          frac = frac - NINT(frac - 0.5_DP)
          i0 = MOD(NINT(frac(1) * REAL(dfft%nr1, KIND=DP)), dfft%nr1)
          j0 = MOD(NINT(frac(2) * REAL(dfft%nr2, KIND=DP)), dfft%nr2)
          k0 = MOD(NINT(frac(3) * REAL(dfft%nr3, KIND=DP)), dfft%nr3)
          IF (i0 < 0) i0 = i0 + dfft%nr1
          IF (j0 < 0) j0 = j0 + dfft%nr2
          IF (k0 < 0) k0 = k0 + dfft%nr3
       ENDIF
    ENDIF
    CALL mp_bcast(i0, ionode_id, intra_image_comm)
    CALL mp_bcast(j0, ionode_id, intra_image_comm)
    CALL mp_bcast(k0, ionode_id, intra_image_comm)
  END SUBROUTINE rdmft_select_reference_ijk
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_build_particle_density()
    !---------------------------------------------------------------
    USE rdmft_energy, ONLY : rdmft_set_wg_from_n
    USE wvfct,        ONLY : wg, nbnd
    USE klist,        ONLY : nks
    REAL(DP), ALLOCATABLE :: wg_save(:,:)
    ALLOCATE(wg_save(nbnd, nks))
    wg_save = wg
    CALL rdmft_set_wg_from_n(1)
    CALL sum_band()
    wg = wg_save
    DEALLOCATE(wg_save)
  END SUBROUTINE rdmft_build_particle_density
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_total_rho_on_dfft(rho_tot)
    !---------------------------------------------------------------
    USE fft_base,         ONLY : dfftp
    USE scf,              ONLY : rho
    USE lsda_mod,         ONLY : nspin, lsda
    USE noncollin_module, ONLY : noncolin
    REAL(DP), INTENT(OUT) :: rho_tot(:)
    INTEGER :: is
    rho_tot = 0.0_DP
    IF (noncolin) THEN
       rho_tot(1:dfftp%nnr) = rho%of_r(1:dfftp%nnr, 1)
    ELSE IF (lsda .AND. nspin == 2) THEN
       DO is = 1, nspin
          rho_tot(1:dfftp%nnr) = rho_tot(1:dfftp%nnr) + rho%of_r(1:dfftp%nnr, is)
       ENDDO
    ELSE
       rho_tot(1:dfftp%nnr) = rho%of_r(1:dfftp%nnr, 1)
    ENDIF
  END SUBROUTINE rdmft_total_rho_on_dfft
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_symmetrize_smooth_field(gamma_s, gamma_p)
    !---------------------------------------------------------------
    USE fft_base, ONLY : dffts, dfftp
    USE gvect,    ONLY : ngm
    USE fft_rho,  ONLY : rho_r2g, rho_g2r
    USE symme,    ONLY : sym_rho
    REAL(DP), INTENT(IN)  :: gamma_s(:)
    REAL(DP), INTENT(OUT) :: gamma_p(:)
    COMPLEX(DP), ALLOCATABLE :: rhog(:,:)
    ALLOCATE(rhog(ngm, 1))
    CALL rho_r2g(dffts, gamma_s, rhog)
    CALL sym_rho(1, rhog)
    CALL rho_g2r(dfftp, rhog(:, 1), gamma_p)
    DEALLOCATE(rhog)
  END SUBROUTINE rdmft_symmetrize_smooth_field
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_build_xhole_fields(i0, j0, k0, g_ref, g_onsite)
    !---------------------------------------------------------------
    USE io_global,         ONLY : stdout, ionode
    USE wvfct,             ONLY : nbnd, npwx, current_k
    USE klist,             ONLY : nks, ngk, igk_k, xk, wk
    USE wavefunctions,     ONLY : evc, psic, psic_nc
    USE lsda_mod,          ONLY : lsda, isk, current_spin
    USE noncollin_module,  ONLY : noncolin, npol
    USE control_flags,     ONLY : gamma_only
    USE io_files,          ONLY : nwordwfc, iunwfc
    USE buffers,           ONLY : get_buffer
    USE uspp,              ONLY : okvan, vkb, nkb
    USE uspp_init,         ONLY : init_us_2
    USE becmod,            ONLY : becp
    USE fft_base,          ONLY : dffts, dfftp
    USE fft_wave,          ONLY : wave_g2r
    USE mp,                ONLY : mp_sum
    USE mp_pools,          ONLY : inter_pool_comm
    USE mp_bands,          ONLY : inter_bgrp_comm
    USE mp_images,         ONLY : intra_image_comm
    USE rdmft_module,      ONLY : rdmft_n
    USE fft_types,         ONLY : fft_index_to_3d
    INTEGER, INTENT(IN)  :: i0, j0, k0
    REAL(DP), INTENT(OUT) :: g_ref(:), g_onsite(:)
    REAL(DP), ALLOCATABLE :: gamma_s(:), gamma_p(:), rho_tot(:)
    REAL(DP) :: weight, rho0
    COMPLEX(DP), ALLOCATABLE :: psicd(:)
    INTEGER :: ik, ib, npw, ir
    LOGICAL :: off
    INTEGER :: i, j, k
    ALLOCATE(gamma_s(dffts%nnr), gamma_p(dfftp%nnr), rho_tot(dfftp%nnr))
    gamma_s = 0.0_DP
    DO ik = 1, nks
       npw = ngk(ik)
       current_k = ik
       IF (lsda) current_spin = isk(ik)
       IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       IF (nkb > 0) CALL init_us_2(npw, igk_k(1, ik), xk(1, ik), vkb, .TRUE.)
       DO ib = 1, nbnd
          weight = wk(ik) * rdmft_n(ib, ik)
          IF (ABS(weight) < 1.0e-30_DP) CYCLE
          IF (noncolin) THEN
             CALL wave_g2r(evc(1:npw, ib:ib), psic_nc(:, 1), dffts, igk=igk_k(:, ik))
             CALL wave_g2r(evc(npwx + 1:npwx + npw, ib:ib), psic_nc(:, 2), &
                  dffts, igk=igk_k(:, ik))
             CALL rdmft_add_gamma_contribution_nc(i0, j0, k0, dffts, weight, &
                  psic_nc(:, 1), psic_nc(:, 2), gamma_s)
          ELSE IF (gamma_only) THEN
             CALL wave_g2r(evc(1:npw, ib:ib), psic, dffts)
             CALL rdmft_add_gamma_contribution(i0, j0, k0, dffts, weight, psic, gamma_s)
          ELSE
             ALLOCATE(psicd(dffts%nnr))
             CALL wave_g2r(evc(1:npw, ib:ib), psicd, dffts, igk=igk_k(:, ik))
             CALL rdmft_add_gamma_contribution(i0, j0, k0, dffts, weight, psicd, gamma_s)
             DEALLOCATE(psicd)
          ENDIF
       ENDDO
    ENDDO
    CALL mp_sum(gamma_s, inter_pool_comm)
    CALL mp_sum(gamma_s, inter_bgrp_comm)
    CALL rdmft_symmetrize_smooth_field(gamma_s, gamma_p)
    CALL rdmft_total_rho_on_dfft(rho_tot)
    rho0 = 0.0_DP
    DO ir = 1, dfftp%nnr
       CALL fft_index_to_3d(ir, dfftp, i, j, k, off)
       IF (.NOT. off .AND. i == i0 .AND. j == j0 .AND. k == k0) rho0 = rho_tot(ir)
    ENDDO
    CALL mp_sum(rho0, intra_image_comm)
    DO ir = 1, dfftp%nnr
       g_ref(ir) = gamma_p(ir) - rho0 * rho_tot(ir)
       g_onsite(ir) = rho_tot(ir) - rho_tot(ir)**2
    ENDDO
    IF (ionode) THEN
       WRITE(stdout, '(5X,A,3(I0,1X))') 'RDMFT coherency slice reference grid (i,j,k) = ', &
            i0, j0, k0
       WRITE(stdout, '(5X,A,F14.8)') '   rho(r0) = ', rho0
    ENDIF
    DEALLOCATE(gamma_s, gamma_p, rho_tot)
  END SUBROUTINE rdmft_build_xhole_fields
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_ensure_xhole_cache()
    !---------------------------------------------------------------
    USE fft_base, ONLY : dfftp
    REAL(DP), ALLOCATABLE :: rho_loc(:)
    IF (xhole_fields_built) RETURN
    IF (.NOT. ALLOCATED(g_ref_cache)) ALLOCATE(g_ref_cache(dfftp%nnr))
    IF (.NOT. ALLOCATED(g_onsite_cache)) ALLOCATE(g_onsite_cache(dfftp%nnr))
    ALLOCATE(rho_loc(dfftp%nnr))
    CALL rdmft_total_rho_on_dfft(rho_loc)
    CALL rdmft_select_reference_ijk(dfftp, rho_loc, ref_i0, ref_j0, ref_k0)
    CALL rdmft_build_xhole_fields(ref_i0, ref_j0, ref_k0, g_ref_cache, g_onsite_cache)
    DEALLOCATE(rho_loc)
    xhole_fields_built = .TRUE.
  END SUBROUTINE rdmft_ensure_xhole_cache
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_fill_raux_rho(raux)
    !---------------------------------------------------------------
    USE fft_base, ONLY : dfftp
    REAL(DP), INTENT(OUT) :: raux(:)
    REAL(DP), ALLOCATABLE :: rho_loc(:)
    ALLOCATE(rho_loc(dfftp%nnr))
    CALL rdmft_total_rho_on_dfft(rho_loc)
    raux(1:dfftp%nnr) = rho_loc(1:dfftp%nnr)
    DEALLOCATE(rho_loc)
  END SUBROUTINE rdmft_fill_raux_rho
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_fill_raux_xhole_ref(raux, title)
    !---------------------------------------------------------------
    USE fft_base, ONLY : dfftp
    CHARACTER(LEN=*), INTENT(OUT) :: title
    REAL(DP), INTENT(OUT) :: raux(:)
    CALL rdmft_ensure_xhole_cache()
    raux(1:dfftp%nnr) = g_ref_cache(1:dfftp%nnr)
     WRITE(title, '(A,3(I0,1X),A)') 'RDMFT 1-RDM coherency h(r0,r), ref grid (i,j,k)= ', &
         ref_i0, ref_j0, ref_k0, ' '
  END SUBROUTINE rdmft_fill_raux_xhole_ref
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_fill_raux_xhole_onsite(raux)
    !---------------------------------------------------------------
    USE fft_base, ONLY : dfftp
    REAL(DP), INTENT(OUT) :: raux(:)
    CALL rdmft_ensure_xhole_cache()
    raux(1:dfftp%nnr) = g_onsite_cache(1:dfftp%nnr)
  END SUBROUTINE rdmft_fill_raux_xhole_onsite
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_print_density_integrals()
    !---------------------------------------------------------------
    USE io_global, ONLY : stdout, ionode
    USE fft_base,  ONLY : dfftp
    USE klist,     ONLY : nelec
    REAL(DP), ALLOCATABLE :: rho_loc(:), rho_g(:)
    REAL(DP) :: charge_rho
    IF (.NOT. ionode) RETURN
    ALLOCATE(rho_loc(dfftp%nnr), rho_g(dfftp%nr1x * dfftp%nr2x * dfftp%nr3x))
    CALL rdmft_total_rho_on_dfft(rho_loc)
    CALL rdmft_gather_real_field(dfftp, rho_loc, rho_g)
    charge_rho = rdmft_integrate_field(dfftp, rho_g)
    WRITE(stdout, '(5X,A,F12.6,A,F12.6)') 'Integrated rho = ', charge_rho, &
         '   Ne target = ', nelec
    DEALLOCATE(rho_loc, rho_g)
  END SUBROUTINE rdmft_print_density_integrals
  !
END MODULE rdmft_density
