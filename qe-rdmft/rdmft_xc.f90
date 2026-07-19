!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_xc
  !------------------------------------------------------------------
  !! RDMFT exchange-correlation functionals.
  !!
  !! Muller / Power / GEO / HybOpt small-$n$ handling mirrors ABACUS
  !! \texttt{rdmft\_xc\_functional.h} (\texttt{pow\_reg} / \texttt{dpow\_reg}):
  !!
  !!   rdmft_g(n)  = n^\alpha  for n \ge \texttt{rdmft\_reg\_eps}
  !!         = \varepsilon^\alpha + \alpha\varepsilon^{\alpha-1}(n-\varepsilon)
  !!           for n < \varepsilon
  !!   g'(n) = \alpha \max(n,\varepsilon)^{\alpha-1}
  !!
  !! so $rdmft_g(0)=(1-\alpha)\varepsilon^\alpha$, $rdmft_g(n)\to n^\alpha$ for
  !! $n\gg\varepsilon$, and $|g'(0)|=\alpha\,\varepsilon^{\alpha-1}$ is
  !! bounded when $\alpha<1$.
  !!
  !! For a separable functional ``f(n_i, n_j) = rdmft_g(n_i) rdmft_g(n_j)`` we
  !! expose the coupling function ``rdmft_g(n)`` and its derivative.  For the
  !! non-separable kernels (GU, CHF, CGA, GEO, HybOpt) we expose the
  !! channel decomposition that PWscf needs to reuse the existing
  !! EXX-via-ACE machinery: each channel is a (coef_t, w_t(n), w_t'(n))
  !! triplet such that
  !!
  !!   E_xc[gamma] = sum_t coef_t * E_x^HF[gamma_t]
  !!
  !! where \texttt{gamma\_t} is the modified density matrix built with
  !! per-band weight \texttt{w\_t(n\_ik)}.  ABACUS' production code
  !! does exactly this (see rdmft\_derivation.md \S2.4).
  !!
  !! **Adding a new functional** (energy layer unchanged):
  !!
  !! 1. Register the name in :func:`rdmft_xc_set_type` (+ default ``rdmft_power_alpha``).
  !! 2. Choose eval kind: ``XC_EVAL_CHANNELS`` or ``XC_EVAL_PAIR``.
  !! 3. Implement kernel(s): extend :func:`rdmft_xc_channel` / :func:`rdmft_xc_n_channels`
  !!    **or** ``pair_kernel_*``.
  !! 4. If an on-site extra is needed (like GU): extend :func:`rdmft_xc_has_extra_diag`
  !!    and the ``xc_extra_*`` hooks.
  !
  USE kinds, ONLY : DP
  USE rdmft_module, ONLY : rdmft_xc_id, rdmft_power_alpha, rdmft_reg_eps, rdmft_n, &
                        rdmft_exxbuff_stale, rdmft_exxbuff_nbnd_occ, &
                        rdmft_xi_cache_valid, rdmft_xi_cache_channel, &
                        rdmft_xi_cache_n
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER, PRIVATE :: XC_EVAL_CHANNELS = 1
  INTEGER, PARAMETER, PRIVATE :: XC_EVAL_PAIR     = 2
  !
  ! XC type integer codes (kept in sync with rdmft_xc_string_to_id).
  INTEGER, PARAMETER :: RDMFT_XC_HF      = 1
  INTEGER, PARAMETER :: RDMFT_XC_MULLER  = 2
  INTEGER, PARAMETER :: RDMFT_XC_POWER   = 3
  INTEGER, PARAMETER :: RDMFT_XC_GU      = 4
  INTEGER, PARAMETER :: RDMFT_XC_CHF     = 5
  INTEGER, PARAMETER :: RDMFT_XC_CGA     = 6
  INTEGER, PARAMETER :: RDMFT_XC_GEO     = 7
  INTEGER, PARAMETER :: RDMFT_XC_HYBOPT  = 8
  INTEGER, PARAMETER :: RDMFT_XC_BOW      = 9
  INTEGER, PARAMETER :: RDMFT_XC_BOWMOD   = 10
  !
  REAL(DP), PARAMETER :: HYBOPT_POWER_W   = 0.938328_DP
  REAL(DP), PARAMETER :: HYBOPT_POWER_EXP = 0.541076_DP
  REAL(DP), PARAMETER :: BOW_DEFAULT_ALPHA = 0.61_DP
  REAL(DP), PARAMETER :: BOWMOD_DEFAULT_BETA = 0.375_DP
  REAL(DP), PARAMETER :: BOWMOD_COEF_HF = 1.0_DP
  REAL(DP), PARAMETER :: BOWMOD_COEF_ALPHA = -0.14648044160527188_DP
  REAL(DP), PARAMETER :: BOWMOD_COEF_MIX = 0.666761947146244_DP
  !
  PUBLIC :: rdmft_xc_set_type, rdmft_xc_alpha, rdmft_g, rdmft_dg, rdmft_xc_n_channels
  PUBLIC :: rdmft_xc_eval_kind, rdmft_xc_is_separable
  PUBLIC :: rdmft_binary_entropy_f, rdmft_binary_entropy_dfdn
  PUBLIC :: rdmft_xc_set_exx_wg, rdmft_xc_refresh_x_occupation, rdmft_xc_reinit_exxbuff
  PUBLIC :: rdmft_xc_exxinit_once, rdmft_xc_ensure_exxbuff_tsm_capacity
  PUBLIC :: rdmft_xc_set_wg_for_channel
  PUBLIC :: rdmft_xc_ensure_min_exx_occ_support, rdmft_xc_tsm_trim_wg_for_ace
  PUBLIC :: rdmft_xc_channel_energy_term
  PUBLIC :: rdmft_xc_total_energy, rdmft_xc_add_occ_gradient, rdmft_xc_add_occ_gradient_from_diag
  PUBLIC :: rdmft_xc_compute_exchange, rdmft_xc_occ_grad_add, rdmft_xc_occ_grad_add_from_diag
  PUBLIC :: rdmft_xc_extra_energy, rdmft_xc_extra_grad_add
  PUBLIC :: rdmft_xc_orb_exchange, rdmft_xc_probe_exchange_add, rdmft_xc_diag_vx_weight
  PUBLIC :: RDMFT_XC_HF, RDMFT_XC_MULLER, RDMFT_XC_POWER, RDMFT_XC_GU, &
            RDMFT_XC_CHF, RDMFT_XC_CGA, RDMFT_XC_GEO, RDMFT_XC_HYBOPT, &
            RDMFT_XC_BOW, RDMFT_XC_BOWMOD
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_set_type(name)
    !---------------------------------------------------------------
    !! Translate the ``functional`` string into the integer
    !! code \texttt{rdmft\_xc\_id} and force the appropriate value of
    !! \texttt{rdmft\_power\_alpha} when the functional fixes it
    !! (matches the ABACUS XCFunctional constructor).
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=16) :: lname
    INTEGER :: i
    !
    lname = ADJUSTL(name)
    DO i = 1, LEN_TRIM(lname)
       IF (lname(i:i) >= 'A' .AND. lname(i:i) <= 'Z') &
            lname(i:i) = ACHAR(IACHAR(lname(i:i)) + 32)
    ENDDO
    !
    SELECT CASE (TRIM(lname))
    CASE ('hf')
       rdmft_xc_id = RDMFT_XC_HF
       rdmft_power_alpha = 1.0_DP
    CASE ('muller')
       rdmft_xc_id = RDMFT_XC_MULLER
       rdmft_power_alpha = 0.5_DP
    CASE ('power')
       rdmft_xc_id = RDMFT_XC_POWER
       ! keep the user-provided alpha (or default 0.656)
    CASE ('gu')
       rdmft_xc_id = RDMFT_XC_GU
       rdmft_power_alpha = 0.5_DP
    CASE ('chf')
       rdmft_xc_id = RDMFT_XC_CHF
       rdmft_power_alpha = 1.0_DP
    CASE ('cga')
       rdmft_xc_id = RDMFT_XC_CGA
       rdmft_power_alpha = 1.0_DP
    CASE ('geo')
       rdmft_xc_id = RDMFT_XC_GEO
       rdmft_power_alpha = 0.75_DP
    CASE ('hybopt')
       rdmft_xc_id = RDMFT_XC_HYBOPT
       rdmft_power_alpha = HYBOPT_POWER_EXP
    CASE ('bow')
       rdmft_xc_id = RDMFT_XC_BOW
       rdmft_power_alpha = BOW_DEFAULT_ALPHA
    CASE ('bowmod')
       rdmft_xc_id = RDMFT_XC_BOWMOD
       rdmft_power_alpha = BOW_DEFAULT_ALPHA
    CASE DEFAULT
       CALL errore('rdmft_xc_set_type', &
            'Unknown functional='//TRIM(lname)//'. Allowed: hf, muller, power, gu, chf, cga, geo, hybopt, bow, bowmod.', 1)
    END SELECT
    !
  END SUBROUTINE rdmft_xc_set_type
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_xc_alpha()
    !---------------------------------------------------------------
    rdmft_xc_alpha = rdmft_power_alpha
  END FUNCTION rdmft_xc_alpha
  !
  !-----------------------------------------------------------------
  LOGICAL FUNCTION rdmft_xc_is_separable()
    !---------------------------------------------------------------
    !! ``f(n_i, n_j) = rdmft_g(n_i) rdmft_g(n_j)`` for single-channel forms.
    rdmft_xc_is_separable = (rdmft_xc_id == RDMFT_XC_HF      .OR. &
                              rdmft_xc_id == RDMFT_XC_MULLER  .OR. &
                              rdmft_xc_id == RDMFT_XC_POWER   .OR. &
                              rdmft_xc_id == RDMFT_XC_BOWMOD)
  END FUNCTION rdmft_xc_is_separable
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_binary_entropy_f(n)
    !---------------------------------------------------------------
    !! \(f_{\mathrm{bin}}(n)=n\ln n+(1-n)\ln(1-n)\), clamped for HF entropy.
    REAL(DP), INTENT(IN) :: n
    REAL(DP), PARAMETER :: occ_entropy_n_eps = 1.0e-12_DP
    REAL(DP) :: nc
    nc = MAX(occ_entropy_n_eps, MIN(1.0_DP - occ_entropy_n_eps, n))
    rdmft_binary_entropy_f = nc * LOG(nc) + (1.0_DP - nc) * LOG(1.0_DP - nc)
  END FUNCTION rdmft_binary_entropy_f
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_binary_entropy_dfdn(n)
    !---------------------------------------------------------------
    REAL(DP), INTENT(IN) :: n
    REAL(DP), PARAMETER :: occ_entropy_n_eps = 1.0e-12_DP
    REAL(DP) :: nc
    nc = MAX(occ_entropy_n_eps, MIN(1.0_DP - occ_entropy_n_eps, n))
    rdmft_binary_entropy_dfdn = LOG(nc / (1.0_DP - nc))
  END FUNCTION rdmft_binary_entropy_dfdn
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_gu_diag_factor(n)
    !---------------------------------------------------------------
    !! GU on-site diagonal factor \(n^2-n\) (ABACUS \texttt{gu\_diag\_factor}).
    REAL(DP), INTENT(IN) :: n
    REAL(DP) :: nc
    nc = MAX(0.0_DP, MIN(1.0_DP, n))
    rdmft_gu_diag_factor = nc * nc - nc
  END FUNCTION rdmft_gu_diag_factor
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_gu_diag_factor_deriv(n)
    !---------------------------------------------------------------
    REAL(DP), INTENT(IN) :: n
    REAL(DP) :: nc
    nc = MAX(0.0_DP, MIN(1.0_DP, n))
    rdmft_gu_diag_factor_deriv = 2.0_DP * nc - 1.0_DP
  END FUNCTION rdmft_gu_diag_factor_deriv
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION bow_one_minus_u(u, alpha)
    !---------------------------------------------------------------
    !! Regularised \((1-u)^{1/\alpha}\) on \(u\in[0,1]\); value and
    !! derivative vanish at \(u=1\) (full occupation).
    REAL(DP), INTENT(IN) :: u, alpha
    REAL(DP) :: uc, eps, om, g_eps, dg_eps
    !
    uc = MAX(0.0_DP, MIN(1.0_DP, u))
    eps = rdmft_reg_eps
    IF (uc >= 1.0_DP - eps) THEN
       bow_one_minus_u = 0.0_DP
       RETURN
    ENDIF
    om = 1.0_DP - uc
    IF (alpha >= 1.0_DP - 1.0e-12_DP) THEN
       bow_one_minus_u = om
    ELSEIF (om >= eps) THEN
       bow_one_minus_u = om**(1.0_DP / alpha)
    ELSE
       g_eps  = eps**(1.0_DP / alpha)
       dg_eps = (1.0_DP / alpha) * eps**(1.0_DP / alpha - 1.0_DP)
       bow_one_minus_u = g_eps + dg_eps * (om - eps)
    ENDIF
  END FUNCTION bow_one_minus_u
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION bow_done_minus_u(u, alpha)
    !---------------------------------------------------------------
    !! \(d/du\,(1-u)^{1/\alpha}\) with the same clamp at \(u=1\).
    REAL(DP), INTENT(IN) :: u, alpha
    REAL(DP) :: uc, eps, om, omin
    !
    uc = MAX(0.0_DP, MIN(1.0_DP, u))
    eps = rdmft_reg_eps
    IF (uc >= 1.0_DP - eps) THEN
       bow_done_minus_u = 0.0_DP
       RETURN
    ENDIF
    om = 1.0_DP - uc
    IF (alpha >= 1.0_DP - 1.0e-12_DP) THEN
       bow_done_minus_u = -1.0_DP
       RETURN
    ENDIF
    IF (om >= eps) THEN
       bow_done_minus_u = -(1.0_DP / alpha) * om**(1.0_DP / alpha - 1.0_DP)
    ELSE
       omin = MAX(om, eps)
       bow_done_minus_u = -(1.0_DP / alpha) * omin**(1.0_DP / alpha - 1.0_DP)
    ENDIF
  END FUNCTION bow_done_minus_u
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION bow_df_du(u, alpha)
    !---------------------------------------------------------------
    !! \(df^{\mathrm{BOW}}/du\) for \(f=(u^\alpha-\alpha u+\alpha)-\alpha(1-u)^{1/\alpha}\).
    REAL(DP), INTENT(IN) :: u, alpha
    REAL(DP) :: uc
    !
    uc = MAX(0.0_DP, MIN(1.0_DP, u))
    bow_df_du = pow_dg(uc, alpha) - alpha - alpha * bow_done_minus_u(uc, alpha)
  END FUNCTION bow_df_du
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION bow_f(ni, nj)
    !---------------------------------------------------------------
    !! Baldsiefen BOW kernel \(f^{\mathrm{BOW}}(n_i,n_j;\alpha)\).
    REAL(DP), INTENT(IN) :: ni, nj
    REAL(DP) :: nic, njc, u, alpha
    !
    nic = MAX(0.0_DP, MIN(1.0_DP, ni))
    njc = MAX(0.0_DP, MIN(1.0_DP, nj))
    u = nic * njc
    alpha = rdmft_power_alpha
    bow_f = pow_g(u, alpha) - alpha * u + alpha &
         - alpha * bow_one_minus_u(u, alpha)
  END FUNCTION bow_f
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION bow_df_dni(ni, nj)
    !---------------------------------------------------------------
    !! \(\partial f^{\mathrm{BOW}}/\partial n_i = (df/du)\, n_j\).
    REAL(DP), INTENT(IN) :: ni, nj
    REAL(DP) :: nic, njc, u
    !
    nic = MAX(0.0_DP, MIN(1.0_DP, ni))
    njc = MAX(0.0_DP, MIN(1.0_DP, nj))
    u = nic * njc
    bow_df_dni = bow_df_du(u, rdmft_power_alpha) * njc
  END FUNCTION bow_df_dni
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION bowmod_mix_w(n, beta)
    !---------------------------------------------------------------
    !! \(q(n)=[n(1-n)]^\beta\), regularised through \texttt{pow_g}.
    REAL(DP), INTENT(IN) :: n, beta
    REAL(DP) :: nc, s
    !
    nc = MAX(0.0_DP, MIN(1.0_DP, n))
    s = MAX(0.0_DP, nc * (1.0_DP - nc))
    bowmod_mix_w = pow_g(s, beta)
  END FUNCTION bowmod_mix_w
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION bowmod_mix_dw(n, beta)
    !---------------------------------------------------------------
    !! \(dq/dn\) for \(q(n)=[n(1-n)]^\beta\), regularised through
    !! \texttt{pow_dg}.
    REAL(DP), INTENT(IN) :: n, beta
    REAL(DP) :: nc, s
    !
    nc = MAX(0.0_DP, MIN(1.0_DP, n))
    s = MAX(0.0_DP, nc * (1.0_DP - nc))
    bowmod_mix_dw = pow_dg(s, beta) * (1.0_DP - 2.0_DP * nc)
  END FUNCTION bowmod_mix_dw
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION pow_g(n, alpha)
    !---------------------------------------------------------------
    !! Piecewise power coupling on $[0,1]$ (ABACUS \texttt{pow\_reg}).
    REAL(DP), INTENT(IN) :: n, alpha
    REAL(DP) :: nc, eps, g_eps, dg_eps
    !
    nc = MAX(0.0_DP, MIN(1.0_DP, n))
    IF (alpha >= 1.0_DP - 1.0e-12_DP) THEN
       pow_g = nc
    ELSE
       eps = rdmft_reg_eps
       IF (nc >= eps) THEN
          pow_g = nc**alpha
       ELSE
          g_eps  = eps**alpha
          dg_eps = alpha * eps**(alpha - 1.0_DP)
          pow_g  = g_eps + dg_eps * (nc - eps)
       ENDIF
    ENDIF
  END FUNCTION pow_g
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION pow_dg(n, alpha)
    !---------------------------------------------------------------
    !! Derivative of :func:`pow_g` w.r.t.\ $n$ (ABACUS \texttt{dpow\_reg}).
    REAL(DP), INTENT(IN) :: n, alpha
    REAL(DP) :: nc, eps, nmin
    !
    nc = MAX(0.0_DP, MIN(1.0_DP, n))
    IF (alpha >= 1.0_DP - 1.0e-12_DP) THEN
       pow_dg = 1.0_DP
       RETURN
    ENDIF
    eps  = rdmft_reg_eps
    nmin = MAX(nc, eps)
    pow_dg = alpha * nmin**(alpha - 1.0_DP)
  END FUNCTION pow_dg
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_g(n)
    !---------------------------------------------------------------
    !! Coupling function rdmft_g(n) for the separable / single-channel form.
    REAL(DP), INTENT(IN) :: n
    !
    SELECT CASE (rdmft_xc_id)
    CASE (RDMFT_XC_HF)
       rdmft_g = MAX(0.0_DP, n)
    CASE DEFAULT
       rdmft_g = pow_g(n, rdmft_power_alpha)
    END SELECT
  END FUNCTION rdmft_g
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_dg(n)
    !---------------------------------------------------------------
    REAL(DP), INTENT(IN) :: n
    !
    SELECT CASE (rdmft_xc_id)
    CASE (RDMFT_XC_HF)
       rdmft_dg = 1.0_DP
    CASE DEFAULT
       rdmft_dg = pow_dg(n, rdmft_power_alpha)
    END SELECT
  END FUNCTION rdmft_dg
  !
  !-----------------------------------------------------------------
  INTEGER FUNCTION rdmft_xc_n_channels()
    !---------------------------------------------------------------
    !! Number of separable channels into which the current functional
    !! decomposes (1 for HF / Muller / Power / GU; multiple for the
    !! mixed kernels CHF / CGA / GEO / HybOpt).
    SELECT CASE (rdmft_xc_id)
    CASE (RDMFT_XC_HF, RDMFT_XC_MULLER, RDMFT_XC_POWER, RDMFT_XC_GU)
       rdmft_xc_n_channels = 1
    CASE (RDMFT_XC_CHF, RDMFT_XC_CGA)
       rdmft_xc_n_channels = 2
    CASE (RDMFT_XC_GEO)
       rdmft_xc_n_channels = 3
    CASE (RDMFT_XC_HYBOPT)
       rdmft_xc_n_channels = 2
    CASE (RDMFT_XC_BOWMOD)
       rdmft_xc_n_channels = 4
    CASE (RDMFT_XC_BOW)
       rdmft_xc_n_channels = 0
    CASE DEFAULT
       rdmft_xc_n_channels = 1
    END SELECT
  END FUNCTION rdmft_xc_n_channels
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_channel(it, n, coef, w, dw)
    !---------------------------------------------------------------
    !! Evaluate channel \texttt{it} (1-based) at occupation \texttt{n}.
    !! \texttt{coef} is the constant prefactor in the energy sum,
    !! \texttt{w(n)} is the per-band weight used to build the modified
    !! density matrix \(\gamma_t\), and \texttt{dw(n)} is its derivative.
    INTEGER,  INTENT(IN)  :: it
    REAL(DP), INTENT(IN)  :: n
    REAL(DP), INTENT(OUT) :: coef, w, dw
    !
    REAL(DP) :: nc, a, s
    !
    nc = MAX(0.0_DP, MIN(1.0_DP, n))
    !
    SELECT CASE (rdmft_xc_id)
    CASE (RDMFT_XC_HF, RDMFT_XC_MULLER, RDMFT_XC_POWER, RDMFT_XC_GU)
       coef = 1.0_DP
       w    = rdmft_g(nc)
       dw   = rdmft_dg(nc)
    CASE (RDMFT_XC_CHF)
       SELECT CASE (it)
       CASE (1)
          coef = 0.5_DP
          w    = nc
          dw   = 1.0_DP
       CASE (2)
          coef = 0.5_DP
          w    = SQRT(MAX(0.0_DP, nc * (1.0_DP - nc)))
          s    = MAX(0.0_DP, nc * (1.0_DP - nc))
          IF (s <= rdmft_reg_eps * rdmft_reg_eps) THEN
             dw = 0.0_DP
          ELSE
             dw = 0.5_DP * (1.0_DP - 2.0_DP*nc) / w
          ENDIF
       END SELECT
    CASE (RDMFT_XC_CGA)
       SELECT CASE (it)
       CASE (1)
          coef = 0.25_DP
          w    = nc
          dw   = 1.0_DP
       CASE (2)
          coef = 0.25_DP
          w    = SQRT(MAX(0.0_DP, nc * (2.0_DP - nc)))
          s    = MAX(0.0_DP, nc * (2.0_DP - nc))
          IF (s <= rdmft_reg_eps * rdmft_reg_eps) THEN
             dw = 0.0_DP
          ELSE
             dw = (1.0_DP - nc) / w
          ENDIF
       END SELECT
    CASE (RDMFT_XC_GEO)
       SELECT CASE (it)
       CASE (1) ; coef = 0.25_DP ; a = 1.0_DP
       CASE (2) ; coef = 0.25_DP ; a = 0.5_DP
       CASE (3) ; coef = 0.5_DP  ; a = 0.75_DP
       CASE DEFAULT
          coef = 0.0_DP ; a = 1.0_DP
       END SELECT
       w  = pow_g(nc, a)
       dw = pow_dg(nc, a)
    CASE (RDMFT_XC_HYBOPT)
       SELECT CASE (it)
       CASE (1)
          coef = 1.0_DP - HYBOPT_POWER_W
          w    = nc
          dw   = 1.0_DP
       CASE (2)
          coef = HYBOPT_POWER_W
          w    = pow_g(nc, HYBOPT_POWER_EXP)
          dw   = pow_dg(nc, HYBOPT_POWER_EXP)
       END SELECT
    CASE (RDMFT_XC_BOWMOD)
       SELECT CASE (it)
       CASE (1)
          ! HF-like bridge channel; keeps idempotent occupancies close
          ! to the exact-exchange limit while preserving separability.
          coef = BOWMOD_COEF_HF
          w    = nc
          dw   = 1.0_DP
       CASE (2)
          coef = BOWMOD_COEF_ALPHA
          w    = pow_g(nc, rdmft_power_alpha)
          dw   = pow_dg(nc, rdmft_power_alpha)
       CASE (3)
          coef = BOWMOD_COEF_ALPHA
          w    = pow_g(1.0_DP - nc, rdmft_power_alpha)
          dw   = -pow_dg(1.0_DP - nc, rdmft_power_alpha)
       CASE (4)
          coef = BOWMOD_COEF_MIX
          w    = bowmod_mix_w(nc, BOWMOD_DEFAULT_BETA)
          dw   = bowmod_mix_dw(nc, BOWMOD_DEFAULT_BETA)
       CASE DEFAULT
          coef = 0.0_DP
          w    = 0.0_DP
          dw   = 0.0_DP
       END SELECT
    END SELECT
  END SUBROUTINE rdmft_xc_channel
  !
  !-----------------------------------------------------------------
  INTEGER FUNCTION rdmft_xc_eval_kind()
    !---------------------------------------------------------------
    SELECT CASE (rdmft_xc_id)
    CASE (RDMFT_XC_BOW)
       rdmft_xc_eval_kind = XC_EVAL_PAIR
    CASE DEFAULT
       rdmft_xc_eval_kind = XC_EVAL_CHANNELS
    END SELECT
  END FUNCTION rdmft_xc_eval_kind
  !
  !-----------------------------------------------------------------
  LOGICAL FUNCTION rdmft_xc_has_extra_diag()
    !---------------------------------------------------------------
    rdmft_xc_has_extra_diag = (rdmft_xc_id == RDMFT_XC_GU)
  END FUNCTION rdmft_xc_has_extra_diag
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION pair_kernel_f(ni, nj)
    !---------------------------------------------------------------
    REAL(DP), INTENT(IN) :: ni, nj
    SELECT CASE (rdmft_xc_id)
    CASE (RDMFT_XC_BOW)
       pair_kernel_f = bow_f(ni, nj)
    CASE DEFAULT
       pair_kernel_f = 0.0_DP
    END SELECT
  END FUNCTION pair_kernel_f
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION pair_kernel_df_dni(ni, nj)
    !---------------------------------------------------------------
    REAL(DP), INTENT(IN) :: ni, nj
    SELECT CASE (rdmft_xc_id)
    CASE (RDMFT_XC_BOW)
       pair_kernel_df_dni = bow_df_dni(ni, nj)
    CASE DEFAULT
       pair_kernel_df_dni = 0.0_DP
    END SELECT
  END FUNCTION pair_kernel_df_dni
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_set_exx_wg()
    !---------------------------------------------------------------
    USE klist, ONLY : nks, wk
    USE wvfct, ONLY : nbnd, wg
    INTEGER :: it, ib, ik, nch
    REAL(DP) :: coef, w, dw, wch
    !
    IF (rdmft_xc_eval_kind() == XC_EVAL_PAIR) THEN
       DO ik = 1, nks
          DO ib = 1, nbnd
             wg(ib, ik) = wk(ik) * rdmft_n(ib, ik)
          ENDDO
       ENDDO
       RETURN
    ENDIF
    nch = rdmft_xc_n_channels()
    DO ik = 1, nks
       DO ib = 1, nbnd
          wg(ib, ik) = 0.0_DP
          DO it = 1, nch
             CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
             wch = wk(ik) * w
             IF (wch > wg(ib, ik)) wg(ib, ik) = wch
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE rdmft_xc_set_exx_wg
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_set_wg_for_channel(it)
    !---------------------------------------------------------------
    USE klist, ONLY : nks, wk
    USE wvfct, ONLY : nbnd, wg
    INTEGER, INTENT(IN) :: it
    REAL(DP) :: coef, w, dw
    INTEGER :: ib, ik
    !
    DO ik = 1, nks
       DO ib = 1, nbnd
          CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
          wg(ib, ik) = wk(ik) * w
       ENDDO
    ENDDO
  END SUBROUTINE rdmft_xc_set_wg_for_channel
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_xc_channel_energy_term(it, vx_diag, wg)
    !---------------------------------------------------------------
    !! ``coef_it * 0.5 * sum(wg * vx_diag)`` for channel ``it``.
    USE wvfct, ONLY : nbnd
    USE klist, ONLY : nks
    INTEGER, INTENT(IN) :: it
    REAL(DP), INTENT(IN) :: vx_diag(:,:), wg(:,:)
    REAL(DP) :: coef, w, dw
    !
    CALL rdmft_xc_channel(it, 0.5_DP, coef, w, dw)
    rdmft_xc_channel_energy_term = coef * 0.5_DP &
         * SUM(wg(1:nbnd, 1:nks) * vx_diag(1:nbnd, 1:nks))
  END FUNCTION rdmft_xc_channel_energy_term
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_refresh_x_occupation()
    !---------------------------------------------------------------
    USE klist,    ONLY : nks, nkstot, wk
    USE wvfct,    ONLY : nbnd, wg
    USE exx_base, ONLY : x_occupation, x_nbnd_occ
    REAL(DP), ALLOCATABLE :: occ(:,:)
    INTEGER :: ik, ibnd
    REAL(DP), PARAMETER :: eps_occ = 1.0e-8_DP
    !
    IF (.NOT. ALLOCATED(x_occupation)) THEN
       ALLOCATE(x_occupation(nbnd, nkstot))
       x_occupation = 0.0_DP
    ENDIF
    ALLOCATE(occ(nbnd, nks))
    DO ik = 1, nks
       IF (ABS(wk(ik)) > eps_occ) THEN
          occ(1:nbnd, ik) = wg(1:nbnd, ik) / wk(ik)
       ELSE
          occ(1:nbnd, ik) = 0.0_DP
       ENDIF
    ENDDO
    CALL poolcollect(nbnd, nks, occ, nkstot, x_occupation)
    DEALLOCATE(occ)
    x_nbnd_occ = 0
    DO ik = 1, nkstot
       DO ibnd = MAX(1, x_nbnd_occ), nbnd
          IF (ABS(x_occupation(ibnd, ik)) > eps_occ) x_nbnd_occ = ibnd
       ENDDO
    ENDDO
  END SUBROUTINE rdmft_xc_refresh_x_occupation
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_reinit_exxbuff()
    !---------------------------------------------------------------
    USE exx, ONLY : exxinit, exxbuff, exxbuff_d
    USE exx_base, ONLY : x_nbnd_occ
    !
    CALL rdmft_xc_set_exx_wg()
    IF (ALLOCATED(exxbuff))   DEALLOCATE(exxbuff)
    IF (ALLOCATED(exxbuff_d)) DEALLOCATE(exxbuff_d)
    CALL exxinit(.FALSE.)
    rdmft_exxbuff_nbnd_occ = x_nbnd_occ
    rdmft_exxbuff_stale = .FALSE.
    rdmft_xi_cache_valid = .FALSE.
    rdmft_xi_cache_channel = 0
  END SUBROUTINE rdmft_xc_reinit_exxbuff
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_occ_grad_add(grad, vx_diag, ib, ik, wk, n, it)
    !---------------------------------------------------------------
    REAL(DP), INTENT(INOUT) :: grad(:,:)
    REAL(DP), INTENT(IN)    :: vx_diag(:,:)
    INTEGER, INTENT(IN)     :: ib, ik
    REAL(DP), INTENT(IN)    :: wk, n
    INTEGER, INTENT(IN), OPTIONAL :: it
    REAL(DP) :: coef, w, dw
    INTEGER :: ich
    !
    IF (rdmft_xc_eval_kind() == XC_EVAL_PAIR) THEN
       grad(ib, ik) = grad(ib, ik) + wk * vx_diag(ib, ik)
    ELSE
       ich = 1
       IF (PRESENT(it)) ich = it
       CALL rdmft_xc_channel(ich, n, coef, w, dw)
       grad(ib, ik) = grad(ib, ik) + coef * wk * dw * vx_diag(ib, ik)
    ENDIF
  END SUBROUTINE rdmft_xc_occ_grad_add
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_occ_grad_add_from_diag(grad, diag_vx, ib, ik, wk, n, it)
    !---------------------------------------------------------------
    !! Occupation gradient when ``diag_vx`` already stores
    !! ``coef * <psi|Vx|psi>`` (cached band diagnostics).
    REAL(DP), INTENT(INOUT) :: grad(:,:)
    REAL(DP), INTENT(IN)    :: diag_vx(:,:)
    INTEGER, INTENT(IN)     :: ib, ik
    REAL(DP), INTENT(IN)    :: wk, n
    INTEGER, INTENT(IN), OPTIONAL :: it
    REAL(DP) :: coef, w, dw
    INTEGER :: ich
    !
    ich = 1
    IF (PRESENT(it)) ich = it
    CALL rdmft_xc_channel(ich, n, coef, w, dw)
    grad(ib, ik) = grad(ib, ik) + wk * dw * diag_vx(ib, ik)
  END SUBROUTINE rdmft_xc_occ_grad_add_from_diag
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_extra_energy(e_add, vx_diag)
    !---------------------------------------------------------------
    USE wvfct,    ONLY : nbnd
    USE klist,    ONLY : nks, wk
    USE mp,       ONLY : mp_sum
    USE mp_pools, ONLY : inter_pool_comm
    REAL(DP), INTENT(OUT) :: e_add
    REAL(DP), INTENT(IN)  :: vx_diag(:,:)
    REAL(DP) :: fac_gu
    INTEGER :: ik, ib
    !
    e_add = 0.0_DP
    IF (.NOT. rdmft_xc_has_extra_diag()) RETURN
    fac_gu = 0.0_DP
    DO ik = 1, nks
       DO ib = 1, nbnd
          fac_gu = fac_gu + wk(ik)**2 * rdmft_gu_diag_factor(rdmft_n(ib, ik)) &
                   * vx_diag(ib, ik)
       ENDDO
    ENDDO
    CALL mp_sum(fac_gu, inter_pool_comm)
    e_add = 0.5_DP * fac_gu
  END SUBROUTINE rdmft_xc_extra_energy
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_extra_grad_add(grad, vx_diag, ib, ik, wk, n)
    !---------------------------------------------------------------
    REAL(DP), INTENT(INOUT) :: grad(:,:)
    REAL(DP), INTENT(IN)    :: vx_diag(:,:)
    INTEGER, INTENT(IN)     :: ib, ik
    REAL(DP), INTENT(IN)    :: wk, n
    !
    IF (.NOT. rdmft_xc_has_extra_diag()) RETURN
    grad(ib, ik) = grad(ib, ik) + 0.5_DP * wk**2 &
         * rdmft_gu_diag_factor_deriv(n) * vx_diag(ib, ik)
  END SUBROUTINE rdmft_xc_extra_grad_add
  !
  !-----------------------------------------------------------------
  REAL(DP) FUNCTION rdmft_xc_diag_vx_weight(ib, ik, n, it)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN) :: ib, ik
    REAL(DP), INTENT(IN) :: n
    INTEGER, INTENT(IN), OPTIONAL :: it
    REAL(DP) :: coef, w, dw
    INTEGER :: ich
    !
    IF (rdmft_xc_eval_kind() == XC_EVAL_PAIR) THEN
       rdmft_xc_diag_vx_weight = 1.0_DP
    ELSE
       ich = 1
       IF (PRESENT(it)) ich = it
       CALL rdmft_xc_channel(ich, n, coef, w, dw)
       rdmft_xc_diag_vx_weight = coef
    ENDIF
  END FUNCTION rdmft_xc_diag_vx_weight
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_probe_exchange_add(contrib, vx_diag, ib, ik, n_probe, mode, it)
    !---------------------------------------------------------------
    USE klist, ONLY : wk
    REAL(DP), INTENT(INOUT) :: contrib
    REAL(DP), INTENT(IN)    :: vx_diag(:,:)
    INTEGER, INTENT(IN)     :: ib, ik
    REAL(DP), INTENT(IN)    :: n_probe
    CHARACTER(LEN=*), INTENT(IN) :: mode
    INTEGER, INTENT(IN), OPTIONAL :: it
    REAL(DP) :: coef, w, dw
    INTEGER :: ich, nch
    !
    IF (rdmft_xc_eval_kind() == XC_EVAL_PAIR) THEN
       IF (TRIM(mode) == 'tsm') THEN
          contrib = contrib + wk(ik) * vx_diag(ib, ik)
       ELSEIF (TRIM(mode) == 'koop') THEN
          contrib = contrib + vx_diag(ib, ik)
       ENDIF
       IF (rdmft_xc_has_extra_diag() .AND. TRIM(mode) == 'tsm') THEN
          contrib = contrib + 0.5_DP * wk(ik)**2 &
               * rdmft_gu_diag_factor_deriv(n_probe) * vx_diag(ib, ik)
       ENDIF
       RETURN
    ENDIF
    IF (PRESENT(it)) THEN
       CALL rdmft_xc_channel(it, n_probe, coef, w, dw)
       IF (TRIM(mode) == 'tsm') THEN
          contrib = contrib + coef * wk(ik) * dw * vx_diag(ib, ik)
       ELSEIF (TRIM(mode) == 'koop') THEN
          contrib = contrib + coef * vx_diag(ib, ik)
       ENDIF
       IF (rdmft_xc_has_extra_diag() .AND. TRIM(mode) == 'tsm' .AND. it == 1) THEN
          contrib = contrib + 0.5_DP * wk(ik)**2 &
               * rdmft_gu_diag_factor_deriv(n_probe) * vx_diag(ib, ik)
       ENDIF
    ELSE
       nch = rdmft_xc_n_channels()
       DO ich = 1, nch
          CALL rdmft_xc_channel(ich, n_probe, coef, w, dw)
          IF (TRIM(mode) == 'tsm') THEN
             contrib = contrib + coef * wk(ik) * dw * vx_diag(ib, ik)
          ELSEIF (TRIM(mode) == 'koop') THEN
             contrib = contrib + coef * vx_diag(ib, ik)
          ENDIF
       ENDDO
       IF (rdmft_xc_has_extra_diag() .AND. TRIM(mode) == 'tsm') THEN
          contrib = contrib + 0.5_DP * wk(ik)**2 &
               * rdmft_gu_diag_factor_deriv(n_probe) * vx_diag(ib, ik)
       ENDIF
    ENDIF
  END SUBROUTINE rdmft_xc_probe_exchange_add
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_compute_exchange_pair(e_t, vx_diag, vxpsi, only_k, ib_only)
    !---------------------------------------------------------------
    USE wvfct,             ONLY : nbnd, npwx, current_k
    USE klist,             ONLY : nks, nkstot, ngk, igk_k, xk, wk
    USE wavefunctions,     ONLY : evc
    USE noncollin_module,  ONLY : npol
    USE control_flags,     ONLY : gamma_only
    USE gvect,             ONLY : gstart
    USE io_files,          ONLY : nwordwfc, iunwfc
    USE buffers,           ONLY : get_buffer
    USE lsda_mod,          ONLY : lsda, current_spin, isk
    USE uspp,              ONLY : nkb, vkb, okvan
    USE uspp_init,         ONLY : init_us_2
    USE becmod,            ONLY : bec_type, allocate_bec_type, deallocate_bec_type, calbec
    USE exx,               ONLY : vexx
    USE exx_base,          ONLY : x_nbnd_occ, x_occupation
    USE mp,                ONLY : mp_sum
    USE mp_bands,          ONLY : intra_bgrp_comm
    USE mp_pools,          ONLY : inter_pool_comm
    REAL(DP), INTENT(OUT) :: e_t
    REAL(DP), INTENT(OUT) :: vx_diag(:,:)
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: vxpsi(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: only_k, ib_only
    COMPLEX(DP), ALLOCATABLE :: vx_loc(:,:)
    TYPE(bec_type) :: becpsi
    REAL(DP), ALLOCATABLE :: n_global(:,:), wk_global(:)
    REAL(DP) :: fpair, wsrc, wtgt, vdot
    INTEGER :: ikt, ib, jb, iks, npw, klo, khi, jb_lo, jb_hi, loc, m_vexx
    !
    IF (PRESENT(ib_only) .AND. .NOT. PRESENT(only_k)) &
         CALL errore('rdmft_xc_compute_exchange_pair', 'ib_only requires only_k', 1)
    IF (PRESENT(ib_only) .AND. PRESENT(vxpsi)) &
         CALL errore('rdmft_xc_compute_exchange_pair', 'ib_only incompatible with vxpsi', 2)
    IF (PRESENT(only_k)) THEN
       klo = only_k ; khi = only_k
    ELSE
       klo = 1 ; khi = nks
    ENDIF
    IF (PRESENT(ib_only)) THEN
       jb_lo = ib_only ; jb_hi = ib_only ; m_vexx = 1
    ELSE
       jb_lo = 1 ; jb_hi = nbnd ; m_vexx = nbnd
    ENDIF
    vx_diag = 0.0_DP
    IF (PRESENT(vxpsi)) vxpsi = (0.0_DP, 0.0_DP)
    e_t = 0.0_DP
    ALLOCATE(n_global(nbnd, nkstot), wk_global(nkstot))
    CALL poolcollect(nbnd, nks, rdmft_n, nkstot, n_global)
    CALL poolcollect(1, nks, wk, nkstot, wk_global)
    IF (.NOT. ALLOCATED(x_occupation)) ALLOCATE(x_occupation(nbnd, nkstot))
    x_occupation = 0.0_DP
    ALLOCATE(vx_loc(npwx*npol, m_vexx))
    DO ikt = klo, khi
       npw = ngk(ikt)
       current_k = ikt
       IF (lsda) current_spin = isk(ikt)
       IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ikt)
       CALL g2_kin(ikt)
       IF (okvan) THEN
          CALL init_us_2(npw, igk_k(1, ikt), xk(:, ikt), vkb)
          CALL allocate_bec_type(nkb, m_vexx, becpsi)
          CALL calbec(npw, vkb, evc(:, jb_lo:jb_hi), becpsi, m_vexx)
       ENDIF
       DO iks = 1, nkstot
          DO ib = 1, nbnd
             x_occupation(ib, iks) = 1.0_DP
             x_nbnd_occ = ib
             vx_loc = (0.0_DP, 0.0_DP)
             IF (okvan) THEN
                CALL vexx(npwx, npw, m_vexx, evc(:, jb_lo:jb_hi), vx_loc, becpsi)
             ELSE
                CALL vexx(npwx, npw, m_vexx, evc(:, jb_lo:jb_hi), vx_loc)
             ENDIF
             x_occupation(ib, iks) = 0.0_DP
             wsrc = wk_global(iks)
             wtgt = wk(ikt)
             DO jb = jb_lo, jb_hi
                loc = jb - jb_lo + 1
                IF (gamma_only) THEN
                   vdot = 2.0_DP * REAL(SUM(CONJG(evc(1:npw, jb)) &
                        * vx_loc(1:npw, loc)), KIND=DP)
                   IF (gstart == 2) vdot = vdot - &
                        REAL(evc(1, jb), KIND=DP) * REAL(vx_loc(1, loc), KIND=DP)
                ELSE
                   vdot = REAL(SUM(CONJG(evc(1:npw, jb)) &
                        * vx_loc(1:npw, loc)), KIND=DP)
                   IF (npol == 2) vdot = vdot + &
                        REAL(SUM(CONJG(evc(npwx+1:npwx+npw, jb)) &
                                   * vx_loc(npwx+1:npwx+npw, loc)), KIND=DP)
                ENDIF
                fpair = pair_kernel_f(n_global(jb, ikt), n_global(ib, iks))
                vx_diag(jb, ikt) = vx_diag(jb, ikt) &
                     + wsrc * pair_kernel_df_dni(n_global(jb, ikt), n_global(ib, iks)) * vdot
                IF (.NOT. PRESENT(ib_only)) &
                     e_t = e_t + 0.5_DP * wtgt * wsrc * fpair * vdot
                IF (PRESENT(vxpsi)) &
                     vxpsi(:, jb, ikt) = vxpsi(:, jb, ikt) + 0.5_DP * wsrc * fpair * vx_loc(:, loc)
             ENDDO
          ENDDO
       ENDDO
       IF (okvan) CALL deallocate_bec_type(becpsi)
    ENDDO
    DEALLOCATE(vx_loc, n_global, wk_global)
    CALL mp_sum(vx_diag, intra_bgrp_comm)
    IF (.NOT. PRESENT(ib_only)) THEN
       CALL mp_sum(e_t, inter_pool_comm)
    ELSE
       e_t = 0.0_DP
    ENDIF
    rdmft_xi_cache_valid = .FALSE.
    rdmft_xi_cache_channel = 0
  END SUBROUTINE rdmft_xc_compute_exchange_pair
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_compute_exchange_channels(it, e_t, vx_diag, vxpsi, only_k, ib_only)
    !---------------------------------------------------------------
    USE wvfct,             ONLY : nbnd, npwx, current_k, wg
    USE klist,             ONLY : nks, ngk, igk_k, xk
    USE wavefunctions,     ONLY : evc
    USE noncollin_module,  ONLY : npol
    USE control_flags,     ONLY : gamma_only
    USE gvect,             ONLY : gstart
    USE io_files,          ONLY : nwordwfc, iunwfc
    USE buffers,           ONLY : get_buffer
    USE lsda_mod,          ONLY : lsda, current_spin, isk
    USE uspp,              ONLY : nkb, vkb, okvan
    USE uspp_init,         ONLY : init_us_2
    USE becmod,            ONLY : bec_type, allocate_bec_type, deallocate_bec_type, calbec
    USE exx,               ONLY : aceinit, vexxace_gamma, vexxace_k, vexx, &
                                  use_ace, domat, xi
    USE mp,                ONLY : mp_sum, mp_min
    USE mp_bands,          ONLY : intra_bgrp_comm
    USE mp_images,         ONLY : intra_image_comm
    USE exx_base,          ONLY : x_nbnd_occ
    INTEGER, INTENT(IN)   :: it
    REAL(DP), INTENT(OUT) :: e_t
    REAL(DP), INTENT(OUT) :: vx_diag(:,:)
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: vxpsi(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: only_k, ib_only
    COMPLEX(DP), ALLOCATABLE :: vx_loc(:,:)
    INTEGER :: ik, npw, ib, klo, khi, ib_lo, ib_hi, m_vexx, xi_reuse_flag
    REAL(DP) :: dummy
    TYPE(bec_type) :: becpsi
    LOGICAL :: xi_reuse, direct_vexx_path
    REAL(DP), PARAMETER :: xi_cache_tol = 1.0e-12_DP
    !
    CALL rdmft_xc_set_wg_for_channel(it)
    ! Keep ``x_occupation`` globally consistent with the current channel
    ! weights before any ``vexx``/ACE call.  The direct ``only_k`` path
    ! still needs all k-points because the exchange kernel reads the
    ! full occupation table across ``inter_pool_comm``.
    CALL rdmft_xc_refresh_x_occupation()
    IF (.NOT. PRESENT(only_k)) THEN
       IF (rdmft_exxbuff_stale .OR. x_nbnd_occ > rdmft_exxbuff_nbnd_occ) THEN
          CALL rdmft_xc_reinit_exxbuff()
          CALL rdmft_xc_set_wg_for_channel(it)
          CALL rdmft_xc_refresh_x_occupation()
       ENDIF
    ENDIF
    IF (PRESENT(only_k)) THEN
       klo = only_k ; khi = only_k
    ELSE
       klo = 1 ; khi = nks
    ENDIF
    direct_vexx_path = PRESENT(only_k)
    IF (PRESENT(ib_only) .AND. .NOT. PRESENT(only_k)) &
         CALL errore('rdmft_xc_compute_exchange_channels', 'ib_only requires only_k', 1)
    IF (PRESENT(ib_only) .AND. PRESENT(vxpsi)) &
         CALL errore('rdmft_xc_compute_exchange_channels', 'ib_only incompatible with vxpsi', 2)
    IF (direct_vexx_path) THEN
       ik = only_k
       npw = ngk(ik)
       current_k = ik
       IF (lsda) current_spin = isk(ik)
       IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       CALL g2_kin(ik)
       IF (PRESENT(ib_only)) THEN
          ib_lo = ib_only ; ib_hi = ib_only ; m_vexx = 1
       ELSE
          ib_lo = 1 ; ib_hi = nbnd ; m_vexx = nbnd
       ENDIF
       IF (okvan) THEN
          CALL init_us_2(npw, igk_k(1, ik), xk(:, ik), vkb)
          CALL allocate_bec_type(nkb, m_vexx, becpsi)
          CALL calbec(npw, vkb, evc(:, ib_lo:ib_hi), becpsi, m_vexx)
       ENDIF
       vx_diag = 0.0_DP
       IF (PRESENT(vxpsi)) vxpsi = (0.0_DP, 0.0_DP)
       ALLOCATE(vx_loc(npwx*npol, m_vexx))
       vx_loc = (0.0_DP, 0.0_DP)
       IF (okvan) THEN
          CALL vexx(npwx, npw, m_vexx, evc(:, ib_lo:ib_hi), vx_loc, becpsi)
       ELSE
          CALL vexx(npwx, npw, m_vexx, evc(:, ib_lo:ib_hi), vx_loc)
       ENDIF
       DO ib = ib_lo, ib_hi
          IF (gamma_only) THEN
             vx_diag(ib, ik) = 2.0_DP * REAL(SUM(CONJG(evc(1:npw, ib)) * &
                  vx_loc(1:npw, ib - ib_lo + 1)), KIND=DP)
             IF (gstart == 2) THEN
                vx_diag(ib, ik) = vx_diag(ib, ik) - &
                     REAL(evc(1, ib), KIND=DP) * &
                     REAL(vx_loc(1, ib - ib_lo + 1), KIND=DP)
             ENDIF
          ELSE
             vx_diag(ib, ik) = REAL(SUM(CONJG(evc(1:npw, ib)) * &
                  vx_loc(1:npw, ib - ib_lo + 1)), KIND=DP)
             IF (npol == 2) THEN
                vx_diag(ib, ik) = vx_diag(ib, ik) + &
                     REAL(SUM(CONJG(evc(npwx+1:npwx+npw, ib)) * &
                                vx_loc(npwx+1:npwx+npw, ib - ib_lo + 1)), KIND=DP)
             ENDIF
          ENDIF
       ENDDO
       IF (PRESENT(vxpsi)) vxpsi(:, ib_lo:ib_hi, ik) = vx_loc(:, 1:m_vexx)
       CALL mp_sum(vx_diag, intra_bgrp_comm)
       IF (okvan) CALL deallocate_bec_type(becpsi)
       DEALLOCATE(vx_loc)
       rdmft_xi_cache_valid = .FALSE.
       rdmft_xi_cache_channel = 0
       e_t = 0.0_DP
       RETURN
    ELSE
       xi_reuse = .FALSE.
       IF (ALLOCATED(xi) .AND. rdmft_xi_cache_valid &
           .AND. rdmft_xi_cache_channel == it &
           .AND. ALLOCATED(rdmft_xi_cache_n) &
           .AND. SIZE(rdmft_xi_cache_n, 1) == nbnd &
           .AND. SIZE(rdmft_xi_cache_n, 2) == SIZE(rdmft_n, 2)) THEN
          IF (ALL(ABS(rdmft_n - rdmft_xi_cache_n) < xi_cache_tol)) xi_reuse = .TRUE.
       ENDIF
       ! ``xi_reuse`` is decided from the POOL-LOCAL occupation block
       ! ``rdmft_n(:,1:nks)``, so two pools can disagree (e.g. near
       ! convergence one pool's k-points stop moving while a
       ! Fermi-surface pool still updates).  The reuse branch skips
       ! ``aceinit`` -- which performs an ``mp_sum`` over
       ! ``inter_pool_comm`` -- so a disagreement makes some ranks call
       ! the collective while others skip it and the run DEADLOCKS
       ! (this is the ``mpirun -nk > 1`` hang).  Force a single global
       ! decision: rebuild the ACE projector everywhere unless EVERY
       ! rank in the image can reuse its cache.
       xi_reuse_flag = 0
       IF (xi_reuse) xi_reuse_flag = 1
       CALL mp_min(xi_reuse_flag, intra_image_comm)
       xi_reuse = (xi_reuse_flag == 1)
       IF (xi_reuse) THEN
          e_t = 0.0_DP
          domat = .FALSE.
       ELSE
          CALL aceinit(.FALSE., e_t)
          IF (.NOT. ALLOCATED(rdmft_xi_cache_n)) THEN
             ALLOCATE(rdmft_xi_cache_n(nbnd, SIZE(rdmft_n, 2)))
          ELSE IF (SIZE(rdmft_xi_cache_n, 1) /= nbnd .OR. &
                   SIZE(rdmft_xi_cache_n, 2) /= SIZE(rdmft_n, 2)) THEN
             DEALLOCATE(rdmft_xi_cache_n)
             ALLOCATE(rdmft_xi_cache_n(nbnd, SIZE(rdmft_n, 2)))
          ENDIF
          rdmft_xi_cache_n = rdmft_n
          rdmft_xi_cache_channel = it
          rdmft_xi_cache_valid = .TRUE.
       ENDIF
    ENDIF
    vx_diag = 0.0_DP
    IF (PRESENT(vxpsi)) vxpsi = (0.0_DP, 0.0_DP)
    ALLOCATE(vx_loc(npwx*npol, nbnd))
    DO ik = klo, khi
       npw = ngk(ik)
       current_k = ik
       IF (lsda) current_spin = isk(ik)
       IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       CALL g2_kin(ik)
       IF (okvan) CALL init_us_2(npw, igk_k(1, ik), xk(:, ik), vkb)
       vx_loc = (0.0_DP, 0.0_DP)
       IF (gamma_only) THEN
          CALL vexxace_gamma(npw, nbnd, evc, dummy, vx_loc)
       ELSE
          CALL vexxace_k(npw, nbnd, evc, dummy, vx_loc)
       ENDIF
!$omp parallel do default(shared) private(ib)
       DO ib = 1, nbnd
          IF (gamma_only) THEN
             vx_diag(ib, ik) = 2.0_DP * REAL(SUM(CONJG(evc(1:npw, ib)) * vx_loc(1:npw, ib)), KIND=DP)
             IF (gstart == 2) THEN
                vx_diag(ib, ik) = vx_diag(ib, ik) - &
                     REAL(evc(1, ib), KIND=DP) * REAL(vx_loc(1, ib), KIND=DP)
             ENDIF
          ELSE
             vx_diag(ib, ik) = REAL(SUM(CONJG(evc(1:npw, ib)) * vx_loc(1:npw, ib)), KIND=DP)
             IF (npol == 2) THEN
                vx_diag(ib, ik) = vx_diag(ib, ik) + &
                     REAL(SUM(CONJG(evc(npwx+1:npwx+npw, ib)) * &
                                vx_loc(npwx+1:npwx+npw, ib)), KIND=DP)
             ENDIF
          ENDIF
       ENDDO
!$omp end parallel do
       IF (PRESENT(vxpsi)) vxpsi(:, :, ik) = vx_loc(:, :)
    ENDDO
    CALL mp_sum(vx_diag, intra_bgrp_comm)
    DEALLOCATE(vx_loc)
  END SUBROUTINE rdmft_xc_compute_exchange_channels
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_compute_exchange(it, e_t, vx_diag, vxpsi, only_k, ib_only)
    !---------------------------------------------------------------
    INTEGER, INTENT(IN)   :: it
    REAL(DP), INTENT(OUT) :: e_t
    REAL(DP), INTENT(OUT) :: vx_diag(:,:)
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: vxpsi(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: only_k, ib_only
    !
    IF (rdmft_xc_eval_kind() == XC_EVAL_PAIR) THEN
       CALL rdmft_xc_compute_exchange_pair(e_t, vx_diag, vxpsi, only_k, ib_only)
    ELSE
       CALL rdmft_xc_compute_exchange_channels(it, e_t, vx_diag, vxpsi, only_k, ib_only)
    ENDIF
  END SUBROUTINE rdmft_xc_compute_exchange
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_compute_exchange_kopt(it, e_t, vx_diag, only_k, ib_only)
    !---------------------------------------------------------------
    !! Forward ``only_k`` / ``ib_only`` to :func:`rdmft_xc_compute_exchange`.
    INTEGER, INTENT(IN) :: it
    REAL(DP), INTENT(OUT) :: e_t
    REAL(DP), INTENT(OUT) :: vx_diag(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: only_k, ib_only
    !
    IF (PRESENT(only_k)) THEN
       IF (PRESENT(ib_only)) THEN
          CALL rdmft_xc_compute_exchange(it, e_t, vx_diag, only_k=only_k, ib_only=ib_only)
       ELSE
          CALL rdmft_xc_compute_exchange(it, e_t, vx_diag, only_k=only_k)
       ENDIF
    ELSE
       CALL rdmft_xc_compute_exchange(it, e_t, vx_diag)
    ENDIF
  END SUBROUTINE rdmft_xc_compute_exchange_kopt
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_add_occ_gradient_exact_kopt(grad_n, only_k, ib_only)
    !---------------------------------------------------------------
    !! Forward ``only_k`` / ``ib_only`` to
    !! :func:`rdmft_xc_add_occ_gradient_exact`.
    REAL(DP), INTENT(INOUT) :: grad_n(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: only_k, ib_only
    !
    IF (PRESENT(only_k)) THEN
       IF (PRESENT(ib_only)) THEN
          CALL rdmft_xc_add_occ_gradient_exact(grad_n, only_k=only_k, ib_only=ib_only)
       ELSE
          CALL rdmft_xc_add_occ_gradient_exact(grad_n, only_k=only_k)
       ENDIF
    ELSE
       CALL rdmft_xc_add_occ_gradient_exact(grad_n)
    ENDIF
  END SUBROUTINE rdmft_xc_add_occ_gradient_exact_kopt
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_orb_exchange(e_xc, vxpsi_total, vx_diag, vxpsi_buf)
    !---------------------------------------------------------------
    USE wvfct,    ONLY : nbnd, npwx, wg
    USE klist,    ONLY : nks
    USE mp,       ONLY : mp_sum
    USE mp_pools, ONLY : inter_pool_comm
    REAL(DP), INTENT(OUT) :: e_xc
    COMPLEX(DP), INTENT(OUT) :: vxpsi_total(:,:,:)
    REAL(DP), INTENT(OUT) :: vx_diag(:,:)
    COMPLEX(DP), INTENT(INOUT) :: vxpsi_buf(:,:,:)
    REAL(DP) :: coef, w, dw, e_t, e_xc_partial, e_add
    INTEGER :: it, nch, ik, ib
    !
    vxpsi_total = (0.0_DP, 0.0_DP)
    e_xc = 0.0_DP
    IF (rdmft_xc_eval_kind() == XC_EVAL_PAIR) THEN
       CALL rdmft_xc_compute_exchange_pair(e_t, vx_diag, vxpsi_buf)
       e_xc = e_t
       vxpsi_total = vxpsi_buf
    ELSE
       nch = rdmft_xc_n_channels()
       DO it = 1, nch
          CALL rdmft_xc_compute_exchange_channels(it, e_t, vx_diag, vxpsi_buf)
          CALL rdmft_xc_channel(it, 0.5_DP, coef, w, dw)
          e_xc_partial = coef * 0.5_DP * SUM(wg(1:nbnd, 1:nks) * vx_diag(1:nbnd, 1:nks))
          CALL mp_sum(e_xc_partial, inter_pool_comm)
          e_xc = e_xc + e_xc_partial
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
    CALL rdmft_xc_extra_energy(e_add, vx_diag)
    e_xc = e_xc + e_add
  END SUBROUTINE rdmft_xc_orb_exchange
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_exxinit_once()
    !---------------------------------------------------------------
    !! Run \texttt{exxinit} once before the RDMFT loop (sizes
    !! \texttt{exxbuff} from the widest channel support).
    CALL rdmft_xc_reinit_exxbuff()
  END SUBROUTINE rdmft_xc_exxinit_once
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_ensure_exxbuff_tsm_capacity()
    !---------------------------------------------------------------
    !! Grow ``exxbuff`` if the TSM probe pattern (\(n_{ik}=1/2\) on
    !! every band) needs a wider occupied-band bound than the current
    !! buffer.
    USE wvfct,        ONLY : nbnd
    USE klist,        ONLY : nks
    USE exx_base,     ONLY : x_nbnd_occ
    !
    REAL(DP), ALLOCATABLE :: n_save(:,:)
    INTEGER :: ik, ib
    REAL(DP), PARAMETER :: n_probe = 0.5_DP
    !
    ALLOCATE(n_save(nbnd, nks))
    n_save = rdmft_n
    DO ik = 1, nks
       DO ib = 1, nbnd
          rdmft_n(ib, ik) = n_probe
       ENDDO
    ENDDO
    CALL rdmft_xc_set_exx_wg()
    CALL rdmft_xc_refresh_x_occupation()
    IF (x_nbnd_occ > rdmft_exxbuff_nbnd_occ) THEN
       rdmft_exxbuff_stale = .TRUE.
       CALL rdmft_xc_reinit_exxbuff()
    ENDIF
    rdmft_n = n_save
    DEALLOCATE(n_save)
    CALL rdmft_xc_set_exx_wg()
    CALL rdmft_xc_refresh_x_occupation()
  END SUBROUTINE rdmft_xc_ensure_exxbuff_tsm_capacity
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_ensure_min_exx_occ_support(it, ik)
    !---------------------------------------------------------------
    !! Inflate ``wg`` at ``ik`` so ACE has at least two active bands.
    USE klist,        ONLY : wk
    USE wvfct,        ONLY : nbnd, wg
    !
    INTEGER, INTENT(IN) :: it, ik
    !
    INTEGER :: ib, n_active
    REAL(DP) :: coef, w, dw, w_probe, w_min
    REAL(DP), PARAMETER :: occ_eps = 1.0e-8_DP
    INTEGER, PARAMETER :: min_active = 2
    !
    CALL rdmft_xc_channel(it, 0.5_DP, coef, w_probe, dw)
    w_min = wk(ik) * w_probe
    n_active = 0
    DO ib = 1, nbnd
       CALL rdmft_xc_channel(it, rdmft_n(ib, ik), coef, w, dw)
       IF (wk(ik) * w > occ_eps) n_active = n_active + 1
    ENDDO
    IF (n_active >= min_active) RETURN
    DO ib = 1, nbnd
       IF (wg(ib, ik) < w_min) wg(ib, ik) = w_min
    ENDDO
    CALL rdmft_xc_refresh_x_occupation()
  END SUBROUTINE rdmft_xc_ensure_min_exx_occ_support
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_tsm_trim_wg_for_ace(it, ik, ib_probe)
    !---------------------------------------------------------------
    !! Drop nearly empty bands from ``wg`` at ``ik`` for TSM ACE builds.
    USE klist,        ONLY : wk
    USE wvfct,        ONLY : nbnd, wg
    !
    INTEGER, INTENT(IN) :: it, ik, ib_probe
    !
    INTEGER :: jb
    REAL(DP) :: n_thresh
    !
    IF (ib_probe < 1 .OR. ib_probe > nbnd) RETURN
    IF (ik < 1) RETURN
    IF (ABS(wk(ik)) <= 0.0_DP) RETURN
    !
    n_thresh = rdmft_reg_eps
    DO jb = 1, nbnd
       IF (jb == ib_probe) CYCLE
       IF (rdmft_n(jb, ik) < n_thresh) wg(jb, ik) = 0.0_DP
    ENDDO
    CALL rdmft_xc_ensure_min_exx_occ_support(it, ik)
  END SUBROUTINE rdmft_xc_tsm_trim_wg_for_ace
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_total_energy(e_xc_total)
    !---------------------------------------------------------------
    !! RDMFT XC exchange energy (+ GU on-site extra when active).
    USE wvfct,    ONLY : nbnd, wg
    USE klist,    ONLY : nks
    USE mp,       ONLY : mp_sum
    USE mp_pools, ONLY : inter_pool_comm
    !
    REAL(DP), INTENT(OUT) :: e_xc_total
    !
    REAL(DP) :: e_t, e_xc_partial, e_add
    REAL(DP), ALLOCATABLE :: vx_diag(:,:)
    INTEGER :: it, nch
    !
    ALLOCATE(vx_diag(nbnd, nks))
    e_xc_total = 0.0_DP
    nch = rdmft_xc_n_channels()
    IF (nch == 0) THEN
       CALL rdmft_xc_compute_exchange(1, e_t, vx_diag)
       e_xc_total = e_t
       ! ``rdmft_xc_compute_exchange_pair`` already performs the
       ! inter-pool reduction for ``e_t`` on the pair-kernel path.
    ELSE
       DO it = 1, nch
          CALL rdmft_xc_compute_exchange(it, e_t, vx_diag)
          e_xc_partial = rdmft_xc_channel_energy_term(it, vx_diag, wg)
          CALL mp_sum(e_xc_partial, inter_pool_comm)
          e_xc_total = e_xc_total + e_xc_partial
       ENDDO
    ENDIF
    CALL rdmft_xc_extra_energy(e_add, vx_diag)
    e_xc_total = e_xc_total + e_add
    DEALLOCATE(vx_diag)
  END SUBROUTINE rdmft_xc_total_energy
  !
  !-----------------------------------------------------------------
  FUNCTION rdmft_xc_vdot_band(npw, ib, vx_loc) RESULT(vdot)
    !---------------------------------------------------------------
    !! \(\Re\langle\psi_{ib}|V_x|\psi_{ib}\rangle\) from ``vexx`` output.
    USE wvfct,            ONLY : npwx
    USE wavefunctions,    ONLY : evc
    USE noncollin_module, ONLY : npol
    USE control_flags,    ONLY : gamma_only
    USE gvect,            ONLY : gstart
    INTEGER, INTENT(IN) :: npw, ib
    COMPLEX(DP), INTENT(IN) :: vx_loc(:,:)
    REAL(DP) :: vdot
    !
    IF (gamma_only) THEN
       vdot = 2.0_DP * REAL(SUM(CONJG(evc(1:npw, ib)) * vx_loc(1:npw, ib)), KIND=DP)
       IF (gstart == 2) vdot = vdot - &
            REAL(evc(1, ib), KIND=DP) * REAL(vx_loc(1, ib), KIND=DP)
    ELSE
       vdot = REAL(SUM(CONJG(evc(1:npw, ib)) * vx_loc(1:npw, ib)), KIND=DP)
       IF (npol == 2) vdot = vdot + &
            REAL(SUM(CONJG(evc(npwx+1:npwx+npw, ib)) &
                       * vx_loc(npwx+1:npwx+npw, ib)), KIND=DP)
    ENDIF
  END FUNCTION rdmft_xc_vdot_band
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_add_occ_gradient_exact(grad_n, only_k, ib_only)
    !---------------------------------------------------------------
    !! \textbf{DEPRECATED / BUGGY -- no longer called.}  Retained only
    !! as historical reference and for any out-of-tree callers; left
    !! unchanged on purpose so that ``git blame`` keeps pointing at
    !! the issue.  The current driver
    !! (\texttt{rdmft\_xc\_add\_occ\_gradient}) ignores this routine
    !! and uses the closed-form contraction
    !! \(c_t w_{\mathbf{k}} w_t'(n)\,D_{ii}\), which matches the
    !! ABACUS reference and the derivation in
    !! \texttt{doc/rdmft\_calculation.md} \S9.4.
    !!
    !! \textbf{Original (incorrect) claim.}  The routine was meant to
    !! evaluate
    !!
    !! \(\partial E_x / \partial n_{j\mathbf{k}} = \tfrac{1}{2} c_t w_{\mathbf{k}}
    !! w_t'(n_{j\mathbf{k}}) \bigl[ D_{jj} + \sum_{\mathbf{k}',i}
    !! w_{\mathbf{k}'} w_t(n_{i\mathbf{k}'}) K_{i\mathbf{k}',j\mathbf{k}} \bigr]\)
    !!
    !! and combine an ``explicit'' \(D_{jj}/2\) with an ``implicit''
    !! response built from rank-one ``vexx'' calls.  Two problems:
    !!
    !! (a) The ``rank-one'' ``vexx'' diagonal at fixed unit
    !! \texttt{x\_occupation} of a single \((i_s,\mathbf{k}_s)\) source
    !! folds in the BZ-star weight \(w_{\mathbf{k}_s}/s_d\) coming from
    !! the per-q ``\(\times 1/n_{qs}\)'' factor in
    !! \texttt{exx.f90}, so the ``response'' sum equals
    !! \(w_{\mathbf{k}_s} D_{i_s i_s}\), not the dimensionless
    !! \(D_{i_s i_s}\) the formula assumes.
    !! (b) Even with a hypothetical ``\(D_{i_s i_s}\)'' response, the
    !! quadratic form \(E_x = -\tfrac{1}{2} g^T K g\) with symmetric
    !! \(K\) has gradient \(D_\alpha = -(Kg)_\alpha\) directly: there
    !! is no separate implicit-chain-rule piece to add.
    !!
    !! Net effect: this path overestimates the XC occupation gradient
    !! by a factor \((1+w_{\mathbf{k}})/2\) per source k-point
    !! (\(\times 3/2\) on \(\Gamma\)-only \(n_\mathrm{spin}=1\) cells
    !! where \(w_{\mathbf{k}}=2\); \(\times 0.625\) on a cubic
    !! \(2\times 2\times 2\) mesh where \(w_{\mathbf{k}}=0.25\); etc).
    !! The bug went unnoticed because ``rdmft\_grad\_check`` was
    !! off by default in every shipped example and HF runs that start
    !! at the KS-HF fixed point converge regardless (the line search
    !! absorbs the magnitude error).
    USE rdmft_module,      ONLY : rdmft_n, rdmft_exxbuff_stale
    USE wvfct,             ONLY : nbnd, npwx, current_k
    USE klist,             ONLY : nks, nkstot, ngk, igk_k, xk, wk
    USE wavefunctions,     ONLY : evc
    USE noncollin_module,  ONLY : npol
    USE control_flags,     ONLY : gamma_only
    USE gvect,             ONLY : gstart
    USE io_files,          ONLY : nwordwfc, iunwfc
    USE buffers,           ONLY : get_buffer
    USE lsda_mod,          ONLY : lsda, current_spin, isk
    USE uspp,              ONLY : nkb, vkb, okvan
    USE uspp_init,         ONLY : init_us_2
    USE becmod,            ONLY : bec_type, allocate_bec_type, deallocate_bec_type, calbec
    USE exx,               ONLY : vexx
    USE exx_base,          ONLY : x_nbnd_occ, x_occupation
    USE mp,                ONLY : mp_sum
    USE mp_bands,          ONLY : intra_bgrp_comm
    USE mp_pools,          ONLY : inter_pool_comm
    !
    REAL(DP), INTENT(INOUT) :: grad_n(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: only_k, ib_only
    !
    REAL(DP), ALLOCATABLE :: vx_diag(:,:), n_global(:,:), wk_global(:), x_occ_save(:,:)
    COMPLEX(DP), ALLOCATABLE :: vx_loc(:,:)
    REAL(DP) :: e_t, coef, w_src, dw_src, w_tgt, dw_dummy, pref, response, vdot
    TYPE(bec_type) :: becpsi
    INTEGER :: it, nch, iks, iks_lo, iks_hi, ibs, ibs_lo, ibs_hi
    INTEGER :: ikt, ikt_g, ik_src, npw, jb
    INTEGER, EXTERNAL :: global_kpoint_index, local_kpoint_index
    REAL(DP), PARAMETER :: dw_eps = 1.0e-30_DP
    !
    IF (PRESENT(ib_only) .AND. .NOT. PRESENT(only_k)) &
         CALL errore('rdmft_xc_add_occ_gradient_exact', 'ib_only requires only_k', 1)
    !
    nch = rdmft_xc_n_channels()
    IF (nch <= 0) RETURN
    !
    IF (rdmft_exxbuff_stale) CALL rdmft_xc_reinit_exxbuff()
    !
    ALLOCATE(n_global(nbnd, nkstot), wk_global(nkstot))
    CALL poolcollect(nbnd, nks, rdmft_n, nkstot, n_global)
    CALL poolcollect(1, nks, wk, nkstot, wk_global)
    IF (.NOT. ALLOCATED(x_occupation)) ALLOCATE(x_occupation(nbnd, nkstot))
    ALLOCATE(x_occ_save(nbnd, nkstot))
    x_occ_save = x_occupation
    !
    IF (PRESENT(only_k)) THEN
       iks_lo = global_kpoint_index(nkstot, only_k)
       iks_hi = iks_lo
    ELSE
       iks_lo = 1
       iks_hi = nkstot
    ENDIF
    IF (PRESENT(ib_only)) THEN
       ibs_lo = ib_only
       ibs_hi = ib_only
    ELSE
       ibs_lo = 1
       ibs_hi = nbnd
    ENDIF
    !
    ALLOCATE(vx_diag(nbnd, nks), vx_loc(npwx*npol, nbnd))
    DO it = 1, nch
       CALL rdmft_xc_compute_exchange_kopt(it, e_t, vx_diag, only_k, ib_only)
       !
       DO iks = iks_lo, iks_hi
          DO ibs = ibs_lo, ibs_hi
             CALL rdmft_xc_channel(it, n_global(ibs, iks), coef, w_src, dw_src)
             IF (ABS(dw_src) < dw_eps) CYCLE
             response = 0.0_DP
             DO ikt = 1, nks
                ikt_g = global_kpoint_index(nkstot, ikt)
                npw = ngk(ikt)
                current_k = ikt
                IF (lsda) current_spin = isk(ikt)
                IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ikt)
                CALL g2_kin(ikt)
                x_occupation = 0.0_DP
                x_occupation(ibs, iks) = 1.0_DP
                x_nbnd_occ = ibs
                vx_loc = (0.0_DP, 0.0_DP)
                IF (okvan) THEN
                   CALL init_us_2(npw, igk_k(1, ikt), xk(:, ikt), vkb)
                   CALL allocate_bec_type(nkb, nbnd, becpsi)
                   CALL calbec(npw, vkb, evc, becpsi, nbnd)
                   CALL vexx(npwx, npw, nbnd, evc, vx_loc, becpsi)
                   CALL deallocate_bec_type(becpsi)
                ELSE
                   CALL vexx(npwx, npw, nbnd, evc, vx_loc)
                ENDIF
                DO jb = 1, nbnd
                   CALL rdmft_xc_channel(it, n_global(jb, ikt_g), coef, w_tgt, dw_dummy)
                   vdot = rdmft_xc_vdot_band(npw, jb, vx_loc)
                   CALL mp_sum(vdot, intra_bgrp_comm)
                   response = response + wk(ikt) * w_tgt * vdot
                ENDDO
             ENDDO
             CALL mp_sum(response, inter_pool_comm)
             pref = 0.5_DP * coef * wk_global(iks) * dw_src
             ik_src = local_kpoint_index(nkstot, iks)
             IF (ik_src > 0) THEN
                grad_n(ibs, ik_src) = grad_n(ibs, ik_src) &
                     + pref * (vx_diag(ibs, ik_src) + response)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !
    x_occupation(1:nbnd, 1:nkstot) = x_occ_save
    CALL rdmft_xc_refresh_x_occupation()
    rdmft_xi_cache_valid = .FALSE.
    rdmft_xi_cache_channel = 0
    DEALLOCATE(vx_diag, vx_loc, n_global, wk_global, x_occ_save)
  END SUBROUTINE rdmft_xc_add_occ_gradient_exact
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_add_occ_gradient(grad_n, only_k, ib_only)
    !---------------------------------------------------------------
    !! Add the exchange part of \(\partial E/\partial n_{i\mathbf{k}}\)
    !! to ``grad_n``.
    !!
    !! \textbf{Derivation.} For separable channel functionals the
    !! exchange energy is the quadratic form
    !!
    !! \[ E_x[\gamma_t] = -\tfrac{1}{2}\sum_{\alpha\beta} g_\alpha
    !!    g_\beta K_{\alpha\beta} = \tfrac{1}{2}\sum_\alpha g_\alpha
    !!    D_\alpha,\qquad g_\alpha = w_{\mathbf{k}_\alpha}\,
    !!    w_t(n_\alpha), \quad D_\alpha = -\sum_\beta K_{\alpha\beta}
    !!    g_\beta = \langle\alpha|V_x[\gamma_t]|\alpha\rangle. \]
    !!
    !! Because \(K\) is symmetric, the gradient of the quadratic
    !! \(-\tfrac{1}{2} g^T K g\) reads \(\partial E_x/\partial g_\alpha
    !! = -(K g)_\alpha = D_\alpha\) -- the implicit chain rule through
    !! the Fock operator is \emph{automatically} absorbed into the
    !! single contraction.  Chain ruling to \(n\),
    !!
    !! \[ \frac{\partial E_x^{(t)}}{\partial n_{i\mathbf{k}}}
    !!    = w_{\mathbf{k}}\, w_t'(n_{i\mathbf{k}})\,
    !!      D_{ii}^{(t,\mathbf{k})}. \]
    !!
    !! This is exactly the ABACUS reference formula
    !! (\texttt{rdmft\_energy\_gradient.cpp::compute}) and the formula
    !! collected in \texttt{doc/rdmft\_calculation.md} \S9.4.  It is
    !! \emph{not} an approximation: the ``frozen-Fock'' label that
    !! used to be attached to this code path was a misnomer -- there
    !! is no separate ``response'' term to add.
    !!
    !! Previous revisions of this routine took a longer
    !! ``\texttt{exact}'' path that built a per-band rank-one
    !! \texttt{vexx} response and combined it as
    !! ``\(\tfrac{1}{2} c_t w_{\mathbf{k}} w_t'(D_{ii} +
    !! \mathrm{response})\)''.  The response there evaluates to
    !! \(w_{\mathbf{k}} D_{ii}\) (see derivation note below), so that
    !! combination produced \(\tfrac{1}{2} c_t w_{\mathbf{k}}
    !! w_t'\,(1+w_{\mathbf{k}})\,D_{ii}\) -- a factor
    !! \((1+w_{\mathbf{k}})/2\) too large for \(w_{\mathbf{k}}\neq 1\)
    !! (e.g. \(\times 3/2\) on a \(\Gamma\)-only \(n_\mathrm{spin}=1\)
    !! cell where \(w_{\mathbf{k}} = 2\)).  The current implementation
    !! uses the correct closed-form contraction; ``rdmft\_frozen\_fock''
    !! now only controls the cached-coupling fast path in
    !! \texttt{rdmft\_grad\_n\_occ\_fast} (it has no effect on the
    !! formula itself).
    !!
    !! \texttt{bow}/\texttt{bowmod} (\(nch=0\)) follow the pair-kernel
    !! derivative in \texttt{rdmft\_xc\_compute\_exchange\_pair}.
    USE wvfct, ONLY : nbnd
    USE klist, ONLY : nks, wk
    USE rdmft_module, ONLY : rdmft_n
    !
    REAL(DP), INTENT(INOUT) :: grad_n(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: only_k, ib_only
    !
    REAL(DP), ALLOCATABLE :: vx_diag(:,:)
    REAL(DP) :: e_t
    INTEGER :: ik, ib, it, nch, klo, khi, ib_lo, ib_hi
    !
    IF (PRESENT(ib_only) .AND. .NOT. PRESENT(only_k)) &
         CALL errore('rdmft_xc_add_occ_gradient', 'ib_only requires only_k', 1)
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
    ALLOCATE(vx_diag(nbnd, nks))
    nch = rdmft_xc_n_channels()
    IF (nch == 0) THEN
       CALL rdmft_xc_compute_exchange_kopt(1, e_t, vx_diag, only_k, ib_only)
       DO ik = klo, khi
          DO ib = ib_lo, ib_hi
             CALL rdmft_xc_occ_grad_add(grad_n, vx_diag, ib, ik, wk(ik), rdmft_n(ib, ik))
          ENDDO
       ENDDO
    ELSE
       DO it = 1, nch
          CALL rdmft_xc_compute_exchange_kopt(it, e_t, vx_diag, only_k, ib_only)
          DO ik = klo, khi
             DO ib = ib_lo, ib_hi
                CALL rdmft_xc_occ_grad_add(grad_n, vx_diag, ib, ik, wk(ik), &
                     rdmft_n(ib, ik), it=it)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    DO ik = klo, khi
       DO ib = ib_lo, ib_hi
          CALL rdmft_xc_extra_grad_add(grad_n, vx_diag, ib, ik, wk(ik), rdmft_n(ib, ik))
       ENDDO
    ENDDO
    DEALLOCATE(vx_diag)
  END SUBROUTINE rdmft_xc_add_occ_gradient
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_xc_add_occ_gradient_from_diag(grad_n, diag_vx)
    !---------------------------------------------------------------
    !! Occupation gradient from cached XC diagonals (no ACE rebuild).
    USE wvfct, ONLY : nbnd
    USE klist, ONLY : nks, wk
    !
    REAL(DP), INTENT(INOUT) :: grad_n(:,:)
    REAL(DP), INTENT(IN)    :: diag_vx(:,:)
    !
    INTEGER :: ik, ib, it
    !
    it = 1
    DO ik = 1, nks
       DO ib = 1, nbnd
          CALL rdmft_xc_occ_grad_add_from_diag(grad_n, diag_vx, ib, ik, wk(ik), &
               rdmft_n(ib, ik), it=it)
       ENDDO
    ENDDO
    DO ik = 1, nks
       DO ib = 1, nbnd
          CALL rdmft_xc_extra_grad_add(grad_n, diag_vx, ib, ik, wk(ik), rdmft_n(ib, ik))
       ENDDO
    ENDDO
  END SUBROUTINE rdmft_xc_add_occ_gradient_from_diag
  !
END MODULE rdmft_xc
