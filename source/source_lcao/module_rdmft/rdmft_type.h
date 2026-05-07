//==========================================================
// Author: Kai Luo
// Email: kluo@njust.edu.cn
// DATE : April 2026
//==========================================================
#ifndef RDMFT_TYPE_H
#define RDMFT_TYPE_H

#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>

namespace rdmft
{

enum class XCFunctionalType
{
    HF,
    Muller,
    Power,
    GU,
    /// Corrected Hartree-Fock (Csanyi-Arias PRB 61, 7348):
    /// f(n_i,n_j) = 1/2 n_i n_j + 1/2 sqrt(n_i(1-n_i)) sqrt(n_j(1-n_j))
    CHF,
    /// Csanyi-Goedecker-Arias (PRA 65, 032510):
    /// f(n_i,n_j) = 1/4 n_i n_j + 1/4 sqrt(n_i(2-n_i)) sqrt(n_j(2-n_j))
    CGA,
    /// GEO functional:
    /// f(n_i, n_j) = [n_i n_j + (n_i n_j)^(1/2) + 2 (n_i n_j)^(3/4)] / 4
    GEO,
    /// OptGM functional: convex combination of HF and Power(alpha)
    /// kernels with fixed literature parameters.
    OptGM
};

enum class OccParamType
{
    CosineSq,
    Logistic,
    /// Unconstrained parameterization: n_{ik} = σ(z_{ik} + λ) where λ is
    /// solved each step by bisection so that Σ_k w_k Σ_i n_{ik} = N_e.
    /// The electron-count constraint is satisfied exactly at every iteration;
    /// no augmented-Lagrangian penalty is needed in the joint strategy.
    SigmaShift
};

enum class OccInitMode
{
    /// Keep occupations from the preceding KS calculation.
    KS,
    /// Around the Fermi window, add +delta above and -delta below.
    Perturbed,
    /// Around the Fermi window, set to 1-delta (occ>=0.5) or delta (occ<0.5).
    Binary,
    /// Around the Fermi window, replace by uniform average occupation.
    Uniform
};

enum class ConstraintMethod
{
    AugmentedLagrangian,
    ProjectedGradient,
    ActiveSet
};

enum class OptimizerType
{
    SteepestDescent,
    ConjugateGradient,
    LBFGS,
    Adam
};

enum class SolverStrategy
{
    Alternating,
    /// Joint (product-manifold) optimisation: orbitals and occupations are
    /// packed into a single point on the product manifold
    ///     M = R^{Nk x Nb}  x  prod_k St(Nb, Nbasis; S^k)
    /// and stepped simultaneously. Previously called "ProductManifold".
    Joint
};

enum class BBStepMode
{
    BB1,
    BB2,
    Alternate
};

/// Initial trial-step policy for occupation Armijo line search.
enum class LineSearchInitStep
{
    Fixed,          // alpha0 = occ_ls_stepsize (INPUT rdmft_occ_ls_stepsize)
    BarzilaiBorwein,   // alpha0 from BB estimate (fallback to occ_ls_stepsize)
    Quadratic          // alpha0 from previous-step quadratic model (fallback when unavailable)
};

/// Default initial Armijo trial step and BB/quadratic fallback scale (no INPUT keyword).
constexpr double RDMFT_DEFAULT_LS_ALPHA_INIT = 1.0;

/// Line-search family for occupation/orbital inner iterations (INPUT `rdmft_occ_ls_type`,
/// `rdmft_orb_ls_type`). For **occupations**, `Auto` is Armijo. For **orbitals** (alternating
/// `optimize_orbitals` only), `Auto` uses strong Wolfe with CG, weak Wolfe with L-BFGS, and
/// Armijo with SD/Adam. Use `sw` / `wolfe` / `armijo` to force a preset regardless of optimiser.
enum class RdmftLineSearchPreset
{
    Auto,
    Armijo,
    StrongWolfe,
    WeakWolfe
};

/// Selects which occupation-number weighting is applied by occNum_func /
/// occNum_MulPsi / occNum_Mul_wfcHwfc.
///
///  For a given occupation η and XC coupling function g(η):
///   Occupation        → η
///   HalfOccupation    → 0.5 η
///   Coupling          → g(η)
///   HalfCoupling      → 0.5 g(η)
///   CouplingDerivative→ g'(η) = d g(η)/d η
///   Unity             → 1.0
enum class OccWeightMode
{
    Occupation = 0,
    HalfOccupation = 1,
    Coupling = 2,
    HalfCoupling = 3,
    CouplingDerivative = 4,
    Unity = 5
};

inline XCFunctionalType parse_xc_type(const std::string& name)
{
    if (name == "hf") return XCFunctionalType::HF;
    if (name == "muller") return XCFunctionalType::Muller;
    if (name == "power") return XCFunctionalType::Power;
    if (name == "gu") return XCFunctionalType::GU;
    if (name == "chf") return XCFunctionalType::CHF;
    if (name == "cga") return XCFunctionalType::CGA;
    if (name == "geo") return XCFunctionalType::GEO;
    if (name == "optgm") return XCFunctionalType::OptGM;
    throw std::invalid_argument("Unknown RDMFT XC functional: " + name);
}

inline std::string xc_type_to_string(XCFunctionalType type)
{
    switch (type)
    {
        case XCFunctionalType::HF: return "hf";
        case XCFunctionalType::Muller: return "muller";
        case XCFunctionalType::Power: return "power";
        case XCFunctionalType::GU: return "gu";
        case XCFunctionalType::CHF: return "chf";
        case XCFunctionalType::CGA: return "cga";
        case XCFunctionalType::GEO: return "geo";
        case XCFunctionalType::OptGM: return "optgm";
    }
    return "unknown";
}

/// How the RDMFT electron equality target N_e was built (for logging and diagnostics).
struct RDMFTNelectronTargetMeta
{
    /// Sum of KS occupation weights sum_{ik,ib} wg(ik, ib) before RDMFT.
    double sum_initial_wg = 0.;
    /// PARAM.inp.nelec used elsewhere in the LCAO driver.
    double input_nelec = 0.;
    /// If true, N_e base is input_nelec; else base is sum_initial_wg.
    bool use_input_nelec = false;
    /// Additive RDMFT-only offset: N_e = base + rdmft_nelec_delta.
    double rdmft_nelec_delta = 0.;
};

struct RDMFTConfig
{
    XCFunctionalType xc_type = XCFunctionalType::Power;
    double alpha_power = 0.656;

    OccParamType occ_param = OccParamType::CosineSq;
    ConstraintMethod constraint_method = ConstraintMethod::AugmentedLagrangian;

    /// Occupation optimiser (line search defaults to Armijo unless `rdmft_occ_ls_type` overrides).
    OptimizerType occ_optimizer = OptimizerType::ConjugateGradient;
    /// Orbital optimiser. With `rdmft_orb_ls_type` = auto: CG uses strong Wolfe, L-BFGS weak Wolfe,
    /// SD/Adam Armijo; explicit `sw` / `wolfe` / `armijo` overrides that pairing.
    OptimizerType orb_optimizer = OptimizerType::ConjugateGradient;
    /// Single unified optimiser used by SolverStrategy::Joint. The joint
    /// strategy packs (occupation parameters, orbital coefficients) into one
    /// point on the product manifold and applies a single optimiser of this
    /// type to the packed gradient (dE/dp, Riemannian dE/dC).
    /// Joint strategy uses Armijo backtracking on the packed (occ, orb) vector.
    OptimizerType joint_optimizer = OptimizerType::LBFGS;
    /// Occupation inner: `auto` / `armijo` → Armijo; `sw` / `wolfe` → Wolfe families.
    RdmftLineSearchPreset occ_ls_preset = RdmftLineSearchPreset::Auto;
    /// Orbital inner: auto pairs Wolfe type with `orb_optimizer` (CG/sw, LBFGS/ww, else Armijo).
    RdmftLineSearchPreset orb_ls_preset = RdmftLineSearchPreset::Auto;

    SolverStrategy strategy = SolverStrategy::Alternating;

    /// Outer alternating / joint cycles (occ then orb per cycle, or one
    /// joint product-manifold step for SolverStrategy::Joint).
    int outer_maxiter = 200;
    /// Inner iterations for optimize_orbitals (fixed occupations).
    int orb_maxiter = 50;
    /// Orbital SD only: one line-search trial at `orb_ls_stepsize` (INPUT
    /// `rdmft_orb_ls_stepsize`) with Armijo acceptance; no rho backtracking. Ignored
    /// for cg / lbfgs / adam.
    bool orb_ls_fixed_step = false;
    /// Orbital line search: initial trial α₀ for sd/cg/lbfgs (Adam uses 1 internally). Also the sole trial when
    /// `orb_ls_fixed_step` (INPUT `rdmft_orb_ls_fixed_step`) is true with sd.
    double orb_ls_stepsize = 1.0;
    /// Occupation SD only (ALM Armijo, PG monotone, AS monotone): one trial at
    /// `occ_ls_stepsize` with the same acceptance as the usual monotone path; no
    /// backtracking. Ignored for non-SD optimisers or Wolfe line searches.
    bool occ_ls_fixed_step = false;
    /// First trial α₀ when `occ_line_search_init_step` is Fixed (`rdmft_occ_ls_init_step` fixed); also used for
    /// `occ_ls_fixed_step` single-trial mode. Default 1.0.
    double occ_ls_stepsize = 1.0;
    /// Inner iterations for optimize_occupations (fixed orbitals). 0 = skip occupation
    /// optimization each outer cycle (alternating strategy only).
    int occ_maxiter = 50;
    /// Alternating / joint outer: when >0, require |dE| < this between outer iters
    /// **and** occupation+orbital inner `converged` flags. When <=0, outer stops on inner flags only
    /// (after first outer iter). Same field as INPUT `rdmft_energy_tol`.
    double energy_tol = 1e-8;
    /// Inner orbital step: Riemannian gradient norm ||G_R|| must fall below this
    /// for convergence (and, when orb_energy_tol > 0, the energy change criterion
    /// must also be satisfied).
    double orb_grad_tol = 1e-6;
    /// Alternating orbital inner loop: when > 0, convergence additionally requires
    /// |E_k - E_{k-1}| before the step and |E_new - E| after an accepted line search
    /// to stay below this (Ry). Set <= 0 to disable the energy criterion (gradient-only).
    double orb_energy_tol = 1e-8;
    /// Augmented Lagrangian (and similar): stop an inner step when sum_i |Δn_i| in one iteration is below this.
    double rdmft_occ_tol = 1e-8;
    /// Reserved / unused for projected_gradient (PG uses occ_grad_tol on ||g_proj|| only; kept for INPUT compat).
    double occ_energy_tol = 1e-8;
    /// PG: ||g_proj||_inf < this at post-step (g_proj = (n - P(n - τ∇E))/τ; τ from occupation line search:
    /// initial trial α₀ pre-step, accepted Armijo α post-step, or RDMFT_DEFAULT_LS_ALPHA_INIT when τ is invalid).
    /// Also used for ALM first-inner gradient norm, active set, joint, and other checks as in the solver.
    double occ_grad_tol = 1e-6;
    /// HF-only occupation entropy prefactor γ (binary entropy); 0 disables
    double occ_entropy_gamma = 0.0;

    double aug_lag_mu_init = 1.0;
    double aug_lag_mu_factor = 2.0;
    double aug_lag_mu_max = 1e6;
    double aug_lag_lambda_init = 0.0;

    /// PG with occupation optimiser CG: cap the occupation line-search initial trial
    ///   α₀ = min(seed, cap)  where seed comes from rdmft_occ_ls_init_step / BB / quad.
    /// Set <= 0 to disable capping (previous behaviour).
    double pg_occ_cg_ls_alpha_cap = 1e-3;
    /// PG: after failed non-monotone Wolfe (occ CG), try SD from this initial α (monotone
    /// projected backtracking). Set <= 0 to skip this recovery.
    double pg_occ_ls_recovery_alpha = 1e-8;
    double line_search_c1 = 1e-4;
    double line_search_rho = 0.5;
    int line_search_max_iter = 20;
    /// Curvature coefficient for Wolfe conditions:
    ///   strong Wolfe (CG): |g(alpha)| <= c2 * |g0|
    ///   weak Wolfe (LBFGS): g(alpha) >= c2 * g0
    /// Typical: 0.9 (lbfgs), smaller for nonlinear CG (e.g. 0.1).
    double line_search_c2 = 0.9;
    /// Zhang–Hager memory η for non-monotone Strong Wolfe (PG/AS occupation + CG
    /// in alternating mode only; ALM keeps monotone Wolfe). Must be in (0, 1].
    double occ_cg_nonmonotone_eta = 0.85;
    /// Zoom iteration cap in `strong_wolfe_line_search` (Nocedal & Wright zoom).
    int line_search_max_zoom = 20;
    /// If true, backtracking after a failed Armijo trial uses a quadratic
    /// model on the first failure and a cubic on later failures; if false, use
    /// geometric reduction (multiply by `line_search_rho` only). Polynomial
    /// suggestions are safeguarded inside the line search (fixed 0.1--0.5 of
    /// the last failed step; see `armijo_line_search` in `rdmft_optimizer.h`).
    bool line_search_polynomial = true;
    LineSearchInitStep occ_line_search_init_step = LineSearchInitStep::BarzilaiBorwein;

    /// Barzilai–Borwein line-search seed (ALM params, joint packed SD/CG).
    /// When enabled, the first trial α uses BB (spectral ratio ||x||/||g|| when
    /// no prior step exists, else BB1/BB2 per alm_bb_mode), clamped to
    /// [alm_bb_alpha_min, alm_bb_alpha_max], with fallback to
    /// RDMFT_DEFAULT_LS_ALPHA_INIT when the estimate is unusable.
    bool alm_bb_enabled = true;
    BBStepMode alm_bb_mode = BBStepMode::Alternate;
    double alm_bb_alpha_min = 1e-8;
    double alm_bb_alpha_max = 10.0;

    /// Relative scaling between the orbital and occupation parameter blocks
    /// in the joint (product-manifold) strategy. The packed optimisation
    /// variable is (p, C) where p are the unconstrained occupation
    /// parameters (in radians for cosine^2, dimensionless for logistic)
    /// and C are the orbital coefficients. The natural scale of dE/dp
    /// depends on the Jacobian dn/dp (which can range from 0 to 1 across
    /// the Brillouin zone and band index), while dE/dC scales with the
    /// Hamiltonian matrix elements. In a single-alpha Armijo line search
    /// this block-scale mismatch manifests as either (a) well-behaved
    /// occupation steps together with far too aggressive orbital steps,
    /// or (b) vice versa. Multiplying the orbital gradient (as fed to
    /// the Euclidean optimiser) and the orbital direction (as applied by
    /// retraction) by joint_orb_scale is equivalent to changing the
    /// units of the orbital block inside the packed vector:
    ///   C_internal = C / joint_orb_scale
    ///   dE/dC_internal = joint_orb_scale * dE/dC
    ///   alpha * dir_internal = alpha * joint_orb_scale * (-dE/dC)
    /// so the *physical* orbital step, alpha * joint_orb_scale * dE/dC,
    /// is attenuated or amplified by joint_orb_scale^2 (which makes this
    /// knob behave like a per-block pre-conditioner).
    ///
    /// Default = 1.0 (no rescaling). For systems where the initial joint
    /// step would move the orbitals too aggressively a value < 1 damps
    /// the orbital block; for systems where the orbital sub-problem is
    /// hard a value > 1 gives it more weight.
    double joint_orb_scale = 1.0;

    double adam_lr = 0.001;
    double adam_beta1 = 0.9;
    double adam_beta2 = 0.999;
    double adam_eps = 1e-8;

    int lbfgs_memory = 10;

    bool use_roptlite = false;

    double fd_epsilon = 1e-5;

    /// Optional additive perturbation magnitude applied on top of the
    /// selected initial-occupation mode before the first optimisation step.
    OccInitMode occ_init_mode = OccInitMode::KS;
    /// Optional additive perturbation magnitude (delta) used by
    /// OccInitMode::Perturbed.
    double occ_init_perturb = 0.0;
    /// Number of bands in the Fermi window per side (above/below Fermi),
    /// used by OccInitMode::Perturbed and OccInitMode::Uniform.
    int occ_init_nbands_top = 0;

    /// Log per-k Stiefel Gram residual (alternating outer loop); expensive, default off.
    bool print_stiefel_gram = false;

    /// After solve: print ELK rdmeval-style ε_ik (probe n=0.5, O(nk*nbands) compute calls).
    bool print_evals = false;

    /// Each `EnergyGradient::compute`: log one-body / Hartree / EXX orbital gradient norms
    /// (ambient Frobenius in C- and X-space, canonical after per-term projection).
    bool print_orb_grad_decomp = false;

    /// If true, `RDMFTSolver::solve` runs finite-difference gradient checks after
    /// `precompute_cholesky_S()` and `wfc_C_to_X` (must not run earlier: orbitals are
    /// still in C-space before that).
    bool grad_check = false;
};

} // namespace rdmft

#endif // RDMFT_TYPE_H
