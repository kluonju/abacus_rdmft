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
    FixedOne,          // always alpha0 = 1.0
    BarzilaiBorwein,   // alpha0 from BB estimate (fallback to line_search_alpha_init)
    Quadratic          // alpha0 from previous-step quadratic model (fallback when unavailable)
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

    /// ALM occupation inner: Strong Wolfe line search if occ_optimizer is LBFGS, else Armijo.
    OptimizerType occ_optimizer = OptimizerType::ConjugateGradient;
    /// Alternating orbital inner: Armijo only (all optimiser types).
    OptimizerType orb_optimizer = OptimizerType::ConjugateGradient;
    /// Single unified optimiser used by SolverStrategy::Joint. The joint
    /// strategy packs (occupation parameters, orbital coefficients) into one
    /// point on the product manifold and applies a single optimiser of this
    /// type to the packed gradient (dE/dp, Riemannian dE/dC).
    /// Joint line search: Strong Wolfe if joint_optimizer is LBFGS, else Armijo.
    OptimizerType joint_optimizer = OptimizerType::LBFGS;

    SolverStrategy strategy = SolverStrategy::Alternating;

    /// Outer alternating / joint cycles (occ then orb per cycle, or one
    /// joint product-manifold step for SolverStrategy::Joint).
    int outer_maxiter = 200;
    /// Inner iterations for optimize_orbitals (fixed occupations).
    int orb_maxiter = 50;
    /// Inner iterations for optimize_occupations (fixed orbitals).
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
    /// initial trial α₀ pre-step, accepted Armijo α post-step, line_search_alpha_init on SD fallback).
    /// Also used for ALM first-inner gradient norm, active set, joint, and other checks as in the solver.
    double occ_grad_tol = 1e-6;
    /// HF-only occupation entropy prefactor γ (binary entropy); 0 disables
    double occ_entropy_gamma = 0.0;

    double aug_lag_mu_init = 1.0;
    double aug_lag_mu_factor = 2.0;
    double aug_lag_mu_max = 1e6;
    double aug_lag_lambda_init = 0.0;

    double line_search_alpha_init = 1.0;
    double line_search_c1 = 1e-4;
    double line_search_rho = 0.5;
    int line_search_max_iter = 20;
    /// Curvature coefficient for Strong Wolfe: |g(alpha)| <= c2 * |g0|.
    /// Typical: 0.9 (lbfgs), smaller for nonlinear CG (e.g. 0.1).
    double line_search_c2 = 0.9;
    /// Zoom iteration cap in `strong_wolfe_line_search` (Nocedal & Wright zoom).
    int line_search_max_zoom = 20;
    /// If true, backtracking after a failed Armijo trial uses a quadratic
    /// model on the first failure and a cubic on later failures; if false, use
    /// geometric reduction (multiply by `line_search_rho` only). Polynomial
    /// suggestions are safeguarded inside the line search (fixed 0.1--0.5 of
    /// the last failed step; see `armijo_line_search` in `rdmft_optimizer.h`).
    bool line_search_polynomial = true;
    LineSearchInitStep occ_line_search_init_step = LineSearchInitStep::BarzilaiBorwein;

    /// ALM occupation line-search seed policy.
    /// When enabled, Armijo starts from a Barzilai-Borwein step estimate
    /// computed in occupation-parameter space (fallback to
    /// line_search_alpha_init when unavailable).
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

    /// If true, `RDMFTSolver::solve` runs finite-difference gradient checks after
    /// `precompute_cholesky_S()` and `wfc_C_to_X` (must not run earlier: orbitals are
    /// still in C-space before that).
    bool grad_check = false;
};

} // namespace rdmft

#endif // RDMFT_TYPE_H
