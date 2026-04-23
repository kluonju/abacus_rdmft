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
    GU
};

enum class OccParamType
{
    CosineSq,
    Logistic
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

inline XCFunctionalType parse_xc_type(const std::string& name)
{
    if (name == "hf") return XCFunctionalType::HF;
    if (name == "muller") return XCFunctionalType::Muller;
    if (name == "power") return XCFunctionalType::Power;
    if (name == "gu") return XCFunctionalType::GU;
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
    }
    return "unknown";
}

struct RDMFTConfig
{
    XCFunctionalType xc_type = XCFunctionalType::Power;
    double alpha_power = 0.656;

    OccParamType occ_param = OccParamType::CosineSq;
    ConstraintMethod constraint_method = ConstraintMethod::AugmentedLagrangian;

    OptimizerType occ_optimizer = OptimizerType::ConjugateGradient;
    OptimizerType orb_optimizer = OptimizerType::ConjugateGradient;
    /// Single unified optimiser used by SolverStrategy::Joint. The joint
    /// strategy packs (occupation parameters, orbital coefficients) into one
    /// point on the product manifold and applies a single optimiser of this
    /// type to the packed gradient (dE/dp, Riemannian dE/dC).
    OptimizerType joint_optimizer = OptimizerType::LBFGS;

    SolverStrategy strategy = SolverStrategy::Alternating;

    /// Outer alternating / joint cycles (occ then orb per cycle, or one
    /// joint product-manifold step for SolverStrategy::Joint).
    int outer_maxiter = 200;
    /// Inner iterations for optimize_orbitals (fixed occupations).
    int orb_maxiter = 50;
    /// Inner iterations for optimize_occupations (fixed orbitals).
    int occ_maxiter = 50;
    double energy_tol = 1e-8;
    /// Inner orbital step: stop when Riemannian gradient norm ||G_R|| falls below this.
    double orb_grad_tol = 1e-6;
    /// Inner occupation step: stop when sum_i |Δn_i| in one iteration falls below this
    /// as the occupation convergence criterion.
    double rdmft_occ_tol = 1e-8;

    double aug_lag_mu_init = 1.0;
    double aug_lag_mu_factor = 2.0;
    double aug_lag_mu_max = 1e6;
    double aug_lag_lambda_init = 0.0;

    double line_search_alpha_init = 0.1;
    double line_search_c1 = 1e-4;
    double line_search_rho = 0.5;
    int line_search_max_iter = 30;

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
};

} // namespace rdmft

#endif // RDMFT_TYPE_H
