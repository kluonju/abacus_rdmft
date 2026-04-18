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
    ProductManifold
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

    SolverStrategy strategy = SolverStrategy::Alternating;

    int orb_maxiter = 200;
    int occ_maxiter = 50;
    double energy_tol = 1e-8;
    double grad_tol = 1e-6;

    double aug_lag_mu_init = 1.0;
    double aug_lag_mu_factor = 2.0;
    double aug_lag_mu_max = 1e6;

    double line_search_alpha_init = 0.1;
    double line_search_c1 = 1e-4;
    double line_search_rho = 0.5;
    int line_search_max_iter = 30;

    double adam_lr = 0.001;
    double adam_beta1 = 0.9;
    double adam_beta2 = 0.999;
    double adam_eps = 1e-8;

    int lbfgs_memory = 10;

    bool use_roptlite = false;

    double fd_epsilon = 1e-5;

    /// Distance inside [0,1] that the initial occupations are clamped away
    /// from the boundary before being passed to the optimiser. For the
    /// cosine^2 / logistic parameterisations the Jacobian dn/dp vanishes at
    /// 0 and 1, which stalls the parameter-space optimiser on a KS seed
    /// with integer occupations. Setting occ_init_margin > 0 maps
    ///     n -> clamp(n, m, 1-m)
    /// (keeping sum preserved by normalising afterwards) so the parameter
    /// gradient is non-zero at step 1.
    /// Use 0.0 to disable; default 1e-3 is small enough to leave the
    /// physical answer unchanged once the optimiser converges.
    double occ_init_margin = 1e-3;
};

} // namespace rdmft

#endif // RDMFT_TYPE_H
