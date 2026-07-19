#ifndef SOURCE_RDMFT_RDMFT_PARAMS_H
#define SOURCE_RDMFT_RDMFT_PARAMS_H

#include <string>
#include <algorithm>
#include <cctype>

//! \file rdmft_params.h
//! \brief Basis-independent configuration for the redesigned RDMFT module.
//!
//! This module implements Reduced Density Matrix Functional Theory (RDMFT)
//! as a standalone, basis-agnostic optimisation core.  It is intentionally
//! decoupled from the plane-wave / LCAO Hamiltonian machinery: everything
//! the optimisers need is expressed through the abstract
//! rdmft::RdmftBackend oracle (see rdmft_backend.h).  The optimisation
//! algorithms (SPG2 occupation block, EBI occupation block and the Stiefel
//! orbital block) are ported from the Quantum ESPRESSO RDMFT reference so
//! that the same natural-orbital functionals can be driven from either the
//! plane-wave or the LCAO ESolver_RDMFT backends.

namespace rdmft
{

//! Natural-orbital exchange-correlation functional.
enum class XcType
{
    HF,      //!< Hartree-Fock (g(n)=n).
    Muller,  //!< Mueller / BB (alpha = 1/2).
    Power,   //!< Power functional (user alpha, default 0.656).
    GU,      //!< Goedecker-Umrigar (alpha = 1/2 + on-site correction).
    CHF,     //!< Corrected Hartree-Fock (2 channels).
    CGA,     //!< Corrected Gilbert-Aryasetiawan (2 channels).
    GEO,     //!< Geometric mean functional (3 channels).
    HybOpt,  //!< Optimised power/HF hybrid (2 channels).
    BOWMOD   //!< Modified Baldsiefen-Ohlrogge-Weimer functional (4 channels).
};

//! Occupation-number constraint / parameterisation strategy.
enum class OccOptimizerType
{
    SPG2,  //!< Birgin-Martinez-Raydan spectral projected gradient (default).
    EBI    //!< Explicit-by-implicit erf parameterisation (Yao et al. 2022).
};

//! Orbital (Stiefel) descent method.
enum class OrbOptimizerType
{
    SD,    //!< Steepest descent.
    CG,    //!< Polak-Ribiere conjugate gradient (default).
    LBFGS  //!< Limited-memory BFGS on the packed tangent space.
};

//! Top-level solver strategy.
enum class SolverStrategy
{
    Alternating, //!< Alternate occupation and orbital inner blocks (default).
    OccOnly,     //!< Only optimise occupations (frozen orbitals).
    OrbOnly      //!< Only optimise orbitals (frozen occupations).
};

//! Line-search flavour used by the inner blocks.
enum class LineSearchType
{
    Armijo,      //!< Monotone Armijo backtracking (default for occupations).
    StrongWolfe  //!< Strong-Wolfe zoom (default for orbitals).
};

//! Complete configuration for a RDMFT run.  All defaults mirror the
//! Quantum ESPRESSO RDMFT reference (see qe-rdmft/rdmft_module.f90).
struct RdmftParams
{
    // Functional -----------------------------------------------------------
    XcType xc = XcType::Muller;   //!< Natural-orbital functional.
    double power_alpha = 0.656;   //!< Power exponent (overwritten by fixed functionals).
    double reg_eps = 1.0e-8;      //!< Small-occupation regularisation.

    // Strategy -------------------------------------------------------------
    SolverStrategy strategy = SolverStrategy::Alternating;
    OccOptimizerType occ_optimizer = OccOptimizerType::SPG2;
    OrbOptimizerType orb_optimizer = OrbOptimizerType::CG;

    // Iteration limits -----------------------------------------------------
    int outer_maxiter = 50;   //!< Outer alternating cycles.
    int occ_maxiter = 20;     //!< Inner occupation iterations per block.
    int orb_maxiter = 20;     //!< Inner orbital iterations per block.
    int lbfgs_memory = 10;    //!< L-BFGS history depth.

    // Tolerances -----------------------------------------------------------
    double energy_tol = 1.0e-7;    //!< Outer |dE| convergence.
    double occ_grad_tol = 1.0e-5;  //!< SPG ||g1||_inf / EBI ||grad_x||.
    double occ_tol = 1.0e-6;       //!< EBI sum|dn| convergence.
    double orb_grad_tol = 1.0e-5;  //!< ||G_R|| Riemannian gradient.

    // Line search ----------------------------------------------------------
    LineSearchType occ_ls = LineSearchType::Armijo;
    LineSearchType orb_ls = LineSearchType::StrongWolfe;
    double ls_c1 = 1.0e-4;   //!< Armijo / Wolfe sufficient-decrease constant.
    double ls_c2 = 0.9;      //!< Wolfe curvature constant.
    double ls_rho = 0.5;     //!< Backtracking contraction factor.
    int ls_max_iter = 40;    //!< Max backtracking / bracket iterations.
    int ls_max_zoom = 20;    //!< Max Wolfe zoom iterations.

    // Barzilai-Borwein spectral steplength (SPG / initial guesses) ---------
    double bb_alpha_min = 1.0e-8;
    double bb_alpha_max = 1.0e2;

    // Electron-count constraint -------------------------------------------
    double nelec = 0.0;             //!< Total electrons (weighted-sum target).
    bool fix_magnetization = false; //!< Constrain N_up / N_down separately.
    double nelec_up = 0.0;          //!< Spin-up target (fix_magnetization).
    double nelec_down = 0.0;        //!< Spin-down target (fix_magnetization).

    int verbosity = 1;              //!< 0 silent, 1 outer/inner, 2 detailed.
};

//! Case-insensitive functional-name parsing.  Returns true on success and
//! also fixes power_alpha for functionals that pin the exponent.
inline bool parse_xc_type(const std::string& name_in, XcType& out, double& power_alpha)
{
    std::string s = name_in;
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    if (s == "hf")
    {
        out = XcType::HF;
        power_alpha = 1.0;
    }
    else if (s == "muller" || s == "mueller" || s == "bb")
    {
        out = XcType::Muller;
        power_alpha = 0.5;
    }
    else if (s == "power")
    {
        out = XcType::Power;
        // keep the user-provided alpha
    }
    else if (s == "gu")
    {
        out = XcType::GU;
        power_alpha = 0.5;
    }
    else if (s == "chf")
    {
        out = XcType::CHF;
        power_alpha = 1.0;
    }
    else if (s == "cga")
    {
        out = XcType::CGA;
        power_alpha = 1.0;
    }
    else if (s == "geo")
    {
        out = XcType::GEO;
        power_alpha = 0.75;
    }
    else if (s == "hybopt")
    {
        out = XcType::HybOpt;
        power_alpha = 0.541076;
    }
    else if (s == "bowmod")
    {
        out = XcType::BOWMOD;
        power_alpha = 0.61;
    }
    else
    {
        return false;
    }
    return true;
}

inline bool parse_occ_optimizer(const std::string& name_in, OccOptimizerType& out)
{
    std::string s = name_in;
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    if (s == "spg2" || s == "spg" || s == "sd" || s.empty())
    {
        out = OccOptimizerType::SPG2;
    }
    else if (s == "ebi")
    {
        out = OccOptimizerType::EBI;
    }
    else
    {
        return false;
    }
    return true;
}

inline bool parse_orb_optimizer(const std::string& name_in, OrbOptimizerType& out)
{
    std::string s = name_in;
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    if (s == "sd")
    {
        out = OrbOptimizerType::SD;
    }
    else if (s == "cg" || s.empty())
    {
        out = OrbOptimizerType::CG;
    }
    else if (s == "lbfgs" || s == "l-bfgs" || s == "bfgs")
    {
        out = OrbOptimizerType::LBFGS;
    }
    else
    {
        return false;
    }
    return true;
}

} // namespace rdmft

#endif // SOURCE_RDMFT_RDMFT_PARAMS_H
