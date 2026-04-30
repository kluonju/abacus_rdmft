#ifndef RDMFT_XC_FUNCTIONAL_H
#define RDMFT_XC_FUNCTIONAL_H

#include "rdmft_type.h"
#include <cmath>
#include <cassert>
#include <algorithm>

namespace rdmft
{

/// Coupling function g(n) and its derivative for RDMFT exchange-correlation functionals.
/// For separable functionals: f(n_i, n_j) = g(n_i) * g(n_j).
///
/// Regularization near n = 0 (Muller/Power/GU)
/// -------------------------------------------
/// For alpha < 1 the analytic derivative
///     g'(n) = alpha * n^(alpha-1)
/// diverges as n -> 0+. In practice this blows up the occupation gradient
/// for unoccupied bands and destroys the Armijo line-search / CG history
/// in the occupation loop.
///
/// We regularise by smoothly replacing the singular power with its Taylor
/// expansion around a small cutoff eps:
///     g(n)  = n^alpha                                     for n >= eps
///           = eps^alpha + alpha * eps^(alpha-1) * (n-eps)  for n < eps
/// so that g(0) = eps^alpha - alpha * eps^alpha = (1-alpha) * eps^alpha,
/// g'(0) = alpha * eps^(alpha-1) (bounded), and g,g' are continuous at eps.
/// The default cutoff eps = 1e-8 is many orders of magnitude below the
/// target tolerances, so it does not affect the final converged energy
/// for the ZnO-size problems we target, while keeping derivatives bounded
/// away from the 0/1 boundaries.
class XCFunctional
{
  public:
    XCFunctional(XCFunctionalType type,
                 double alpha = 0.656,
                 double reg_eps = 1e-8)
        : type_(type), alpha_(alpha), reg_eps_(reg_eps)
    {
        if (type == XCFunctionalType::HF) alpha_ = 1.0;
        else if (type == XCFunctionalType::Muller) alpha_ = 0.5;
        else if (type == XCFunctionalType::GU) alpha_ = 0.5; // Müller sqrt(n) for RI DM; GU SIC added in engine
        else if (type == XCFunctionalType::BBC3) alpha_ = 0.5; // Müller-like power for fallback g
        // GEO: alpha_ unused for the energy / DM construction; the GEO functional decomposes
        // as a weighted sum of three separable terms with exponents (1, 1/2, 3/4). The
        // single-g entry points (g/dg/d2g) are still provided using the dominant 3/4 power
        // for callers that probe a single coupling (e.g. is_separable() == false handlers
        // ignore them).
        else if (type == XCFunctionalType::GEO) alpha_ = 0.75;
    }

    XCFunctionalType type() const { return type_; }
    double alpha() const { return alpha_; }
    double reg_eps() const { return reg_eps_; }

    /// g(n): coupling function for the separable part (regularised near 0)
    double g(double n) const
    {
        n = std::max(0.0, std::min(1.0, n));
        if (alpha_ >= 1.0 - 1e-12) return n; // HF: no regularisation needed
        if (n >= reg_eps_) return std::pow(n, alpha_);
        // Linear extrapolation from eps down to 0
        const double g_eps  = std::pow(reg_eps_, alpha_);
        const double dg_eps = alpha_ * std::pow(reg_eps_, alpha_ - 1.0);
        return g_eps + dg_eps * (n - reg_eps_);
    }

    /// g'(n): derivative of g(n) w.r.t. n (regularised to a bounded value on [0,eps])
    double dg(double n) const
    {
        n = std::max(0.0, std::min(1.0, n));
        if (alpha_ >= 1.0 - 1e-12) return 1.0; // HF
        const double nmin = std::max(n, reg_eps_);
        return alpha_ * std::pow(nmin, alpha_ - 1.0);
    }

    /// g''(n): second derivative (0 on the regularised interval)
    double d2g(double n) const
    {
        n = std::max(0.0, std::min(1.0, n));
        if (n < reg_eps_) return 0.0;
        return alpha_ * (alpha_ - 1.0) * std::pow(n, alpha_ - 2.0);
    }

    /// Whether the functional is separable: f(n_i, n_j) = g(n_i)*g(n_j) for all i,j
    bool is_separable() const
    {
        return type_ != XCFunctionalType::GU && type_ != XCFunctionalType::BBC3
               && type_ != XCFunctionalType::GEO;
    }

    /// GEO functional decomposition.
    ///   f^GEO(n_p, n_q) = (1/4) g_1(n_p) g_1(n_q)
    ///                   + (1/4) g_{1/2}(n_p) g_{1/2}(n_q)
    ///                   + (1/2) g_{3/4}(n_p) g_{3/4}(n_q)
    /// where g_alpha(n) = n^alpha (regularised near 0 the same way as g()).
    /// Returns the number of GEO separable terms (3); coefficients and exponents are
    /// indexed in [0, num_geo_terms()).
    static constexpr int num_geo_terms() { return 3; }
    static double geo_coef(int idx)
    {
        // 0.25 * (n_p n_q) + 0.25 * sqrt(n_p n_q) + 0.5 * (n_p n_q)^{3/4}
        switch (idx)
        {
            case 0: return 0.25; // n^1
            case 1: return 0.25; // n^{1/2}
            case 2: return 0.50; // n^{3/4}
            default: return 0.0;
        }
    }
    static double geo_alpha(int idx)
    {
        switch (idx)
        {
            case 0: return 1.0;
            case 1: return 0.5;
            case 2: return 0.75;
            default: return 1.0;
        }
    }
    /// Eps-regularised n^alpha with the same C^1 cutoff used for g/dg.
    double pow_reg(double n, double alpha) const
    {
        n = std::max(0.0, std::min(1.0, n));
        if (alpha >= 1.0 - 1e-12) return n;
        if (n >= reg_eps_) return std::pow(n, alpha);
        const double g_eps = std::pow(reg_eps_, alpha);
        const double dg_eps = alpha * std::pow(reg_eps_, alpha - 1.0);
        return g_eps + dg_eps * (n - reg_eps_);
    }
    /// Eps-regularised d(n^alpha)/dn (bounded near 0).
    double dpow_reg(double n, double alpha) const
    {
        n = std::max(0.0, std::min(1.0, n));
        if (alpha >= 1.0 - 1e-12) return 1.0;
        const double nmin = std::max(n, reg_eps_);
        return alpha * std::pow(nmin, alpha - 1.0);
    }

    /// Full coupling f(n_i, n_j) for the GU functional (non-separable).
    /// For separable functionals, this returns g(n_i)*g(n_j).
    /// The GU off-diagonal piece reuses the same eps-regularised sqrt(n)
    /// used for Muller so we don't hit a divergent derivative at n_i=0.
    double f(double ni, double nj, bool same_orbital = false) const
    {
        if (type_ == XCFunctionalType::GU)
        {
            const double nic = std::max(0.0, std::min(1.0, ni));
            const double njc = std::max(0.0, std::min(1.0, nj));
            if (same_orbital) return nic * nic;
            // eps-regularised sqrt: matches the Muller branch of g()
            auto sqrt_reg = [&](double x) {
                if (x >= reg_eps_) return std::sqrt(x);
                const double g_eps  = std::sqrt(reg_eps_);
                const double dg_eps = 0.5 / std::sqrt(reg_eps_);
                return g_eps + dg_eps * (x - reg_eps_);
            };
            return sqrt_reg(nic) * sqrt_reg(njc);
        }
        if (type_ == XCFunctionalType::GEO)
        {
            double sum = 0.0;
            for (int t = 0; t < num_geo_terms(); ++t)
            {
                sum += geo_coef(t) * pow_reg(ni, geo_alpha(t)) * pow_reg(nj, geo_alpha(t));
            }
            return sum;
        }
        return g(ni) * g(nj);
    }

    /// df/dn_i for the GU functional (eps-regularised)
    double df_dni(double ni, double nj, bool same_orbital = false) const
    {
        if (type_ == XCFunctionalType::GU)
        {
            const double nic = std::max(0.0, std::min(1.0, ni));
            const double njc = std::max(0.0, std::min(1.0, nj));
            if (same_orbital) return 2.0 * nic;
            const double ni_safe = std::max(nic, reg_eps_);
            const double nj_safe = std::max(njc, reg_eps_);
            return 0.5 / std::sqrt(ni_safe) * std::sqrt(nj_safe);
        }
        if (type_ == XCFunctionalType::GEO)
        {
            double sum = 0.0;
            for (int t = 0; t < num_geo_terms(); ++t)
            {
                const double a = geo_alpha(t);
                sum += geo_coef(t) * dpow_reg(ni, a) * pow_reg(nj, a);
            }
            return sum;
        }
        return dg(ni) * g(nj);
    }

    /// For the GU functional, the diagonal correction energy per orbital.
    /// E_GU = E_Muller + 0.5 * sum_{ik} w_k^2 * (n_ik^2 - n_ik) * J_{ii}
    /// Returns (n^2 - n) which is the diagonal factor difference.
    double gu_diag_factor(double n) const
    {
        assert(type_ == XCFunctionalType::GU);
        n = std::max(0.0, std::min(1.0, n));
        return n * n - n;
    }

    /// Derivative of the GU diagonal factor w.r.t. n: 2n - 1
    double gu_diag_factor_deriv(double n) const
    {
        assert(type_ == XCFunctionalType::GU);
        n = std::max(0.0, std::min(1.0, n));
        return 2.0 * n - 1.0;
    }

  private:
    XCFunctionalType type_;
    double alpha_;
    double reg_eps_;
};

} // namespace rdmft

#endif // RDMFT_XC_FUNCTIONAL_H
