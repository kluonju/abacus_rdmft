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
/// We regularise by replacing the singular branch on [0, eps] with a C^1 cubic
/// Hermite polynomial in x = n/eps:
///     p(x) = (alpha-2)x^3 + (3-alpha)x^2
/// and
///     g(n)  = n^alpha                    for n >= eps
///           = eps^alpha * p(n/eps)       for n < eps
/// This enforces g(0)=0, g'(0)=0, g(eps)=eps^alpha, g'(eps)=alpha*eps^(alpha-1).
/// So g and g' remain continuous at eps while empty bands keep zero coupling.
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
        else if (type == XCFunctionalType::OptGM) alpha_ = 0.75;
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
        return regularized_power_value(n, alpha_);
    }

    /// g'(n): derivative of g(n) w.r.t. n (regularised to a bounded value on [0,eps])
    double dg(double n) const
    {
        n = std::max(0.0, std::min(1.0, n));
        if (alpha_ >= 1.0 - 1e-12) return 1.0; // HF
        if (n >= reg_eps_) return alpha_ * std::pow(n, alpha_ - 1.0);
        return regularized_power_derivative(n, alpha_);
    }

    /// g''(n): second derivative (finite on the regularised interval)
    double d2g(double n) const
    {
        n = std::max(0.0, std::min(1.0, n));
        if (alpha_ >= 1.0 - 1e-12) return 0.0; // HF branch
        if (n < reg_eps_) return regularized_power_second_derivative(n, alpha_);
        return alpha_ * (alpha_ - 1.0) * std::pow(n, alpha_ - 2.0);
    }

    /// Whether the functional is separable: f(n_i, n_j) = g(n_i)*g(n_j) for all i,j
    bool is_separable() const
    {
        return type_ != XCFunctionalType::GU && type_ != XCFunctionalType::BBC3
               && type_ != XCFunctionalType::GEO && type_ != XCFunctionalType::OptGM;
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
    /// optGM mixture weights for the same (alpha_0, alpha_1, alpha_2) as GEO.
    static double optgm_coef(int idx)
    {
        switch (idx)
        {
            case 0: return 0.00675; // n^1
            case 1: return 0.64213; // n^{1/2}
            case 2: return 0.35112; // n^{3/4}
            default: return 0.0;
        }
    }
    /// Coefficient c_t for GEO or optGM mixture term `idx`.
    double mixture_coef(int idx) const
    {
        if (type_ == XCFunctionalType::GEO) return geo_coef(idx);
        if (type_ == XCFunctionalType::OptGM) return optgm_coef(idx);
        assert(false && "mixture_coef only for GEO / optGM");
        return 0.0;
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
        return regularized_power_value(n, alpha);
    }
    /// Eps-regularised d(n^alpha)/dn (bounded near 0).
    double dpow_reg(double n, double alpha) const
    {
        n = std::max(0.0, std::min(1.0, n));
        if (alpha >= 1.0 - 1e-12) return 1.0;
        if (n >= reg_eps_) return alpha * std::pow(n, alpha - 1.0);
        return regularized_power_derivative(n, alpha);
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
            // eps-regularised sqrt branch consistent with pow_reg(alpha=0.5).
            return pow_reg(nic, 0.5) * pow_reg(njc, 0.5);
        }
        if (type_ == XCFunctionalType::GEO || type_ == XCFunctionalType::OptGM)
        {
            double sum = 0.0;
            for (int t = 0; t < num_geo_terms(); ++t)
            {
                sum += mixture_coef(t) * pow_reg(ni, geo_alpha(t)) * pow_reg(nj, geo_alpha(t));
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
            return dpow_reg(nic, 0.5) * pow_reg(njc, 0.5);
        }
        if (type_ == XCFunctionalType::GEO || type_ == XCFunctionalType::OptGM)
        {
            double sum = 0.0;
            for (int t = 0; t < num_geo_terms(); ++t)
            {
                const double a = geo_alpha(t);
                sum += mixture_coef(t) * dpow_reg(ni, a) * pow_reg(nj, a);
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
    double regularized_power_value(double n, double alpha) const
    {
        // C1 cubic Hermite on [0, reg_eps] matching value and slope at reg_eps,
        // while enforcing g(0)=0 and g'(0)=0 to keep empty states inactive.
        const double x = n / reg_eps_;
        const double p = (alpha - 2.0) * x * x * x + (3.0 - alpha) * x * x;
        return std::pow(reg_eps_, alpha) * p;
    }

    double regularized_power_derivative(double n, double alpha) const
    {
        const double x = n / reg_eps_;
        const double dpdx = 3.0 * (alpha - 2.0) * x * x + 2.0 * (3.0 - alpha) * x;
        return std::pow(reg_eps_, alpha - 1.0) * dpdx;
    }

    double regularized_power_second_derivative(double n, double alpha) const
    {
        const double x = n / reg_eps_;
        const double d2pdx2 = 6.0 * (alpha - 2.0) * x + 2.0 * (3.0 - alpha);
        return std::pow(reg_eps_, alpha - 2.0) * d2pdx2;
    }

    XCFunctionalType type_;
    double alpha_;
    double reg_eps_;
};

} // namespace rdmft

#endif // RDMFT_XC_FUNCTIONAL_H
