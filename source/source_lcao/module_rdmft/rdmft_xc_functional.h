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
class XCFunctional
{
  public:
    XCFunctional(XCFunctionalType type, double alpha = 0.656)
        : type_(type), alpha_(alpha)
    {
        if (type == XCFunctionalType::HF) alpha_ = 1.0;
        else if (type == XCFunctionalType::Muller) alpha_ = 0.5;
    }

    XCFunctionalType type() const { return type_; }
    double alpha() const { return alpha_; }

    /// g(n): coupling function for the separable part
    double g(double n) const
    {
        n = std::max(0.0, std::min(1.0, n));
        return std::pow(n, alpha_);
    }

    /// g'(n): derivative of g(n) w.r.t. n
    double dg(double n) const
    {
        n = std::max(1e-30, std::min(1.0, n));
        return alpha_ * std::pow(n, alpha_ - 1.0);
    }

    /// g''(n): second derivative (needed for some optimizers)
    double d2g(double n) const
    {
        n = std::max(1e-30, std::min(1.0, n));
        return alpha_ * (alpha_ - 1.0) * std::pow(n, alpha_ - 2.0);
    }

    /// Whether the functional is separable: f(n_i, n_j) = g(n_i)*g(n_j) for all i,j
    bool is_separable() const
    {
        return type_ != XCFunctionalType::GU;
    }

    /// Full coupling f(n_i, n_j) for the GU functional (non-separable)
    /// For separable functionals, this returns g(n_i)*g(n_j).
    double f(double ni, double nj, bool same_orbital = false) const
    {
        if (type_ == XCFunctionalType::GU)
        {
            if (same_orbital)
                return std::max(0.0, std::min(1.0, ni)) * std::max(0.0, std::min(1.0, ni));
            else
                return std::sqrt(std::max(0.0, ni)) * std::sqrt(std::max(0.0, nj));
        }
        return g(ni) * g(nj);
    }

    /// df/dn_i for the GU functional
    double df_dni(double ni, double nj, bool same_orbital = false) const
    {
        if (type_ == XCFunctionalType::GU)
        {
            ni = std::max(1e-30, std::min(1.0, ni));
            nj = std::max(1e-30, std::min(1.0, nj));
            if (same_orbital)
                return 2.0 * ni;
            else
                return 0.5 / std::sqrt(ni) * std::sqrt(nj);
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
};

} // namespace rdmft

#endif // RDMFT_XC_FUNCTIONAL_H
