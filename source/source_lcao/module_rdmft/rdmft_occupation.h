//==========================================================
// Author: Kai Luo
// Email: kluo@njust.edu.cn
// DATE : April 2026
//==========================================================
#ifndef RDMFT_OCCUPATION_H
#define RDMFT_OCCUPATION_H

#include "rdmft_type.h"
#include "source_base/global_variable.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <numeric>

namespace rdmft
{

class OccupationConstraint;

/// Parameterization of occupation numbers to enforce box constraint n in [0,1].
/// Provides mapping from unconstrained parameter to n, and chain-rule Jacobian.
class OccupationParam
{
  public:
    explicit OccupationParam(OccParamType type = OccParamType::CosineSq)
        : type_(type) {}

    OccParamType type() const { return type_; }

    /// Map unconstrained parameter p -> occupation n in [0,1]
    double to_occ(double p) const
    {
        if (type_ == OccParamType::CosineSq)
        {
            double c = std::cos(p);
            return c * c;
        }
        else // Logistic or SigmaShift (SigmaShift uses logistic formula; shift handled separately)
        {
            return 1.0 / (1.0 + std::exp(-p));
        }
    }

    /// Map occupation n in [0,1] -> unconstrained parameter p
    double to_param(double n) const
    {
        n = std::max(1e-12, std::min(1.0 - 1e-12, n));
        if (type_ == OccParamType::CosineSq)
        {
            return std::acos(std::sqrt(n));
        }
        else // Logistic or SigmaShift
        {
            return std::log(n / (1.0 - n));
        }
    }

    /// Jacobian dn/dp
    double jacobian(double p) const
    {
        if (type_ == OccParamType::CosineSq)
        {
            return -std::sin(2.0 * p);
        }
        else // Logistic or SigmaShift
        {
            double n = to_occ(p);
            return n * (1.0 - n);
        }
    }

    /// Transform gradient from dE/dn to dE/dp = (dE/dn) * (dn/dp)
    double transform_gradient(double dE_dn, double p) const
    {
        return dE_dn * jacobian(p);
    }

    /// Batch operations on vectors
    void params_to_occ(const std::vector<double>& params, std::vector<double>& occ) const
    {
        occ.resize(params.size());
        for (size_t i = 0; i < params.size(); ++i)
            occ[i] = to_occ(params[i]);
    }

    void occ_to_params(const std::vector<double>& occ, std::vector<double>& params) const
    {
        params.resize(occ.size());
        for (size_t i = 0; i < occ.size(); ++i)
            params[i] = to_param(occ[i]);
    }

    void transform_gradient_batch(const std::vector<double>& dE_dn,
                                  const std::vector<double>& params,
                                  std::vector<double>& dE_dp) const
    {
        dE_dp.resize(dE_dn.size());
        for (size_t i = 0; i < dE_dn.size(); ++i)
            dE_dp[i] = transform_gradient(dE_dn[i], params[i]);
    }

  private:
    OccParamType type_;
};


/// Handles the electron number constraint: sum_k w_k sum_i n_ik = N_e.
/// Supports augmented Lagrangian, projected gradient, and active set methods.
class OccupationConstraint
{
  public:
    OccupationConstraint(ConstraintMethod method, double n_electrons,
                         const std::vector<double>& kweights,
                         int nbands)
        : method_(method), n_electrons_(n_electrons),
          kweights_(kweights), nbands_(nbands),
          lambda_(0.0), mu_(1.0)
    {
        nk_ = kweights_.size();
    }

    ConstraintMethod method() const { return method_; }

    /// Weighted occupation total: sum_k w_k sum_i n_ik (same sum as in the constraint).
    double weighted_occupation_sum(const std::vector<double>& occ) const
    {
        assert(occ.size() == static_cast<size_t>(nk_ * nbands_));
        double sum = 0.0;
        for (int ik = 0; ik < nk_; ++ik)
            for (int i = 0; i < nbands_; ++i)
                sum += kweights_[ik] * occ[ik * nbands_ + i];
        return sum;
    }

    /// Constraint violation: c(n) = sum_k w_k sum_i n_ik - N_e
    double constraint_violation(const std::vector<double>& occ) const
    {
        return weighted_occupation_sum(occ) - n_electrons_;
    }

    /// Augmented Lagrangian penalty: lambda*c + mu/2*c^2
    double augmented_lagrangian_penalty(const std::vector<double>& occ) const
    {
        double c = constraint_violation(occ);
        return lambda_ * c + 0.5 * mu_ * c * c;
    }

    /// Gradient of the augmented Lagrangian penalty w.r.t. n_ik
    /// Returns (lambda + mu*c) * w_k for each (ik, i)
    void augmented_lagrangian_gradient(const std::vector<double>& occ,
                                       std::vector<double>& grad_penalty) const
    {
        double c = constraint_violation(occ);
        double factor = lambda_ + mu_ * c;
        grad_penalty.resize(occ.size());
        for (int ik = 0; ik < nk_; ++ik)
            for (int i = 0; i < nbands_; ++i)
                grad_penalty[ik * nbands_ + i] = factor * kweights_[ik];
    }

    /// Update Lagrange multiplier after inner optimization
    void update_multiplier(const std::vector<double>& occ)
    {
        double c = constraint_violation(occ);
        lambda_ += mu_ * c;
    }

    /// Increase penalty parameter
    void increase_penalty(double factor = 2.0, double max_mu = 1e6)
    {
        mu_ = std::min(mu_ * factor, max_mu);
    }

    /// Project occupations onto feasible set [0,1] with electron number constraint
    void project(std::vector<double>& occ) const
    {
        for (auto& n : occ)
            n = std::max(0.0, std::min(1.0, n));
        rescale_to_nel(occ);
    }

    /// Active set method: identify active constraints and solve reduced problem
    struct ActiveSetInfo
    {
        std::vector<bool> at_lower; // n = 0
        std::vector<bool> at_upper; // n = 1
        std::vector<bool> is_free;  // 0 < n < 1
        double lagrange_mult = 0.0; // for the equality constraint
    };

    ActiveSetInfo identify_active_set(const std::vector<double>& occ,
                                      const std::vector<double>& grad,
                                      double tol = 1e-8) const
    {
        ActiveSetInfo info;
        int N = occ.size();
        info.at_lower.resize(N, false);
        info.at_upper.resize(N, false);
        info.is_free.resize(N, false);

        for (int idx = 0; idx < N; ++idx)
        {
            int ik = idx / nbands_;
            if (occ[idx] <= tol && grad[idx] > 0)
                info.at_lower[idx] = true;
            else if (occ[idx] >= 1.0 - tol && grad[idx] < 0)
                info.at_upper[idx] = true;
            else
                info.is_free[idx] = true;
        }

        double sum_grad_w = 0.0;
        double sum_w2 = 0.0;
        for (int idx = 0; idx < N; ++idx)
        {
            if (info.is_free[idx])
            {
                int ik = idx / nbands_;
                sum_grad_w += grad[idx] * kweights_[ik];
                sum_w2 += kweights_[ik] * kweights_[ik];
            }
        }
        if (sum_w2 > 0)
            info.lagrange_mult = -sum_grad_w / sum_w2;

        return info;
    }

    /// Apply active set step: set gradient to zero for active constraints,
    /// subtract Lagrange multiplier projection for free variables
    void apply_active_set(const ActiveSetInfo& info,
                          std::vector<double>& grad) const
    {
        for (size_t idx = 0; idx < grad.size(); ++idx)
        {
            if (!info.is_free[idx])
            {
                grad[idx] = 0.0;
            }
            else
            {
                int ik = idx / nbands_;
                grad[idx] += info.lagrange_mult * kweights_[ik];
            }
        }
    }

    double lambda() const { return lambda_; }
    double mu() const { return mu_; }
    void set_mu(double mu) { mu_ = mu; }
    void set_lambda(double lambda) { lambda_ = lambda; }
    int nk() const { return nk_; }
    int nbands() const { return nbands_; }
    double n_electrons() const { return n_electrons_; }
    const std::vector<double>& kweights() const { return kweights_; }

  private:
    // Project occupations onto
    //   sum_k w_k * sum_i n_ki = N_e,  0 <= n_ki <= 1
    // by solving the 1D dual variable lambda:
    //   n_ki(lambda) = clip(n_ki^tmp - lambda * w_k, 0, 1).
    // The weighted sum is monotone decreasing in lambda, so we bracket
    // and solve by bisection.
    void rescale_to_nel(std::vector<double>& occ) const
    {
        assert(occ.size() == static_cast<size_t>(nk_ * nbands_));
        for (auto& n : occ)
            n = std::max(0.0, std::min(1.0, n));

        const double target = n_electrons_;
        const double tol_sum = 1e-12;
        const std::vector<double> occ_tmp(occ);

        auto weighted_sum_from_lambda = [&](double lambda)
        {
            double sum = 0.0;
            for (int ik = 0; ik < nk_; ++ik)
                for (int i = 0; i < nbands_; ++i)
                {
                    const double w = kweights_[ik];
                    const double n = occ_tmp[ik * nbands_ + i];
                    const double y = std::max(0.0, std::min(1.0, n - lambda * w));
                    sum += w * y;
                }
            return sum;
        };

        const double sum0 = weighted_sum_from_lambda(0.0);
        if (std::abs(sum0 - target) < tol_sum)
            return;

        double lam_lo = 0.0;
        double lam_hi = 0.0;
        double sum_lo = sum0;
        double sum_hi = sum0;
        const double expand_max = 1e12;

        if (sum0 > target)
        {
            lam_hi = 1.0;
            sum_hi = weighted_sum_from_lambda(lam_hi);
            while (sum_hi > target && lam_hi < expand_max)
            {
                lam_hi *= 2.0;
                sum_hi = weighted_sum_from_lambda(lam_hi);
            }
        }
        else
        {
            lam_lo = -1.0;
            sum_lo = weighted_sum_from_lambda(lam_lo);
            while (sum_lo < target && std::abs(lam_lo) < expand_max)
            {
                lam_lo *= 2.0;
                sum_lo = weighted_sum_from_lambda(lam_lo);
            }
        }

        if (!(sum_lo >= target && sum_hi <= target))
        {
            const double infeas_hi = weighted_sum_from_lambda(-expand_max);
            const double infeas_lo = weighted_sum_from_lambda(+expand_max);
            GlobalV::ofs_running << "WARNING: RDMFT rescale_to_nel: cannot bracket lambda for projection."
                                 << " target=" << target
                                 << " reachable_range=[" << infeas_lo << ", " << infeas_hi << "]"
                                 << std::endl;
            const double lambda = (std::abs(target - infeas_lo) < std::abs(target - infeas_hi)) ? expand_max : -expand_max;
            for (int ik = 0; ik < nk_; ++ik)
                for (int i = 0; i < nbands_; ++i)
                {
                    const double w = kweights_[ik];
                    occ[ik * nbands_ + i] = std::max(0.0, std::min(1.0, occ_tmp[ik * nbands_ + i] - lambda * w));
                }
            return;
        }

        for (int it = 0; it < 100; ++it)
        {
            const double lam_mid = 0.5 * (lam_lo + lam_hi);
            const double sum_mid = weighted_sum_from_lambda(lam_mid);
            if (std::abs(sum_mid - target) < tol_sum)
            {
                lam_lo = lam_mid;
                lam_hi = lam_mid;
                break;
            }
            if (sum_mid > target)
            {
                lam_lo = lam_mid;
            }
            else
            {
                lam_hi = lam_mid;
            }
        }

        const double lambda = 0.5 * (lam_lo + lam_hi);
        for (int ik = 0; ik < nk_; ++ik)
            for (int i = 0; i < nbands_; ++i)
            {
                const double w = kweights_[ik];
                occ[ik * nbands_ + i] = std::max(0.0, std::min(1.0, occ_tmp[ik * nbands_ + i] - lambda * w));
            }
    }

    ConstraintMethod method_;
    double n_electrons_;
    std::vector<double> kweights_;
    int nbands_;
    int nk_;
    double lambda_;
    double mu_;
};


/// Unconstrained occupation parameterization via sigmoid with adaptive shift.
///
/// Given unconstrained z ∈ R^{Nk×Nb}, the shift λ ∈ R is determined each
/// step by bisection from:
///   Σ_k w_k Σ_i σ(z_{ik} + λ) = N_e,   σ(x) = 1 / (1 + e^{-x})
/// then n_{ik} = σ(z_{ik} + λ) automatically satisfies 0 < n < 1 and the
/// electron-count constraint.  No augmented-Lagrangian penalty is needed in
/// the joint (product-manifold) strategy.
///
/// λ is an implicit function of z from the electron-count equation.  With
/// h_{ik} = ∂E/∂n_{ik} (ABACUS `grad_occ`, already including k-point weights
/// in the usual RDMFT energy derivative),
///   ∂E/∂z_{ik} = σ'_{ik} ( h_{ik} - w_k · (Σ_{jq} h_{jq} σ'_{jq}) / (Σ_{jq} w_j σ'_{jq}) ),
/// σ'_{ik} = σ(z_{ik}+λ)(1−σ(...)).  If Σ w σ' ≈ 0 (saturated occupations), ∂E/∂z is set to 0.
class SigmaShiftOccParam
{
  public:
    SigmaShiftOccParam(double n_electrons,
                       const std::vector<double>& kweights,
                       int nbands)
        : n_electrons_(n_electrons), kweights_(kweights),
          nbands_(nbands), nk_(static_cast<int>(kweights.size())),
          lambda_(0.0)
    {}

    /// Numerically stable σ(x); avoids exp overflow for |x| large.
    static double stable_sigmoid(double x)
    {
        if (x >= 0.0)
        {
            const double z = std::exp(-x);
            return 1.0 / (1.0 + z);
        }
        const double z = std::exp(x);
        return z / (1.0 + z);
    }

    /// σ'(x) using stable σ.
    static double stable_sigmoid_prime(double x)
    {
        const double s = stable_sigmoid(x);
        return s * (1.0 - s);
    }

    /// Find λ by bisection such that Σ_k w_k Σ_i σ(z_{ik}+λ) = N_e.
    /// Stores and returns the computed λ.
    double compute_lambda(const std::vector<double>& z_params)
    {
        assert(static_cast<int>(z_params.size()) == nk_ * nbands_);
        // Convergence tolerance: electron-count error smaller than ~1e-12 electrons.
        const double tol = 1e-12;
        // Safety cap: |λ| > 1e6 would push all σ(z+λ) to 0 or 1, so the
        // weighted sum is bounded and the bisection must have converged by then.
        const double lambda_bound = 1e6;
        // 100 bisection steps give 2^{-100} ≈ 1e-30 accuracy in the bracket,
        // which is far better than the tolerance above.
        const int max_bisect_iter = 100;

        auto weighted_sum = [&](double lam) -> double {
            double s = 0.0;
            for (int ik = 0; ik < nk_; ++ik)
                for (int ib = 0; ib < nbands_; ++ib)
                {
                    const double x = z_params[ik * nbands_ + ib] + lam;
                    s += kweights_[ik] * stable_sigmoid(x);
                }
            return s;
        };

        const double s0 = weighted_sum(0.0);
        if (std::abs(s0 - n_electrons_) < tol)
        {
            lambda_ = 0.0;
            return lambda_;
        }

        double lam_lo, lam_hi;
        if (s0 < n_electrons_)
        {
            // Need to shift sigmoid right (increase all occupations).
            lam_lo = 0.0;
            lam_hi = 1.0;
            while (weighted_sum(lam_hi) < n_electrons_ && lam_hi < lambda_bound)
                lam_hi *= 2.0;
        }
        else
        {
            // Need to shift sigmoid left (decrease all occupations).
            lam_lo = -1.0;
            lam_hi = 0.0;
            while (weighted_sum(lam_lo) > n_electrons_ && lam_lo > -lambda_bound)
                lam_lo *= 2.0;
        }

        for (int it = 0; it < max_bisect_iter; ++it)
        {
            const double lam_mid = 0.5 * (lam_lo + lam_hi);
            const double s_mid = weighted_sum(lam_mid);
            if (std::abs(s_mid - n_electrons_) < tol)
            {
                lam_lo = lam_hi = lam_mid;
                break;
            }
            if (s_mid < n_electrons_)
                lam_lo = lam_mid;
            else
                lam_hi = lam_mid;
        }

        lambda_ = 0.5 * (lam_lo + lam_hi);
        return lambda_;
    }

    /// Map z → n using given lambda: n_i = σ(z_i + lambda).
    void params_to_occ(const std::vector<double>& z,
                       double lambda,
                       std::vector<double>& occ) const
    {
        occ.resize(z.size());
        for (size_t i = 0; i < z.size(); ++i)
            occ[i] = stable_sigmoid(z[i] + lambda);
    }

    /// Map n → z (initial parameterization): z_i = logit(n_i) = log(n/(1-n)).
    void occ_to_params(const std::vector<double>& occ,
                       std::vector<double>& z) const
    {
        z.resize(occ.size());
        for (size_t i = 0; i < occ.size(); ++i)
        {
            // Clamp strictly away from {0,1} to keep logit finite.
            const double n = std::max(1e-12, std::min(1.0 - 1e-12, occ[i]));
            z[i] = std::log(n / (1.0 - n));
        }
    }

    /// ∂E/∂z with λ(z) from Σ_k w_k Σ_i σ(z_{ik}+λ) = N_e (implicit differentiation).
    void transform_gradient_batch(const std::vector<double>& dE_dn,
                                  const std::vector<double>& z,
                                  double lambda,
                                  std::vector<double>& dE_dz) const
    {
        assert(static_cast<int>(dE_dn.size()) == nk_ * nbands_);
        assert(static_cast<int>(z.size()) == nk_ * nbands_);
        dE_dz.resize(dE_dn.size());

        double sum_s = 0.0;
        double sum_hsp = 0.0;
        for (int ik = 0; ik < nk_; ++ik)
        {
            const double wk = kweights_[ik];
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const int idx = ik * nbands_ + ib;
                const double sp = stable_sigmoid_prime(z[idx] + lambda);
                sum_s += wk * sp;
                sum_hsp += dE_dn[idx] * sp;
            }
        }

        constexpr double sum_s_floor = 1e-12;
        if (std::abs(sum_s) < sum_s_floor)
        {
            std::fill(dE_dz.begin(), dE_dz.end(), 0.0);
            return;
        }

        const double ratio = sum_hsp / sum_s;
        for (int ik = 0; ik < nk_; ++ik)
        {
            const double wk = kweights_[ik];
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const int idx = ik * nbands_ + ib;
                const double sp = stable_sigmoid_prime(z[idx] + lambda);
                dE_dz[idx] = sp * (dE_dn[idx] - wk * ratio);
            }
        }
    }

    double lambda() const { return lambda_; }

  private:
    double n_electrons_;
    std::vector<double> kweights_;
    int nbands_;
    int nk_;
    double lambda_;
};

/// Add augmented-Lagrangian penalty gradient in-place:
///   grad_occ <- grad_occ + d/dn [lambda*c + 0.5*mu*c^2]
/// where c = sum_k w_k sum_i n_ki - N_e.
inline void add_augmented_lagrangian_occ_gradient(const OccupationConstraint& constraint,
                                                  const std::vector<double>& occ,
                                                  std::vector<double>& grad_occ)
{
    assert(occ.size() == grad_occ.size());
    std::vector<double> grad_pen;
    constraint.augmented_lagrangian_gradient(occ, grad_pen);
    assert(grad_pen.size() == grad_occ.size());
    for (size_t i = 0; i < grad_occ.size(); ++i)
        grad_occ[i] += grad_pen[i];
}

} // namespace rdmft

#endif // RDMFT_OCCUPATION_H
