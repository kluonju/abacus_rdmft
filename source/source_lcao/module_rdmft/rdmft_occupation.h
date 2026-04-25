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
        else
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
        else
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
        else
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

} // namespace rdmft

#endif // RDMFT_OCCUPATION_H
