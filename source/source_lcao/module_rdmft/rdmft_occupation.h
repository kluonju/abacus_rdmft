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
///
/// When constructed with a spin map (`isk` of length `nk`, with values in
/// [0, nspin)), the class enforces nspin **separate** equality constraints:
/// sum_k w_k sum_i n_iks = N_s for s = 0..nspin-1.  This is equivalent to
/// having two effective chemical potentials (one per spin) when nspin=2 with
/// fixed N_↑ and N_↓, mirroring ABACUS KS's split-Fermi treatment when
/// `nupdown` is set (see docs/advanced/scf/spin.md).  The single-constraint
/// behavior is recovered when `isk` is empty / nspin==1.
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

    /// Spin-resolved constructor: nspin separate equality constraints.
    /// `isk[ik]` is the spin index in [0, nspin) for k-spin row `ik`, and
    /// `n_electrons_per_spin[s]` is the per-spin equality target. The
    /// "total" target n_electrons (used for diagnostics / fallbacks) is the
    /// sum of `n_electrons_per_spin`.
    OccupationConstraint(ConstraintMethod method,
                         const std::vector<double>& kweights,
                         int nbands,
                         int nspin,
                         const std::vector<int>& isk,
                         const std::vector<double>& n_electrons_per_spin)
        : method_(method), kweights_(kweights), nbands_(nbands),
          lambda_(0.0), mu_(1.0),
          spin_resolved_(nspin > 1),
          nspin_(nspin), isk_(isk),
          n_electrons_per_spin_(n_electrons_per_spin),
          lambda_per_spin_(nspin, 0.0),
          mu_per_spin_(nspin, 1.0)
    {
        nk_ = kweights_.size();
        n_electrons_ = 0.0;
        for (double v : n_electrons_per_spin_) n_electrons_ += v;
        // For nspin==1 fall back to the single-constraint code path.
        if (nspin <= 1)
        {
            spin_resolved_ = false;
            isk_.clear();
            n_electrons_per_spin_.clear();
        }
    }

    ConstraintMethod method() const { return method_; }
    bool spin_resolved() const { return spin_resolved_; }
    int nspin() const { return nspin_; }
    const std::vector<int>& isk() const { return isk_; }
    const std::vector<double>& n_electrons_per_spin() const { return n_electrons_per_spin_; }

    /// Weighted occupation total: sum_k w_k sum_i n_ik (sum over all (k, spin)).
    double weighted_occupation_sum(const std::vector<double>& occ) const
    {
        assert(occ.size() == static_cast<size_t>(nk_ * nbands_));
        double sum = 0.0;
        for (int ik = 0; ik < nk_; ++ik)
            for (int i = 0; i < nbands_; ++i)
                sum += kweights_[ik] * occ[ik * nbands_ + i];
        return sum;
    }

    /// Per-spin weighted sum (size nspin); same as weighted_occupation_sum
    /// with a single entry when not spin-resolved.
    std::vector<double> weighted_occupation_sum_per_spin(const std::vector<double>& occ) const
    {
        assert(occ.size() == static_cast<size_t>(nk_ * nbands_));
        if (!spin_resolved_)
        {
            return std::vector<double>{weighted_occupation_sum(occ)};
        }
        std::vector<double> sums(nspin_, 0.0);
        for (int ik = 0; ik < nk_; ++ik)
        {
            const int is = isk_[ik];
            for (int i = 0; i < nbands_; ++i)
                sums[is] += kweights_[ik] * occ[ik * nbands_ + i];
        }
        return sums;
    }

    /// Per-spin equality residual c_s(n) = sum_k(s) w_k sum_i n_iks - N_s.
    /// For non-spin-resolved instances returns a single-element vector.
    std::vector<double> constraint_violation_per_spin(const std::vector<double>& occ) const
    {
        std::vector<double> sums = weighted_occupation_sum_per_spin(occ);
        if (!spin_resolved_)
        {
            sums[0] -= n_electrons_;
            return sums;
        }
        for (int s = 0; s < nspin_; ++s)
            sums[s] -= n_electrons_per_spin_[s];
        return sums;
    }

    /// Aggregated equality residual:
    /// - non-spin-resolved: c(n) = sum_k w_k sum_i n_ik - N_e (scalar);
    /// - spin-resolved: sum_s c_s (preserves the legacy behavior of
    ///   monitoring a single number, e.g. for the |c| diagnostic line). Use
    ///   `constraint_violation_per_spin` to inspect each spin separately.
    double constraint_violation(const std::vector<double>& occ) const
    {
        if (!spin_resolved_)
            return weighted_occupation_sum(occ) - n_electrons_;
        const std::vector<double> cs = constraint_violation_per_spin(occ);
        double c = 0.0;
        for (double v : cs) c += v;
        return c;
    }

    /// Augmented Lagrangian penalty: sum_s [lambda_s c_s + 0.5 mu_s c_s^2].
    double augmented_lagrangian_penalty(const std::vector<double>& occ) const
    {
        if (!spin_resolved_)
        {
            double c = constraint_violation(occ);
            return lambda_ * c + 0.5 * mu_ * c * c;
        }
        const std::vector<double> cs = constraint_violation_per_spin(occ);
        double pen = 0.0;
        for (int s = 0; s < nspin_; ++s)
            pen += lambda_per_spin_[s] * cs[s] + 0.5 * mu_per_spin_[s] * cs[s] * cs[s];
        return pen;
    }

    /// Gradient of the augmented Lagrangian penalty w.r.t. n_ik.
    /// Single-constraint: (lambda + mu*c) * w_k for every (ik, i).
    /// Per-spin: (lambda_s + mu_s * c_s) * w_k where s = isk_[ik].
    void augmented_lagrangian_gradient(const std::vector<double>& occ,
                                       std::vector<double>& grad_penalty) const
    {
        grad_penalty.resize(occ.size());
        if (!spin_resolved_)
        {
            const double c = constraint_violation(occ);
            const double factor = lambda_ + mu_ * c;
            for (int ik = 0; ik < nk_; ++ik)
                for (int i = 0; i < nbands_; ++i)
                    grad_penalty[ik * nbands_ + i] = factor * kweights_[ik];
            return;
        }
        const std::vector<double> cs = constraint_violation_per_spin(occ);
        std::vector<double> factor(nspin_, 0.0);
        for (int s = 0; s < nspin_; ++s)
            factor[s] = lambda_per_spin_[s] + mu_per_spin_[s] * cs[s];
        for (int ik = 0; ik < nk_; ++ik)
        {
            const int is = isk_[ik];
            const double f = factor[is];
            for (int i = 0; i < nbands_; ++i)
                grad_penalty[ik * nbands_ + i] = f * kweights_[ik];
        }
    }

    /// Update Lagrange multiplier(s) after an inner optimization. For the
    /// spin-resolved case, λ_s ← λ_s + μ_s c_s independently per spin.
    void update_multiplier(const std::vector<double>& occ)
    {
        if (!spin_resolved_)
        {
            const double c = constraint_violation(occ);
            lambda_ += mu_ * c;
            return;
        }
        const std::vector<double> cs = constraint_violation_per_spin(occ);
        for (int s = 0; s < nspin_; ++s)
            lambda_per_spin_[s] += mu_per_spin_[s] * cs[s];
    }

    /// Increase penalty parameter(s)
    void increase_penalty(double factor = 2.0, double max_mu = 1e6)
    {
        mu_ = std::min(mu_ * factor, max_mu);
        for (auto& m : mu_per_spin_) m = std::min(m * factor, max_mu);
    }

    /// Project occupations onto feasible set [0,1] with electron number constraint.
    /// Pre-clips the input to [0,1] before the bisection-based rescale, suitable
    /// for cleaning up an already-feasible occupation that has drifted slightly
    /// out of range. **Not** the L2 proximal projection of an arbitrary point:
    /// the pre-clip discards how far n_ki is below 0 or above 1, which the
    /// proximal projection needs.  Use `proximal_project` for SPG steps.
    void project(std::vector<double>& occ) const
    {
        for (auto& n : occ)
            n = std::max(0.0, std::min(1.0, n));
        rescale_to_nel_dispatch(occ);
    }

    /// Same as `project` but uses a uniform (k-weight-independent) shift.
    /// **Not** the L2 proximal projection of an arbitrary point.  Use
    /// `proximal_project_uniform` for SPG steps.
    void project_uniform(std::vector<double>& occ) const
    {
        for (auto& n : occ)
            n = std::max(0.0, std::min(1.0, n));
        rescale_to_nel_uniform_dispatch(occ);
    }

    /// L2 proximal projection (wk-proportional shift) onto
    ///   { y in [0,1]^N : sum_k w_k sum_i y_ki = N_e }.
    ///
    ///   y = argmin_y ||y - x||^2  s.t.  y in [0,1], sum_k w_k y = N_e
    ///
    /// Solved via the dual: y_ki(λ) = clip(x_ki - λ w_k, 0, 1), bisecting λ to
    /// satisfy the equality constraint.  Unlike `project`, this does NOT
    /// pre-clip x; the bisection sees the un-clipped input so the convex hull
    /// of the box constraint is respected.  This is the projection that makes
    /// the SPG step `y = P(x - α g)` a descent direction (Fejer property),
    /// and at a KKT point of E with this metric, P(x - α g) = x exactly.
    ///
    /// In spin-resolved mode (nspin == 2 with `isk` set), the projection is
    /// separable across spin blocks (the two equality constraints touch
    /// disjoint coordinates), so we run an independent 1D bisection per
    /// spin and the result is still the L2 nearest feasible point.
    void proximal_project(std::vector<double>& occ) const
    {
        rescale_to_nel_dispatch(occ);
    }

    /// L2 proximal projection with a uniform shift (k-weight-independent).
    ///   y_ki = clip(x_ki - mu, 0, 1),  mu chosen so sum_k w_k sum_i y_ki = N_e.
    /// This is the proximal projection in the metric ||y||^2 = sum_ki y^2
    /// (no w_k weighting on the diagonal); it pairs naturally with SPG using
    /// the unweighted gradient eps_ki = (∂E/∂n_ki) / w_k.  Unlike
    /// `project_uniform`, this does NOT pre-clip, preserving the proximal
    /// geometry required for SPG descent and KKT preservation.
    ///
    /// Spin-resolved variant: a separate uniform shift is solved per spin
    /// block (independent dual problems).
    void proximal_project_uniform(std::vector<double>& occ) const
    {
        rescale_to_nel_uniform_dispatch(occ);
    }

    /// Active set method: identify active constraints and solve reduced problem
    struct ActiveSetInfo
    {
        std::vector<bool> at_lower; // n = 0
        std::vector<bool> at_upper; // n = 1
        std::vector<bool> is_free;  // 0 < n < 1
        double lagrange_mult = 0.0; // for the equality constraint (single)
        /// Per-spin multipliers when the constraint is spin-resolved.
        /// Empty otherwise; `lagrange_mult` then holds max abs as a single
        /// number for diagnostic logging.
        std::vector<double> lagrange_mult_per_spin;
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
            if (occ[idx] <= tol && grad[idx] > 0)
                info.at_lower[idx] = true;
            else if (occ[idx] >= 1.0 - tol && grad[idx] < 0)
                info.at_upper[idx] = true;
            else
                info.is_free[idx] = true;
        }

        if (!spin_resolved_)
        {
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
        }
        else
        {
            // Two equality constraints (one per spin); independent multipliers.
            std::vector<double> sum_grad_w(nspin_, 0.0);
            std::vector<double> sum_w2(nspin_, 0.0);
            for (int idx = 0; idx < N; ++idx)
            {
                if (info.is_free[idx])
                {
                    int ik = idx / nbands_;
                    int is = isk_[ik];
                    sum_grad_w[is] += grad[idx] * kweights_[ik];
                    sum_w2[is] += kweights_[ik] * kweights_[ik];
                }
            }
            info.lagrange_mult_per_spin.assign(nspin_, 0.0);
            for (int s = 0; s < nspin_; ++s)
                if (sum_w2[s] > 0.0)
                    info.lagrange_mult_per_spin[s] = -sum_grad_w[s] / sum_w2[s];
            // Keep the scalar field for back-compat logging (max abs).
            info.lagrange_mult = 0.0;
            for (int s = 0; s < nspin_; ++s)
                info.lagrange_mult = std::max(info.lagrange_mult,
                                              std::abs(info.lagrange_mult_per_spin[s]));
        }

        return info;
    }

    /// Apply active set step: set gradient to zero for active constraints,
    /// subtract Lagrange multiplier projection for free variables (per-spin
    /// when applicable).
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
                int ik = static_cast<int>(idx) / nbands_;
                if (!spin_resolved_)
                {
                    grad[idx] += info.lagrange_mult * kweights_[ik];
                }
                else
                {
                    const int is = isk_[ik];
                    grad[idx] += info.lagrange_mult_per_spin[is] * kweights_[ik];
                }
            }
        }
    }

    double lambda() const { return lambda_; }
    double mu() const { return mu_; }
    void set_mu(double mu)
    {
        mu_ = mu;
        std::fill(mu_per_spin_.begin(), mu_per_spin_.end(), mu);
    }
    void set_lambda(double lambda)
    {
        lambda_ = lambda;
        std::fill(lambda_per_spin_.begin(), lambda_per_spin_.end(), lambda);
    }
    /// Per-spin multipliers (empty when not spin-resolved).
    const std::vector<double>& lambda_per_spin() const { return lambda_per_spin_; }
    const std::vector<double>& mu_per_spin() const { return mu_per_spin_; }
    int nk() const { return nk_; }
    int nbands() const { return nbands_; }
    double n_electrons() const { return n_electrons_; }
    const std::vector<double>& kweights() const { return kweights_; }

  private:
    /// Dispatcher: route to single- or per-spin rescale depending on mode.
    void rescale_to_nel_dispatch(std::vector<double>& occ) const
    {
        if (!spin_resolved_)
        {
            rescale_to_nel(occ);
            return;
        }
        // Independent per-spin bisection.
        for (int s = 0; s < nspin_; ++s)
            rescale_to_nel_one_spin(occ, s, /*uniform_shift=*/false);
    }

    void rescale_to_nel_uniform_dispatch(std::vector<double>& occ) const
    {
        if (!spin_resolved_)
        {
            rescale_to_nel_uniform(occ);
            return;
        }
        for (int s = 0; s < nspin_; ++s)
            rescale_to_nel_one_spin(occ, s, /*uniform_shift=*/true);
    }

    /// Per-spin projection: for every (ik with isk[ik]==spin, ib) entry,
    /// solve a 1D dual to satisfy
    ///   sum_{ik in spin} w_k sum_i clip(n_iks - shift * (uniform ? 1 : w_k), 0, 1) = N_s.
    /// Implementation mirrors `rescale_to_nel` / `rescale_to_nel_uniform`,
    /// restricted to indices belonging to the given spin block.
    void rescale_to_nel_one_spin(std::vector<double>& occ,
                                 int spin,
                                 bool uniform_shift) const
    {
        assert(occ.size() == static_cast<size_t>(nk_ * nbands_));
        const double target = n_electrons_per_spin_[spin];
        const double tol_sum = 1e-12;
        const std::vector<double> occ_tmp(occ);

        auto weighted_sum = [&](double shift) -> double {
            double sum = 0.0;
            for (int ik = 0; ik < nk_; ++ik)
            {
                if (isk_[ik] != spin) continue;
                const double w = kweights_[ik];
                for (int i = 0; i < nbands_; ++i)
                {
                    const double n = occ_tmp[ik * nbands_ + i];
                    const double sub = uniform_shift ? shift : shift * w;
                    const double y = std::max(0.0, std::min(1.0, n - sub));
                    sum += w * y;
                }
            }
            return sum;
        };

        const double sum0 = weighted_sum(0.0);
        if (std::abs(sum0 - target) < tol_sum)
            return;

        double lo = 0.0, hi = 0.0;
        double sum_lo = sum0, sum_hi = sum0;
        const double expand_max = 1.0e12;
        if (sum0 > target)
        {
            hi = 1.0;
            sum_hi = weighted_sum(hi);
            while (sum_hi > target && hi < expand_max)
            {
                hi *= 2.0;
                sum_hi = weighted_sum(hi);
            }
        }
        else
        {
            lo = -1.0;
            sum_lo = weighted_sum(lo);
            while (sum_lo < target && std::abs(lo) < expand_max)
            {
                lo *= 2.0;
                sum_lo = weighted_sum(lo);
            }
        }

        double shift_final = 0.0;
        if (!(sum_lo >= target && sum_hi <= target))
        {
            // Saturated; clamp to whichever side is closer to target.
            shift_final = (std::abs(target - weighted_sum(-expand_max))
                            < std::abs(target - weighted_sum(+expand_max)))
                ? -expand_max : expand_max;
        }
        else
        {
            for (int it = 0; it < 100; ++it)
            {
                const double mid = 0.5 * (lo + hi);
                const double sm = weighted_sum(mid);
                if (std::abs(sm - target) < tol_sum)
                {
                    lo = hi = mid;
                    break;
                }
                if (sm > target) lo = mid; else hi = mid;
            }
            shift_final = 0.5 * (lo + hi);
        }

        for (int ik = 0; ik < nk_; ++ik)
        {
            if (isk_[ik] != spin) continue;
            const double w = kweights_[ik];
            for (int i = 0; i < nbands_; ++i)
            {
                const double n = occ_tmp[ik * nbands_ + i];
                const double sub = uniform_shift ? shift_final : shift_final * w;
                occ[ik * nbands_ + i] = std::max(0.0, std::min(1.0, n - sub));
            }
        }
    }

    // Project occupations onto
    //   sum_k w_k * sum_i n_ki = N_e,  0 <= n_ki <= 1
    // by solving the 1D dual variable lambda:
    //   n_ki(lambda) = clip(n_ki^tmp - lambda * w_k, 0, 1).
    // The weighted sum is monotone decreasing in lambda, so we bracket
    // and solve by bisection.
    void rescale_to_nel(std::vector<double>& occ) const
    {
        assert(occ.size() == static_cast<size_t>(nk_ * nbands_));
        // Note: do NOT pre-clip occ here.  The bisection below applies the
        // box clip [0,1] inside `weighted_sum_from_lambda`, so passing the
        // un-clipped input is what makes this an L2 proximal projection.
        // Callers that want to pre-clip (e.g. `project()` cleaning up a
        // drifted feasible occupation) do so before calling this routine.

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

    void rescale_to_nel_uniform(std::vector<double>& occ) const
    {
        assert(occ.size() == static_cast<size_t>(nk_ * nbands_));
        // No internal pre-clip; see comment in `rescale_to_nel`.

        const double target = n_electrons_;
        const double tol_sum = 1e-12;
        const std::vector<double> occ_tmp(occ);

        auto weighted_sum_from_mu = [&](double mu)
        {
            double sum = 0.0;
            for (int ik = 0; ik < nk_; ++ik)
                for (int i = 0; i < nbands_; ++i)
                {
                    const double y = std::max(0.0, std::min(1.0, occ_tmp[ik * nbands_ + i] - mu));
                    sum += kweights_[ik] * y;
                }
            return sum;
        };

        const double sum0 = weighted_sum_from_mu(0.0);
        if (std::abs(sum0 - target) < tol_sum)
            return;

        double mu_lo = 0.0, mu_hi = 0.0;
        double sum_lo = sum0, sum_hi = sum0;
        // Allow large expansion because the un-clipped input to the proximal
        // projection can be far outside [0,1] for SPG steps with
        // α·ε > O(1) on some bands.
        const double expand_max = 1.0e12;

        if (sum0 > target)
        {
            mu_hi = 1.0;
            sum_hi = weighted_sum_from_mu(mu_hi);
            while (sum_hi > target && mu_hi < expand_max)
            {
                mu_hi *= 2.0;
                sum_hi = weighted_sum_from_mu(mu_hi);
            }
        }
        else
        {
            mu_lo = -1.0;
            sum_lo = weighted_sum_from_mu(mu_lo);
            while (sum_lo < target && std::abs(mu_lo) < expand_max)
            {
                mu_lo *= 2.0;
                sum_lo = weighted_sum_from_mu(mu_lo);
            }
        }

        if (!(sum_lo >= target && sum_hi <= target))
        {
            // Saturated: clamp mu to the side closer to feasibility.
            const double mu_sat = (sum_lo < target) ? mu_lo : mu_hi;
            for (int ik = 0; ik < nk_; ++ik)
                for (int i = 0; i < nbands_; ++i)
                    occ[ik * nbands_ + i] = std::max(0.0, std::min(1.0, occ_tmp[ik * nbands_ + i] - mu_sat));
            return;
        }

        for (int it = 0; it < 100; ++it)
        {
            const double mu_mid = 0.5 * (mu_lo + mu_hi);
            const double sum_mid = weighted_sum_from_mu(mu_mid);
            if (std::abs(sum_mid - target) < tol_sum)
            {
                mu_lo = mu_mid;
                mu_hi = mu_mid;
                break;
            }
            if (sum_mid > target)
                mu_lo = mu_mid;
            else
                mu_hi = mu_mid;
        }

        const double mu = 0.5 * (mu_lo + mu_hi);
        for (int ik = 0; ik < nk_; ++ik)
            for (int i = 0; i < nbands_; ++i)
                occ[ik * nbands_ + i] = std::max(0.0, std::min(1.0, occ_tmp[ik * nbands_ + i] - mu));
    }

    ConstraintMethod method_;
    double n_electrons_;
    std::vector<double> kweights_;
    int nbands_;
    int nk_;
    double lambda_;
    double mu_;

    // Spin-resolved equality-constraint state (active when nspin > 1).
    bool spin_resolved_ = false;
    int nspin_ = 1;
    std::vector<int> isk_;                     // length nk_, values in [0, nspin_)
    std::vector<double> n_electrons_per_spin_; // size nspin_
    std::vector<double> lambda_per_spin_;      // size nspin_
    std::vector<double> mu_per_spin_;          // size nspin_
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

    /// Spin-resolved sigma-shift: a separate shift λ_s solves
    ///   sum_k(s) w_k sum_i σ(z_iks + λ_s) = N_s
    /// independently per spin block (s = 0..nspin-1). When nspin == 1 this
    /// reduces to the single-shift constructor above.
    SigmaShiftOccParam(const std::vector<double>& kweights,
                       int nbands,
                       int nspin,
                       const std::vector<int>& isk,
                       const std::vector<double>& n_electrons_per_spin)
        : n_electrons_(0.0), kweights_(kweights),
          nbands_(nbands), nk_(static_cast<int>(kweights.size())),
          lambda_(0.0),
          spin_resolved_(nspin > 1),
          nspin_(nspin), isk_(isk),
          n_electrons_per_spin_(n_electrons_per_spin),
          lambda_per_spin_(nspin, 0.0)
    {
        for (double v : n_electrons_per_spin_) n_electrons_ += v;
        if (nspin <= 1)
        {
            spin_resolved_ = false;
            isk_.clear();
            n_electrons_per_spin_.clear();
            lambda_per_spin_.clear();
        }
    }

    bool spin_resolved() const { return spin_resolved_; }
    int nspin() const { return nspin_; }
    const std::vector<double>& lambda_per_spin() const { return lambda_per_spin_; }

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
    /// Stores and returns the computed λ. In spin-resolved mode, runs a
    /// separate bisection per spin block (filling `lambda_per_spin_`) and
    /// returns λ_↑ (with λ_↓ in `lambda_per_spin_[1]`).
    double compute_lambda(const std::vector<double>& z_params)
    {
        assert(static_cast<int>(z_params.size()) == nk_ * nbands_);
        if (!spin_resolved_)
        {
            lambda_ = solve_one_lambda(z_params, /*spin=*/-1, n_electrons_);
            return lambda_;
        }
        for (int s = 0; s < nspin_; ++s)
        {
            lambda_per_spin_[s] = solve_one_lambda(z_params, s, n_electrons_per_spin_[s]);
        }
        // Use λ_↑ for the legacy scalar return; downstream code that needs
        // both should query lambda_per_spin().
        lambda_ = lambda_per_spin_[0];
        return lambda_;
    }

  private:
    /// Bisection helper: solve `weighted_sum(λ) = target`. When spin = -1,
    /// the sum runs over all (ik, ib); otherwise only over `isk_[ik] == spin`.
    double solve_one_lambda(const std::vector<double>& z_params,
                            int spin,
                            double target) const
    {
        const double tol = 1e-12;
        const double lambda_bound = 1e6;
        const int max_bisect_iter = 100;

        auto weighted_sum = [&](double lam) -> double {
            double s = 0.0;
            for (int ik = 0; ik < nk_; ++ik)
            {
                if (spin >= 0 && isk_[ik] != spin) continue;
                for (int ib = 0; ib < nbands_; ++ib)
                {
                    const double x = z_params[ik * nbands_ + ib] + lam;
                    s += kweights_[ik] * stable_sigmoid(x);
                }
            }
            return s;
        };

        const double s0 = weighted_sum(0.0);
        if (std::abs(s0 - target) < tol)
            return 0.0;

        double lam_lo, lam_hi;
        if (s0 < target)
        {
            lam_lo = 0.0;
            lam_hi = 1.0;
            while (weighted_sum(lam_hi) < target && lam_hi < lambda_bound)
                lam_hi *= 2.0;
        }
        else
        {
            lam_lo = -1.0;
            lam_hi = 0.0;
            while (weighted_sum(lam_lo) > target && lam_lo > -lambda_bound)
                lam_lo *= 2.0;
        }

        for (int it = 0; it < max_bisect_iter; ++it)
        {
            const double lam_mid = 0.5 * (lam_lo + lam_hi);
            const double s_mid = weighted_sum(lam_mid);
            if (std::abs(s_mid - target) < tol)
            {
                lam_lo = lam_hi = lam_mid;
                break;
            }
            if (s_mid < target)
                lam_lo = lam_mid;
            else
                lam_hi = lam_mid;
        }
        return 0.5 * (lam_lo + lam_hi);
    }

  public:

    /// Map z → n using given lambda: n_i = σ(z_i + lambda).
    /// In spin-resolved mode the per-spin λ_s in `lambda_per_spin_` is used
    /// regardless of the `lambda` argument (which is left for back-compat).
    void params_to_occ(const std::vector<double>& z,
                       double lambda,
                       std::vector<double>& occ) const
    {
        occ.resize(z.size());
        if (!spin_resolved_)
        {
            for (size_t i = 0; i < z.size(); ++i)
                occ[i] = stable_sigmoid(z[i] + lambda);
            return;
        }
        for (int ik = 0; ik < nk_; ++ik)
        {
            const double lam_s = lambda_per_spin_[isk_[ik]];
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const int idx = ik * nbands_ + ib;
                occ[idx] = stable_sigmoid(z[idx] + lam_s);
            }
        }
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
    /// In spin-resolved mode each spin block has its own implicit constraint,
    /// so the chain-rule subtraction `w_k * (Σ h σ' / Σ w σ')` is computed
    /// per spin and only applied within that spin's block.
    void transform_gradient_batch(const std::vector<double>& dE_dn,
                                  const std::vector<double>& z,
                                  double lambda,
                                  std::vector<double>& dE_dz) const
    {
        assert(static_cast<int>(dE_dn.size()) == nk_ * nbands_);
        assert(static_cast<int>(z.size()) == nk_ * nbands_);
        dE_dz.resize(dE_dn.size());

        constexpr double sum_s_floor = 1e-12;

        if (!spin_resolved_)
        {
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
            return;
        }

        // Per-spin: separate implicit-derivative ratios.
        std::vector<double> sum_s(nspin_, 0.0);
        std::vector<double> sum_hsp(nspin_, 0.0);
        for (int ik = 0; ik < nk_; ++ik)
        {
            const int is = isk_[ik];
            const double wk = kweights_[ik];
            const double lam_s = lambda_per_spin_[is];
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const int idx = ik * nbands_ + ib;
                const double sp = stable_sigmoid_prime(z[idx] + lam_s);
                sum_s[is] += wk * sp;
                sum_hsp[is] += dE_dn[idx] * sp;
            }
        }
        std::vector<double> ratio(nspin_, 0.0);
        for (int s = 0; s < nspin_; ++s)
            if (std::abs(sum_s[s]) >= sum_s_floor)
                ratio[s] = sum_hsp[s] / sum_s[s];
        for (int ik = 0; ik < nk_; ++ik)
        {
            const int is = isk_[ik];
            const double wk = kweights_[ik];
            const double lam_s = lambda_per_spin_[is];
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const int idx = ik * nbands_ + ib;
                const double sp = stable_sigmoid_prime(z[idx] + lam_s);
                dE_dz[idx] = sp * (dE_dn[idx] - wk * ratio[is]);
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

    // Spin-resolved sigma-shift state (active when nspin > 1).
    bool spin_resolved_ = false;
    int nspin_ = 1;
    std::vector<int> isk_;
    std::vector<double> n_electrons_per_spin_;
    std::vector<double> lambda_per_spin_;
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
