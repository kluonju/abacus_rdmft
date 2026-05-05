#ifndef RDMFT_OPTIMIZER_H
#define RDMFT_OPTIMIZER_H

#include "rdmft_type.h"
#include <vector>
#include <functional>
#include <utility>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cassert>

namespace rdmft
{

/// Result of a single optimization run
struct OptResult
{
    int iterations = 0;
    double final_energy = 0.0;
    double grad_norm = 0.0;
    bool converged = false;
};

/// Result of a line search (Strong Wolfe, etc.)
struct LineSearchResult
{
    /// Accepted or last-tried step size
    double step = 0.0;
    /// Objective at `step` (or last energy evaluation on failure)
    double f_new = 0.0;
    bool success = false;
    /// Initial trial step passed into the line search
    double alpha_init = 0.0;
    /// Number of scalar objective evaluations `phi(alpha)` plus derivative
    /// evaluations `phi'(alpha)` (each counted once).
    int n_feval = 0;
};

/// Standalone Barzilai-Borwein step-size estimator.
///
/// Usage pattern:
/// 1) set_mode()/set_bounds() once.
/// 2) each iteration, call suggest(current_x, current_grad, fallback_alpha).
/// 3) after an accepted iterate is committed, call record_state(new_x, new_grad).
class BarzilaiBorweinStep
{
  public:
    BarzilaiBorweinStep() = default;

    void set_mode(BBStepMode mode)
    {
        mode_ = mode;
    }

    void set_bounds(double alpha_min, double alpha_max)
    {
        alpha_min_ = std::max(1e-16, alpha_min);
        alpha_max_ = std::max(alpha_min_, alpha_max);
    }

    void reset()
    {
        have_prev_ = false;
        next_use_bb1_ = true;
        prev_x_.clear();
        prev_g_.clear();
    }

    void record_state(std::vector<double> x,
                      std::vector<double> grad)
    {
        prev_x_ = std::move(x);
        prev_g_ = std::move(grad);
        have_prev_ = true;
    }

    double suggest(const std::vector<double>& x,
                   const std::vector<double>& grad,
                   const double fallback_alpha)
    {
        if (!have_prev_ || prev_x_.size() != x.size() || prev_g_.size() != grad.size())
        {
            return clamp_alpha(fallback_alpha);
        }

        double ss = 0.0;
        double sy = 0.0;
        double yy = 0.0;
        for (size_t i = 0; i < x.size(); ++i)
        {
            const double s = x[i] - prev_x_[i];
            const double y = grad[i] - prev_g_[i];
            ss += s * s;
            sy += s * y;
            yy += y * y;
        }

        const auto valid = [](double v) {
            return std::isfinite(v) && v > 0.0;
        };

        double alpha = fallback_alpha;
        if (mode_ == BBStepMode::BB1)
        {
            if (valid(ss) && valid(sy)) alpha = ss / sy;
        }
        else if (mode_ == BBStepMode::BB2)
        {
            if (valid(sy) && valid(yy)) alpha = sy / yy;
        }
        else
        {
            const bool use_bb1 = next_use_bb1_;
            next_use_bb1_ = !next_use_bb1_;
            if (use_bb1)
            {
                if (valid(ss) && valid(sy)) alpha = ss / sy;
            }
            else
            {
                if (valid(sy) && valid(yy)) alpha = sy / yy;
            }
        }

        return clamp_alpha(alpha);
    }

  private:
    double clamp_alpha(double alpha) const
    {
        if (!std::isfinite(alpha) || alpha <= 0.0)
        {
            alpha = alpha_min_;
        }
        return std::min(alpha_max_, std::max(alpha_min_, alpha));
    }

    BBStepMode mode_ = BBStepMode::Alternate;
    double alpha_min_ = 1e-8;
    double alpha_max_ = 10.0;
    bool have_prev_ = false;
    bool next_use_bb1_ = true;
    std::vector<double> prev_x_;
    std::vector<double> prev_g_;
};

namespace detail
{

/// Quadratic \f$p(t) = f_0 + t\phi'_0 + c t^2\f$ with \f$p(a)=f_a\f$. Returns
/// the unconstrained critical point, or quiet NaN if it is ill-defined.
inline double quadratic_unconstrained_min_t(double f0, double dphi0, double a, double fa)
{
    if (!(std::isfinite(a) && a > 0.0) || !std::isfinite(f0) || !std::isfinite(fa) || !std::isfinite(dphi0))
        return std::numeric_limits<double>::quiet_NaN();
    const double rhs = fa - f0 - a * dphi0;
    const double denom = 2.0 * rhs;
    if (std::abs(denom) < 1.0e-20 * (1.0 + std::max(std::abs(f0), std::abs(fa))))
        return std::numeric_limits<double>::quiet_NaN();
    // p(t) = f0 + t dphi0 + t^2 * rhs / a^2,  p'(t) = 0  =>  t* = -dphi0 a^2 /(2*rhs)
    const double t = -dphi0 * a * a / denom;
    if (!std::isfinite(t))
        return std::numeric_limits<double>::quiet_NaN();
    return t;
}

/// Cubic \f$p(t) = f_0 + t\phi'_0 + c_2 t^2 + c_3 t^3\f$ through \f$(a_1,f_1)\f$, \f$(a_2,f_2)\f$.
/// Returns a local minimizer in \f$(0, \text{step\_cap})\f$ with \f$p''>0\f$, or quiet NaN.
inline double cubic_local_min_t(
    double f0, double dphi0, double a1, double f1, double a2, double f2, double step_cap)
{
    if (!(std::isfinite(a1) && a1 > 0.0) || !(std::isfinite(a2) && a2 > 0.0) || !std::isfinite(step_cap)
        || step_cap <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    if (std::abs(a1 - a2) < 1.0e-20 * (1.0 + std::abs(a1) + std::abs(a2)))
        return std::numeric_limits<double>::quiet_NaN();
    const double r1 = f1 - f0 - a1 * dphi0;
    const double r2 = f2 - f0 - a2 * dphi0;
    // [a1^2  a1^3][c2]   [r1]     D = a1^2 a2^3 - a1^3 a2^2
    // [a2^2  a2^3][c3] = [r2]
    const double D = a1 * a1 * a2 * a2 * a2 - a1 * a1 * a1 * a2 * a2;
    if (std::abs(D) < 1.0e-20)
        return std::numeric_limits<double>::quiet_NaN();
    const double c2 = (r1 * a2 * a2 * a2 - r2 * a1 * a1 * a1) / D;
    const double c3 = (a1 * a1 * r2 - a2 * a2 * r1) / D;
    if (!std::isfinite(c2) || !std::isfinite(c3))
        return std::numeric_limits<double>::quiet_NaN();

    // p'(t) = dphi0 + 2 c2 t + 3 c3 t^2,  p''(t) = 2 c2 + 6 c3 t
    const auto p = [&](double t) { return f0 + t * dphi0 + t * t * c2 + t * t * t * c3; };
    const auto pp2 = [&](double t) { return 2.0 * c2 + 6.0 * c3 * t; };

    const double t_small = 1.0e-15;
    const double hi = std::min(step_cap, std::min(a1, a2)) * (1.0 - 1.0e-10);

    double t_best = std::numeric_limits<double>::quiet_NaN();
    double f_best = std::numeric_limits<double>::infinity();

    if (std::abs(c3) < 1.0e-20 * (1.0 + std::abs(c2)) + 1.0e-20)
    {
        if (std::abs(c2) < 1.0e-20)
            return std::numeric_limits<double>::quiet_NaN();
        const double t = -dphi0 / (2.0 * c2);
        if (t > t_small && t < hi && pp2(t) > 0.0 && p(t) < f_best)
        {
            t_best = t;
            f_best = p(t);
        }
    }
    else
    {
        const double disc = 4.0 * c2 * c2 - 12.0 * c3 * dphi0;
        if (disc >= 0.0)
        {
            const double s = std::sqrt(disc);
            const double t_a = (-2.0 * c2 + s) / (6.0 * c3);
            const double t_b = (-2.0 * c2 - s) / (6.0 * c3);
            for (double t : {t_a, t_b})
            {
                if (t > t_small && t < hi && pp2(t) > 0.0 && p(t) < f_best)
                {
                    t_best = t;
                    f_best = p(t);
                }
            }
        }
    }

    if (!std::isfinite(f_best) || f_best == std::numeric_limits<double>::infinity())
        return std::numeric_limits<double>::quiet_NaN();
    return t_best;
}

} // namespace detail

/// Strong Wolfe line search (Nocedal & Wright, Algorithms 3.5–3.6),
/// with sufficient decrease
/// \f$\varphi(\alpha)\le f_{\text{ref}}+c_1\alpha\varphi'(0)\f$.
/// Pass \f$f_{\text{ref}}=\varphi(0)\f$ for the standard monotone Armijo
/// condition. Requires \f$\varphi'(0)<0\f$, \f$0<c_1<c_2<1\f$.
///
/// `phi` evaluates \f$\varphi(\alpha)\f$; `deriv_phi` evaluates
/// \f$\varphi'(\alpha)\f$ (directional derivative along the search ray).
inline LineSearchResult nonmonotone_strong_wolfe_line_search(
    const std::function<double(double)>& phi,
    const std::function<double(double)>& deriv_phi,
    double phi0,
    double derphi0,
    double f_ref_nm,
    double alpha_init,
    double c1,
    double c2,
    int max_bracket_iter,
    int max_zoom_iter,
    double alpha_expand_cap = 1e10)
{
    LineSearchResult result;
    result.alpha_init = alpha_init;
    result.f_new = phi0;
    result.step = 0.0;
    result.success = false;
    result.n_feval = 0;

    if (!(derphi0 < 0.0) || !std::isfinite(derphi0) || !std::isfinite(phi0) || !std::isfinite(f_ref_nm)
        || !(alpha_init > 0.0) || !(c1 > 0.0 && c1 < 1.0) || !(c2 > c1 && c2 < 1.0)
        || max_bracket_iter <= 0 || max_zoom_iter <= 0)
    {
        return result;
    }

    const auto armijo_nm = [&](double a, double fa) {
        return fa <= f_ref_nm + c1 * a * derphi0;
    };
    const auto curvature_ok = [&](double gp) {
        return std::abs(gp) <= -c2 * derphi0;
    };

    // Precondition: alpha_lo < alpha_hi and phi_lo_known = phi(alpha_lo).
    auto zoom = [&](double alpha_lo,
                    double alpha_hi,
                    double phi_lo_known,
                    double& alpha_out,
                    double& f_out,
                    bool& ok_zoom) {
        ok_zoom = false;
        alpha_out = alpha_lo;
        f_out = phi_lo_known;
        double lo = alpha_lo;
        double hi = alpha_hi;
        double flo = phi_lo_known;
        if (!(lo < hi) || !std::isfinite(lo) || !std::isfinite(hi) || !std::isfinite(flo))
        {
            return;
        }
        constexpr double tol = 1e-12;
        for (int z = 0; z < max_zoom_iter; ++z)
        {
            if (hi - lo <= tol * (1.0 + std::abs(hi)))
            {
                break;
            }
            const double alpha_j = 0.5 * (lo + hi);
            const double phi_j = phi(alpha_j);
            ++result.n_feval;
            if (!std::isfinite(phi_j))
            {
                hi = alpha_j;
                continue;
            }
            if (!armijo_nm(alpha_j, phi_j) || phi_j >= flo)
            {
                hi = alpha_j;
                continue;
            }
            const double der_j = deriv_phi(alpha_j);
            ++result.n_feval;
            if (!std::isfinite(der_j))
            {
                hi = alpha_j;
                continue;
            }
            if (curvature_ok(der_j))
            {
                alpha_out = alpha_j;
                f_out = phi_j;
                ok_zoom = true;
                return;
            }
            if (der_j * (hi - lo) >= 0.0)
            {
                hi = lo;
            }
            lo = alpha_j;
            flo = phi_j;
        }
    };

    double alpha_prev = 0.0;
    double phi_prev = phi0;
    double alpha = alpha_init;

    for (int i = 0; i < max_bracket_iter; ++i)
    {
        const double phi_a = phi(alpha);
        ++result.n_feval;
        if (!std::isfinite(phi_a))
        {
            alpha *= 0.5;
            if (alpha < 1e-30 * std::max(alpha_init, 1.0))
            {
                return result;
            }
            continue;
        }

        if (!armijo_nm(alpha, phi_a) || (phi_a >= phi_prev && alpha_prev > 0.0))
        {
            double z_alpha = 0.0;
            double z_f = 0.0;
            bool zok = false;
            if (alpha_prev < alpha)
            {
                zoom(alpha_prev, alpha, phi_prev, z_alpha, z_f, zok);
            }
            else if (alpha_prev > alpha)
            {
                zoom(alpha, alpha_prev, phi_a, z_alpha, z_f, zok);
            }
            else
            {
                zok = false;
            }
            result.step = z_alpha;
            result.f_new = z_f;
            result.success = zok;
            return result;
        }

        const double der_a = deriv_phi(alpha);
        ++result.n_feval;
        if (!std::isfinite(der_a))
        {
            alpha *= 0.5;
            if (alpha < 1e-30 * std::max(alpha_init, 1.0))
            {
                return result;
            }
            continue;
        }

        if (curvature_ok(der_a))
        {
            result.step = alpha;
            result.f_new = phi_a;
            result.success = true;
            return result;
        }

        if (der_a >= 0.0)
        {
            double z_alpha = 0.0;
            double z_f = 0.0;
            bool zok = false;
            if (alpha_prev < alpha)
            {
                zoom(alpha_prev, alpha, phi_prev, z_alpha, z_f, zok);
            }
            else if (alpha_prev > alpha)
            {
                zoom(alpha, alpha_prev, phi_a, z_alpha, z_f, zok);
            }
            else
            {
                zok = false;
            }
            result.step = z_alpha;
            result.f_new = z_f;
            result.success = zok;
            return result;
        }

        alpha_prev = alpha;
        phi_prev = phi_a;
        const double alpha_next = std::min(2.0 * alpha, alpha_expand_cap);
        if (alpha_next <= alpha * (1.0 + 1e-14))
        {
            return result;
        }
        alpha = alpha_next;
    }

    result.step = alpha_prev;
    result.f_new = phi_prev;
    result.success = false;
    return result;
}

/// Euclidean optimizer for occupation numbers (after parameterization).
/// Works in unconstrained parameter space.
class EuclideanOptimizer
{
  public:
    explicit EuclideanOptimizer(OptimizerType type, const RDMFTConfig& config = RDMFTConfig())
        : type_(type), config_(config) {}

    OptimizerType type() const { return type_; }

    /// Initialize the optimizer state for n_params parameters
    void init(int n_params)
    {
        n_ = n_params;
        step_ = 0;
        prev_grad_.assign(n_, 0.0);
        prev_dir_.assign(n_, 0.0);
    }

    /// Compute search direction from current gradient.
    /// Returns the direction d (caller should do x += step * d).
    void compute_direction(const std::vector<double>& grad, std::vector<double>& dir)
    {
        dir.resize(n_);
        switch (type_)
        {
            case OptimizerType::SteepestDescent:
                for (int i = 0; i < n_; ++i)
                    dir[i] = -grad[i];
                break;

            case OptimizerType::ConjugateGradient:
                compute_cg_direction(grad, dir);
                break;
        }
        prev_grad_ = grad;
    }

    /// Update state after a step (advances CG counters and stores the new gradient).
    void update(const std::vector<double>& new_grad, const std::vector<double>& /*step_vec*/)
    {
        prev_grad_ = new_grad;
        step_++;
    }

  private:
    void compute_cg_direction(const std::vector<double>& grad, std::vector<double>& dir)
    {
        if (step_ == 0)
        {
            for (int i = 0; i < n_; ++i)
                dir[i] = -grad[i];
        }
        else
        {
            // Polak-Ribière with restart
            double num = 0.0, den = 0.0;
            for (int i = 0; i < n_; ++i)
            {
                num += grad[i] * (grad[i] - prev_grad_[i]);
                den += prev_grad_[i] * prev_grad_[i];
            }
            double beta = (den > 1e-30) ? std::max(0.0, num / den) : 0.0;
            for (int i = 0; i < n_; ++i)
                dir[i] = -grad[i] + beta * prev_dir_[i];

            // Check descent condition
            double dd = 0.0;
            for (int i = 0; i < n_; ++i)
                dd += dir[i] * grad[i];
            if (dd >= 0)
                for (int i = 0; i < n_; ++i)
                    dir[i] = -grad[i];
        }
        prev_dir_ = dir;
    }

    OptimizerType type_;
    RDMFTConfig config_;
    int n_ = 0;
    int step_ = 0;
    std::vector<double> prev_grad_;
    std::vector<double> prev_dir_;
};

} // namespace rdmft

#endif // RDMFT_OPTIMIZER_H
