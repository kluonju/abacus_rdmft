//==========================================================
// Author: Kai Luo
// Email: kluo@njust.edu.cn
// DATE : April 2026
//==========================================================
#ifndef RDMFT_OPTIMIZER_H
#define RDMFT_OPTIMIZER_H

#include "rdmft_type.h"
#include <vector>
#include <deque>
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

/// Result of Armijo backtracking (and similar) line search
struct LineSearchResult
{
    /// Accepted or last-tried step size
    double step = 0.0;
    /// Objective at `step` (or last energy evaluation on failure)
    double f_new = 0.0;
    bool success = false;
    /// Initial trial step (the `alpha_init` argument to `armijo_line_search`)
    double alpha_init = 0.0;
    /// Number of objective evaluations in `f_at_step` (one per trial, plus
    /// one on the failure finalisation path when `success` is false).
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

/// Armijo backtracking line search. When `use_polynomial` is true, a failed
/// trial is followed by: (1) a quadratic model using \f$(0,f_0,\phi'_0)\f$ and
/// the first failed \f$(\alpha,f(\alpha))\f$, and (2) on later failures, a
/// cubic through \f$(0,f_0,\phi'_0)\f$ and the last two \f$(\alpha,f(\alpha))\f$
/// pairs. The polynomial suggestion is then clamped to \f$[0.1,0.5]\f$ times the
/// last failed step length (Nocedal & Wright style), then capped by
/// \f$0.99\,\alpha_{\text{fail}}\f$. Otherwise use pure geometric backtracking:
/// multiply by `rho` only.
inline LineSearchResult armijo_line_search(
    std::function<double(double step)> f_at_step,
    double f0,
    double directional_deriv,
    double alpha_init = 0.1,
    double c1 = 1e-4,
    double rho = 0.5,
    int max_iter = 30,
    bool use_polynomial = true)
{
    // Safeguard multipliers for the polynomial trial (fixed; not user-tunable).
    constexpr double k_poly_clamp_lo = 0.1;
    constexpr double k_poly_clamp_hi = 0.5;
    LineSearchResult result;
    result.alpha_init = alpha_init;
    double alpha = alpha_init;
    double prev_fail_alpha = -1.0;
    double prev_fail_f = 0.0;

    for (int i = 0; i < max_iter; ++i)
    {
        const double fail_alpha = alpha;
        const double f_new = f_at_step(fail_alpha);
        ++result.n_feval;
        if (f_new <= f0 + c1 * fail_alpha * directional_deriv)
        {
            result.step = fail_alpha;
            result.f_new = f_new;
            result.success = true;
            return result;
        }

        double alpha_next = 0.0;
        if (use_polynomial)
        {
            double t_s = std::numeric_limits<double>::quiet_NaN();
            if (prev_fail_alpha > 0.0 && std::isfinite(prev_fail_f)
                && std::abs(fail_alpha - prev_fail_alpha)
                       > 1.0e-20 * (1.0 + std::abs(fail_alpha) + std::abs(prev_fail_alpha)))
            {
                const double cap = std::min(fail_alpha, prev_fail_alpha);
                t_s = detail::cubic_local_min_t(
                    f0, directional_deriv, prev_fail_alpha, prev_fail_f, fail_alpha, f_new, cap);
            }
            else
            {
                t_s = detail::quadratic_unconstrained_min_t(
                    f0, directional_deriv, fail_alpha, f_new);
            }
            if (!std::isfinite(t_s) || t_s <= 0.0)
            {
                alpha_next = rho * fail_alpha;
            }
            else
            {
                const double lo = k_poly_clamp_lo * fail_alpha;
                const double hi = k_poly_clamp_hi * fail_alpha;
                if (lo <= hi)
                    t_s = std::min(hi, std::max(lo, t_s));
                else
                    t_s = std::min(lo, std::max(hi, t_s));
                t_s = std::min(t_s, 0.99 * fail_alpha);
                if (t_s <= 0.0 || t_s >= fail_alpha * (1.0 - 1.0e-10))
                    alpha_next = rho * fail_alpha;
                else
                    alpha_next = t_s;
            }
        }
        else
        {
            alpha_next = rho * fail_alpha;
        }

        prev_fail_alpha = fail_alpha;
        prev_fail_f = f_new;
        alpha = alpha_next;
        const double alpha_floor = 1.0e-16 * std::max(alpha_init, 1.0);
        if (alpha < alpha_floor)
        {
            alpha = alpha_floor;
        }
    }

    result.step = alpha;
    result.f_new = f_at_step(alpha);
    ++result.n_feval;
    result.success = false;
    return result;
}

/// Strong Wolfe line search (Nocedal & Wright, Algorithm 3.5 + 3.6).
///
/// `phi(alpha)` returns `{ f(alpha), g(alpha) }` where `g(alpha) = f'(0)` along
/// the line at `x + alpha d` = ((grad f at trial) . direction).
/// `f0`, `g0` are the values at `alpha = 0` (`g0` must be < 0 for a descent
/// direction). Returns a step satisfying strong Wolfe, or a best-effort
/// `success == false` when the budget is exceeded.
inline LineSearchResult strong_wolfe_line_search(
    std::function<std::pair<double, double>(double)> phi,
    double f0,
    double g0,
    double alpha_init = 1.0,
    double c1 = 1e-4,
    double c2 = 0.9,
    int max_iter = 20,
    int max_zoom = 20)
{
    LineSearchResult result;
    result.alpha_init = alpha_init;
    int n_feval = 0;
    auto phi_wrap = [&](double a) -> std::pair<double, double> {
        ++n_feval;
        return phi(a);
    };

    if (!(g0 < 0.0))
    {
        result.step = 0.0;
        result.f_new = f0;
        result.n_feval = 0;
        result.success = false;
        return result;
    }

    // Cubic interpolation: bracket (a, fa, ga) and (b, fb, gb) -> local minimizer.
    auto cubic_min = [](double a, double fa, double ga, double b, double fb, double gb) -> double
    {
        if (std::abs(b - a) <= std::numeric_limits<double>::epsilon()
            * (std::abs(a) + std::abs(b) + 1.0))
        {
            return 0.5 * (a + b);
        }
        const double d1 = ga + gb - 3.0 * (fb - fa) / (b - a);
        const double d2_sq = d1 * d1 - ga * gb;
        if (d2_sq < 0.0)
        {
            return 0.5 * (a + b);
        }
        const double d2 = std::sqrt(d2_sq);
        const double alpha_star = b - (b - a) * (gb + d2 - d1) / (gb - ga + 2.0 * d2);
        const double lo = std::min(a, b);
        const double hi = std::max(a, b);
        const double margin = 0.1 * (hi - lo);
        return std::min(std::max(alpha_star, lo + margin), hi - margin);
    };

    auto zoom = [&](double alpha_lo, double f_lo, double g_lo, double alpha_hi, double f_hi, double g_hi) -> LineSearchResult
    {
        LineSearchResult z_result;
        for (int j = 0; j < max_zoom; ++j)
        {
            const double alpha_j = cubic_min(alpha_lo, f_lo, g_lo, alpha_hi, f_hi, g_hi);

            std::pair<double, double> fg_j = phi_wrap(alpha_j);
            double f_j = fg_j.first;
            double g_j = fg_j.second;

            if (f_j > f0 + c1 * alpha_j * g0 || f_j >= f_lo)
            {
                alpha_hi = alpha_j;
                f_hi = f_j;
                g_hi = g_j;
            }
            else
            {
                if (std::abs(g_j) <= c2 * std::abs(g0))
                {
                    z_result.step = alpha_j;
                    z_result.f_new = f_j;
                    z_result.n_feval = n_feval;
                    z_result.alpha_init = alpha_init;
                    z_result.success = true;
                    return z_result;
                }
                if (g_j * (alpha_hi - alpha_lo) >= 0.0)
                {
                    alpha_hi = alpha_lo;
                    f_hi = f_lo;
                    g_hi = g_lo;
                }
                alpha_lo = alpha_j;
                f_lo = f_j;
                g_lo = g_j;
            }
        }
        z_result.step = alpha_lo;
        z_result.f_new = f_lo;
        z_result.n_feval = n_feval;
        z_result.alpha_init = alpha_init;
        z_result.success = false;
        return z_result;
    };

    const double alpha_max = alpha_init * 100.0;
    double alpha_prev = 0.0;
    double f_prev = f0;
    double g_prev = g0;
    double alpha = alpha_init;

    for (int i = 0; i < max_iter; ++i)
    {
        std::pair<double, double> fg_i = phi_wrap(alpha);
        double f_i = fg_i.first;
        double g_i = fg_i.second;

        if (f_i > f0 + c1 * alpha * g0 || (i > 0 && f_i >= f_prev))
        {
            LineSearchResult zr
                = zoom(alpha_prev, f_prev, g_prev, alpha, f_i, g_i);
            return zr;
        }

        if (std::abs(g_i) <= c2 * std::abs(g0))
        {
            result.step = alpha;
            result.f_new = f_i;
            result.n_feval = n_feval;
            result.success = true;
            return result;
        }

        if (g_i >= 0.0)
        {
            return zoom(alpha, f_i, g_i, alpha_prev, f_prev, g_prev);
        }

        double alpha_new = std::min(2.0 * alpha, alpha_max);
        if (alpha_new <= alpha)
        {
            result.step = alpha;
            result.f_new = f_i;
            result.n_feval = n_feval;
            result.success = false;
            return result;
        }
        alpha_prev = alpha;
        f_prev = f_i;
        g_prev = g_i;
        alpha = alpha_new;
    }

    result.step = alpha_prev;
    result.f_new = f_prev;
    result.n_feval = n_feval;
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

        if (type_ == OptimizerType::Adam)
        {
            m_.assign(n_, 0.0);
            v_.assign(n_, 0.0);
        }
        if (type_ == OptimizerType::LBFGS)
        {
            s_history_.clear();
            y_history_.clear();
        }
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

            case OptimizerType::LBFGS:
                compute_lbfgs_direction(grad, dir);
                break;

            case OptimizerType::Adam:
                compute_adam_direction(grad, dir);
                break;
        }
        prev_grad_ = grad;
    }

    /// Update state after a step (for lbfgs history, etc.)
    void update(const std::vector<double>& new_grad, const std::vector<double>& step_vec)
    {
        if (type_ == OptimizerType::LBFGS)
        {
            std::vector<double> y(n_);
            for (int i = 0; i < n_; ++i)
                y[i] = new_grad[i] - prev_grad_[i];

            double sy = 0.0;
            for (int i = 0; i < n_; ++i)
                sy += step_vec[i] * y[i];

            if (sy > 1e-12)
            {
                s_history_.push_back(step_vec);
                y_history_.push_back(y);
                if (static_cast<int>(s_history_.size()) > config_.lbfgs_memory)
                {
                    s_history_.pop_front();
                    y_history_.pop_front();
                }
            }
        }
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

    void compute_lbfgs_direction(const std::vector<double>& grad, std::vector<double>& dir)
    {
        dir = grad; // q = grad
        int m = s_history_.size();

        std::vector<double> alpha_hist(m);
        std::vector<double> rho_hist(m);

        // First loop
        for (int i = m - 1; i >= 0; --i)
        {
            double sy = 0.0, yy = 0.0;
            for (int j = 0; j < n_; ++j)
            {
                sy += s_history_[i][j] * y_history_[i][j];
                yy += y_history_[i][j] * y_history_[i][j];
            }
            rho_hist[i] = (sy > 1e-30) ? 1.0 / sy : 0.0;

            double a = 0.0;
            for (int j = 0; j < n_; ++j)
                a += s_history_[i][j] * dir[j];
            alpha_hist[i] = rho_hist[i] * a;

            for (int j = 0; j < n_; ++j)
                dir[j] -= alpha_hist[i] * y_history_[i][j];
        }

        // Scale by gamma = s_k^T y_k / y_k^T y_k
        if (m > 0)
        {
            double sy = 0.0, yy = 0.0;
            for (int j = 0; j < n_; ++j)
            {
                sy += s_history_.back()[j] * y_history_.back()[j];
                yy += y_history_.back()[j] * y_history_.back()[j];
            }
            double gamma = (yy > 1e-30) ? sy / yy : 1.0;
            for (int j = 0; j < n_; ++j)
                dir[j] *= gamma;
        }

        // Second loop
        for (int i = 0; i < m; ++i)
        {
            double b = 0.0;
            for (int j = 0; j < n_; ++j)
                b += y_history_[i][j] * dir[j];
            b *= rho_hist[i];
            for (int j = 0; j < n_; ++j)
                dir[j] += (alpha_hist[i] - b) * s_history_[i][j];
        }

        // Negate for descent
        for (int j = 0; j < n_; ++j)
            dir[j] = -dir[j];
    }

    void compute_adam_direction(const std::vector<double>& grad, std::vector<double>& dir)
    {
        step_++;
        double beta1 = config_.adam_beta1;
        double beta2 = config_.adam_beta2;
        double eps = config_.adam_eps;

        for (int i = 0; i < n_; ++i)
        {
            m_[i] = beta1 * m_[i] + (1.0 - beta1) * grad[i];
            v_[i] = beta2 * v_[i] + (1.0 - beta2) * grad[i] * grad[i];
        }

        double bc1 = 1.0 - std::pow(beta1, step_);
        double bc2 = 1.0 - std::pow(beta2, step_);

        for (int i = 0; i < n_; ++i)
        {
            double m_hat = m_[i] / bc1;
            double v_hat = v_[i] / bc2;
            dir[i] = -config_.adam_lr * m_hat / (std::sqrt(v_hat) + eps);
        }
    }

    OptimizerType type_;
    RDMFTConfig config_;
    int n_ = 0;
    int step_ = 0;
    std::vector<double> prev_grad_;
    std::vector<double> prev_dir_;

    // Adam state
    std::vector<double> m_;
    std::vector<double> v_;

    // lbfgs history
    std::deque<std::vector<double>> s_history_;
    std::deque<std::vector<double>> y_history_;
};

} // namespace rdmft

#endif // RDMFT_OPTIMIZER_H
