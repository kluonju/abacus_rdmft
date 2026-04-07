#ifndef RDMFT_OPTIMIZER_H
#define RDMFT_OPTIMIZER_H

#include "rdmft_type.h"
#include <vector>
#include <deque>
#include <functional>
#include <cmath>
#include <algorithm>
#include <numeric>
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

/// Abstract base for line search
struct LineSearchResult
{
    double step = 0.0;
    double f_new = 0.0;
    bool success = false;
};

/// Armijo backtracking line search
/// f_and_grad: evaluates function and gradient at point x + step * d
inline LineSearchResult armijo_line_search(
    std::function<double(double step)> f_at_step,
    double f0,
    double directional_deriv,
    double alpha_init = 0.1,
    double c1 = 1e-4,
    double rho = 0.5,
    int max_iter = 30)
{
    LineSearchResult result;
    double alpha = alpha_init;

    for (int i = 0; i < max_iter; ++i)
    {
        double f_new = f_at_step(alpha);
        if (f_new <= f0 + c1 * alpha * directional_deriv)
        {
            result.step = alpha;
            result.f_new = f_new;
            result.success = true;
            return result;
        }
        alpha *= rho;
    }

    result.step = alpha;
    result.f_new = f_at_step(alpha);
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

    /// Update state after a step (for L-BFGS history, etc.)
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

    // L-BFGS history
    std::deque<std::vector<double>> s_history_;
    std::deque<std::vector<double>> y_history_;
};

} // namespace rdmft

#endif // RDMFT_OPTIMIZER_H
