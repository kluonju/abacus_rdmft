#include "rdmft_orbital_optimizer.h"

#include "rdmft_line_search.h"

#include <algorithm>
#include <cmath>
#include <deque>

namespace rdmft_core
{

namespace
{
//! Limited-memory BFGS history on the packed ambient tangent representation.
//! Directions are re-projected onto the tangent space by the caller.
class PackedLBFGS
{
  public:
    void init(int memory)
    {
        memory_ = std::max(1, memory);
        s_.clear();
        y_.clear();
        rho_.clear();
    }
    void reset()
    {
        s_.clear();
        y_.clear();
        rho_.clear();
    }
    void push(const std::vector<double>& s, const std::vector<double>& y)
    {
        double sy = 0.0;
        for (size_t i = 0; i < s.size(); ++i)
        {
            sy += s[i] * y[i];
        }
        if (sy <= 1.0e-30)
        {
            return; // skip non-curvature pairs
        }
        s_.push_back(s);
        y_.push_back(y);
        rho_.push_back(1.0 / sy);
        if ((int)s_.size() > memory_)
        {
            s_.pop_front();
            y_.pop_front();
            rho_.pop_front();
        }
    }
    //! Two-loop recursion returning d = -H * g.
    void direction(const std::vector<double>& g, std::vector<double>& d) const
    {
        const int n = (int)g.size();
        d.assign(g.begin(), g.end());
        const int m = (int)s_.size();
        std::vector<double> alpha(m, 0.0);
        for (int i = m - 1; i >= 0; --i)
        {
            double a = 0.0;
            for (int k = 0; k < n; ++k)
            {
                a += s_[i][k] * d[k];
            }
            a *= rho_[i];
            alpha[i] = a;
            for (int k = 0; k < n; ++k)
            {
                d[k] -= a * y_[i][k];
            }
        }
        double gamma = 1.0;
        if (m > 0)
        {
            double ys = 0.0;
            double yy = 0.0;
            for (int k = 0; k < n; ++k)
            {
                ys += s_[m - 1][k] * y_[m - 1][k];
                yy += y_[m - 1][k] * y_[m - 1][k];
            }
            if (yy > 1.0e-30)
            {
                gamma = ys / yy;
            }
        }
        for (int k = 0; k < n; ++k)
        {
            d[k] *= gamma;
        }
        for (int i = 0; i < m; ++i)
        {
            double b = 0.0;
            for (int k = 0; k < n; ++k)
            {
                b += y_[i][k] * d[k];
            }
            b *= rho_[i];
            const double coeff = alpha[i] - b;
            for (int k = 0; k < n; ++k)
            {
                d[k] += coeff * s_[i][k];
            }
        }
        for (int k = 0; k < n; ++k)
        {
            d[k] = -d[k];
        }
    }
    bool empty() const { return s_.empty(); }

  private:
    int memory_ = 10;
    std::deque<std::vector<double>> s_;
    std::deque<std::vector<double>> y_;
    std::deque<double> rho_;
};
} // namespace

OrbBlockResult OrbitalOptimizer::run(RdmftBackend& backend, const RdmftParams& params,
                                     const std::vector<double>& occ, double& etot)
{
    OrbBlockResult res;
    if (!backend.has_orbital_optimization() || params.orb_maxiter <= 0)
    {
        res.converged = true;
        res.energy = etot;
        return res;
    }
    const int ndim = backend.orb_dim();
    if (ndim <= 0)
    {
        res.converged = true;
        res.energy = etot;
        return res;
    }

    etot = backend.total_energy(occ);

    std::vector<double> gR(ndim);
    std::vector<double> gR_prev(ndim);
    std::vector<double> xi(ndim);
    std::vector<double> xi_prev(ndim, 0.0);
    bool have_prev = false;

    PackedLBFGS lbfgs;
    if (params.orb_optimizer == OrbOptimizerType::LBFGS)
    {
        lbfgs.init(params.lbfgs_memory);
    }

    for (int inner = 0; inner < params.orb_maxiter; ++inner)
    {
        res.iterations = inner + 1;
        backend.orb_save();
        const double gnorm2 = backend.riemannian_gradient(occ, gR);
        const double gnorm = std::sqrt(std::max(0.0, gnorm2));
        res.grad_norm = gnorm;
        if (gnorm < params.orb_grad_tol)
        {
            res.converged = true;
            break;
        }

        // --- search direction -------------------------------------------------
        if (params.orb_optimizer == OrbOptimizerType::SD || !have_prev)
        {
            for (int i = 0; i < ndim; ++i)
            {
                xi[i] = -gR[i];
            }
        }
        else if (params.orb_optimizer == OrbOptimizerType::CG)
        {
            // Polak-Ribiere with Hager safeguard, product Stiefel metric.
            std::vector<double> diff(ndim);
            for (int i = 0; i < ndim; ++i)
            {
                diff[i] = gR[i] - gR_prev[i];
            }
            const double denom = backend.orb_inner(gR_prev, gR_prev);
            double beta = 0.0;
            if (denom > 1.0e-30)
            {
                beta = std::max(0.0, backend.orb_inner(gR, diff) / denom);
            }
            for (int i = 0; i < ndim; ++i)
            {
                xi[i] = -gR[i] + beta * xi_prev[i];
            }
            backend.orb_project_tangent(xi);
        }
        else // LBFGS
        {
            std::vector<double> d;
            lbfgs.direction(gR, d);
            xi = d;
            backend.orb_project_tangent(xi);
        }

        // Descent safeguard.
        double slope = backend.orb_inner(gR, xi);
        if (slope >= 0.0)
        {
            for (int i = 0; i < ndim; ++i)
            {
                xi[i] = -gR[i];
            }
            if (params.orb_optimizer == OrbOptimizerType::LBFGS)
            {
                lbfgs.reset();
            }
            slope = backend.orb_inner(gR, xi);
            if (slope >= 0.0)
            {
                res.converged = true;
                break;
            }
        }

        // --- line search along the retraction ---------------------------------
        const double E0 = etot;
        const bool need_grad = (params.orb_ls == LineSearchType::StrongWolfe);
        std::vector<double> gR_trial(ndim);
        LineSearchEval eval = [&](double alpha, double& f, double& g, int& ierr) {
            ierr = 0;
            backend.orb_restore();
            backend.orb_retract(xi, alpha);
            f = backend.total_energy(occ);
            if (need_grad)
            {
                backend.riemannian_gradient(occ, gR_trial);
                g = backend.orb_inner(gR_trial, xi);
            }
            else
            {
                g = slope;
            }
        };

        LineSearchResult ls;
        const double alpha0 = 1.0;
        if (params.orb_ls == LineSearchType::StrongWolfe)
        {
            ls = strong_wolfe_line_search(eval, E0, slope, alpha0, params.ls_c1, params.ls_c2,
                                          params.ls_max_iter, params.ls_max_zoom, E0);
        }
        if (!ls.success)
        {
            ls = armijo_line_search(eval, E0, slope, alpha0, params.ls_c1, params.ls_rho,
                                    params.ls_max_iter, E0, params.bb_alpha_min);
        }

        if (!ls.success || ls.step <= 0.0)
        {
            backend.orb_restore();
            etot = backend.total_energy(occ);
            break;
        }

        // Commit the accepted retraction.
        backend.orb_restore();
        backend.orb_retract(xi, ls.step);
        etot = ls.f_new;

        // History updates for CG / L-BFGS.
        std::vector<double> gR_new(ndim);
        backend.riemannian_gradient(occ, gR_new);
        if (params.orb_optimizer == OrbOptimizerType::LBFGS)
        {
            std::vector<double> s(ndim);
            std::vector<double> y(ndim);
            for (int i = 0; i < ndim; ++i)
            {
                s[i] = ls.step * xi[i];
                y[i] = gR_new[i] - gR[i];
            }
            lbfgs.push(s, y);
        }
        gR_prev = gR;
        xi_prev = xi;
        gR = gR_new;
        have_prev = true;

        const double gnorm_new = std::sqrt(std::max(0.0, backend.orb_inner(gR_new, gR_new)));
        res.grad_norm = gnorm_new;
        if (gnorm_new < params.orb_grad_tol)
        {
            res.converged = true;
            break;
        }
    }

    res.energy = etot;
    return res;
}

} // namespace rdmft_core
