#include "rdmft_occ_optimizer.h"

#include "rdmft_line_search.h"

#include <algorithm>
#include <cmath>

namespace rdmft
{

namespace
{
double dot(const std::vector<double>& a, const std::vector<double>& b)
{
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        s += a[i] * b[i];
    }
    return s;
}

double sum_abs_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        s += std::fabs(a[i] - b[i]);
    }
    return s;
}
} // namespace

OccBlockResult OccOptimizer::run(RdmftBackend& backend, const RdmftParams& params,
                                 std::vector<double>& occ, double& etot)
{
    if (params.occ_optimizer == OccOptimizerType::EBI)
    {
        return run_ebi(backend, params, occ, etot);
    }
    return run_spg2(backend, params, occ, etot);
}

// ----------------------------------------------------------------------------
// SPG2 (Birgin-Martinez-Raydan) occupation block
// ----------------------------------------------------------------------------
OccBlockResult OccOptimizer::run_spg2(RdmftBackend& backend, const RdmftParams& params,
                                      std::vector<double>& occ, double& etot)
{
    OccBlockResult res;
    const OccConstraints& con = backend.occ_constraints();
    const int ndim = con.size();
    if (ndim <= 0 || params.occ_maxiter <= 0)
    {
        res.converged = true;
        return res;
    }

    // Ensure feasibility of the starting point.
    con.proximal_project(occ);
    etot = backend.total_energy(occ);

    BarzilaiBorwein bb;
    bb.init(params.bb_alpha_min, params.bb_alpha_max);

    std::vector<double> grad_n(ndim);
    std::vector<double> n0(ndim);
    std::vector<double> dir(ndim);
    std::vector<double> grad_new(ndim);

    bool grad_ready = false;
    double g1_inf = 0.0;

    for (int inner = 0; inner < params.occ_maxiter; ++inner)
    {
        if (!grad_ready)
        {
            backend.grad_occ(occ, grad_n);
            g1_inf = con.pg_map_grad_inf(occ, grad_n, 1.0);
        }
        grad_ready = false;
        res.iterations = inner + 1;
        res.grad_norm = g1_inf;

        if (params.occ_grad_tol > 0.0 && g1_inf <= params.occ_grad_tol)
        {
            res.converged = true;
            break;
        }

        // Spectral steplength: first step alpha0 = 1/||g1||_inf, later steps
        // from the safeguarded inverse Rayleigh quotient on accepted pairs.
        double alpha_init = (g1_inf > 1.0e-30) ? 1.0 / g1_inf : params.bb_alpha_max;
        alpha_init = std::min(params.bb_alpha_max, std::max(params.bb_alpha_min, alpha_init));
        const double alpha_spectral = bb.spectral_step(occ, grad_n, alpha_init);

        n0 = occ;
        const double E0 = etot;

        // Build the SPG chord d = P_w(n0 - alpha*grad_n) - n0.
        std::vector<double> trial(ndim);
        for (int i = 0; i < ndim; ++i)
        {
            trial[i] = n0[i] - alpha_spectral * grad_n[i];
        }
        con.proximal_project(trial);
        for (int i = 0; i < ndim; ++i)
        {
            dir[i] = trial[i] - n0[i];
        }
        const double phi0 = dot(grad_n, dir);
        if (phi0 >= 0.0)
        {
            // Chord uphill (projection clipped the pre-image): stop the block.
            res.converged = (g1_inf <= std::max(params.occ_grad_tol, 1.0e-12));
            break;
        }

        // Straight-chord line search: n(lambda) = n0 + lambda*d, feasible for
        // lambda in [0,1] (convex set), project only if extrapolated.
        const bool need_grad = (params.occ_ls == LineSearchType::StrongWolfe);
        LineSearchEval eval = [&](double alpha, double& f, double& g, int& ierr) {
            ierr = 0;
            std::vector<double> nt(ndim);
            for (int i = 0; i < ndim; ++i)
            {
                nt[i] = n0[i] + alpha * dir[i];
            }
            if (alpha > 1.0)
            {
                con.proximal_project(nt);
                if (con.constraint_violation(nt) > 1.0e-6)
                {
                    ierr = 1;
                    f = 1.0e300;
                    g = 0.0;
                    return;
                }
            }
            f = backend.total_energy(nt);
            if (need_grad)
            {
                std::vector<double> gt(ndim);
                backend.grad_occ(nt, gt);
                g = dot(gt, dir);
            }
            else
            {
                g = phi0;
            }
        };

        LineSearchResult ls;
        if (params.occ_ls == LineSearchType::StrongWolfe)
        {
            ls = strong_wolfe_line_search(eval, E0, phi0, 1.0, params.ls_c1, params.ls_c2,
                                          params.ls_max_iter, params.ls_max_zoom, E0);
        }
        if (!ls.success)
        {
            ls = armijo_line_search(eval, E0, phi0, 1.0, params.ls_c1, params.ls_rho,
                                    params.ls_max_iter, E0, params.bb_alpha_min);
        }

        if (!ls.success || ls.step <= 0.0)
        {
            // Keep current occupations; block stalled.
            occ = n0;
            etot = backend.total_energy(occ);
            break;
        }

        for (int i = 0; i < ndim; ++i)
        {
            occ[i] = n0[i] + ls.step * dir[i];
        }
        if (ls.step > 1.0)
        {
            con.proximal_project(occ);
        }
        etot = ls.f_new;

        // Gradient at the accepted iterate: BB pair + reuse next iteration.
        backend.grad_occ(occ, grad_new);
        bb.record(occ, grad_new);
        grad_n = grad_new;
        g1_inf = con.pg_map_grad_inf(occ, grad_n, 1.0);
        grad_ready = true;
        res.grad_norm = g1_inf;

        if (params.occ_grad_tol > 0.0 && g1_inf <= params.occ_grad_tol)
        {
            res.converged = true;
            break;
        }
    }

    res.energy = etot;
    return res;
}

// ----------------------------------------------------------------------------
// EBI@GD occupation block
// ----------------------------------------------------------------------------
OccBlockResult OccOptimizer::run_ebi(RdmftBackend& backend, const RdmftParams& params,
                                     std::vector<double>& occ, double& etot)
{
    OccBlockResult res;
    const OccConstraints& con = backend.occ_constraints();
    const int ndim = con.size();
    if (ndim <= 0 || params.occ_maxiter <= 0)
    {
        res.converged = true;
        return res;
    }

    std::vector<double> x(ndim);
    con.proximal_project(occ);
    con.ebi_occ_to_params(occ, x);
    con.ebi_fold_mu(x);
    con.ebi_sync_occ_from_x(x, occ);
    etot = backend.total_energy(occ);

    BarzilaiBorwein bb;
    bb.init(params.bb_alpha_min, params.bb_alpha_max);

    std::vector<double> grad_n(ndim);
    std::vector<double> grad_x(ndim);
    std::vector<double> dir(ndim);
    std::vector<double> x0(ndim);
    std::vector<double> n_save(ndim);

    for (int inner = 0; inner < params.occ_maxiter; ++inner)
    {
        res.iterations = inner + 1;
        backend.grad_occ(occ, grad_n);

        double mu_up, mu_dw = 0.0;
        if (con.fix_magnetization)
        {
            mu_up = con.ebi_solve_mu(x, con.n_target_up, 1, true);
            mu_dw = con.ebi_solve_mu(x, con.n_target_down, 2, true);
        }
        else
        {
            mu_up = con.ebi_solve_mu(x, con.n_target, 0, false);
        }
        con.ebi_transform_gradient(grad_n, x, mu_up, mu_dw, grad_x);

        for (int i = 0; i < ndim; ++i)
        {
            dir[i] = -grad_x[i];
        }
        const double gnorm = std::sqrt(dot(grad_x, grad_x));
        res.grad_norm = gnorm;
        if ((params.occ_grad_tol > 0.0 && gnorm <= params.occ_grad_tol) || gnorm < 1.0e-30)
        {
            res.converged = true;
            break;
        }

        x0 = x;
        n_save = occ;
        const double E0 = etot;
        const double phi0 = dot(grad_x, dir); // = -||grad_x||^2 <= 0
        if (phi0 >= 0.0)
        {
            res.converged = true;
            break;
        }

        double alpha0 = bb.spectral_step(x, grad_x, 1.0);
        if (alpha0 <= 0.0)
        {
            alpha0 = 1.0;
        }
        bb.record(x, grad_x);

        const bool need_grad = (params.occ_ls == LineSearchType::StrongWolfe);
        LineSearchEval eval = [&](double alpha, double& f, double& g, int& ierr) {
            ierr = 0;
            std::vector<double> xt(ndim);
            for (int i = 0; i < ndim; ++i)
            {
                xt[i] = x0[i] + alpha * dir[i];
            }
            std::vector<double> nt;
            con.ebi_sync_occ_from_x(xt, nt);
            if (con.constraint_violation(nt) > std::max(10.0 * params.occ_tol, 1.0e-10))
            {
                ierr = 1;
                f = 1.0e300;
                g = 0.0;
                return;
            }
            f = backend.total_energy(nt);
            if (need_grad)
            {
                std::vector<double> gn(ndim);
                backend.grad_occ(nt, gn);
                double mu_u, mu_d = 0.0;
                if (con.fix_magnetization)
                {
                    mu_u = con.ebi_solve_mu(xt, con.n_target_up, 1, true);
                    mu_d = con.ebi_solve_mu(xt, con.n_target_down, 2, true);
                }
                else
                {
                    mu_u = con.ebi_solve_mu(xt, con.n_target, 0, false);
                }
                std::vector<double> gx(ndim);
                con.ebi_transform_gradient(gn, xt, mu_u, mu_d, gx);
                g = dot(gx, dir);
            }
            else
            {
                g = phi0;
            }
        };

        LineSearchResult ls;
        if (params.occ_ls == LineSearchType::StrongWolfe)
        {
            ls = strong_wolfe_line_search(eval, E0, phi0, alpha0, params.ls_c1, params.ls_c2,
                                          params.ls_max_iter, params.ls_max_zoom, E0);
        }
        if (!ls.success)
        {
            ls = armijo_line_search(eval, E0, phi0, alpha0, params.ls_c1, params.ls_rho,
                                    params.ls_max_iter, E0, params.bb_alpha_min);
        }

        if (!ls.success || ls.step <= 0.0)
        {
            x = x0;
            con.ebi_sync_occ_from_x(x, occ);
            etot = backend.total_energy(occ);
            break;
        }

        for (int i = 0; i < ndim; ++i)
        {
            x[i] = x0[i] + ls.step * dir[i];
        }
        con.ebi_fold_mu(x);
        con.ebi_sync_occ_from_x(x, occ);
        etot = backend.total_energy(occ);

        const double sum_dn = sum_abs_diff(occ, n_save);
        if (params.occ_tol > 0.0 && sum_dn < params.occ_tol)
        {
            res.converged = true;
            break;
        }
    }

    res.energy = etot;
    return res;
}

} // namespace rdmft
