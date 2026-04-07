#include "rdmft_solver.h"
#include "source_base/timer.h"
#include "source_base/parallel_reduce.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <memory>

namespace rdmft
{

// Helper to get |x|^2 for both real and complex types
inline double abs2(double x) { return x * x; }
inline double abs2(std::complex<double> x) { return std::norm(x); }

template <typename TK, typename TR>
void RDMFTSolver<TK, TR>::init(
    const RDMFTConfig& config,
    EnergyGradient<TK, TR>& energy_grad,
    const K_Vectors* kv,
    int nbands,
    double n_electrons)
{
    config_ = config;
    energy_grad_ = &energy_grad;
    kv_ = kv;
    nbands_ = nbands;
    n_electrons_ = n_electrons;
    nk_ = energy_grad.nk();

    // Initialize occupation parameterization
    occ_param_.reset(new OccupationParam(config_.occ_param));

    // Build k-point weight vector
    std::vector<double> kweights(nk_);
    for (int ik = 0; ik < nk_; ++ik)
        kweights[ik] = kv_->wk[ik];

    occ_constraint_.reset(new OccupationConstraint(
        config_.constraint_method, n_electrons_, kweights, nbands_));

    occ_optimizer_.reset(new EuclideanOptimizer(config_.occ_optimizer, config_));
    orb_optimizer_.reset(new EuclideanOptimizer(config_.orb_optimizer, config_));

    occ_constraint_->set_mu(config_.aug_lag_mu_init);
}

template <typename TK, typename TR>
double RDMFTSolver<TK, TR>::solve(
    std::vector<double>& occ_flat,
    psi::Psi<TK>& wfc)
{
    ModuleBase::timer::tick("RDMFT", "solve");

    double E = 0.0;
    switch (config_.strategy)
    {
        case SolverStrategy::Alternating:
            E = solve_alternating(occ_flat, wfc);
            break;
        case SolverStrategy::ProductManifold:
            E = solve_product_manifold(occ_flat, wfc);
            break;
    }

    ModuleBase::timer::tick("RDMFT", "solve");
    return E;
}

template <typename TK, typename TR>
double RDMFTSolver<TK, TR>::solve_alternating(
    std::vector<double>& occ_flat,
    psi::Psi<TK>& wfc)
{
    std::cout << "\n===== RDMFT Alternating Optimization =====" << std::endl;
    GlobalV::ofs_running << "\n===== RDMFT Alternating Optimization =====" << std::endl;

    double E_prev = 1e30;
    double E = 0.0;

    for (int iter = 0; iter < config_.max_iter; ++iter)
    {
        // 1. Optimize occupations with orbitals fixed
        auto occ_result = optimize_occupations(occ_flat, wfc);
        GlobalV::ofs_running << "    occ inner: " << occ_result.iterations << " iters, gnorm="
            << std::scientific << occ_result.grad_norm
            << "  E=" << std::fixed << std::setprecision(10) << occ_result.final_energy
            << (occ_result.converged ? "  (converged)" : "") << std::endl;

        // 2. Optimize orbitals with occupations fixed
        auto orb_result = optimize_orbitals(occ_flat, wfc);
        GlobalV::ofs_running << "    orb inner: " << orb_result.iterations << " iters, gnorm="
            << std::scientific << orb_result.grad_norm
            << "  E=" << std::fixed << std::setprecision(10) << orb_result.final_energy
            << (orb_result.converged ? "  (converged)" : "") << std::endl;

        // 3. Evaluate full energy
        std::vector<double> grad_occ;
        psi::Psi<TK> grad_wfc;
        E = energy_grad_->compute(occ_flat, wfc, grad_occ, grad_wfc);

        double dE = std::abs(E - E_prev);
        double constraint_viol = occ_constraint_->constraint_violation(occ_flat);

        std::cout << std::fixed << std::setprecision(10)
            << "  RDMFT iter " << iter + 1
            << "  E = " << E
            << "  dE = " << std::scientific << dE
            << "  |c| = " << std::abs(constraint_viol)
            << std::endl;

        GlobalV::ofs_running << std::fixed << std::setprecision(10)
            << "  RDMFT iter " << iter + 1
            << "  E = " << E
            << "  dE = " << std::scientific << dE
            << "  |c| = " << std::abs(constraint_viol)
            << "  occ_gnorm = " << occ_result.grad_norm
            << "  orb_gnorm = " << orb_result.grad_norm
            << std::endl;

        if (dE < config_.energy_tol && std::abs(constraint_viol) < 1e-8)
        {
            last_result_.converged = true;
            last_result_.iterations = iter + 1;
            last_result_.final_energy = E;
            last_result_.grad_norm = std::max(occ_result.grad_norm, orb_result.grad_norm);
            break;
        }

        E_prev = E;
    }

    last_result_.final_energy = E;
    return E;
}

template <typename TK, typename TR>
double RDMFTSolver<TK, TR>::solve_product_manifold(
    std::vector<double>& occ_flat,
    psi::Psi<TK>& wfc)
{
    std::cout << "\n===== RDMFT Product Manifold Optimization =====" << std::endl;
    GlobalV::ofs_running << "\n===== RDMFT Product Manifold Optimization =====" << std::endl;

    // Convert occupations to unconstrained parameters
    std::vector<double> params(occ_flat.size());
    occ_param_->occ_to_params(occ_flat, params);

    int n_occ_params = params.size();
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();

    EuclideanOptimizer prod_occ_opt(config_.occ_optimizer, config_);
    prod_occ_opt.init(n_occ_params);

    double E_prev = 1e30;
    double E = 0.0;

    for (int iter = 0; iter < config_.max_iter; ++iter)
    {
        // Convert params -> occupations
        occ_param_->params_to_occ(params, occ_flat);

        // Compute energy and gradients
        std::vector<double> grad_occ;
        psi::Psi<TK> grad_wfc;
        E = energy_grad_->compute(occ_flat, wfc, grad_occ, grad_wfc);

        // Add augmented Lagrangian penalty
        if (config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
        {
            E += occ_constraint_->augmented_lagrangian_penalty(occ_flat);
            std::vector<double> penalty_grad;
            occ_constraint_->augmented_lagrangian_gradient(occ_flat, penalty_grad);
            for (size_t i = 0; i < grad_occ.size(); ++i)
                grad_occ[i] += penalty_grad[i];
        }

        // Transform occupation gradient to parameter space
        std::vector<double> grad_params;
        occ_param_->transform_gradient_batch(grad_occ, params, grad_params);

        // Compute Riemannian gradient for orbitals (project onto tangent space)
        // For now, we do a simple projected gradient step on orbitals
        // and a Euclidean step on occupation parameters

        // Occupation parameter step
        std::vector<double> occ_dir;
        prod_occ_opt.compute_direction(grad_params, occ_dir);

        // Line search: try different step sizes
        double alpha = config_.line_search_alpha_init;
        double dd = 0.0;
        for (size_t i = 0; i < grad_params.size(); ++i)
            dd += occ_dir[i] * grad_params[i];

        auto f_at_step = [&](double step) -> double {
            std::vector<double> params_trial(params);
            for (size_t i = 0; i < params_trial.size(); ++i)
                params_trial[i] += step * occ_dir[i];
            std::vector<double> occ_trial;
            occ_param_->params_to_occ(params_trial, occ_trial);
            double E_trial = energy_grad_->compute_energy(occ_trial, wfc);
            if (config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
                E_trial += occ_constraint_->augmented_lagrangian_penalty(occ_trial);
            return E_trial;
        };

        auto ls_result = armijo_line_search(f_at_step, E, dd, alpha,
                                             config_.line_search_c1,
                                             config_.line_search_rho,
                                             config_.line_search_max_iter);

        // Update params
        std::vector<double> step_vec(params.size());
        for (size_t i = 0; i < params.size(); ++i)
        {
            step_vec[i] = ls_result.step * occ_dir[i];
            params[i] += step_vec[i];
        }
        prod_occ_opt.update(grad_params, step_vec);

        // Orbital retraction step
        double orb_step = config_.line_search_alpha_init * 0.1;
        const int nk = wfc.get_nk();
        for (int ik = 0; ik < nk; ++ik)
        {
            // Simple gradient descent with reorthogonalization
            TK* C = &wfc(ik, 0, 0);
            const TK* G = &grad_wfc(ik, 0, 0);
            for (int i = 0; i < nb_local * nbs_local; ++i)
                C[i] -= TK(orb_step) * G[i];
        }
        // TODO: proper Stiefel retraction via roptlite

        // Update augmented Lagrangian multiplier periodically
        if (config_.constraint_method == ConstraintMethod::AugmentedLagrangian
            && (iter + 1) % 10 == 0)
        {
            occ_param_->params_to_occ(params, occ_flat);
            occ_constraint_->update_multiplier(occ_flat);
            occ_constraint_->increase_penalty(config_.aug_lag_mu_factor, config_.aug_lag_mu_max);
        }

        double dE = std::abs(E - E_prev);
        double gnorm = 0.0;
        for (auto g : grad_params) gnorm += g * g;
        gnorm = std::sqrt(gnorm);

        std::cout << std::fixed << std::setprecision(10)
            << "  RDMFT prod-iter " << iter + 1
            << "  E = " << E
            << "  dE = " << std::scientific << dE
            << "  gnorm = " << gnorm
            << std::endl;

        GlobalV::ofs_running << std::fixed << std::setprecision(10)
            << "  RDMFT prod-iter " << iter + 1
            << "  E = " << E
            << "  dE = " << std::scientific << dE
            << "  gnorm = " << gnorm
            << std::endl;

        if (dE < config_.energy_tol && gnorm < config_.grad_tol)
        {
            last_result_.converged = true;
            last_result_.iterations = iter + 1;
            last_result_.final_energy = E;
            last_result_.grad_norm = gnorm;
            break;
        }
        E_prev = E;
    }

    // Final conversion
    occ_param_->params_to_occ(params, occ_flat);
    last_result_.final_energy = E;
    return E;
}

template <typename TK, typename TR>
OptResult RDMFTSolver<TK, TR>::optimize_occupations(
    std::vector<double>& occ_flat,
    const psi::Psi<TK>& wfc)
{
    OptResult result;

    switch (config_.constraint_method)
    {
        case ConstraintMethod::AugmentedLagrangian:
        {
            // Convert to unconstrained parameters
            std::vector<double> params;
            occ_param_->occ_to_params(occ_flat, params);

            EuclideanOptimizer opt(config_.occ_optimizer, config_);
            opt.init(params.size());

            for (int inner = 0; inner < config_.max_inner_iter; ++inner)
            {
                occ_param_->params_to_occ(params, occ_flat);

                std::vector<double> grad_occ;
                psi::Psi<TK> grad_wfc_dummy;
                double E = energy_grad_->compute(occ_flat, const_cast<psi::Psi<TK>&>(wfc),
                                                  grad_occ, grad_wfc_dummy);
                E += occ_constraint_->augmented_lagrangian_penalty(occ_flat);

                std::vector<double> penalty_grad;
                occ_constraint_->augmented_lagrangian_gradient(occ_flat, penalty_grad);
                for (size_t i = 0; i < grad_occ.size(); ++i)
                    grad_occ[i] += penalty_grad[i];

                std::vector<double> grad_params;
                occ_param_->transform_gradient_batch(grad_occ, params, grad_params);

                result.grad_norm = 0.0;
                for (auto g : grad_params) result.grad_norm += g * g;
                result.grad_norm = std::sqrt(result.grad_norm);

                GlobalV::ofs_running << "      occ inner " << inner + 1
                    << "  E=" << std::fixed << std::setprecision(10) << E
                    << "  gnorm=" << std::scientific << result.grad_norm << std::endl;

                if (result.grad_norm < config_.grad_tol)
                {
                    result.converged = true;
                    result.iterations = inner + 1;
                    result.final_energy = E;
                    break;
                }

                std::vector<double> dir;
                opt.compute_direction(grad_params, dir);

                double dd = 0.0;
                for (size_t i = 0; i < dir.size(); ++i)
                    dd += dir[i] * grad_params[i];

                auto f_at_step = [&](double step) -> double {
                    std::vector<double> p_trial(params);
                    for (size_t i = 0; i < p_trial.size(); ++i)
                        p_trial[i] += step * dir[i];
                    std::vector<double> occ_trial;
                    occ_param_->params_to_occ(p_trial, occ_trial);
                    double Et = energy_grad_->compute_energy(occ_trial,
                                    const_cast<psi::Psi<TK>&>(wfc));
                    Et += occ_constraint_->augmented_lagrangian_penalty(occ_trial);
                    return Et;
                };

                auto ls = armijo_line_search(f_at_step, E, dd,
                    config_.line_search_alpha_init, config_.line_search_c1,
                    config_.line_search_rho, config_.line_search_max_iter);

                std::vector<double> step_vec(params.size());
                for (size_t i = 0; i < params.size(); ++i)
                {
                    step_vec[i] = ls.step * dir[i];
                    params[i] += step_vec[i];
                }

                occ_param_->params_to_occ(params, occ_flat);
                std::vector<double> new_grad_occ;
                energy_grad_->compute(occ_flat, const_cast<psi::Psi<TK>&>(wfc),
                                       new_grad_occ, grad_wfc_dummy);
                std::vector<double> new_grad_params;
                occ_param_->transform_gradient_batch(new_grad_occ, params, new_grad_params);
                opt.update(new_grad_params, step_vec);

                result.iterations = inner + 1;
                result.final_energy = ls.f_new;
            }

            // Update Lagrange multiplier
            occ_param_->params_to_occ(params, occ_flat);
            occ_constraint_->update_multiplier(occ_flat);
            if (std::abs(occ_constraint_->constraint_violation(occ_flat)) > 1e-6)
                occ_constraint_->increase_penalty(config_.aug_lag_mu_factor, config_.aug_lag_mu_max);
            break;
        }

        case ConstraintMethod::ProjectedGradient:
        {
            for (int inner = 0; inner < config_.max_inner_iter; ++inner)
            {
                std::vector<double> grad_occ;
                psi::Psi<TK> grad_wfc_dummy;
                double E = energy_grad_->compute(occ_flat, const_cast<psi::Psi<TK>&>(wfc),
                                                  grad_occ, grad_wfc_dummy);

                result.grad_norm = 0.0;
                for (auto g : grad_occ) result.grad_norm += g * g;
                result.grad_norm = std::sqrt(result.grad_norm);

                if (result.grad_norm < config_.grad_tol) break;

                double step = config_.line_search_alpha_init;
                for (size_t i = 0; i < occ_flat.size(); ++i)
                    occ_flat[i] -= step * grad_occ[i];

                occ_constraint_->project(occ_flat);
                result.iterations = inner + 1;
                result.final_energy = E;
            }
            break;
        }

        case ConstraintMethod::ActiveSet:
        {
            for (int inner = 0; inner < config_.max_inner_iter; ++inner)
            {
                std::vector<double> grad_occ;
                psi::Psi<TK> grad_wfc_dummy;
                double E = energy_grad_->compute(occ_flat, const_cast<psi::Psi<TK>&>(wfc),
                                                  grad_occ, grad_wfc_dummy);

                auto as_info = occ_constraint_->identify_active_set(occ_flat, grad_occ);
                occ_constraint_->apply_active_set(as_info, grad_occ);

                result.grad_norm = 0.0;
                for (auto g : grad_occ) result.grad_norm += g * g;
                result.grad_norm = std::sqrt(result.grad_norm);

                if (result.grad_norm < config_.grad_tol) break;

                double step = config_.line_search_alpha_init;
                for (size_t i = 0; i < occ_flat.size(); ++i)
                    occ_flat[i] -= step * grad_occ[i];

                // Clip to [0,1]
                for (auto& n : occ_flat)
                    n = std::max(0.0, std::min(1.0, n));

                result.iterations = inner + 1;
                result.final_energy = E;
            }
            break;
        }
    }

    return result;
}

template <typename TK, typename TR>
OptResult RDMFTSolver<TK, TR>::optimize_orbitals(
    const std::vector<double>& occ_flat,
    psi::Psi<TK>& wfc)
{
    OptResult result;

    const int nk = wfc.get_nk();
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();

    energy_grad_->invalidate_hone_cache();

    for (int inner = 0; inner < config_.max_inner_iter; ++inner)
    {
        std::vector<double> grad_occ;
        psi::Psi<TK> grad_wfc;
        double E = energy_grad_->compute(const_cast<std::vector<double>&>(occ_flat),
                                          wfc, grad_occ, grad_wfc);

        result.grad_norm = 0.0;
        for (int ik = 0; ik < nk; ++ik)
            for (int ib = 0; ib < nb_local; ++ib)
                for (int mu = 0; mu < nbs_local; ++mu)
                    result.grad_norm += abs2(grad_wfc(ik, ib, mu));
        result.grad_norm = std::sqrt(result.grad_norm);

        GlobalV::ofs_running << "      orb inner " << inner + 1
            << "  E=" << std::fixed << std::setprecision(10) << E
            << "  gnorm=" << std::scientific << result.grad_norm << std::endl;

        if (result.grad_norm < config_.grad_tol)
        {
            result.converged = true;
            result.iterations = inner + 1;
            result.final_energy = E;
            break;
        }

        double step = config_.line_search_alpha_init * 0.1;
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int ib = 0; ib < nb_local; ++ib)
                for (int mu = 0; mu < nbs_local; ++mu)
                    wfc(ik, ib, mu) -= TK(step) * grad_wfc(ik, ib, mu);
        }

        energy_grad_->invalidate_hone_cache();

        result.iterations = inner + 1;
        result.final_energy = E;
    }

    return result;
}

template <typename TK, typename TR>
bool RDMFTSolver<TK, TR>::check_gradient_consistency(
    const std::vector<double>& occ_flat,
    const psi::Psi<TK>& wfc,
    double epsilon,
    double tolerance)
{
    GlobalV::ofs_running << "\n===== RDMFT Gradient Consistency Check =====" << std::endl;

    // Compute analytic gradients
    std::vector<double> grad_occ;
    psi::Psi<TK> grad_wfc;
    double E0 = energy_grad_->compute(
        const_cast<std::vector<double>&>(occ_flat),
        const_cast<psi::Psi<TK>&>(wfc),
        grad_occ, grad_wfc);

    bool all_pass = true;

    // Check occupation gradients by finite difference
    GlobalV::ofs_running << "\n-- Occupation gradient check --" << std::endl;
    int n_occ_check = std::min(static_cast<int>(occ_flat.size()), 10);

    for (int idx = 0; idx < n_occ_check; ++idx)
    {
        std::vector<double> occ_plus(occ_flat);
        std::vector<double> occ_minus(occ_flat);

        double n_orig = occ_flat[idx];
        if (n_orig + epsilon > 1.0 || n_orig - epsilon < 0.0)
            continue;

        occ_plus[idx] += epsilon;
        occ_minus[idx] -= epsilon;

        double E_plus = energy_grad_->compute_energy(occ_plus,
                            const_cast<psi::Psi<TK>&>(wfc));
        double E_minus = energy_grad_->compute_energy(occ_minus,
                            const_cast<psi::Psi<TK>&>(wfc));

        double fd_grad = (E_plus - E_minus) / (2.0 * epsilon);
        double analytic = grad_occ[idx];
        double rel_err = (std::abs(analytic) > 1e-10)
                         ? std::abs(fd_grad - analytic) / std::abs(analytic)
                         : std::abs(fd_grad - analytic);

        bool pass = rel_err < tolerance;
        if (!pass) all_pass = false;

        GlobalV::ofs_running << std::fixed << std::setprecision(8)
            << "  occ[" << idx << "]: analytic=" << analytic
            << "  fd=" << fd_grad
            << "  rel_err=" << rel_err
            << (pass ? "  PASS" : "  FAIL")
            << std::endl;
    }

    GlobalV::ofs_running << "\nGradient check " << (all_pass ? "PASSED" : "FAILED") << std::endl;
    return all_pass;
}

template <typename TK, typename TR>
double RDMFTSolver<TK, TR>::occ_grad_norm(const std::vector<double>& grad) const
{
    double norm = 0.0;
    for (auto g : grad) norm += g * g;
    return std::sqrt(norm);
}

template <typename TK, typename TR>
double RDMFTSolver<TK, TR>::orb_grad_norm(const psi::Psi<TK>& rgrad) const
{
    double norm = 0.0;
    for (int ik = 0; ik < rgrad.get_nk(); ++ik)
        for (int ib = 0; ib < rgrad.get_nbands(); ++ib)
            for (int mu = 0; mu < rgrad.get_nbasis(); ++mu)
                norm += abs2(rgrad(ik, ib, mu));
    return std::sqrt(norm);
}

// Explicit template instantiations
template class RDMFTSolver<double, double>;
template class RDMFTSolver<std::complex<double>, double>;
template class RDMFTSolver<std::complex<double>, std::complex<double>>;

} // namespace rdmft
