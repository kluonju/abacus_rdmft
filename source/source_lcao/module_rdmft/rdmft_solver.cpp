#include "rdmft_solver.h"
#include "source_base/timer.h"
#include "source_base/parallel_reduce.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <memory>
#include <type_traits>

namespace rdmft
{
namespace
{
double sum_abs_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    assert(a.size() == b.size());
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        s += std::abs(a[i] - b[i]);
    return s;
}

/// Stop occupation inner loop if the parameter-space gradient is small, or if the
/// total change in occupations sum_i |Δn_i| is small. Exact-zero Δn with a large
/// gradient (e.g. failed line search) does not count as converged.
bool occ_inner_should_stop(double sum_abs_dn,
                           double grad_norm,
                           double dn_tol,
                           double occ_grad_tol)
{
    if (grad_norm < occ_grad_tol)
        return true;
    const double abs_floor = 1e-20;
    return (sum_abs_dn < dn_tol && sum_abs_dn > abs_floor);
}

void log_occ_inner_line(const std::string& line)
{
    GlobalV::ofs_running << line << std::endl;
}

/// One summary line, then per-ik lines listing each n(ik,ib) and step change dn from the
/// previous inner-iteration start (occ_prev_for_dn empty => first step, dn = 0).
void log_occ_inner_summary_and_nik(const std::string& summary_first_line,
                                     const std::vector<double>& occ_flat,
                                     const std::vector<double>& occ_prev_for_dn,
                                     const int nk,
                                     const int nbands)
{
    log_occ_inner_line(summary_first_line);
    const bool have_prev = (occ_prev_for_dn.size() == occ_flat.size());
    for (int ik = 0; ik < nk; ++ik)
    {
        std::ostringstream row;
        row << "        ik=" << ik;
        for (int ib = 0; ib < nbands; ++ib)
        {
            const int idx = ik * nbands + ib;
            const double n = occ_flat[idx];
            const double dn = have_prev ? (n - occ_prev_for_dn[idx]) : 0.0;
            row << "  n(" << ik << "," << ib << ")=" << std::fixed << std::setprecision(8) << n
                << " dn=" << dn;
        }
        log_occ_inner_line(row.str());
    }
}
} // namespace

// Helper to get |x|^2 for both real and complex types
inline double abs2(double x) { return x * x; }
inline double abs2(std::complex<double> x) { return std::norm(x); }

// ------------------------------------------------------------------------
// Flatten / unflatten helpers for the Stiefel-manifold orbital optimiser.
// To reuse the Euclidean L-BFGS / Adam optimiser working on a flat
// std::vector<double>, we map psi::Psi<TK> onto a real-valued flat array:
//   * TK = double              -> size = N, 1 double per entry
//   * TK = complex<double>     -> size = 2N, (Re, Im) per entry
// With this layout the Euclidean inner product on the flattened vector
// equals Re Tr(X^H Y) in Psi space, which is the Frobenius inner product
// the EuclideanOptimizer expects (up to overlap S -- see note below).
//
// Note on the Riemannian metric: the Stiefel manifold uses the S-weighted
// inner product  <X, Y>_S = Re Tr(X^H S Y).  The flattened Euclidean dot
// product coincides with the Riemannian one only when S = I (PW basis or
// already S-orthonormalised LCAO columns).  For S != I we still get a
// valid descent direction after projecting back onto the tangent space,
// but the L-BFGS preconditioning is an approximation.  This matches the
// standard "Riemannian L-BFGS by projection" prescription used by the
// ROPTLITE library and is adequate for the current applications.
// ------------------------------------------------------------------------
inline void psi_to_flat(const psi::Psi<double>& P, std::vector<double>& flat)
{
    const int nk = P.get_nk();
    const int nb = P.get_nbands();
    const int nbs = P.get_nbasis();
    const int n = nk * nb * nbs;
    flat.resize(n);
    const double* p = &P(0, 0, 0);
    std::copy(p, p + n, flat.begin());
}

inline void psi_to_flat(const psi::Psi<std::complex<double>>& P,
                        std::vector<double>& flat)
{
    const int nk = P.get_nk();
    const int nb = P.get_nbands();
    const int nbs = P.get_nbasis();
    const int n = nk * nb * nbs;
    flat.resize(2 * n);
    const std::complex<double>* p = &P(0, 0, 0);
    for (int i = 0; i < n; ++i)
    {
        flat[2 * i]     = p[i].real();
        flat[2 * i + 1] = p[i].imag();
    }
}

inline void flat_to_psi(const std::vector<double>& flat, psi::Psi<double>& P)
{
    const int nk = P.get_nk();
    const int nb = P.get_nbands();
    const int nbs = P.get_nbasis();
    const int n = nk * nb * nbs;
    double* p = &P(0, 0, 0);
    for (int i = 0; i < n; ++i) p[i] = flat[i];
}

inline void flat_to_psi(const std::vector<double>& flat,
                        psi::Psi<std::complex<double>>& P)
{
    const int nk = P.get_nk();
    const int nb = P.get_nbands();
    const int nbs = P.get_nbasis();
    const int n = nk * nb * nbs;
    std::complex<double>* p = &P(0, 0, 0);
    for (int i = 0; i < n; ++i)
        p[i] = std::complex<double>(flat[2 * i], flat[2 * i + 1]);
}

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
    ModuleBase::timer::start("RDMFT", "solve");

    // Clamp the initial occupations away from the [0, 1] boundary. With the
    // cosine^2 / logistic parameterisations dn/dp vanishes at n = 0 and n = 1,
    // so a KS seed (which typically has integer occupations) stalls the
    // parameter-space optimiser on iteration 1 even though the analytic dE/dn
    // is non-zero. Clamping by a small margin and rescaling to preserve the
    // electron count restores a non-zero dp-gradient without perturbing the
    // final converged energy (the margin vanishes when the optimum has
    // non-integer occupations).
    const double m = config_.occ_init_margin;
    if (m > 0.0 && m < 0.5)
    {
        // occ_init_nbands_top == 0: all bands (legacy). K > 0: only top K bands per k.
        const int K_cfg = config_.occ_init_nbands_top;
        const int K_eff = (K_cfg == 0) ? nbands_ : std::min(K_cfg, nbands_);
        const int ib_min_target = nbands_ - K_eff;

        for (int ik = 0; ik < nk_; ++ik)
        {
            for (int ib = ib_min_target; ib < nbands_; ++ib)
            {
                double& n = occ_flat[ik * nbands_ + ib];
                n = std::max(m, std::min(1.0 - m, n));
            }
        }

        // Re-scale to restore the electron-number constraint  sum wk n = Ne.
        // Only scale target-band values in (m, 1-m); values at the clamped boundaries are
        // held fixed so they stay feasible under a uniform rescale.
        double current = 0.0;
        double free_sum = 0.0;
        for (int ik = 0; ik < nk_; ++ik)
        {
            const double wk = kv_->wk[ik];
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const double n = occ_flat[ik * nbands_ + ib];
                current += wk * n;
                if (ib >= ib_min_target && n > m && n < 1.0 - m)
                {
                    free_sum += wk * n;
                }
            }
        }
        const double delta = n_electrons_ - current;
        if (std::abs(delta) > 1e-12 && free_sum > 1e-12)
        {
            const double scale = (free_sum + delta) / free_sum;
            for (int ik = 0; ik < nk_; ++ik)
            {
                for (int ib = ib_min_target; ib < nbands_; ++ib)
                {
                    double& n = occ_flat[ik * nbands_ + ib];
                    if (n > m && n < 1.0 - m)
                    {
                        n *= scale;
                        n = std::max(m, std::min(1.0 - m, n));
                    }
                }
            }
        }
    }

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

    ModuleBase::timer::end("RDMFT", "solve");
    return E;
}

template <typename TK, typename TR>
double RDMFTSolver<TK, TR>::solve_alternating(
    std::vector<double>& occ_flat,
    psi::Psi<TK>& wfc)
{
    GlobalV::ofs_running << "\n===== RDMFT Alternating Optimization =====" << std::endl;

    double E_prev = 1e30;
    double E = 0.0;

    for (int iter = 0; iter < config_.outer_maxiter; ++iter)
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

    for (int iter = 0; iter < config_.outer_maxiter; ++iter)
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

        // Project orbital gradient onto the Stiefel tangent space so the orbital
        // step moves within the manifold and converges at the correct critical point.
        energy_grad_->project_orbital_gradient(wfc, grad_wfc);

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

        // Orbital retraction step using the Riemannian gradient
        double orb_step = config_.line_search_alpha_init * 0.1;
        const int nk = wfc.get_nk();
        for (int ik = 0; ik < nk; ++ik)
        {
            TK* C = &wfc(ik, 0, 0);
            const TK* G = &grad_wfc(ik, 0, 0);
            for (int i = 0; i < nb_local * nbs_local; ++i)
                C[i] -= TK(orb_step) * G[i];
        }

        // Update the augmented-Lagrangian multiplier every outer iteration so
        // that lambda tracks the constraint violation closely. This prevents
        // the transient |c| excursions seen when lambda was only refreshed
        // every 10 iterations. The penalty parameter mu is only doubled when
        // the constraint violation is still above a threshold, to avoid the
        // penalty term dominating the gradient and causing erratic steps.
        if (config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
        {
            occ_param_->params_to_occ(params, occ_flat);
            occ_constraint_->update_multiplier(occ_flat);
            if (std::abs(occ_constraint_->constraint_violation(occ_flat)) > 1e-6)
            {
                occ_constraint_->increase_penalty(
                    config_.aug_lag_mu_factor, config_.aug_lag_mu_max);
            }
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

            double L_prev = 0.0;
            bool have_L_prev = false;
            std::vector<double> occ_snap_start;

            for (int inner = 0; inner < config_.occ_maxiter; ++inner)
            {
                occ_param_->params_to_occ(params, occ_flat);

                const std::vector<double> occ_prev_for_dn = occ_snap_start;
                occ_snap_start = occ_flat;
                const std::vector<double> occ_at_step_start(occ_flat);

                std::vector<double> grad_occ;
                psi::Psi<TK> grad_wfc_dummy;
                const double E_phys = energy_grad_->compute(occ_flat, const_cast<psi::Psi<TK>&>(wfc),
                                                              grad_occ, grad_wfc_dummy);
                const double c = occ_constraint_->constraint_violation(occ_flat);
                const double pen = occ_constraint_->augmented_lagrangian_penalty(occ_flat);
                const double L = E_phys + pen;
                const double dL = have_L_prev ? (L - L_prev) : 0.0;
                L_prev = L;
                have_L_prev = true;

                std::vector<double> penalty_grad;
                occ_constraint_->augmented_lagrangian_gradient(occ_flat, penalty_grad);
                for (size_t i = 0; i < grad_occ.size(); ++i)
                    grad_occ[i] += penalty_grad[i];

                std::vector<double> grad_params;
                occ_param_->transform_gradient_batch(grad_occ, params, grad_params);

                result.grad_norm = 0.0;
                for (auto g : grad_params) result.grad_norm += g * g;
                result.grad_norm = std::sqrt(result.grad_norm);

                {
                    std::ostringstream os;
                    os << "      occ inner " << (inner + 1) << "  E_phys=" << std::fixed
                       << std::setprecision(10) << E_phys << "  L_aug=" << L << "  dL=" << std::scientific
                       << dL << "  |c|=" << std::abs(c) << "  gnorm=" << result.grad_norm;
                    log_occ_inner_summary_and_nik(os.str(), occ_flat, occ_prev_for_dn, nk_, nbands_);
                }

                if (result.grad_norm < config_.occ_grad_tol)
                {
                    result.converged = true;
                    result.iterations = inner + 1;
                    result.final_energy = L;
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

                auto ls = armijo_line_search(f_at_step, L, dd,
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

                const double sum_abs_dn = sum_abs_diff(occ_flat, occ_at_step_start);
                GlobalV::ofs_running << "      sum|dn|=" << std::scientific << sum_abs_dn
                    << std::endl;
                if (occ_inner_should_stop(sum_abs_dn,
                        result.grad_norm,
                        config_.occ_dn_sum_tol,
                        config_.occ_grad_tol))
                {
                    result.converged = true;
                    break;
                }
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
            std::vector<double> occ_snap_start;
            double E_prev = 0.0;
            bool have_E_prev = false;

            for (int inner = 0; inner < config_.occ_maxiter; ++inner)
            {
                const std::vector<double> occ_prev_for_dn = occ_snap_start;
                occ_snap_start = occ_flat;

                std::vector<double> grad_occ;
                psi::Psi<TK> grad_wfc_dummy;
                const double E = energy_grad_->compute(occ_flat, const_cast<psi::Psi<TK>&>(wfc),
                                                       grad_occ, grad_wfc_dummy);

                const double c_abs = std::abs(occ_constraint_->constraint_violation(occ_flat));
                const double dE = have_E_prev ? (E - E_prev) : 0.0;
                E_prev = E;
                have_E_prev = true;

                result.grad_norm = 0.0;
                for (auto g : grad_occ) result.grad_norm += g * g;
                result.grad_norm = std::sqrt(result.grad_norm);

                {
                    std::ostringstream os;
                    os << "      occ inner " << (inner + 1) << "  E=" << std::fixed
                       << std::setprecision(10) << E << "  dE=" << std::scientific << dE
                       << "  |c|=" << c_abs << "  gnorm=" << result.grad_norm;
                    log_occ_inner_summary_and_nik(os.str(), occ_flat, occ_prev_for_dn, nk_, nbands_);
                }

                if (result.grad_norm < config_.occ_grad_tol)
                {
                    result.converged = true;
                    result.iterations = inner + 1;
                    result.final_energy = E;
                    break;
                }

                const std::vector<double> occ_before_step(occ_flat);
                double step = config_.line_search_alpha_init;
                for (size_t i = 0; i < occ_flat.size(); ++i)
                    occ_flat[i] -= step * grad_occ[i];

                occ_constraint_->project(occ_flat);
                const double sum_abs_dn = sum_abs_diff(occ_flat, occ_before_step);
                GlobalV::ofs_running << "      sum|dn|=" << std::scientific << sum_abs_dn
                    << std::endl;
                if (occ_inner_should_stop(sum_abs_dn,
                        result.grad_norm,
                        config_.occ_dn_sum_tol,
                        config_.occ_grad_tol))
                {
                    result.converged = true;
                    result.iterations = inner + 1;
                    result.final_energy = E;
                    break;
                }
                result.iterations = inner + 1;
                result.final_energy = E;
            }
            break;
        }

        case ConstraintMethod::ActiveSet:
        {
            std::vector<double> occ_snap_start;
            double E_prev = 0.0;
            bool have_E_prev = false;

            for (int inner = 0; inner < config_.occ_maxiter; ++inner)
            {
                const std::vector<double> occ_prev_for_dn = occ_snap_start;
                occ_snap_start = occ_flat;

                std::vector<double> grad_occ;
                psi::Psi<TK> grad_wfc_dummy;
                const double E = energy_grad_->compute(occ_flat, const_cast<psi::Psi<TK>&>(wfc),
                                                       grad_occ, grad_wfc_dummy);

                auto as_info = occ_constraint_->identify_active_set(occ_flat, grad_occ);
                occ_constraint_->apply_active_set(as_info, grad_occ);

                const double c_abs = std::abs(occ_constraint_->constraint_violation(occ_flat));
                const double dE = have_E_prev ? (E - E_prev) : 0.0;
                E_prev = E;
                have_E_prev = true;

                result.grad_norm = 0.0;
                for (auto g : grad_occ) result.grad_norm += g * g;
                result.grad_norm = std::sqrt(result.grad_norm);

                {
                    std::ostringstream os;
                    os << "      occ inner " << (inner + 1) << "  E=" << std::fixed
                       << std::setprecision(10) << E << "  dE=" << std::scientific << dE
                       << "  |c|=" << c_abs << "  gnorm=" << result.grad_norm;
                    log_occ_inner_summary_and_nik(os.str(), occ_flat, occ_prev_for_dn, nk_, nbands_);
                }

                if (result.grad_norm < config_.occ_grad_tol)
                {
                    result.converged = true;
                    result.iterations = inner + 1;
                    result.final_energy = E;
                    break;
                }

                const std::vector<double> occ_before_step(occ_flat);
                double step = config_.line_search_alpha_init;
                for (size_t i = 0; i < occ_flat.size(); ++i)
                    occ_flat[i] -= step * grad_occ[i];

                // Clip to [0,1]
                for (auto& n : occ_flat)
                    n = std::max(0.0, std::min(1.0, n));

                const double sum_abs_dn = sum_abs_diff(occ_flat, occ_before_step);
                GlobalV::ofs_running << "      sum|dn|=" << std::scientific << sum_abs_dn
                    << std::endl;
                if (occ_inner_should_stop(sum_abs_dn,
                        result.grad_norm,
                        config_.occ_dn_sum_tol,
                        config_.occ_grad_tol))
                {
                    result.converged = true;
                    result.iterations = inner + 1;
                    result.final_energy = E;
                    break;
                }
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
    const int total_size = nk * nb_local * nbs_local;

    const OptimizerType opt_type = config_.orb_optimizer;
    const bool use_cg    = (opt_type == OptimizerType::ConjugateGradient);
    const bool use_lbfgs = (opt_type == OptimizerType::LBFGS);
    const bool use_adam  = (opt_type == OptimizerType::Adam);

    energy_grad_->invalidate_hone_cache();

    // EuclideanOptimizer for L-BFGS / Adam.  SD and CG are handled by the
    // existing bespoke Riemannian code paths below because they need
    // manifold-aware vector transport / restart heuristics.
    EuclideanOptimizer eucl_opt(opt_type, config_);
    int flat_size = 0;
    if (use_lbfgs || use_adam)
    {
        flat_size = total_size;
        if constexpr (!std::is_same<TK, double>::value) flat_size *= 2;
        eucl_opt.init(flat_size);
    }

    // Work buffers: snapshot of the current wfc used for line-search rollback,
    // previous Riemannian gradient and previous search direction (both needed
    // for conjugate-gradient). prev_grad / prev_dir remain uninitialised when
    // use_cg is false.
    psi::Psi<TK> wfc_save(wfc);
    psi::Psi<TK> prev_grad;
    psi::Psi<TK> prev_dir;
    if (use_cg)
    {
        prev_grad = wfc;
        prev_dir = wfc;
        prev_grad.zero_out();
        prev_dir.zero_out();
    }

    double prev_gnorm2 = 0.0;

    for (int inner = 0; inner < config_.orb_maxiter; ++inner)
    {
        std::vector<double> grad_occ;
        psi::Psi<TK> grad_wfc;
        double E = energy_grad_->compute(const_cast<std::vector<double>&>(occ_flat),
                                          wfc, grad_occ, grad_wfc);

        // Project onto the tangent space of the generalised Stiefel manifold
        // with overlap S:  G_R = G - S C sym(C^H G).
        // At any S-orthonormal critical point, G_R == 0.
        energy_grad_->project_orbital_gradient(wfc, grad_wfc);

        // Compute ||G_R||^2 in the S-weighted metric so the descent / CG
        // quantities are consistent with the Stiefel geometry.
        double gnorm2 = energy_grad_->s_inner_product(grad_wfc, grad_wfc);
        result.grad_norm = std::sqrt(std::max(0.0, gnorm2));

        GlobalV::ofs_running << "      orb inner " << inner + 1
            << "  E=" << std::fixed << std::setprecision(10) << E
            << "  gnorm=" << std::scientific << result.grad_norm << std::endl;

        // Orbital sub-problem: exit when Riemannian gradient norm is below threshold.
        if (result.grad_norm < config_.grad_tol)
        {
            result.converged = true;
            result.iterations = inner + 1;
            result.final_energy = E;
            break;
        }

        // ---- Build search direction ----
        //   SD:       d = -G
        //   CG (FR):  d = -G + beta * T(d_prev),  beta = ||G||^2 / ||G_prev||^2
        //             Safeguards: beta >= 0 (FR+), restart on first step and
        //             whenever <G, d> >= 0 (not a descent direction) or the
        //             inner product <G_prev, G> / ||G_prev||^2 is too large
        //             (Powell restart).
        //   L-BFGS:   d computed by EuclideanOptimizer on the flattened
        //             gradient, then projected onto the Stiefel tangent
        //             space.  The update (s, y) pair is computed from the
        //             Riemannian gradient before/after the retraction.
        //   Adam:     same as L-BFGS; Adam already carries its own
        //             learning rate so we step with unit alpha in the
        //             Armijo line search (with backtracking as a safety
        //             net in case the manifold curvature invalidates the
        //             Euclidean step scale).
        psi::Psi<TK> dir(wfc);
        // dir = -grad_wfc
        for (int ik = 0; ik < nk; ++ik)
            for (int ib = 0; ib < nb_local; ++ib)
                for (int mu = 0; mu < nbs_local; ++mu)
                    dir(ik, ib, mu) = -grad_wfc(ik, ib, mu);

        bool restart = true;

        if (use_lbfgs || use_adam)
        {
            // Flatten the Riemannian gradient into a real-valued vector and
            // let the Euclidean optimiser produce a direction.
            std::vector<double> grad_flat, dir_flat;
            psi_to_flat(grad_wfc, grad_flat);
            eucl_opt.compute_direction(grad_flat, dir_flat);

            // Unflatten the direction back into a psi::Psi and project it
            // onto the current tangent space (vector transport by
            // projection).  This is standard for Riemannian L-BFGS and is
            // numerically stable here because the flattened inner product
            // matches the Euclidean ambient metric.
            flat_to_psi(dir_flat, dir);
            energy_grad_->project_orbital_gradient(wfc, dir);

            // Check descent.  If not, fall back to steepest descent (this
            // can happen in the first few iterations of L-BFGS before the
            // Hessian approximation has been built up, or when an Adam
            // momentum term points uphill along the projection).
            double dd_check = energy_grad_->s_inner_product(grad_wfc, dir);
            if (dd_check >= 0.0)
            {
                for (int ik = 0; ik < nk; ++ik)
                    for (int ib = 0; ib < nb_local; ++ib)
                        for (int mu = 0; mu < nbs_local; ++mu)
                            dir(ik, ib, mu) = -grad_wfc(ik, ib, mu);
            }
            restart = false;
        }
        if (use_cg && inner > 0 && prev_gnorm2 > 1e-30)
        {
            // Vector-transport prev_dir to the current tangent space. We use the
            // simplest transport: project onto the new tangent space (this is
            // consistent with our descent-direction projection).
            energy_grad_->project_orbital_gradient(wfc, prev_dir);

            // Powell restart: if <grad, transported grad_prev> / ||grad||^2 is
            // too large, the gradients have lost orthogonality and CG memory is
            // unreliable.
            energy_grad_->project_orbital_gradient(wfc, prev_grad);
            double prev_grad_dot = energy_grad_->s_inner_product(grad_wfc, prev_grad);
            const double powell_thr = 0.1;
            if (std::abs(prev_grad_dot) <= powell_thr * gnorm2)
            {
                double beta = gnorm2 / prev_gnorm2;
                // Polak-Ribiere+ style: enforce beta >= 0 for FR too.
                beta = std::max(0.0, beta);
                // dir += beta * prev_dir
                for (int ik = 0; ik < nk; ++ik)
                    for (int ib = 0; ib < nb_local; ++ib)
                        for (int mu = 0; mu < nbs_local; ++mu)
                            dir(ik, ib, mu) += TK(beta) * prev_dir(ik, ib, mu);

                // Re-project to kill the small tangent-space drift introduced by
                // the linear combination.
                energy_grad_->project_orbital_gradient(wfc, dir);

                // Check descent condition  <grad, dir> < 0.
                double dd_check = energy_grad_->s_inner_product(grad_wfc, dir);
                if (dd_check < 0.0)
                {
                    restart = false;
                }
            }
        }
        if (restart)
        {
            // Steepest descent fallback (also used when use_cg is false)
            for (int ik = 0; ik < nk; ++ik)
                for (int ib = 0; ib < nb_local; ++ib)
                    for (int mu = 0; mu < nbs_local; ++mu)
                        dir(ik, ib, mu) = -grad_wfc(ik, ib, mu);
        }

        // Directional derivative  dd = <grad, dir>  (must be negative)
        double dd = energy_grad_->s_inner_product(grad_wfc, dir);

        // ---- Armijo line search along the retracted direction ----
        // For step size alpha, we want  wfc_new = R_wfc(alpha * dir).
        // retract_orbitals(wfc, step, alpha) internally does
        //     wfc <- wfc - alpha * step,
        // so we pass  step = -dir  to get the effective update  wfc + alpha*dir.
        psi::Psi<TK> neg_dir(dir);
        {
            TK* q = &neg_dir(0, 0, 0);
            const TK* p = &dir(0, 0, 0);
            for (int i = 0; i < total_size; ++i) q[i] = -p[i];
        }

        const double c1 = config_.line_search_c1;
        const double rho = config_.line_search_rho;
        // Initial line-search step size:
        //   - SD/CG: use the configured Armijo start (small, because the
        //            gradient has arbitrary scale).
        //   - L-BFGS: start at 1.0 (the quasi-Newton step is already
        //             properly scaled) and backtrack if needed.
        //   - Adam:   start at 1.0 (Adam incorporates its own learning
        //             rate into the direction).
        double alpha = (use_lbfgs || use_adam) ? 1.0 : config_.line_search_alpha_init;
        double E_new = E;
        bool ls_success = false;

        // Save current wfc for line-search rollback
        for (int ik = 0; ik < nk; ++ik)
            for (int ib = 0; ib < nb_local; ++ib)
                for (int mu = 0; mu < nbs_local; ++mu)
                    wfc_save(ik, ib, mu) = wfc(ik, ib, mu);

        for (int ls = 0; ls < config_.line_search_max_iter; ++ls)
        {
            // Reset wfc from the saved copy
            for (int ik = 0; ik < nk; ++ik)
                for (int ib = 0; ib < nb_local; ++ib)
                    for (int mu = 0; mu < nbs_local; ++mu)
                        wfc(ik, ib, mu) = wfc_save(ik, ib, mu);

            // Retract along +alpha * dir  (i.e. -alpha * neg_dir)
            energy_grad_->retract_orbitals(wfc, neg_dir, alpha);
            energy_grad_->invalidate_hone_cache();

            E_new = energy_grad_->compute_energy(occ_flat, wfc);
            if (E_new <= E + c1 * alpha * dd)
            {
                ls_success = true;
                break;
            }
            alpha *= rho;
        }

        if (!ls_success)
        {
            // Line search failed; revert to previous orbitals.
            for (int ik = 0; ik < nk; ++ik)
                for (int ib = 0; ib < nb_local; ++ib)
                    for (int mu = 0; mu < nbs_local; ++mu)
                        wfc(ik, ib, mu) = wfc_save(ik, ib, mu);
            energy_grad_->invalidate_hone_cache();

            if (use_cg && !restart)
            {
                // Retry with a fresh steepest-descent step before giving up:
                // the CG direction may simply be poorly-conditioned. We force
                // a restart by zeroing prev_dir and continuing.
                prev_gnorm2 = 0.0;
                for (int ik = 0; ik < nk; ++ik)
                    for (int ib = 0; ib < nb_local; ++ib)
                        for (int mu = 0; mu < nbs_local; ++mu)
                            prev_dir(ik, ib, mu) = TK(0);
                result.iterations = inner + 1;
                result.final_energy = E;
                continue;
            }

            if (use_lbfgs || use_adam)
            {
                // Reset the Euclidean optimiser state so the next step
                // starts from steepest descent.  This is the Riemannian
                // analogue of the CG Powell-restart safeguard above.
                eucl_opt.init(flat_size);
                result.iterations = inner + 1;
                result.final_energy = E;
                continue;
            }

            result.iterations = inner + 1;
            result.final_energy = E;
            break;
        }

        // Save this step's gradient and direction for the next CG iteration.
        if (use_cg)
        {
            for (int ik = 0; ik < nk; ++ik)
                for (int ib = 0; ib < nb_local; ++ib)
                    for (int mu = 0; mu < nbs_local; ++mu)
                    {
                        prev_grad(ik, ib, mu) = grad_wfc(ik, ib, mu);
                        prev_dir(ik, ib, mu)  = dir(ik, ib, mu);
                    }
            prev_gnorm2 = gnorm2;
        }

        // Update the L-BFGS history.  For L-BFGS we need the Riemannian
        // gradient at the retracted point; that costs one extra
        // energy+gradient evaluation per step but is essential for the
        // (s, y) curvature pair.  Adam does not need this call because its
        // internal moments are fully updated inside compute_direction() and
        // the pre-existing EuclideanOptimizer::update() increments step_
        // a second time, which would corrupt Adam's bias correction.
        if (use_lbfgs)
        {
            std::vector<double> grad_occ_new;
            psi::Psi<TK> grad_wfc_new;
            energy_grad_->compute(const_cast<std::vector<double>&>(occ_flat),
                                   wfc, grad_occ_new, grad_wfc_new);
            energy_grad_->project_orbital_gradient(wfc, grad_wfc_new);

            // step_vec = alpha * dir in the flattened real-vector space.
            // (For the retraction the effective tangent step is also
            // alpha * dir; the curvature introduced by reorthonormalisation
            // is O(alpha^2) and can safely be absorbed into the L-BFGS
            // approximation error.)
            std::vector<double> new_grad_flat, step_flat, dir_flat_step;
            psi_to_flat(grad_wfc_new, new_grad_flat);
            psi_to_flat(dir, dir_flat_step);
            step_flat.resize(dir_flat_step.size());
            for (size_t i = 0; i < step_flat.size(); ++i)
                step_flat[i] = alpha * dir_flat_step[i];

            eucl_opt.update(new_grad_flat, step_flat);
        }

        energy_grad_->invalidate_hone_cache();
        result.iterations = inner + 1;
        result.final_energy = E_new;
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
        double Ep_one = energy_grad_->E_one_body();
        double Ep_h = energy_grad_->E_hartree();
        double Ep_x = energy_grad_->E_xc();

        double E_minus = energy_grad_->compute_energy(occ_minus,
                            const_cast<psi::Psi<TK>&>(wfc));
        double Em_one = energy_grad_->E_one_body();
        double Em_h = energy_grad_->E_hartree();
        double Em_x = energy_grad_->E_xc();

        double fd_grad = (E_plus - E_minus) / (2.0 * epsilon);
        double fd_one = (Ep_one - Em_one) / (2.0 * epsilon);
        double fd_h   = (Ep_h - Em_h) / (2.0 * epsilon);
        double fd_x   = (Ep_x - Em_x) / (2.0 * epsilon);

        double analytic = grad_occ[idx];
        double rel_err = (std::abs(analytic) > 1e-10)
                         ? std::abs(fd_grad - analytic) / std::abs(analytic)
                         : std::abs(fd_grad - analytic);

        bool pass = rel_err < tolerance;
        if (!pass) all_pass = false;

        GlobalV::ofs_running << std::fixed << std::setprecision(8)
            << "  occ[" << idx << "]: analytic=" << analytic
            << "  fd=" << fd_grad
            << "  (fd_one=" << fd_one
            << "  fd_h=" << fd_h
            << "  fd_x=" << fd_x << ")"
            << "  rel_err=" << rel_err
            << (pass ? "  PASS" : "  FAIL")
            << std::endl;
    }

    // ---- Orbital gradient check (directional derivative along Riemannian G) ----
    GlobalV::ofs_running << "\n-- Orbital gradient check --" << std::endl;
    {
        // Project Euclidean gradient onto tangent space; this is the direction
        // we take the step along.
        psi::Psi<TK> G_tan(grad_wfc);
        energy_grad_->project_orbital_gradient(const_cast<psi::Psi<TK>&>(wfc), G_tan);

        // <G, G_tan> = directional derivative along direction = G_tan
        double dd = 0.0;
        const int nk = wfc.get_nk();
        const int nb_local = wfc.get_nbands();
        const int nbs_local = wfc.get_nbasis();
        for (int ik = 0; ik < nk; ++ik)
            for (int ib = 0; ib < nb_local; ++ib)
                for (int mu = 0; mu < nbs_local; ++mu)
                {
                    // Real part of conj(grad_wfc) * G_tan
                    TK a = grad_wfc(ik, ib, mu);
                    TK b = G_tan(ik, ib, mu);
                    if constexpr (std::is_same<TK, double>::value)
                        dd += a * b;
                    else
                        dd += (std::conj(a) * b).real();
                }
#ifdef __MPI
        Parallel_Reduce::reduce_all(dd);
#endif

        double gnorm2 = 0.0;
        for (int ik = 0; ik < nk; ++ik)
            for (int ib = 0; ib < nb_local; ++ib)
                for (int mu = 0; mu < nbs_local; ++mu)
                    gnorm2 += abs2(G_tan(ik, ib, mu));
#ifdef __MPI
        Parallel_Reduce::reduce_all(gnorm2);
#endif

        // f(t) = E(occ, R_C(-t * G_tan)).  f'(0) = -<G, G_tan> via Armijo convention.
        auto f_at = [&](double t) -> double {
            psi::Psi<TK> wfc_trial(wfc);
            energy_grad_->retract_orbitals(wfc_trial, G_tan, t);
            energy_grad_->invalidate_hone_cache();
            return energy_grad_->compute_energy(occ_flat, wfc_trial);
        };
        const double t = epsilon;
        double f_plus  = f_at(+t);
        double f_minus = f_at(-t);
        double fd = (f_plus - f_minus) / (2.0 * t);
        // Along the steepest-descent direction the analytic directional derivative
        // is  -<grad, G_tan> = -||G_tan||^2 up to manifold curvature.
        double analytic_dd = -dd;
        double rel_err = (std::abs(analytic_dd) > 1e-10)
                         ? std::abs(fd - analytic_dd) / std::abs(analytic_dd)
                         : std::abs(fd - analytic_dd);
        bool pass = rel_err < tolerance;
        if (!pass) all_pass = false;
        GlobalV::ofs_running << std::fixed << std::setprecision(8)
            << "  orb dir deriv: analytic=" << analytic_dd
            << "  fd=" << fd
            << "  ||G_R||^2=" << gnorm2
            << "  rel_err=" << rel_err
            << (pass ? "  PASS" : "  FAIL") << std::endl;

        // Restore state: invalidate cache so next compute() recomputes correctly.
        energy_grad_->invalidate_hone_cache();
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
