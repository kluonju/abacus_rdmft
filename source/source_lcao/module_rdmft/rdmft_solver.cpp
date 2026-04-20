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

void print_rdmft_outer_energy_stdout(bool converged, double E)
{
    std::cout << std::fixed << std::setprecision(10);
    if (converged)
    {
        std::cout << "  RDMFT outer: converged  E = " << E << std::defaultfloat << std::endl;
    }
    else
    {
        std::cout << "  RDMFT outer: NOT converged  E = " << E << std::defaultfloat << std::endl;
    }
}

void print_occ_table_stdout(const std::vector<double>& occ_flat, int nk, int nbands)
{
    std::cout << "  WARNING: RDMFT outer loop did not converge. Occupations n(ik, ib):" << std::endl;
    std::cout << std::left << std::setw(6) << "ik" << std::setw(8) << "ib" << std::setw(18) << "n"
              << std::endl;
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            const int idx = ik * nbands + ib;
            std::cout << std::left << std::setw(6) << ik << std::setw(8) << ib << std::fixed
                      << std::setprecision(10) << std::setw(18) << occ_flat[idx] << std::defaultfloat
                      << std::endl;
        }
    }
}

void print_nonconverged_report_running_alternating(double E,
                                                   double dE,
                                                   double abs_c,
                                                   double occ_gnorm,
                                                   double orb_gnorm,
                                                   const std::vector<double>& occ_final,
                                                   const std::vector<double>& occ_start_last_outer,
                                                   int nk,
                                                   int nbands)
{
    GlobalV::ofs_running << "\n  WARNING: RDMFT outer loop did not converge." << std::endl;
    GlobalV::ofs_running << std::fixed << std::setprecision(10) << "  Last outer: E = " << E
                         << "  dE = " << std::scientific << dE << std::fixed << "  |c| = " << abs_c
                         << "  occ_gnorm = " << occ_gnorm << "  orb_gnorm = " << orb_gnorm
                         << std::endl;
    const bool have_prev = (occ_start_last_outer.size() == occ_final.size());
    double sum_abs_dn = 0.0;
    GlobalV::ofs_running << std::left << std::setw(6) << "ik" << std::setw(8) << "ib" << std::setw(18)
                         << "n" << std::setw(18) << "dn" << std::endl;
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            const int idx = ik * nbands + ib;
            const double n = occ_final[idx];
            const double dn = have_prev ? (n - occ_start_last_outer[idx]) : 0.0;
            if (have_prev)
            {
                sum_abs_dn += std::abs(dn);
            }
            GlobalV::ofs_running << std::left << std::setw(6) << ik << std::setw(8) << ib << std::fixed
                                 << std::setprecision(10) << std::setw(18) << n << std::setw(18) << dn
                                 << std::endl;
        }
    }
    if (have_prev)
    {
        GlobalV::ofs_running << "  sum|dn| (last outer cycle) = " << std::scientific << sum_abs_dn
                             << std::endl;
    }
}

void print_nonconverged_report_running_joint(double E,
                                             double dE,
                                             double abs_c,
                                             double gnorm_occ,
                                             double gnorm_orb,
                                             double gnorm_total,
                                             const std::vector<double>& occ_final,
                                             const std::vector<double>& occ_start_last_outer,
                                             int nk,
                                             int nbands)
{
    GlobalV::ofs_running << "\n  WARNING: RDMFT outer loop did not converge." << std::endl;
    GlobalV::ofs_running << std::fixed << std::setprecision(10) << "  Last outer: E = " << E
                         << "  dE = " << std::scientific << dE << std::fixed << "  |c| = " << abs_c
                         << "  |grad_occ| = " << gnorm_occ << "  |grad_orb| = " << gnorm_orb
                         << "  |grad_total| = " << gnorm_total << std::endl;
    const bool have_prev = (occ_start_last_outer.size() == occ_final.size());
    double sum_abs_dn = 0.0;
    GlobalV::ofs_running << std::left << std::setw(6) << "ik" << std::setw(8) << "ib" << std::setw(18)
                         << "n" << std::setw(18) << "dn" << std::endl;
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            const int idx = ik * nbands + ib;
            const double n = occ_final[idx];
            const double dn = have_prev ? (n - occ_start_last_outer[idx]) : 0.0;
            if (have_prev)
            {
                sum_abs_dn += std::abs(dn);
            }
            GlobalV::ofs_running << std::left << std::setw(6) << ik << std::setw(8) << ib << std::fixed
                                 << std::setprecision(10) << std::setw(18) << n << std::setw(18) << dn
                                 << std::endl;
        }
    }
    if (have_prev)
    {
        GlobalV::ofs_running << "  sum|dn| (last outer cycle) = " << std::scientific << sum_abs_dn
                             << std::endl;
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
    last_result_ = {};

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

    // S-orthonormalize the initial orbitals so that all subsequent gradient
    // projections and Cholesky retractions start from a valid point on the
    // generalised Stiefel manifold { C : C^H S C = I }.  Calling
    // retract_orbitals with alpha = 0 applies only the Cholesky-QR
    // re-orthonormalisation without any directional step (the gradient
    // argument is irrelevant when alpha = 0, so we reuse wfc itself).
    energy_grad_->retract_orbitals(wfc, wfc, 0.0);

    // The X-variable path (X = U*C with S = U^H U) is currently disabled in
    // MPI runs because it can trigger heap corruption on layouts where some
    // ranks own zero local blocks. Keep the stable C-space manifold updates.
#ifdef __MPI
    const bool enable_cholesky_x_variable = false;
#else
    const bool enable_cholesky_x_variable = true;
#endif

    energy_grad_->disable_X_variable();
    if (enable_cholesky_x_variable)
    {
        // Precompute the Cholesky factorisation S_k = U_k^H U_k and switch to
        // the X_k = U_k C_k variable for all subsequent manifold operations.
        // In X-space X^H X = I, so projection, retraction, and inner product
        // are the standard (S = I) Stiefel forms.
        energy_grad_->precompute_cholesky_S();
        if (energy_grad_->use_X_variable())
        {
            energy_grad_->wfc_C_to_X(wfc);
        }
    }

    double E = 0.0;
    switch (config_.strategy)
    {
        case SolverStrategy::Alternating:
            E = solve_alternating(occ_flat, wfc);
            break;
        case SolverStrategy::Joint:
            E = solve_joint(occ_flat, wfc);
            break;
    }

    // Transform back X -> C before returning to the caller.
    if (enable_cholesky_x_variable && energy_grad_->use_X_variable())
    {
        energy_grad_->wfc_X_to_C(wfc);
        energy_grad_->disable_X_variable();
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

    std::vector<double> occ_at_outer_start;
    double last_dE = 0.0;
    double last_abs_c = 0.0;
    double last_occ_gnorm = 0.0;
    double last_orb_gnorm = 0.0;
    int outer_iters_done = 0;

    for (int iter = 0; iter < config_.outer_maxiter; ++iter)
    {
        occ_at_outer_start.assign(occ_flat.begin(), occ_flat.end());

        if (config_.print_stiefel_gram)
        {
            std::vector<double> stiefel_gram_frob_per_ik;
            energy_grad_->stiefel_gram_residual_frobenius_per_k(wfc, stiefel_gram_frob_per_ik);
            GlobalV::ofs_running << "    Stiefel Gram residual ||G_k - I||_F per k-point (G_k = "
                << (energy_grad_->use_X_variable() ? "X_k^H X_k" : "C_k^H S_k C_k") << "):" << std::endl;
            for (int ik = 0; ik < wfc.get_nk(); ++ik)
            {
                GlobalV::ofs_running << std::fixed << std::setprecision(10) << "      ik=" << ik
                    << "  ||G_k - I||_F = " << stiefel_gram_frob_per_ik[ik] << std::endl;
            }
        }

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

        last_dE = dE;
        last_abs_c = std::abs(constraint_viol);
        last_occ_gnorm = occ_result.grad_norm;
        last_orb_gnorm = orb_result.grad_norm;
        outer_iters_done = iter + 1;

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
    if (!last_result_.converged)
    {
        last_result_.iterations = outer_iters_done;
        last_result_.grad_norm = std::max(last_occ_gnorm, last_orb_gnorm);
    }

    print_rdmft_outer_energy_stdout(last_result_.converged, E);
    if (!last_result_.converged)
    {
        print_occ_table_stdout(occ_flat, nk_, nbands_);
        print_nonconverged_report_running_alternating(E,
                                                      last_dE,
                                                      last_abs_c,
                                                      last_occ_gnorm,
                                                      last_orb_gnorm,
                                                      occ_flat,
                                                      occ_at_outer_start,
                                                      nk_,
                                                      nbands_);
    }
    return E;
}

// ----------------------------------------------------------------------------
// Joint (product-manifold) optimisation
//
// We treat the point
//     x = (p, C^1, ..., C^{Nk})
// where p are the unconstrained occupation parameters (cosine^2 / logistic
// already handles the box constraint n in [0,1]) and each C^k lives on the
// generalised Stiefel manifold St(Nb, Nbasis; S^k), as a single point on the
// product manifold
//     M = R^{Nk x Nb}  x  prod_k St(Nb, Nbasis; S^k).
//
// Key design: occupations and orbitals are packed into ONE flat vector
//     z = (p_1, ..., p_{Np},   flat(C^1), ..., flat(C^{Nk}))
// and a SINGLE Euclidean optimiser (configured via rdmft_joint_optimizer)
// consumes the packed gradient
//     g = (dE/dp_1, ..., dE/dp_{Np},   flat(G_R^1), ..., flat(G_R^{Nk}))
// where G_R^k is the Riemannian (tangent-space) projection of the Euclidean
// orbital gradient dE/dC^k at the current C^k. The produced search direction
// is split back into an occupation block (linear update on parameters p) and
// an orbital block (projected and retracted onto each Stiefel fibre) before a
// single Armijo line search along the packed direction commits the step. The
// augmented-Lagrangian multiplier is refreshed once per outer iteration.
// ----------------------------------------------------------------------------
template <typename TK, typename TR>
double RDMFTSolver<TK, TR>::solve_joint(
    std::vector<double>& occ_flat,
    psi::Psi<TK>& wfc)
{
    GlobalV::ofs_running << "\n===== RDMFT Joint (Product-Manifold) Optimization ====="
                         << std::endl;

    // Convert occupations to unconstrained parameters.
    std::vector<double> params(occ_flat.size());
    occ_param_->occ_to_params(occ_flat, params);

    const int n_occ_params = static_cast<int>(params.size());
    const int nk = wfc.get_nk();
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();
    const int orb_total_size = nk * nb_local * nbs_local;

    // Flat-size of the orbital block in the packed vector: one real per
    // orbital coefficient when TK = double, two reals (Re, Im) when TK is
    // complex.
    int orb_flat_size = orb_total_size;
    if constexpr (!std::is_same<TK, double>::value) orb_flat_size *= 2;

    const int packed_size = n_occ_params + orb_flat_size;

    // Single unified optimiser on the packed variable (p, C_flat).
    const OptimizerType joint_type = config_.joint_optimizer;
    const bool joint_is_lbfgs = (joint_type == OptimizerType::LBFGS);
    const bool joint_is_adam  = (joint_type == OptimizerType::Adam);
    EuclideanOptimizer joint_opt(joint_type, config_);
    joint_opt.init(packed_size);

    // Work buffer for line-search rollback.
    psi::Psi<TK> wfc_save(wfc);

    double E_prev = 1e30;
    double E = 0.0;

    std::vector<double> occ_at_outer_start;
    double joint_last_E = 0.0;
    double joint_last_dE = 0.0;
    double joint_last_abs_c = 0.0;
    double joint_gn_occ = 0.0;
    double joint_gn_orb = 0.0;
    double joint_gn_tot = 0.0;
    int joint_outer_done = 0;

    auto total_energy = [&](const std::vector<double>& occ_in,
                            const psi::Psi<TK>& wfc_in) -> double {
        double E_val = energy_grad_->compute_energy(occ_in,
                                const_cast<psi::Psi<TK>&>(wfc_in));
        if (config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
            E_val += occ_constraint_->augmented_lagrangian_penalty(occ_in);
        return E_val;
    };

    // Scratch buffers reused across outer iterations.
    std::vector<double> packed_grad(packed_size, 0.0);
    std::vector<double> packed_dir(packed_size, 0.0);
    std::vector<double> packed_step(packed_size, 0.0);

    for (int iter = 0; iter < config_.outer_maxiter; ++iter)
    {
        // Refresh occupations from params so downstream code always sees a
        // consistent (p, n) pair.
        occ_param_->params_to_occ(params, occ_flat);
        occ_at_outer_start.assign(occ_flat.begin(), occ_flat.end());

        // Full energy + Euclidean gradients at the current point.
        std::vector<double> grad_occ;
        psi::Psi<TK> grad_wfc;
        E = energy_grad_->compute(occ_flat, wfc, grad_occ, grad_wfc);

        // Add augmented Lagrangian penalty.
        if (config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
        {
            E += occ_constraint_->augmented_lagrangian_penalty(occ_flat);
            std::vector<double> penalty_grad;
            occ_constraint_->augmented_lagrangian_gradient(occ_flat, penalty_grad);
            for (size_t i = 0; i < grad_occ.size(); ++i)
                grad_occ[i] += penalty_grad[i];
        }

        // Transform occupation gradient dE/dn -> dE/dp (chain rule through the
        // cosine^2 / logistic parameterisation).
        std::vector<double> grad_params;
        occ_param_->transform_gradient_batch(grad_occ, params, grad_params);

        // Project orbital gradient onto the Stiefel tangent space at C.
        energy_grad_->project_orbital_gradient(wfc, grad_wfc);

        // ---- Pack gradients into the unified vector ----
        //
        // We multiply the orbital block by joint_orb_scale both on the way
        // in (gradient) and on the way out (direction) so that the packed
        // optimiser effectively works with a rescaled orbital coefficient
        //   C_internal = C / joint_orb_scale
        // while alpha remains a single scalar. This is exactly a diagonal
        // preconditioner on the packed vector; it does not change the
        // search direction's sign, only its relative magnitude per block.
        // For a SD direction at iteration 1 this gives
        //   alpha * orb_step = alpha * joint_orb_scale^2 * grad_wfc
        // so a value of joint_orb_scale < 1 reduces the physical step in
        // the orbital block while leaving the occupation block unaffected.
        const double orb_scale = config_.joint_orb_scale;
        std::vector<double> grad_orb_flat;
        psi_to_flat(grad_wfc, grad_orb_flat);
        for (int i = 0; i < n_occ_params; ++i)
            packed_grad[i] = grad_params[i];
        for (int i = 0; i < orb_flat_size; ++i)
            packed_grad[n_occ_params + i] = orb_scale * grad_orb_flat[i];

        // Ask the single unified optimiser for a packed descent direction.
        joint_opt.compute_direction(packed_grad, packed_dir);

        // Split the packed direction into occupation and orbital blocks.
        // Apply the orbital-block scale on the direction as well so that
        // alpha * orb_dir (physical) = alpha * orb_scale * orb_dir_internal.
        std::vector<double> occ_dir(n_occ_params);
        for (int i = 0; i < n_occ_params; ++i) occ_dir[i] = packed_dir[i];

        psi::Psi<TK> orb_dir(wfc);
        {
            std::vector<double> orb_dir_flat(orb_flat_size);
            for (int i = 0; i < orb_flat_size; ++i)
                orb_dir_flat[i] = orb_scale * packed_dir[n_occ_params + i];
            flat_to_psi(orb_dir_flat, orb_dir);
        }

        // Orbital gradient norm (Riemannian) for diagnostics and for the
        // steepest-descent fallback below.
        const double orb_gnorm2 = energy_grad_->s_inner_product(grad_wfc, grad_wfc);

        // Project the orbital search direction back onto the tangent space
        // at C. The optimiser's Euclidean update (especially for L-BFGS /
        // Adam, which precondition the gradient) generally leaves the tangent
        // space; projection restores a valid Riemannian direction without
        // changing its component on the tangent space.
        energy_grad_->project_orbital_gradient(wfc, orb_dir);

        // Directional derivative of the augmented energy along the packed
        // direction:  dd_total = <grad_p, occ_dir> + <G_R, orb_dir>_S.
        double occ_dd = 0.0;
        for (int i = 0; i < n_occ_params; ++i)
            occ_dd += occ_dir[i] * grad_params[i];
        double orb_dd = energy_grad_->s_inner_product(grad_wfc, orb_dir);
        double dd_total = occ_dd + orb_dd;

        // Descent safeguard: if the optimiser-produced direction is not
        // descent on the product manifold, fall back to the full packed
        // steepest-descent direction d = -g.
        if (dd_total >= 0.0)
        {
            for (int i = 0; i < n_occ_params; ++i)
                occ_dir[i] = -grad_params[i];
            for (int ik = 0; ik < nk; ++ik)
                for (int ib = 0; ib < nb_local; ++ib)
                    for (int mu = 0; mu < nbs_local; ++mu)
                        orb_dir(ik, ib, mu) = -grad_wfc(ik, ib, mu);
            occ_dd = 0.0;
            for (int i = 0; i < n_occ_params; ++i)
                occ_dd += occ_dir[i] * grad_params[i];
            orb_dd = -orb_gnorm2;
            dd_total = occ_dd + orb_dd;
        }

        // neg_orb_dir = -orb_dir (so retract_orbitals(wfc, neg_orb_dir, alpha)
        // applies wfc <- wfc - alpha * (-orb_dir) = wfc + alpha * orb_dir).
        psi::Psi<TK> neg_orb_dir(orb_dir);
        {
            TK* q = &neg_orb_dir(0, 0, 0);
            const TK* p = &orb_dir(0, 0, 0);
            for (int i = 0; i < orb_total_size; ++i) q[i] = -p[i];
        }

        // Decide whether the orbital block actually needs to be stepped in
        // this iteration. Mirroring the alternating strategy's
        // optimize_orbitals (which early-exits on ||G_R|| < grad_tol and
        // therefore never invokes retract_orbitals with a zero / tiny
        // direction), we skip the orbital retraction when the Riemannian
        // orbital gradient is below grad_tol. Without this guard the
        // Cholesky-QR S-orthonormalisation inside retract_orbitals is
        // executed every iteration, which is a numerically non-trivial
        // O(nbasis^3) update whose round-off amplifies when the trial
        // step is otherwise supposed to leave the orbitals unchanged; on
        // H2 (where the KS seed is already a fixed point of the orbital
        // sub-problem, so ||G_R|| == 0) this spurious re-orthonormalisation
        // perturbs C by a few percent at the start and destabilises the
        // occupation line search.
        const double orb_gnorm = std::sqrt(std::max(0.0, orb_gnorm2));
        const bool step_orbitals = (orb_gnorm >= config_.grad_tol);

        // Snapshot the current orbitals for line-search rollback (needed
        // only when we will actually retract).
        if (step_orbitals)
        {
            for (int ik = 0; ik < nk; ++ik)
                for (int ib = 0; ib < nb_local; ++ib)
                    for (int mu = 0; mu < nbs_local; ++mu)
                        wfc_save(ik, ib, mu) = wfc(ik, ib, mu);
        }

        // Joint Armijo line search along the packed direction. Initial step
        // is 1.0 for L-BFGS / Adam (already well-scaled) and the configured
        // alpha_init for SD / CG.
        const double alpha_init = (joint_is_lbfgs || joint_is_adam)
                                      ? 1.0
                                      : config_.line_search_alpha_init;
        double alpha = alpha_init;
        const double c1 = config_.line_search_c1;
        const double rho = config_.line_search_rho;
        double E_new = E;
        bool ls_success = false;

        std::vector<double> occ_flat_new;
        std::vector<double> params_new(params.size());

        for (int ls = 0; ls < config_.line_search_max_iter; ++ls)
        {
            if (step_orbitals)
            {
                // Roll back orbitals and retract along +alpha * orb_dir.
                for (int ik = 0; ik < nk; ++ik)
                    for (int ib = 0; ib < nb_local; ++ib)
                        for (int mu = 0; mu < nbs_local; ++mu)
                            wfc(ik, ib, mu) = wfc_save(ik, ib, mu);
                energy_grad_->retract_orbitals(wfc, neg_orb_dir, alpha);
                energy_grad_->invalidate_hone_cache();
            }

            for (int i = 0; i < n_occ_params; ++i)
                params_new[i] = params[i] + alpha * occ_dir[i];
            occ_param_->params_to_occ(params_new, occ_flat_new);

            E_new = total_energy(occ_flat_new, wfc);
            if (E_new <= E + c1 * alpha * dd_total)
            {
                ls_success = true;
                break;
            }
            alpha *= rho;
        }

        if (!ls_success)
        {
            // Line search failed: roll back, restart the optimiser state and
            // continue with steepest descent at the next iteration.
            if (step_orbitals)
            {
                for (int ik = 0; ik < nk; ++ik)
                    for (int ib = 0; ib < nb_local; ++ib)
                        for (int mu = 0; mu < nbs_local; ++mu)
                            wfc(ik, ib, mu) = wfc_save(ik, ib, mu);
                energy_grad_->invalidate_hone_cache();
            }

            joint_opt.init(packed_size);

            std::cout << "  RDMFT joint-iter " << iter + 1
                      << "  line search failed, restarting optimiser state"
                      << std::endl;
            GlobalV::ofs_running << "  RDMFT joint-iter " << iter + 1
                << "  line search failed, restarting optimiser state" << std::endl;
            E_prev = E;
            continue;
        }

        // Commit occupation step (orbitals are already retracted by the
        // final, successful line-search trial above).
        for (int i = 0; i < n_occ_params; ++i)
        {
            packed_step[i] = alpha * occ_dir[i];
            params[i] += packed_step[i];
        }
        occ_param_->params_to_occ(params, occ_flat);

        // Build the packed step for the orbital block. Store the step in
        // the *internal* (scaled) units that match packed_dir:
        //   step_internal = alpha * dir_internal = alpha * orb_dir / orb_scale
        // so that the optimiser's (s, y) history stays consistent with the
        // rescaled gradient feed.
        {
            std::vector<double> orb_dir_flat;
            psi_to_flat(orb_dir, orb_dir_flat);
            const double inv_scale = (orb_scale != 0.0) ? (1.0 / orb_scale) : 1.0;
            for (int i = 0; i < orb_flat_size; ++i)
                packed_step[n_occ_params + i] = alpha * orb_dir_flat[i] * inv_scale;
        }

        // Update the unified optimiser's history with the new packed gradient.
        if (joint_is_lbfgs || joint_is_adam)
        {
            std::vector<double> new_grad_occ;
            psi::Psi<TK> new_grad_wfc;
            energy_grad_->compute(occ_flat, wfc, new_grad_occ, new_grad_wfc);
            if (config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
            {
                std::vector<double> pen_grad;
                occ_constraint_->augmented_lagrangian_gradient(occ_flat, pen_grad);
                for (size_t i = 0; i < new_grad_occ.size(); ++i)
                    new_grad_occ[i] += pen_grad[i];
            }
            std::vector<double> new_grad_params;
            occ_param_->transform_gradient_batch(new_grad_occ, params, new_grad_params);
            energy_grad_->project_orbital_gradient(wfc, new_grad_wfc);

            std::vector<double> new_grad_orb_flat;
            psi_to_flat(new_grad_wfc, new_grad_orb_flat);

            std::vector<double> new_packed_grad(packed_size);
            for (int i = 0; i < n_occ_params; ++i)
                new_packed_grad[i] = new_grad_params[i];
            for (int i = 0; i < orb_flat_size; ++i)
                new_packed_grad[n_occ_params + i] = orb_scale * new_grad_orb_flat[i];

            joint_opt.update(new_packed_grad, packed_step);
        }

        // Augmented-Lagrangian multiplier refresh.
        if (config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
        {
            occ_constraint_->update_multiplier(occ_flat);
            if (std::abs(occ_constraint_->constraint_violation(occ_flat)) > 1e-6)
            {
                occ_constraint_->increase_penalty(
                    config_.aug_lag_mu_factor, config_.aug_lag_mu_max);
            }
        }

        const double dE = std::abs(E_new - E_prev);
        double gnorm2_occ = 0.0;
        for (auto g : grad_params) gnorm2_occ += g * g;
        const double gnorm_total = std::sqrt(gnorm2_occ + std::max(0.0, orb_gnorm2));

        std::cout << std::fixed << std::setprecision(10)
            << "  RDMFT joint-iter " << iter + 1
            << "  E = " << E_new
            << "  dE = " << std::scientific << dE
            << "  alpha = " << alpha
            << "  |grad| = " << gnorm_total
            << std::endl;

        GlobalV::ofs_running << std::fixed << std::setprecision(10)
            << "  RDMFT joint-iter " << iter + 1
            << "  E = " << E_new
            << "  dE = " << std::scientific << dE
            << "  alpha = " << alpha
            << "  |grad_occ| = " << std::sqrt(gnorm2_occ)
            << "  |grad_orb| = " << std::sqrt(std::max(0.0, orb_gnorm2))
            << "  |grad_total| = " << gnorm_total
            << std::endl;

        const double abs_c_joint = std::abs(occ_constraint_->constraint_violation(occ_flat));
        joint_last_E = E_new;
        joint_last_dE = dE;
        joint_last_abs_c = abs_c_joint;
        joint_gn_occ = std::sqrt(gnorm2_occ);
        joint_gn_orb = std::sqrt(std::max(0.0, orb_gnorm2));
        joint_gn_tot = gnorm_total;
        joint_outer_done = iter + 1;

        if (dE < config_.energy_tol && gnorm_total < config_.grad_tol)
        {
            last_result_.converged = true;
            last_result_.iterations = iter + 1;
            last_result_.final_energy = E_new;
            last_result_.grad_norm = gnorm_total;
            E = E_new;
            break;
        }
        E_prev = E_new;
        E = E_new;
    }

    occ_param_->params_to_occ(params, occ_flat);
    last_result_.final_energy = E;
    if (!last_result_.converged)
    {
        last_result_.iterations = joint_outer_done;
        last_result_.grad_norm = joint_gn_tot;
    }

    print_rdmft_outer_energy_stdout(last_result_.converged, E);
    if (!last_result_.converged)
    {
        print_occ_table_stdout(occ_flat, nk_, nbands_);
        print_nonconverged_report_running_joint(joint_last_E,
                                                joint_last_dE,
                                                joint_last_abs_c,
                                                joint_gn_occ,
                                                joint_gn_orb,
                                                joint_gn_tot,
                                                occ_flat,
                                                occ_at_outer_start,
                                                nk_,
                                                nbands_);
    }
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

        // Project onto the tangent space of the Stiefel manifold.
        // In C-space (generalised Stiefel with overlap S):
        //     G_R = G - C sym(C^H S G)
        // In X-space (standard Stiefel, use_X_variable_ = true):
        //     G_R = G - X sym(X^H G)
        energy_grad_->project_orbital_gradient(wfc, grad_wfc);

        // Compute ||G_R||^2 in the appropriate metric (S-weighted or
        // Euclidean when in X-space) for descent / CG consistency.
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
