#include "rdmft_solver.h"
#include "source_base/constants.h"
#include "source_base/formatter.h"
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
#include <complex>
#include <chrono>
#include <limits>
#include <stdexcept>

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

/// Line search in RDMFT: Strong Wolfe for LBFGS *and* nonlinear CG; Armijo
/// (backtracking) for SD / Adam.
///
/// Nonlinear CG on a non-convex objective (e.g. RDMFT Müller / Power
/// exchange on the Stiefel manifold) requires Strong Wolfe (Nocedal &
/// Wright §3.5): the curvature condition `|g(alpha)| <= c2 |g(0)|`
/// rejects line-search trials that overshoot the linear regime, which
/// pure Armijo cannot detect. Pure Armijo accepts any trial that beats
/// `f(0) + c1*alpha*g(0)` with no upper bound on the size of the
/// "downhill" overshoot, and on this landscape the Power / Müller
/// exchange has spurious low-energy regions far from the physical
/// basin. The bug shows up immediately on the user's NiO test (clean
/// KS init, nelec_delta=0): the second orbital line-search trial
/// accepts alpha=1.0 with one feval and drops E from -338 Ry to
/// -806 Ry on a `dd ~ -6`, then PR+ amplifies the gradient mismatch
/// by 6 orders of magnitude on the next iter, and the iterate is in
/// the wrong basin permanently.
inline bool line_search_uses_strong_wolfe(const OptimizerType opt)
{
    return opt == OptimizerType::LBFGS
        || opt == OptimizerType::ConjugateGradient;
}

/// Stop augmented-Lagrangian occupation inner loop when the total occupation change is small.
bool occ_inner_should_stop(double sum_abs_dn,
                           double dn_tol)
{
    return (sum_abs_dn < dn_tol);
}

inline double euclidean_norm(const std::vector<double>& v)
{
    double s = 0.0;
    for (double x : v)
        s += x * x;
    return std::sqrt(s);
}

/// Stationarity measure r = n - P(n - g) for the Spectral Projected Gradient
/// (Bertsekas projected-gradient map at unit τ, the textbook choice; the
/// resulting `r` is τ-free and well-defined regardless of the spectral step).
/// Returns ‖r‖_2 and ‖r‖_∞.
void spg_stationarity_norms(const std::vector<double>& occ,
                            const std::vector<double>& grad,
                            const OccupationConstraint& constraint,
                            double& l2,
                            double& linf)
{
    std::vector<double> trial(occ.size());
    for (size_t i = 0; i < occ.size(); ++i)
    {
        trial[i] = occ[i] - grad[i];
    }
    constraint.project(trial);

    double n2 = 0.0;
    linf = 0.0;
    for (size_t i = 0; i < occ.size(); ++i)
    {
        const double d = occ[i] - trial[i];
        n2 += d * d;
        linf = std::max(linf, std::abs(d));
    }
    l2 = std::sqrt(std::max(0.0, n2));
}

double active_set_dual_complementarity_violation(const std::vector<double>& grad,
                                                 const OccupationConstraint::ActiveSetInfo& as_info,
                                                 const OccupationConstraint& constraint)
{
    const int nbands = constraint.nbands();
    const auto& kweights = constraint.kweights();
    double max_violation = 0.0;
    for (size_t idx = 0; idx < grad.size(); ++idx)
    {
        const int ik = static_cast<int>(idx) / nbands;
        const double wk = kweights[ik];
        const double gkkt = grad[idx] + as_info.lagrange_mult * wk;

        // KKT sign conditions on active bounds.
        if (as_info.at_lower[idx])
        {
            max_violation = std::max(max_violation, std::max(0.0, -gkkt));
        }
        else if (as_info.at_upper[idx])
        {
            max_violation = std::max(max_violation, std::max(0.0, gkkt));
        }
    }
    return max_violation;
}

/// Spectral (Barzilai-Borwein) step length for the SPG occupation inner.
///
/// Standard formula (BB1) `α = (s,s)/(s,y)` with `s = n − n_prev`,
/// `y = g − g_prev`.  Returns `fallback` when no previous iterate is available
/// or the curvature pair is degenerate.  Clamped to `[alpha_min, alpha_max]`
/// (Birgin–Martínez–Raydan safeguarding).
inline double spectral_bb_step(const std::vector<double>& occ,
                               const std::vector<double>& grad,
                               const std::vector<double>& occ_prev,
                               const std::vector<double>& grad_prev,
                               bool have_prev,
                               double fallback,
                               double alpha_min,
                               double alpha_max)
{
    if (!have_prev || occ_prev.size() != occ.size() || grad_prev.size() != grad.size())
    {
        return std::min(alpha_max, std::max(alpha_min, fallback));
    }
    double ss = 0.0;
    double sy = 0.0;
    for (size_t i = 0; i < occ.size(); ++i)
    {
        const double s = occ[i] - occ_prev[i];
        const double y = grad[i] - grad_prev[i];
        ss += s * s;
        sy += s * y;
    }
    double alpha = fallback;
    if (sy > 1e-30 && std::isfinite(ss) && std::isfinite(sy))
    {
        alpha = ss / sy;
    }
    if (!std::isfinite(alpha) || alpha <= 0.0)
    {
        alpha = fallback;
    }
    return std::min(alpha_max, std::max(alpha_min, alpha));
}

int find_fermi_boundary_index(const std::vector<double>& occ_flat,
                              const int ik,
                              const int nbands)
{
    // Occupations are in [0,1] in RDMFT spin-channel convention.
    // The Fermi boundary is approximated by the first band with n < 0.5.
    for (int ib = 0; ib < nbands; ++ib)
    {
        if (occ_flat[ik * nbands + ib] < 0.5)
        {
            return ib - 1;
        }
    }

    // Fallback for atypical ordering: last non-negligibly occupied band.
    int last_occ = -1;
    for (int ib = 0; ib < nbands; ++ib)
    {
        if (occ_flat[ik * nbands + ib] > 1e-8)
        {
            last_occ = ib;
        }
    }
    return last_occ;
}

/// Last band index ib with cumulative \(\sum_{jb\le ib} w_k n_{\mathrm{KS}}(jb) < n_{\mathrm{target},k}\).
/// Used when all KS occupations at this k are \(\ge 0.5\) so the 0.5 rule gives no gap (e.g. nelec_delta).
int find_fermi_boundary_from_electron_target(const std::vector<double>& occ_ks,
                                             int ik,
                                             int nbands,
                                             double wk,
                                             double n_target_k)
{
    if (n_target_k <= 0.0)
    {
        return -1;
    }
    double acc = 0.0;
    int ib_last = -1;
    for (int ib = 0; ib < nbands; ++ib)
    {
        const double n = std::max(0.0, std::min(1.0, occ_ks[ik * nbands + ib]));
        const double contrib = wk * n;
        if (acc + contrib < n_target_k - 1e-10)
        {
            ib_last = ib;
            acc += contrib;
        }
        else
        {
            break;
        }
    }
    return ib_last;
}

/// Prefer the 0.5-gap Fermi index for metallic KS; otherwise use the electron-target boundary.
int find_init_fermi_boundary(const std::vector<double>& occ_ks,
                           int ik,
                           int nbands,
                           double wk,
                           double n_target_k)
{
    bool has_below_half = false;
    for (int ib = 0; ib < nbands; ++ib)
    {
        if (occ_ks[ik * nbands + ib] < 0.5)
        {
            has_below_half = true;
            break;
        }
    }
    if (has_below_half)
    {
        return find_fermi_boundary_index(occ_ks, ik, nbands);
    }
    return find_fermi_boundary_from_electron_target(occ_ks, ik, nbands, wk, n_target_k);
}

void log_occ_inner_line(const std::string& line)
{
    GlobalV::ofs_running << line << std::endl;
}

/// One summary line, then per-ik lines listing each n(ik,ib) and dn vs the *previous*
/// inner iteration's start (occ_prev_for_dn empty on first inner => table dn = 0 by convention).
void log_occ_inner_summary_and_nik(const std::string& summary_first_line,
                                   const std::vector<double>& occ_flat,
                                   const std::vector<double>& occ_prev_for_dn,
                                   const int nk,
                                   const int nbands)
{
    log_occ_inner_line(summary_first_line);

    const bool have_prev = (occ_prev_for_dn.size() == occ_flat.size());
    if (!have_prev)
    {
        log_occ_inner_line("      note: table dn compares n to the start of the *previous* inner iteration; "
                           "the first inner lists dn=0 by convention (not the step change). "
                           "The line sum|dn|_step is the L1 size of the occupation change that step, "
                           "not the electron total (use sum(w*n) vs N_e on that line for the count).");
    }
    const size_t nrows = static_cast<size_t>(nk) * static_cast<size_t>(nbands);
    if (nrows == 0)
    {
        log_occ_inner_line("        (empty occupations)");
        return;
    }

    std::vector<std::string> col_ik(nrows);
    std::vector<std::string> col_ib(nrows);
    std::vector<double> col_occ(nrows);
    std::vector<double> col_dn(nrows);

    size_t row = 0;
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            const int idx = ik * nbands + ib;
            col_ik[row] = std::to_string(ik);
            col_ib[row] = std::to_string(ib);
            col_occ[row] = occ_flat[idx];
            col_dn[row] = have_prev ? (occ_flat[idx] - occ_prev_for_dn[idx]) : 0.0;
            ++row;
        }
    }

    FmtTable table(/*titles=*/{"ik", "ib", "n(ik,ib)", "dn"},
                   /*nrows=*/nrows,
                   /*formats=*/{"%-6s", "%-6s", "%20.8f", "%20.8f"},
                   /*indents=*/4,
                   /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::CENTER});
    table << col_ik << col_ib << col_occ << col_dn;
    GlobalV::ofs_running << table.str() << std::endl;
}

void print_rdmft_outer_energy_stdout(bool converged, double E)
{
    GlobalV::ofs_running << std::fixed << std::setprecision(10);
    if (converged)
    {
        GlobalV::ofs_running << "  RDMFT outer: converged  E = " << E << std::defaultfloat << std::endl;
    }
    else
    {
        GlobalV::ofs_running << "  RDMFT outer: NOT converged  E = " << E << std::defaultfloat << std::endl;
    }
}

void print_rdmft_energy_table_running(const std::vector<std::string>& titles,
                                      const std::vector<double>& energies_ry)
{
    std::vector<double> energies_ev(energies_ry.size());
    std::transform(energies_ry.begin(), energies_ry.end(), energies_ev.begin(), [](const double e) {
        return e * ModuleBase::Ry_to_eV;
    });

    FmtTable table(/*titles=*/{"Energy", "Rydberg", "eV"},
                   /*nrows=*/titles.size(),
                   /*formats=*/{"%-14s", "%20.10f", "%20.10f"},
                   /*indents=*/1,
                   /*align=*/{/*value*/FmtTable::Align::LEFT, /*title*/FmtTable::Align::CENTER});
    table << titles << energies_ry << energies_ev;
    GlobalV::ofs_running << table.str() << std::endl;
}

void print_occ_table_to_stream(std::ostream& os,
                               const std::vector<double>& occ_flat,
                               int nk,
                               int nbands,
                               const std::string& header,
                               int indents)
{
    os << header << std::endl;

    const size_t nrows = static_cast<size_t>(nk) * static_cast<size_t>(nbands);
    if (nrows == 0)
    {
        os << "  (empty)" << std::endl;
        return;
    }

    std::vector<std::string> col_ik(nrows);
    std::vector<std::string> col_ib(nrows);
    std::vector<double> col_occ(nrows);

    size_t row = 0;
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            col_ik[row] = std::to_string(ik);
            col_ib[row] = std::to_string(ib);
            col_occ[row] = occ_flat[ik * nbands + ib];
            ++row;
        }
    }

    FmtTable table(/*titles=*/{"ik", "ib", "n(ik,ib)"},
                   /*nrows=*/nrows,
                   /*formats=*/{"%-6s", "%-6s", "%20.10f"},
                   /*indents=*/indents,
                   /*align=*/{/*value*/FmtTable::Align::RIGHT, /*title*/FmtTable::Align::CENTER});
    table << col_ik << col_ib << col_occ;
    os << table.str() << std::endl;
}

void print_occ_table_running(const std::vector<double>& occ_flat, int nk, int nbands)
{
    print_occ_table_to_stream(GlobalV::ofs_running, occ_flat, nk, nbands,
                              "  Occupations n(ik, ib):", 2);
}

/// Inner occ/orb iterations: log to running file only (no stdout).
/// Orbital inner iterations pass print_occupations=false to avoid repeating n(ik,ib).
void print_inner_loop_stdout(const std::string& label,
                             int inner_iter,
                             double end_energy,
                             const std::vector<double>& occ_flat,
                             int nk,
                             int nbands,
                             const bool print_occupations = true)
{
    GlobalV::ofs_running << std::fixed << std::setprecision(10)
                         << "  " << label << " " << inner_iter << " end: E = " << end_energy
                         << std::defaultfloat << std::endl;
    if (print_occupations)
    {
        print_occ_table_to_stream(GlobalV::ofs_running, occ_flat, nk, nbands,
                                  "  Occupations n(ik, ib):", 2);
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

std::string occ_param_to_string(const OccParamType t)
{
    switch (t)
    {
        case OccParamType::CosineSq: return "cosine_sq";
        case OccParamType::Logistic: return "logistic";
        case OccParamType::SigmaShift: return "sigma_shift";
    }
    return "unknown";
}

std::string occ_init_mode_to_string(const OccInitMode m)
{
    switch (m)
    {
        case OccInitMode::KS: return "ks";
        case OccInitMode::Perturbed: return "perturbed";
        case OccInitMode::Binary: return "binary";
        case OccInitMode::Uniform: return "uniform";
    }
    return "unknown";
}

std::string constraint_method_to_string(const ConstraintMethod m)
{
    switch (m)
    {
        case ConstraintMethod::AugmentedLagrangian: return "augmented_lagrangian";
        case ConstraintMethod::ProjectedGradient: return "projected_gradient";
        case ConstraintMethod::ActiveSet: return "active_set";
    }
    return "unknown";
}

std::string optimizer_to_string(const OptimizerType t)
{
    switch (t)
    {
        case OptimizerType::SteepestDescent: return "sd";
        case OptimizerType::ConjugateGradient: return "cg";
        case OptimizerType::LBFGS: return "lbfgs";
        case OptimizerType::Adam: return "adam";
    }
    return "unknown";
}

std::string strategy_to_string(const SolverStrategy s)
{
    switch (s)
    {
        case SolverStrategy::Alternating: return "alternating";
        case SolverStrategy::Joint: return "joint";
    }
    return "unknown";
}

std::string bb_mode_to_string(const BBStepMode m)
{
    switch (m)
    {
        case BBStepMode::BB1: return "bb1";
        case BBStepMode::BB2: return "bb2";
        case BBStepMode::Alternate: return "alternate";
    }
    return "unknown";
}

void emit_rdmft_config_kv_table(const char* section_title,
                                const std::vector<std::string>& keys,
                                const std::vector<std::string>& vals)
{
    if (keys.empty())
    {
        return;
    }
    GlobalV::ofs_running << "\n  --- " << section_title << " ---" << std::endl;
    FmtTable table(/*titles=*/{"RDMFT option", "value"},
                   /*nrows=*/keys.size(),
                   /*formats=*/{"%-34s", "%-40s"},
                   /*indents=*/1,
                   /*align=*/{/*value*/FmtTable::Align::LEFT, /*title*/FmtTable::Align::CENTER});
    table << keys << vals;
    GlobalV::ofs_running << table.str() << std::endl;
}

void print_rdmft_run_config(const RDMFTConfig& cfg,
                            const int nk,
                            const int nbands,
                            const double n_electrons,
                            const RDMFTNelectronTargetMeta& meta)
{
    GlobalV::ofs_running << "\n===== RDMFT Run Configuration =====" << std::endl;
    GlobalV::ofs_running << "  Path: KS -> RDMFT [strategy=" << strategy_to_string(cfg.strategy)
                         << ", constraint=" << constraint_method_to_string(cfg.constraint_method)
                         << "]" << std::endl;

    auto as_sci = [](const double v) {
        std::ostringstream os;
        os << std::scientific << std::setprecision(6) << v;
        return os.str();
    };
    auto add_kv = [](std::vector<std::string>& keys,
                     std::vector<std::string>& vals,
                     const std::string& k,
                     const std::string& v) {
        keys.push_back(k);
        vals.push_back(v);
    };

    // --- Section A: problem size and electron target ---
    {
        std::vector<std::string> keys;
        std::vector<std::string> vals;
        add_kv(keys, vals, "nk", std::to_string(nk));
        add_kv(keys, vals, "nbands", std::to_string(nbands));
        add_kv(keys, vals, "rdmft_nelec_use_input", meta.use_input_nelec ? "true" : "false");
        add_kv(keys, vals, "sum_initial_wg", as_sci(meta.sum_initial_wg));
        add_kv(keys, vals, "input_nelec", as_sci(meta.input_nelec));
        add_kv(keys, vals, "rdmft_nelec_delta", as_sci(meta.rdmft_nelec_delta));
        add_kv(keys, vals, "n_electrons_constraint_Ne", as_sci(n_electrons));
        emit_rdmft_config_kv_table("Problem size and electron constraint", keys, vals);
    }
    if (meta.use_input_nelec && std::abs(meta.sum_initial_wg - meta.input_nelec) > 1.0e-6)
    {
        GlobalV::ofs_running << "  Note: sum_initial_wg differs from input_nelec; initial KS occupations may "
                                 "violate the equality target until the optimiser updates n."
                             << std::endl;
    }

    // --- Section B: strategy ---
    {
        std::vector<std::string> keys;
        std::vector<std::string> vals;
        add_kv(keys, vals, "rdmft_solver_strategy", strategy_to_string(cfg.strategy));
        if (cfg.strategy == SolverStrategy::Alternating)
        {
            add_kv(keys, vals, "rdmft_occ_optimizer", optimizer_to_string(cfg.occ_optimizer));
            add_kv(keys, vals, "rdmft_orb_optimizer", optimizer_to_string(cfg.orb_optimizer));
            add_kv(keys, vals, "rdmft_outer_maxiter", std::to_string(cfg.outer_maxiter));
            add_kv(keys, vals, "rdmft_occ_maxiter", std::to_string(cfg.occ_maxiter));
            add_kv(keys, vals, "rdmft_orb_maxiter", std::to_string(cfg.orb_maxiter));
            add_kv(keys, vals, "rdmft_print_stiefel_gram", cfg.print_stiefel_gram ? "true" : "false");
        }
        else
        {
            add_kv(keys, vals, "rdmft_joint_optimizer", optimizer_to_string(cfg.joint_optimizer));
            add_kv(keys, vals, "joint_orb_scale", as_sci(cfg.joint_orb_scale));
            add_kv(keys, vals, "rdmft_outer_maxiter", std::to_string(cfg.outer_maxiter));
        }
        add_kv(keys, vals, "rdmft_orb_retraction", orb_retraction_to_string(cfg.orb_retraction));
        emit_rdmft_config_kv_table("Solver strategy", keys, vals);
    }

    // --- Section C: constraint method ---
    {
        std::vector<std::string> keys;
        std::vector<std::string> vals;
        add_kv(keys, vals, "rdmft_constraint", constraint_method_to_string(cfg.constraint_method));
        if (cfg.constraint_method == ConstraintMethod::AugmentedLagrangian)
        {
            add_kv(keys, vals, "rdmft_occ_tol", as_sci(cfg.rdmft_occ_tol));
            add_kv(keys, vals, "rdmft_occ_energy_tol", as_sci(cfg.occ_energy_tol));
            add_kv(keys, vals, "alm_lambda_init", as_sci(cfg.aug_lag_lambda_init));
            add_kv(keys, vals, "alm_mu_init", as_sci(cfg.aug_lag_mu_init));
            add_kv(keys, vals, "alm_mu_factor", as_sci(cfg.aug_lag_mu_factor));
            add_kv(keys, vals, "alm_mu_max", as_sci(cfg.aug_lag_mu_max));
        }
        else
        {
            add_kv(keys, vals, "rdmft_occ_grad_tol", as_sci(cfg.occ_grad_tol));
        }
        emit_rdmft_config_kv_table("Electron-number constraint method", keys, vals);
    }

    // --- Occupation line search / BB (ALM only; SPG manages its own BB step) ---
    {
        std::vector<std::string> keys;
        std::vector<std::string> vals;
        add_kv(keys, vals, "rdmft_alpha_step", as_sci(cfg.line_search_alpha_init));
        add_kv(keys, vals, "alm_bb_enabled", cfg.alm_bb_enabled ? "true" : "false");
        add_kv(keys, vals, "alm_bb_mode", bb_mode_to_string(cfg.alm_bb_mode));
        add_kv(keys, vals, "alm_bb_alpha_min", as_sci(cfg.alm_bb_alpha_min));
        add_kv(keys, vals, "alm_bb_alpha_max", as_sci(cfg.alm_bb_alpha_max));
        emit_rdmft_config_kv_table("Occupation line search and BB seed", keys, vals);
    }

    // --- Section D: occupation model ---
    {
        std::vector<std::string> keys;
        std::vector<std::string> vals;
        add_kv(keys, vals, "rdmft_occ_param", occ_param_to_string(cfg.occ_param));
        add_kv(keys, vals, "rdmft_occ_init_mode", occ_init_mode_to_string(cfg.occ_init_mode));
        add_kv(keys, vals, "rdmft_occ_init_nbands_top", std::to_string(cfg.occ_init_nbands_top));
        add_kv(keys, vals, "rdmft_occ_init_perturb", as_sci(cfg.occ_init_perturb));
        if (cfg.xc_type == XCFunctionalType::HF)
        {
            add_kv(keys, vals, "rdmft_occ_entropy_gamma", as_sci(cfg.occ_entropy_gamma));
        }
        emit_rdmft_config_kv_table("Occupation model and initialisation", keys, vals);
    }

    // --- Section E: global tolerances and line search ---
    {
        std::vector<std::string> keys;
        std::vector<std::string> vals;
        add_kv(keys, vals, "rdmft_energy_tol", as_sci(cfg.energy_tol));
        add_kv(keys, vals, "rdmft_orb_grad_tol", as_sci(cfg.orb_grad_tol));
        add_kv(keys, vals, "rdmft_orb_energy_tol", as_sci(cfg.orb_energy_tol));
        add_kv(keys, vals, "line_search_c1", as_sci(cfg.line_search_c1));
        add_kv(keys, vals, "line_search_rho", as_sci(cfg.line_search_rho));
        add_kv(keys, vals, "line_search_max_iter", std::to_string(cfg.line_search_max_iter));
        add_kv(keys, vals, "line_search_polynomial", cfg.line_search_polynomial ? "true" : "false");
        add_kv(keys, vals, "line_search_c2", as_sci(cfg.line_search_c2));
        add_kv(keys, vals, "line_search_max_zoom", std::to_string(cfg.line_search_max_zoom));
        add_kv(keys, vals, "rdmft_grad_check", cfg.grad_check ? "true" : "false");
        emit_rdmft_config_kv_table("Global tolerances and line search", keys, vals);
    }

    // --- Section F: optimiser hyperparameters (active optimisers only) ---
    {
        std::vector<std::string> keys;
        std::vector<std::string> vals;
        bool need_lbfgs = false;
        bool need_adam = false;
        if (cfg.strategy == SolverStrategy::Joint)
        {
            need_lbfgs = (cfg.joint_optimizer == OptimizerType::LBFGS);
            need_adam = (cfg.joint_optimizer == OptimizerType::Adam);
        }
        else
        {
            need_lbfgs = (cfg.occ_optimizer == OptimizerType::LBFGS) || (cfg.orb_optimizer == OptimizerType::LBFGS);
            need_adam = (cfg.occ_optimizer == OptimizerType::Adam) || (cfg.orb_optimizer == OptimizerType::Adam);
        }
        if (need_lbfgs)
        {
            add_kv(keys, vals, "lbfgs_memory", std::to_string(cfg.lbfgs_memory));
        }
        if (need_adam)
        {
            add_kv(keys, vals, "adam_lr", as_sci(cfg.adam_lr));
        }
        emit_rdmft_config_kv_table("Optimiser hyperparameters", keys, vals);
    }
}

void print_rdmft_optimization_summary_alternating(const int outer_iters,
                                                  const int occ_calls,
                                                  const int orb_calls,
                                                  const int occ_inner_total,
                                                  const int orb_inner_total,
                                                  const double occ_time_sec,
                                                  const double orb_time_sec,
                                                  const double total_time_sec)
{
    GlobalV::ofs_running << "\n===== RDMFT Optimization Summary =====" << std::endl;
    std::vector<std::string> block;
    std::vector<std::string> calls;
    std::vector<std::string> inner_iters;
    std::vector<double> wall_time_sec;
    std::vector<double> avg_time_per_call;

    block.push_back("occupation");
    calls.push_back(std::to_string(occ_calls));
    inner_iters.push_back(std::to_string(occ_inner_total));
    wall_time_sec.push_back(occ_time_sec);
    avg_time_per_call.push_back(occ_calls > 0 ? occ_time_sec / static_cast<double>(occ_calls) : 0.0);

    block.push_back("orbital");
    calls.push_back(std::to_string(orb_calls));
    inner_iters.push_back(std::to_string(orb_inner_total));
    wall_time_sec.push_back(orb_time_sec);
    avg_time_per_call.push_back(orb_calls > 0 ? orb_time_sec / static_cast<double>(orb_calls) : 0.0);

    block.push_back("total");
    calls.push_back(std::to_string(outer_iters));
    inner_iters.push_back("-");
    wall_time_sec.push_back(total_time_sec);
    avg_time_per_call.push_back(outer_iters > 0 ? total_time_sec / static_cast<double>(outer_iters) : 0.0);

    FmtTable table(/*titles=*/{"Block", "Calls", "Inner iters", "Wall time (s)", "Avg/call (s)"},
                   /*nrows=*/block.size(),
                   /*formats=*/{"%-14s", "%-10s", "%-14s", "%16.6f", "%16.6f"},
                   /*indents=*/1,
                   /*align=*/{/*value*/FmtTable::Align::LEFT, /*title*/FmtTable::Align::CENTER});
    table << block << calls << inner_iters << wall_time_sec << avg_time_per_call;
    GlobalV::ofs_running << table.str() << std::endl;
}

void print_rdmft_optimization_summary_joint(const int outer_iters,
                                            const double total_time_sec)
{
    GlobalV::ofs_running << "\n===== RDMFT Optimization Summary =====" << std::endl;
    GlobalV::ofs_running << "  Strategy=joint: occupations and orbitals are optimized simultaneously." << std::endl;
    std::vector<std::string> keys = {"joint_outer_iters", "joint_total_time_s", "joint_avg_time_per_outer_s"};
    std::vector<std::string> vals;
    vals.push_back(std::to_string(outer_iters));
    {
        std::ostringstream os;
        os << std::fixed << std::setprecision(6) << total_time_sec;
        vals.push_back(os.str());
    }
    {
        std::ostringstream os;
        os << std::fixed << std::setprecision(6)
           << (outer_iters > 0 ? total_time_sec / static_cast<double>(outer_iters) : 0.0);
        vals.push_back(os.str());
    }
    FmtTable table(/*titles=*/{"Metric", "value"},
                   /*nrows=*/keys.size(),
                   /*formats=*/{"%-28s", "%-24s"},
                   /*indents=*/1,
                   /*align=*/{/*value*/FmtTable::Align::LEFT, /*title*/FmtTable::Align::CENTER});
    table << keys << vals;
    GlobalV::ofs_running << table.str() << std::endl;
}

} // namespace

// Helper to get |x|^2 for both real and complex types
inline double abs2(double x) { return x * x; }
inline double abs2(std::complex<double> x) { return std::norm(x); }

// ------------------------------------------------------------------------
// Flatten / unflatten helpers for the Stiefel-manifold orbital optimiser.
// To reuse the Euclidean lbfgs / Adam optimiser working on a flat
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
// but the lbfgs preconditioning is an approximation.  This matches the
// standard "Riemannian lbfgs by projection" prescription used by the
// ROPTLITE library and is adequate for the current applications.
// ------------------------------------------------------------------------
inline void psi_to_flat(const psi::Psi<double>& P, std::vector<double>& flat)
{
    const int nk = P.get_nk();
    const int nb = P.get_nbands();
    const int nbs = P.get_nbasis();
    const int n = nk * nb * nbs;
    flat.resize(n);
    if (n == 0)
    {
        return;
    }
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
    if (n == 0)
    {
        return;
    }
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
    if (n == 0)
    {
        return;
    }
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
    if (n == 0)
    {
        return;
    }
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
    double n_electrons,
    const RDMFTNelectronTargetMeta& nelec_meta)
{
    if (kv == nullptr)
    {
        throw std::invalid_argument("RDMFTSolver::init: kv pointer is null");
    }
    if (nbands <= 0)
    {
        throw std::invalid_argument("RDMFTSolver::init: nbands must be positive, got "
                                    + std::to_string(nbands));
    }
    if (n_electrons <= 0.0)
    {
        throw std::invalid_argument("RDMFTSolver::init: n_electrons must be positive, got "
                                    + std::to_string(n_electrons));
    }

    config_ = config;
    energy_grad_ = &energy_grad;
    kv_ = kv;
    nbands_ = nbands;
    n_electrons_ = n_electrons;
    nelec_meta_ = nelec_meta;
    if (nelec_meta_.sum_initial_wg == 0.0 && nelec_meta_.input_nelec == 0.0 && !nelec_meta_.use_input_nelec
        && nelec_meta_.rdmft_nelec_delta == 0.0)
    {
        nelec_meta_.sum_initial_wg = n_electrons;
        nelec_meta_.input_nelec = n_electrons;
    }
    nk_ = energy_grad.nk();

    if (nk_ <= 0)
    {
        throw std::invalid_argument("RDMFTSolver::init: energy_grad.nk() must be positive, got "
                                    + std::to_string(nk_));
    }

    if (config_.occ_param == OccParamType::SigmaShift && config_.strategy == SolverStrategy::Alternating)
    {
        throw std::invalid_argument(
            "RDMFTSolver::init: rdmft_occ_param sigma_shift requires rdmft_solver_strategy joint "
            "(alternating occupation optimisation does not apply the sigma-shift lambda solve).");
    }

    // Initialize occupation parameterization
    occ_param_ = std::make_unique<OccupationParam>(config_.occ_param);

    // Build k-point weight vector
    std::vector<double> kweights(nk_);
    for (int ik = 0; ik < nk_; ++ik)
        kweights[ik] = kv_->wk[ik];

    occ_constraint_ = std::make_unique<OccupationConstraint>(
        config_.constraint_method, n_electrons_, kweights, nbands_);

    occ_optimizer_ = std::make_unique<EuclideanOptimizer>(config_.occ_optimizer, config_);
    orb_optimizer_ = std::make_unique<EuclideanOptimizer>(config_.orb_optimizer, config_);

    // Forward the orbital retraction selector to EnergyGradient. Default
    // Polar reproduces the historical behaviour byte-for-byte; QR / Cayley
    // are honoured by retract_orbitals (serial-only; MPI builds fall back
    // to Polar with a one-time warning, see EnergyGradient::retract_orbitals).
    energy_grad_->set_orb_retraction(config_.orb_retraction);

    occ_constraint_->set_lambda(config_.aug_lag_lambda_init);
    occ_constraint_->set_mu(config_.aug_lag_mu_init);

    // Sigma-shift parameterization: electron constraint is satisfied exactly
    // each step via bisection; no augmented-Lagrangian penalty needed.
    if (config_.occ_param == OccParamType::SigmaShift)
    {
        sigma_shift_param_ = std::make_unique<SigmaShiftOccParam>(
            n_electrons_, kweights, nbands_);
    }
}

template <typename TK, typename TR>
double RDMFTSolver<TK, TR>::solve(
    std::vector<double>& occ_flat,
    psi::Psi<TK>& wfc)
{
    ModuleBase::timer::start("RDMFT", "solve");
    last_result_ = {};
    logged_initial_occ_for_first_occ_gradient_ = false;

    print_rdmft_run_config(config_, nk_, nbands_, n_electrons_, nelec_meta_);

    // occ_flat on entry is the KS occupation seed from pelec->wg.
    const std::vector<double> occ_ks_seed = occ_flat;
    double sum_k_weights = 0.0;
    for (int jk = 0; jk < nk_; ++jk)
    {
        sum_k_weights += occ_constraint_->kweights()[jk];
    }
    sum_k_weights = std::max(sum_k_weights, 1e-20);
    const double c_ks_seed = occ_constraint_->constraint_violation(occ_ks_seed);

    auto log_init_common = [&](const std::string& mode_name,
                               const std::string& details) {
        GlobalV::ofs_running << "RDMFT init occupations: mode=" << mode_name
                             << ", constraint_ks=" << std::scientific << c_ks_seed
                             << ", " << details << std::defaultfloat << std::endl;
    };

    if (config_.occ_init_mode == OccInitMode::KS)
    {
        occ_flat = occ_ks_seed;
        log_init_common("ks", "note=using KS occupations directly");
    }
    else if (config_.occ_init_mode == OccInitMode::Perturbed)
    {
        occ_flat = occ_ks_seed;
        const int K = std::min(config_.occ_init_nbands_top, nbands_);
        const double delta = config_.occ_init_perturb;

        if (K <= 0 || delta <= 0.0)
        {
            log_init_common("perturbed",
                            "note=K<=0 or delta<=0, fallback=ks");
        }
        else
        {
            int n_plus = 0;
            int n_minus = 0;
            std::vector<bool> touched(static_cast<size_t>(nk_ * nbands_), false);
            for (int ik = 0; ik < nk_; ++ik)
            {
                const double wk = occ_constraint_->kweights()[ik];
                const double n_target_k = n_electrons_ * wk / sum_k_weights;
                const int ib_fermi = find_init_fermi_boundary(occ_ks_seed, ik, nbands_, wk, n_target_k);
                for (int t = 0; t < K; ++t)
                {
                    const int ib_above = ib_fermi + 1 + t;
                    if (ib_above >= 0 && ib_above < nbands_)
                    {
                        const int idx = ik * nbands_ + ib_above;
                        occ_flat[idx] = occ_ks_seed[idx] + delta;
                        touched[idx] = true;
                        ++n_plus;
                    }

                    const int ib_below = ib_fermi - t;
                    if (ib_below >= 0 && ib_below < nbands_)
                    {
                        const int idx = ik * nbands_ + ib_below;
                        occ_flat[idx] = occ_ks_seed[idx] - delta;
                        touched[idx] = true;
                        ++n_minus;
                    }
                }
            }

            const double c_before_project = occ_constraint_->constraint_violation(occ_flat);
            occ_constraint_->project_preserving_ks_on_fixed(occ_flat, occ_ks_seed, touched);
            const double c_after_project = occ_constraint_->constraint_violation(occ_flat);
            const double sum_abs_init_change = sum_abs_diff(occ_flat, occ_ks_seed);

            std::ostringstream os;
            os << "K=" << K
               << ", delta=" << delta
               << ", n_plus=" << n_plus
               << ", n_minus=" << n_minus
               << ", sum|n_init-ks|=" << sum_abs_init_change
               << ", constraint_before_project=" << c_before_project
               << ", constraint_after_project=" << c_after_project;
            log_init_common("perturbed", os.str());
        }
    }
    else if (config_.occ_init_mode == OccInitMode::Binary)
    {
        occ_flat = occ_ks_seed;
        const int K = std::min(config_.occ_init_nbands_top, nbands_);
        const double delta = config_.occ_init_perturb;

        if (K <= 0 || delta <= 0.0)
        {
            log_init_common("binary",
                            "note=K<=0 or delta<=0, fallback=ks");
        }
        else
        {
            int n_high = 0;
            int n_low = 0;
            std::vector<bool> touched(static_cast<size_t>(nk_ * nbands_), false);
            for (int ik = 0; ik < nk_; ++ik)
            {
                const double wk = occ_constraint_->kweights()[ik];
                const double n_target_k = n_electrons_ * wk / sum_k_weights;
                const int ib_fermi = find_init_fermi_boundary(occ_ks_seed, ik, nbands_, wk, n_target_k);
                for (int t = 0; t < K; ++t)
                {
                    const int ib_above = ib_fermi + 1 + t;
                    if (ib_above >= 0 && ib_above < nbands_)
                    {
                        const int idx = ik * nbands_ + ib_above;
                        occ_flat[idx] = delta;
                        touched[idx] = true;
                        ++n_low;
                    }

                    const int ib_below = ib_fermi - t;
                    if (ib_below >= 0 && ib_below < nbands_)
                    {
                        const int idx = ik * nbands_ + ib_below;
                        occ_flat[idx] = 1.0 - delta;
                        touched[idx] = true;
                        ++n_high;
                    }
                }
            }

            const double c_before_project = occ_constraint_->constraint_violation(occ_flat);
            occ_constraint_->project_preserving_ks_on_fixed(occ_flat, occ_ks_seed, touched);
            const double c_after_project = occ_constraint_->constraint_violation(occ_flat);
            const double sum_abs_init_change = sum_abs_diff(occ_flat, occ_ks_seed);

            std::ostringstream os;
            os << "K=" << K
               << ", delta=" << delta
               << ", n_high=" << n_high
               << ", n_low=" << n_low
               << ", sum|n_init-ks|=" << sum_abs_init_change
               << ", constraint_before_project=" << c_before_project
               << ", constraint_after_project=" << c_after_project;
            log_init_common("binary", os.str());
        }
    }
    else
    {
        occ_flat = occ_ks_seed;
        const int K = std::min(config_.occ_init_nbands_top, nbands_);

        if (K <= 0)
        {
            log_init_common("uniform",
                            "note=K<=0, fallback=ks");
        }
        else
        {
            int n_uniformized = 0;
            std::vector<bool> touched(static_cast<size_t>(nk_ * nbands_), false);
            for (int ik = 0; ik < nk_; ++ik)
            {
                const double wk = occ_constraint_->kweights()[ik];
                const double n_target_k = n_electrons_ * wk / sum_k_weights;
                const int ib_fermi = find_init_fermi_boundary(occ_ks_seed, ik, nbands_, wk, n_target_k);
                std::vector<int> selected;
                selected.reserve(2 * K);

                for (int t = 0; t < K; ++t)
                {
                    const int ib_above = ib_fermi + 1 + t;
                    if (ib_above >= 0 && ib_above < nbands_)
                    {
                        selected.push_back(ib_above);
                    }
                    const int ib_below = ib_fermi - t;
                    if (ib_below >= 0 && ib_below < nbands_)
                    {
                        selected.push_back(ib_below);
                    }
                }

                if (selected.empty())
                {
                    continue;
                }

                double n_top = 0.0;
                for (const int ib : selected)
                {
                    n_top += occ_flat[ik * nbands_ + ib];
                }
                const double n_uniform = n_top / static_cast<double>(selected.size());
                for (const int ib : selected)
                {
                    const int idx = ik * nbands_ + ib;
                    occ_flat[idx] = n_uniform;
                    touched[idx] = true;
                }
                n_uniformized += static_cast<int>(selected.size());
            }

            const double c_before_project = occ_constraint_->constraint_violation(occ_flat);
            occ_constraint_->project_preserving_ks_on_fixed(occ_flat, occ_ks_seed, touched);
            const double c_after_project = occ_constraint_->constraint_violation(occ_flat);
            const double sum_abs_init_change = sum_abs_diff(occ_flat, occ_ks_seed);

            std::ostringstream os;
            os << "K=" << K
               << ", n_uniformized=" << n_uniformized
               << ", sum|n_init-ks|=" << sum_abs_init_change
               << ", constraint_before_project=" << c_before_project
               << ", constraint_after_project=" << c_after_project;
            log_init_common("uniform", os.str());
        }
    }

    // Precompute S_k = U_k^H U_k and switch to the X_k = U_k C_k variable.
    // All RDMFT manifold operations are performed only in X-space.
    energy_grad_->precompute_cholesky_S();

    energy_grad_->wfc_C_to_X(wfc);

    // Orthonormalize in X-space so all subsequent manifold operations start
    // from a valid point on the standard Stiefel manifold X^H X = I.
    energy_grad_->retract_orbitals(wfc, wfc, 0.0);

    // Sigma-shift: map KS occupations to n = σ(z+λ) with λ enforcing Σ w n = N_e
    // before gradient check or strategy entry (KS wg/wk may not satisfy N_e exactly).
    if (config_.occ_param == OccParamType::SigmaShift && sigma_shift_param_ != nullptr)
    {
        std::vector<double> z_pre(occ_flat.size());
        sigma_shift_param_->occ_to_params(occ_flat, z_pre);
        const double lam_pre = sigma_shift_param_->compute_lambda(z_pre);
        sigma_shift_param_->params_to_occ(z_pre, lam_pre, occ_flat);
    }

    // Gradient check must run only after X-space setup above; compute()/retract
    // paths require precompute_cholesky_S() (see EnergyGradient::wfc_X_to_C).
    if (config_.grad_check)
    {
        (void)check_gradient_consistency(occ_flat, wfc);
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
    energy_grad_->wfc_X_to_C(wfc);

    ModuleBase::timer::end("RDMFT", "solve");
    return E;
}

template <typename TK, typename TR>
double RDMFTSolver<TK, TR>::solve_alternating(
    std::vector<double>& occ_flat,
    psi::Psi<TK>& wfc)
{
    GlobalV::ofs_running << "\n===== RDMFT Alternating Optimization =====" << std::endl;
    const auto t_alternating_start = std::chrono::steady_clock::now();

    double E_prev = 1e30;
    double E = 0.0;

    std::vector<double> occ_at_outer_start;
    double last_dE = 0.0;
    double last_abs_c = 0.0;
    double last_occ_gnorm = 0.0;
    double last_orb_gnorm = 0.0;
    int outer_iters_done = 0;
    int occ_calls = 0;
    int orb_calls = 0;
    int occ_inner_total = 0;
    int orb_inner_total = 0;
    double occ_time_sec = 0.0;
    double orb_time_sec = 0.0;

    for (int iter = 0; iter < config_.outer_maxiter; ++iter)
    {
        occ_at_outer_start.assign(occ_flat.begin(), occ_flat.end());

        if (config_.print_stiefel_gram)
        {
            std::vector<double> stiefel_gram_frob_per_ik;
            energy_grad_->stiefel_gram_residual_frobenius_per_k(wfc, stiefel_gram_frob_per_ik);
            GlobalV::ofs_running << "    Stiefel Gram residual ||G_k - I||_F per k-point (G_k = X_k^H X_k):"
                << std::endl;
            for (int ik = 0; ik < wfc.get_nk(); ++ik)
            {
                GlobalV::ofs_running << std::fixed << std::setprecision(10) << "      ik=" << ik
                    << "  ||G_k - I||_F = " << stiefel_gram_frob_per_ik[ik] << std::endl;
            }
        }

        // 1. Optimize occupations with orbitals fixed (every outer iteration).
        OptResult occ_result;
        const auto t_occ0 = std::chrono::steady_clock::now();
        occ_result = optimize_occupations(occ_flat, wfc);
        const auto t_occ1 = std::chrono::steady_clock::now();
        occ_time_sec += std::chrono::duration<double>(t_occ1 - t_occ0).count();
        ++occ_calls;
        occ_inner_total += occ_result.iterations;
        const double occ_outer_dn_sum = sum_abs_diff(occ_flat, occ_at_outer_start);
        GlobalV::ofs_running << "    occ inner: " << occ_result.iterations << " iters, gnorm="
            << std::scientific << occ_result.grad_norm
            << "  E=" << std::fixed << std::setprecision(10) << occ_result.final_energy
            << (occ_result.converged ? "  (converged)" : "")
            << "  sum|dn|_outer=" << std::scientific << occ_outer_dn_sum
            << std::defaultfloat << std::endl;

        // 2. Optimize orbitals with occupations fixed (every outer iteration).
        const auto t_orb0 = std::chrono::steady_clock::now();
        OptResult orb_result = optimize_orbitals(occ_flat, wfc);
        const auto t_orb1 = std::chrono::steady_clock::now();
        orb_time_sec += std::chrono::duration<double>(t_orb1 - t_orb0).count();
        ++orb_calls;
        orb_inner_total += orb_result.iterations;
        GlobalV::ofs_running << "    orb inner: " << orb_result.iterations << " iters, gnorm="
            << std::scientific << orb_result.grad_norm
            << "  E=" << std::fixed << std::setprecision(10) << orb_result.final_energy
            << (orb_result.converged ? "  (converged)" : "") << std::endl;

        // 3. Evaluate full energy
        std::vector<double> grad_occ;
        psi::Psi<TK> grad_wfc;
        E = energy_grad_->compute(occ_flat, wfc, grad_occ, grad_wfc);
        const double E_one = energy_grad_->E_one_body();
        const double E_hartree = energy_grad_->E_hartree();
        const double E_xc = energy_grad_->E_xc();
        const double E_entropy = energy_grad_->E_entropy();
        const double E_ewald = energy_grad_->E_ewald();
        const double E_total = energy_grad_->E_total();

        double dE = std::abs(E - E_prev);
        double constraint_viol = occ_constraint_->constraint_violation(occ_flat);

        GlobalV::ofs_running << std::fixed << std::setprecision(10)
            << "  RDMFT iter " << iter + 1
            << "  E = " << E
            << "  dE = " << std::scientific << dE
            << "  |c| = " << std::abs(constraint_viol)
            << "  occ_gnorm = " << occ_result.grad_norm
            << "  orb_gnorm = " << orb_result.grad_norm
            << "  occ_conv = " << (occ_result.converged ? "Y" : "N")
            << "  orb_conv = " << (orb_result.converged ? "Y" : "N")
            << std::endl;

        if (config_.occ_entropy_gamma > 0.0 && config_.xc_type == XCFunctionalType::HF)
        {
            print_rdmft_energy_table_running({"E_one_elec", "E_Hartree", "E_xc", "E_entropy", "E_Ewald", "E_total"},
                                             {E_one, E_hartree, E_xc, E_entropy, E_ewald, E_total});
        }
        else
        {
            print_rdmft_energy_table_running({"E_one_elec", "E_Hartree", "E_xc", "E_Ewald", "E_total"},
                                             {E_one, E_hartree, E_xc, E_ewald, E_total});
        }

        // stdout: one block per outer iteration (inner loops stay off stdout).
        std::cout << std::fixed << std::setprecision(10)
                  << "  RDMFT outer iter " << (iter + 1) << "  E_total=" << E_total << "  (Ry)";
        if (iter > 0)
        {
            std::cout << std::scientific << "  dE=" << (E - E_prev) << std::fixed;
        }
        std::cout << std::defaultfloat << std::endl;

        last_dE = dE;
        last_abs_c = std::abs(constraint_viol);
        last_occ_gnorm = occ_result.grad_norm;
        last_orb_gnorm = orb_result.grad_norm;
        outer_iters_done = iter + 1;

        // Outer stop: energy change below tol (if tol > 0) or both inner sub-problems converged.
        const bool inner_both = occ_result.converged && orb_result.converged;
        const bool energy_ok
            = (config_.energy_tol > 0.0) && (dE < config_.energy_tol);
        const bool outer_converged = (iter > 0) && (energy_ok || inner_both);
        if (outer_converged)
        {
            last_result_.converged = true;
            last_result_.iterations = iter + 1;
            last_result_.final_energy = E;
            last_result_.grad_norm = std::max(occ_result.grad_norm, orb_result.grad_norm);
            if (energy_ok && inner_both)
            {
                GlobalV::ofs_running << "  RDMFT alternating: outer loop stopped (|dE| < rdmft_energy_tol and "
                                        "OCC&ORB inner converged)"
                                     << std::endl;
            }
            else if (energy_ok)
            {
                GlobalV::ofs_running << "  RDMFT alternating: outer loop stopped (|dE| < rdmft_energy_tol)"
                                     << std::endl;
            }
            else
            {
                GlobalV::ofs_running << "  RDMFT alternating: outer loop stopped (OCC&ORB inner converged)"
                                     << std::endl;
            }
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
    // Always print final occupations in tabular form.
    print_occ_table_running(occ_flat, nk_, nbands_);
    const auto t_alternating_end = std::chrono::steady_clock::now();
    const double total_time_sec = std::chrono::duration<double>(t_alternating_end - t_alternating_start).count();
    print_rdmft_optimization_summary_alternating(outer_iters_done,
                                                 occ_calls,
                                                 orb_calls,
                                                 occ_inner_total,
                                                 orb_inner_total,
                                                 occ_time_sec,
                                                 orb_time_sec,
                                                 total_time_sec);
    if (!last_result_.converged)
    {
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
    const auto t_joint_start = std::chrono::steady_clock::now();

    // Sigma-shift active flag and mutable shift λ (updated each outer iteration).
    // If the enum is SigmaShift, sigma_shift_param_ must have been created in init().
    assert(config_.occ_param != OccParamType::SigmaShift || sigma_shift_param_ != nullptr);
    const bool use_sigma_shift = (config_.occ_param == OccParamType::SigmaShift);
    double sigma_lambda = 0.0;

    // Convert occupations to unconstrained parameters.
    std::vector<double> params(occ_flat.size());
    if (use_sigma_shift)
    {
        sigma_shift_param_->occ_to_params(occ_flat, params);
        // Fold the first implicit λ into z: with z' = z + λ, n_i = σ(z'_i) and
        // Σ_k w_k Σ_i n_i = N_e still holds, while the next λ solve is ~0.  This
        // avoids a large shift from logit(KS) vs the target electron sum and
        // keeps early joint iterations numerically on the electron constraint.
        const double lam_fold = sigma_shift_param_->compute_lambda(params);
        for (double& z : params)
            z += lam_fold;
        const double lam_canon = sigma_shift_param_->compute_lambda(params);
        sigma_shift_param_->params_to_occ(params, lam_canon, occ_flat);
    }
    else
    {
        occ_param_->occ_to_params(occ_flat, params);
        occ_param_->params_to_occ(params, occ_flat);
    }

    const int n_occ_params = static_cast<int>(params.size());
    const int nk = wfc.get_nk();
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();
    const int orb_total_size = nk * nb_local * nbs_local;

    // Flat-size of the orbital block in the packed vector: one real per
    // orbital coefficient when TK = double, two reals (Re, Im) when TK is
    // complex.
    int orb_flat_size = orb_total_size;
    if (!std::is_same<TK, double>::value) orb_flat_size *= 2;

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
    // After a joint line-search failure we reset the optimiser (next iter is
    // effectively steepest descent). If line search fails again on that
    // recovery iteration, stop the outer loop.
    bool joint_ls_failed_prev = false;

    auto total_energy = [&](const std::vector<double>& occ_in,
                            const psi::Psi<TK>& wfc_in) -> double {
        double E_val = energy_grad_->compute_energy(occ_in,
                                const_cast<psi::Psi<TK>&>(wfc_in));
        if (!use_sigma_shift && config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
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
        if (use_sigma_shift)
        {
            // Sigma-shift: determine λ by bisection, then n = σ(z+λ).
            sigma_lambda = sigma_shift_param_->compute_lambda(params);
            sigma_shift_param_->params_to_occ(params, sigma_lambda, occ_flat);
        }
        else
        {
            occ_param_->params_to_occ(params, occ_flat);
        }
        const std::vector<double> occ_before_joint_step(occ_flat);
        occ_at_outer_start.assign(occ_flat.begin(), occ_flat.end());

        // SCF stdout / running: occupations at the beginning of this joint outer iteration.
        {
            std::ostringstream hdr_start;
            if (iter == 0 && !logged_initial_occ_for_first_occ_gradient_)
            {
                GlobalV::ofs_running << "\n  RDMFT joint: first ∂E/∂n is evaluated at the occupations in this table "
                                        "(fractional n after init / param sync)."
                                     << std::endl;
                logged_initial_occ_for_first_occ_gradient_ = true;
            }
            hdr_start << std::fixed << std::setprecision(10)
                      << "  RDMFT joint iter " << (iter + 1)
                      << "  occupations n(ik,ib) at iteration START"
                      << "  sum(w*n)=" << occ_constraint_->weighted_occupation_sum(occ_flat)
                      << "  N_e=" << n_electrons_;
            if (use_sigma_shift)
            {
                hdr_start << "  sigma_lambda=" << sigma_lambda;
            }
            print_occ_table_to_stream(GlobalV::ofs_running, occ_flat, nk_, nbands_, hdr_start.str(), 2);
        }

        // Full energy + Euclidean gradients at the current point.
        std::vector<double> grad_occ;
        psi::Psi<TK> grad_wfc;
        E = energy_grad_->compute(occ_flat, wfc, grad_occ, grad_wfc);

        // Add augmented Lagrangian penalty (skipped for sigma-shift).
        if (!use_sigma_shift && config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
        {
            E += occ_constraint_->augmented_lagrangian_penalty(occ_flat);
            std::vector<double> penalty_grad;
            occ_constraint_->augmented_lagrangian_gradient(occ_flat, penalty_grad);
            for (size_t i = 0; i < grad_occ.size(); ++i)
                grad_occ[i] += penalty_grad[i];
        }

        // Transform occupation gradient dE/dn -> dE/dp (chain rule through the
        // parameterisation: cosine^2, logistic, or sigma-shift).
        std::vector<double> grad_params;
        if (use_sigma_shift)
            sigma_shift_param_->transform_gradient_batch(grad_occ, params, sigma_lambda, grad_params);
        else
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
        // Use Euclidean (plain) norm of the projected gradient G_R
        double orb_gnorm2 = 0.0;
        {
            const int nk_local = grad_wfc.get_nk();
            const int nb_local = grad_wfc.get_nbands();
            const int nbs_local = grad_wfc.get_nbasis();
            for (int ik = 0; ik < nk_local; ++ik)
            {
                const int nelem = nb_local * nbs_local;
                if (nelem == 0)
                {
                    continue;
                }
                const TK* gk = &grad_wfc(ik, 0, 0);
                for (int i = 0; i < nelem; ++i)
                    orb_gnorm2 += std::real(std::conj(gk[i]) * gk[i]);
            }
            Parallel_Reduce::reduce_all(orb_gnorm2);
        }

        // Project the orbital search direction back onto the tangent space
        // at C. The optimiser's Euclidean update (especially for lbfgs /
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
            if (orb_total_size > 0)
            {
                TK* q = &neg_orb_dir(0, 0, 0);
                const TK* p = &orb_dir(0, 0, 0);
                for (int i = 0; i < orb_total_size; ++i) q[i] = -p[i];
            }
        }

        // Decide whether the orbital block actually needs to be stepped in
        // this iteration. Mirroring the alternating strategy's
        // optimize_orbitals (which early-exits on ||G_R|| < orb_grad_tol and
        // therefore never invokes retract_orbitals with a zero / tiny
        // direction), we skip the orbital retraction when the Riemannian
        // orbital gradient is below orb_grad_tol. Without this guard the
        // Cholesky-QR S-orthonormalisation inside retract_orbitals is
        // executed every iteration, which is a numerically non-trivial
        // O(nbasis^3) update whose round-off amplifies when the trial
        // step is otherwise supposed to leave the orbitals unchanged; on
        // H2 (where the KS seed is already a fixed point of the orbital
        // sub-problem, so ||G_R|| == 0) this spurious re-orthonormalisation
        // perturbs C by a few percent at the start and destabilises the
        // occupation line search.
        const double orb_gnorm = std::sqrt(std::max(0.0, orb_gnorm2));
        const bool step_orbitals = (orb_gnorm >= config_.orb_grad_tol);

        // Snapshot the current orbitals for line-search rollback (needed
        // only when we will actually retract).
        if (step_orbitals)
        {
            for (int ik = 0; ik < nk; ++ik)
                for (int ib = 0; ib < nb_local; ++ib)
                    for (int mu = 0; mu < nbs_local; ++mu)
                        wfc_save(ik, ib, mu) = wfc(ik, ib, mu);
        }

        // Joint line search: Strong Wolfe for lbfgs only; Armijo for sd / cg / adam
        // (see line_search_uses_strong_wolfe).
        const double alpha_init = (joint_is_lbfgs || joint_is_adam)
                                      ? 1.0
                                      : config_.line_search_alpha_init;
        double alpha = alpha_init;
        const double joint_ls_alpha0 = alpha_init;
        const double c1 = config_.line_search_c1;
        const double rho = config_.line_search_rho;
        double E_new = E;
        bool ls_success = false;
        int joint_ls_trials = 0;
        int joint_n_feval = 0;
        const bool joint_use_sw = line_search_uses_strong_wolfe(joint_type);

        std::vector<double> occ_flat_new;
        std::vector<double> params_new(params.size());

        auto eval_at_step = [&](double step,
                                std::vector<double>* grad_params_out,
                                psi::Psi<TK>* grad_wfc_out,
                                double* dd_out) -> double
        {
            if (step_orbitals)
            {
                for (int ik = 0; ik < nk; ++ik)
                    for (int ib = 0; ib < nb_local; ++ib)
                        for (int mu = 0; mu < nbs_local; ++mu)
                            wfc(ik, ib, mu) = wfc_save(ik, ib, mu);
                energy_grad_->retract_orbitals(wfc, neg_orb_dir, step);
                energy_grad_->invalidate_hone_cache();
            }

            for (int i = 0; i < n_occ_params; ++i)
                params_new[i] = params[i] + step * occ_dir[i];

            // Compute occupations at trial point.
            double lambda_trial = 0.0;
            if (use_sigma_shift)
            {
                lambda_trial = sigma_shift_param_->compute_lambda(params_new);
                sigma_shift_param_->params_to_occ(params_new, lambda_trial, occ_flat_new);
            }
            else
            {
                occ_param_->params_to_occ(params_new, occ_flat_new);
            }

            if (grad_params_out && grad_wfc_out && dd_out)
            {
                std::vector<double> g_occ_trial;
                psi::Psi<TK> g_wfc_trial;
                double E_val = energy_grad_->compute(occ_flat_new, wfc, g_occ_trial, g_wfc_trial);
                if (use_sigma_shift)
                {
                    sigma_shift_param_->transform_gradient_batch(
                        g_occ_trial, params_new, lambda_trial, *grad_params_out);
                }
                else if (config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
                {
                    E_val += occ_constraint_->augmented_lagrangian_penalty(occ_flat_new);
                    std::vector<double> pen_g;
                    occ_constraint_->augmented_lagrangian_gradient(occ_flat_new, pen_g);
                    for (size_t i = 0; i < g_occ_trial.size(); ++i)
                        g_occ_trial[i] += pen_g[i];
                    occ_param_->transform_gradient_batch(g_occ_trial, params_new, *grad_params_out);
                }
                else
                {
                    occ_param_->transform_gradient_batch(g_occ_trial, params_new, *grad_params_out);
                }
                energy_grad_->project_orbital_gradient(wfc, g_wfc_trial);
                *grad_wfc_out = g_wfc_trial;

                double occ_dd_t = 0.0;
                for (int i = 0; i < n_occ_params; ++i)
                    occ_dd_t += occ_dir[i] * (*grad_params_out)[i];
                const double orb_dd_t = energy_grad_->s_inner_product(*grad_wfc_out, orb_dir);
                *dd_out = occ_dd_t + orb_dd_t;
                return E_val;
            }
            return total_energy(occ_flat_new, wfc);
        };

        if (joint_use_sw)
        {
            auto phi_joint = [&](double step) -> std::pair<double, double> {
                std::vector<double> g_params_t;
                psi::Psi<TK> g_wfc_t;
                double dd_t = 0.0;
                const double E_t = eval_at_step(step, &g_params_t, &g_wfc_t, &dd_t);
                return {E_t, dd_t};
            };

            const LineSearchResult sw_result
                = strong_wolfe_line_search(phi_joint, E, dd_total, alpha_init, c1, config_.line_search_c2,
                    config_.line_search_max_iter, config_.line_search_max_zoom);
            alpha = sw_result.step;
            E_new = sw_result.f_new;
            ls_success = sw_result.success;
            joint_n_feval = sw_result.n_feval;

            if (step_orbitals)
            {
                for (int ik = 0; ik < nk; ++ik)
                    for (int ib = 0; ib < nb_local; ++ib)
                        for (int mu = 0; mu < nbs_local; ++mu)
                            wfc(ik, ib, mu) = wfc_save(ik, ib, mu);
                energy_grad_->retract_orbitals(wfc, neg_orb_dir, alpha);
                energy_grad_->invalidate_hone_cache();
            }
            for (int i = 0; i < n_occ_params; ++i)
                params_new[i] = params[i] + alpha * occ_dir[i];
            if (use_sigma_shift)
            {
                const double lam_sw = sigma_shift_param_->compute_lambda(params_new);
                sigma_shift_param_->params_to_occ(params_new, lam_sw, occ_flat_new);
            }
            else
            {
                occ_param_->params_to_occ(params_new, occ_flat_new);
            }
        }
        else
        {
            for (int ls = 0; ls < config_.line_search_max_iter; ++ls)
            {
                joint_ls_trials = ls + 1;
                E_new = eval_at_step(alpha, nullptr, nullptr, nullptr);
                if (E_new <= E + c1 * alpha * dd_total)
                {
                    ls_success = true;
                    break;
                }
                alpha *= rho;
            }
        }

        if (ls_success)
        {
            if (joint_use_sw)
            {
                GlobalV::ofs_running << "  RDMFT joint-iter " << (iter + 1)
                    << "  line search (joint StrongWolfe): alpha_init=" << std::scientific << joint_ls_alpha0
                    << " step=" << alpha << " n_feval=" << joint_n_feval << "  E0=" << E
                    << " E1=" << E_new << "  dd=" << dd_total << "  c1=" << std::defaultfloat << c1
                    << " c2=" << config_.line_search_c2
                    << "  step_orb=" << (step_orbitals ? 1 : 0) << std::endl;
            }
            else
            {
                GlobalV::ofs_running << "  RDMFT joint-iter " << (iter + 1)
                    << "  line search (joint Armijo): alpha_init=" << std::scientific << joint_ls_alpha0
                    << " step=" << alpha << " n_trial=" << joint_ls_trials << "  E0=" << E
                    << " E1=" << E_new << "  dd=" << dd_total << "  c1=" << std::defaultfloat << c1
                    << " rho=" << rho << "  step_orb=" << (step_orbitals ? 1 : 0) << std::endl;
            }
        }

        if (!ls_success)
        {
            if (joint_ls_failed_prev)
            {
                GlobalV::ofs_running << "  RDMFT joint-iter " << (iter + 1)
                    << "  line search failed again after optimiser reset (SD recovery); stopping outer loop."
                    << std::endl;
                last_result_.converged = false;
                break;
            }
            joint_ls_failed_prev = true;
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

            GlobalV::ofs_running << "  RDMFT joint-iter " << (iter + 1)
                << "  line search (joint " << (joint_use_sw ? "StrongWolfe" : "Armijo") << ") failed: alpha_init="
                << std::scientific << joint_ls_alpha0 << " last_alpha=" << alpha;
            if (joint_use_sw)
            {
                GlobalV::ofs_running << " n_feval=" << joint_n_feval;
            }
            else
            {
                GlobalV::ofs_running << " n_trial=" << joint_ls_trials;
            }
            GlobalV::ofs_running << "  E_last=" << E_new << std::defaultfloat
                << ", restarting optimiser state (next outer: steepest-descent step)" << std::endl;
            E_prev = E;
            continue;
        }

        joint_ls_failed_prev = false;

        // Commit occupation step (orbitals are already retracted by the
        // final, successful line-search trial above).
        for (int i = 0; i < n_occ_params; ++i)
        {
            packed_step[i] = alpha * occ_dir[i];
            params[i] += packed_step[i];
        }
        if (use_sigma_shift)
        {
            // Recompute λ for the new params, then update occupations.
            sigma_lambda = sigma_shift_param_->compute_lambda(params);
            sigma_shift_param_->params_to_occ(params, sigma_lambda, occ_flat);
        }
        else
        {
            occ_param_->params_to_occ(params, occ_flat);
        }

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

            std::vector<double> new_grad_params;
            if (use_sigma_shift)
            {
                // sigma_lambda was already updated when committing the step above.
                sigma_shift_param_->transform_gradient_batch(
                    new_grad_occ, params, sigma_lambda, new_grad_params);
            }
            else
            {
                if (config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
                {
                    std::vector<double> pen_grad;
                    occ_constraint_->augmented_lagrangian_gradient(occ_flat, pen_grad);
                    for (size_t i = 0; i < new_grad_occ.size(); ++i)
                        new_grad_occ[i] += pen_grad[i];
                }
                occ_param_->transform_gradient_batch(new_grad_occ, params, new_grad_params);
            }
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

        // Augmented-Lagrangian multiplier refresh (not needed for sigma-shift).
        if (!use_sigma_shift && config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
        {
            occ_constraint_->update_multiplier(occ_flat);
            if (std::abs(occ_constraint_->constraint_violation(occ_flat)) > 1e-6)
            {
                occ_constraint_->increase_penalty(
                    config_.aug_lag_mu_factor, config_.aug_lag_mu_max);
            }
        }

        const double dE = std::abs(E_new - E_prev);
        const double sum_abs_dn_joint = sum_abs_diff(occ_flat, occ_before_joint_step);
        double gnorm2_occ = 0.0;
        for (auto g : grad_params) gnorm2_occ += g * g;
        const double gnorm_occ = std::sqrt(gnorm2_occ);
        const double gnorm_orb = std::sqrt(std::max(0.0, orb_gnorm2));
        const double gnorm_total = std::sqrt(gnorm2_occ + std::max(0.0, orb_gnorm2));
        const bool occ_conv_joint = (sum_abs_dn_joint < config_.rdmft_occ_tol);
        const bool orb_conv_joint = (gnorm_orb < config_.orb_grad_tol);
        const double E_one = energy_grad_->E_one_body();
        const double E_hartree = energy_grad_->E_hartree();
        const double E_xc = energy_grad_->E_xc();
        const double E_entropy = energy_grad_->E_entropy();
        const double E_ewald = energy_grad_->E_ewald();
        const double E_total = energy_grad_->E_total();
        const double E_penalty = (config_.constraint_method == ConstraintMethod::AugmentedLagrangian)
            ? occ_constraint_->augmented_lagrangian_penalty(occ_flat)
            : 0.0;

        GlobalV::ofs_running << std::fixed << std::setprecision(10)
            << "  RDMFT joint-iter " << iter + 1
            << "  E = " << E_new
            << "  dE = " << std::scientific << dE
            << "  alpha = " << alpha
            << "  sum|dn|_step (L1) = " << sum_abs_dn_joint
            << "  sum(w*n) = " << occ_constraint_->weighted_occupation_sum(occ_flat)
            << "  N_e = " << n_electrons_
            << "  |grad_occ| = " << gnorm_occ
            << "  |grad_orb| = " << gnorm_orb
            << "  |grad_total| = " << gnorm_total
            << "  occ_conv = " << (occ_conv_joint ? "Y" : "N")
            << "  orb_conv = " << (orb_conv_joint ? "Y" : "N")
            << std::endl;

        if (config_.occ_entropy_gamma > 0.0 && config_.xc_type == XCFunctionalType::HF)
        {
            print_rdmft_energy_table_running({"E_one_elec", "E_Hartree", "E_xc", "E_entropy", "E_Ewald", "E_total",
                                              "E_penalty", "E_aug"},
                {E_one, E_hartree, E_xc, E_entropy, E_ewald, E_total, E_penalty, E_total + E_penalty});
        }
        else
        {
            print_rdmft_energy_table_running(
                {"E_one_elec", "E_Hartree", "E_xc", "E_Ewald", "E_total", "E_penalty", "E_aug"},
                {E_one, E_hartree, E_xc, E_ewald, E_total, E_penalty, E_total + E_penalty});
        }

        // stdout: one block per joint outer iteration (inner line search stays off stdout).
        std::cout << std::fixed << std::setprecision(10)
                  << "  RDMFT joint iter " << (iter + 1) << "  E_total=" << E_total << "  (Ry)";
        if (iter > 0)
        {
            std::cout << std::scientific << "  dE=" << (E_new - E_prev) << std::fixed;
        }
        std::cout << std::defaultfloat << std::endl;
        {
            std::ostringstream hdr_end;
            hdr_end << std::fixed << std::setprecision(10)
                    << "  RDMFT joint iter " << (iter + 1)
                    << "  occupations n(ik,ib) after step"
                    << "  sum(w*n)=" << occ_constraint_->weighted_occupation_sum(occ_flat)
                    << "  N_e=" << n_electrons_;
            if (use_sigma_shift)
            {
                hdr_end << "  sigma_lambda=" << sigma_lambda;
            }
            print_occ_table_to_stream(GlobalV::ofs_running, occ_flat, nk_, nbands_, hdr_end.str(), 2);
        }

        const double abs_c_joint = std::abs(occ_constraint_->constraint_violation(occ_flat));
        joint_last_E = E_new;
        joint_last_dE = dE;
        joint_last_abs_c = abs_c_joint;
        joint_gn_occ = gnorm_occ;
        joint_gn_orb = gnorm_orb;
        joint_gn_tot = gnorm_total;
        joint_outer_done = iter + 1;

        // Outer stop: energy change below tol (if tol > 0) or both occ/orb stationarity flags.
        const bool inner_both_joint = occ_conv_joint && orb_conv_joint;
        const bool energy_ok_joint
            = (config_.energy_tol > 0.0) && (dE < config_.energy_tol);
        const bool outer_converged_joint = (iter > 0) && (energy_ok_joint || inner_both_joint);
        if (outer_converged_joint)
        {
            last_result_.converged = true;
            last_result_.iterations = iter + 1;
            last_result_.final_energy = E_new;
            last_result_.grad_norm = gnorm_total;
            E = E_new;
            if (energy_ok_joint && inner_both_joint)
            {
                GlobalV::ofs_running << "  RDMFT joint: outer loop stopped (|dE| < rdmft_energy_tol and OCC&ORB "
                                        "flags)"
                                     << std::endl;
            }
            else if (energy_ok_joint)
            {
                GlobalV::ofs_running << "  RDMFT joint: outer loop stopped (|dE| < rdmft_energy_tol)"
                                     << std::endl;
            }
            else
            {
                GlobalV::ofs_running << "  RDMFT joint: outer loop stopped (OCC&ORB flags)" << std::endl;
            }
            break;
        }
        E_prev = E_new;
        E = E_new;
    }

    if (use_sigma_shift)
    {
        sigma_lambda = sigma_shift_param_->compute_lambda(params);
        sigma_shift_param_->params_to_occ(params, sigma_lambda, occ_flat);
    }
    else
    {
        occ_param_->params_to_occ(params, occ_flat);
    }
    last_result_.final_energy = E;
    if (!last_result_.converged)
    {
        last_result_.iterations = joint_outer_done;
        last_result_.grad_norm = joint_gn_tot;
    }

    print_rdmft_outer_energy_stdout(last_result_.converged, E);
    // Always print final occupations in tabular form.
    print_occ_table_running(occ_flat, nk_, nbands_);
    const auto t_joint_end = std::chrono::steady_clock::now();
    const double total_time_sec = std::chrono::duration<double>(t_joint_end - t_joint_start).count();
    print_rdmft_optimization_summary_joint(joint_outer_done, total_time_sec);
    if (!last_result_.converged)
    {
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

    if (!logged_initial_occ_for_first_occ_gradient_)
    {
        GlobalV::ofs_running << "\n  RDMFT: initial occupations for the first occupation-gradient (∂E/∂n) step"
                             << "\n     (fractional n after init mode / projection; energy uses this same n)"
                             << "\n     sum(w*n)=" << occ_constraint_->weighted_occupation_sum(occ_flat)
                             << "  N_e=" << n_electrons_ << std::endl;
        print_occ_table_to_stream(GlobalV::ofs_running, occ_flat, nk_, nbands_, "     n(ik,ib):", 3);
        logged_initial_occ_for_first_occ_gradient_ = true;
    }

    switch (config_.constraint_method)
    {
        case ConstraintMethod::AugmentedLagrangian:
        {
            // Convert to unconstrained parameters
            std::vector<double> params;
            occ_param_->occ_to_params(occ_flat, params);
            // Commit n ≡ chart(params) so the first ∂E/∂n is evaluated at the same fractional occupations
            // as the unconstrained parameters (avoids navigating from a mismatched (p, n) pair).
            occ_param_->params_to_occ(params, occ_flat);

            BarzilaiBorweinStep bb_step;
            bb_step.set_mode(config_.alm_bb_mode);
            bb_step.set_bounds(config_.alm_bb_alpha_min, config_.alm_bb_alpha_max);
            bb_step.reset();

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

                if (inner == 0 && result.grad_norm < config_.occ_grad_tol)
                {
                    result.converged = true;
                    result.iterations = 1;
                    result.final_energy = L;
                    print_inner_loop_stdout("RDMFT occ inner", 1, result.final_energy, occ_flat, nk_, nbands_);
                    GlobalV::ofs_running << "      occ inner ALM: converged at first inner (||dL/dp|| < occ_grad_tol); "
                                             "skipping line search."
                                         << std::endl;
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

                auto phi_at_step = [&](double step) -> std::pair<double, double> {
                    std::vector<double> p_trial(params);
                    for (size_t i = 0; i < p_trial.size(); ++i)
                        p_trial[i] += step * dir[i];
                    std::vector<double> occ_trial;
                    occ_param_->params_to_occ(p_trial, occ_trial);
                    std::vector<double> g_occ_trial;
                    psi::Psi<TK> g_wfc_trial;
                    double Et = energy_grad_->compute(occ_trial,
                                    const_cast<psi::Psi<TK>&>(wfc), g_occ_trial, g_wfc_trial);
                    Et += occ_constraint_->augmented_lagrangian_penalty(occ_trial);
                    std::vector<double> pen_g;
                    occ_constraint_->augmented_lagrangian_gradient(occ_trial, pen_g);
                    for (size_t i = 0; i < g_occ_trial.size(); ++i)
                        g_occ_trial[i] += pen_g[i];
                    std::vector<double> g_params_trial;
                    occ_param_->transform_gradient_batch(g_occ_trial, p_trial, g_params_trial);
                    double dd_trial = 0.0;
                    for (size_t i = 0; i < dir.size(); ++i)
                        dd_trial += dir[i] * g_params_trial[i];
                    return {Et, dd_trial};
                };

                const bool use_strong_wolfe = line_search_uses_strong_wolfe(config_.occ_optimizer);
                LineSearchResult ls;
                if (use_strong_wolfe)
                {
                    const double sw_alpha0 = config_.alm_bb_enabled
                        ? bb_step.suggest(params, grad_params, config_.line_search_alpha_init)
                        : config_.line_search_alpha_init;
                    ls = strong_wolfe_line_search(phi_at_step, L, dd, sw_alpha0, config_.line_search_c1,
                        config_.line_search_c2, config_.line_search_max_iter, config_.line_search_max_zoom);
                }
                else
                {
                    ls = armijo_line_search(f_at_step, L, dd,
                        config_.alm_bb_enabled
                            ? bb_step.suggest(params, grad_params, config_.line_search_alpha_init)
                            : config_.line_search_alpha_init,
                        config_.line_search_c1,
                        config_.line_search_rho,
                        config_.line_search_max_iter,
                        config_.line_search_polynomial);
                }

                {
                    const char* const ls_name = use_strong_wolfe ? "StrongWolfe" : "Armijo";
                    GlobalV::ofs_running << "      occ line search (ALM " << ls_name << "): alpha_init=" << std::scientific
                        << ls.alpha_init << " step=" << ls.step << " n_feval=" << ls.n_feval << "  L0=" << L
                        << " L1=" << ls.f_new << " dd=" << dd;
                    if (!use_strong_wolfe)
                    {
                        GlobalV::ofs_running << "  c1=" << std::defaultfloat << config_.line_search_c1
                            << " rho=" << config_.line_search_rho;
                    }
                    else
                    {
                        GlobalV::ofs_running << "  c1=" << std::defaultfloat << config_.line_search_c1
                            << " c2=" << config_.line_search_c2;
                    }
                    GlobalV::ofs_running << (ls.success ? "  ok" : "  fail") << std::endl;
                }

                std::vector<double> step_vec(params.size());
                for (size_t i = 0; i < params.size(); ++i)
                {
                    step_vec[i] = ls.step * dir[i];
                    params[i] += step_vec[i];
                }

                if (!ls.success)
                {
                    GlobalV::ofs_running << "      ALM occ inner: "
                                         << (use_strong_wolfe ? "Strong Wolfe" : "Armijo")
                                         << " line search failed at inner=" << (inner + 1)
                                         << "; accepting best-effort step (step=" << std::scientific
                                         << ls.step << ")" << std::defaultfloat << std::endl;
                }

                occ_param_->params_to_occ(params, occ_flat);
                std::vector<double> new_grad_occ;
                energy_grad_->compute(occ_flat, const_cast<psi::Psi<TK>&>(wfc),
                                       new_grad_occ, grad_wfc_dummy);
                std::vector<double> new_grad_params;
                occ_param_->transform_gradient_batch(new_grad_occ, params, new_grad_params);
                opt.update(new_grad_params, step_vec);
                if (config_.alm_bb_enabled)
                {
                    bb_step.record_state(params, new_grad_params);
                }

                result.iterations = inner + 1;
                result.final_energy = ls.f_new;
                print_inner_loop_stdout("RDMFT occ inner", inner + 1, result.final_energy,
                                        occ_flat, nk_, nbands_);

                const double sum_abs_dn = sum_abs_diff(occ_flat, occ_at_step_start);
                GlobalV::ofs_running << "      sum|dn|_step (L1 move)=" << std::scientific
                    << sum_abs_dn << "  sum(w*n)=" << occ_constraint_->weighted_occupation_sum(occ_flat)
                    << "  N_e=" << n_electrons_ << std::endl;
                if (occ_inner_should_stop(sum_abs_dn,
                        config_.rdmft_occ_tol))
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
        case ConstraintMethod::ActiveSet:
        {
            // Spectral Projected Gradient (SPG, Birgin-Martinez-Raydan,
            // SIAM J. Optim. 10 (2000) 1196) with non-monotone Armijo line
            // search (Grippo-Lampariello-Lucidi, SIAM J. Numer. Anal. 23
            // (1986) 707).
            //
            // Per inner iteration k:
            //   1. Compute g_k = ∂E/∂n at current n_k.
            //   2. Spectral step α_k = (s,s)/(s,y) clamped to [α_min, α_max]
            //      with s = n_k - n_{k-1}, y = g_k - g_{k-1} (BB1; first
            //      iteration uses α_0 = 1 / max(1, ||g_k||_∞)).
            //   3. Spectral projected direction
            //          d_k = P_C(n_k − α_k g_k) − n_k.
            //      `d_k` is provably a descent direction with
            //          g_k · d_k ≤ −‖d_k‖² / α_k ≤ 0
            //      (Birgin-Martinez-Raydan, Lemma 2.1) for any closed convex
            //      set C.  This removes the need for any descent-direction
            //      safeguard or "trust-region cap" on α: the spectral step is
            //      already bounded by the BB clamp.
            //   4. Non-monotone Armijo: find the smallest j ≥ 0 such that
            //          E(n_k + λ_j d_k) ≤ f_max + c_1 λ_j (g_k · d_k)
            //      with λ_j = ρ^j and f_max = max_{0 ≤ i ≤ min(k, M-1)} f(n_{k-i}).
            //      Since C is convex and `n_k + λ d_k = (1−λ) n_k + λ P_C(...)`
            //      is a convex combination, every trial point is feasible
            //      *without* re-projection or rejection.
            //
            // Convergence: stationarity is measured by the τ-free Bertsekas
            // residual r = n − P_C(n − g) (textbook normalisation; the SPG
            // paper uses the same residual to check first-order optimality).
            // Stop when ‖r‖_∞ < `rdmft_occ_grad_tol` *or* the L1 occupation
            // move per step drops below `rdmft_occ_tol`.
            //
            // The legacy `ProjectedGradient` and `ActiveSet` constraint
            // methods both route to this SPG implementation.  The only
            // difference between the two used to be that AS additionally
            // zeroed the gradient at pinned occupations; that bookkeeping
            // is now handled implicitly by the projection P_C inside the
            // spectral step.  `rdmft_occ_optimizer` (sd / cg / lbfgs / adam)
            // is *not* consulted on this path: SPG already has a cleanly
            // proven globally convergent step and any preconditioner would
            // need to be tied to the projection in a non-trivial way.
            const char* const method_name
                = (config_.constraint_method == ConstraintMethod::ProjectedGradient) ? "SPG" : "SPG-AS";
            GlobalV::ofs_running << "      " << method_name
                                 << ": rdmft_occ_maxiter=" << config_.occ_maxiter << std::endl;

            // SPG safeguarding constants (Birgin-Martinez-Raydan, Algorithm 2.1).
            //   alpha_min, alpha_max  -- BB clamp; standard values.
            //   nm_memory M           -- non-monotone history length; M = 10
            //                            is the value used in the original SPG paper.
            //   rho                   -- backtracking ratio (re-uses the user's
            //                            line_search_rho so the fallback cadence
            //                            matches existing logs).
            const double alpha_min = 1e-10;
            const double alpha_max = 1e10;
            const int nm_memory = 10;
            const double rho = (config_.line_search_rho > 0.0 && config_.line_search_rho < 1.0)
                                   ? config_.line_search_rho
                                   : 0.5;
            const double c1 = config_.line_search_c1;

            std::deque<double> f_history;

            std::vector<double> occ_prev;
            std::vector<double> grad_prev;
            bool have_prev = false;

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

                double r_l2 = 0.0;
                double r_inf = 0.0;
                spg_stationarity_norms(occ_flat, grad_occ, *occ_constraint_, r_l2, r_inf);

                double grad_l2 = 0.0;
                for (auto g : grad_occ) grad_l2 += g * g;
                grad_l2 = std::sqrt(grad_l2);

                {
                    std::ostringstream os;
                    os << "      occ inner " << method_name << " " << (inner + 1) << "  E=" << std::fixed
                       << std::setprecision(10) << E << "  dE=" << std::scientific << dE
                       << "  |c|=" << c_abs
                       << "  ||r||_inf=||n-P(n-∇E)||_inf=" << r_inf
                       << "  ||r||_2=" << r_l2
                       << "  ||grad_n E||_2=" << grad_l2;
                    log_occ_inner_summary_and_nik(os.str(), occ_flat, occ_prev_for_dn, nk_, nbands_);
                }

                if (r_inf < config_.occ_grad_tol)
                {
                    result.converged = true;
                    result.iterations = inner + 1;
                    result.final_energy = E;
                    result.grad_norm = r_inf;
                    print_inner_loop_stdout("RDMFT occ inner", inner + 1, result.final_energy,
                                            occ_flat, nk_, nbands_);
                    GlobalV::ofs_running << "      occ inner " << method_name
                                         << ": converged (||r||_inf=" << std::scientific << r_inf
                                         << " < rdmft_occ_grad_tol=" << config_.occ_grad_tol << ")"
                                         << std::defaultfloat << std::endl;
                    break;
                }

                // Spectral (BB1) step.  First inner: α₀ = 1 / max(1, ||g||_∞)
                // (BMR's standard initial step; bounds the first trial direction
                // ‖α₀ g‖_∞ ≤ 1, i.e. at most a unit move in occupation space).
                double alpha_bb;
                if (have_prev)
                {
                    alpha_bb = spectral_bb_step(occ_flat, grad_occ, occ_prev, grad_prev,
                                                /*have_prev=*/true,
                                                /*fallback=*/1.0, alpha_min, alpha_max);
                }
                else
                {
                    double g_inf = 0.0;
                    for (double g : grad_occ) g_inf = std::max(g_inf, std::abs(g));
                    const double init = 1.0 / std::max(1.0, g_inf);
                    alpha_bb = std::min(alpha_max, std::max(alpha_min, init));
                }

                // Spectral projected direction d = P(n - α_BB g) - n.
                std::vector<double> trial_proj(occ_flat.size());
                for (size_t i = 0; i < occ_flat.size(); ++i)
                {
                    trial_proj[i] = occ_flat[i] - alpha_bb * grad_occ[i];
                }
                occ_constraint_->project(trial_proj);

                std::vector<double> dir(occ_flat.size());
                double dd = 0.0;
                for (size_t i = 0; i < occ_flat.size(); ++i)
                {
                    dir[i] = trial_proj[i] - occ_flat[i];
                    dd += grad_occ[i] * dir[i];
                }

                // Non-monotone reference (Grippo-Lampariello-Lucidi).
                if (static_cast<int>(f_history.size()) >= nm_memory)
                {
                    f_history.pop_front();
                }
                f_history.push_back(E);
                double f_max = E;
                for (double f : f_history)
                {
                    if (f > f_max) f_max = f;
                }

                // Non-monotone Armijo backtracking along the convex segment
                // n + λ d.  Every trial point is feasible (convex combination
                // of two feasible points), so no per-trial re-projection is needed.
                const std::vector<double> occ_before_step(occ_flat);
                double lambda = 1.0;
                int n_trial = 0;
                bool ls_success = false;
                std::vector<double> occ_trial(occ_flat.size());
                double E_last_trial = E;
                double lambda_acc = 0.0;
                for (int ls = 0; ls < config_.line_search_max_iter; ++ls)
                {
                    n_trial = ls + 1;
                    for (size_t i = 0; i < occ_flat.size(); ++i)
                    {
                        occ_trial[i] = occ_flat[i] + lambda * dir[i];
                    }
                    const double E_trial = energy_grad_->compute_energy(
                        occ_trial, const_cast<psi::Psi<TK>&>(wfc));
                    E_last_trial = E_trial;
                    if (E_trial <= f_max + c1 * lambda * dd)
                    {
                        ls_success = true;
                        lambda_acc = lambda;
                        break;
                    }
                    lambda *= rho;
                }

                {
                    std::ostringstream pg_ls;
                    pg_ls << "      occ line search (" << method_name << " non-monotone Armijo)"
                          << ": alpha_BB=" << std::scientific << alpha_bb
                          << " lambda=" << (ls_success ? lambda_acc : 0.0)
                          << " n_trial=" << n_trial
                          << " E0=" << E
                          << " E_trial=" << E_last_trial
                          << " f_max=" << f_max
                          << " dd=g·d=" << dd
                          << "  c1=" << std::defaultfloat << c1
                          << "  rho=" << rho
                          << (ls_success ? "  ok" : "  fail");
                    log_occ_inner_line(pg_ls.str());
                }

                // Save (n_k, g_k) for the next BB step *before* committing the move.
                occ_prev = occ_flat;
                grad_prev = grad_occ;
                have_prev = true;

                if (ls_success)
                {
                    occ_flat = occ_trial;
                }
                else
                {
                    // SPG line search cannot fail in exact arithmetic when d is
                    // a descent direction (which it provably is for the spectral
                    // projected step).  In floating point this can happen if dd
                    // is degenerate (≈ 0); accept zero step so the outer loop
                    // can move on.
                    GlobalV::ofs_running
                        << "      occ inner " << method_name
                        << ": non-monotone Armijo failed (dd=" << std::scientific << dd
                        << "); accepting zero step" << std::defaultfloat << std::endl;
                    occ_flat = occ_before_step;
                }

                // Diagnostics + convergence test on the new iterate.
                std::vector<double> new_grad_occ;
                const double E_post = energy_grad_->compute(occ_flat,
                                                            const_cast<psi::Psi<TK>&>(wfc),
                                                            new_grad_occ, grad_wfc_dummy);
                double r_l2_post = 0.0;
                double r_inf_post = 0.0;
                spg_stationarity_norms(occ_flat, new_grad_occ, *occ_constraint_,
                                       r_l2_post, r_inf_post);
                result.grad_norm = r_inf_post;
                result.iterations = inner + 1;
                result.final_energy = E_post;

                print_inner_loop_stdout("RDMFT occ inner", inner + 1, result.final_energy,
                                        occ_flat, nk_, nbands_);

                const double sum_abs_dn = sum_abs_diff(occ_flat, occ_before_step);
                GlobalV::ofs_running << "      sum|dn|_step (L1 move)=" << std::scientific
                    << sum_abs_dn
                    << "  sum(w*n)=" << occ_constraint_->weighted_occupation_sum(occ_flat)
                    << "  N_e=" << n_electrons_ << std::endl;
                GlobalV::ofs_running << "      " << method_name
                                     << " ||r||_inf (post-step)=" << std::scientific << r_inf_post
                                     << "  ||r||_2=" << r_l2_post
                                     << "  rdmft_occ_grad_tol=" << config_.occ_grad_tol
                                     << std::defaultfloat << std::endl;

                if (r_inf_post < config_.occ_grad_tol)
                {
                    result.converged = true;
                    break;
                }
                if (occ_inner_should_stop(sum_abs_dn, config_.rdmft_occ_tol))
                {
                    result.converged = true;
                    GlobalV::ofs_running
                        << "      occ inner " << method_name << ": converged on sum|dn|="
                        << std::scientific << sum_abs_dn
                        << " < rdmft_occ_tol=" << config_.rdmft_occ_tol
                        << std::defaultfloat << std::endl;
                    break;
                }
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
    // Plain Riemannian gradient method on the Stiefel manifold
    // (Absil-Mahony-Sepulchre, "Optimization Algorithms on Matrix Manifolds",
    // Princeton 2008, §4.2), specialised for each of the four optimiser
    // families exposed via `rdmft_orb_optimizer`.
    //
    // Notation. At inner iteration k:
    //   X_k       current orbital iterate on St(N_b, N) (X_k^H X_k = I).
    //   G_k       Euclidean gradient dE/dX evaluated at X_k.
    //   G_R       Riemannian gradient at X_k:
    //                G_R = G_k - X_k sym(X_k^H G_k)
    //             (`EnergyGradient::project_orbital_gradient`, with the
    //             generalised-Stiefel S already absorbed by the X-space
    //             change of variable `wfc_C_to_X`).
    //   R_X(eta)  retraction selected by `rdmft_orb_retraction` (Polar by
    //             default; QR or Cayley serial-only). `retract_orbitals(W, S, a)`
    //             implements W <- R_W(-a * S), so to follow direction D we pass
    //             step S = -D.
    //   <A, B>_F  Re Tr(A^H B) summed over k-points (S-inner_product, since
    //             after wfc_C_to_X we work with X^H X = I).
    //
    // Per inner iteration:
    //   1. Compute E and the Euclidean gradient G_k. Project to G_R.
    //      Convergence test (gradient and optional energy criterion) on
    //      pre-step values, identical to the ALM occupation inner loop.
    //   2. Build the search direction D in the tangent space at X_k:
    //        SD    : D = -G_R.
    //        CG    : Polak-Ribiere+ with vector transport by tangent-space
    //                projection of the previous (G_R, D) onto T_{X_k} St:
    //                  beta = max(0, <G_R, G_R - T(G_R^{prev})> / <T(G_R^{prev}), T(G_R^{prev})>)
    //                  D    = -G_R + beta * T(D_{prev})
    //                Re-project D to be safe against numerical drift.
    //        LBFGS : two-loop recursion in flat Euclidean coordinates
    //                produces a quasi-Newton direction; project onto
    //                T_{X_k} St. (Riemannian L-BFGS by projection.)
    //        Adam  : Euclidean Adam direction (m, v moments updated from
    //                the Euclidean gradient G_k), projected onto T_{X_k} St.
    //      Descent safeguard: if <G_R, D>_F >= 0, reset D = -G_R.
    //   3. Line search along the retracted curve X(alpha) = R_{X_k}(alpha * D):
    //        sd / cg / adam : monotone Armijo (rdmft_optimizer.h).
    //        lbfgs          : Strong Wolfe (Nocedal & Wright Algorithm 3.5).
    //      Initial step alpha_init: 1.0 for lbfgs / adam (the optimiser
    //      sets the natural scale); user-tunable line_search_alpha_init for
    //      sd / cg.
    //   4. Commit X_{k+1}. For lbfgs, update the (s, y) history with
    //      s = alpha * D (in flat coords) and y = G_R(X_{k+1}) - G_R(X_k)
    //      (transported by projection at X_{k+1}). Snapshot (G_R, D) for CG.
    //
    // No Barzilai-Borwein step, no non-monotone history, no trust radius,
    // no suspicious-descent guard. Convergence criterion: ||G_R||_F below
    // `rdmft_orb_grad_tol` (and, when `rdmft_orb_energy_tol > 0`, also
    // |E - E_prev| below `rdmft_orb_energy_tol`).

    OptResult result;

    const int nk = wfc.get_nk();
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();
    const int total_size = nk * nb_local * nbs_local;

    const OptimizerType opt_type = config_.orb_optimizer;
    const bool use_sd    = (opt_type == OptimizerType::SteepestDescent);
    const bool use_cg    = (opt_type == OptimizerType::ConjugateGradient);
    const bool use_lbfgs = (opt_type == OptimizerType::LBFGS);
    const bool use_adam  = (opt_type == OptimizerType::Adam);
    const char* orb_opt_name = use_lbfgs ? "lbfgs" : (use_adam ? "adam" : (use_cg ? "cg" : "sd"));
    const bool use_strong_wolfe = line_search_uses_strong_wolfe(opt_type);

    energy_grad_->invalidate_hone_cache();

    EuclideanOptimizer eucl_opt(opt_type, config_);
    int flat_size = 0;
    if (use_lbfgs || use_adam)
    {
        flat_size = total_size;
        if (!std::is_same<TK, double>::value) flat_size *= 2;
        eucl_opt.init(flat_size);
    }

    auto frob_inner = [&](const psi::Psi<TK>& A, const psi::Psi<TK>& B) {
        double s = 0.0;
        const int nkA = A.get_nk();
        const int nbA = A.get_nbands();
        const int nbsA = A.get_nbasis();
        for (int ik = 0; ik < nkA; ++ik)
        {
            const int nelem = nbA * nbsA;
            if (nelem == 0) continue;
            const TK* a = &A(ik, 0, 0);
            const TK* b = &B(ik, 0, 0);
            for (int i = 0; i < nelem; ++i) s += std::real(std::conj(a[i]) * b[i]);
        }
        Parallel_Reduce::reduce_all(s);
        return s;
    };

    auto frob2 = [&](const psi::Psi<TK>& A) {
        return frob_inner(A, A);
    };

    psi::Psi<TK> wfc_save(wfc);
    psi::Psi<TK> prev_grad(wfc);
    psi::Psi<TK> prev_dir(wfc);
    bool have_prev = false;

    double prev_orb_E = 0.0;

    for (int inner = 0; inner < config_.orb_maxiter; ++inner)
    {
        // ---- 1. Energy and Riemannian gradient at the current iterate ----
        std::vector<double> grad_occ;
        psi::Psi<TK> grad_wfc;
        const double E = energy_grad_->compute(
            const_cast<std::vector<double>&>(occ_flat), wfc, grad_occ, grad_wfc);

        energy_grad_->project_orbital_gradient(wfc, grad_wfc);
        const double gnorm2 = frob2(grad_wfc);
        result.grad_norm = std::sqrt(std::max(0.0, gnorm2));

        const bool orb_e_enabled = (config_.orb_energy_tol > 0.0);
        const double dE_abs = (inner > 0 && orb_e_enabled)
                                  ? std::abs(E - prev_orb_E)
                                  : std::numeric_limits<double>::infinity();

        {
            std::ostringstream os;
            os << "      orb inner " << inner + 1
               << "  E=" << std::fixed << std::setprecision(10) << E
               << "  gnorm=" << std::scientific << result.grad_norm;
            if (inner > 0)
                os << "  dE=" << std::scientific << (E - prev_orb_E);
            GlobalV::ofs_running << os.str() << std::endl;
        }
        prev_orb_E = E;

        const bool grad_conv = (result.grad_norm < config_.orb_grad_tol);
        const bool energy_conv = orb_e_enabled && (dE_abs < config_.orb_energy_tol);
        const bool inner_energy_ok = !orb_e_enabled || energy_conv;
        if (grad_conv && inner_energy_ok)
        {
            result.converged = true;
            result.iterations = inner + 1;
            result.final_energy = E;
            print_inner_loop_stdout("RDMFT orb inner", inner + 1, result.final_energy,
                                    occ_flat, nk_, nbands_, false);
            break;
        }

        // ---- 2. Build search direction D in T_{X_k} St ----
        psi::Psi<TK> dir(wfc);
        for (int ik = 0; ik < nk; ++ik)
            for (int ib = 0; ib < nb_local; ++ib)
                for (int mu = 0; mu < nbs_local; ++mu)
                    dir(ik, ib, mu) = -grad_wfc(ik, ib, mu);

        if (use_lbfgs || use_adam)
        {
            // Riemannian L-BFGS / Adam by projection: take the Euclidean
            // direction produced by the optimiser and project it onto
            // T_{X_k} St. The (s, y) history for L-BFGS is built from
            // tangent-space increments at each iterate, so the two-loop
            // recursion already operates "morally" on the manifold; the
            // projection at the end keeps the result tangent.
            std::vector<double> grad_flat, dir_flat;
            psi_to_flat(grad_wfc, grad_flat);
            eucl_opt.compute_direction(grad_flat, dir_flat);
            flat_to_psi(dir_flat, dir);
            energy_grad_->project_orbital_gradient(wfc, dir);
        }
        else if (use_cg && have_prev)
        {
            // Polak-Ribiere+ CG with vector transport by tangent-space
            // projection of (G_R^{prev}, D_{prev}) onto T_{X_k}.
            psi::Psi<TK> g_prev_T(prev_grad);
            energy_grad_->project_orbital_gradient(wfc, g_prev_T);
            const double gn2_prev = frob2(g_prev_T);
            if (gn2_prev > 1e-30)
            {
                psi::Psi<TK> y(grad_wfc);
                for (int ik = 0; ik < nk; ++ik)
                    for (int ib = 0; ib < nb_local; ++ib)
                        for (int mu = 0; mu < nbs_local; ++mu)
                            y(ik, ib, mu) = grad_wfc(ik, ib, mu) - g_prev_T(ik, ib, mu);
                double beta_pr = frob_inner(grad_wfc, y) / gn2_prev;
                beta_pr = std::max(0.0, beta_pr);
                psi::Psi<TK> d_prev_T(prev_dir);
                energy_grad_->project_orbital_gradient(wfc, d_prev_T);
                for (int ik = 0; ik < nk; ++ik)
                    for (int ib = 0; ib < nb_local; ++ib)
                        for (int mu = 0; mu < nbs_local; ++mu)
                            dir(ik, ib, mu) = -grad_wfc(ik, ib, mu)
                                              + TK(beta_pr) * d_prev_T(ik, ib, mu);
                energy_grad_->project_orbital_gradient(wfc, dir);
            }
        }

        // Descent safeguard: with our convention `retract_orbitals(W, S, a)`
        // applies W <- R_W(-a S), the Armijo slope along D = -S is
        // -<G_R, S> at a = 0. Equivalently <G_R, D> must be < 0 for the
        // step to descend. Reset to SD on failure.
        double dd = energy_grad_->s_inner_product(grad_wfc, dir);
        if (!(dd < 0.0))
        {
            for (int ik = 0; ik < nk; ++ik)
                for (int ib = 0; ib < nb_local; ++ib)
                    for (int mu = 0; mu < nbs_local; ++mu)
                        dir(ik, ib, mu) = -grad_wfc(ik, ib, mu);
            dd = energy_grad_->s_inner_product(grad_wfc, dir);
        }
        (void)use_sd;

        // Snapshot for line-search rollback / repeat trials.
        for (int ik = 0; ik < nk; ++ik)
            for (int ib = 0; ib < nb_local; ++ib)
                for (int mu = 0; mu < nbs_local; ++mu)
                    wfc_save(ik, ib, mu) = wfc(ik, ib, mu);

        // retract_orbitals applies wfc <- R_{wfc}(-lambda * step). To follow +D
        // we pass step = -D.
        psi::Psi<TK> neg_dir(dir);
        if (total_size > 0)
        {
            TK* q = &neg_dir(0, 0, 0);
            const TK* p = &dir(0, 0, 0);
            for (int i = 0; i < total_size; ++i) q[i] = -p[i];
        }

        // ---- 3. Line search along the retracted curve ----
        // L-BFGS / Adam set their own scale; their natural initial step is
        // alpha_init = 1.0. SD / CG use the user's line_search_alpha_init.
        const double alpha_init_cfg = (config_.line_search_alpha_init > 0.0)
                                          ? config_.line_search_alpha_init : 1.0;
        const double alpha_init = (use_lbfgs || use_adam) ? 1.0 : alpha_init_cfg;
        const double c1  = config_.line_search_c1;
        const double rho = (config_.line_search_rho > 0.0 && config_.line_search_rho < 1.0)
                               ? config_.line_search_rho : 0.5;

        auto eval_at_step = [&](double a) -> double {
            for (int ik = 0; ik < nk; ++ik)
                for (int ib = 0; ib < nb_local; ++ib)
                    for (int mu = 0; mu < nbs_local; ++mu)
                        wfc(ik, ib, mu) = wfc_save(ik, ib, mu);
            energy_grad_->retract_orbitals(wfc, neg_dir, a);
            energy_grad_->invalidate_hone_cache();
            return energy_grad_->compute_energy(occ_flat, wfc);
        };

        LineSearchResult ls;
        if (use_strong_wolfe)
        {
            // Strong Wolfe: phi(a) returns (E(a), <G_R(X(a)), D>_F). The
            // gradient slope is computed by retracting, evaluating
            // grad_wfc at the trial point, projecting onto T_{X(a)}, and
            // taking <G_R, D>_F (D fixed in the ambient space). This is the
            // standard Wolfe along the retracted curve when D is treated
            // as a constant ambient direction.
            psi::Psi<TK> grad_trial(wfc);
            std::vector<double> grad_occ_trial;
            auto phi = [&](double a) -> std::pair<double, double> {
                for (int ik = 0; ik < nk; ++ik)
                    for (int ib = 0; ib < nb_local; ++ib)
                        for (int mu = 0; mu < nbs_local; ++mu)
                            wfc(ik, ib, mu) = wfc_save(ik, ib, mu);
                energy_grad_->retract_orbitals(wfc, neg_dir, a);
                energy_grad_->invalidate_hone_cache();
                grad_occ_trial.clear();
                const double E_a = energy_grad_->compute(
                    const_cast<std::vector<double>&>(occ_flat), wfc, grad_occ_trial, grad_trial);
                energy_grad_->project_orbital_gradient(wfc, grad_trial);
                const double slope = energy_grad_->s_inner_product(grad_trial, dir);
                return {E_a, slope};
            };
            ls = strong_wolfe_line_search(
                phi, E, dd, alpha_init, c1, config_.line_search_c2,
                config_.line_search_max_iter, config_.line_search_max_zoom);
        }
        else
        {
            ls = armijo_line_search(
                eval_at_step, E, dd, alpha_init, c1, rho,
                config_.line_search_max_iter, config_.line_search_polynomial);
        }

        const double lambda = ls.success ? ls.step : 0.0;
        const double E_new = ls.success ? ls.f_new : E;

        {
            std::ostringstream orb_ls;
            orb_ls << "      orb line search ("
                   << (use_strong_wolfe ? "Strong Wolfe" : "Armijo")
                   << "): optim=" << orb_opt_name
                   << " alpha_init=" << std::scientific << alpha_init
                   << " lambda=" << lambda
                   << " n_feval=" << ls.n_feval
                   << " E0=" << E << " E_trial=" << E_new
                   << "  dd=<G_R,dir>=" << dd
                   << "  c1=" << std::defaultfloat << c1 << "  rho=" << rho
                   << (ls.success ? "  ok" : "  fail");
            GlobalV::ofs_running << orb_ls.str() << std::endl;
        }

        // Save (G_R^k, D_k) BEFORE the move so prev_grad / prev_dir refer
        // to the *previous* iterate when CG transports them next round.
        for (int ik = 0; ik < nk; ++ik)
            for (int ib = 0; ib < nb_local; ++ib)
                for (int mu = 0; mu < nbs_local; ++mu)
                {
                    prev_grad(ik, ib, mu) = grad_wfc(ik, ib, mu);
                    prev_dir(ik, ib, mu) = dir(ik, ib, mu);
                }
        have_prev = true;

        if (!ls.success)
        {
            // Roll back; line search failed (rare, usually a stationary
            // iterate). Let the outer loop notice via gnorm.
            for (int ik = 0; ik < nk; ++ik)
                for (int ib = 0; ib < nb_local; ++ib)
                    for (int mu = 0; mu < nbs_local; ++mu)
                        wfc(ik, ib, mu) = wfc_save(ik, ib, mu);
            energy_grad_->invalidate_hone_cache();
            result.iterations = inner + 1;
            result.final_energy = E;
            print_inner_loop_stdout("RDMFT orb inner", inner + 1, result.final_energy,
                                    occ_flat, nk_, nbands_, false);
            break;
        }

        // The line search above may have been called repeatedly with
        // intermediate trial steps; restore the iterate corresponding to
        // the accepted step length lambda.
        for (int ik = 0; ik < nk; ++ik)
            for (int ib = 0; ib < nb_local; ++ib)
                for (int mu = 0; mu < nbs_local; ++mu)
                    wfc(ik, ib, mu) = wfc_save(ik, ib, mu);
        energy_grad_->retract_orbitals(wfc, neg_dir, lambda);
        energy_grad_->invalidate_hone_cache();

        // ---- 4. Post-step gradient + L-BFGS history update ----
        std::vector<double> grad_occ_new;
        psi::Psi<TK> grad_wfc_new;
        energy_grad_->compute(
            const_cast<std::vector<double>&>(occ_flat), wfc, grad_occ_new, grad_wfc_new);
        energy_grad_->project_orbital_gradient(wfc, grad_wfc_new);
        const double post_gnorm2 = frob2(grad_wfc_new);
        result.grad_norm = std::sqrt(std::max(0.0, post_gnorm2));

        if (use_lbfgs)
        {
            // Riemannian L-BFGS by projection: store s = lambda * D (in
            // flat coords) and y = G_R(X_{k+1}) - G_R(X_k) so the next
            // two-loop recursion sees a curvature pair that approximates
            // the Stiefel Hessian along the connecting curve.
            std::vector<double> new_grad_flat, dir_flat_step, step_flat;
            psi_to_flat(grad_wfc_new, new_grad_flat);
            psi_to_flat(dir, dir_flat_step);
            step_flat.resize(dir_flat_step.size());
            for (size_t i = 0; i < step_flat.size(); ++i)
                step_flat[i] = lambda * dir_flat_step[i];
            eucl_opt.update(new_grad_flat, step_flat);
        }

        const double dE_step_abs = std::abs(E_new - E);
        const bool orb_e_post = (config_.orb_energy_tol > 0.0);
        const bool post_energy_ok = !orb_e_post || (dE_step_abs < config_.orb_energy_tol);
        const bool post_grad_conv = (result.grad_norm < config_.orb_grad_tol);
        if (post_grad_conv && post_energy_ok)
        {
            result.converged = true;
            result.iterations = inner + 1;
            result.final_energy = E_new;
            print_inner_loop_stdout("RDMFT orb inner", inner + 1, result.final_energy,
                                    occ_flat, nk_, nbands_, false);
            energy_grad_->invalidate_hone_cache();
            break;
        }

        energy_grad_->invalidate_hone_cache();
        result.iterations = inner + 1;
        result.final_energy = E_new;
        print_inner_loop_stdout("RDMFT orb inner", inner + 1, result.final_energy,
                                occ_flat, nk_, nbands_, false);
    }

    GlobalV::ofs_running << "      orb sub-problem: " << result.iterations << " iters"
        << "  optim=" << orb_opt_name
        << "  final_gnorm=" << std::scientific << result.grad_norm
        << "  final_E=" << std::fixed << std::setprecision(10) << result.final_energy
        << (result.converged ? "  CONVERGED" : "  not converged")
        << std::endl;

    return result;
}

// Run an occupation- and orbital-gradient FD check at a single given (occ, wfc)
// point. Logs to GlobalV::ofs_running and returns true iff every per-band
// occupation FD and the orbital directional-derivative FD pass `tolerance`.
template <typename TK, typename TR>
static bool grad_check_at_point(
    rdmft::EnergyGradient<TK, TR>& energy_grad,
    const std::vector<double>& occ_flat,
    const psi::Psi<TK>& wfc,
    const std::string& label,
    double epsilon,
    double tolerance)
{
    GlobalV::ofs_running << "\n----- Probe: " << label << " -----" << std::endl;
    bool all_pass = true;

    // Print a short occupation summary at this probe so the log shows what
    // configuration we are testing.
    {
        const std::size_t n = occ_flat.size();
        const std::size_t n_print = std::min<std::size_t>(n, 12);
        std::ostringstream os;
        os << "  occ summary [first " << n_print << " of " << n << "]: ";
        os << std::fixed << std::setprecision(4);
        for (std::size_t i = 0; i < n_print; ++i) os << occ_flat[i] << ' ';
        double s = 0.0, mn = (n > 0 ? occ_flat[0] : 0.0), mx = mn;
        for (std::size_t i = 0; i < n; ++i)
        {
            s += occ_flat[i];
            mn = std::min(mn, occ_flat[i]);
            mx = std::max(mx, occ_flat[i]);
        }
        os << " sum=" << s << "  min=" << mn << "  max=" << mx;
        GlobalV::ofs_running << os.str() << std::endl;
    }

    // Analytic gradients at the probe point.
    std::vector<double> grad_occ;
    psi::Psi<TK> grad_wfc;
    double E0 = energy_grad.compute(
        const_cast<std::vector<double>&>(occ_flat),
        const_cast<psi::Psi<TK>&>(wfc),
        grad_occ, grad_wfc);

    // ---- Occupation gradient FD ----
    GlobalV::ofs_running << "  -- Occupation gradient --" << std::endl;
    const std::size_t n_occ_check = std::min<std::size_t>(occ_flat.size(), 10);
    for (std::size_t idx = 0; idx < n_occ_check; ++idx)
    {
        std::vector<double> occ_plus(occ_flat);
        std::vector<double> occ_minus(occ_flat);
        const double n_orig = occ_flat[idx];
        if (n_orig + epsilon > 1.0 || n_orig - epsilon < 0.0) continue;
        occ_plus[idx] += epsilon;
        occ_minus[idx] -= epsilon;

        const double E_plus
            = energy_grad.compute_energy(occ_plus, const_cast<psi::Psi<TK>&>(wfc));
        const double Ep_one = energy_grad.E_one_body();
        const double Ep_h = energy_grad.E_hartree();
        const double Ep_x = energy_grad.E_xc();
        const double E_minus
            = energy_grad.compute_energy(occ_minus, const_cast<psi::Psi<TK>&>(wfc));
        const double Em_one = energy_grad.E_one_body();
        const double Em_h = energy_grad.E_hartree();
        const double Em_x = energy_grad.E_xc();

        const double fd_grad = (E_plus - E_minus) / (2.0 * epsilon);
        const double fd_one = (Ep_one - Em_one) / (2.0 * epsilon);
        const double fd_h = (Ep_h - Em_h) / (2.0 * epsilon);
        const double fd_x = (Ep_x - Em_x) / (2.0 * epsilon);

        const double analytic = grad_occ[idx];
        const double rel_err = (std::abs(analytic) > 1e-10)
                                   ? std::abs(fd_grad - analytic) / std::abs(analytic)
                                   : std::abs(fd_grad - analytic);
        const bool pass = rel_err < tolerance;
        if (!pass) all_pass = false;

        GlobalV::ofs_running << std::fixed << std::setprecision(8)
            << "    occ[" << idx << "]"
            << "  n=" << n_orig
            << "  analytic=" << analytic
            << "  fd=" << fd_grad
            << "  (fd_one=" << fd_one
            << "  fd_h=" << fd_h
            << "  fd_x=" << fd_x << ")"
            << "  rel_err=" << rel_err
            << (pass ? "  PASS" : "  FAIL") << std::endl;
    }

    // ---- Orbital directional-derivative FD ----
    // Project the Euclidean gradient onto the Stiefel tangent (matching
    // optimize_orbitals), then forward-differentiate along R(- t * G_R).
    GlobalV::ofs_running << "  -- Orbital directional derivative --" << std::endl;
    energy_grad.project_orbital_gradient(const_cast<psi::Psi<TK>&>(wfc), grad_wfc);
    psi::Psi<TK> G_dir(grad_wfc);
    const double gnorm2 = energy_grad.s_inner_product(grad_wfc, grad_wfc);
    const double analytic_dd = -gnorm2;

    auto f_at = [&](double step) -> double {
        psi::Psi<TK> wfc_trial(wfc);
        energy_grad.retract_orbitals(wfc_trial, G_dir, step);
        energy_grad.invalidate_hone_cache();
        return energy_grad.compute_energy(occ_flat, wfc_trial);
    };
    const double fd = (f_at(epsilon) - E0) / epsilon;

    const double orb_rel_err = (std::abs(analytic_dd) > 1e-10)
                                   ? std::abs(fd - analytic_dd) / std::abs(analytic_dd)
                                   : std::abs(fd - analytic_dd);
    const bool orb_pass = orb_rel_err < tolerance;
    if (!orb_pass) all_pass = false;

    GlobalV::ofs_running << std::fixed << std::setprecision(8)
        << "    analytic_dd=" << analytic_dd
        << "  fd_fwd=" << fd
        << "  ||G_R||^2=" << gnorm2
        << "  rel_err=" << orb_rel_err
        << (orb_pass ? "  PASS" : "  FAIL") << std::endl;

    energy_grad.invalidate_hone_cache();
    return all_pass;
}

template <typename TK, typename TR>
bool RDMFTSolver<TK, TR>::check_gradient_consistency(
    const std::vector<double>& occ_flat,
    const psi::Psi<TK>& wfc,
    double epsilon,
    double tolerance)
{
    GlobalV::ofs_running << "\n===== RDMFT Gradient Consistency Check =====" << std::endl;
    bool all_pass = true;

    // 1. Probe at the input (KS-derived) occupations
    all_pass &= grad_check_at_point<TK, TR>(*energy_grad_, occ_flat, wfc,
                                             "input occupations", epsilon, tolerance);

    // Build several "strange" occupation configurations on top of the input
    // (so the orbital iterate / k-weights / nbands match), each renormalised
    // to the input weighted electron count. They probe regions of n that
    // RDMFT optimisation can visit (uniform spread, near-zero, near-one,
    // randomly perturbed, inverted), where the n^alpha coupling and its
    // regularisation matter.
    const int nk = wfc.get_nk();
    const int nb = nbands_;
    const std::size_t total = static_cast<std::size_t>(nk) * static_cast<std::size_t>(nb);
    if (occ_flat.size() == total && nk > 0 && nb > 0)
    {
        const auto& kw = (kv_ != nullptr) ? kv_->wk : std::vector<double>(nk, 1.0);
        // Reference electron count Σ_k w_k Σ_b n_kb (matches the constraint).
        auto weighted_sum = [&](const std::vector<double>& occ) {
            double s = 0.0;
            for (int ik = 0; ik < nk; ++ik)
            {
                const double w = (ik < static_cast<int>(kw.size())) ? kw[ik] : 1.0;
                for (int ib = 0; ib < nb; ++ib)
                    s += w * occ[ik * nb + ib];
            }
            return s;
        };
        const double Ne_target = weighted_sum(occ_flat);

        // Affinely rescale `occ` (clipped to [eps, 1-eps]) to a target weighted sum.
        // Uses one bisection on a uniform shift x s.t. Σ w_k Σ_b clip(occ_kb + x) = N_e.
        auto rescale_to_Ne = [&](std::vector<double>& occ, double Ne) {
            const double clip_eps = 1e-6;
            auto clip = [&](double n) {
                return std::min(1.0 - clip_eps, std::max(clip_eps, n));
            };
            auto S = [&](double x) {
                double s = 0.0;
                for (int ik = 0; ik < nk; ++ik)
                {
                    const double w = (ik < static_cast<int>(kw.size())) ? kw[ik] : 1.0;
                    for (int ib = 0; ib < nb; ++ib)
                        s += w * clip(occ[ik * nb + ib] + x);
                }
                return s;
            };
            double lo = -1.0, hi = 1.0;
            for (int it = 0; it < 100; ++it)
            {
                const double mid = 0.5 * (lo + hi);
                ((S(mid) > Ne) ? hi : lo) = mid;
                if (hi - lo < 1e-12) break;
            }
            const double xstar = 0.5 * (lo + hi);
            for (auto& n : occ) n = clip(n + xstar);
        };

        // Probe (a): uniform smear (n = N_e / (Σ w · nbands) clipped per-band).
        {
            std::vector<double> occ(total, 0.0);
            const double w_sum_nb = [&]() {
                double s = 0.0;
                for (int ik = 0; ik < nk; ++ik)
                    s += ((ik < static_cast<int>(kw.size())) ? kw[ik] : 1.0) * nb;
                return s;
            }();
            const double n_uniform = (w_sum_nb > 0.0) ? (Ne_target / w_sum_nb) : 0.5;
            std::fill(occ.begin(), occ.end(), n_uniform);
            rescale_to_Ne(occ, Ne_target);
            all_pass &= grad_check_at_point<TK, TR>(*energy_grad_, occ, wfc,
                                                     "uniform smear", epsilon, tolerance);
        }

        // Probe (b): near-zero with one band per k near 1 (sparsest configuration).
        {
            std::vector<double> occ(total, 1e-3);
            for (int ik = 0; ik < nk; ++ik) occ[ik * nb + 0] = 1.0 - 1e-3;
            rescale_to_Ne(occ, Ne_target);
            all_pass &= grad_check_at_point<TK, TR>(*energy_grad_, occ, wfc,
                                                     "near-empty (one nearly-full per k)",
                                                     epsilon, tolerance);
        }

        // Probe (c): inverted (1 - n) - swap occupied / virtual, renormalised.
        {
            std::vector<double> occ(occ_flat);
            for (auto& n : occ) n = 1.0 - n;
            rescale_to_Ne(occ, Ne_target);
            all_pass &= grad_check_at_point<TK, TR>(*energy_grad_, occ, wfc,
                                                     "inverted (1 - n)", epsilon, tolerance);
        }

        // Probe (d): pseudo-random perturbation (deterministic, reproducible).
        {
            std::vector<double> occ(occ_flat);
            // LCG-based pseudo-random in [0, 1); deterministic & MPI-independent.
            std::uint64_t s = 0x9E3779B97F4A7C15ULL;
            for (std::size_t i = 0; i < occ.size(); ++i)
            {
                s = s * 6364136223846793005ULL + 1442695040888963407ULL;
                const double r = static_cast<double>((s >> 32) & 0xFFFFFFFFULL)
                                 / static_cast<double>(0x100000000ULL);
                occ[i] += 0.2 * (r - 0.5); // +/- 0.1 perturbation
            }
            rescale_to_Ne(occ, Ne_target);
            all_pass &= grad_check_at_point<TK, TR>(*energy_grad_, occ, wfc,
                                                     "random perturbation +/- 0.1",
                                                     epsilon, tolerance);
        }

        // Probe (e): half-filled fractional (every band at 0.5, rescaled).
        {
            std::vector<double> occ(total, 0.5);
            rescale_to_Ne(occ, Ne_target);
            all_pass &= grad_check_at_point<TK, TR>(*energy_grad_, occ, wfc,
                                                     "half-filled (n = 0.5)",
                                                     epsilon, tolerance);
        }
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
