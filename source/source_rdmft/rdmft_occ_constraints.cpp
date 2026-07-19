#include "rdmft_occ_constraints.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace rdmft
{

namespace
{
inline double clip01(double v)
{
    return std::max(0.0, std::min(1.0, v));
}
const double INV_SQRT_PI = 0.5641895835477562869480794515607725858440506293289988;
} // namespace

double OccConstraints::weighted_sum(const std::vector<double>& occ) const
{
    double s = 0.0;
    for (int ik = 0; ik < nks; ++ik)
    {
        const double w = wk[ik];
        const int base = ik * nbnd;
        for (int ib = 0; ib < nbnd; ++ib)
        {
            s += w * occ[base + ib];
        }
    }
    return s;
}

double OccConstraints::weighted_sum_ispin(const std::vector<double>& occ, int ispin) const
{
    double s = 0.0;
    for (int ik = 0; ik < nks; ++ik)
    {
        if (isk[ik] != ispin)
        {
            continue;
        }
        const double w = wk[ik];
        const int base = ik * nbnd;
        for (int ib = 0; ib < nbnd; ++ib)
        {
            s += w * occ[base + ib];
        }
    }
    return s;
}

double OccConstraints::magnetization_electrons(const std::vector<double>& occ) const
{
    return weighted_sum_ispin(occ, 1) - weighted_sum_ispin(occ, 2);
}

double OccConstraints::constraint_violation(const std::vector<double>& occ) const
{
    if (fix_magnetization)
    {
        return std::fabs(weighted_sum_ispin(occ, 1) - n_target_up)
               + std::fabs(weighted_sum_ispin(occ, 2) - n_target_down);
    }
    return std::fabs(weighted_sum(occ) - n_target);
}

double OccConstraints::eval_weighted_shift_sum(const std::vector<double>& x, double shift,
                                               int ispin, bool use_filter) const
{
    double ssum = 0.0;
    for (int ik = 0; ik < nks; ++ik)
    {
        if (use_filter && isk[ik] != ispin)
        {
            continue;
        }
        const double w = wk[ik];
        const double sub = shift * w;
        const int base = ik * nbnd;
        for (int ib = 0; ib < nbnd; ++ib)
        {
            ssum += w * clip01(x[base + ib] - sub);
        }
    }
    return ssum;
}

double OccConstraints::eval_weighted_uniform_sum(const std::vector<double>& x, double mu,
                                                 int ispin, bool use_filter) const
{
    double ssum = 0.0;
    for (int ik = 0; ik < nks; ++ik)
    {
        if (use_filter && isk[ik] != ispin)
        {
            continue;
        }
        const double w = wk[ik];
        const int base = ik * nbnd;
        for (int ib = 0; ib < nbnd; ++ib)
        {
            ssum += w * clip01(x[base + ib] - mu);
        }
    }
    return ssum;
}

void OccConstraints::proximal_project_core(std::vector<double>& occ, double target_ne, int ispin,
                                           bool use_filter) const
{
    const double tol_sum = 1.0e-12;
    const double expand_max = 1.0e12;
    std::vector<double> x = occ;

    bool already_box = true;
    for (int ik = 0; ik < nks && already_box; ++ik)
    {
        if (use_filter && isk[ik] != ispin)
        {
            continue;
        }
        const int base = ik * nbnd;
        for (int ib = 0; ib < nbnd; ++ib)
        {
            const double v = x[base + ib];
            if (v < -1.0e-15 || v > 1.0 + 1.0e-15)
            {
                already_box = false;
                break;
            }
        }
    }

    const double s0 = eval_weighted_shift_sum(x, 0.0, ispin, use_filter);
    if (already_box && std::fabs(s0 - target_ne) < tol_sum)
    {
        for (int ik = 0; ik < nks; ++ik)
        {
            if (use_filter && isk[ik] != ispin)
            {
                continue;
            }
            const int base = ik * nbnd;
            for (int ib = 0; ib < nbnd; ++ib)
            {
                occ[base + ib] = clip01(x[base + ib]);
            }
        }
        return;
    }

    double lam_lo, lam_hi, sum_lo, sum_hi;
    if (s0 > target_ne)
    {
        lam_lo = 0.0;
        lam_hi = 1.0;
        sum_lo = s0;
        sum_hi = eval_weighted_shift_sum(x, lam_hi, ispin, use_filter);
        while (sum_hi > target_ne && lam_hi < expand_max)
        {
            lam_hi *= 2.0;
            sum_hi = eval_weighted_shift_sum(x, lam_hi, ispin, use_filter);
        }
    }
    else
    {
        lam_hi = 0.0;
        lam_lo = -1.0;
        sum_hi = s0;
        sum_lo = eval_weighted_shift_sum(x, lam_lo, ispin, use_filter);
        while (sum_lo < target_ne && std::fabs(lam_lo) < expand_max)
        {
            lam_lo *= 2.0;
            sum_lo = eval_weighted_shift_sum(x, lam_lo, ispin, use_filter);
        }
    }

    double shift;
    if (!(sum_lo >= target_ne && sum_hi <= target_ne))
    {
        shift = (std::fabs(target_ne - sum_lo) < std::fabs(target_ne - sum_hi)) ? lam_lo : lam_hi;
    }
    else
    {
        for (int it = 0; it < 100; ++it)
        {
            const double lam_mid = 0.5 * (lam_lo + lam_hi);
            const double sum_mid = eval_weighted_shift_sum(x, lam_mid, ispin, use_filter);
            if (std::fabs(sum_mid - target_ne) < tol_sum)
            {
                lam_lo = lam_mid;
                lam_hi = lam_mid;
                break;
            }
            if (sum_mid > target_ne)
            {
                lam_lo = lam_mid;
            }
            else
            {
                lam_hi = lam_mid;
            }
        }
        shift = 0.5 * (lam_lo + lam_hi);
    }

    for (int ik = 0; ik < nks; ++ik)
    {
        if (use_filter && isk[ik] != ispin)
        {
            continue;
        }
        const double sub = shift * wk[ik];
        const int base = ik * nbnd;
        for (int ib = 0; ib < nbnd; ++ib)
        {
            occ[base + ib] = clip01(x[base + ib] - sub);
        }
    }
}

void OccConstraints::proximal_project(std::vector<double>& occ) const
{
    if (fix_magnetization)
    {
        proximal_project_core(occ, n_target_up, 1, true);
        proximal_project_core(occ, n_target_down, 2, true);
    }
    else
    {
        proximal_project_core(occ, n_target, 0, false);
    }
}

void OccConstraints::proximal_project_uniform_core(std::vector<double>& occ, double target_ne,
                                                   int ispin, bool use_filter) const
{
    const double tol_sum = 1.0e-12;
    const double expand_max = 1.0e12;
    std::vector<double> x = occ;

    const double s0 = eval_weighted_uniform_sum(x, 0.0, ispin, use_filter);
    if (std::fabs(s0 - target_ne) < tol_sum)
    {
        for (int ik = 0; ik < nks; ++ik)
        {
            if (use_filter && isk[ik] != ispin)
            {
                continue;
            }
            const int base = ik * nbnd;
            for (int ib = 0; ib < nbnd; ++ib)
            {
                occ[base + ib] = clip01(x[base + ib]);
            }
        }
        return;
    }

    double mu_lo, mu_hi, sum_lo, sum_hi;
    if (s0 > target_ne)
    {
        mu_lo = 0.0;
        mu_hi = 1.0;
        sum_lo = s0;
        sum_hi = eval_weighted_uniform_sum(x, mu_hi, ispin, use_filter);
        while (sum_hi > target_ne && mu_hi < expand_max)
        {
            mu_hi *= 2.0;
            sum_hi = eval_weighted_uniform_sum(x, mu_hi, ispin, use_filter);
        }
    }
    else
    {
        mu_hi = 0.0;
        mu_lo = -1.0;
        sum_hi = s0;
        sum_lo = eval_weighted_uniform_sum(x, mu_lo, ispin, use_filter);
        while (sum_lo < target_ne && std::fabs(mu_lo) < expand_max)
        {
            mu_lo *= 2.0;
            sum_lo = eval_weighted_uniform_sum(x, mu_lo, ispin, use_filter);
        }
    }

    double mu;
    if (!(sum_lo >= target_ne && sum_hi <= target_ne))
    {
        mu = (std::fabs(target_ne - sum_lo) < std::fabs(target_ne - sum_hi)) ? mu_lo : mu_hi;
    }
    else
    {
        for (int it = 0; it < 100; ++it)
        {
            const double mu_mid = 0.5 * (mu_lo + mu_hi);
            const double sum_mid = eval_weighted_uniform_sum(x, mu_mid, ispin, use_filter);
            if (std::fabs(sum_mid - target_ne) < tol_sum)
            {
                mu_lo = mu_mid;
                mu_hi = mu_mid;
                break;
            }
            if (sum_mid > target_ne)
            {
                mu_lo = mu_mid;
            }
            else
            {
                mu_hi = mu_mid;
            }
        }
        mu = 0.5 * (mu_lo + mu_hi);
    }

    for (int ik = 0; ik < nks; ++ik)
    {
        if (use_filter && isk[ik] != ispin)
        {
            continue;
        }
        const int base = ik * nbnd;
        for (int ib = 0; ib < nbnd; ++ib)
        {
            occ[base + ib] = clip01(x[base + ib] - mu);
        }
    }
}

void OccConstraints::proximal_project_uniform(std::vector<double>& occ) const
{
    if (fix_magnetization)
    {
        proximal_project_uniform_core(occ, n_target_up, 1, true);
        proximal_project_uniform_core(occ, n_target_down, 2, true);
    }
    else
    {
        proximal_project_uniform_core(occ, n_target, 0, false);
    }
}

double OccConstraints::pg_map_grad_inf(const std::vector<double>& occ,
                                       const std::vector<double>& grad, double alpha) const
{
    std::vector<double> trial(occ.size());
    for (size_t i = 0; i < occ.size(); ++i)
    {
        trial[i] = occ[i] - alpha * grad[i];
    }
    proximal_project(trial);
    double gp = 0.0;
    for (size_t i = 0; i < occ.size(); ++i)
    {
        gp = std::max(gp, std::fabs(occ[i] - trial[i]));
    }
    return gp;
}

double OccConstraints::pg_kkt_residual(const std::vector<double>& occ,
                                       const std::vector<double>& grad, double boundary_tol) const
{
    double V = -std::numeric_limits<double>::max();
    double W = std::numeric_limits<double>::max();
    for (int i = 0; i < size(); ++i)
    {
        const double e = grad[i];
        const bool at_upper = (occ[i] >= 1.0 - boundary_tol);
        const bool at_lower = (occ[i] <= boundary_tol);
        if (at_upper)
        {
            V = std::max(V, e);
        }
        if (at_lower)
        {
            W = std::min(W, e);
        }
        if (!at_upper && !at_lower)
        {
            V = std::max(V, e);
            W = std::min(W, e);
        }
    }
    if (V < -0.5 * std::numeric_limits<double>::max() || W > 0.5 * std::numeric_limits<double>::max())
    {
        return 0.0;
    }
    return std::max(0.0, V - W);
}

// ----------------------------------------------------------------------------
// EBI erf parameterisation
// ----------------------------------------------------------------------------

double OccConstraints::ebi_occ_from_arg(double arg)
{
    return 0.5 * (std::erf(arg) + 1.0);
}

double OccConstraints::ebi_occ_prime(double arg)
{
    return INV_SQRT_PI * std::exp(-arg * arg);
}

double OccConstraints::ebi_erfinv(double y)
{
    const double eps = 1.0e-12;
    const double tol = 1.0e-14;
    const double yclip = std::max(-1.0 + eps, std::min(1.0 - eps, y));
    if (std::fabs(yclip) < 1.0e-30)
    {
        return 0.0;
    }
    double x = std::copysign(std::sqrt(-std::log(0.5 * (1.0 - std::fabs(yclip)))), yclip);
    for (int it = 0; it < 20; ++it)
    {
        const double dy = std::erf(x) - yclip;
        if (std::fabs(dy) < tol)
        {
            break;
        }
        const double dfdx = 2.0 * INV_SQRT_PI * std::exp(-x * x);
        x -= dy / dfdx;
    }
    return x;
}

double OccConstraints::ebi_weighted_occ_sum(const std::vector<double>& x, double mu, int ispin,
                                            bool use_filter) const
{
    double ssum = 0.0;
    for (int ik = 0; ik < nks; ++ik)
    {
        if (use_filter && isk[ik] != ispin)
        {
            continue;
        }
        const double w = wk[ik];
        const int base = ik * nbnd;
        for (int ib = 0; ib < nbnd; ++ib)
        {
            ssum += w * ebi_occ_from_arg(x[base + ib] + mu);
        }
    }
    return ssum;
}

double OccConstraints::ebi_solve_mu(const std::vector<double>& x, double target, int ispin,
                                    bool use_filter) const
{
    const double tol = 1.0e-12;
    const double mu_bound = 1.0e6;
    const int max_bisect = 100;

    const double s0 = ebi_weighted_occ_sum(x, 0.0, ispin, use_filter);
    if (std::fabs(s0 - target) < tol)
    {
        return 0.0;
    }

    double mu_lo, mu_hi, s_lo, s_hi;
    if (s0 < target)
    {
        mu_lo = 0.0;
        mu_hi = 1.0;
        s_hi = ebi_weighted_occ_sum(x, mu_hi, ispin, use_filter);
        while (s_hi < target && mu_hi < mu_bound)
        {
            mu_hi *= 2.0;
            s_hi = ebi_weighted_occ_sum(x, mu_hi, ispin, use_filter);
        }
    }
    else
    {
        mu_lo = -1.0;
        mu_hi = 0.0;
        s_lo = ebi_weighted_occ_sum(x, mu_lo, ispin, use_filter);
        while (s_lo > target && mu_lo > -mu_bound)
        {
            mu_lo *= 2.0;
            s_lo = ebi_weighted_occ_sum(x, mu_lo, ispin, use_filter);
        }
    }

    double mu = 0.5 * (mu_lo + mu_hi);
    for (int it = 0; it < max_bisect; ++it)
    {
        mu = 0.5 * (mu_lo + mu_hi);
        const double s_mid = ebi_weighted_occ_sum(x, mu, ispin, use_filter);
        if (std::fabs(s_mid - target) < tol)
        {
            return mu;
        }
        if (s_mid < target)
        {
            mu_lo = mu;
        }
        else
        {
            mu_hi = mu;
        }
    }
    return mu;
}

void OccConstraints::ebi_params_to_occ(const std::vector<double>& x, double mu_up, double mu_dw,
                                       std::vector<double>& occ) const
{
    occ.resize(size());
    const bool split = fix_magnetization;
    for (int ik = 0; ik < nks; ++ik)
    {
        const double mu = split ? (isk[ik] == 1 ? mu_up : mu_dw) : mu_up;
        const int base = ik * nbnd;
        for (int ib = 0; ib < nbnd; ++ib)
        {
            occ[base + ib] = ebi_occ_from_arg(x[base + ib] + mu);
        }
    }
}

void OccConstraints::ebi_occ_to_params(const std::vector<double>& occ, std::vector<double>& x) const
{
    const double eps = 1.0e-12;
    x.resize(size());
    for (int i = 0; i < size(); ++i)
    {
        const double nclip = std::min(std::max(eps, occ[i]), 1.0 - eps);
        x[i] = ebi_erfinv(2.0 * nclip - 1.0);
    }
}

void OccConstraints::ebi_fold_mu(std::vector<double>& x) const
{
    if (fix_magnetization)
    {
        const double mu_up = ebi_solve_mu(x, n_target_up, 1, true);
        const double mu_dw = ebi_solve_mu(x, n_target_down, 2, true);
        for (int ik = 0; ik < nks; ++ik)
        {
            const double mu = (isk[ik] == 1) ? mu_up : mu_dw;
            const int base = ik * nbnd;
            for (int ib = 0; ib < nbnd; ++ib)
            {
                x[base + ib] += mu;
            }
        }
    }
    else
    {
        const double mu = ebi_solve_mu(x, n_target, 0, false);
        for (int i = 0; i < size(); ++i)
        {
            x[i] += mu;
        }
    }
}

void OccConstraints::ebi_transform_gradient(const std::vector<double>& grad_n,
                                            const std::vector<double>& x, double mu_up, double mu_dw,
                                            std::vector<double>& grad_x) const
{
    const double sum_s_floor = 1.0e-30;
    grad_x.assign(size(), 0.0);

    if (!fix_magnetization)
    {
        double sum_s = 0.0;
        double sum_hsp = 0.0;
        for (int ik = 0; ik < nks; ++ik)
        {
            const double w = wk[ik];
            const int base = ik * nbnd;
            for (int ib = 0; ib < nbnd; ++ib)
            {
                const double sp = ebi_occ_prime(x[base + ib] + mu_up);
                sum_s += w * sp;
                sum_hsp += grad_n[base + ib] * sp;
            }
        }
        if (std::fabs(sum_s) < sum_s_floor)
        {
            return;
        }
        const double ratio = sum_hsp / sum_s;
        for (int ik = 0; ik < nks; ++ik)
        {
            const double w = wk[ik];
            const int base = ik * nbnd;
            for (int ib = 0; ib < nbnd; ++ib)
            {
                const double sp = ebi_occ_prime(x[base + ib] + mu_up);
                grad_x[base + ib] = sp * (grad_n[base + ib] - w * ratio);
            }
        }
        return;
    }

    // Spin-resolved: independent ratios per channel.
    for (int spin = 1; spin <= 2; ++spin)
    {
        const double mu = (spin == 1) ? mu_up : mu_dw;
        double sum_s = 0.0;
        double sum_hsp = 0.0;
        for (int ik = 0; ik < nks; ++ik)
        {
            if (isk[ik] != spin)
            {
                continue;
            }
            const double w = wk[ik];
            const int base = ik * nbnd;
            for (int ib = 0; ib < nbnd; ++ib)
            {
                const double sp = ebi_occ_prime(x[base + ib] + mu);
                sum_s += w * sp;
                sum_hsp += grad_n[base + ib] * sp;
            }
        }
        if (std::fabs(sum_s) < sum_s_floor)
        {
            continue;
        }
        const double ratio = sum_hsp / sum_s;
        for (int ik = 0; ik < nks; ++ik)
        {
            if (isk[ik] != spin)
            {
                continue;
            }
            const double w = wk[ik];
            const int base = ik * nbnd;
            for (int ib = 0; ib < nbnd; ++ib)
            {
                const double sp = ebi_occ_prime(x[base + ib] + mu);
                grad_x[base + ib] = sp * (grad_n[base + ib] - w * ratio);
            }
        }
    }
}

void OccConstraints::ebi_sync_occ_from_x(const std::vector<double>& x, std::vector<double>& occ) const
{
    double mu_up, mu_dw = 0.0;
    if (fix_magnetization)
    {
        mu_up = ebi_solve_mu(x, n_target_up, 1, true);
        mu_dw = ebi_solve_mu(x, n_target_down, 2, true);
    }
    else
    {
        mu_up = ebi_solve_mu(x, n_target, 0, false);
    }
    ebi_params_to_occ(x, mu_up, mu_dw, occ);
}

} // namespace rdmft
