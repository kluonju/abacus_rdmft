#include "rdmft_line_search.h"

#include <algorithm>
#include <cmath>

namespace rdmft_core
{

LineSearchResult armijo_line_search(const LineSearchEval& eval, double f0, double phi0,
                                    double alpha0, double c1, double rho, int max_iter,
                                    double f_ref, double alpha_floor)
{
    LineSearchResult res;
    if (phi0 >= 0.0)
    {
        return res; // not a descent direction
    }
    double alpha = alpha0 > 0.0 ? alpha0 : 1.0;
    const double ref = (f_ref > f0) ? f_ref : f0;

    for (int it = 0; it < max_iter; ++it)
    {
        double f = 0.0;
        double g = 0.0;
        int ierr = 0;
        eval(alpha, f, g, ierr);
        ++res.n_eval;
        if (ierr == 0 && f <= ref + c1 * alpha * phi0)
        {
            res.success = true;
            res.step = alpha;
            res.f_new = f;
            res.g_new = g;
            return res;
        }
        alpha *= rho;
        if (alpha < alpha_floor)
        {
            break;
        }
    }
    return res;
}

namespace
{
// Cubic/quadratic-safeguarded zoom for the strong-Wolfe search.
LineSearchResult zoom(const LineSearchEval& eval, double f0, double phi0, double c1, double c2,
                      double f_ref, double a_lo, double f_lo, double phi_lo, double a_hi,
                      int max_zoom, int& n_eval)
{
    (void)f0;
    (void)phi_lo;
    LineSearchResult res;
    for (int it = 0; it < max_zoom; ++it)
    {
        const double a_j = 0.5 * (a_lo + a_hi);
        double f_j = 0.0;
        double g_j = 0.0;
        int ierr = 0;
        eval(a_j, f_j, g_j, ierr);
        ++n_eval;
        if (ierr != 0)
        {
            a_hi = a_j;
            continue;
        }
        if (f_j > f_ref + c1 * a_j * phi0 || f_j >= f_lo)
        {
            a_hi = a_j;
        }
        else
        {
            if (std::fabs(g_j) <= -c2 * phi0)
            {
                res.success = true;
                res.step = a_j;
                res.f_new = f_j;
                res.g_new = g_j;
                return res;
            }
            if (g_j * (a_hi - a_lo) >= 0.0)
            {
                a_hi = a_lo;
            }
            a_lo = a_j;
            f_lo = f_j;
            phi_lo = g_j;
        }
    }
    return res;
}
} // namespace

LineSearchResult strong_wolfe_line_search(const LineSearchEval& eval, double f0, double phi0,
                                          double alpha0, double c1, double c2, int max_iter,
                                          int max_zoom, double f_ref)
{
    LineSearchResult res;
    if (phi0 >= 0.0)
    {
        return res;
    }
    const double ref = (f_ref > f0) ? f_ref : f0;
    double a_prev = 0.0;
    double f_prev = f0;
    double phi_prev = phi0;
    double alpha = alpha0 > 0.0 ? alpha0 : 1.0;

    for (int it = 0; it < max_iter; ++it)
    {
        double f = 0.0;
        double g = 0.0;
        int ierr = 0;
        eval(alpha, f, g, ierr);
        ++res.n_eval;
        if (ierr != 0)
        {
            // Trial infeasible: shrink towards the previous point.
            alpha = 0.5 * (a_prev + alpha);
            if (alpha <= a_prev)
            {
                break;
            }
            continue;
        }
        if (f > ref + c1 * alpha * phi0 || (it > 0 && f >= f_prev))
        {
            int nz = 0;
            LineSearchResult z = zoom(eval, f0, phi0, c1, c2, ref, a_prev, f_prev, phi_prev,
                                      alpha, max_zoom, nz);
            res.n_eval += nz;
            if (z.success)
            {
                res.success = true;
                res.step = z.step;
                res.f_new = z.f_new;
                res.g_new = z.g_new;
            }
            return res;
        }
        if (std::fabs(g) <= -c2 * phi0)
        {
            res.success = true;
            res.step = alpha;
            res.f_new = f;
            res.g_new = g;
            return res;
        }
        if (g >= 0.0)
        {
            int nz = 0;
            LineSearchResult z = zoom(eval, f0, phi0, c1, c2, ref, alpha, f, g, a_prev, max_zoom, nz);
            res.n_eval += nz;
            if (z.success)
            {
                res.success = true;
                res.step = z.step;
                res.f_new = z.f_new;
                res.g_new = z.g_new;
            }
            return res;
        }
        a_prev = alpha;
        f_prev = f;
        phi_prev = g;
        alpha *= 2.0;
    }
    return res;
}

// ----------------------------------------------------------------------------
// Barzilai-Borwein
// ----------------------------------------------------------------------------

void BarzilaiBorwein::init(double alpha_min, double alpha_max)
{
    alpha_min_ = alpha_min;
    alpha_max_ = alpha_max;
    reset();
}

void BarzilaiBorwein::reset()
{
    have_prev_ = false;
    x_prev_.clear();
    g_prev_.clear();
}

void BarzilaiBorwein::record(const std::vector<double>& x, const std::vector<double>& g)
{
    x_prev_ = x;
    g_prev_ = g;
    have_prev_ = true;
}

double BarzilaiBorwein::spectral_step(const std::vector<double>& x, const std::vector<double>& g,
                                      double alpha_init) const
{
    const double fallback = std::min(alpha_max_, std::max(alpha_min_, alpha_init));
    if (!have_prev_ || x_prev_.size() != x.size() || g_prev_.size() != g.size())
    {
        return fallback;
    }
    double sy = 0.0;
    double ss = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        const double s = x[i] - x_prev_[i];
        const double y = g[i] - g_prev_[i];
        sy += s * y;
        ss += s * s;
    }
    if (sy <= 0.0)
    {
        return fallback;
    }
    return std::min(alpha_max_, std::max(alpha_min_, ss / sy));
}

} // namespace rdmft_core
