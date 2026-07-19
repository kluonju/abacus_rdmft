#include "rdmft_xc.h"

#include <algorithm>
#include <cmath>

namespace rdmft_core
{

namespace
{
const double HYBOPT_POWER_W = 0.938328;
const double HYBOPT_POWER_EXP = 0.541076;
const double BOWMOD_DEFAULT_BETA = 0.375;
const double BOWMOD_COEF_HF = 1.0;
const double BOWMOD_COEF_ALPHA = -0.14648044160527188;
const double BOWMOD_COEF_MIX = 0.666761947146244;

inline double clip01(double n)
{
    return std::max(0.0, std::min(1.0, n));
}
} // namespace

RdmftXC::RdmftXC(XcType type, double alpha, double reg_eps)
{
    set(type, alpha, reg_eps);
}

void RdmftXC::set(XcType type, double alpha, double reg_eps)
{
    type_ = type;
    alpha_ = alpha;
    reg_eps_ = reg_eps;
}

// Piecewise power coupling on [0,1] (ABACUS pow_reg / QE pow_g).
double RdmftXC::pow_g(double n, double a) const
{
    const double nc = clip01(n);
    if (a >= 1.0 - 1.0e-12)
    {
        return nc;
    }
    const double eps = reg_eps_;
    if (nc >= eps)
    {
        return std::pow(nc, a);
    }
    const double g_eps = std::pow(eps, a);
    const double dg_eps = a * std::pow(eps, a - 1.0);
    return g_eps + dg_eps * (nc - eps);
}

// Derivative of pow_g (ABACUS dpow_reg / QE pow_dg).
double RdmftXC::pow_dg(double n, double a) const
{
    const double nc = clip01(n);
    if (a >= 1.0 - 1.0e-12)
    {
        return 1.0;
    }
    const double eps = reg_eps_;
    const double nmin = std::max(nc, eps);
    return a * std::pow(nmin, a - 1.0);
}

double RdmftXC::g(double n) const
{
    if (type_ == XcType::HF)
    {
        return std::max(0.0, n);
    }
    return pow_g(n, alpha_);
}

double RdmftXC::dg(double n) const
{
    if (type_ == XcType::HF)
    {
        return 1.0;
    }
    return pow_dg(n, alpha_);
}

bool RdmftXC::is_separable() const
{
    return type_ == XcType::HF || type_ == XcType::Muller || type_ == XcType::Power;
}

int RdmftXC::n_channels() const
{
    switch (type_)
    {
    case XcType::HF:
    case XcType::Muller:
    case XcType::Power:
    case XcType::GU:
        return 1;
    case XcType::CHF:
    case XcType::CGA:
    case XcType::HybOpt:
        return 2;
    case XcType::GEO:
        return 3;
    case XcType::BOWMOD:
        return 4;
    default:
        return 1;
    }
}

XcChannel RdmftXC::channel(int it, double n) const
{
    XcChannel ch;
    const double nc = clip01(n);

    switch (type_)
    {
    case XcType::HF:
    case XcType::Muller:
    case XcType::Power:
    case XcType::GU:
        ch.coef = 1.0;
        ch.w = g(nc);
        ch.dw = dg(nc);
        break;
    case XcType::CHF:
        if (it == 1)
        {
            ch.coef = 0.5;
            ch.w = nc;
            ch.dw = 1.0;
        }
        else
        {
            ch.coef = 0.5;
            const double s = std::max(0.0, nc * (1.0 - nc));
            ch.w = std::sqrt(s);
            ch.dw = (s <= reg_eps_ * reg_eps_) ? 0.0 : 0.5 * (1.0 - 2.0 * nc) / ch.w;
        }
        break;
    case XcType::CGA:
        if (it == 1)
        {
            ch.coef = 0.25;
            ch.w = nc;
            ch.dw = 1.0;
        }
        else
        {
            ch.coef = 0.25;
            const double s = std::max(0.0, nc * (2.0 - nc));
            ch.w = std::sqrt(s);
            ch.dw = (s <= reg_eps_ * reg_eps_) ? 0.0 : (1.0 - nc) / ch.w;
        }
        break;
    case XcType::GEO:
    {
        double a = 1.0;
        if (it == 1)
        {
            ch.coef = 0.25;
            a = 1.0;
        }
        else if (it == 2)
        {
            ch.coef = 0.25;
            a = 0.5;
        }
        else if (it == 3)
        {
            ch.coef = 0.5;
            a = 0.75;
        }
        else
        {
            ch.coef = 0.0;
            a = 1.0;
        }
        ch.w = pow_g(nc, a);
        ch.dw = pow_dg(nc, a);
        break;
    }
    case XcType::HybOpt:
        if (it == 1)
        {
            ch.coef = 1.0 - HYBOPT_POWER_W;
            ch.w = nc;
            ch.dw = 1.0;
        }
        else
        {
            ch.coef = HYBOPT_POWER_W;
            ch.w = pow_g(nc, HYBOPT_POWER_EXP);
            ch.dw = pow_dg(nc, HYBOPT_POWER_EXP);
        }
        break;
    case XcType::BOWMOD:
        if (it == 1)
        {
            ch.coef = BOWMOD_COEF_HF;
            ch.w = nc;
            ch.dw = 1.0;
        }
        else if (it == 2)
        {
            ch.coef = BOWMOD_COEF_ALPHA;
            ch.w = pow_g(nc, alpha_);
            ch.dw = pow_dg(nc, alpha_);
        }
        else if (it == 3)
        {
            ch.coef = BOWMOD_COEF_ALPHA;
            ch.w = pow_g(1.0 - nc, alpha_);
            ch.dw = -pow_dg(1.0 - nc, alpha_);
        }
        else
        {
            ch.coef = BOWMOD_COEF_MIX;
            const double s = std::max(0.0, nc * (1.0 - nc));
            ch.w = pow_g(s, BOWMOD_DEFAULT_BETA);
            ch.dw = pow_dg(s, BOWMOD_DEFAULT_BETA) * (1.0 - 2.0 * nc);
        }
        break;
    default:
        ch.coef = 1.0;
        ch.w = g(nc);
        ch.dw = dg(nc);
        break;
    }
    return ch;
}

bool RdmftXC::has_extra_diag() const
{
    return type_ == XcType::GU;
}

double RdmftXC::gu_diag_factor(double n) const
{
    const double nc = clip01(n);
    return nc * nc - nc;
}

double RdmftXC::gu_diag_factor_deriv(double n) const
{
    const double nc = clip01(n);
    return 2.0 * nc - 1.0;
}

double RdmftXC::binary_entropy(double n)
{
    const double eps = 1.0e-12;
    const double nc = std::max(eps, std::min(1.0 - eps, n));
    return nc * std::log(nc) + (1.0 - nc) * std::log(1.0 - nc);
}

double RdmftXC::binary_entropy_deriv(double n)
{
    const double eps = 1.0e-12;
    const double nc = std::max(eps, std::min(1.0 - eps, n));
    return std::log(nc / (1.0 - nc));
}

} // namespace rdmft_core
