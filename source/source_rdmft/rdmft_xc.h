#ifndef SOURCE_RDMFT_RDMFT_XC_H
#define SOURCE_RDMFT_RDMFT_XC_H

#include "rdmft_params.h"

//! \file rdmft_xc.h
//! \brief Natural-orbital exchange-correlation coupling kernels.
//!
//! Ported from qe-rdmft/rdmft_xc.f90.  The functionals are expressed as a
//! sum of separable exchange channels
//!
//!   E_xc[gamma] = sum_t coef_t * E_x^HF[gamma_t],
//!
//! where gamma_t is a "modified" one-body density matrix built with the
//! per-band weight w_t(n_ik).  Each channel therefore exposes a
//! (coef_t, w_t(n), w_t'(n)) triplet so the basis-specific backend can
//! reuse its exact-exchange builder.  Both the plane-wave and LCAO backends
//! consume exactly this interface, which is what makes the functional layer
//! basis independent.

namespace rdmft
{

//! Channel descriptor: E contribution is coef * 0.5 * sum_ik wk*w * <psi|Vx[gamma_t]|psi>,
//! and dE/dn contribution is coef * wk * dw * <psi|Vx[gamma_t]|psi>.
struct XcChannel
{
    double coef = 0.0; //!< constant prefactor in the energy sum
    double w = 0.0;    //!< modified-DM weight w_t(n)
    double dw = 0.0;   //!< derivative w_t'(n)
};

//! Analytic natural-orbital XC coupling functions.
class RdmftXC
{
  public:
    RdmftXC() = default;
    RdmftXC(XcType type, double alpha, double reg_eps);

    void set(XcType type, double alpha, double reg_eps);

    XcType type() const { return type_; }
    double alpha() const { return alpha_; }

    //! Single-channel separable coupling g(n) and derivative dg(n).
    double g(double n) const;
    double dg(double n) const;

    //! Whether the functional is a single separable channel f = g(ni)g(nj).
    bool is_separable() const;

    //! Number of separable exchange channels the functional decomposes into.
    int n_channels() const;

    //! Evaluate channel it (1-based) at occupation n.
    XcChannel channel(int it, double n) const;

    //! GU-type on-site extra: energy factor (n^2-n) and its derivative (2n-1).
    bool has_extra_diag() const;
    double gu_diag_factor(double n) const;
    double gu_diag_factor_deriv(double n) const;

    //! Binary entropy helpers used by the optional finite-temperature term.
    static double binary_entropy(double n);
    static double binary_entropy_deriv(double n);

  private:
    double pow_g(double n, double a) const;
    double pow_dg(double n, double a) const;

    XcType type_ = XcType::Muller;
    double alpha_ = 0.5;
    double reg_eps_ = 1.0e-8;
};

} // namespace rdmft

#endif // SOURCE_RDMFT_RDMFT_XC_H
