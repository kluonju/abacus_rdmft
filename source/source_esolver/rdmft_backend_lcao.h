#ifndef RDMFT_BACKEND_LCAO_H
#define RDMFT_BACKEND_LCAO_H

#include "source_rdmft/rdmft_backend.h"
#include "source_lcao/module_rdmft/rdmft_energy_gradient.h"
#include "source_psi/psi.h"
#include "source_cell/klist.h"

#include <complex>
#include <vector>

//! \file rdmft_backend_lcao.h
//! \brief LCAO adapter for the basis-independent RDMFT optimisation core.
//!
//! Implements rdmft::RdmftBackend by delegating the energy/gradient physics to
//! the existing (validated) LCAO engine rdmft::EnergyGradient<TK,TR>, while the
//! optimisation itself (SPG2/EBI occupations, Stiefel orbitals) is driven by
//! the new modular core.  The engine works in Cholesky X-space
//! (X_k = U_k C_k, S_k = U_k^H U_k), so the orbital manifold reduces to the
//! standard Stiefel manifold X^H X = I with the canonical metric provided by
//! EnergyGradient::stiefel_canonical_inner_product.
//!
//! Occupation convention.  The core enforces per-spin natural occupations in
//! [0,1] with the spin degeneracy folded into the BZ weights (QE convention).
//! ABACUS stores wg = wk * f with f in [0,2] for nspin=1.  The adapter folds a
//! factor spin_deg (=2 for nspin=1, else 1) into wk so n = wg/(spin_deg*wk) in
//! [0,1], and converts occupations/gradients back to the engine's convention
//! (occ_flat = spin_deg * n, dE/dn = spin_deg * dE/docc_flat).

namespace ModuleESolver
{

namespace rdmft_detail
{
// Flatten / restore psi coefficient buffers as a real double vector.
inline void psi_to_flat(const psi::Psi<double>& p, std::vector<double>& f)
{
    f.assign(p.get_pointer(), p.get_pointer() + p.size());
}
inline void psi_to_flat(const psi::Psi<std::complex<double>>& p, std::vector<double>& f)
{
    const std::complex<double>* d = p.get_pointer();
    const size_t n = p.size();
    f.resize(2 * n);
    for (size_t i = 0; i < n; ++i)
    {
        f[2 * i] = d[i].real();
        f[2 * i + 1] = d[i].imag();
    }
}
inline void flat_to_psi(const std::vector<double>& f, psi::Psi<double>& p)
{
    double* d = p.get_pointer();
    for (size_t i = 0; i < p.size(); ++i)
    {
        d[i] = f[i];
    }
}
inline void flat_to_psi(const std::vector<double>& f, psi::Psi<std::complex<double>>& p)
{
    std::complex<double>* d = p.get_pointer();
    for (size_t i = 0; i < p.size(); ++i)
    {
        d[i] = std::complex<double>(f[2 * i], f[2 * i + 1]);
    }
}
} // namespace rdmft_detail

template <typename TK, typename TR>
class RdmftBackendLCAO : public rdmft_core::RdmftBackend
{
  public:
    //! eg must already be init()'d and update_ion()'d.  psi holds the KS
    //! orbitals in C-space; this adapter transforms them to X-space in place.
    RdmftBackendLCAO(rdmft::EnergyGradient<TK, TR>& eg, psi::Psi<TK>& psi, const K_Vectors& kv,
                     int nspin, double nelec, bool fix_mag, double nelec_up, double nelec_down)
        : eg_(eg), psi_(psi), grad_wfc_(psi), scratch_(psi), saved_(psi)
    {
        const int nk = eg_.nk();
        const int nbands = eg_.nbands();
        const int nks = kv.get_nks();
        (void)nspin;

        con_.nbnd = nbands;
        con_.nks = nk;
        con_.wk.resize(nk);
        con_.isk.resize(nk);
        for (int ik = 0; ik < nk; ++ik)
        {
            const int kslot = (ik < nks) ? ik : ik - nks;
            // kv.wk already folds the spin degeneracy (for nspin=1 it carries the
            // factor 2), so n = wg/kv.wk is the per-spin natural occupation in
            // [0,1] and equals the occupation EnergyGradient expects (wg/wk).
            con_.wk[ik] = kv.wk[kslot];
            con_.isk[ik] = (nspin == 2 && ik >= nks) ? 2 : 1;
        }
        con_.n_target = nelec;
        con_.fix_magnetization = fix_mag;
        con_.n_target_up = nelec_up;
        con_.n_target_down = nelec_down;

        // Move to X-space: the manifold becomes X^H X = I.
        eg_.precompute_cholesky_S();
        eg_.wfc_C_to_X(psi_);
        occ_scratch_.assign(con_.size(), 0.0);
        grad_scratch_.assign(con_.size(), 0.0);
    }

    const rdmft_core::OccConstraints& occ_constraints() const override { return con_; }

    double total_energy(const std::vector<double>& occ) override
    {
        to_engine_occ(occ, occ_scratch_);
        return eg_.compute_energy(occ_scratch_, psi_);
    }

    void grad_occ(const std::vector<double>& occ, std::vector<double>& grad) override
    {
        to_engine_occ(occ, occ_scratch_);
        eg_.compute(occ_scratch_, psi_, grad_scratch_, grad_wfc_);
        grad.resize(con_.size());
        for (int i = 0; i < con_.size(); ++i)
        {
            // occ_flat == n (see to_engine_occ), so dE/dn = dE/docc_flat.
            grad[i] = grad_scratch_[i];
        }
    }

    bool has_orbital_optimization() const override { return true; }

    int orb_dim() const override
    {
        return static_cast<int>(orb_flat_size(psi_));
    }

    double riemannian_gradient(const std::vector<double>& occ, std::vector<double>& gR) override
    {
        to_engine_occ(occ, occ_scratch_);
        eg_.compute(occ_scratch_, psi_, grad_scratch_, grad_wfc_);
        eg_.project_orbital_gradient(psi_, grad_wfc_);
        rdmft_detail::psi_to_flat(grad_wfc_, gR);
        return eg_.stiefel_canonical_inner_product(psi_, grad_wfc_, grad_wfc_);
    }

    double orb_inner(const std::vector<double>& a, const std::vector<double>& b) override
    {
        // Use dedicated scratch buffers so orb_save()/orb_restore() state is
        // never clobbered by an inner product evaluated during a line search.
        rdmft_detail::flat_to_psi(a, grad_wfc_);
        rdmft_detail::flat_to_psi(b, scratch_);
        return eg_.stiefel_canonical_inner_product(psi_, grad_wfc_, scratch_);
    }

    void orb_project_tangent(std::vector<double>& v) override
    {
        rdmft_detail::flat_to_psi(v, grad_wfc_);
        eg_.project_orbital_gradient(psi_, grad_wfc_);
        rdmft_detail::psi_to_flat(grad_wfc_, v);
    }

    void orb_retract(const std::vector<double>& dir, double alpha) override
    {
        rdmft_detail::flat_to_psi(dir, grad_wfc_);
        // retract_orbitals moves C <- C - alpha*G; pass -alpha to advance along +dir.
        eg_.retract_orbitals(psi_, grad_wfc_, -alpha);
        eg_.invalidate_hone_cache();
    }

    void orb_save() override { copy_psi(psi_, saved_); }
    void orb_restore() override
    {
        copy_psi(saved_, psi_);
        eg_.invalidate_hone_cache();
    }

    //! Convert the optimised X-space orbitals back to C-space (call once done).
    void finalize_orbitals() { eg_.wfc_X_to_C(psi_); }

  private:
    static size_t orb_flat_size(const psi::Psi<double>& p) { return p.size(); }
    static size_t orb_flat_size(const psi::Psi<std::complex<double>>& p) { return 2 * p.size(); }

    static void copy_psi(const psi::Psi<TK>& src, psi::Psi<TK>& dst)
    {
        dst.set_all_psi(src.get_pointer(), src.size());
    }

    void to_engine_occ(const std::vector<double>& n, std::vector<double>& occ_flat) const
    {
        occ_flat.resize(con_.size());
        for (int i = 0; i < con_.size(); ++i)
        {
            // n is already the per-spin occupation wg/kv.wk expected by EnergyGradient.
            occ_flat[i] = n[i];
        }
    }

    rdmft::EnergyGradient<TK, TR>& eg_;
    psi::Psi<TK>& psi_;
    psi::Psi<TK> grad_wfc_;
    psi::Psi<TK> scratch_;
    psi::Psi<TK> saved_;
    rdmft_core::OccConstraints con_;
    std::vector<double> occ_scratch_;
    std::vector<double> grad_scratch_;
};

} // namespace ModuleESolver

#endif // RDMFT_BACKEND_LCAO_H
