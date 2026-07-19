#include "esolver_rdmft_pw.h"

#include "source_base/global_variable.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"

#ifdef __RDMFT
#include "source_base/parallel_reduce.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_pw/module_pwdft/exx_helper.h"
#include "source_pw/module_pwdft/op_pw_exx.h"
#include "source_rdmft/rdmft_backend.h"
#include "source_rdmft/rdmft_driver.h"
#include "source_rdmft/rdmft_params.h"
#include "source_rdmft/rdmft_xc.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>
#endif

namespace ModuleESolver
{

template <typename T, typename Device>
ESolver_RDMFT_PW<T, Device>::ESolver_RDMFT_PW()
{
    this->classname = "ESolver_RDMFT_PW";
}

template <typename T, typename Device>
ESolver_RDMFT_PW<T, Device>::~ESolver_RDMFT_PW()
{
}

#ifdef __RDMFT
namespace
{
//! Plane-wave RDMFT backend: occupation optimisation with frozen natural
//! orbitals, exact exchange via the ACE operator.  Energy is assembled by
//! reusing the KS density/potential/energy machinery:
//!   E = E_bandlike - E_Hartree + E_xc + E_ewald,
//! where E_bandlike = sum_ik w_k n_ik <psi|H_KS(noEXX)|psi> (H_KS = T+V_nl+
//! V_loc+V_H for the HF XC context), and E_xc is the ACE exchange energy with
//! the modified occupation weights w_t(n).  The double-counting bookkeeping
//! follows E_one + E_H = E_bandlike - E_H (since <psi|V_H|psi> summed = 2 E_H).
template <typename T, typename Device>
class RdmftBackendPW : public rdmft_core::RdmftBackend
{
  public:
    RdmftBackendPW(hamilt::Hamilt<T, Device>* p_hamilt, Exx_HelperBase* exx_helper,
                   psi::Psi<T, Device>* psi, elecstate::ElecState* pelec, const K_Vectors& kv,
                   UnitCell& ucell, const rdmft_core::RdmftXC& xc, double hybrid_alpha, int nspin,
                   double nelec, bool fix_mag, double nelec_up, double nelec_down)
        : p_hamilt_(p_hamilt), exx_helper_(exx_helper), psi_(psi), pelec_(pelec), kv_(kv),
          ucell_(ucell), xc_(xc), hybrid_alpha_(hybrid_alpha)
    {
        nbands_ = psi_->get_nbands();
        nks_ = kv.get_nks();
        npol_ = (nspin == 4) ? 2 : 1;

        con_.nbnd = nbands_;
        con_.nks = nks_;
        con_.wk.resize(nks_);
        con_.isk.resize(nks_);
        for (int ik = 0; ik < nks_; ++ik)
        {
            // kv.wk already folds the spin degeneracy (for nspin=1 it carries the
            // factor 2), so the natural occupation n = wg/kv.wk lies in [0,1] and
            // the modified-DM weight is kv.wk * w_t(n).  Multiplying by an extra
            // spin factor would double-count it and break nonlinear functionals.
            con_.wk[ik] = kv.wk[ik];
            // kv.isk is 0-based (0=up, 1=down); OccConstraints uses 1/2.
            con_.isk[ik] = (ik < static_cast<int>(kv.isk.size())) ? kv.isk[ik] + 1 : 1;
        }
        con_.n_target = nelec;
        con_.fix_magnetization = fix_mag;
        con_.n_target_up = nelec_up;
        con_.n_target_down = nelec_down;

        wg_.create(nks_, nbands_);
        wg_x_.create(nks_, nbands_);
        diag_.assign(nks_ * nbands_, 0.0);
        vx_diag_.assign(nks_ * nbands_, 0.0);
        const size_t full = static_cast<size_t>(nks_) * nbands_ * psi_->get_nbasis();
        hpsi_all_.assign(full, T(0.0, 0.0));
        vxpsi_all_.assign(full, T(0.0, 0.0));
    }

    ~RdmftBackendPW() override
    {
        // Restore the EXX operator's weight pointer (we repointed it at wg_x_).
        Exx_Helper<T, Device>* helper = static_cast<Exx_Helper<T, Device>*>(exx_helper_);
        if (helper != nullptr && helper->op_exx != nullptr)
        {
            helper->op_exx->set_wg(&pelec_->wg);
        }
    }

    const rdmft_core::OccConstraints& occ_constraints() const override { return con_; }

    double total_energy(const std::vector<double>& occ) override
    {
        refresh_state(occ);
        // E_bandlike = sum w_k n_ik <psi|H_KS(noEXX)|psi>.
        double e_bandlike = 0.0;
        for (int ik = 0; ik < nks_; ++ik)
        {
            for (int ib = 0; ib < nbands_; ++ib)
            {
                e_bandlike += con_.wk[ik] * occ[ik * nbands_ + ib] * diag_[ik * nbands_ + ib];
            }
        }
        const double e_h = pelec_->f_en.hartree_energy;
        const double e_ewald = pelec_->f_en.ewald_energy;
        // E_bandlike = E_one + 2 E_H + vtxc (Veff = V_loc + V_H + V_xc), so the
        // RDMFT energy (which replaces semilocal XC by the natural-orbital
        // functional) is E_one + E_H + E_x + E_ewald = E_bandlike - E_H - vtxc
        // + E_x + E_ewald.  vtxc = 0 in the pure-HF XC context used here, so this
        // reduces to the HF assembly but stays correct if any semilocal V_xc is
        // present.
        const double vtxc = pelec_->f_en.vtxc;
        const double e_x = exchange_energy(occ);
        return e_bandlike - e_h - vtxc + e_x + e_ewald;
    }

    void grad_occ(const std::vector<double>& occ, std::vector<double>& grad) override
    {
        refresh_state(occ);
        compute_exchange_action(occ); // fills vx_diag_ + vxpsi_all_
        grad.assign(con_.size(), 0.0);
        const int nch = xc_.n_channels();
        for (int ik = 0; ik < nks_; ++ik)
        {
            const double wk = con_.wk[ik];
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const int i = ik * nbands_ + ib;
                double gx = 0.0;
                for (int it = 1; it <= nch; ++it)
                {
                    const rdmft_core::XcChannel ch = xc_.channel(it, occ[i]);
                    gx += ch.coef * ch.dw;
                }
                grad[i] = wk * (diag_[i] + gx * vx_diag_[i]);
            }
        }
    }

    // --- orbital (Stiefel, S = I) oracle ------------------------------------
    bool has_orbital_optimization() const override { return true; }

    int orb_dim() const override { return 2 * static_cast<int>(psi_->size()); }

    double riemannian_gradient(const std::vector<double>& occ, std::vector<double>& gR) override
    {
        refresh_state(occ);            // fills hpsi_all_
        compute_exchange_action(occ);  // fills vxpsi_all_
        gR.assign(orb_dim(), 0.0);
        const int nbasis = psi_->get_nbasis();
        const int npwx = nbasis / npol_;
        std::vector<T> Gk(static_cast<size_t>(nbands_) * nbasis);
        for (int ik = 0; ik < nks_; ++ik)
        {
            psi_->fix_k(ik);
            T* Cbuf = psi_->get_pointer();
            const int npwk = psi_->get_current_nbas();
            const double wk = con_.wk[ik];
            std::fill(Gk.begin(), Gk.end(), T(0.0, 0.0));
            const size_t koff = static_cast<size_t>(ik) * nbands_ * nbasis;
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const double n = occ[ik * nbands_ + ib];
                const double gc = xc_.g(n); // single-channel exchange coupling
                const size_t boff = static_cast<size_t>(ib) * nbasis;
                for (int ipol = 0; ipol < npol_; ++ipol)
                {
                    for (int g = 0; g < npwk; ++g)
                    {
                        const size_t idx = boff + static_cast<size_t>(ipol) * npwx + g;
                        Gk[idx] = wk * (n * hpsi_all_[koff + idx] + gc * vxpsi_all_[koff + idx]);
                    }
                }
            }
            project_tangent_k(Cbuf, Gk.data(), npwk);
            // pack Gk (this k-block) into gR
            for (size_t j = 0; j < static_cast<size_t>(nbands_) * nbasis; ++j)
            {
                const size_t full = koff + j;
                gR[2 * full] = Gk[j].real();
                gR[2 * full + 1] = Gk[j].imag();
            }
        }
        return orb_inner(gR, gR);
    }

    double orb_inner(const std::vector<double>& a, const std::vector<double>& b) override
    {
        double s = 0.0;
        for (size_t i = 0; i < a.size(); ++i)
        {
            s += a[i] * b[i];
        }
#ifdef __MPI
        double sg = 0.0;
        MPI_Allreduce(&s, &sg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        s = sg;
#endif
        return s;
    }

    void orb_project_tangent(std::vector<double>& v) override
    {
        const int nbasis = psi_->get_nbasis();
        std::vector<T> Vk(static_cast<size_t>(nbands_) * nbasis);
        for (int ik = 0; ik < nks_; ++ik)
        {
            psi_->fix_k(ik);
            T* Cbuf = psi_->get_pointer();
            const int npwk = psi_->get_current_nbas();
            const size_t koff = static_cast<size_t>(ik) * nbands_ * nbasis;
            for (size_t j = 0; j < static_cast<size_t>(nbands_) * nbasis; ++j)
            {
                Vk[j] = T(v[2 * (koff + j)], v[2 * (koff + j) + 1]);
            }
            project_tangent_k(Cbuf, Vk.data(), npwk);
            for (size_t j = 0; j < static_cast<size_t>(nbands_) * nbasis; ++j)
            {
                v[2 * (koff + j)] = Vk[j].real();
                v[2 * (koff + j) + 1] = Vk[j].imag();
            }
        }
    }

    void orb_retract(const std::vector<double>& dir, double alpha) override
    {
        const int nbasis = psi_->get_nbasis();
        const int npwx = nbasis / npol_;
        for (int ik = 0; ik < nks_; ++ik)
        {
            psi_->fix_k(ik);
            T* Cbuf = psi_->get_pointer();
            const int npwk = psi_->get_current_nbas();
            const size_t koff = static_cast<size_t>(ik) * nbands_ * nbasis;
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const size_t boff = static_cast<size_t>(ib) * nbasis;
                for (int ipol = 0; ipol < npol_; ++ipol)
                {
                    for (int g = 0; g < npwk; ++g)
                    {
                        const size_t idx = boff + static_cast<size_t>(ipol) * npwx + g;
                        Cbuf[idx] += T(alpha * dir[2 * (koff + idx)], alpha * dir[2 * (koff + idx) + 1]);
                    }
                }
            }
            orthonormalize_k(Cbuf, npwk);
        }
    }

    void orb_save() override
    {
        // The per-k loops above leave psi's cursor (current_k / psi_bias) at the
        // last k-point, so get_pointer() no longer points at the buffer start.
        // Reset it before copying the full [0, size()) buffer, otherwise this
        // reads/writes past the allocation for nks_ > 1 (multi-k or nspin=2).
        psi_->fix_k(0);
        saved_psi_.assign(psi_->get_pointer(), psi_->get_pointer() + psi_->size());
    }
    void orb_restore() override
    {
        psi_->fix_k(0);
        std::copy(saved_psi_.begin(), saved_psi_.end(), psi_->get_pointer());
    }

  private:
    // Sync density/potential and recompute one-body action/diagonals for occ.
    void refresh_state(const std::vector<double>& occ)
    {
        for (int ik = 0; ik < nks_; ++ik)
        {
            for (int ib = 0; ib < nbands_; ++ib)
            {
                wg_(ik, ib) = con_.wk[ik] * occ[ik * nbands_ + ib];
            }
        }
        pelec_->wg = wg_;
        pelec_->psiToRho(*psi_);
        pelec_->pot->update_from_charge(pelec_->charge, &ucell_);
        pelec_->cal_energies(2);
        compute_onebody();
    }

    // hpsi_all_[koff+boff+idx] = (H_KS(no EXX) psi)_ig; diag_ = Re <psi|H|psi>.
    void compute_onebody()
    {
        const bool prev_first = exx_helper_->get_op_first_iter();
        exx_helper_->set_op_first_iter(true); // disable EXX in the operator chain
        const int nbasis = psi_->get_nbasis();
        std::vector<T> hpsi(static_cast<size_t>(nbasis) * nbands_);
        for (int ik = 0; ik < nks_; ++ik)
        {
            p_hamilt_->updateHk(ik);
            psi_->fix_k(ik);
            T* psi_k = psi_->get_pointer();
            const int npw = psi_->get_current_nbas();
            std::fill(hpsi.begin(), hpsi.end(), T(0.0, 0.0));
            psi::Psi<T, Device> wrapper(psi_k, 1, nbands_, nbasis, npw);
            psi::Range bands(true, 0, 0, nbands_ - 1);
            typename hamilt::Operator<T, Device>::hpsi_info info(&wrapper, bands, hpsi.data());
            p_hamilt_->ops->hPsi(info);
            const size_t koff = static_cast<size_t>(ik) * nbands_ * nbasis;
            store_action(hpsi, psi_k, npw, koff, hpsi_all_, diag_, ik);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(diag_.data(), static_cast<int>(diag_.size()));
#endif
        exx_helper_->set_op_first_iter(prev_first);
    }

    // Point the EXX operator at the modified-DM weights wg_x = wk * w_t(n) and
    // rebuild the ACE projectors for the current orbitals.  Exx_Helper::set_wg
    // only updates the helper's own pointer, so the operator's weight pointer
    // (used by construct_ace / cal_exx_energy) must be set explicitly.
    hamilt::OperatorEXXPW<T, Device>* prepare_exx(const std::vector<double>& occ)
    {
        set_exchange_weights(occ);
        Exx_Helper<T, Device>* helper = static_cast<Exx_Helper<T, Device>*>(exx_helper_);
        hamilt::OperatorEXXPW<T, Device>* op = helper->op_exx;
        if (op != nullptr)
        {
            op->set_wg(&wg_x_);
            op->first_iter = false;
        }
        exx_helper_->set_wg(&wg_x_);
        exx_helper_->set_op_first_iter(false);
        exx_helper_->set_psi(psi_); // op->set_psi + construct_ace(wg_x, current orbitals)
        return op;
    }

    // Build ACE exchange with weights w_t(n) and return the exchange energy.
    double exchange_energy(const std::vector<double>& occ)
    {
        prepare_exx(occ);
        return hybrid_alpha_ * exx_helper_->cal_exx_energy(psi_);
    }

    // vxpsi_all_ = Vx[w_t(n)] psi (hybrid_alpha-scaled); vx_diag_ its diagonal.
    void compute_exchange_action(const std::vector<double>& occ)
    {
        hamilt::OperatorEXXPW<T, Device>* op = prepare_exx(occ);
        std::fill(vx_diag_.begin(), vx_diag_.end(), 0.0);
        std::fill(vxpsi_all_.begin(), vxpsi_all_.end(), T(0.0, 0.0));
        if (op == nullptr)
        {
            return;
        }
        const int nbasis = psi_->get_nbasis();
        std::vector<T> vxpsi(static_cast<size_t>(nbasis) * nbands_);
        for (int ik = 0; ik < nks_; ++ik)
        {
            p_hamilt_->updateHk(ik);
            psi_->fix_k(ik);
            T* psi_k = psi_->get_pointer();
            const int npw = psi_->get_current_nbas();
            std::fill(vxpsi.begin(), vxpsi.end(), T(0.0, 0.0));
            op->act(nbands_, nbasis, 1, psi_k, vxpsi.data(), npw, true);
            const size_t koff = static_cast<size_t>(ik) * nbands_ * nbasis;
            store_action(vxpsi, psi_k, npw, koff, vxpsi_all_, vx_diag_, ik);
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(vx_diag_.data(), static_cast<int>(vx_diag_.size()));
#endif
        // op->act already applies hybrid_alpha internally.
    }

    // Store the operator action into the full buffer (zeroing padding) and the
    // per-band diagonal Re <psi|A|psi> over the active plane waves.
    void store_action(const std::vector<T>& apsi, const T* psi_k, int npw, size_t koff,
                      std::vector<T>& dst, std::vector<double>& diag, int ik)
    {
        const int nbasis = psi_->get_nbasis();
        const int npwx = nbasis / npol_;
        for (int ib = 0; ib < nbands_; ++ib)
        {
            const size_t boff = static_cast<size_t>(ib) * nbasis;
            const T* pc = psi_k + boff;
            const T* ac = apsi.data() + boff;
            double d = 0.0;
            for (int ipol = 0; ipol < npol_; ++ipol)
            {
                for (int g = 0; g < npw; ++g)
                {
                    const size_t idx = boff + static_cast<size_t>(ipol) * npwx + g;
                    dst[koff + idx] = apsi[idx];
                    d += (std::conj(pc[static_cast<size_t>(ipol) * npwx + g])
                          * ac[static_cast<size_t>(ipol) * npwx + g])
                             .real();
                }
            }
            diag[ik * nbands_ + ib] = d;
        }
    }

    void set_exchange_weights(const std::vector<double>& occ)
    {
        // Modified-DM exchange weights wg_x = wk * max_t w_t(n), exactly as
        // qe-rdmft's rdmft_xc_set_exx_wg (max over channels of the coupling
        // weight w_t(n)).  Following qe, the main energy/gradient ACE build is
        // NOT trimmed: the regularised small-n tail of w_t(n) (w_t(0) != 0 for
        // Muller/Power) is kept for every band, so the PW/ACE exchange matches
        // qe's plane-wave ACE.  (qe only trims near-empty bands in its TSM/DOS
        // probe path, rdmft_xc_tsm_trim_wg_for_ace, not here.)
        const int nch = xc_.n_channels();
        for (int ik = 0; ik < nks_; ++ik)
        {
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const double n = occ[ik * nbands_ + ib];
                double w = 0.0;
                for (int it = 1; it <= nch; ++it)
                {
                    w = std::max(w, xc_.channel(it, n).w);
                }
                wg_x_(ik, ib) = con_.wk[ik] * w;
            }
        }
    }

    // Gram matrix A[a*nb+b] = sum_active conj(C[a]) D[b], reduced over the pool.
    void gram(const T* Cbuf, const T* Dbuf, int npwk, std::vector<T>& A)
    {
        const int nbasis = psi_->get_nbasis();
        const int npwx = nbasis / npol_;
        A.assign(static_cast<size_t>(nbands_) * nbands_, T(0.0, 0.0));
        for (int a = 0; a < nbands_; ++a)
        {
            const T* ca = Cbuf + static_cast<size_t>(a) * nbasis;
            for (int b = 0; b < nbands_; ++b)
            {
                const T* db = Dbuf + static_cast<size_t>(b) * nbasis;
                T s(0.0, 0.0);
                for (int ipol = 0; ipol < npol_; ++ipol)
                {
                    for (int g = 0; g < npwk; ++g)
                    {
                        const size_t idx = static_cast<size_t>(ipol) * npwx + g;
                        s += std::conj(ca[idx]) * db[idx];
                    }
                }
                A[static_cast<size_t>(a) * nbands_ + b] = s;
            }
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(reinterpret_cast<double*>(A.data()), 2 * nbands_ * nbands_);
#endif
    }

    // eta = G - C sym(C^H G), sym(A) = (A + A^H)/2.
    void project_tangent_k(const T* Cbuf, T* Gbuf, int npwk)
    {
        std::vector<T> A;
        gram(Cbuf, Gbuf, npwk, A);
        std::vector<T> As(static_cast<size_t>(nbands_) * nbands_);
        for (int a = 0; a < nbands_; ++a)
        {
            for (int b = 0; b < nbands_; ++b)
            {
                As[a * nbands_ + b]
                    = 0.5 * (A[a * nbands_ + b] + std::conj(A[b * nbands_ + a]));
            }
        }
        const int nbasis = psi_->get_nbasis();
        const int npwx = nbasis / npol_;
        for (int b = 0; b < nbands_; ++b)
        {
            for (int ipol = 0; ipol < npol_; ++ipol)
            {
                for (int g = 0; g < npwk; ++g)
                {
                    const size_t idx = static_cast<size_t>(ipol) * npwx + g;
                    T acc(0.0, 0.0);
                    for (int a = 0; a < nbands_; ++a)
                    {
                        acc += Cbuf[static_cast<size_t>(a) * nbasis + idx] * As[a * nbands_ + b];
                    }
                    Gbuf[static_cast<size_t>(b) * nbasis + idx] -= acc;
                }
            }
        }
    }

    // In-place Cholesky orthonormalisation Y <- Y (L^H)^{-1}, Y^H Y = L L^H.
    void orthonormalize_k(T* Ybuf, int npwk)
    {
        std::vector<T> M;
        gram(Ybuf, Ybuf, npwk, M);
        // Cholesky lower L (M = L L^H).
        std::vector<T> L(static_cast<size_t>(nbands_) * nbands_, T(0.0, 0.0));
        for (int j = 0; j < nbands_; ++j)
        {
            T s = M[j * nbands_ + j];
            for (int k = 0; k < j; ++k)
            {
                s -= L[j * nbands_ + k] * std::conj(L[j * nbands_ + k]);
            }
            double djj = std::sqrt(std::max(1.0e-300, s.real()));
            L[j * nbands_ + j] = T(djj, 0.0);
            for (int i = j + 1; i < nbands_; ++i)
            {
                T t = M[i * nbands_ + j];
                for (int k = 0; k < j; ++k)
                {
                    t -= L[i * nbands_ + k] * std::conj(L[j * nbands_ + k]);
                }
                L[i * nbands_ + j] = t / djj;
            }
        }
        // Lower-triangular inverse Linv (L Linv = I).
        std::vector<T> Linv(static_cast<size_t>(nbands_) * nbands_, T(0.0, 0.0));
        for (int i = 0; i < nbands_; ++i)
        {
            Linv[i * nbands_ + i] = T(1.0, 0.0) / L[i * nbands_ + i];
            for (int j = 0; j < i; ++j)
            {
                T acc(0.0, 0.0);
                for (int k = j; k < i; ++k)
                {
                    acc += L[i * nbands_ + k] * Linv[k * nbands_ + j];
                }
                Linv[i * nbands_ + j] = -acc / L[i * nbands_ + i];
            }
        }
        // Y_new[:,b] = sum_a Y[:,a] conj(Linv[b,a]).
        const int nbasis = psi_->get_nbasis();
        const int npwx = nbasis / npol_;
        std::vector<T> col(nbands_);
        for (int ipol = 0; ipol < npol_; ++ipol)
        {
            for (int g = 0; g < npwk; ++g)
            {
                const size_t idx = static_cast<size_t>(ipol) * npwx + g;
                for (int b = 0; b < nbands_; ++b)
                {
                    T acc(0.0, 0.0);
                    for (int a = 0; a < nbands_; ++a)
                    {
                        acc += Ybuf[static_cast<size_t>(a) * nbasis + idx]
                               * std::conj(Linv[static_cast<size_t>(b) * nbands_ + a]);
                    }
                    col[b] = acc;
                }
                for (int b = 0; b < nbands_; ++b)
                {
                    Ybuf[static_cast<size_t>(b) * nbasis + idx] = col[b];
                }
            }
        }
    }

    hamilt::Hamilt<T, Device>* p_hamilt_;
    Exx_HelperBase* exx_helper_;
    psi::Psi<T, Device>* psi_;
    elecstate::ElecState* pelec_;
    const K_Vectors& kv_;
    UnitCell& ucell_;
    rdmft_core::RdmftXC xc_;
    double hybrid_alpha_;
    int nbands_ = 0;
    int nks_ = 0;
    int npol_ = 1;
    rdmft_core::OccConstraints con_;
    ModuleBase::matrix wg_;
    ModuleBase::matrix wg_x_;
    std::vector<double> diag_;
    std::vector<double> vx_diag_;
    std::vector<T> hpsi_all_;
    std::vector<T> vxpsi_all_;
    std::vector<T> saved_psi_;
};

rdmft_core::OccOptimizerType map_occ_optimizer_pw(const std::string& s)
{
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
    if (lower.find("ebi") != std::string::npos)
    {
        return rdmft_core::OccOptimizerType::EBI;
    }
    if (lower.find("bgd") != std::string::npos || lower == "gd")
    {
        return rdmft_core::OccOptimizerType::BGD;
    }
    return rdmft_core::OccOptimizerType::SPG2;
}

rdmft_core::OrbOptimizerType map_orb_optimizer_pw(const std::string& s)
{
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
    if (lower == "sd")
    {
        return rdmft_core::OrbOptimizerType::SD;
    }
    if (lower.find("lbfgs") != std::string::npos || lower.find("bfgs") != std::string::npos)
    {
        return rdmft_core::OrbOptimizerType::LBFGS;
    }
    return rdmft_core::OrbOptimizerType::CG;
}
} // namespace
#endif

template <typename T, typename Device>
void ESolver_RDMFT_PW<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ESolver_KS_PW<T, Device>::after_scf(ucell, istep, conv_esolver);

#ifdef __RDMFT
    if (this->exx_helper == nullptr || !GlobalC::exx_info.info_global.cal_exx)
    {
        GlobalV::ofs_running << "RDMFT(PW): exact exchange is not enabled; skipping RDMFT."
                             << std::endl;
        return;
    }
    const Input_para& inp = PARAM.inp;
    psi::Psi<T, Device>* psi = this->stp.template get_psi_t<T, Device>();
    if (psi == nullptr)
    {
        return;
    }

    if (PARAM.inp.nspin == 4)
    {
        ModuleBase::WARNING_QUIT("ESolver_RDMFT_PW",
                                 "non-collinear (nspin=4) RDMFT with the plane-wave ACE backend is "
                                 "under development.");
    }

    rdmft_core::XcType xc_type = rdmft_core::XcType::HF;
    double alpha = inp.rdmft_power_alpha;
    rdmft_core::parse_xc_type(inp.rdmft_functional, xc_type, alpha);
    const rdmft_core::RdmftXC xc(xc_type, alpha, 1.0e-8);

    const int nspin = PARAM.inp.nspin;
    const double nelec = PARAM.inp.nelec;
    double nelec_up = 0.5 * nelec;
    double nelec_down = 0.5 * nelec;
    bool fix_mag = false;
    if (nspin == 2)
    {
        // Collinear RDMFT: the exchange functional is spin-diagonal, so the
        // per-spin electron numbers are conserved independently.  Derive the two
        // targets from the converged KS reference occupations rather than from
        // inp.nupdown (which may be auto-filled internally).
        double nup = 0.0;
        double ndw = 0.0;
        const int nkloc = this->pelec->wg.nr;
        const int nbloc = this->pelec->wg.nc;
        for (int ik = 0; ik < nkloc; ++ik)
        {
            const int is = (ik < static_cast<int>(this->kv.isk.size())) ? this->kv.isk[ik] : 0;
            double s = 0.0;
            for (int ib = 0; ib < nbloc; ++ib)
            {
                s += this->pelec->wg(ik, ib);
            }
            if (is == 0)
            {
                nup += s;
            }
            else
            {
                ndw += s;
            }
        }
        nelec_up = nup;
        nelec_down = ndw;
        fix_mag = true;
    }
    const double hybrid_alpha = GlobalC::exx_info.info_global.hybrid_alpha;

    RdmftBackendPW<T, Device> backend(static_cast<hamilt::Hamilt<T, Device>*>(this->p_hamilt),
                                      this->exx_helper, psi, this->pelec, this->kv, ucell, xc,
                                      hybrid_alpha, nspin, nelec, fix_mag, nelec_up, nelec_down);
    const rdmft_core::OccConstraints& con = backend.occ_constraints();

    const int nk = this->pelec->wg.nr;
    const int nbands = this->pelec->wg.nc;
    std::vector<double> occ(con.size(), 0.0);
    for (int ik = 0; ik < nk; ++ik)
    {
        const double wkeff = con.wk[ik];
        for (int ib = 0; ib < nbands; ++ib)
        {
            occ[ik * nbands + ib] = (wkeff > 0.0) ? this->pelec->wg(ik, ib) / wkeff : 0.0;
        }
    }

    rdmft_core::RdmftParams params;
    params.xc = xc_type;
    params.power_alpha = alpha;
    params.occ_optimizer = map_occ_optimizer_pw(inp.rdmft_occ_optimizer);
    params.orb_optimizer = map_orb_optimizer_pw(inp.rdmft_orb_optimizer);
    params.outer_maxiter = inp.rdmft_outer_maxiter;
    params.occ_maxiter = std::max(0, inp.rdmft_occ_maxiter);
    params.orb_maxiter = std::max(0, inp.rdmft_orb_maxiter);
    if (params.occ_maxiter > 0 && params.orb_maxiter > 0)
    {
        params.strategy = rdmft_core::SolverStrategy::Alternating;
    }
    else if (params.orb_maxiter > 0)
    {
        params.strategy = rdmft_core::SolverStrategy::OrbOnly;
    }
    else
    {
        params.strategy = rdmft_core::SolverStrategy::OccOnly;
    }
    params.energy_tol = inp.rdmft_energy_tol;
    params.occ_grad_tol = inp.rdmft_occ_grad_tol;
    params.occ_tol = inp.rdmft_occ_tol;
    params.ls_c1 = inp.rdmft_line_search_c1;
    params.ls_c2 = inp.rdmft_line_search_c2;
    params.nelec = nelec;
    params.fix_magnetization = fix_mag;
    params.nelec_up = nelec_up;
    params.nelec_down = nelec_down;
    if (inp.rdmft_occ_init_mode == "perturbed")
    {
        params.occ_init_mode = rdmft_core::OccInitMode::Perturbed;
    }
    else if (inp.rdmft_occ_init_mode == "binary")
    {
        params.occ_init_mode = rdmft_core::OccInitMode::Binary;
    }
    else if (inp.rdmft_occ_init_mode == "uniform")
    {
        params.occ_init_mode = rdmft_core::OccInitMode::Uniform;
    }
    params.occ_init_perturb = inp.rdmft_occ_init_perturb;
    params.occ_init_nbands_top = inp.rdmft_occ_init_nbands_top;

    double etot = 0.0;
    rdmft_core::RdmftDriver driver;
    rdmft_core::DriverResult result = driver.solve(backend, params, occ, etot);

    for (int ik = 0; ik < nk; ++ik)
    {
        const double wkeff = con.wk[ik];
        for (int ib = 0; ib < nbands; ++ib)
        {
            this->pelec->wg(ik, ib) = wkeff * occ[ik * nbands + ib];
        }
    }
    this->pelec->f_en.etot = etot;
    GlobalV::ofs_running << "\n RDMFT (PW/ACE esolver) finished: E = " << etot
                         << " Ry, converged = " << (result.converged ? "T" : "F") << std::endl;
#endif
}

template class ESolver_RDMFT_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_RDMFT_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace ModuleESolver
