#include "esolver_rdmft_pw.h"

#include "source_base/global_variable.h"
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
        const double spin_deg = (nspin == 1) ? 2.0 : 1.0;
        spin_deg_ = spin_deg;

        con_.nbnd = nbands_;
        con_.nks = nks_;
        con_.wk.resize(nks_);
        con_.isk.resize(nks_);
        for (int ik = 0; ik < nks_; ++ik)
        {
            con_.wk[ik] = spin_deg * kv.wk[ik];
            con_.isk[ik] = (nspin == 2 && ik >= nks_ / 2) ? 2 : 1;
        }
        con_.n_target = nelec;
        con_.fix_magnetization = fix_mag;
        con_.n_target_up = nelec_up;
        con_.n_target_down = nelec_down;

        wg_.create(nks_, nbands_);
        wg_x_.create(nks_, nbands_);
        diag_.assign(nks_ * nbands_, 0.0);
        vx_diag_.assign(nks_ * nbands_, 0.0);
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
        const double e_x = exchange_energy(occ);
        return e_bandlike - e_h + e_x + e_ewald;
    }

    void grad_occ(const std::vector<double>& occ, std::vector<double>& grad) override
    {
        refresh_state(occ);
        exchange_diag(occ); // fills vx_diag_ (already hybrid_alpha-scaled)
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

  private:
    // Sync density/potential and recompute one-body diagonals for occ.
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
        compute_onebody_diag();
    }

    // diag_[ik*nbands+ib] = Re <psi|H_KS(no EXX)|psi>.
    void compute_onebody_diag()
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
            psi::Psi<T, Device> wrapper(psi_k, 1, nbands_, nbasis, npw);
            psi::Range bands(true, 0, 0, nbands_ - 1);
            typename hamilt::Operator<T, Device>::hpsi_info info(&wrapper, bands, hpsi.data());
            p_hamilt_->ops->hPsi(info);
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const T* pc = psi_k + static_cast<size_t>(ib) * nbasis;
                const T* hc = hpsi.data() + static_cast<size_t>(ib) * nbasis;
                double d = 0.0;
                for (int ig = 0; ig < npw; ++ig)
                {
                    d += (std::conj(pc[ig]) * hc[ig]).real();
                }
                diag_[ik * nbands_ + ib] = d;
            }
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(diag_.data(), static_cast<int>(diag_.size()));
#endif
        exx_helper_->set_op_first_iter(prev_first);
    }

    // Build ACE exchange with weights w_t(n) and return the exchange energy.
    double exchange_energy(const std::vector<double>& occ)
    {
        set_exchange_weights(occ);
        exx_helper_->set_wg(&wg_x_);
        exx_helper_->set_psi(psi_); // rebuilds ACE projectors (exxace && separate_loop)
        exx_helper_->set_op_first_iter(false);
        return hybrid_alpha_ * exx_helper_->cal_exx_energy(psi_);
    }

    // vx_diag_[i] = hybrid_alpha * <psi|Vx[w_t(n)]|psi> summed over channels.
    void exchange_diag(const std::vector<double>& occ)
    {
        set_exchange_weights(occ);
        exx_helper_->set_wg(&wg_x_);
        exx_helper_->set_psi(psi_);
        exx_helper_->set_op_first_iter(false);
        Exx_Helper<T, Device>* helper = static_cast<Exx_Helper<T, Device>*>(exx_helper_);
        hamilt::OperatorEXXPW<T, Device>* op = helper->op_exx;
        std::fill(vx_diag_.begin(), vx_diag_.end(), 0.0);
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
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const T* pc = psi_k + static_cast<size_t>(ib) * nbasis;
                const T* vc = vxpsi.data() + static_cast<size_t>(ib) * nbasis;
                double d = 0.0;
                for (int ig = 0; ig < npw; ++ig)
                {
                    d += (std::conj(pc[ig]) * vc[ig]).real();
                }
                vx_diag_[ik * nbands_ + ib] = d;
            }
        }
#ifdef __MPI
        Parallel_Reduce::reduce_pool(vx_diag_.data(), static_cast<int>(vx_diag_.size()));
#endif
        // op->act already applies hybrid_alpha internally.
    }

    void set_exchange_weights(const std::vector<double>& occ)
    {
        const int nch = xc_.n_channels();
        for (int ik = 0; ik < nks_; ++ik)
        {
            for (int ib = 0; ib < nbands_; ++ib)
            {
                double w = 0.0;
                for (int it = 1; it <= nch; ++it)
                {
                    w = std::max(w, xc_.channel(it, occ[ik * nbands_ + ib]).w);
                }
                wg_x_(ik, ib) = con_.wk[ik] * w;
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
    double spin_deg_ = 1.0;
    int nbands_ = 0;
    int nks_ = 0;
    rdmft_core::OccConstraints con_;
    ModuleBase::matrix wg_;
    ModuleBase::matrix wg_x_;
    std::vector<double> diag_;
    std::vector<double> vx_diag_;
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

    rdmft_core::XcType xc_type = rdmft_core::XcType::HF;
    double alpha = inp.rdmft_power_alpha;
    rdmft_core::parse_xc_type(inp.rdmft_functional, xc_type, alpha);
    const rdmft_core::RdmftXC xc(xc_type, alpha, 1.0e-8);

    const int nspin = PARAM.inp.nspin;
    const bool fix_mag = (nspin == 2) && (std::fabs(inp.nupdown) > 1.0e-12);
    const double nelec = PARAM.inp.nelec;
    const double nelec_up = 0.5 * (nelec + inp.nupdown);
    const double nelec_down = 0.5 * (nelec - inp.nupdown);
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
    params.strategy = rdmft_core::SolverStrategy::OccOnly; // orbitals frozen at hybrid-KS
    params.outer_maxiter = inp.rdmft_outer_maxiter;
    params.occ_maxiter = std::max(1, inp.rdmft_occ_maxiter);
    params.orb_maxiter = 0;
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
