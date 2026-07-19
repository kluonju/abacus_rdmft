#include "esolver_rdmft.h"

#include "source_base/global_variable.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"

#include <complex>

#ifdef __RDMFT
#include "rdmft_backend_lcao.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_lcao/module_rdmft/rdmft_input_parse.h"
#include "source_lcao/module_rdmft/rdmft_xc_functional.h"
#include "source_rdmft/rdmft_driver.h"
#include "source_rdmft/rdmft_params.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>
#endif

namespace ModuleESolver
{

template <typename TK, typename TR>
ESolver_RDMFT<TK, TR>::ESolver_RDMFT()
{
    this->classname = "ESolver_RDMFT";
}

template <typename TK, typename TR>
ESolver_RDMFT<TK, TR>::~ESolver_RDMFT()
{
}

#ifdef __RDMFT
namespace
{
rdmft_core::OccOptimizerType map_occ_optimizer(const std::string& s)
{
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
    if (lower.find("ebi") != std::string::npos)
    {
        return rdmft_core::OccOptimizerType::EBI;
    }
    // sd / cg / lbfgs / spg2 all use the canonical SPG2 occupation block.
    return rdmft_core::OccOptimizerType::SPG2;
}

rdmft_core::OrbOptimizerType map_orb_optimizer(const std::string& s)
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

template <typename TK, typename TR>
void ESolver_RDMFT<TK, TR>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    // Finish the Kohn-Sham post-processing first (density, outputs, ...).
    ESolver_KS_LCAO<TK, TR>::after_scf(ucell, istep, conv_esolver);

#ifdef __RDMFT
    if (this->psi == nullptr)
    {
        return;
    }

    if (PARAM.inp.nspin == 4)
    {
        ModuleBase::WARNING_QUIT("ESolver_RDMFT",
                                 "non-collinear (nspin=4) RDMFT is under development.");
    }

    const Input_para& inp = PARAM.inp;

    // --- functional / XC context ------------------------------------------
    const rdmft::XCFunctionalType xc_type
        = rdmft::parse_xc_type_or_quit(inp.rdmft_functional, "ESolver_RDMFT::after_scf");
    const rdmft::XCFunctional xc_func(xc_type, inp.rdmft_power_alpha);
    // RDMFT exchange kernels are evaluated in a HF XC context (see the legacy
    // path in ESolver_KS_LCAO::after_scf); the occupation coupling is carried
    // by xc_type inside EnergyGradient.
    XC_Functional::set_xc_type("hf");

    // EnergyGradient enables its own exact-exchange infrastructure only when
    // PARAM.inp.rdmft is set (rdmft_energy_gradient.cpp).  This esolver is
    // selected through esolver_type=rdmft rather than that flag, so enable it
    // around the RDMFT stage and restore it afterwards (the base KS after_scf
    // above has already run with the original flag, so its legacy embedded
    // RDMFT path is not triggered).
    const bool rdmft_flag_saved = PARAM.inp.rdmft;
    const_cast<Input_para&>(PARAM.inp).rdmft = true;

    if (!this->eg_ready_)
    {
        this->eg_.init(&this->pv, &ucell, &this->gd, &this->kv, this->pelec, &this->orb_,
                       &this->two_center_bundle_, xc_func);
        this->eg_ready_ = true;
    }
    this->eg_.update_ion(ucell, *(this->pw_rho), this->locpp.vloc, this->sf.strucFac);

    // --- electron-number targets ------------------------------------------
    const int nspin = PARAM.inp.nspin;
    const double nelec = PARAM.inp.nelec;
    double nelec_up = 0.5 * nelec;
    double nelec_down = 0.5 * nelec;
    bool fix_mag = false;
    if (nspin == 2)
    {
        // Collinear RDMFT: the exact-exchange (and power) functional is
        // spin-diagonal, so the per-spin electron numbers are conserved
        // independently during the occupation optimisation.  Derive the two
        // targets from the converged KS reference occupations rather than from
        // inp.nupdown, which may be resolved internally in a way that does not
        // reflect the actual spin populations (e.g. auto-filled defaults).
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

    // The KS solve leaves psi's internal cursor (current_k / psi_bias) pointing
    // at the last processed spin/k block.  The RDMFT engine and backend index
    // k-blocks explicitly from the buffer start (psi(ik,0,0)) but the backend's
    // flat save/restore helpers use get_pointer()+size(); a non-zero psi_bias
    // would make those read one block past the buffer (harmless for nspin=1
    // where current_k stays 0, but an out-of-bounds access for nspin>=2).  Reset
    // the cursor so get_pointer() == buffer start before any backend copies.
    this->psi->fix_k(0);

    // --- LCAO backend around the validated energy/gradient engine ----------
    ModuleESolver::RdmftBackendLCAO<TK, TR> backend(this->eg_, *(this->psi), this->kv, nspin, nelec,
                                                    fix_mag, nelec_up, nelec_down);
    const rdmft_core::OccConstraints& con = backend.occ_constraints();

    // --- initial occupations n_ik = wg / (spin_deg * wk) -------------------
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

    // --- configuration ----------------------------------------------------
    rdmft_core::RdmftParams params;
    params.xc = xc_type == rdmft::XCFunctionalType::HF ? rdmft_core::XcType::HF : rdmft_core::XcType::Power;
    params.power_alpha = inp.rdmft_power_alpha;
    params.occ_optimizer = map_occ_optimizer(inp.rdmft_occ_optimizer);
    params.orb_optimizer = map_orb_optimizer(inp.rdmft_orb_optimizer);
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
    else
    {
        params.occ_init_mode = rdmft_core::OccInitMode::KS;
    }
    params.occ_init_perturb = inp.rdmft_occ_init_perturb;
    params.occ_init_nbands_top = inp.rdmft_occ_init_nbands_top;
    params.outer_maxiter = inp.rdmft_outer_maxiter;
    params.occ_maxiter = std::max(0, inp.rdmft_occ_maxiter);
    params.orb_maxiter = std::max(0, inp.rdmft_orb_maxiter);
    params.energy_tol = inp.rdmft_energy_tol;
    params.occ_grad_tol = inp.rdmft_occ_grad_tol;
    params.orb_grad_tol = inp.rdmft_orb_grad_tol;
    params.occ_tol = inp.rdmft_occ_tol;
    params.ls_c1 = inp.rdmft_line_search_c1;
    params.ls_c2 = inp.rdmft_line_search_c2;
    params.ls_max_zoom = std::max(1, inp.rdmft_line_search_max_zoom);
    params.nelec = nelec;
    params.fix_magnetization = fix_mag;
    params.nelec_up = nelec_up;
    params.nelec_down = nelec_down;
    if (inp.rdmft_occ_maxiter <= 0)
    {
        params.strategy = rdmft_core::SolverStrategy::OrbOnly;
    }
    else if (inp.rdmft_orb_maxiter <= 0)
    {
        params.strategy = rdmft_core::SolverStrategy::OccOnly;
    }
    else
    {
        params.strategy = rdmft_core::SolverStrategy::Alternating;
    }

    // --- solve ------------------------------------------------------------
    double etot = 0.0;
    rdmft_core::RdmftDriver driver;
    rdmft_core::DriverResult result = driver.solve(backend, params, occ, etot);

    // Diagnostic: report the optimised natural occupations at the first k-point
    // and the Stiefel gram residual ||X^H X - I|| (X-space; the generalised
    // C^H S C = I constraint is enforced through the Cholesky S transform).
    {
        std::vector<double> gram;
        this->eg_.stiefel_gram_residual_frobenius_per_k(*(this->psi), gram);
        double gmax = 0.0;
        for (double g : gram) { gmax = std::max(gmax, g); }
        GlobalV::ofs_running << " RDMFT(LCAO) natural occ (ik=0):";
        for (int ib = 0; ib < std::min(nbands, 8); ++ib)
        {
            GlobalV::ofs_running << " " << occ[ib];
        }
        GlobalV::ofs_running << " | max||X^H X - I||_F=" << gmax << std::endl;
    }

    backend.finalize_orbitals(); // X-space -> C-space natural orbitals

    // --- write back optimised occupations and energy ----------------------
    for (int ik = 0; ik < nk; ++ik)
    {
        const double wkeff = con.wk[ik];
        for (int ib = 0; ib < nbands; ++ib)
        {
            this->pelec->wg(ik, ib) = wkeff * occ[ik * nbands + ib];
        }
    }
    this->pelec->f_en.etot = etot;

    GlobalV::ofs_running << "\n RDMFT (modular esolver) finished: E = " << etot
                         << " Ry, converged = " << (result.converged ? "T" : "F")
                         << ", outer iterations = " << result.outer_iterations << std::endl;

    const_cast<Input_para&>(PARAM.inp).rdmft = rdmft_flag_saved;
#endif
}

template class ESolver_RDMFT<double, double>;
template class ESolver_RDMFT<std::complex<double>, double>;
template class ESolver_RDMFT<std::complex<double>, std::complex<double>>;

} // namespace ModuleESolver
