#include "esolver_rdmft.h"

#include "source_base/global_variable.h"
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
rdmft::OccOptimizerType map_occ_optimizer(const std::string& s)
{
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
    if (lower.find("ebi") != std::string::npos)
    {
        return rdmft::OccOptimizerType::EBI;
    }
    // sd / cg / lbfgs / spg2 all use the canonical SPG2 occupation block.
    return rdmft::OccOptimizerType::SPG2;
}

rdmft::OrbOptimizerType map_orb_optimizer(const std::string& s)
{
    std::string lower = s;
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
    if (lower == "sd")
    {
        return rdmft::OrbOptimizerType::SD;
    }
    if (lower.find("lbfgs") != std::string::npos || lower.find("bfgs") != std::string::npos)
    {
        return rdmft::OrbOptimizerType::LBFGS;
    }
    return rdmft::OrbOptimizerType::CG;
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

    const Input_para& inp = PARAM.inp;

    // --- functional / XC context ------------------------------------------
    const rdmft::XCFunctionalType xc_type
        = rdmft::parse_xc_type_or_quit(inp.rdmft_functional, "ESolver_RDMFT::after_scf");
    const rdmft::XCFunctional xc_func(xc_type, inp.rdmft_power_alpha);
    // RDMFT exchange kernels are evaluated in a HF XC context (see the legacy
    // path in ESolver_KS_LCAO::after_scf); the occupation coupling is carried
    // by xc_type inside EnergyGradient.
    XC_Functional::set_xc_type("hf");

    if (!this->eg_ready_)
    {
        this->eg_.init(&this->pv, &ucell, &this->gd, &this->kv, this->pelec, &this->orb_,
                       &this->two_center_bundle_, xc_func);
        this->eg_ready_ = true;
    }
    this->eg_.update_ion(ucell, *(this->pw_rho), this->locpp.vloc, this->sf.strucFac);

    // --- electron-number targets ------------------------------------------
    const int nspin = PARAM.inp.nspin;
    const bool fix_mag = (nspin == 2) && (std::fabs(inp.nupdown) > 1.0e-12);
    const double nelec = PARAM.inp.nelec;
    const double nelec_up = 0.5 * (nelec + inp.nupdown);
    const double nelec_down = 0.5 * (nelec - inp.nupdown);

    // --- LCAO backend around the validated energy/gradient engine ----------
    ModuleESolver::RdmftBackendLCAO<TK, TR> backend(this->eg_, *(this->psi), this->kv, nspin, nelec,
                                                    fix_mag, nelec_up, nelec_down);
    const rdmft::OccConstraints& con = backend.occ_constraints();

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
    rdmft::RdmftParams params;
    params.xc = xc_type == rdmft::XCFunctionalType::HF ? rdmft::XcType::HF : rdmft::XcType::Power;
    params.power_alpha = inp.rdmft_power_alpha;
    params.occ_optimizer = map_occ_optimizer(inp.rdmft_occ_optimizer);
    params.orb_optimizer = map_orb_optimizer(inp.rdmft_orb_optimizer);
    if (inp.rdmft_occ_init_mode == "perturbed")
    {
        params.occ_init_mode = rdmft::OccInitMode::Perturbed;
    }
    else if (inp.rdmft_occ_init_mode == "binary")
    {
        params.occ_init_mode = rdmft::OccInitMode::Binary;
    }
    else if (inp.rdmft_occ_init_mode == "uniform")
    {
        params.occ_init_mode = rdmft::OccInitMode::Uniform;
    }
    else
    {
        params.occ_init_mode = rdmft::OccInitMode::KS;
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
        params.strategy = rdmft::SolverStrategy::OrbOnly;
    }
    else if (inp.rdmft_orb_maxiter <= 0)
    {
        params.strategy = rdmft::SolverStrategy::OccOnly;
    }
    else
    {
        params.strategy = rdmft::SolverStrategy::Alternating;
    }

    // --- solve ------------------------------------------------------------
    double etot = 0.0;
    rdmft::RdmftDriver driver;
    rdmft::DriverResult result = driver.solve(backend, params, occ, etot);
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
#endif
}

template class ESolver_RDMFT<double, double>;
template class ESolver_RDMFT<std::complex<double>, double>;
template class ESolver_RDMFT<std::complex<double>, std::complex<double>>;

} // namespace ModuleESolver
