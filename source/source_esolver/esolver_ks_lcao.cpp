#include "esolver_ks_lcao.h"
#include "source_base/module_external/blacs_connector.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_estate/elecstate_tools.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_lcao/module_deltaspin/deltaspin_lcao.h"
#include "source_lcao/dftu_lcao.h"
#include "source_lcao/hs_matrix_k.hpp" // there may be multiple definitions if using hpp
#include "source_estate/module_charge/symmetry_rho.h"
#include "source_lcao/LCAO_domain.h" // need DeePKS_init
#include "source_lcao/FORCE_STRESS.h"
#include "source_lcao/module_gint/gint.h"
#include "source_estate/elecstate_lcao.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_hsolver/hsolver_lcao.h"
#ifdef __EXX
#include "../source_lcao/module_ri/exx_opt_orb.h"
#endif
#ifdef __RDMFT
#include "source_lcao/module_rdmft/rdmft.h"
#include "source_lcao/module_rdmft/rdmft_input_parse.h"
#include "source_lcao/module_rdmft/rdmft_solver.h"
#endif
#include "source_estate/module_charge/chgmixing.h" // use charge mixing, mohan add 20251006
#include "source_estate/module_dm/init_dm.h" // init dm from electronic wave functions
#include "source_io/module_ctrl/ctrl_runner_lcao.h" // use ctrl_runner_lcao() 
#include "source_io/module_ctrl/ctrl_iter_lcao.h" // use ctrl_iter_lcao() 
#include "source_io/module_ctrl/ctrl_scf_lcao.h" // use ctrl_scf_lcao()
#include "source_io/module_output/print_info.h"
#include "source_lcao/rho_tau_lcao.h" // mohan add 20251024
#include "source_lcao/LCAO_set.h" // mohan add 20251111
#include "source_psi/setup_psi.h" // use Setup_Psi for deallocate_psi
namespace ModuleESolver
{

template <typename TK, typename TR>
ESolver_KS_LCAO<TK, TR>::ESolver_KS_LCAO()
{
    this->classname = "ESolver_KS_LCAO";
    this->basisname = "LCAO";
}

template <typename TK, typename TR>
ESolver_KS_LCAO<TK, TR>::~ESolver_KS_LCAO()
{
    //****************************************************
    // do not add any codes in this deconstructor funcion
    //****************************************************
    Setup_Psi<TK>::deallocate_psi(this->psi);
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_all_runners");
    ModuleBase::timer::start("ESolver_KS_LCAO", "before_all_runners");

    // 0) init EXX - moved from constructor to ensure GlobalC::exx_info.info_global is already set
    this->exx_nao.init();

    // 1) before_all_runners in ESolver_KS
    ESolver_KS::before_all_runners(ucell, inp);

    // 2) autoset nbands in ElecState before init_basis (for Psi 2d division)
    if (this->pelec == nullptr)
    {
        // TK stands for double and std::complex<double>?
        this->pelec = new elecstate::ElecStateLCAO<TK>(&(this->chr), &(this->kv),
          this->kv.get_nks(), this->pw_big);
    }

    // 3) read LCAO orbitals/projectors and construct the interpolation tables.
    LCAO_domain::init_basis_lcao(this->pv, inp.onsite_radius, inp.lcao_ecut,
      inp.lcao_dk, inp.lcao_dr, inp.lcao_rmax, ucell, two_center_bundle_, orb_);

    // 4) setup EXX calculations
    if (inp.calculation == "gen_opt_abfs")
    {
#ifdef __EXX
        Exx_Opt_Orb exx_opt_orb;
        exx_opt_orb.generate_matrix(GlobalC::exx_info.info_opt_abfs, this->kv, ucell, this->orb_);
#else
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::before_all_runners", "calculation=gen_opt_abfs must compile __EXX");
#endif
        return;
    }

    LCAO_domain::set_psi_occ_dm_chg<TK>(this->kv, this->psi, this->pv, this->pelec,
      this->dmat, this->chr, inp);

    LCAO_domain::set_pot<TK>(ucell, this->kv, this->sf, *this->pw_rho, *this->pw_rhod,
      this->pelec, this->orb_, this->pv, this->locpp, this->dftu,
      this->solvent, this->exx_nao, this->deepks, inp);

    //! if kpar is not divisible by nks, print a warning
    ModuleIO::print_kpar(this->kv.get_nks(), PARAM.globalv.kpar_lcao);

    ModuleBase::timer::end("ESolver_KS_LCAO", "before_all_runners");
    return;
}


template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "before_scf");
    ModuleBase::timer::start("ESolver_KS_LCAO", "before_scf");

    //! 1) call before_scf() of ESolver_KS.
    ESolver_KS::before_scf(ucell, istep);

    //! 2) find search radius
    double search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
      PARAM.inp.out_level, orb_.get_rcutmax_Phi(), ucell.infoNL.get_rcutmax_Beta(),
      PARAM.globalv.gamma_only_local);

    //! 3) use search_radius to search adj atoms
    atom_arrange::search(PARAM.globalv.search_pbc, GlobalV::ofs_running,
      this->gd, ucell, search_radius, PARAM.inp.test_atom_input);

    //! 4) initialize NAO basis set
    // here new is a unique pointer, which will be deleted automatically
    gint_info_.reset(
        new ModuleGint::GintInfo(
        this->pw_big->nbx, this->pw_big->nby, this->pw_big->nbz,
        this->pw_rho->nx, this->pw_rho->ny, this->pw_rho->nz,
        0, 0, this->pw_big->nbzp_start,
        this->pw_big->nbx, this->pw_big->nby, this->pw_big->nbzp,
        orb_.Phi, ucell, this->gd));
    ModuleGint::Gint::set_gint_info(gint_info_.get());

    // 7) For each atom, calculate the adjacent atoms in different cells
    // and allocate the space for H(R) and S(R).
    // If k point is used here, allocate HlocR after atom_arrange.
    this->RA.for_2d(ucell, this->gd, this->pv, PARAM.globalv.gamma_only_local, orb_.cutoffs());

    // 8) initialize the Hamiltonian operators
    // if atom moves, then delete old pointer and add a new one
    if (this->p_hamilt != nullptr)
    {
        delete this->p_hamilt;
        this->p_hamilt = nullptr;
    }
    if (this->p_hamilt == nullptr)
    {
        this->p_hamilt = new hamilt::HamiltLCAO<TK, TR>(
            ucell, this->gd, &this->pv, this->pelec->pot, this->kv,
            two_center_bundle_, orb_, this->dmat.dm, &this->dftu, this->deepks, istep, exx_nao);
    }

    // 9) for each ionic step, the overlap <phi|alpha> must be rebuilt
    // since it depends on ionic positions.
    // overlap_orb_alpha is only built when DeePKS is enabled (descriptor
    // orbitals); guard the dereference so non-DeePKS runs don't form a
    // reference from a null unique_ptr (undefined behaviour).
    if (two_center_bundle_.overlap_orb_alpha)
    {
        this->deepks.build_overlap(ucell, orb_, pv, gd, *(two_center_bundle_.overlap_orb_alpha), PARAM.inp);
    }

    // 10) prepare sc calculation
    init_deltaspin_lcao<TK>(ucell, PARAM.inp, &(this->pv), this->kv, this->p_hamilt, this->psi, this->dmat.dm, this->pelec);

    // 11) set xc type before the first cal of xc in pelec->init_scf, Peize Lin add 2016-12-03
    this->exx_nao.before_scf(ucell, this->kv, orb_, this->p_chgmix, istep, PARAM.inp);

    // 12) initalize DM(R), which has the same size with Hamiltonian(R)
    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);

    if(!hamilt_lcao)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::before_scf","p_hamilt does not exist");
    }
    this->dmat.dm->init_DMR(*hamilt_lcao->getHR());

    // 13.1) decide the strategy for initializing DMR and HR
    if(istep == 0)//if the first scf step, readin DMR from file,
    {
        //calculate or readin the density matrix DMR
        if(PARAM.inp.init_chg == "dm" || PARAM.inp.init_chg == "dm_no_renormalize")
        {
            //! 13.1.1) init charge density from density matrix file
            LCAO_domain::init_chg_dm<TK>(PARAM.globalv.global_readin_dir, PARAM.inp.nspin,
                this->dmat, ucell, &(this->pv), this->pelec->charge);
        }
        if(PARAM.inp.init_chg == "hr")
        {
            //! 13.1.2) init charge density from Hamiltonian matrix file
            LCAO_domain::init_chg_hr<TK, TR>(PARAM.globalv.global_readin_dir, PARAM.inp.nspin,
                static_cast<hamilt::Hamilt<TK>*>(this->p_hamilt), ucell, &(this->pv), this->psi[0], this->pelec, *this->dmat.dm,
                this->chr, PARAM.inp.ks_solver);
        }
    }
    else if(PARAM.inp.esolver_type!="tddft")//if not, use the DMR calculated from last step
    {
        // 13.1.2) two cases are considered:
        // 1. DMK in DensityMatrix is not empty (istep > 0), then DMR is initialized by DMK
        // 2. DMK in DensityMatrix is empty (istep == 0), then DMR is initialized by zeros
        this->dmat.dm->cal_DMR();
    }
    // 13.2) init_scf, should be before_scf? mohan add 2025-03-10
    elecstate::init_scf(ucell, this->Pgrid, this->sf.strucFac, this->locpp.numeric,
                          istep, PARAM.globalv.global_out_dir, PARAM.inp, this->pelec);

#ifdef __MLALGO
    // 14) initialize DM2(R) of DeePKS, the DM2(R) is different from DM(R)
    this->deepks.ld.init_DMR(ucell, orb_, this->pv, this->gd);
#endif

    // 16) the electron charge density should be symmetrized,
    Symmetry_rho::symmetrize_rho(PARAM.inp.nspin, this->chr, this->pw_rho, ucell.symm);

#ifdef __RDMFT
    // 17) update of RDMFT (only after first init in after_scf, post-KS)
    if (PARAM.inp.rdmft == true && this->rdmft_module_initialized)
    {
        const bool use_new_rdmft_engine = !PARAM.inp.rdmft_functional.empty();
        if (use_new_rdmft_engine)
        {
            if (rdmft_eg_initialized)
            {
                rdmft_eg.update_ion(ucell, *(this->pw_rho), this->locpp.vloc, this->sf.strucFac);
            }
        }
        else
        {
            rdmft_solver.update_ion(ucell, *(this->pw_rho), this->locpp.vloc, this->sf.strucFac);
        }
    }
#endif

    ModuleBase::timer::end("ESolver_KS_LCAO", "before_scf");
    return;
}


template <typename TK, typename TR>
double ESolver_KS_LCAO<TK, TR>::cal_energy()
{
    return this->pelec->f_en.etot;
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_force(UnitCell& ucell, ModuleBase::matrix& force)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_force");
    ModuleBase::timer::start("ESolver_KS_LCAO", "cal_force");

    Force_Stress_LCAO<TK> fsl(this->RA, ucell.nat);

    deepks.dpks_out_type = "tot";  // for deepks method

    fsl.getForceStress(ucell, PARAM.inp.cal_force, PARAM.inp.cal_stress, 
                       PARAM.inp.test_force, PARAM.inp.test_stress,
                       this->gd, this->pv, this->pelec, this->dmat, this->psi,
                       two_center_bundle_, orb_, force, this->scs,
                       this->locpp, this->sf, this->kv,
                       this->pw_rho, this->solvent, this->dftu, this->deepks,
                       this->exx_nao, &ucell.symm, PARAM.inp.td_stype,
                       static_cast<hamilt::Hamilt<TK>*>(this->p_hamilt));

    // delete RA after cal_force
    this->RA.delete_grid();

    this->have_force = true;

    ModuleBase::timer::end("ESolver_KS_LCAO", "cal_force");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::cal_stress(UnitCell& ucell, ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "cal_stress");
    ModuleBase::timer::start("ESolver_KS_LCAO", "cal_stress");

    if (!this->have_force)
    {
        ModuleBase::matrix fcs;
        this->cal_force(ucell, fcs);
    }

    // the stress has been calculated in 'cal_force'
    stress = this->scs;
    this->have_force = false;

    ModuleBase::timer::end("ESolver_KS_LCAO", "cal_stress");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_all_runners(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_all_runners");
    ModuleBase::timer::start("ESolver_KS_LCAO", "after_all_runners");

    ESolver_KS::after_all_runners(ucell);

    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);
    if(!hamilt_lcao)
    {
	    ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_all_runners","p_hamilt does not exist");
    }

    ModuleIO::ctrl_runner_lcao<TK, TR>(ucell,
		    PARAM.inp, this->kv, this->pelec, this->dmat, this->pv, this->Pgrid, 
		    this->gd, this->psi, this->chr, hamilt_lcao,
		    this->two_center_bundle_,
		    this->orb_, this->pw_rho, this->pw_rhod,
		    this->sf, this->locpp.vloc, this->exx_nao, this->solvent);


#ifdef __MPI
#ifdef __LCAO
    // Exit BLACS environment for LCAO calculations
    Cblacs_exit(1);
#endif
#endif

    ModuleBase::timer::end("ESolver_KS_LCAO", "after_all_runners");
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_init");

    // call iter_init() of ESolver_KS
    ESolver_KS::iter_init(ucell, istep, iter);

    module_charge::chgmixing_ks_lcao(iter, this->p_chgmix, this->dftu, 
      this->dmat.dm->get_DMR_pointer(1)->get_nnr(), PARAM.inp); 

    if (iter == 1)
    {
        this->gint_precision_controller_.set_mode(PARAM.inp.gint_precision);
        this->gint_precision_controller_.reset_for_new_scf();
        this->gint_info_->set_exec_precision(this->gint_precision_controller_.current_precision());
        if (PARAM.inp.gint_precision == "mix")
        {
            GlobalV::ofs_running << "\n >> Gint mixed-precision mode: starting SCF with fp32"
                                 << " (will switch to fp64 when drho is small enough)" << std::endl;
            std::cout << " >> NOTICE: Gint grid-integration starts with fp32 (mixed-precision mode)" << std::endl;
        }
        else if (PARAM.inp.gint_precision == "single")
        {
            GlobalV::ofs_running << "\n >> Gint single-precision mode: using fp32 throughout SCF" << std::endl;
            std::cout << " >> NOTICE: Gint grid-integration uses fp32 throughout SCF (single-precision mode)" << std::endl;
        }
    }

    // mohan update 2012-06-05
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband(ucell);

    if (istep == 0 && PARAM.inp.init_wfc == "file")
	{
		int exx_two_level_step = 0;
#ifdef __EXX
		if (GlobalC::exx_info.info_global.cal_exx)
		{
			// the following steps are only needed in the first outer exx loop
			exx_two_level_step
				= GlobalC::exx_info.info_ri.real_number ? 
                  this->exx_nao.exd->two_level_step : this->exx_nao.exc->two_level_step;
		}
#endif
		elecstate::init_dm<TK>(ucell, this->pelec, this->dmat, this->psi, this->chr, iter, exx_two_level_step);
	}

#ifdef __EXX
    // calculate exact-exchange
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exx_nao.exd->exx_eachiterinit(istep, ucell, *this->dmat.dm, this->kv, iter);
        }
        else
        {
            this->exx_nao.exc->exx_eachiterinit(istep, ucell, *this->dmat.dm, this->kv, iter);
        }
    }
#endif

    init_dftu_lcao<TK>(istep, iter, PARAM.inp, &(this->dftu), this->dmat.dm, ucell, this->chr.rho, this->pw_rho->nrxx);

#ifdef __MLALGO
    // the density matrixes of DeePKS have been updated in each iter
    this->deepks.ld.set_hr_cal(true);

    // HR in HamiltLCAO should be recalculate
    if (PARAM.inp.deepks_scf)
    {
        this->p_hamilt->refresh();
    }
#endif

    if (PARAM.inp.vl_in_h)
    {
        // update real space Hamiltonian
        this->p_hamilt->refresh();
    }

    // save density matrix DMR for mixing
    if (PARAM.inp.mixing_restart > 0 && PARAM.inp.mixing_dmr && this->p_chgmix->mixing_restart_count > 0)
    {
        this->dmat.dm->save_DMR();
    }
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::hamilt2rho_single(UnitCell& ucell, int istep, int iter, double ethr)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "hamilt2rho_single");

    // 1) reset energy
    this->pelec->f_en.eband = 0.0;
    this->pelec->f_en.demet = 0.0;
    bool skip_charge = PARAM.inp.calculation == "nscf" ? true : false;

    // 2) run the inner lambda loop to contrain atomic moments with the DeltaSpin method
    bool skip_solve = false;
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
        if (PARAM.inp.sc_lambda_strategy == "linear_scan")
        {
            sc.run_lambda_linear_scan(iter - 1);
            skip_solve = true;
        }
        else if (!sc.mag_converged() && this->drho > 0 && this->drho < PARAM.inp.sc_scf_thr)
        {
            sc.run_lambda_loop(iter - 1);
            sc.set_mag_converged(true);
            skip_solve = true;
        }
        else if (sc.mag_converged())
        {
            sc.run_lambda_loop(iter - 1);
            skip_solve = true;
        }
    }

    // 3) run Hsolver
    if (!skip_solve)
    {
        hsolver::HSolverLCAO<TK> hsolver_lcao_obj(&(this->pv), PARAM.inp.ks_solver);
        hsolver_lcao_obj.solve(static_cast<hamilt::Hamilt<TK>*>(this->p_hamilt), this->psi[0], this->pelec, *this->dmat.dm, 
          this->chr, PARAM.inp.nspin, skip_charge);
    }
    else
    {
        // Lambda loop updated the density matrix (DM) but not the real-space charge density.
        // HSolver was skipped, so we need to sync rho from DM manually.
        LCAO_domain::dm2rho(this->dmat.dm->get_DMR_vector(), PARAM.inp.nspin, &this->chr);
    }

    // 4) EXX
#ifdef __EXX
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_ri.real_number)
        {
            this->exx_nao.exd->exx_hamilt2rho(*this->pelec, this->pv, iter);
        }
        else
        {
            this->exx_nao.exc->exx_hamilt2rho(*this->pelec, this->pv, iter);
        }
    }
#endif

    // 5) symmetrize the charge density
    Symmetry_rho::symmetrize_rho(PARAM.inp.nspin, this->chr, this->pw_rho, ucell.symm);

    // 6) calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);
}


template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "iter_finish");

    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);

    if(!hamilt_lcao)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::iter_finish","p_hamilt does not exist");
    }

	const std::vector<std::vector<TK>>& dm_vec = this->dmat.dm->get_DMK_vector();

    // 1) calculate the local occupation number matrix and energy correction in DFT+U
    finish_dftu_lcao<TK>(iter, conv_esolver, PARAM.inp, &(this->dftu), ucell, dm_vec, this->kv, this->p_chgmix->get_mixing_beta(), hamilt_lcao);

    // 2) for deepks, calculate delta_e, output labels during electronic steps
    this->deepks.delta_e(ucell, this->kv, this->orb_, this->pv, this->gd, dm_vec, this->pelec->f_en, PARAM.inp);

    // 3) for delta spin
    cal_mi_lcao_wrapper<TK>(iter, PARAM.inp);

    // call iter_finish() of ESolver_KS, where band gap is printed,
    // eig and occ are printed, magnetization is calculated,
    // charge mixing is performed, potential is updated, 
    // HF and kS energies are computed, meta-GGA, Jason and restart
    ESolver_KS::iter_finish(ucell, istep, iter, conv_esolver);
    const bool precision_switched = this->gint_precision_controller_.update_after_iteration(this->drho, this->scf_thr);
    this->gint_info_->set_exec_precision(this->gint_precision_controller_.current_precision());
    if (precision_switched)
    {
        GlobalV::ofs_running << "\n >> Gint precision switched: fp32 -> fp64 (drho = "
                             << this->drho << ")" << std::endl;
        std::cout << " >> NOTICE: Gint grid-integration precision switched from fp32 to fp64" << std::endl;
    }

    // mix density matrix if mixing_restart + mixing_dmr + not first
    // mixing_restart at every iter except the last iter
    if(iter != PARAM.inp.scf_nmax && !conv_esolver)
    {
        if (PARAM.inp.mixing_restart > 0 && this->p_chgmix->mixing_restart_count > 0 && PARAM.inp.mixing_dmr)
        {
            this->p_chgmix->mix_dmr(this->dmat.dm);
        }
    }

    // control the output related to the finished iteration
    ModuleIO::ctrl_iter_lcao<TK, TR>(ucell, PARAM.inp, this->kv, this->pelec, *this->dmat.dm,
      this->pv, this->gd, this->psi, this->chr, this->p_chgmix, 
      hamilt_lcao, this->orb_, this->deepks, 
      this->exx_nao, iter, istep, conv_esolver, this->scf_ene_thr);
}

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_scf");
    ModuleBase::timer::start("ESolver_KS_LCAO", "after_scf");

    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);

    if(!hamilt_lcao)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_scf","p_hamilt does not exist");
    }

    if (PARAM.inp.out_elf[0] > 0)
	{
		LCAO_domain::dm2tau(this->dmat.dm->get_DMR_vector(), PARAM.inp.nspin, this->pelec->charge);
	}

    //! 1) call after_scf() of ESolver_KS
    ESolver_KS::after_scf(ucell, istep, conv_esolver);

#ifdef __RDMFT
    //! Lazy-init RDMFT after KS has finished for this SCF (not in before_all_runners).
    if (PARAM.inp.rdmft == true && !this->rdmft_module_initialized)
    {
        const Input_para& inp_rdmft = PARAM.inp;
        const bool use_new_rdmft_engine = !inp_rdmft.rdmft_functional.empty();

        if (use_new_rdmft_engine)
        {
            const rdmft::XCFunctionalType xc_type
                = rdmft::parse_xc_type_or_quit(inp_rdmft.rdmft_functional, "ESolver_KS_LCAO::after_scf");
            const rdmft::XCFunctional xc_func(xc_type, inp_rdmft.rdmft_power_alpha);
            // Ensure global XC context is already HF before EnergyGradient::init,
            // because EXX backends snapshot XC-dependent settings at init time.
            XC_Functional::set_xc_type("hf");
            this->rdmft_eg.init(&this->pv, &ucell, &this->gd, &this->kv, this->pelec, &this->orb_,
                                &two_center_bundle_, xc_func);
            this->rdmft_eg_initialized = true;

            // Build ion-dependent one-body terms and set EnergyGradient::ion_initialized_.
            this->rdmft_eg.update_ion(ucell, *(this->pw_rho), this->locpp.vloc, this->sf.strucFac);
        }
        else
        {
            this->rdmft_solver.init(this->pv, ucell, this->gd, this->kv, *(this->pelec), this->orb_,
                                    two_center_bundle_, inp_rdmft.dft_functional, inp_rdmft.rdmft_power_alpha);
            this->rdmft_solver.update_ion(ucell, *(this->pw_rho), this->locpp.vloc, this->sf.strucFac);
        }

        this->rdmft_module_initialized = true;
    }
#endif

    //! 1.5) Run RDMFT optimization when the new engine is active
#ifdef __RDMFT
    if (rdmft_eg_initialized && this->psi != nullptr)
    {
        ModuleBase::timer::start("ESolver_KS_LCAO", "rdmft_solve");

        const Input_para& inp = PARAM.inp;
        const rdmft::XCFunctionalType rdmft_xc_type
            = rdmft::parse_xc_type_or_quit(inp.rdmft_functional, "ESolver_KS_LCAO::after_scf");
        const std::string xc_func_restore = inp.dft_functional;
        bool rdmft_xc_context_switched = false;
        if (!inp.rdmft_functional.empty())
        {
            // Decouple the RDMFT objective/gradient path from KS `dft_functional`.
            // RDMFT exchange kernels are evaluated in a HF XC context, while the
            // occupation coupling (hf/muller/power/...) is handled by rdmft_xc_type.
            XC_Functional::set_xc_type("hf");
            rdmft_xc_context_switched = true;
            if (inp.dft_functional != "hf")
            {
                GlobalV::ofs_running
                    << "RDMFT: using internal HF XC context for RDMFT solve (independent of dft_functional="
                    << inp.dft_functional << ")."
                    << std::endl;
            }
        }
        const int nk = this->pelec->wg.nr;
        const int nbands = this->pelec->wg.nc;
        const int nks = this->kv.get_nks(); // kv.wk has nks entries; ik >= nks wraps (spin-down)

        // Helper to map k-spin index ik to the kv.wk slot
        auto kv_wk = [&](int ik) -> double {
            return this->kv.wk[ik < nks ? ik : ik - nks];
        };

        // Build flat occupation vector n_{ik} = wg(ik, ib) / wk[ik]
        std::vector<double> occ_flat(nk * nbands, 0.0);
        for (int ik = 0; ik < nk; ++ik)
        {
            const double wk = kv_wk(ik);
            for (int ib = 0; ib < nbands; ++ib)
            {
                occ_flat[ik * nbands + ib] = (wk > 0.0)
                    ? this->pelec->wg(ik, ib) / wk : 0.0;
            }
        }

        // Build RDMFTConfig from input parameters
        rdmft::RDMFTConfig rdmft_config;
        rdmft_config.xc_type = rdmft_xc_type;
        rdmft_config.alpha_power = inp.rdmft_power_alpha;
        rdmft_config.outer_maxiter = inp.rdmft_outer_maxiter;
        rdmft_config.orb_maxiter = inp.rdmft_orb_maxiter;
        // occ_maxiter == 0 skips occupation optimization in alternating RDMFT (orbitals-only inner).
        rdmft_config.occ_maxiter = std::max(0, inp.rdmft_occ_maxiter);
        if (inp.rdmft_occ_init_mode == "perturbed")
            rdmft_config.occ_init_mode = rdmft::OccInitMode::Perturbed;
        else if (inp.rdmft_occ_init_mode == "binary")
            rdmft_config.occ_init_mode = rdmft::OccInitMode::Binary;
        else if (inp.rdmft_occ_init_mode == "uniform")
            rdmft_config.occ_init_mode = rdmft::OccInitMode::Uniform;
        else
            rdmft_config.occ_init_mode = rdmft::OccInitMode::KS;
        rdmft_config.occ_init_perturb = inp.rdmft_occ_init_perturb;
        rdmft_config.occ_init_nbands_top = inp.rdmft_occ_init_nbands_top;
        rdmft_config.energy_tol = inp.rdmft_energy_tol;
        rdmft_config.orb_grad_tol = inp.rdmft_orb_grad_tol;
        rdmft_config.orb_energy_tol = inp.rdmft_orb_energy_tol;
        rdmft_config.orb_ls_fixed_step = inp.rdmft_orb_ls_fixed_step;
        rdmft_config.orb_ls_stepsize = inp.rdmft_orb_ls_stepsize;
        rdmft_config.occ_ls_fixed_step = inp.rdmft_occ_ls_fixed_step;
        rdmft_config.occ_ls_stepsize = inp.rdmft_occ_ls_stepsize;
        rdmft_config.rdmft_occ_tol = inp.rdmft_occ_tol;
        rdmft_config.occ_energy_tol = inp.rdmft_occ_energy_tol;
        rdmft_config.occ_grad_tol = inp.rdmft_occ_grad_tol;
        rdmft_config.aug_lag_lambda_init = inp.rdmft_alm_lambda_init;
        rdmft_config.aug_lag_mu_init = inp.rdmft_alm_mu_init;
        rdmft_config.aug_lag_mu_factor = inp.rdmft_alm_mu_factor;
        rdmft_config.pg_occ_cg_ls_alpha_cap = inp.rdmft_pg_occ_cg_ls_alpha_cap;
        rdmft_config.pg_occ_ls_recovery_alpha = inp.rdmft_pg_occ_ls_recovery_alpha;
        if (inp.rdmft_occ_ls_init_step == "fixed")
            rdmft_config.occ_line_search_init_step = rdmft::LineSearchInitStep::Fixed;
        else if (inp.rdmft_occ_ls_init_step == "quad")
            rdmft_config.occ_line_search_init_step = rdmft::LineSearchInitStep::Quadratic;
        else
            rdmft_config.occ_line_search_init_step = rdmft::LineSearchInitStep::BarzilaiBorwein;
        rdmft_config.line_search_polynomial = inp.rdmft_line_search_polynomial;
        rdmft_config.line_search_c1 = inp.rdmft_line_search_c1;
        rdmft_config.line_search_c2 = inp.rdmft_line_search_c2;
        rdmft_config.line_search_max_zoom = inp.rdmft_line_search_max_zoom;
        rdmft_config.alm_bb_enabled = inp.rdmft_alm_bb_enabled;
        if (inp.rdmft_alm_bb_mode == "bb1")
            rdmft_config.alm_bb_mode = rdmft::BBStepMode::BB1;
        else if (inp.rdmft_alm_bb_mode == "bb2")
            rdmft_config.alm_bb_mode = rdmft::BBStepMode::BB2;
        else
            rdmft_config.alm_bb_mode = rdmft::BBStepMode::Alternate;
        rdmft_config.alm_bb_alpha_min = inp.rdmft_alm_bb_alpha_min;
        rdmft_config.alm_bb_alpha_max = inp.rdmft_alm_bb_alpha_max;
        rdmft_config.lbfgs_memory = inp.rdmft_lbfgs_memory;
        rdmft_config.adam_lr = inp.rdmft_adam_lr;
        rdmft_config.joint_orb_scale = inp.rdmft_joint_orb_scale;
        rdmft_config.print_stiefel_gram = inp.rdmft_print_stiefel_gram;
        rdmft_config.print_evals = inp.rdmft_print_evals;
        rdmft_config.print_orb_grad_decomp = inp.rdmft_print_orb_grad_decomp;
        rdmft_config.occ_entropy_gamma = inp.rdmft_occ_entropy_gamma;
        if (rdmft_config.xc_type != rdmft::XCFunctionalType::HF)
        {
            rdmft_config.occ_entropy_gamma = 0.0;
        }

        // Parse strategy. "joint" is the current name for the simultaneous
        // product-manifold optimisation; "product_manifold" is accepted as a
        // backwards-compatible alias.
        if (inp.rdmft_solver_strategy == "alternating")
            rdmft_config.strategy = rdmft::SolverStrategy::Alternating;
        else
            rdmft_config.strategy = rdmft::SolverStrategy::Joint;

        if (inp.rdmft_constraint == "projected_gradient")
            rdmft_config.constraint_method = rdmft::ConstraintMethod::ProjectedGradient;
        else if (inp.rdmft_constraint == "active_set")
            rdmft_config.constraint_method = rdmft::ConstraintMethod::ActiveSet;
        else
            rdmft_config.constraint_method = rdmft::ConstraintMethod::AugmentedLagrangian;

        if (inp.rdmft_occ_param == "logistic")
            rdmft_config.occ_param = rdmft::OccParamType::Logistic;
        else if (inp.rdmft_occ_param == "sigma_shift")
            rdmft_config.occ_param = rdmft::OccParamType::SigmaShift;
        else
            rdmft_config.occ_param = rdmft::OccParamType::CosineSq;

        rdmft_config.occ_optimizer = rdmft::parse_optimizer_input_or_quit(inp.rdmft_occ_optimizer,
                                                                            "rdmft_occ_optimizer",
                                                                            "ESolver_KS_LCAO::after_scf");
        rdmft_config.orb_optimizer = rdmft::parse_optimizer_input_or_quit(inp.rdmft_orb_optimizer,
                                                                          "rdmft_orb_optimizer",
                                                                          "ESolver_KS_LCAO::after_scf");
        rdmft_config.orb_cg_precond = inp.rdmft_orb_precond;
        rdmft_config.orb_cg_precond_delta = inp.rdmft_orb_precond_delta;
        rdmft_config.occ_ls_preset = rdmft::parse_line_search_preset_or_quit(inp.rdmft_occ_ls_type,
                                                                             "rdmft_occ_ls_type",
                                                                             "ESolver_KS_LCAO::after_scf");
        rdmft_config.orb_ls_preset = rdmft::parse_line_search_preset_or_quit(inp.rdmft_orb_ls_type,
                                                                             "rdmft_orb_ls_type",
                                                                             "ESolver_KS_LCAO::after_scf");
        rdmft_config.joint_optimizer = rdmft::parse_optimizer_input_or_quit(inp.rdmft_joint_optimizer,
                                                                            "rdmft_joint_optimizer",
                                                                            "ESolver_KS_LCAO::after_scf");
        rdmft_config.grad_check = inp.rdmft_grad_check;

        // RDMFT equality target N_e: default base is sum wg (matches loaded occupations);
        // optional base PARAM.inp.nelec plus rdmft_nelec_delta (see INPUT).
        // When sum(wg) ≈ 0 (e.g. scf_nmax=0 so calculate_weights never ran), fall
        // back to inp.nelec so RDMFT can still proceed with the correct electron count.
        double sum_wg = 0.0;
        for (int ik = 0; ik < nk; ++ik)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                sum_wg += this->pelec->wg(ik, ib);
            }
        }
        double nelec_base;
        if (inp.rdmft_nelec_use_input)
        {
            nelec_base = inp.nelec;
        }
        else if (sum_wg > 1.0e-12)
        {
            nelec_base = sum_wg;
        }
        else
        {
            nelec_base = inp.nelec;
            GlobalV::ofs_running
                << "RDMFT: sum(wg) ~ 0 (no KS occupations available, e.g. scf_nmax=0); "
                << "falling back to inp.nelec=" << inp.nelec << " for N_e target."
                << std::endl;
        }
        const double n_electrons = nelec_base + inp.rdmft_nelec_delta;
        if (n_electrons <= 0.0)
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_scf",
                                     "RDMFT electron target N_e must be positive. Check rdmft_nelec_use_input, "
                                     "nelec, sum(wg), and rdmft_nelec_delta.");
        }
        rdmft::RDMFTNelectronTargetMeta nelec_meta;
        nelec_meta.sum_initial_wg = sum_wg;
        nelec_meta.input_nelec = inp.nelec;
        nelec_meta.use_input_nelec = inp.rdmft_nelec_use_input;
        nelec_meta.rdmft_nelec_delta = inp.rdmft_nelec_delta;
        // Per-spin equality targets for nspin=2: split N_e via nupdown.  This
        // mirrors KS-LCAO `nelec_spin` and the split-Fermi treatment of
        // `nupdown` (`PARAM.globalv.two_fermi`); when nupdown is zero we
        // still feed two paired N_s = N_e/2 targets so the RDMFT solver
        // enforces the same number of equality constraints (two) as the
        // collinear KS reference.
        nelec_meta.nupdown = inp.nupdown;
        if (PARAM.inp.nspin == 2)
        {
            nelec_meta.two_fermi_active = true;
            // Prefer the per-spin sums of the KS occupation seed when they
            // are populated (rdmft_nelec_use_input == false): this preserves
            // the spin polarization KS converged to.  Otherwise fall back to
            // a symmetric split using PARAM.inp.nupdown (matching ABACUS
            // KS's `nelec_spin` initialiser).
            double n_up = 0.0, n_dn = 0.0;
            const int nks_kpt = this->kv.get_nks() / 2; // per-spin k-points
            if (!inp.rdmft_nelec_use_input && sum_wg > 1.0e-12 && nks_kpt > 0)
            {
                for (int ik = 0; ik < nk; ++ik)
                {
                    const int isk = (ik < nks_kpt) ? 0 : 1;
                    for (int ib = 0; ib < nbands; ++ib)
                    {
                        if (isk == 0) n_up += this->pelec->wg(ik, ib);
                        else          n_dn += this->pelec->wg(ik, ib);
                    }
                }
                // Apply rdmft_nelec_delta proportionally between spins so
                // total stays consistent with n_electrons.
                const double total_wg = n_up + n_dn;
                if (total_wg > 1.0e-12)
                {
                    const double scale = n_electrons / total_wg;
                    n_up *= scale;
                    n_dn *= scale;
                }
            }
            else
            {
                n_up = 0.5 * (n_electrons + inp.nupdown);
                n_dn = 0.5 * (n_electrons - inp.nupdown);
            }
            nelec_meta.n_electrons_per_spin = {n_up, n_dn};
        }
        rdmft::RDMFTSolver<TK, TR> rdmft_new_solver;
        rdmft_new_solver.init(rdmft_config, rdmft_eg, &this->kv, nbands, n_electrons, nelec_meta);
        this->rdmft_eg.set_occ_entropy_gamma(rdmft_config.occ_entropy_gamma);
        this->rdmft_eg.set_print_orb_grad_decomp(rdmft_config.print_orb_grad_decomp);

        // Run the optimization (`rdmft_grad_check` runs inside solve() after X-space setup)
        double etot_rdmft = rdmft_new_solver.solve(occ_flat, *this->psi);

        // Update pelec->wg from optimized occupations
        for (int ik = 0; ik < nk; ++ik)
        {
            const double wk = kv_wk(ik);
            for (int ib = 0; ib < nbands; ++ib)
            {
                this->pelec->wg(ik, ib) = occ_flat[ik * nbands + ib] * wk;
            }
        }

        // Update the total energy record
        this->pelec->f_en.etot = etot_rdmft;
        if (rdmft_xc_context_switched)
        {
            XC_Functional::set_xc_type(xc_func_restore);
        }

        ModuleBase::timer::end("ESolver_KS_LCAO", "rdmft_solve");
    }
#endif

    //! 2) output of lcao every few ionic steps
    ModuleIO::ctrl_scf_lcao<TK, TR>(ucell,
            PARAM.inp, this->kv, this->pelec, this->dmat.dm, this->pv,
            this->gd, this->psi, hamilt_lcao, this->dftu, this->two_center_bundle_,
            this->orb_, this->pw_wfc, this->pw_rho, this->pw_big, this->sf,
            this->pw_rhod, this->locpp.vloc, this->solvent,
            this->deepks, this->exx_nao,
            this->conv_esolver, this->scf_nmax_flag, istep);

    //! 3) Clean up RA, which is used to serach for adjacent atoms
    if (!PARAM.inp.cal_force && !PARAM.inp.cal_stress)
    {
        this->RA.delete_grid();
    }

    ModuleBase::timer::end("ESolver_KS_LCAO", "after_scf");
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver
