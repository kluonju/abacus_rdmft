//==========================================================
// Author: Jingang Han
// DATE : 2024-03-11
//==========================================================
//==========================================================
// Author: Kai Luo
// Email: kluo@njust.edu.cn
// DATE : April 2026
//==========================================================
#include "source_lcao/module_rdmft/rdmft_tools.h"
// used by class Veff_rdmft
#include "source_base/tool_title.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_estate/module_pot/H_Hartree_pw.h"
#include "source_estate/module_pot/pot_local.h"
#include "source_estate/module_pot/pot_xc.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_lcao/module_gint/gint_interface.h"

#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <cassert>

namespace rdmft
{

template <>
void conj_psi<double>(psi::Psi<double>& wfc) {}

template <>
void HkPsi<double>(const Parallel_Orbitals* ParaV, 
                   const double& HK, 
                   const double& wfc, 
                   double& H_wfc)
{
    const int one_int = 1;
    const double one_double = 1.0;
    const double zero_double = 0.0;
    const char N_char = 'N';
    const char C_char = 'C';

#ifdef __MPI
    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    //because wfc(bands, basis'), H(basis, basis'), we do wfc*H^T(in the perspective of cpp, not in fortran). And get H_wfc(bands, basis) is correct.
    ScalapackConnector::gemm( C_char, N_char, nbasis, nbands, nbasis, one_double, &HK, 1, 1, ParaV->desc,
        &wfc, 1, 1, ParaV->desc_wfc, zero_double, &H_wfc, 1, 1, ParaV->desc_wfc );
#endif

}


template <>
void cal_bra_op_ket<double>(const Parallel_Orbitals* ParaV, 
                            const Parallel_2D& para_Eij_in,
                            const double& wfc, 
                            const double& H_wfc, 
                            std::vector<double>& Dmn)
{
    const int one_int = 1;
    const double one_double = 1.0;
    const double zero_double = 0.0;
    const char N_char = 'N';
    const char T_char = 'T';

    const int nrow_bands = para_Eij_in.get_row_size();
    const int ncol_bands = para_Eij_in.get_col_size();

#ifdef __MPI
    const int nbasis = ParaV->desc[2];
    const int nbands = ParaV->desc_wfc[3];

    ScalapackConnector::gemm( T_char, N_char, nbands, nbands, nbasis, one_double, &wfc, 1, 1, ParaV->desc_wfc,
            &H_wfc, 1, 1, ParaV->desc_wfc, zero_double, &Dmn[0], 1, 1, para_Eij_in.desc );
#endif
}


// occNum_wfcHwfc = occNum*wfcHwfc + occNum_wfcHwfc
// See OccWeightMode for meaning of weight_mode. Default OccWeightMode::Occupation.
void occNum_Mul_wfcHwfc(const ModuleBase::matrix& occ_number, 
                            const ModuleBase::matrix& wfcHwfc, 
                            ModuleBase::matrix& occNum_wfcHwfc,
                            OccWeightMode weight_mode, 
                            const std::string XC_func_rdmft, 
                            const double alpha)
{
    for(int ir=0; ir<occ_number.nr; ++ ir)
    {
        for(int ic=0; ic<occ_number.nc; ++ic) 
        { 
            occNum_wfcHwfc(ir, ic) += occNum_func(occ_number(ir, ic), weight_mode, XC_func_rdmft, alpha) * wfcHwfc(ir, ic);
        }
    } 
}


// for the gradient of Etotal with respect to occupation numbers
void add_occNum(const K_Vectors& kv, 
                    const ModuleBase::matrix& occ_number, 
                    const ModuleBase::matrix& wfcHwfc_TV_in, 
                    const ModuleBase::matrix& wfcHwfc_hartree_in,
                    const ModuleBase::matrix& wfcHwfc_dft_XC_in, 
                    const ModuleBase::matrix& wfcHwfc_exx_XC_in, 
                    ModuleBase::matrix& occNum_wfcHwfc, 
                    const std::string XC_func_rdmft, 
                    const double alpha)
{ 
    occNum_wfcHwfc.zero_out();
    occNum_Mul_wfcHwfc(occ_number, wfcHwfc_exx_XC_in, occNum_wfcHwfc, OccWeightMode::CouplingDerivative, XC_func_rdmft, alpha);
    occNum_wfcHwfc+=(wfcHwfc_TV_in);
    occNum_wfcHwfc+=(wfcHwfc_hartree_in);
    occNum_wfcHwfc+=(wfcHwfc_dft_XC_in);

    // consider W_k for dE/d_occNum
    for(int ik=0; ik<occ_number.nr; ++ik)
    {
        for(int inb=0; inb<occ_number.nc; ++inb) 
        { 
            occNum_wfcHwfc(ik, inb) *= kv.wk[ik];
        }
    } 
}


//! give certain occNum_wfcHwfc, get the corresponding energy
double getEnergy(const ModuleBase::matrix& occNum_wfcHwfc)
{
    double energy = 0.0;
    for(int ir=0; ir<occNum_wfcHwfc.nr; ++ ir)
    {
        for(int ic=0; ic<occNum_wfcHwfc.nc; ++ic) 
        { 
            energy += occNum_wfcHwfc(ir, ic);
        }
    }
    return energy;
}


//! for HF, Muller and power functional, g(eta) = eta, eta^0.5, eta^alpha respectively.
//! weight_mode controls the return value:
//!   Occupation        → eta
//!   HalfOccupation    → 0.5*eta
//!   Coupling          → g(eta) = eta^alpha
//!   HalfCoupling      → 0.5*g(eta)
//!   CouplingDerivative→ d_g(eta)/d_eta = alpha*eta^(alpha-1)
//!   Unity             → 1.0
//! Default weight_mode=OccWeightMode::Occupation, XC_func_rdmft="hf", alpha=0.656
double occNum_func(const double eta, const OccWeightMode weight_mode, const std::string XC_func_rdmft, double alpha)
{
    if( XC_func_rdmft == "power" || XC_func_rdmft == "wp22" || XC_func_rdmft == "cwp22" ) { ; }
    else if( XC_func_rdmft == "muller" ) { alpha = 0.5; }
    else { alpha = 1.0; }

    switch (weight_mode)
    {
        case OccWeightMode::Occupation:          return eta;
        case OccWeightMode::HalfOccupation:      return 0.5 * eta;
        case OccWeightMode::Coupling:            return std::pow(eta, alpha);
        case OccWeightMode::HalfCoupling:        return 0.5 * std::pow(eta, alpha);
        case OccWeightMode::CouplingDerivative:  return alpha * std::pow(eta, alpha - 1.0);
        case OccWeightMode::Unity:               return 1.0;
    }

    // unreachable — default to raw occupation
    return eta;
    
}

// this part of the code is copying from class Veff
// initialize_HR()
template <typename TK, typename TR>
void Veff_rdmft<TK, TR>::initialize_HR(const UnitCell* ucell_in, const Grid_Driver* GridD)
{
    ModuleBase::TITLE("Veff", "initialize_HR");
    ModuleBase::timer::start("Veff", "initialize_HR");

    this->nspin = PARAM.inp.nspin;
    auto* paraV = this->hR->get_paraV();// get parallel orbitals from HR
    // TODO: if paraV is nullptr, AtomPair can not use paraV for constructor, I will repair it in the future.

    for (int iat1 = 0; iat1 < ucell_in->nat; iat1++)
    {
        auto tau1 = ucell_in->get_tau(iat1);
        int T1 = 0;
        int I1 = 0;
        ucell_in->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell_in, tau1, T1, I1, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T2 = adjs.ntype[ad1];
            const int I2 = adjs.natom[ad1];
            const int iat2 = ucell_in->itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad1];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (ucell_in->cal_dtau(iat1, iat2, R_index2).norm() * ucell_in->lat0
                < orb_cutoff_[T1] + orb_cutoff_[T2])
            {
                hamilt::AtomPair<TR> tmp(iat1, iat2, R_index2, paraV);
                this->hR->insert_pair(tmp);
            }
        }
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    this->hR->allocate(nullptr, true);

    ModuleBase::timer::end("Veff", "initialize_HR");
}


namespace {
// Unified HR builder shared between the gamma-only (double) and multi-k
// (complex) specialisations of Veff_rdmft::contributeHR(). The loop over
// spin and the branch on `potential_` is identical; only the template
// arguments of the enclosing class differ.
template <typename TR_>
void build_HR_gint(const std::string& potential,
                   const int nspin,
                   const Charge* charge,
                   const UnitCell* ucell,
                   const ModulePW::PW_Basis* rho_basis,
                   const ModuleBase::matrix* vloc,
                   const ModuleBase::ComplexMatrix* sf,
                   double* etxc,
                   double* vtxc,
                   hamilt::HContainer<TR_>* hR)
{
    double* vr_eff = nullptr;

    // RDMFT's HContainer carries a single spin slot, and V_H / V_xc are the
    // same for both spin components in the collinear nspin=2 case
    // (PotHartree returns v(0)=v(1)). Looping over `is` and calling
    // cal_gint_vl(v(is), hR) for each spin doubles the matrix elements,
    // breaking the RDMFT energy / gradient identities for nspin=2. We add
    // exactly once instead.
    if (potential == "hartree")
    {
        ModuleBase::matrix v_matrix(nspin, charge->nrxx);
        elecstate::PotHartree potH(rho_basis);
        potH.cal_v_eff(charge, ucell, v_matrix);
        vr_eff = &v_matrix(0, 0);
        ModuleGint::cal_gint_vl(vr_eff, hR);
    }
    else if (potential == "local")
    {
        double vlocal_of_0 = 0.0;
        ModuleBase::matrix v_matrix(1, charge->nrxx);
        elecstate::PotLocal potL(vloc, sf, rho_basis, vlocal_of_0);
        potL.cal_fixed_v(&v_matrix(0, 0));
        vr_eff = &v_matrix(0, 0);
        ModuleGint::cal_gint_vl(vr_eff, hR);
    }
    else if (potential == "xc")
    {
        // meta-GGA not supported yet for Veff_rdmft.  Same single-slot
        // single-add convention as the "hartree" branch above.
        ModuleBase::matrix vofk = *vloc;
        vofk.zero_out();
        ModuleBase::matrix v_matrix(nspin, charge->nrxx);
        elecstate::PotXC potXC(rho_basis, etxc, vtxc, &vofk);
        potXC.cal_v_eff(charge, ucell, v_matrix);
        vr_eff = &v_matrix(0, 0);
        ModuleGint::cal_gint_vl(vr_eff, hR);
    }
    else
    {
        ModuleBase::WARNING("Veff_rdmft", "unknown potential type " + potential);
    }
}
} // anonymous namespace

template <>
void Veff_rdmft<std::complex<double>, double>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::start("Veff", "contributeHR");

    build_HR_gint(potential_, this->nspin, charge_, ucell, rho_basis_, vloc_, sf_,
                  etxc, vtxc, this->hR);

    if (this->nspin == 2) { this->current_spin = 1 - this->current_spin; }

    ModuleBase::timer::end("Veff", "contributeHR");
}

template <>
void Veff_rdmft<std::complex<double>, std::complex<double>>::contributeHR()
{
    // nspin = 4 case not implemented currently.
}

template <>
void Veff_rdmft<double, double>::contributeHR()
{
    ModuleBase::TITLE("Veff", "contributeHR");
    ModuleBase::timer::start("Veff", "contributeHR");

    build_HR_gint(potential_, this->nspin, charge_, ucell, rho_basis_, vloc_, sf_,
                  etxc, vtxc, this->hR);

    this->new_e_iteration = false;
    if (this->nspin == 2) { this->current_spin = 1 - this->current_spin; }

    ModuleBase::timer::end("Veff", "contributeHR");
}

}
template class rdmft::Veff_rdmft<double, double>;

template class rdmft::Veff_rdmft<std::complex<double>, double>;

template class rdmft::Veff_rdmft<std::complex<double>, std::complex<double>>;


