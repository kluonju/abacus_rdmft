#include "rdmft_energy_gradient.h"
#include "source_base/timer.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_estate/module_charge/symmetry_rho.h"
#include "source_estate/module_pot/H_Hartree_pw.h"
#include "source_estate/module_pot/pot_local.h"
#include "source_estate/module_pot/pot_xc.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_lcao/module_gint/gint_interface.h"
#include "source_lcao/module_operator_lcao/ekinetic.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"

#ifdef __EXX
#include "source_lcao/module_ri/RI_2D_Comm.h"
#include "source_lcao/module_operator_lcao/op_exx_lcao.h"
#endif

#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <type_traits>

namespace rdmft
{

namespace detail {
// ScaLAPACK dispatch helpers for real and complex types
inline void pgemm_wrapper(char ta, char tb, int m, int n, int k,
    double alpha, const double* A, int ia, int ja, const int* descA,
    const double* B, int ib, int jb, const int* descB,
    double beta, double* C, int ic, int jc, const int* descC)
{
    ScalapackConnector::gemm(ta, tb, m, n, k, alpha, A, ia, ja, descA,
                              B, ib, jb, descB, beta, C, ic, jc, descC);
}
inline void pgemm_wrapper(char ta, char tb, int m, int n, int k,
    std::complex<double> alpha, const std::complex<double>* A, int ia, int ja, const int* descA,
    const std::complex<double>* B, int ib, int jb, const int* descB,
    std::complex<double> beta, std::complex<double>* C, int ic, int jc, const int* descC)
{
    ScalapackConnector::gemm(ta, tb, m, n, k, alpha, A, ia, ja, descA,
                              B, ib, jb, descB, beta, C, ic, jc, descC);
}

inline double real_of(double v) { return v; }
inline double real_of(std::complex<double> v) { return v.real(); }

// Transpose character for real/complex
inline char trans_char(double) { return 'T'; }
inline char trans_char(std::complex<double>) { return 'C'; }
} // namespace detail

// ---- Forward declarations for Veff_rdmft (local to this TU) ----
template <typename TK, typename TR>
class Veff_rdmft_local : public hamilt::OperatorLCAO<TK, TR>
{
  public:
    Veff_rdmft_local(hamilt::HS_Matrix_K<TK>* hsk_in,
                     const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                     elecstate::Potential* pot_in,
                     hamilt::HContainer<TR>* hR_in,
                     const UnitCell* ucell_in,
                     const std::vector<double>& orb_cutoff,
                     const Grid_Driver* GridD_in,
                     int nspin,
                     const Charge* charge_in,
                     const ModulePW::PW_Basis* rho_basis_in,
                     const ModuleBase::matrix* vloc_in,
                     const ModuleBase::ComplexMatrix* sf_in,
                     const std::string& potential_in,
                     double* etxc_in = nullptr,
                     double* vtxc_in = nullptr);

    ~Veff_rdmft_local() = default;

    void contributeHR() override;

  private:
    void initialize_HR(const UnitCell* ucell_in, const Grid_Driver* GridD);

    std::vector<double> orb_cutoff_;
    elecstate::Potential* pot_ = nullptr;
    const UnitCell* ucell_ = nullptr;
    const Grid_Driver* gd_ = nullptr;
    int nspin_ = 1;
    const Charge* charge_ = nullptr;
    const ModulePW::PW_Basis* rho_basis_ = nullptr;
    const ModuleBase::matrix* vloc_ = nullptr;
    const ModuleBase::ComplexMatrix* sf_ = nullptr;
    std::string potential_;
    double* etxc_ = nullptr;
    double* vtxc_ = nullptr;
};

// ---- Veff_rdmft_local implementation ----

template <typename TK, typename TR>
Veff_rdmft_local<TK, TR>::Veff_rdmft_local(
    hamilt::HS_Matrix_K<TK>* hsk_in,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
    elecstate::Potential* pot_in,
    hamilt::HContainer<TR>* hR_in,
    const UnitCell* ucell_in,
    const std::vector<double>& orb_cutoff,
    const Grid_Driver* GridD_in,
    int nspin,
    const Charge* charge_in,
    const ModulePW::PW_Basis* rho_basis_in,
    const ModuleBase::matrix* vloc_in,
    const ModuleBase::ComplexMatrix* sf_in,
    const std::string& potential_in,
    double* etxc_in,
    double* vtxc_in)
    : orb_cutoff_(orb_cutoff), pot_(pot_in), ucell_(ucell_in), gd_(GridD_in),
      hamilt::OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in),
      nspin_(nspin), charge_(charge_in), rho_basis_(rho_basis_in),
      vloc_(vloc_in), sf_(sf_in), potential_(potential_in),
      etxc_(etxc_in), vtxc_(vtxc_in)
{
    this->cal_type = hamilt::calculation_type::lcao_gint;
    this->initialize_HR(ucell_in, GridD_in);
}

template <typename TK, typename TR>
void Veff_rdmft_local<TK, TR>::initialize_HR(const UnitCell* ucell_in, const Grid_Driver* GridD)
{
    auto* paraV = this->hR->get_paraV();
    for (int iat1 = 0; iat1 < ucell_in->nat; iat1++)
    {
        auto tau1 = ucell_in->get_tau(iat1);
        int T1 = 0, I1 = 0;
        ucell_in->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell_in, tau1, T1, I1, &adjs);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T2 = adjs.ntype[ad1];
            const int I2 = adjs.natom[ad1];
            const int iat2 = ucell_in->itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
                continue;
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad1];
            if (ucell_in->cal_dtau(iat1, iat2, R_index2).norm() * ucell_in->lat0
                < orb_cutoff_[T1] + orb_cutoff_[T2])
            {
                hamilt::AtomPair<TR> tmp(iat1, iat2, R_index2, paraV);
                this->hR->insert_pair(tmp);
            }
        }
    }
    this->hR->allocate(nullptr, true);
}

// Explicit specialization for complex<double>, double (nspin=1,2, multi-k)
template <>
void Veff_rdmft_local<std::complex<double>, double>::contributeHR()
{
    double* vr_eff = nullptr;
    if (potential_ == "hartree")
    {
        ModuleBase::matrix v_hartree(nspin_, charge_->nrxx);
        elecstate::PotHartree potH(rho_basis_);
        potH.cal_v_eff(charge_, ucell_, v_hartree);
        for (int is = 0; is < nspin_; ++is)
        {
            vr_eff = &v_hartree(is, 0);
            ModuleGint::cal_gint_vl(vr_eff, this->hR);
        }
    }
    else if (potential_ == "local")
    {
        double vlocal_of_0 = 0.0;
        ModuleBase::matrix v_local(1, charge_->nrxx);
        elecstate::PotLocal potL(vloc_, sf_, rho_basis_, vlocal_of_0);
        potL.cal_fixed_v(&v_local(0, 0));
        vr_eff = &v_local(0, 0);
        ModuleGint::cal_gint_vl(vr_eff, this->hR);
    }
    else if (potential_ == "xc")
    {
        ModuleBase::matrix vofk = *vloc_;
        vofk.zero_out();
        ModuleBase::matrix v_xc(nspin_, charge_->nrxx);
        elecstate::PotXC potXC(rho_basis_, etxc_, vtxc_, &vofk);
        potXC.cal_v_eff(charge_, ucell_, v_xc);
        for (int is = 0; is < nspin_; ++is)
        {
            vr_eff = &v_xc(is, 0);
            ModuleGint::cal_gint_vl(vr_eff, this->hR);
        }
    }
}

template <>
void Veff_rdmft_local<double, double>::contributeHR()
{
    double* vr_eff = nullptr;
    if (potential_ == "hartree")
    {
        ModuleBase::matrix v_hartree(nspin_, charge_->nrxx);
        elecstate::PotHartree potH(rho_basis_);
        potH.cal_v_eff(charge_, ucell_, v_hartree);
        for (int is = 0; is < nspin_; ++is)
        {
            vr_eff = &v_hartree(is, 0);
            ModuleGint::cal_gint_vl(vr_eff, this->hR);
        }
    }
    else if (potential_ == "local")
    {
        double vlocal_of_0 = 0.0;
        ModuleBase::matrix v_local(1, charge_->nrxx);
        elecstate::PotLocal potL(vloc_, sf_, rho_basis_, vlocal_of_0);
        potL.cal_fixed_v(&v_local(0, 0));
        vr_eff = &v_local(0, 0);
        ModuleGint::cal_gint_vl(vr_eff, this->hR);
    }
    else if (potential_ == "xc")
    {
        ModuleBase::matrix vofk = *vloc_;
        vofk.zero_out();
        ModuleBase::matrix v_xc(nspin_, charge_->nrxx);
        elecstate::PotXC potXC(rho_basis_, etxc_, vtxc_, &vofk);
        potXC.cal_v_eff(charge_, ucell_, v_xc);
        for (int is = 0; is < nspin_; ++is)
        {
            vr_eff = &v_xc(is, 0);
            ModuleGint::cal_gint_vl(vr_eff, this->hR);
        }
    }
}

template <>
void Veff_rdmft_local<std::complex<double>, std::complex<double>>::contributeHR()
{
    // nspin=4 not implemented
}

// ---- EnergyGradient implementation ----

template <typename TK, typename TR>
EnergyGradient<TK, TR>::~EnergyGradient()
{
    delete HR_one_;
    delete HR_hartree_;
    delete HR_exx_;
    delete hsk_one_;
    delete hsk_hartree_;
    delete hsk_exx_;
    delete op_ekinetic_;
    delete op_nonlocal_;
    delete op_local_;
    delete op_hartree_;
    delete op_exx_;
#ifdef __EXX
    delete exx_lri_d_;
    delete exx_lri_c_;
#endif
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::init(
    const Parallel_Orbitals* ParaV_in,
    const UnitCell* ucell_in,
    const Grid_Driver* gd_in,
    const K_Vectors* kv_in,
    elecstate::ElecState* pelec_in,
    const LCAO_Orbitals* orb_in,
    const TwoCenterBundle* two_center_bundle_in,
    const XCFunctional& xc_func_in)
{
    ParaV_ = ParaV_in;
    ucell_ = ucell_in;
    gd_ = gd_in;
    kv_ = kv_in;
    pelec_ = pelec_in;
    charge_ = pelec_in->charge;
    orb_ = orb_in;
    two_center_bundle_ = two_center_bundle_in;
    xc_func_ = xc_func_in;

    nspin_ = PARAM.inp.nspin;
    nbands_ = PARAM.inp.nbands;
    nk_ = ModuleSymmetry::Symmetry::symm_flag == -1
              ? kv_->get_nkstot_full()
              : kv_->get_nks();
    nk_ *= nspin_;
    nbasis_local_ = ParaV_->nrow;

#ifdef __MPI
    para_Eij_.set(nbands_, nbands_, ParaV_->nb, ParaV_->blacs_ctxt);
#endif

    // Allocate Hamiltonian containers
    HR_one_ = new hamilt::HContainer<TR>(*ucell_, ParaV_);
    HR_hartree_ = new hamilt::HContainer<TR>(*ucell_, ParaV_);
    HR_exx_ = new hamilt::HContainer<TR>(*ucell_, ParaV_);

    hsk_one_ = new hamilt::HS_Matrix_K<TK>(ParaV_, true);
    hsk_hartree_ = new hamilt::HS_Matrix_K<TK>(ParaV_, true);
    hsk_exx_ = new hamilt::HS_Matrix_K<TK>(ParaV_, true);

    if (PARAM.inp.gamma_only)
    {
        HR_one_->fix_gamma();
        HR_hartree_->fix_gamma();
        HR_exx_->fix_gamma();
    }

#ifdef __EXX
    if (PARAM.inp.rdmft)
    {
        // RDMFT manages its own EXX infrastructure, independent of the
        // global cal_exx flag which controls the KS-SCF hybrid loop.
        // Ensure Coulomb parameters are set for HF-type exchange.
        if (GlobalC::exx_info.info_global.coulomb_param.empty())
        {
            std::string sing_corr = PARAM.inp.exx_singularity_correction;
            if (sing_corr == "default") { sing_corr = "spencer"; }
            GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock] = {{
                {"alpha", "1"},
                {"singularity_correction", sing_corr}
            }};
        }

        exx_spacegroup_symmetry_ = (PARAM.inp.nspin < 4
                                    && ModuleSymmetry::Symmetry::symm_flag == 1);
        if (exx_spacegroup_symmetry_)
        {
            const std::array<int, 3>& period = RI_Util::get_Born_vonKarmen_period(*kv_);
            symrot_exx_.find_irreducible_sector(ucell_->symm, ucell_->atoms, ucell_->st,
                RI_Util::get_Born_von_Karmen_cells(period), period, ucell_->lat);
            symrot_exx_.cal_Ms(*kv_, *ucell_, *ParaV_);
        }

        if (GlobalC::exx_info.info_ri.real_number)
        {
            exx_lri_d_ = new Exx_LRI<double>(GlobalC::exx_info.info_ri);
            exx_lri_d_->init(MPI_COMM_WORLD, const_cast<UnitCell&>(*ucell_), *kv_, *orb_);
        }
        else
        {
            exx_lri_c_ = new Exx_LRI<std::complex<double>>(GlobalC::exx_info.info_ri);
            exx_lri_c_->init(MPI_COMM_WORLD, const_cast<UnitCell&>(*ucell_), *kv_, *orb_);
        }
        exx_enabled_ = true;
    }
#endif
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::update_ion(
    const UnitCell& ucell,
    const ModulePW::PW_Basis& rho_basis,
    const ModuleBase::matrix& vloc,
    const ModuleBase::ComplexMatrix& sf)
{
    ucell_ = &ucell;
    rho_basis_ = &rho_basis;
    vloc_ = &vloc;
    sf_ = &sf;

    HR_one_->set_zero();

    // Build one-body operators (kinetic + nonlocal + local pseudopotential)
    delete op_ekinetic_;
    op_ekinetic_ = new hamilt::EKinetic<hamilt::OperatorLCAO<TK, TR>>(
        hsk_one_, kv_->kvec_d, HR_one_, ucell_,
        orb_->cutoffs(), gd_, two_center_bundle_->kinetic_orb.get());

    delete op_nonlocal_;
    op_nonlocal_ = new hamilt::Nonlocal<hamilt::OperatorLCAO<TK, TR>>(
        hsk_one_, kv_->kvec_d, HR_one_, ucell_,
        orb_->cutoffs(), gd_, two_center_bundle_->overlap_orb_beta.get());

    delete op_local_;
    op_local_ = new Veff_rdmft_local<TK, TR>(
        hsk_one_, kv_->kvec_d, pelec_->pot, HR_one_, ucell_,
        orb_->cutoffs(), gd_, nspin_, charge_, rho_basis_, vloc_, sf_, "local");

    op_ekinetic_->contributeHR();
    op_nonlocal_->contributeHR();
    op_local_->contributeHR();

#ifdef __EXX
    if (exx_enabled_)
    {
        if (GlobalC::exx_info.info_ri.real_number)
            exx_lri_d_->cal_exx_ions(const_cast<UnitCell&>(ucell));
        else
            exx_lri_c_->cal_exx_ions(const_cast<UnitCell&>(ucell));
    }
#endif

    ion_initialized_ = true;
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::build_charge(
    const std::vector<double>& occ_flat,
    const psi::Psi<TK>& wfc)
{
    // Build wg matrix: wg(ik, ib) = wk[ik] * n(ik, ib)
    ModuleBase::matrix wg(nk_, nbands_);
    for (int ik = 0; ik < nk_; ++ik)
        for (int ib = 0; ib < nbands_; ++ib)
            wg(ik, ib) = kv_->wk[ik] * occ_flat[ik * nbands_ + ib];

    if (PARAM.inp.gamma_only)
    {
        elecstate::DensityMatrix<TK, double> DM(ParaV_, nspin_);
        elecstate::cal_dm_psi(ParaV_, wg, wfc, DM);
        DM.init_DMR(gd_, ucell_);
        DM.cal_DMR();
        for (int is = 0; is < nspin_; is++)
            ModuleBase::GlobalFunc::ZEROS(charge_->rho[is], charge_->nrxx);
        ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin_, charge_->rho);
    }
    else
    {
        elecstate::DensityMatrix<TK, double> DM(ParaV_, nspin_, kv_->kvec_d, nk_);
        elecstate::cal_dm_psi(ParaV_, wg, wfc, DM);
        DM.init_DMR(gd_, ucell_);
        DM.cal_DMR();
        for (int is = 0; is < nspin_; is++)
            ModuleBase::GlobalFunc::ZEROS(charge_->rho[is], charge_->nrxx);
        ModuleGint::cal_gint_rho(DM.get_DMR_vector(), nspin_, charge_->rho);
    }

    charge_->renormalize_rho();

    Symmetry_rho srho;
    for (int is = 0; is < nspin_; is++)
        srho.begin(is, *charge_, rho_basis_, const_cast<ModuleSymmetry::Symmetry&>(ucell_->symm));
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::build_DM_xc(
    const std::vector<double>& occ_flat,
    const psi::Psi<TK>& wfc,
    std::vector<std::vector<TK>>& DM_XC)
{
    DM_XC.resize(nk_, std::vector<TK>(ParaV_->nloc, TK(0)));

    // Build wk_g matrix: wk[ik] * g(n(ik, ib))
    ModuleBase::matrix wk_g(nk_, nbands_);
    for (int ik = 0; ik < nk_; ++ik)
        for (int ib = 0; ib < nbands_; ++ib)
            wk_g(ik, ib) = kv_->wk[ik] * xc_func_.g(occ_flat[ik * nbands_ + ib]);

    // DM_XC(k) = sum_i wk*g(n_i) * conj(C_i) * C_i^T
    // Using psiMulPsi with the modified weights
    psi::Psi<TK> wfc_copy(wfc);
    // Scale each orbital by sqrt(wk*g(n)) then use standard DM construction
    // Actually, build it directly:
    for (int ik = 0; ik < nk_; ++ik)
    {
        // Create scaled wfc: conj(C) * wk*g(n)
        const int nb_local = wfc.get_nbands();
        const int nbs_local = wfc.get_nbasis();

        // Use cal_dm_psi with wk_g as weights
        // This gives DM_XC(k) = sum_i wk_g(ik,i) * C(ik,i) * C(ik,i)^H
    }

    // Use ABACUS DM infrastructure
    if (PARAM.inp.gamma_only)
    {
        elecstate::DensityMatrix<TK, double> DM_xc(ParaV_, nspin_);
        elecstate::cal_dm_psi(ParaV_, wk_g, wfc, DM_xc);
        for (int ik = 0; ik < nk_; ++ik)
        {
            TK* dmk = DM_xc.get_DMK_pointer(ik);
            std::copy(dmk, dmk + ParaV_->nloc, DM_XC[ik].begin());
        }
    }
    else
    {
        elecstate::DensityMatrix<TK, double> DM_xc(ParaV_, nspin_, kv_->kvec_d, nk_);
        elecstate::cal_dm_psi(ParaV_, wk_g, wfc, DM_xc);
        for (int ik = 0; ik < nk_; ++ik)
        {
            TK* dmk = DM_xc.get_DMK_pointer(ik);
            std::copy(dmk, dmk + ParaV_->nloc, DM_XC[ik].begin());
        }
    }
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::apply_Hk(const TK* HK, const TK* psi_k, TK* Hpsi_k) const
{
#ifdef __MPI
    const int nbasis = ParaV_->desc[2];
    const int nbands = ParaV_->desc_wfc[3];
    const TK one = TK(1.0);
    const TK zero = TK(0.0);
    char tc = detail::trans_char(TK());

    detail::pgemm_wrapper(tc, 'N', nbasis, nbands, nbasis,
        one, HK, 1, 1, ParaV_->desc,
        psi_k, 1, 1, ParaV_->desc_wfc,
        zero, Hpsi_k, 1, 1, ParaV_->desc_wfc);
#endif
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::compute_diagonal(
    const TK* psi_k, const TK* Hpsi_k, double* diag, int ik) const
{
#ifdef __MPI
    const int nbands = ParaV_->desc_wfc[3];
    const int nbasis = ParaV_->desc[2];
    const TK one = TK(1.0);
    const TK zero = TK(0.0);
    char tc = detail::trans_char(TK());

    // Compute Eij = psi^H * Hpsi
    std::vector<TK> Eij(para_Eij_.get_row_size() * para_Eij_.get_col_size(), TK(0));

    detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, psi_k, 1, 1, ParaV_->desc_wfc,
        Hpsi_k, 1, 1, ParaV_->desc_wfc,
        zero, Eij.data(), 1, 1, para_Eij_.desc);

    // Extract diagonal
    const int nrow = para_Eij_.get_row_size();
    const int ncol = para_Eij_.get_col_size();
    for (int i = 0; i < nrow; ++i)
    {
        int ig = para_Eij_.local2global_row(i);
        for (int j = 0; j < ncol; ++j)
        {
            int jg = para_Eij_.local2global_col(j);
            if (ig == jg)
                diag[jg] = detail::real_of(Eij[i + j * nrow]);
        }
    }
#endif
}

template <typename TK, typename TR>
double EnergyGradient<TK, TR>::compute(
    const std::vector<double>& occ_flat,
    const psi::Psi<TK>& wfc,
    std::vector<double>& grad_occ,
    psi::Psi<TK>& grad_wfc)
{
    ModuleBase::timer::tick("RDMFT_EG", "compute");

    assert(ion_initialized_);

    // 1. Build charge density from occupations and orbitals
    build_charge(occ_flat, wfc);

    // 2. Build Hartree potential
    HR_hartree_->set_zero();
    delete op_hartree_;
    op_hartree_ = new Veff_rdmft_local<TK, TR>(
        hsk_hartree_, kv_->kvec_d, pelec_->pot, HR_hartree_, ucell_,
        orb_->cutoffs(), gd_, nspin_, charge_, rho_basis_, vloc_, sf_, "hartree");
    op_hartree_->contributeHR();

    // 3. Build exchange from modified DM (RDMFT's private EXX)
#ifdef __EXX
    if (exx_enabled_)
    {
        HR_exx_->set_zero();
        std::vector<std::vector<TK>> DM_XC;
        build_DM_xc(occ_flat, wfc, DM_XC);

        if (exx_spacegroup_symmetry_)
            DM_XC = symrot_exx_.restore_dm(*kv_, DM_XC, *ParaV_);

        std::vector<const std::vector<TK>*> DM_XC_ptr(DM_XC.size());
        for (size_t ik = 0; ik < DM_XC.size(); ++ik)
            DM_XC_ptr[ik] = &DM_XC[ik];

        if (GlobalC::exx_info.info_ri.real_number)
        {
            auto Ds = std::is_same<TK, double>::value
                ? RI_2D_Comm::split_m2D_ktoR<double>(*ucell_, *kv_, DM_XC_ptr, *ParaV_, nspin_)
                : RI_2D_Comm::split_m2D_ktoR<double>(*ucell_, *kv_, DM_XC_ptr, *ParaV_, nspin_,
                                                      exx_spacegroup_symmetry_);
            if (exx_spacegroup_symmetry_ && GlobalC::exx_info.info_ri.exx_symmetry_realspace)
                exx_lri_d_->cal_exx_elec(Ds, *ucell_, *ParaV_, &symrot_exx_);
            else
                exx_lri_d_->cal_exx_elec(Ds, *ucell_, *ParaV_);

            delete op_exx_;
            op_exx_ = new hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>>(
                hsk_exx_, HR_exx_, *ucell_, *kv_,
                &exx_lri_d_->Hexxs, nullptr, hamilt::Add_Hexx_Type::k);
        }
        else
        {
            auto Ds = std::is_same<TK, double>::value
                ? RI_2D_Comm::split_m2D_ktoR<std::complex<double>>(*ucell_, *kv_, DM_XC_ptr, *ParaV_, nspin_)
                : RI_2D_Comm::split_m2D_ktoR<std::complex<double>>(*ucell_, *kv_, DM_XC_ptr, *ParaV_, nspin_,
                                                                     exx_spacegroup_symmetry_);
            if (exx_spacegroup_symmetry_ && GlobalC::exx_info.info_ri.exx_symmetry_realspace)
                exx_lri_c_->cal_exx_elec(Ds, *ucell_, *ParaV_, &symrot_exx_);
            else
                exx_lri_c_->cal_exx_elec(Ds, *ucell_, *ParaV_);

            delete op_exx_;
            op_exx_ = new hamilt::OperatorEXX<hamilt::OperatorLCAO<TK, TR>>(
                hsk_exx_, HR_exx_, *ucell_, *kv_,
                nullptr, &exx_lri_c_->Hexxs, hamilt::Add_Hexx_Type::k);
        }
    
    }
#endif

    // 4. For each k-point: compute H*psi, diagonal elements, and assemble gradients
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();

    grad_occ.assign(nk_ * nbands_, 0.0);
    grad_wfc.resize(nk_, nb_local, nbs_local);
    grad_wfc.zero_out();

    E_one_ = 0.0;
    E_hartree_ = 0.0;
    E_xc_ = 0.0;

    cached_h_one_diag_.resize(nk_);
    std::vector<double> h_one_diag(nbands_, 0.0);
    std::vector<double> vh_diag(nbands_, 0.0);
    std::vector<double> vx_diag(nbands_, 0.0);

    std::vector<TK> Hpsi_one(nb_local * nbs_local);
    std::vector<TK> Hpsi_h(nb_local * nbs_local);
    std::vector<TK> Hpsi_x(nb_local * nbs_local);

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* psi_k = &wfc(ik, 0, 0);

        // One-body: H_one * psi
        hsk_one_->set_zero_hk();
        op_local_->contributeHk(ik);
        std::fill(Hpsi_one.begin(), Hpsi_one.end(), TK(0));
        apply_Hk(hsk_one_->get_hk(), psi_k, Hpsi_one.data());
        std::fill(h_one_diag.begin(), h_one_diag.end(), 0.0);
        compute_diagonal(psi_k, Hpsi_one.data(), h_one_diag.data(), ik);

        cached_h_one_diag_[ik].assign(h_one_diag.begin(), h_one_diag.end());
        hone_cache_valid_ = true;

        // Hartree: V_H * psi
        hsk_hartree_->set_zero_hk();
        op_hartree_->contributeHk(ik);
        std::fill(Hpsi_h.begin(), Hpsi_h.end(), TK(0));
        apply_Hk(hsk_hartree_->get_hk(), psi_k, Hpsi_h.data());
        std::fill(vh_diag.begin(), vh_diag.end(), 0.0);
        compute_diagonal(psi_k, Hpsi_h.data(), vh_diag.data(), ik);

        // Exchange: H_exx * psi (RDMFT's private EXX)
        std::fill(Hpsi_x.begin(), Hpsi_x.end(), TK(0));
        std::fill(vx_diag.begin(), vx_diag.end(), 0.0);
#ifdef __EXX
        if (exx_enabled_ && op_exx_)
        {
            hsk_exx_->set_zero_hk();
            op_exx_->contributeHk(ik);
            apply_Hk(hsk_exx_->get_hk(), psi_k, Hpsi_x.data());
            compute_diagonal(psi_k, Hpsi_x.data(), vx_diag.data(), ik);
        }
#endif

        double wk = kv_->wk[ik];

        // Accumulate energies
        for (int ib = 0; ib < nbands_; ++ib)
        {
            double n = occ_flat[ik * nbands_ + ib];
            double gn = xc_func_.g(n);
            E_one_ += wk * n * h_one_diag[ib];
            E_hartree_ += wk * n * vh_diag[ib] * 0.5; // factor 1/2 for Hartree
            E_xc_ += wk * gn * vx_diag[ib] * 0.5;     // factor 1/2 for exchange
        }

        // Occupation gradient: dE/dn_ik
        for (int ib = 0; ib < nbands_; ++ib)
        {
            double n = occ_flat[ik * nbands_ + ib];
            double dgn = xc_func_.dg(n);
            grad_occ[ik * nbands_ + ib] = wk * (h_one_diag[ib] + vh_diag[ib])
                                          + wk * dgn * vx_diag[ib];
        }

        // Orbital gradient: dE/dC*(ik) = wk * [n * (H_one + V_H) * C + g(n) * H_exx * C]
        for (int ib_local = 0; ib_local < nb_local; ++ib_local)
        {
            int ib_global = ParaV_->local2global_col(ib_local);
            if (ib_global >= nbands_) continue;

            double n = occ_flat[ik * nbands_ + ib_global];
            double gn = xc_func_.g(n);

            TK* grad_ptr = &grad_wfc(ik, ib_local, 0);
            const TK* hone_ptr = &Hpsi_one[ib_local * nbs_local];
            const TK* hh_ptr = &Hpsi_h[ib_local * nbs_local];
            const TK* hx_ptr = &Hpsi_x[ib_local * nbs_local];

            for (int mu = 0; mu < nbs_local; ++mu)
            {
                grad_ptr[mu] = wk * (n * (hone_ptr[mu] + hh_ptr[mu]) + gn * hx_ptr[mu]);
            }
        }
    }

    // Reduce energies across MPI
    Parallel_Reduce::reduce_all(E_one_);
    Parallel_Reduce::reduce_all(E_hartree_);
    Parallel_Reduce::reduce_all(E_xc_);

    E_ewald_ = pelec_->f_en.ewald_energy;
    E_total_ = E_one_ + E_hartree_ + E_xc_ + E_ewald_;

    ModuleBase::timer::tick("RDMFT_EG", "compute");
    return E_total_;
}

template <typename TK, typename TR>
double EnergyGradient<TK, TR>::compute_energy(
    const std::vector<double>& occ_flat,
    const psi::Psi<TK>& wfc)
{
    ModuleBase::timer::tick("RDMFT_EG", "compute_energy");

    assert(ion_initialized_);

    // Rebuild charge and Hartree (these always depend on occupations)
    build_charge(occ_flat, wfc);

    HR_hartree_->set_zero();
    delete op_hartree_;
    op_hartree_ = new Veff_rdmft_local<TK, TR>(
        hsk_hartree_, kv_->kvec_d, pelec_->pot, HR_hartree_, ucell_,
        orb_->cutoffs(), gd_, nspin_, charge_, rho_basis_, vloc_, sf_, "hartree");
    op_hartree_->contributeHR();

    // Populate one-body diagonal cache if not valid
    if (!hone_cache_valid_)
    {
        cached_h_one_diag_.resize(nk_);
        const int nb_local = wfc.get_nbands();
        const int nbs_local = wfc.get_nbasis();
        std::vector<TK> Hpsi_buf(nb_local * nbs_local);

        for (int ik = 0; ik < nk_; ++ik)
        {
            cached_h_one_diag_[ik].assign(nbands_, 0.0);
            const TK* psi_k = &wfc(ik, 0, 0);

            hsk_one_->set_zero_hk();
            op_local_->contributeHk(ik);
            std::fill(Hpsi_buf.begin(), Hpsi_buf.end(), TK(0));
            apply_Hk(hsk_one_->get_hk(), psi_k, Hpsi_buf.data());
            compute_diagonal(psi_k, Hpsi_buf.data(), cached_h_one_diag_[ik].data(), ik);
        }
        hone_cache_valid_ = true;
    }

    // For each k-point, compute Hartree diag and accumulate energies
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();
    std::vector<double> vh_diag(nbands_, 0.0);
    std::vector<double> vx_diag(nbands_, 0.0);
    std::vector<TK> Hpsi_buf(nb_local * nbs_local);

    E_one_ = 0.0;
    E_hartree_ = 0.0;
    E_xc_ = 0.0;

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* psi_k = &wfc(ik, 0, 0);

        // Hartree diag
        hsk_hartree_->set_zero_hk();
        op_hartree_->contributeHk(ik);
        std::fill(Hpsi_buf.begin(), Hpsi_buf.end(), TK(0));
        apply_Hk(hsk_hartree_->get_hk(), psi_k, Hpsi_buf.data());
        std::fill(vh_diag.begin(), vh_diag.end(), 0.0);
        compute_diagonal(psi_k, Hpsi_buf.data(), vh_diag.data(), ik);

        // Exchange diag (RDMFT's private EXX)
        std::fill(vx_diag.begin(), vx_diag.end(), 0.0);
#ifdef __EXX
        if (exx_enabled_ && op_exx_)
        {
            hsk_exx_->set_zero_hk();
            op_exx_->contributeHk(ik);
            std::fill(Hpsi_buf.begin(), Hpsi_buf.end(), TK(0));
            apply_Hk(hsk_exx_->get_hk(), psi_k, Hpsi_buf.data());
            compute_diagonal(psi_k, Hpsi_buf.data(), vx_diag.data(), ik);
        }
#endif

        double wk = kv_->wk[ik];
        for (int ib = 0; ib < nbands_; ++ib)
        {
            double n = occ_flat[ik * nbands_ + ib];
            double gn = xc_func_.g(n);
            E_one_ += wk * n * cached_h_one_diag_[ik][ib];
            E_hartree_ += wk * n * vh_diag[ib] * 0.5;
            E_xc_ += wk * gn * vx_diag[ib] * 0.5;
        }
    }

    Parallel_Reduce::reduce_all(E_one_);
    Parallel_Reduce::reduce_all(E_hartree_);
    Parallel_Reduce::reduce_all(E_xc_);

    E_ewald_ = pelec_->f_en.ewald_energy;
    E_total_ = E_one_ + E_hartree_ + E_xc_ + E_ewald_;

    ModuleBase::timer::tick("RDMFT_EG", "compute_energy");
    return E_total_;
}

// Explicit template instantiations
template class Veff_rdmft_local<double, double>;
template class Veff_rdmft_local<std::complex<double>, double>;
template class Veff_rdmft_local<std::complex<double>, std::complex<double>>;

template class EnergyGradient<double, double>;
template class EnergyGradient<std::complex<double>, double>;
template class EnergyGradient<std::complex<double>, std::complex<double>>;

} // namespace rdmft
