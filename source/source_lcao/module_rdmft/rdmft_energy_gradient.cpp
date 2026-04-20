#include "rdmft_energy_gradient.h"
#include "source_base/timer.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/lapack_connector.h"
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

#include "source_lcao/module_operator_lcao/overlap.h"

// ScaLAPACK Cholesky / triangular-solve prototypes not exposed by
// ScalapackConnector. We only need them inside this TU.
extern "C" {
    // Serial LAPACK/BLAS routines for non-MPI path
    void dtrtri_(const char* uplo, const char* diag, const int* n,
                 double* A, const int* lda, int* info);
    void ztrtri_(const char* uplo, const char* diag, const int* n,
                 std::complex<double>* A, const int* lda, int* info);
    void dtrmm_(const char* side, const char* uplo, const char* trans,
                const char* diag, const int* m, const int* n,
                const double* alpha, const double* A, const int* lda,
                double* B, const int* ldb);
    void ztrmm_(const char* side, const char* uplo, const char* trans,
                const char* diag, const int* m, const int* n,
                const std::complex<double>* alpha, const std::complex<double>* A, const int* lda,
                std::complex<double>* B, const int* ldb);

    // ScaLAPACK (MPI) routines not in scalapack_connector.h
    void pdtrsm_(const char* side, const char* uplo, const char* trans,
                 const char* diag, const int* m, const int* n,
                 const double* alpha,
                 const double* A, const int* ia, const int* ja, const int* descA,
                 double* B, const int* ib, const int* jb, const int* descB);
    void pztrsm_(const char* side, const char* uplo, const char* trans,
                 const char* diag, const int* m, const int* n,
                 const std::complex<double>* alpha,
                 const std::complex<double>* A, const int* ia, const int* ja, const int* descA,
                 std::complex<double>* B, const int* ib, const int* jb, const int* descB);
    // pdtrmm_ and pztrmm_ are already declared (without const) in scalapack_connector.h
    void pdtrtri_(const char* uplo, const char* diag, const int* n,
                  double* A, const int* ia, const int* ja, const int* descA,
                  int* info);
    void pztrtri_(const char* uplo, const char* diag, const int* n,
                  std::complex<double>* A, const int* ia, const int* ja, const int* descA,
                  int* info);
}

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <memory>
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

// |v|^2 for Frobenius accumulation (C++14-friendly; avoids if constexpr / is_same_v)
inline double gram_elem_frob_sq(double v) { return v * v; }
inline double gram_elem_frob_sq(const std::complex<double>& v) { return std::norm(v); }

// Transpose character for real/complex
inline char trans_char(double) { return 'T'; }
inline char trans_char(std::complex<double>) { return 'C'; }

// Non-MPI BLAS helpers: plain column-major dense matrix operations.
// C = alpha * op(A) * op(B) + beta * C  (col-major, lda = leading dimension)
inline void gemm_wrapper(char ta, char tb, int m, int n, int k,
    double alpha, const double* A, int lda,
    const double* B, int ldb,
    double beta, double* C, int ldc)
{
    dgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}
inline void gemm_wrapper(char ta, char tb, int m, int n, int k,
    std::complex<double> alpha, const std::complex<double>* A, int lda,
    const std::complex<double>* B, int ldb,
    std::complex<double> beta, std::complex<double>* C, int ldc)
{
    zgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

// Cholesky factorisation: A = L L^H (lower triangular, in-place)
inline int potrf_lower(double* A, int n)
{
    char uplo = 'L';
    int info = 0;
    dpotrf_(&uplo, &n, A, &n, &info);
    return info;
}
inline int potrf_lower(std::complex<double>* A, int n)
{
    char uplo = 'L';
    int info = 0;
    zpotrf_(&uplo, &n, A, &n, &info);
    return info;
}

// Triangular solve: B <- B * (L^H)^{-1}  (side='R', uplo='L', trans='T' for real/'C' for complex, diag='N')
// i.e. solve B * L^H = B_in for B (in-place). Dimensions: B is m x n.
inline void trsm_right_lower_conjt(int m, int n, double* B, int ldb, const double* L, int ldl)
{
    char side = 'R', uplo = 'L', trans = 'T', diag = 'N';
    double alpha = 1.0;
    dtrsm_(&side, &uplo, &trans, &diag, &m, &n, &alpha, L, &ldl, B, &ldb);
}
inline void trsm_right_lower_conjt(int m, int n, std::complex<double>* B, int ldb,
                                    const std::complex<double>* L, int ldl)
{
    char side = 'R', uplo = 'L', trans = 'C', diag = 'N';
    std::complex<double> alpha = {1.0, 0.0};
    ztrsm_(&side, &uplo, &trans, &diag, &m, &n, &alpha, L, &ldl, B, &ldb);
}

// Upper-triangular Cholesky: A = U^H U (upper triangular, in-place)
inline int potrf_upper(double* A, int n)
{
    char uplo = 'U';
    int info = 0;
    dpotrf_(&uplo, &n, A, &n, &info);
    return info;
}
inline int potrf_upper(std::complex<double>* A, int n)
{
    char uplo = 'U';
    int info = 0;
    zpotrf_(&uplo, &n, A, &n, &info);
    return info;
}

// Triangular inverse: compute U^{-1} in-place for upper triangular U.
inline int trtri_upper(double* U, int n)
{
    char uplo = 'U', diag = 'N';
    int info = 0;
    dtrtri_(&uplo, &diag, &n, U, &n, &info);
    return info;
}
inline int trtri_upper(std::complex<double>* U, int n)
{
    char uplo = 'U', diag = 'N';
    int info = 0;
    ztrtri_(&uplo, &diag, &n, U, &n, &info);
    return info;
}

// Triangular matrix multiply: B <- U * B  (side='L', uplo='U', trans='N', diag='N')
// U is upper triangular n x n, B is n x nrhs, column-major.
inline void trmm_left_upper_notr(int n, int nrhs, const double* U, int ldu,
                                  double* B, int ldb)
{
    char side = 'L', uplo = 'U', trans = 'N', diag = 'N';
    double alpha = 1.0;
    dtrmm_(&side, &uplo, &trans, &diag, &n, &nrhs, &alpha, U, &ldu, B, &ldb);
}
inline void trmm_left_upper_notr(int n, int nrhs, const std::complex<double>* U, int ldu,
                                  std::complex<double>* B, int ldb)
{
    char side = 'L', uplo = 'U', trans = 'N', diag = 'N';
    std::complex<double> alpha = {1.0, 0.0};
    ztrmm_(&side, &uplo, &trans, &diag, &n, &nrhs, &alpha, U, &ldu, B, &ldb);
}

// Triangular solve: B <- U^{-1} * B  (side='L', uplo='U', trans='N')
inline void trsm_left_upper_notr(int n, int nrhs, const double* U, int ldu,
                                  double* B, int ldb)
{
    char side = 'L', uplo = 'U', trans = 'N', diag = 'N';
    double alpha = 1.0;
    dtrsm_(&side, &uplo, &trans, &diag, &n, &nrhs, &alpha, U, &ldu, B, &ldb);
}
inline void trsm_left_upper_notr(int n, int nrhs, const std::complex<double>* U, int ldu,
                                  std::complex<double>* B, int ldb)
{
    char side = 'L', uplo = 'U', trans = 'N', diag = 'N';
    std::complex<double> alpha = {1.0, 0.0};
    ztrsm_(&side, &uplo, &trans, &diag, &n, &nrhs, &alpha, U, &ldu, B, &ldb);
}

// Triangular solve: B <- U^{-H} * B  (side='L', uplo='U', trans='C'/'T')
inline void trsm_left_upper_conjt(int n, int nrhs, const double* U, int ldu,
                                   double* B, int ldb)
{
    char side = 'L', uplo = 'U', trans = 'T', diag = 'N';
    double alpha = 1.0;
    dtrsm_(&side, &uplo, &trans, &diag, &n, &nrhs, &alpha, U, &ldu, B, &ldb);
}
inline void trsm_left_upper_conjt(int n, int nrhs, const std::complex<double>* U, int ldu,
                                   std::complex<double>* B, int ldb)
{
    char side = 'L', uplo = 'U', trans = 'C', diag = 'N';
    std::complex<double> alpha = {1.0, 0.0};
    ztrsm_(&side, &uplo, &trans, &diag, &n, &nrhs, &alpha, U, &ldu, B, &ldb);
}

} // namespace detail

namespace {
/// Threshold below which multiplicative weights (n, g(n)) are treated as zero.
constexpr double rdmft_occ_weight_eps = 1e-12;
/// Skip multiplicative contributions when |weight| is negligible (e.g. |n| < eps).
inline bool rdmft_skip_occ_weight(double x)
{
    return std::fabs(x) < rdmft_occ_weight_eps;
}
} // namespace

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

namespace {
template <typename TR_>
void local_build_HR_gint(const std::string& potential,
                          int nspin,
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
    if (potential == "hartree")
    {
        ModuleBase::matrix v(nspin, charge->nrxx);
        elecstate::PotHartree potH(rho_basis);
        potH.cal_v_eff(charge, ucell, v);
        for (int is = 0; is < nspin; ++is)
        {
            vr_eff = &v(is, 0);
            ModuleGint::cal_gint_vl(vr_eff, hR);
        }
    }
    else if (potential == "local")
    {
        double vlocal_of_0 = 0.0;
        ModuleBase::matrix v(1, charge->nrxx);
        elecstate::PotLocal potL(vloc, sf, rho_basis, vlocal_of_0);
        potL.cal_fixed_v(&v(0, 0));
        vr_eff = &v(0, 0);
        ModuleGint::cal_gint_vl(vr_eff, hR);
    }
    else if (potential == "xc")
    {
        ModuleBase::matrix vofk = *vloc;
        vofk.zero_out();
        ModuleBase::matrix v(nspin, charge->nrxx);
        elecstate::PotXC potXC(rho_basis, etxc, vtxc, &vofk);
        potXC.cal_v_eff(charge, ucell, v);
        for (int is = 0; is < nspin; ++is)
        {
            vr_eff = &v(is, 0);
            ModuleGint::cal_gint_vl(vr_eff, hR);
        }
    }
}
} // anonymous namespace

template <>
void Veff_rdmft_local<std::complex<double>, double>::contributeHR()
{
    local_build_HR_gint(potential_, nspin_, charge_, ucell_, rho_basis_,
                         vloc_, sf_, etxc_, vtxc_, this->hR);
}

template <>
void Veff_rdmft_local<double, double>::contributeHR()
{
    local_build_HR_gint(potential_, nspin_, charge_, ucell_, rho_basis_,
                         vloc_, sf_, etxc_, vtxc_, this->hR);
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
    delete SR_;
    delete hsk_one_;
    delete hsk_hartree_;
    delete hsk_exx_;
    delete hsk_overlap_;
    delete op_ekinetic_;
    delete op_nonlocal_;
    delete op_local_;
    delete op_hartree_;
    delete op_exx_;
    delete op_overlap_;
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
    SR_ = new hamilt::HContainer<TR>(*ucell_, ParaV_);

    hsk_one_ = new hamilt::HS_Matrix_K<TK>(ParaV_, true);
    hsk_hartree_ = new hamilt::HS_Matrix_K<TK>(ParaV_, true);
    hsk_exx_ = new hamilt::HS_Matrix_K<TK>(ParaV_, true);
    // The overlap operator writes S(k) into hsk->sk via hsk->get_sk(),
    // so we must allocate sk (second arg false -> no_s==false)
    hsk_overlap_ = new hamilt::HS_Matrix_K<TK>(ParaV_, false);

    if (PARAM.inp.gamma_only)
    {
        HR_one_->fix_gamma();
        HR_hartree_->fix_gamma();
        HR_exx_->fix_gamma();
        SR_->fix_gamma();
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

        // Propagate input RI parameters (input_conv skips this block when
        // cal_exx is off, but RDMFT needs a fully-initialised info_ri to run
        // its own EXX evaluator).
        GlobalC::exx_info.info_ri.real_number = std::stoi(PARAM.inp.exx_real_number);
        GlobalC::exx_info.info_ri.pca_threshold = PARAM.inp.exx_pca_threshold;
        GlobalC::exx_info.info_ri.C_threshold   = PARAM.inp.exx_c_threshold;
        GlobalC::exx_info.info_ri.V_threshold   = PARAM.inp.exx_v_threshold;
        GlobalC::exx_info.info_ri.dm_threshold  = PARAM.inp.exx_dm_threshold;
        GlobalC::exx_info.info_ri.ccp_rmesh_times = std::stod(PARAM.inp.exx_ccp_rmesh_times);
        GlobalC::exx_info.info_ri.exx_symmetry_realspace = PARAM.inp.exx_symmetry_realspace;
        GlobalC::exx_info.info_ri.Cs_inv_thr = PARAM.inp.exx_cs_inv_thr;
        GlobalC::exx_info.info_ri.shrink_abfs_pca_thr = PARAM.inp.shrink_abfs_pca_thr;
        GlobalC::exx_info.info_ri.shrink_LU_inv_thr = PARAM.inp.shrink_LU_inv_thr;
        GlobalC::exx_info.info_ri.coul_moment = PARAM.inp.exx_coul_moment;
        GlobalC::exx_info.info_ri.rotate_abfs = PARAM.inp.exx_rotate_abfs;
        GlobalC::exx_info.info_ri.multip_moments_threshold = PARAM.inp.exx_multip_moments_threshold;
        GlobalC::exx_info.info_opt_abfs.pca_threshold = PARAM.inp.exx_pca_threshold;
        GlobalC::exx_info.info_opt_abfs.abfs_Lmax = PARAM.inp.exx_opt_orb_lmax;
        GlobalC::exx_info.info_opt_abfs.ecut_exx = PARAM.inp.exx_opt_orb_ecut;
        GlobalC::exx_info.info_opt_abfs.tolerence = PARAM.inp.exx_opt_orb_tolerence;
        GlobalC::exx_info.info_global.hybrid_alpha = 1.0; // full Fock for RDMFT

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

    // Overlap operator (builds SR internally)
    delete op_overlap_;
    op_overlap_ = new hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>(
        hsk_overlap_, kv_->kvec_d, SR_, SR_, ucell_,
        orb_->cutoffs(), gd_, two_center_bundle_->overlap_orb.get());
    op_overlap_->contributeHR();

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

    std::unique_ptr<elecstate::DensityMatrix<TK, double>> DM;
    if (PARAM.inp.gamma_only)
    {
        DM.reset(new elecstate::DensityMatrix<TK, double>(ParaV_, nspin_));
    }
    else
    {
        DM.reset(new elecstate::DensityMatrix<TK, double>(ParaV_, nspin_, kv_->kvec_d, nk_));
    }
    elecstate::cal_dm_psi(ParaV_, wg, wfc, *DM);
    DM->init_DMR(gd_, ucell_);
    DM->cal_DMR();
    for (int is = 0; is < nspin_; ++is)
    {
        ModuleBase::GlobalFunc::ZEROS(charge_->rho[is], charge_->nrxx);
    }
    ModuleGint::cal_gint_rho(DM->get_DMR_vector(), nspin_, charge_->rho);

    // NOTE: Do NOT call charge_->renormalize_rho() here.
    // In RDMFT, the density must be exactly rho = sum_k w_k sum_i n_ik |phi_ik|^2
    // for the analytic occupation gradient to be consistent with finite differences
    // (perturbing n_ik by epsilon would otherwise be undone by renormalization).
    // The electron-number constraint is handled independently at the occupation level
    // (augmented Lagrangian / projection). Renormalizing here would hide constraint
    // violations from the energy and break gradient consistency.

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

    // Build wk_g matrix: wk[ik] * g(n(ik, ib)). Bands with g(n)=0 do not enter the XC density matrix.
    ModuleBase::matrix wk_g(nk_, nbands_);
    for (int ik = 0; ik < nk_; ++ik)
    {
        for (int ib = 0; ib < nbands_; ++ib)
        {
            const double gn = xc_func_.g(occ_flat[ik * nbands_ + ib]);
            wk_g(ik, ib) = rdmft_skip_occ_weight(gn) ? 0.0 : kv_->wk[ik] * gn;
        }
    }

    // DM_XC(k) = sum_i wk*g(n_i) * conj(C_i) * C_i^T, built directly through
    // ABACUS' cal_dm_psi with wk_g as the per-band weight.
    std::unique_ptr<elecstate::DensityMatrix<TK, double>> DM_xc;
    if (PARAM.inp.gamma_only)
    {
        DM_xc.reset(new elecstate::DensityMatrix<TK, double>(ParaV_, nspin_));
    }
    else
    {
        DM_xc.reset(new elecstate::DensityMatrix<TK, double>(ParaV_, nspin_, kv_->kvec_d, nk_));
    }
    elecstate::cal_dm_psi(ParaV_, wk_g, wfc, *DM_xc);
    for (int ik = 0; ik < nk_; ++ik)
    {
        TK* dmk = DM_xc->get_DMK_pointer(ik);
        std::copy(dmk, dmk + ParaV_->nloc, DM_XC[ik].begin());
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

    // Compute Eij = psi^H * Hpsi in the 2D BLACS grid
    std::vector<TK> Eij(para_Eij_.get_row_size() * para_Eij_.get_col_size(), TK(0));

    detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, psi_k, 1, 1, ParaV_->desc_wfc,
        Hpsi_k, 1, 1, ParaV_->desc_wfc,
        zero, Eij.data(), 1, 1, para_Eij_.desc);

    // Extract diagonal. Only the one process that owns the (ig,ig) entry of the
    // 2D BLACS grid writes a non-zero value; all other processes write 0. An
    // Allreduce-sum below brings the full diagonal to every process so that all
    // downstream code (energy accumulation, gradient formulas, line-search) can
    // treat 'diag' as a globally-consistent replicated array.
    std::fill(diag, diag + nbands, 0.0);
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
    Parallel_Reduce::reduce_all(diag, nbands);
#endif
}

template <typename TK, typename TR>
const TK* EnergyGradient<TK, TR>::get_SK(int ik)
{
    if (op_overlap_ == nullptr || hsk_overlap_ == nullptr) return nullptr;
    // NOTE: do NOT call hsk_overlap_->set_zero_sk() here.
    // contributeHk() internally calls set_zero_sk() before recomputing, so the
    // external zero is redundant. More importantly, for the gamma-only real case
    // (TK=double), contributeHk() caches the result and returns early when the
    // k-vector hasn't changed. Calling set_zero_sk() before this early return
    // would leave the sk buffer zeroed, causing all subsequent get_SK() calls
    // to return a zero S-matrix and silently breaking orbital orthonormalization.
    op_overlap_->contributeHk(ik);
    return hsk_overlap_->get_sk();
}

namespace {
inline double real_of_conj_prod(double a, double b) { return a * b; }
inline double real_of_conj_prod(std::complex<double> a, std::complex<double> b)
{
    // Re(conj(a) * b)
    return a.real() * b.real() + a.imag() * b.imag();
}

} // namespace

template <typename TK, typename TR>
double EnergyGradient<TK, TR>::s_inner_product(
    const psi::Psi<TK>& X,
    const psi::Psi<TK>& Y)
{
    double result = 0.0;
    const int nb_local = Y.get_nbands();
    const int nbs_local = Y.get_nbasis();

    // In X-variable mode, X lives in the S^{1/2}-transformed space where the
    // metric is Euclidean (S = I), so we skip the S multiplication entirely.
    const bool skip_S = use_X_variable_;

#ifdef __MPI
    const int nbasis = ParaV_->desc[2];
    const int nbands = ParaV_->desc_wfc[3];
    const TK one = TK(1.0);
    const TK zero = TK(0.0);

    const std::int64_t sy_alloc_ip
        = std::max<std::int64_t>(static_cast<std::int64_t>(nb_local * nbs_local), 1);
    std::vector<TK> SY(static_cast<size_t>(sy_alloc_ip), TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* Xk = &X(ik, 0, 0);
        const TK* Yk = &Y(ik, 0, 0);
        const TK* SK = skip_S ? nullptr : get_SK(ik);

        if (SK != nullptr)
        {
            std::fill(SY.begin(), SY.end(), TK(0));
            detail::pgemm_wrapper('N', 'N', nbasis, nbands, nbasis,
                one, SK, 1, 1, ParaV_->desc,
                Yk, 1, 1, ParaV_->desc_wfc,
                zero, SY.data(), 1, 1, ParaV_->desc_wfc);
            for (int i = 0; i < nb_local * nbs_local; ++i)
                result += real_of_conj_prod(Xk[i], SY[i]);
        }
        else
        {
            for (int i = 0; i < nb_local * nbs_local; ++i)
                result += real_of_conj_prod(Xk[i], Yk[i]);
        }
    }
    Parallel_Reduce::reduce_all(result);
#else
    // Non-MPI: nbs_local == nbasis (full matrix local). Use plain BLAS when
    // S is available; fall back to Euclidean inner product when S == I.
    const int nbasis = nbs_local;
    const int nbands = nb_local;
    std::vector<TK> SY(nbasis * nbands, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* Xk = &X(ik, 0, 0);
        const TK* Yk = &Y(ik, 0, 0);
        const TK* SK = skip_S ? nullptr : get_SK(ik);

        if (SK != nullptr)
        {
            // SY = S * Y
            const TK one = TK(1.0);
            const TK zero = TK(0.0);
            detail::gemm_wrapper('N', 'N', nbasis, nbands, nbasis,
                one, SK, nbasis, Yk, nbasis, zero, SY.data(), nbasis);
            for (int i = 0; i < nbasis * nbands; ++i)
                result += real_of_conj_prod(Xk[i], SY[i]);
        }
        else
        {
            for (int i = 0; i < nbasis * nbands; ++i)
                result += real_of_conj_prod(Xk[i], Yk[i]);
        }
    }
#endif
    return result;
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::stiefel_gram_residual_frobenius_per_k(
    const psi::Psi<TK>& wfc,
    std::vector<double>& frob_per_ik)
{
    frob_per_ik.assign(nk_, 0.0);
    const bool skip_S = use_X_variable_;
    const TK one = TK(1.0);
    const TK zero = TK(0.0);
    const char tc = detail::trans_char(TK());

#ifdef __MPI
    const int nbasis = ParaV_->desc[2];
    const int nbands = ParaV_->desc_wfc[3];
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();
    const int nrow = para_Eij_.get_row_size();
    const int ncol = para_Eij_.get_col_size();
    const int eij_nloc = nrow * ncol;

    // ScaLAPACK must not receive a null matrix pointer when local dimensions are zero;
    // std::vector<T>(0).data() may be nullptr and breaks p*gemm / heap on some ranks.
    const std::int64_t sy_alloc
        = std::max<std::int64_t>(static_cast<std::int64_t>(nb_local * nbs_local), 1);
    const std::int64_t m_alloc = std::max<std::int64_t>(static_cast<std::int64_t>(eij_nloc), 1);
    std::vector<TK> SY(static_cast<size_t>(sy_alloc), TK(0));
    std::vector<TK> M(static_cast<size_t>(m_alloc), TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* C = &wfc(ik, 0, 0);
        const TK* SK = skip_S ? nullptr : get_SK(ik);

        std::fill(SY.begin(), SY.end(), TK(0));
        if (SK != nullptr)
        {
            detail::pgemm_wrapper('N', 'N', nbasis, nbands, nbasis,
                one, SK, 1, 1, ParaV_->desc,
                C, 1, 1, ParaV_->desc_wfc,
                zero, SY.data(), 1, 1, ParaV_->desc_wfc);
        }
        else
        {
            for (int i = 0; i < nb_local * nbs_local; ++i)
            {
                SY[i] = C[i];
            }
        }

        std::fill(M.begin(), M.end(), TK(0));
        detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, C, 1, 1, ParaV_->desc_wfc,
            SY.data(), 1, 1, ParaV_->desc_wfc,
            zero, M.data(), 1, 1, para_Eij_.desc);

        double frob_sq = 0.0;
        double tr = 0.0;
        for (int i = 0; i < nrow; ++i)
        {
            const int ig = para_Eij_.local2global_row(i);
            for (int j = 0; j < ncol; ++j)
            {
                const int jg = para_Eij_.local2global_col(j);
                const TK val = M[i + j * nrow];
                frob_sq += detail::gram_elem_frob_sq(val);
                if (ig == jg)
                {
                    tr += detail::real_of(val);
                }
            }
        }
        Parallel_Reduce::reduce_all(frob_sq);
        Parallel_Reduce::reduce_all(tr);

        const double resid_sq = frob_sq - 2.0 * tr + static_cast<double>(nbands);
        frob_per_ik[ik] = std::sqrt(std::max(0.0, resid_sq));
    }
#else
    const int nbasis = wfc.get_nbasis();
    const int nbands = wfc.get_nbands();
    std::vector<TK> SY(nbasis * nbands, TK(0));
    std::vector<TK> M(nbands * nbands, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* C = &wfc(ik, 0, 0);
        const TK* SK = skip_S ? nullptr : get_SK(ik);

        if (SK != nullptr)
        {
            detail::gemm_wrapper('N', 'N', nbasis, nbands, nbasis,
                one, SK, nbasis, C, nbasis, zero, SY.data(), nbasis);
        }
        else
        {
            for (int i = 0; i < nbasis * nbands; ++i)
            {
                SY[i] = C[i];
            }
        }

        detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, C, nbasis, SY.data(), nbasis, zero, M.data(), nbands);

        double frob_sq = 0.0;
        double tr = 0.0;
        for (int j = 0; j < nbands; ++j)
        {
            for (int i = 0; i < nbands; ++i)
            {
                const TK val = M[i + j * nbands];
                frob_sq += detail::gram_elem_frob_sq(val);
                if (i == j)
                {
                    tr += detail::real_of(val);
                }
            }
        }
        const double resid_sq = frob_sq - 2.0 * tr + static_cast<double>(nbands);
        frob_per_ik[ik] = std::sqrt(std::max(0.0, resid_sq));
    }
#endif
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::project_orbital_gradient(
    const psi::Psi<TK>& wfc,
    psi::Psi<TK>& grad_wfc)
{
    // Compute the Riemannian (tangent-space) gradient on the generalised
    // Stiefel manifold  { C : C^H S C = I } with the canonical (trace) metric:
    //     proj_C(G) = G - C * sym(C^H S G)
    // where sym(M) = 0.5*(M + M^H).
    // Derivation: we seek G_R = G - C*K such that C^H S G_R is skew-Hermitian.
    //   C^H S G_R = C^H S G - C^H S C * K = C^H S G - K  (since C^H S C = I)
    // Skew-Hermitian requirement: K + K^H = C^H S G + G^H S C
    //   => K = sym(C^H S G)
    // So: G_R = G - C * sym(C^H S G).
    // Note: when S = I (or when using X-variable mode where X = U C and
    // X^H X = I) the formula reduces to  G_R = G - X * sym(X^H G),
    // which is the standard Stiefel projection for the canonical metric.
    // The resulting G_R vanishes at any S-orthonormal critical point of E,
    // so the orbital inner loop converges immediately there.
    const bool skip_S = use_X_variable_;
    const TK one = TK(1.0);
    const TK zero = TK(0.0);
    const TK neg_one = TK(-1.0);
    char tc = detail::trans_char(TK());
#ifdef __MPI
    const int nbasis = ParaV_->desc[2];
    const int nbands = ParaV_->desc_wfc[3];
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();

    const int eij_nloc = para_Eij_.get_row_size() * para_Eij_.get_col_size();
    const std::int64_t sc_alloc
        = std::max<std::int64_t>(static_cast<std::int64_t>(nb_local * nbs_local), 1);
    const std::int64_t eij_alloc = std::max<std::int64_t>(static_cast<std::int64_t>(eij_nloc), 1);

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* psi_k = &wfc(ik, 0, 0);
        TK* g_k = &grad_wfc(ik, 0, 0);

        // SC = S * C  (if S available and not in X-mode), otherwise SC == C
        std::vector<TK> SC(static_cast<size_t>(sc_alloc), TK(0));
        const TK* SK = skip_S ? nullptr : get_SK(ik);
        if (SK != nullptr)
        {
            detail::pgemm_wrapper('N', 'N', nbasis, nbands, nbasis,
                one, SK, 1, 1, ParaV_->desc,
                psi_k, 1, 1, ParaV_->desc_wfc,
                zero, SC.data(), 1, 1, ParaV_->desc_wfc);
        }
        else
        {
            const int npsi = nbasis * wfc.get_nbands();
            for (int i = 0; i < npsi; ++i) SC[i] = psi_k[i];
        }

        // A = (SC)^H G = C^H S G  (nbands x nbands)
        std::vector<TK> A(static_cast<size_t>(eij_alloc), TK(0));
        detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, SC.data(), 1, 1, ParaV_->desc_wfc,
            g_k, 1, 1, ParaV_->desc_wfc,
            zero, A.data(), 1, 1, para_Eij_.desc);

        // B = G^H (SC) = G^H S C, then symmetric part  A <- 0.5 (A + B) = sym(C^H S G)
        // This is the correct Riemannian symmetrisation: sym(M) = 0.5 (M + M^H).
        std::vector<TK> B(static_cast<size_t>(eij_alloc), TK(0));
        detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, g_k, 1, 1, ParaV_->desc_wfc,
            SC.data(), 1, 1, ParaV_->desc_wfc,
            zero, B.data(), 1, 1, para_Eij_.desc);
        for (int i = 0; i < eij_nloc; ++i) A[i] = TK(0.5) * (A[i] + B[i]);

        // G <- G - C * sym(C^H S G)  (use psi_k = C, not SC)
        detail::pgemm_wrapper('N', 'N', nbasis, nbands, nbands,
            neg_one, psi_k, 1, 1, ParaV_->desc_wfc,
            A.data(), 1, 1, para_Eij_.desc,
            one, g_k, 1, 1, ParaV_->desc_wfc);
    }
#else
    // Non-MPI: nbs_local == nbasis. Use plain BLAS.
    const int nbasis = wfc.get_nbasis();
    const int nbands = wfc.get_nbands();

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* psi_k = &wfc(ik, 0, 0);
        TK* g_k = &grad_wfc(ik, 0, 0);

        // SC = S * C  (if S available and not in X-mode), otherwise SC == C
        std::vector<TK> SC(nbasis * nbands, TK(0));
        const TK* SK = skip_S ? nullptr : get_SK(ik);
        if (SK != nullptr)
        {
            detail::gemm_wrapper('N', 'N', nbasis, nbands, nbasis,
                one, SK, nbasis, psi_k, nbasis, zero, SC.data(), nbasis);
        }
        else
        {
            for (int i = 0; i < nbasis * nbands; ++i) SC[i] = psi_k[i];
        }

        // A = (SC)^H G = C^H S G  (nbands x nbands)
        std::vector<TK> A(nbands * nbands, TK(0));
        detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, SC.data(), nbasis, g_k, nbasis, zero, A.data(), nbands);

        // B = G^H (SC) = G^H S C, then symmetric part  A <- 0.5 (A + B) = sym(C^H S G)
        std::vector<TK> B(nbands * nbands, TK(0));
        detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, g_k, nbasis, SC.data(), nbasis, zero, B.data(), nbands);
        for (int i = 0; i < nbands * nbands; ++i) A[i] = TK(0.5) * (A[i] + B[i]);

        // G <- G - C * sym(C^H S G)
        detail::gemm_wrapper('N', 'N', nbasis, nbands, nbands,
            neg_one, psi_k, nbasis, A.data(), nbands, one, g_k, nbasis);
    }
#endif
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::retract_orbitals(
    psi::Psi<TK>& wfc,
    const psi::Psi<TK>& grad_wfc,
    double alpha)
{
    // One retraction step on the generalised Stiefel manifold:
    //   Y = C - alpha * G
    //   C_new = Y * (Y^H S Y)^{-1/2}
    // We approximate the inverse-square-root via Cholesky of M = Y^H S Y:
    //     M = L L^H   =>   Y_ortho = Y * L^{-H}
    // This is the "Cholesky QR" S-orthonormalisation, much cheaper than a
    // Hermitian eigendecomposition and perfectly adequate as long as alpha
    // is chosen small enough that M stays well-conditioned (which the outer
    // line search guarantees).
    //
    // In X-variable mode, S = I in the transformed space, so M = Y^H Y
    // and the S multiplication is skipped.
    const bool skip_S = use_X_variable_;
    const TK one = TK(1.0);
    const TK zero = TK(0.0);
    char tc = detail::trans_char(TK());
#ifdef __MPI
    const int nbasis = ParaV_->desc[2];
    const int nbands = ParaV_->desc_wfc[3];
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();
    const int eij_nloc = para_Eij_.get_row_size() * para_Eij_.get_col_size();
    const std::int64_t sy_alloc_retract
        = std::max<std::int64_t>(static_cast<std::int64_t>(nb_local * nbs_local), 1);
    const std::int64_t m_alloc_retract
        = std::max<std::int64_t>(static_cast<std::int64_t>(eij_nloc), 1);

    for (int ik = 0; ik < nk_; ++ik)
    {
        TK* C = &wfc(ik, 0, 0);
        const TK* G = &grad_wfc(ik, 0, 0);

        // Y = C - alpha * G (in-place on C)
        for (int i = 0; i < nb_local * nbs_local; ++i)
            C[i] -= TK(alpha) * G[i];

        // SY = S * Y (or SY = Y when skip_S or no overlap)
        std::vector<TK> SY(static_cast<size_t>(sy_alloc_retract), TK(0));
        const TK* SK = skip_S ? nullptr : get_SK(ik);
        if (SK != nullptr)
        {
            detail::pgemm_wrapper('N', 'N', nbasis, nbands, nbasis,
                one, SK, 1, 1, ParaV_->desc,
                C, 1, 1, ParaV_->desc_wfc,
                zero, SY.data(), 1, 1, ParaV_->desc_wfc);
        }
        else
        {
            for (int i = 0; i < nb_local * nbs_local; ++i) SY[i] = C[i];
        }

        // M = Y^H S Y
        std::vector<TK> M(static_cast<size_t>(m_alloc_retract), TK(0));
        detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, C, 1, 1, ParaV_->desc_wfc,
            SY.data(), 1, 1, ParaV_->desc_wfc,
            zero, M.data(), 1, 1, para_Eij_.desc);

        // Cholesky: M = L L^H, then Y <- Y * L^{-H}
        {
            int info = 0;
            int one_int = 1;
            char uplo = 'L';
            if (std::is_same<TK, double>::value)
            {
                pdpotrf_(&uplo, const_cast<int*>(&nbands),
                    reinterpret_cast<double*>(M.data()),
                    &one_int, &one_int, const_cast<int*>(para_Eij_.desc), &info);
            }
            else
            {
                pzpotrf_(&uplo, const_cast<int*>(&nbands),
                    reinterpret_cast<std::complex<double>*>(M.data()),
                    &one_int, &one_int, const_cast<int*>(para_Eij_.desc), &info);
            }
            if (info != 0)
            {
                // Cholesky failed; fall back to no re-orthonormalisation for
                // this step. The outer line search should reject it.
                continue;
            }

            // Y <- Y * (L^H)^{-1}: solve X * L^H = Y  (side = 'R', trans = 'C')
            char side = 'R';
            char trans = detail::trans_char(TK());
            char diag = 'N';
            if (std::is_same<TK, double>::value)
            {
                double alpha_r = 1.0;
                pdtrsm_(&side, &uplo, &trans, &diag, &nbasis, &nbands,
                    &alpha_r,
                    reinterpret_cast<double*>(M.data()),
                    &one_int, &one_int, para_Eij_.desc,
                    reinterpret_cast<double*>(C),
                    &one_int, &one_int,
                    const_cast<int*>(ParaV_->desc_wfc));
            }
            else
            {
                std::complex<double> alpha_c = {1.0, 0.0};
                pztrsm_(&side, &uplo, &trans, &diag, &nbasis, &nbands,
                    &alpha_c,
                    reinterpret_cast<std::complex<double>*>(M.data()),
                    &one_int, &one_int, para_Eij_.desc,
                    reinterpret_cast<std::complex<double>*>(C),
                    &one_int, &one_int,
                    const_cast<int*>(ParaV_->desc_wfc));
            }
        }
    }
#else
    // Non-MPI: nbs_local == nbasis. Use plain BLAS/LAPACK.
    const int nbasis = wfc.get_nbasis();
    const int nbands = wfc.get_nbands();

    for (int ik = 0; ik < nk_; ++ik)
    {
        TK* C = &wfc(ik, 0, 0);
        const TK* G = &grad_wfc(ik, 0, 0);

        // Y = C - alpha * G (in-place on C)
        for (int i = 0; i < nbasis * nbands; ++i)
            C[i] -= TK(alpha) * G[i];

        // SY = S * Y (or SY = Y when skip_S or no overlap)
        std::vector<TK> SY(nbasis * nbands, TK(0));
        const TK* SK = skip_S ? nullptr : get_SK(ik);
        if (SK != nullptr)
        {
            detail::gemm_wrapper('N', 'N', nbasis, nbands, nbasis,
                one, SK, nbasis, C, nbasis, zero, SY.data(), nbasis);
        }
        else
        {
            for (int i = 0; i < nbasis * nbands; ++i) SY[i] = C[i];
        }

        // M = Y^H S Y  (nbands x nbands)
        std::vector<TK> M(nbands * nbands, TK(0));
        detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, C, nbasis, SY.data(), nbasis, zero, M.data(), nbands);

        // Cholesky: M = L L^H
        const int info = detail::potrf_lower(M.data(), nbands);
        if (info != 0)
        {
            // Cholesky failed; skip re-orthonormalisation for this step.
            // The outer line search should reject this trial point.
            GlobalV::ofs_running << "WARNING: Cholesky factorisation of Y^H S Y failed"
                << " at ik=" << ik << " (info=" << info
                << "); skipping re-orthonormalisation." << std::endl;
            continue;
        }

        // C <- C * (L^H)^{-1}
        detail::trsm_right_lower_conjt(nbasis, nbands, C, nbasis, M.data(), nbands);
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
    ModuleBase::timer::start("RDMFT_EG", "compute");

    assert(ion_initialized_);

    // When in X-variable mode, wfc contains X_k. Convert to C_k for
    // energy evaluation.  We use a non-const copy so that wfc_X_to_C
    // can operate in-place on the temporary.
    psi::Psi<TK> wfc_C;
    const psi::Psi<TK>* wfc_ptr = &wfc;
    if (use_X_variable_ && cholesky_precomputed_)
    {
        wfc_C = wfc;  // deep copy
        wfc_X_to_C(wfc_C);
        wfc_ptr = &wfc_C;
    }
    const psi::Psi<TK>& wfc_eval = *wfc_ptr;

    build_charge(occ_flat, wfc_eval);

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
        build_DM_xc(occ_flat, wfc_eval, DM_XC);

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
    const int nb_local = wfc_eval.get_nbands();
    const int nbs_local = wfc_eval.get_nbasis();

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
        const TK* psi_k = &wfc_eval(ik, 0, 0);

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
        // Use RI_2D_Comm::add_Hexx directly because op_exx_->contributeHk() is
        // short-circuited for non-hybrid KS functionals (it checks the global
        // XC func_type and returns early when it is not 4/5). RDMFT runs on top
        // of an LDA/GGA KS and still needs Fock exchange from its modified DM.
        std::fill(Hpsi_x.begin(), Hpsi_x.end(), TK(0));
        std::fill(vx_diag.begin(), vx_diag.end(), 0.0);
#ifdef __EXX
        if (exx_enabled_)
        {
            hsk_exx_->set_zero_hk();
            if (GlobalC::exx_info.info_ri.real_number)
            {
                RI_2D_Comm::add_Hexx(*ucell_, *kv_, ik,
                    1.0, exx_lri_d_->Hexxs, *ParaV_, hsk_exx_->get_hk());
            }
            else
            {
                RI_2D_Comm::add_Hexx(*ucell_, *kv_, ik,
                    1.0, exx_lri_c_->Hexxs, *ParaV_, hsk_exx_->get_hk());
            }
            apply_Hk(hsk_exx_->get_hk(), psi_k, Hpsi_x.data());
            compute_diagonal(psi_k, Hpsi_x.data(), vx_diag.data(), ik);
        }
#endif

        double wk = kv_->wk[ik];

        // Accumulate energies (skip terms that are identically zero when n=0 or g(n)=0)
        for (int ib = 0; ib < nbands_; ++ib)
        {
            const double n = occ_flat[ik * nbands_ + ib];
            const double gn = xc_func_.g(n);
            if (!rdmft_skip_occ_weight(n))
            {
                E_one_ += wk * n * h_one_diag[ib];
                E_hartree_ += wk * n * vh_diag[ib] * 0.5; // factor 1/2 for Hartree
            }
            if (!rdmft_skip_occ_weight(gn))
            {
                E_xc_ += wk * gn * vx_diag[ib] * 0.5; // factor 1/2 for exchange
            }
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
        // Omit terms when n=0 or g(n)=0 so empty / inactive orbitals do not contribute.
        for (int ib_local = 0; ib_local < nb_local; ++ib_local)
        {
            int ib_global = ParaV_->local2global_col(ib_local);
            if (ib_global >= nbands_) continue;

            const double n = occ_flat[ik * nbands_ + ib_global];
            const double gn = xc_func_.g(n);
            const bool use_one_hart = !rdmft_skip_occ_weight(n);
            const bool use_exx = !rdmft_skip_occ_weight(gn);
            if (!use_one_hart && !use_exx)
            {
                continue;
            }

            TK* grad_ptr = &grad_wfc(ik, ib_local, 0);
            const TK* hone_ptr = &Hpsi_one[ib_local * nbs_local];
            const TK* hh_ptr = &Hpsi_h[ib_local * nbs_local];
            const TK* hx_ptr = &Hpsi_x[ib_local * nbs_local];

            for (int mu = 0; mu < nbs_local; ++mu)
            {
                TK acc = TK(0);
                if (use_one_hart)
                {
                    acc += TK(n) * (hone_ptr[mu] + hh_ptr[mu]);
                }
                if (use_exx)
                {
                    acc += TK(gn) * hx_ptr[mu];
                }
                grad_ptr[mu] = wk * acc;
            }
        }
    }

    // h_one_diag / vh_diag / vx_diag are already globally replicated by
    // compute_diagonal(); the per-rank accumulations above all produced the
    // same value, so NO MPI reduction of E_* is needed (that would double-
    // count). The same is true for grad_occ.

    E_ewald_ = pelec_->f_en.ewald_energy;
    E_total_ = E_one_ + E_hartree_ + E_xc_ + E_ewald_;

    // When in X-variable mode, transform the orbital gradient from C-space
    // to X-space: G_X = U^{-H} G_C.
    if (use_X_variable_ && cholesky_precomputed_)
    {
        grad_C_to_X(grad_wfc);
    }

    ModuleBase::timer::end("RDMFT_EG", "compute");
    return E_total_;
}

template <typename TK, typename TR>
double EnergyGradient<TK, TR>::compute_energy(
    const std::vector<double>& occ_flat,
    const psi::Psi<TK>& wfc)
{
    ModuleBase::timer::start("RDMFT_EG", "compute_energy");

    assert(ion_initialized_);

    // When in X-variable mode, convert X -> C for energy evaluation.
    psi::Psi<TK> wfc_C;
    const psi::Psi<TK>* wfc_ptr = &wfc;
    if (use_X_variable_ && cholesky_precomputed_)
    {
        wfc_C = wfc;
        wfc_X_to_C(wfc_C);
        wfc_ptr = &wfc_C;
    }
    const psi::Psi<TK>& wfc_eval = *wfc_ptr;

    // Rebuild charge and Hartree (these always depend on occupations)
    build_charge(occ_flat, wfc_eval);

    HR_hartree_->set_zero();
    delete op_hartree_;
    op_hartree_ = new Veff_rdmft_local<TK, TR>(
        hsk_hartree_, kv_->kvec_d, pelec_->pot, HR_hartree_, ucell_,
        orb_->cutoffs(), gd_, nspin_, charge_, rho_basis_, vloc_, sf_, "hartree");
    op_hartree_->contributeHR();

    // Exchange must also be rebuilt because the modified DM gamma_xc depends
    // on occupations (and on orbitals, but those are fixed for compute_energy
    // use cases). Without this, line-search / finite-difference calls would
    // use a stale H_exx from the last compute() call.
#ifdef __EXX
    if (exx_enabled_)
    {
        HR_exx_->set_zero();
        std::vector<std::vector<TK>> DM_XC;
        build_DM_xc(occ_flat, wfc_eval, DM_XC);

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

    // Populate one-body diagonal cache if not valid
    if (!hone_cache_valid_)
    {
        cached_h_one_diag_.resize(nk_);
        const int nb_local = wfc_eval.get_nbands();
        const int nbs_local = wfc_eval.get_nbasis();
        std::vector<TK> Hpsi_buf(nb_local * nbs_local);

        for (int ik = 0; ik < nk_; ++ik)
        {
            cached_h_one_diag_[ik].assign(nbands_, 0.0);
            const TK* psi_k = &wfc_eval(ik, 0, 0);

            hsk_one_->set_zero_hk();
            op_local_->contributeHk(ik);
            std::fill(Hpsi_buf.begin(), Hpsi_buf.end(), TK(0));
            apply_Hk(hsk_one_->get_hk(), psi_k, Hpsi_buf.data());
            compute_diagonal(psi_k, Hpsi_buf.data(), cached_h_one_diag_[ik].data(), ik);
        }
        hone_cache_valid_ = true;
    }

    // For each k-point, compute Hartree diag and accumulate energies
    const int nb_local = wfc_eval.get_nbands();
    const int nbs_local = wfc_eval.get_nbasis();
    std::vector<double> vh_diag(nbands_, 0.0);
    std::vector<double> vx_diag(nbands_, 0.0);
    std::vector<TK> Hpsi_buf(nb_local * nbs_local);

    E_one_ = 0.0;
    E_hartree_ = 0.0;
    E_xc_ = 0.0;

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* psi_k = &wfc_eval(ik, 0, 0);

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
        if (exx_enabled_)
        {
            hsk_exx_->set_zero_hk();
            if (GlobalC::exx_info.info_ri.real_number)
            {
                RI_2D_Comm::add_Hexx(*ucell_, *kv_, ik,
                    1.0, exx_lri_d_->Hexxs, *ParaV_, hsk_exx_->get_hk());
            }
            else
            {
                RI_2D_Comm::add_Hexx(*ucell_, *kv_, ik,
                    1.0, exx_lri_c_->Hexxs, *ParaV_, hsk_exx_->get_hk());
            }
            std::fill(Hpsi_buf.begin(), Hpsi_buf.end(), TK(0));
            apply_Hk(hsk_exx_->get_hk(), psi_k, Hpsi_buf.data());
            compute_diagonal(psi_k, Hpsi_buf.data(), vx_diag.data(), ik);
        }
#endif

        double wk = kv_->wk[ik];
        for (int ib = 0; ib < nbands_; ++ib)
        {
            const double n = occ_flat[ik * nbands_ + ib];
            const double gn = xc_func_.g(n);
            if (!rdmft_skip_occ_weight(n))
            {
                E_one_ += wk * n * cached_h_one_diag_[ik][ib];
                E_hartree_ += wk * n * vh_diag[ib] * 0.5;
            }
            if (!rdmft_skip_occ_weight(gn))
            {
                E_xc_ += wk * gn * vx_diag[ib] * 0.5;
            }
        }
    }

    // As in compute(), the *_diag arrays are replicated across ranks so no MPI
    // reduction of E_* is needed here.

    E_ewald_ = pelec_->f_en.ewald_energy;
    E_total_ = E_one_ + E_hartree_ + E_xc_ + E_ewald_;

    ModuleBase::timer::end("RDMFT_EG", "compute_energy");
    return E_total_;
}

// ============================================================================
// Cholesky-based S^{1/2} variable transformation
// ============================================================================

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::precompute_cholesky_S()
{
    ModuleBase::timer::start("RDMFT_EG", "precompute_cholesky_S");

    // Check whether we have an overlap operator at all
    if (op_overlap_ == nullptr || hsk_overlap_ == nullptr)
    {
        // PW basis or no overlap: X = C, nothing to do
        use_X_variable_ = false;
        cholesky_precomputed_ = false;
        ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
        return;
    }

    Uk_.resize(nk_);
    Uk_inv_.resize(nk_);

#ifdef __MPI
    const int nbasis = ParaV_->desc[2];
    const std::int64_t nloc = static_cast<std::int64_t>(ParaV_->nloc);
    const std::int64_t nloc_alloc = std::max<std::int64_t>(nloc, 1);

    for (int ik = 0; ik < nk_; ++ik)
    {
        // Get S_k (distributed, column-major)
        const TK* SK = get_SK(ik);
        if (SK == nullptr)
        {
            GlobalV::ofs_running << "WARNING: overlap matrix S_k is unavailable at ik=" << ik
                << "; disabling X-variable mode." << std::endl;
            use_X_variable_ = false;
            cholesky_precomputed_ = false;
            ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
            return;
        }

        // Copy S_k into Uk_[ik] for in-place Cholesky
        Uk_[ik].assign(static_cast<size_t>(nloc_alloc), TK(0));
        if (nloc > 0)
        {
            std::copy(SK, SK + nloc, Uk_[ik].begin());
        }

        // Upper Cholesky: S = U^H U (uplo='U')
        {
            int info = 0;
            int one_int = 1;
            char uplo = 'U';
            if (std::is_same<TK, double>::value)
            {
                pdpotrf_(&uplo, const_cast<int*>(&nbasis),
                    reinterpret_cast<double*>(Uk_[ik].data()),
                    &one_int, &one_int, const_cast<int*>(ParaV_->desc), &info);
            }
            else
            {
                pzpotrf_(&uplo, const_cast<int*>(&nbasis),
                    reinterpret_cast<std::complex<double>*>(Uk_[ik].data()),
                    &one_int, &one_int, const_cast<int*>(ParaV_->desc), &info);
            }
            if (info != 0)
            {
                GlobalV::ofs_running << "WARNING: Cholesky factorisation of S_k failed"
                    << " at ik=" << ik << " (info=" << info
                    << "); disabling X-variable mode." << std::endl;
                use_X_variable_ = false;
                cholesky_precomputed_ = false;
                ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
                return;
            }
        }

        // Compute U^{-1} by inverting the upper triangular factor
        Uk_inv_[ik] = Uk_[ik];
        {
            int info = 0;
            int one_int = 1;
            char uplo = 'U';
            char diag = 'N';
            if (std::is_same<TK, double>::value)
            {
                pdtrtri_(&uplo, &diag, const_cast<int*>(&nbasis),
                    reinterpret_cast<double*>(Uk_inv_[ik].data()),
                    &one_int, &one_int, const_cast<int*>(ParaV_->desc), &info);
            }
            else
            {
                pztrtri_(&uplo, &diag, const_cast<int*>(&nbasis),
                    reinterpret_cast<std::complex<double>*>(Uk_inv_[ik].data()),
                    &one_int, &one_int, const_cast<int*>(ParaV_->desc), &info);
            }
            if (info != 0)
            {
                GlobalV::ofs_running << "WARNING: Triangular inverse of U_k failed"
                    << " at ik=" << ik << " (info=" << info
                    << "); disabling X-variable mode." << std::endl;
                use_X_variable_ = false;
                cholesky_precomputed_ = false;
                ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
                return;
            }
        }
    }
#else
    // Non-MPI: nbs_local == nbasis
    const int nbasis = nbasis_local_;

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* SK = get_SK(ik);
        if (SK == nullptr)
        {
            Uk_[ik].clear();
            Uk_inv_[ik].clear();
            continue;
        }

        // Copy S_k into Uk_ for in-place Cholesky
        Uk_[ik].assign(SK, SK + nbasis * nbasis);

        // Upper Cholesky: S = U^H U
        int info = detail::potrf_upper(Uk_[ik].data(), nbasis);
        if (info != 0)
        {
            GlobalV::ofs_running << "WARNING: Cholesky factorisation of S_k failed"
                << " at ik=" << ik << " (info=" << info
                << "); disabling X-variable mode." << std::endl;
            use_X_variable_ = false;
            cholesky_precomputed_ = false;
            ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
            return;
        }

        // Compute U^{-1}
        Uk_inv_[ik] = Uk_[ik];
        info = detail::trtri_upper(Uk_inv_[ik].data(), nbasis);
        if (info != 0)
        {
            GlobalV::ofs_running << "WARNING: Triangular inverse of U_k failed"
                << " at ik=" << ik << " (info=" << info
                << "); disabling X-variable mode." << std::endl;
            use_X_variable_ = false;
            cholesky_precomputed_ = false;
            ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
            return;
        }
    }
#endif

    cholesky_precomputed_ = true;
    use_X_variable_ = true;
    GlobalV::ofs_running << "  Cholesky S = U^H U precomputed for all k-points; "
        << "using X_k = U_k C_k as Stiefel variable." << std::endl;

    ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::wfc_C_to_X(psi::Psi<TK>& wfc)
{
    // X_k = U_k * C_k  for each k-point
    if (!cholesky_precomputed_) return;

    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();

#ifdef __MPI
    const int nbasis_const = ParaV_->desc[2];
    const int nbands_const = ParaV_->desc_wfc[3];
    int nbasis = nbasis_const;
    int nbands = nbands_const;

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* C = &wfc(ik, 0, 0);

        // X = U * C via pdtrmm_: B <- alpha * A * B (side='L', uplo='U', trans='N')
        // pdtrmm_ is an in-place operation that overwrites B with A*B,
        // so we operate directly on C (which becomes X after the call).
        int one_int = 1;
        char side = 'L', uplo = 'U', trans = 'N', diag = 'N';
        if (std::is_same<TK, double>::value)
        {
            double alpha_r = 1.0;
            pdtrmm_(&side, &uplo, &trans, &diag, &nbasis, &nbands,
                &alpha_r,
                reinterpret_cast<double*>(const_cast<TK*>(Uk_[ik].data())),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc),
                reinterpret_cast<double*>(C),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc_wfc));
        }
        else
        {
            std::complex<double> alpha_c = {1.0, 0.0};
            pztrmm_(&side, &uplo, &trans, &diag, &nbasis, &nbands,
                &alpha_c,
                reinterpret_cast<std::complex<double>*>(const_cast<TK*>(Uk_[ik].data())),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc),
                reinterpret_cast<std::complex<double>*>(C),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc_wfc));
        }
    }
#else
    const int nbasis = nbs_local;
    const int nbands = nb_local;

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* C = &wfc(ik, 0, 0);

        // X = U * C (in-place via trmm)
        detail::trmm_left_upper_notr(nbasis, nbands, Uk_[ik].data(), nbasis, C, nbasis);
    }
#endif
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::wfc_X_to_C(psi::Psi<TK>& wfc)
{
    // C_k = U_k^{-1} * X_k  for each k-point
    if (!cholesky_precomputed_) return;

    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();

#ifdef __MPI
    int nbasis = ParaV_->desc[2];
    int nbands = ParaV_->desc_wfc[3];

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* X = &wfc(ik, 0, 0);

        // C = U^{-1} X: solve U * C = X for C, i.e. trsm side='L', uplo='U', trans='N'
        int one_int = 1;
        char side = 'L', uplo = 'U', trans = 'N', diag = 'N';
        if (std::is_same<TK, double>::value)
        {
            double alpha_r = 1.0;
            pdtrsm_(&side, &uplo, &trans, &diag, &nbasis, &nbands,
                &alpha_r,
                reinterpret_cast<const double*>(Uk_[ik].data()),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc),
                reinterpret_cast<double*>(X),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc_wfc));
        }
        else
        {
            std::complex<double> alpha_c = {1.0, 0.0};
            pztrsm_(&side, &uplo, &trans, &diag, &nbasis, &nbands,
                &alpha_c,
                reinterpret_cast<const std::complex<double>*>(Uk_[ik].data()),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc),
                reinterpret_cast<std::complex<double>*>(X),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc_wfc));
        }
    }
#else
    const int nbasis = nbs_local;
    const int nbands = nb_local;

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* X = &wfc(ik, 0, 0);

        // C = U^{-1} X via trsm
        detail::trsm_left_upper_notr(nbasis, nbands, Uk_[ik].data(), nbasis, X, nbasis);
    }
#endif
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::grad_C_to_X(psi::Psi<TK>& grad_wfc)
{
    // G_X = U_k^{-H} * G_C  for each k-point
    // (chain rule: C = U^{-1} X => dE/dX* = U^{-H} dE/dC*)
    if (!cholesky_precomputed_) return;

    const int nb_local = grad_wfc.get_nbands();
    const int nbs_local = grad_wfc.get_nbasis();

#ifdef __MPI
    int nbasis = ParaV_->desc[2];
    int nbands = ParaV_->desc_wfc[3];

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* G = &grad_wfc(ik, 0, 0);

        // G_X = U^{-H} G_C: solve U^H * G_X = G_C, i.e.
        // trsm side='L', uplo='U', trans='C'/'T'
        int one_int = 1;
        char side = 'L', uplo = 'U';
        char trans = detail::trans_char(TK());
        char diag = 'N';
        if (std::is_same<TK, double>::value)
        {
            double alpha_r = 1.0;
            pdtrsm_(&side, &uplo, &trans, &diag, &nbasis, &nbands,
                &alpha_r,
                reinterpret_cast<const double*>(Uk_[ik].data()),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc),
                reinterpret_cast<double*>(G),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc_wfc));
        }
        else
        {
            std::complex<double> alpha_c = {1.0, 0.0};
            pztrsm_(&side, &uplo, &trans, &diag, &nbasis, &nbands,
                &alpha_c,
                reinterpret_cast<const std::complex<double>*>(Uk_[ik].data()),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc),
                reinterpret_cast<std::complex<double>*>(G),
                &one_int, &one_int, const_cast<int*>(ParaV_->desc_wfc));
        }
    }
#else
    const int nbasis = nbs_local;
    const int nbands = nb_local;

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* G = &grad_wfc(ik, 0, 0);

        // G_X = U^{-H} G_C via trsm
        detail::trsm_left_upper_conjt(nbasis, nbands, Uk_[ik].data(), nbasis, G, nbasis);
    }
#endif
}

// Explicit template instantiations
template class Veff_rdmft_local<double, double>;
template class Veff_rdmft_local<std::complex<double>, double>;
template class Veff_rdmft_local<std::complex<double>, std::complex<double>>;

template class EnergyGradient<double, double>;
template class EnergyGradient<std::complex<double>, double>;
template class EnergyGradient<std::complex<double>, std::complex<double>>;

} // namespace rdmft
