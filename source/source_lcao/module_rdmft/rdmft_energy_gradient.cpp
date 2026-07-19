//==========================================================
// Author: Kai Luo
// Email: kluo@njust.edu.cn
// DATE : April 2026
//==========================================================
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
#include <sstream>
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
            if (paraV->get_nrow_atom(iat1) <= 0 || paraV->get_ncol_atom(iat2) <= 0)
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
    // RDMFT keeps a single spin slot per HContainer (HR_hartree_, HR_one_,
    // HR_exx_): the per-spin matrix elements are always assembled later by
    // `compute()` itself (folding hR(R) -> hk(k) for every (k, sigma) row of
    // `pelec->wg`). For collinear nspin=2, V_H is identical for both spin
    // components (PotHartree returns v(0)=v(1)=V_H), and V_xc is also added
    // through the same single-slot hR, so we must only accumulate the
    // potential ONCE — looping over `is` and calling cal_gint_vl(v(is), hR)
    // for each spin would double-count V_H (and V_xc) and yield band matrix
    // elements that are 2x larger than KS-LCAO Veff, breaking the RDMFT
    // energy and gradient identities for nspin=2.
    double* vr_eff = nullptr;
    if (potential == "hartree")
    {
        ModuleBase::matrix v(nspin, charge->nrxx);
        elecstate::PotHartree potH(rho_basis);
        potH.cal_v_eff(charge, ucell, v);
        vr_eff = &v(0, 0);
        ModuleGint::cal_gint_vl(vr_eff, hR);
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
        // Same single-slot single-add reasoning as the "hartree" branch above.
        vr_eff = &v(0, 0);
        ModuleGint::cal_gint_vl(vr_eff, hR);
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
EnergyGradient<TK, TR>::~EnergyGradient() = default;

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
    // Use the actual k-spin row count used by occupations / psi in this SCF.
    // kv::get_nks() may already include spin in some paths; multiplying by
    // nspin again can over-count and later drive ik loops out of range.
    nk_ = pelec_->wg.nr;
    nbasis_local_ = ParaV_->nrow;

#ifdef __MPI
    para_Eij_.set(nbands_, nbands_, ParaV_->nb, ParaV_->blacs_ctxt);
#endif

    // Allocate Hamiltonian containers
    HR_one_ = std::make_unique<hamilt::HContainer<TR>>(*ucell_, ParaV_);
    HR_hartree_ = std::make_unique<hamilt::HContainer<TR>>(*ucell_, ParaV_);
    HR_exx_ = std::make_unique<hamilt::HContainer<TR>>(*ucell_, ParaV_);
    SR_ = std::make_unique<hamilt::HContainer<TR>>(*ucell_, ParaV_);

    hsk_one_ = std::make_unique<hamilt::HS_Matrix_K<TK>>(ParaV_, true);
    hsk_hartree_ = std::make_unique<hamilt::HS_Matrix_K<TK>>(ParaV_, true);
    hsk_exx_ = std::make_unique<hamilt::HS_Matrix_K<TK>>(ParaV_, true);
    // The overlap operator writes S(k) into hsk->sk via hsk->get_sk(),
    // so we must allocate sk (second arg false -> no_s==false)
    hsk_overlap_ = std::make_unique<hamilt::HS_Matrix_K<TK>>(ParaV_, false);

    if (PARAM.inp.gamma_only)
    {
        // Match KS-LCAO's gamma-only convention exactly (HamiltLCAO::HamiltLCAO,
        // hamilt_lcao.cpp:136): pre-fix_gamma the H-side containers but
        // **leave SR_ alone**.
        //
        // Why: the operators that fill HR (EKinetic / Nonlocal /
        // Veff_rdmft_local / Exx_LRI) all iterate `adjs_all` over every
        // neighbour cell and += accumulate into `find_matrix(iat1,iat2,R)`,
        // which on a gamma-only HContainer returns the single (0,0,0)
        // BaseMatrix for every R query. So pre-fix_gamma is harmless on the
        // H-side and matches KS exactly.
        //
        // BUT `Overlap::calculate_SR` (`overlap.cpp:143-163`) loops the
        // AtomPair's `tmp.get_R_size()` instead of `adjs_all`, computes per-R
        // S(R) values into per-R BaseMatrix entries, and only at the end of
        // `calculate_SR` calls `SR->fix_gamma()` itself when TK==double.  If
        // we pre-fix_gamma SR_ here, `populate_atom_pairs` collapses every
        // inserted AtomPair to a single (0,0,0) entry **before** values are
        // computed, so `calculate_SR` then only iterates that one entry and
        // misses the periodic-image S(R≠0) sum.  The result is an
        // S(Γ) that is **smaller** than the one KS-LCAO uses, and the wfc
        // ABACUS hands to RDMFT (normalised against the KS S(Γ)) gives
        // wrong matrix elements `<ψ|h|ψ>`, `<ψ|V_H|ψ>`, `<ψ|H_exx|ψ>` —
        // exactly the H2_HF-on-gamma-only-with-small-cell symptom we
        // diagnosed earlier (`||S||_F^2` was 10.04 in the buggy gamma-only
        // path vs 11.87 with the multi-k single-Γ path that *did* match KS).
        //
        // Fix: only fix_gamma the H-side; leave SR_ to be populated and then
        // auto-fix_gamma'd by `Overlap::calculate_SR` itself, exactly as KS
        // does (it never pre-fix_gamma's its sR; see hamilt_lcao.cpp).
        HR_one_->fix_gamma();
        HR_hartree_->fix_gamma();
        HR_exx_->fix_gamma();
    }

#ifdef __EXX
    if (PARAM.inp.rdmft)
    {
        // RDMFT manages its own EXX context and must not inherit KS-side
        // `dft_functional` EXX switches.
        GlobalC::exx_info.info_global.cal_exx = true;
        GlobalC::exx_info.info_global.ccp_type = Conv_Coulomb_Pot_K::Ccp_Type::Hf;
        GlobalC::exx_info.info_global.coulomb_param.clear();
        std::string sing_corr = PARAM.inp.exx_singularity_correction;
        if (sing_corr == "default") { sing_corr = "spencer"; }
        GlobalC::exx_info.info_global.coulomb_param[Conv_Coulomb_Pot_K::Coulomb_Type::Fock] = {{
            {"alpha", "1"},
            {"singularity_correction", sing_corr}
        }};
        GlobalC::exx_info.info_global.separate_loop = PARAM.inp.exx_separate_loop;
        GlobalC::exx_info.info_global.hybrid_step = PARAM.inp.exx_hybrid_step;
        GlobalC::exx_info.info_global.mixing_beta_for_loop1 = PARAM.inp.exx_mixing_beta;

        // Propagate input RI parameters (input_conv skips this block when
        // cal_exx is off, but RDMFT needs a fully-initialised info_ri to run
        // its own EXX evaluator).
        GlobalC::exx_info.info_ri.real_number = std::stoi(PARAM.inp.exx_real_number);
        GlobalC::exx_info.info_ri.pca_threshold = PARAM.inp.exx_pca_threshold;
        GlobalC::exx_info.info_ri.C_threshold   = PARAM.inp.exx_c_threshold;
        GlobalC::exx_info.info_ri.V_threshold   = PARAM.inp.exx_v_threshold;
        GlobalC::exx_info.info_ri.dm_threshold  = PARAM.inp.exx_dm_threshold;
        // Keep RDMFT EXX mesh quality independent of KS `dft_functional`-
        // dependent defaults (e.g. pbe default 1 vs hf default 5).
        const double exx_ccp_rmesh_times_input = std::stod(PARAM.inp.exx_ccp_rmesh_times);
        GlobalC::exx_info.info_ri.ccp_rmesh_times = std::max(5.0, exx_ccp_rmesh_times_input);
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
        // RDMFT EXX must always use full Fock (alpha = 1), independent of
        // KS-side dft_functional setup.
        GlobalC::exx_info.info_global.hybrid_alpha = 1.0;
        XC_Functional::set_hybrid_alpha(1.0);

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
            exx_lri_d_ = std::make_unique<Exx_LRI<double>>(GlobalC::exx_info.info_ri);
            exx_lri_d_->init(MPI_COMM_WORLD, const_cast<UnitCell&>(*ucell_), *kv_, *orb_);
        }
        else
        {
            exx_lri_c_ = std::make_unique<Exx_LRI<std::complex<double>>>(GlobalC::exx_info.info_ri);
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
    op_ekinetic_ = std::make_unique<hamilt::EKinetic<hamilt::OperatorLCAO<TK, TR>>>(
        hsk_one_.get(), kv_->kvec_d, HR_one_.get(), ucell_,
        orb_->cutoffs(), gd_, two_center_bundle_->kinetic_orb.get());

    op_nonlocal_ = std::make_unique<hamilt::Nonlocal<hamilt::OperatorLCAO<TK, TR>>>(
        hsk_one_.get(), kv_->kvec_d, HR_one_.get(), ucell_,
        orb_->cutoffs(), gd_, two_center_bundle_->overlap_orb_beta.get());

    op_local_ = std::make_unique<Veff_rdmft_local<TK, TR>>(
        hsk_one_.get(), kv_->kvec_d, pelec_->pot, HR_one_.get(), ucell_,
        orb_->cutoffs(), gd_, nspin_, charge_, rho_basis_, vloc_, sf_, "local");

    // Overlap operator (builds SR internally)
    op_overlap_ = std::make_unique<hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>>(
        hsk_overlap_.get(), kv_->kvec_d, SR_.get(), SR_.get(), ucell_,
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
        // `nk_` (== pelec_->wg.nr) is the FLATTENED (k, spin) row count, so it
        // already includes the spin doubling from K_Vectors::set_kup_and_kdw.
        // The DensityMatrix multi-k ctor expects the per-spin k-point count
        // (DMK is internally sized as nk_per_spin * nspin), so we have to pass
        // nk_/nspin_ here. Otherwise nspin=2 mis-orders DMK → DMR (mixing
        // spin-up bands into the spin-up DMR while leaving spin-down DMR
        // empty).
        const int nk_per_spin = nk_ / nspin_;
        DM.reset(new elecstate::DensityMatrix<TK, double>(ParaV_, nspin_, kv_->kvec_d, nk_per_spin));
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
    std::vector<std::vector<TK>>& DM_XC,
    double alpha_override,
    const std::vector<double>* occ_weight_override)
{
    DM_XC.resize(nk_, std::vector<TK>(ParaV_->nloc, TK(0)));

    // Build wk_g matrix: wk[ik] * g(n(ik, ib)). Bands with g(n)=0 do not enter the XC density matrix.
    ModuleBase::matrix wk_g(nk_, nbands_);
    const bool use_override = (alpha_override > 0.0);
    const bool use_occ_weight_override = (occ_weight_override != nullptr);
    for (int ik = 0; ik < nk_; ++ik)
    {
        for (int ib = 0; ib < nbands_; ++ib)
        {
            const double n = occ_flat[ik * nbands_ + ib];
            double gn = 0.0;
            if (use_occ_weight_override)
            {
                gn = (*occ_weight_override)[ik * nbands_ + ib];
            }
            else if (use_override)
            {
                gn = xc_func_.pow_reg(n, alpha_override);
            }
            else
            {
                gn = xc_func_.g(n);
            }
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
        // See comment in build_charge(): nk_ is flattened over (k, spin), but
        // the DensityMatrix multi-k ctor wants the per-spin k-count.
        const int nk_per_spin = nk_ / nspin_;
        DM_xc.reset(new elecstate::DensityMatrix<TK, double>(ParaV_, nspin_, kv_->kvec_d, nk_per_spin));
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

    // Compute Eij = psi^H * Hpsi in the 2D BLACS grid.
    // Keep a non-null storage even when this rank owns zero local blocks.
    const int eij_nloc = para_Eij_.get_row_size() * para_Eij_.get_col_size();
    const std::int64_t eij_alloc = std::max<std::int64_t>(static_cast<std::int64_t>(eij_nloc), 1);
    std::vector<TK> Eij(static_cast<size_t>(eij_alloc), TK(0));

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
    if (!op_overlap_ || !hsk_overlap_) return nullptr;
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

/// f(n) = n ln n + (1-n) ln(1-n); n clamped like OccupationParam::to_param
inline double binary_entropy_f(double n)
{
    constexpr double eps = 1e-12;
    n = std::max(eps, std::min(1.0 - eps, n));
    return n * std::log(n) + (1.0 - n) * std::log(1.0 - n);
}

inline double binary_entropy_dfdn(double n)
{
    constexpr double eps = 1e-12;
    n = std::max(eps, std::min(1.0 - eps, n));
    return std::log(n / (1.0 - n));
}

template <typename TK>
inline const TK* psi_k_ptr_or_dummy(const psi::Psi<TK>& psi, int ik, std::vector<TK>& dummy)
{
    if (psi.get_nbands() > 0 && psi.get_nbasis() > 0)
    {
        return &psi(ik, 0, 0);
    }
    if (dummy.empty())
    {
        dummy.assign(1, TK(0));
    }
    else
    {
        dummy[0] = TK(0);
    }
    return dummy.data();
}

template <typename TK>
inline TK* psi_k_ptr_or_dummy(psi::Psi<TK>& psi, int ik, std::vector<TK>& dummy)
{
    if (psi.get_nbands() > 0 && psi.get_nbasis() > 0)
    {
        return &psi(ik, 0, 0);
    }
    if (dummy.empty())
    {
        dummy.assign(1, TK(0));
    }
    else
    {
        dummy[0] = TK(0);
    }
    return dummy.data();
}

/// ∑_k ‖P_k‖_F^2 = ∑ |P|² in the LCAO coefficient layout (ambient Euclidean on coeffs).
template <typename TK>
double ambient_frobenius_norm2_psi(const psi::Psi<TK>& P)
{
    double s = 0.0;
    const int nk = P.get_nk();
    const int nb = P.get_nbands();
    const int nbs = P.get_nbasis();
#ifdef __MPI
    std::vector<TK> dummy(1);
    for (int ik = 0; ik < nk; ++ik)
    {
        const TK* pk = psi_k_ptr_or_dummy(P, ik, dummy);
        for (int i = 0; i < nb * nbs; ++i)
        {
            s += real_of_conj_prod(pk[i], pk[i]);
        }
    }
    Parallel_Reduce::reduce_all(s);
#else
    std::vector<TK> dummy(1);
    for (int ik = 0; ik < nk; ++ik)
    {
        const TK* pk = psi_k_ptr_or_dummy(P, ik, dummy);
        for (int i = 0; i < nb * nbs; ++i)
        {
            s += real_of_conj_prod(pk[i], pk[i]);
        }
    }
#endif
    return s;
}

template <typename TK>
double orb_grad_decomp_residual_frob2(const psi::Psi<TK>& full,
    const psi::Psi<TK>& g1,
    const psi::Psi<TK>& gh,
    const psi::Psi<TK>& gx)
{
    double s = 0.0;
    const int nk = full.get_nk();
    const int nb = full.get_nbands();
    const int nbs = full.get_nbasis();
#ifdef __MPI
    std::vector<TK> d0(1), d1(1), d2(1), d3(1);
    for (int ik = 0; ik < nk; ++ik)
    {
        const TK* pf = psi_k_ptr_or_dummy(full, ik, d0);
        const TK* p1 = psi_k_ptr_or_dummy(g1, ik, d1);
        const TK* ph = psi_k_ptr_or_dummy(gh, ik, d2);
        const TK* px = psi_k_ptr_or_dummy(gx, ik, d3);
        for (int i = 0; i < nb * nbs; ++i)
        {
            const TK d = pf[i] - (p1[i] + ph[i] + px[i]);
            s += real_of_conj_prod(d, d);
        }
    }
    Parallel_Reduce::reduce_all(s);
#else
    std::vector<TK> d0(1), d1(1), d2(1), d3(1);
    for (int ik = 0; ik < nk; ++ik)
    {
        const TK* pf = psi_k_ptr_or_dummy(full, ik, d0);
        const TK* p1 = psi_k_ptr_or_dummy(g1, ik, d1);
        const TK* ph = psi_k_ptr_or_dummy(gh, ik, d2);
        const TK* px = psi_k_ptr_or_dummy(gx, ik, d3);
        for (int i = 0; i < nb * nbs; ++i)
        {
            const TK d = pf[i] - (p1[i] + ph[i] + px[i]);
            s += real_of_conj_prod(d, d);
        }
    }
#endif
    return s;
}

} // namespace

namespace {

/// In-place at one k-point: column-major `G` is the ambient Riesz gradient Ḡ
/// (pairing Re Tr(Ḡ^H ·)).  `X` is the Stiefel point (X^H X = I).  On exit, `G`
/// is the canonical-metric Riemannian gradient ξ with
///   ⟨ξ, V⟩_can = Re Tr(Ḡ^H V)  for every tangent V at X,
/// using the same single-k term as `stiefel_canonical_inner_product`:
///   ⟨U,V⟩_can = Re Tr(U^H V) − ½ Re Tr(U^H X X^H V).
///
/// Explicit construction:
///   T = Ḡ − X sym(X^H Ḡ),
///   W = T + X (X^H T) = (I + X X^H) T,
///   Λ = −¼ (X^H W + W^H X),  ξ = W + 2 X Λ.
template <typename TK>
void canonical_stiefel_riemannian_gradient_one_k_dense(const TK* X,
    TK* G,
    const int nbasis,
    const int nbands,
    const char tc,
    TK* work_W,
    TK* SC,
    TK* A,
    TK* B)
{
    const TK one = TK(1.0);
    const TK zero = TK(0.0);
    const TK neg_one = TK(-1.0);
    const TK two = TK(2.0);
    const TK quarter = TK(0.25);
    const int npsi = nbasis * nbands;

    for (int i = 0; i < npsi; ++i)
    {
        SC[i] = X[i];
    }

    detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, SC, nbasis, G, nbasis, zero, A, nbands);
    detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, G, nbasis, SC, nbasis, zero, B, nbands);
    for (int i = 0; i < nbands * nbands; ++i)
    {
        A[i] = TK(0.5) * (A[i] + B[i]);
    }

    detail::gemm_wrapper('N', 'N', nbasis, nbands, nbands,
        neg_one, X, nbasis, A, nbands, one, G, nbasis);

    detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, SC, nbasis, G, nbasis, zero, B, nbands);

    for (int i = 0; i < npsi; ++i)
    {
        work_W[i] = G[i];
    }
    detail::gemm_wrapper('N', 'N', nbasis, nbands, nbands,
        one, X, nbasis, B, nbands, one, work_W, nbasis);

    detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, SC, nbasis, work_W, nbasis, zero, A, nbands);
    detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, work_W, nbasis, SC, nbasis, zero, B, nbands);
    for (int i = 0; i < nbands * nbands; ++i)
    {
        B[i] = -quarter * (A[i] + B[i]);
    }

    for (int i = 0; i < npsi; ++i)
    {
        G[i] = work_W[i];
    }
    detail::gemm_wrapper('N', 'N', nbasis, nbands, nbands,
        two, X, nbasis, B, nbands, one, G, nbasis);
}

#ifdef __MPI
template <typename TK>
void canonical_stiefel_riemannian_gradient_one_k_scalapack(const TK* X,
    TK* G,
    const int nbasis,
    const int nbands,
    const int npsi,
    const int eij_nloc,
    const int* desc_wfc,
    const int* desc_eij,
    const char tc,
    TK* work_W,
    TK* SC,
    TK* A,
    TK* B)
{
    const TK one = TK(1.0);
    const TK zero = TK(0.0);
    const TK neg_one = TK(-1.0);
    const TK two = TK(2.0);
    const TK quarter = TK(0.25);

    for (int i = 0; i < npsi; ++i)
    {
        SC[i] = X[i];
    }

    detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, SC, 1, 1, desc_wfc,
        G, 1, 1, desc_wfc,
        zero, A, 1, 1, desc_eij);
    detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, G, 1, 1, desc_wfc,
        SC, 1, 1, desc_wfc,
        zero, B, 1, 1, desc_eij);
    for (int i = 0; i < eij_nloc; ++i)
    {
        A[i] = TK(0.5) * (A[i] + B[i]);
    }

    detail::pgemm_wrapper('N', 'N', nbasis, nbands, nbands,
        neg_one, X, 1, 1, desc_wfc,
        A, 1, 1, desc_eij,
        one, G, 1, 1, desc_wfc);

    detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, SC, 1, 1, desc_wfc,
        G, 1, 1, desc_wfc,
        zero, B, 1, 1, desc_eij);

    for (int i = 0; i < npsi; ++i)
    {
        work_W[i] = G[i];
    }
    detail::pgemm_wrapper('N', 'N', nbasis, nbands, nbands,
        one, X, 1, 1, desc_wfc,
        B, 1, 1, desc_eij,
        one, work_W, 1, 1, desc_wfc);

    detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, SC, 1, 1, desc_wfc,
        work_W, 1, 1, desc_wfc,
        zero, A, 1, 1, desc_eij);
    detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
        one, work_W, 1, 1, desc_wfc,
        SC, 1, 1, desc_wfc,
        zero, B, 1, 1, desc_eij);
    for (int i = 0; i < eij_nloc; ++i)
    {
        B[i] = -quarter * (A[i] + B[i]);
    }

    for (int i = 0; i < npsi; ++i)
    {
        G[i] = work_W[i];
    }
    detail::pgemm_wrapper('N', 'N', nbasis, nbands, nbands,
        two, X, 1, 1, desc_wfc,
        B, 1, 1, desc_eij,
        one, G, 1, 1, desc_wfc);
}
#endif

} // namespace

template <typename TK, typename TR>
double EnergyGradient<TK, TR>::stiefel_canonical_inner_product(
    const psi::Psi<TK>& X_stiefel,
    const psi::Psi<TK>& U,
    const psi::Psi<TK>& V)
{
    double sum_frob = 0.0;
    double sum_cross = 0.0;
    const TK one = TK(1.0);
    const TK zero = TK(0.0);
    const char tc = detail::trans_char(TK());

#ifdef __MPI
    const int nbasis = ParaV_->desc[2];
    const int nbands = ParaV_->desc_wfc[3];
    const int nb_local = U.get_nbands();
    const int nbs_local = U.get_nbasis();
    const int nrow = para_Eij_.get_row_size();
    const int ncol = para_Eij_.get_col_size();
    const int eij_nloc = nrow * ncol;
    const std::int64_t sc_alloc
        = std::max<std::int64_t>(static_cast<std::int64_t>(nb_local * nbs_local), 1);
    const std::int64_t eij_alloc = std::max<std::int64_t>(static_cast<std::int64_t>(eij_nloc), 1);
    std::vector<TK> x_dummy(1, TK(0));
    std::vector<TK> u_dummy(1, TK(0));
    std::vector<TK> v_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* Xk = psi_k_ptr_or_dummy(X_stiefel, ik, x_dummy);
        const TK* Uk = psi_k_ptr_or_dummy(U, ik, u_dummy);
        const TK* Vk = psi_k_ptr_or_dummy(V, ik, v_dummy);
        const int npsi = nb_local * nbs_local;
        for (int i = 0; i < npsi; ++i)
        {
            sum_frob += real_of_conj_prod(Uk[i], Vk[i]);
        }

        std::vector<TK> SC(static_cast<size_t>(sc_alloc), TK(0));
        for (int i = 0; i < npsi; ++i)
        {
            SC[i] = Xk[i];
        }
        std::vector<TK> XU(static_cast<size_t>(eij_alloc), TK(0));
        std::vector<TK> XV(static_cast<size_t>(eij_alloc), TK(0));
        detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, SC.data(), 1, 1, ParaV_->desc_wfc,
            Uk, 1, 1, ParaV_->desc_wfc,
            zero, XU.data(), 1, 1, para_Eij_.desc);
        detail::pgemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, SC.data(), 1, 1, ParaV_->desc_wfc,
            Vk, 1, 1, ParaV_->desc_wfc,
            zero, XV.data(), 1, 1, para_Eij_.desc);

        double local_cross_k = 0.0;
        for (int j = 0; j < ncol; ++j)
        {
            for (int i = 0; i < nrow; ++i)
            {
                local_cross_k += real_of_conj_prod(XU[i + j * nrow], XV[i + j * nrow]);
            }
        }
        Parallel_Reduce::reduce_all(local_cross_k);
        sum_cross += local_cross_k;
    }
    Parallel_Reduce::reduce_all(sum_frob);
#else
    const int nbasis = U.get_nbasis();
    const int nbands = U.get_nbands();
    std::vector<TK> x_dummy(1, TK(0));
    std::vector<TK> u_dummy(1, TK(0));
    std::vector<TK> v_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* Xk = psi_k_ptr_or_dummy(X_stiefel, ik, x_dummy);
        const TK* Uk = psi_k_ptr_or_dummy(U, ik, u_dummy);
        const TK* Vk = psi_k_ptr_or_dummy(V, ik, v_dummy);
        const int ntot = nbasis * nbands;
        for (int i = 0; i < ntot; ++i)
        {
            sum_frob += real_of_conj_prod(Uk[i], Vk[i]);
        }
        std::vector<TK> SC(static_cast<size_t>(ntot), TK(0));
        for (int i = 0; i < ntot; ++i)
        {
            SC[i] = Xk[i];
        }
        std::vector<TK> XU(static_cast<size_t>(nbands * nbands), TK(0));
        std::vector<TK> XV(static_cast<size_t>(nbands * nbands), TK(0));
        detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, SC.data(), nbasis, Uk, nbasis, zero, XU.data(), nbands);
        detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, SC.data(), nbasis, Vk, nbasis, zero, XV.data(), nbands);
        for (int j = 0; j < nbands; ++j)
        {
            for (int i = 0; i < nbands; ++i)
            {
                sum_cross += real_of_conj_prod(XU[i + j * nbands], XV[i + j * nbands]);
            }
        }
    }
#endif
    return sum_frob - 0.5 * sum_cross;
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::stiefel_gram_residual_frobenius_per_k(
    const psi::Psi<TK>& wfc,
    std::vector<double>& frob_per_ik)
{
    frob_per_ik.assign(nk_, 0.0);
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
    std::vector<TK> c_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* C = psi_k_ptr_or_dummy(wfc, ik, c_dummy);
        for (int i = 0; i < nb_local * nbs_local; ++i)
        {
            SY[i] = C[i];
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
    std::vector<TK> c_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* C = psi_k_ptr_or_dummy(wfc, ik, c_dummy);
        for (int i = 0; i < nbasis * nbands; ++i)
        {
            SY[i] = C[i];
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
    // See `canonical_stiefel_riemannian_gradient_one_k_*`: explicit T, W, Λ, ξ
    // matching `stiefel_canonical_inner_product` for pairing and norms.
    const char tc = detail::trans_char(TK());
#ifdef __MPI
    const int nbasis = ParaV_->desc[2];
    const int nbands = ParaV_->desc_wfc[3];
    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();

    const int eij_nloc = para_Eij_.get_row_size() * para_Eij_.get_col_size();
    const std::int64_t sc_alloc
        = std::max<std::int64_t>(static_cast<std::int64_t>(nb_local * nbs_local), 1);
    const std::int64_t eij_alloc = std::max<std::int64_t>(static_cast<std::int64_t>(eij_nloc), 1);
    std::vector<TK> work_W(static_cast<size_t>(sc_alloc), TK(0));
    std::vector<TK> SC(static_cast<size_t>(sc_alloc), TK(0));
    std::vector<TK> A(static_cast<size_t>(eij_alloc), TK(0));
    std::vector<TK> B(static_cast<size_t>(eij_alloc), TK(0));
    std::vector<TK> psi_dummy(1, TK(0));
    std::vector<TK> grad_dummy(1, TK(0));
    const int npsi = nb_local * nbs_local;

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* psi_k = psi_k_ptr_or_dummy(wfc, ik, psi_dummy);
        TK* g_k = psi_k_ptr_or_dummy(grad_wfc, ik, grad_dummy);
        canonical_stiefel_riemannian_gradient_one_k_scalapack<TK>(psi_k,
            g_k,
            nbasis,
            nbands,
            npsi,
            eij_nloc,
            ParaV_->desc_wfc,
            para_Eij_.desc,
            tc,
            work_W.data(),
            SC.data(),
            A.data(),
            B.data());
    }
#else
    const int nbasis = wfc.get_nbasis();
    const int nbands = wfc.get_nbands();
    const int npsi = nbasis * nbands;
    std::vector<TK> work_W(static_cast<size_t>(npsi), TK(0));
    std::vector<TK> SC(static_cast<size_t>(npsi), TK(0));
    std::vector<TK> A(static_cast<size_t>(nbands * nbands), TK(0));
    std::vector<TK> B(static_cast<size_t>(nbands * nbands), TK(0));
    std::vector<TK> psi_dummy(1, TK(0));
    std::vector<TK> grad_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* psi_k = psi_k_ptr_or_dummy(wfc, ik, psi_dummy);
        TK* g_k = psi_k_ptr_or_dummy(grad_wfc, ik, grad_dummy);
        canonical_stiefel_riemannian_gradient_one_k_dense<TK>(psi_k,
            g_k,
            nbasis,
            nbands,
            tc,
            work_W.data(),
            SC.data(),
            A.data(),
            B.data());
    }
#endif
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::retract_orbitals(
    psi::Psi<TK>& wfc,
    const psi::Psi<TK>& grad_wfc,
    double alpha)
{
    // One retraction step on the standard Stiefel manifold in X-space:
    //   Y = X - alpha * G
    //   X_new = Y * (Y^H Y)^{-1/2}
    // We approximate the inverse-square-root via Cholesky of M = Y^H Y:
    //     M = L L^H   =>   Y_ortho = Y * L^{-H}
    // This is the "Cholesky QR" S-orthonormalisation, much cheaper than a
    // Hermitian eigendecomposition and perfectly adequate as long as alpha
    // is chosen small enough that M stays well-conditioned (which the outer
    // line search guarantees).
    //
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
    std::vector<TK> c_dummy(1, TK(0));
    std::vector<TK> g_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        TK* C = psi_k_ptr_or_dummy(wfc, ik, c_dummy);
        const TK* G = psi_k_ptr_or_dummy(grad_wfc, ik, g_dummy);

        // Y = C - alpha * G (in-place on C)
        for (int i = 0; i < nb_local * nbs_local; ++i)
            C[i] -= TK(alpha) * G[i];

        // In X-space, SY is just Y.
        std::vector<TK> SY(static_cast<size_t>(sy_alloc_retract), TK(0));
        for (int i = 0; i < nb_local * nbs_local; ++i) SY[i] = C[i];

        // M = Y^H Y
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
    std::vector<TK> c_dummy(1, TK(0));
    std::vector<TK> g_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        TK* C = psi_k_ptr_or_dummy(wfc, ik, c_dummy);
        const TK* G = psi_k_ptr_or_dummy(grad_wfc, ik, g_dummy);

        // Y = C - alpha * G (in-place on C)
        for (int i = 0; i < nbasis * nbands; ++i)
            C[i] -= TK(alpha) * G[i];

        // In X-space, SY is just Y.
        std::vector<TK> SY(nbasis * nbands, TK(0));
        for (int i = 0; i < nbasis * nbands; ++i) SY[i] = C[i];

        // M = Y^H Y  (nbands x nbands)
        std::vector<TK> M(nbands * nbands, TK(0));
        detail::gemm_wrapper(tc, 'N', nbands, nbands, nbasis,
            one, C, nbasis, SY.data(), nbasis, zero, M.data(), nbands);

        // Cholesky: M = L L^H
        const int info = detail::potrf_lower(M.data(), nbands);
        if (info != 0)
        {
            // Cholesky failed; skip re-orthonormalisation for this step.
            // The outer line search should reject this trial point.
            GlobalV::ofs_running << "WARNING: Cholesky factorisation of Y^H Y failed"
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

    // RDMFT operates in X-space. Convert X_k -> C_k for energy evaluation.
    psi::Psi<TK> wfc_C(wfc);
    wfc_X_to_C(wfc_C);
    const psi::Psi<TK>& wfc_eval = wfc_C;

    build_charge(occ_flat, wfc_eval);

    // 2. Build Hartree potential
    HR_hartree_->set_zero();
    if (!op_hartree_)
    {
        op_hartree_ = std::make_unique<Veff_rdmft_local<TK, TR>>(
            hsk_hartree_.get(), kv_->kvec_d, pelec_->pot, HR_hartree_.get(), ucell_,
            orb_->cutoffs(), gd_, nspin_, charge_, rho_basis_, vloc_, sf_, "hartree");
    }
    op_hartree_->contributeHR();

    // 3. Build exchange from modified DM (RDMFT's private EXX)
    const int nb_local = wfc_eval.get_nbands();
    const int nbs_local = wfc_eval.get_nbasis();
    const std::int64_t hpsi_alloc
        = std::max<std::int64_t>(static_cast<std::int64_t>(nb_local * nbs_local), 1);
    const bool is_mixed_exx = (xc_func_.type() == XCFunctionalType::GEO
                               || xc_func_.type() == XCFunctionalType::HybOpt
                               || xc_func_.type() == XCFunctionalType::CHF
                               || xc_func_.type() == XCFunctionalType::CGA);
    auto mixed_num_terms = [&]() -> int {
        switch (xc_func_.type())
        {
            case XCFunctionalType::GEO: return XCFunctional::num_geo_terms();
            case XCFunctionalType::HybOpt:
            case XCFunctionalType::CHF:
            case XCFunctionalType::CGA: return 2;
            default: return 0;
        }
    };
    auto mixed_term_coeff = [&](int t) -> double {
        switch (xc_func_.type())
        {
            case XCFunctionalType::GEO: return XCFunctional::geo_coef(t);
            case XCFunctionalType::HybOpt:
                return (t == 0) ? XCFunctional::hybopt_hf_weight() : XCFunctional::hybopt_power_weight();
            case XCFunctionalType::CHF:
                return (t == 0) ? XCFunctional::chf_hf_weight() : XCFunctional::chf_corr_weight();
            case XCFunctionalType::CGA:
                return (t == 0) ? XCFunctional::cga_hf_weight() : XCFunctional::cga_corr_weight();
            default: return 0.0;
        }
    };
    auto mixed_term_weight = [&](int t, double n) -> double {
        switch (xc_func_.type())
        {
            case XCFunctionalType::GEO: return xc_func_.pow_reg(n, XCFunctional::geo_alpha(t));
            case XCFunctionalType::HybOpt:
                return (t == 0) ? n : xc_func_.pow_reg(n, XCFunctional::hybopt_power_exponent());
            case XCFunctionalType::CHF:
                return (t == 0) ? n : xc_func_.chf_corr_term(n);
            case XCFunctionalType::CGA:
                return (t == 0) ? n : xc_func_.cga_corr_term(n);
            default: return 0.0;
        }
    };
    auto mixed_term_weight_deriv = [&](int t, double n) -> double {
        switch (xc_func_.type())
        {
            case XCFunctionalType::GEO: return xc_func_.dpow_reg(n, XCFunctional::geo_alpha(t));
            case XCFunctionalType::HybOpt:
                return (t == 0) ? 1.0 : xc_func_.dpow_reg(n, XCFunctional::hybopt_power_exponent());
            case XCFunctionalType::CHF:
                return (t == 0) ? 1.0 : xc_func_.chf_corr_term_deriv(n);
            case XCFunctionalType::CGA:
                return (t == 0) ? 1.0 : xc_func_.cga_corr_term_deriv(n);
            default: return 0.0;
        }
    };
    std::vector<std::vector<double>> mix_vx_E_acc(nk_, std::vector<double>(nbands_, 0.0));
    std::vector<std::vector<double>> mix_vx_G_acc(nk_, std::vector<double>(nbands_, 0.0));
    std::vector<std::vector<TK>> mix_Hpsi_x_acc(nk_, std::vector<TK>(static_cast<size_t>(hpsi_alloc), TK(0)));
#ifdef __EXX
    if (exx_enabled_ && !is_mixed_exx)
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

        }
    }
    else if (exx_enabled_ && is_mixed_exx)
    {
        const int n_terms = mixed_num_terms();
        std::vector<TK> Hpsi_x_term(static_cast<size_t>(hpsi_alloc), TK(0));
        std::vector<double> vx_diag_term(nbands_, 0.0);
        std::vector<TK> psi_dummy_mix(1, TK(0));
        std::vector<double> term_weights(nk_ * nbands_, 0.0);
        for (int t = 0; t < n_terms; ++t)
        {
            const double coeff = mixed_term_coeff(t);
            for (int ik = 0; ik < nk_; ++ik)
            {
                for (int ib = 0; ib < nbands_; ++ib)
                {
                    const double n = occ_flat[ik * nbands_ + ib];
                    term_weights[ik * nbands_ + ib] = mixed_term_weight(t, n);
                }
            }

            HR_exx_->set_zero();
            std::vector<std::vector<TK>> DM_XC;
            build_DM_xc(occ_flat, wfc_eval, DM_XC, 0.0, &term_weights);

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
            }

            for (int ik = 0; ik < nk_; ++ik)
            {
                const TK* psi_k = psi_k_ptr_or_dummy(wfc_eval, ik, psi_dummy_mix);
                std::fill(Hpsi_x_term.begin(), Hpsi_x_term.end(), TK(0));
                std::fill(vx_diag_term.begin(), vx_diag_term.end(), 0.0);
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
                apply_Hk(hsk_exx_->get_hk(), psi_k, Hpsi_x_term.data());
                compute_diagonal(psi_k, Hpsi_x_term.data(), vx_diag_term.data(), ik);

                for (int ib = 0; ib < nbands_; ++ib)
                {
                    const double n = occ_flat[ik * nbands_ + ib];
                    const double p = mixed_term_weight(t, n);
                    const double dp = mixed_term_weight_deriv(t, n);
                    mix_vx_E_acc[ik][ib] += coeff * p * vx_diag_term[ib];
                    mix_vx_G_acc[ik][ib] += coeff * dp * vx_diag_term[ib];
                }

                for (int ib_local = 0; ib_local < nb_local; ++ib_local)
                {
                    const int ib_global = ParaV_->local2global_col(ib_local);
                    if (ib_global >= nbands_) continue;
                    const double n = occ_flat[ik * nbands_ + ib_global];
                    const double p = mixed_term_weight(t, n);
                    const double scale = coeff * p;
                    TK* acc_ptr = &mix_Hpsi_x_acc[ik][static_cast<size_t>(ib_local) * nbs_local];
                    const TK* term_ptr = &Hpsi_x_term[static_cast<size_t>(ib_local) * nbs_local];
                    for (int mu = 0; mu < nbs_local; ++mu)
                    {
                        acc_ptr[mu] += TK(scale) * term_ptr[mu];
                    }
                }
            }
        }
    }
#endif

    // 4. For each k-point: compute H*psi, diagonal elements, and assemble gradients

    grad_occ.assign(nk_ * nbands_, 0.0);
    grad_wfc.resize(nk_, nb_local, nbs_local);
    grad_wfc.zero_out();

    E_one_ = 0.0;
    E_hartree_ = 0.0;
    E_xc_ = 0.0;
    E_entropy_ = 0.0;

    cached_h_one_diag_.resize(nk_);
    std::vector<double> h_one_diag(nbands_, 0.0);
    std::vector<double> vh_diag(nbands_, 0.0);
    std::vector<double> vx_diag(nbands_, 0.0);

    std::vector<TK> Hpsi_one(static_cast<size_t>(hpsi_alloc), TK(0));
    std::vector<TK> Hpsi_h(static_cast<size_t>(hpsi_alloc), TK(0));
    std::vector<TK> Hpsi_x(static_cast<size_t>(hpsi_alloc), TK(0));
    std::vector<TK> psi_dummy(1, TK(0));

    std::unique_ptr<psi::Psi<TK>> pg_one;
    std::unique_ptr<psi::Psi<TK>> pg_h;
    std::unique_ptr<psi::Psi<TK>> pg_x;
    if (print_orb_grad_decomp_)
    {
        pg_one = std::make_unique<psi::Psi<TK>>();
        pg_one->resize(nk_, nb_local, nbs_local);
        pg_one->zero_out();
        pg_h = std::make_unique<psi::Psi<TK>>();
        pg_h->resize(nk_, nb_local, nbs_local);
        pg_h->zero_out();
        pg_x = std::make_unique<psi::Psi<TK>>();
        pg_x->resize(nk_, nb_local, nbs_local);
        pg_x->zero_out();
    }

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* psi_k = psi_k_ptr_or_dummy(wfc_eval, ik, psi_dummy);

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
        if (exx_enabled_ && !is_mixed_exx)
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
        else if (exx_enabled_ && is_mixed_exx)
        {
            if (!mix_Hpsi_x_acc.empty())
            {
                std::copy(mix_Hpsi_x_acc[ik].begin(), mix_Hpsi_x_acc[ik].end(), Hpsi_x.begin());
            }
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
            if (is_mixed_exx)
            {
                E_xc_ += wk * mix_vx_E_acc[ik][ib] * 0.5;
            }
            else if (!rdmft_skip_occ_weight(gn))
            {
                E_xc_ += wk * gn * vx_diag[ib] * 0.5; // factor 1/2 for exchange
            }
        }

        // Occupation gradient: dE/dn_ik
        for (int ib = 0; ib < nbands_; ++ib)
        {
            double n = occ_flat[ik * nbands_ + ib];
            grad_occ[ik * nbands_ + ib] = wk * (h_one_diag[ib] + vh_diag[ib]);
            if (is_mixed_exx)
            {
                grad_occ[ik * nbands_ + ib] += wk * mix_vx_G_acc[ik][ib];
            }
            else
            {
                const double dgn = xc_func_.dg(n);
                grad_occ[ik * nbands_ + ib] += wk * dgn * vx_diag[ib];
            }
        }

        // Orbital gradient w.r.t. LCAO coefficients (before X transform), as the
        // Riesz representative for the ambient pairing Re ∑ conj(δC) · G (real TK)
        // / Re ∑ conj(δC_μ) G_μ (complex TK). For E_k ∝ n ⟨C|H|C⟩, ∂E/∂C* brings
        // a factor 2 relative to n H C alone; gamma-only real LCAO matches ∂E/∂C = 2 n H C.
        // Omit terms when n=0 or g(n)=0 so empty / inactive orbitals do not contribute.
        for (int ib_local = 0; ib_local < nb_local; ++ib_local)
        {
            int ib_global = ParaV_->local2global_col(ib_local);
            if (ib_global >= nbands_) continue;

            const double n = occ_flat[ik * nbands_ + ib_global];
            const double gn = xc_func_.g(n);
            const bool use_one_hart = !rdmft_skip_occ_weight(n);
            const bool use_exx = is_mixed_exx
                                     ? !rdmft_skip_occ_weight(n)
                                     : !rdmft_skip_occ_weight(gn);
            if (!use_one_hart && !use_exx)
            {
                continue;
            }

            TK* grad_ptr = &grad_wfc(ik, ib_local, 0);
            const TK* hone_ptr = &Hpsi_one[ib_local * nbs_local];
            const TK* hh_ptr = &Hpsi_h[ib_local * nbs_local];
            const TK* hx_ptr = &Hpsi_x[ib_local * nbs_local];

            TK* po = nullptr;
            TK* ph = nullptr;
            TK* px = nullptr;
            if (print_orb_grad_decomp_)
            {
                po = &(*pg_one)(ik, ib_local, 0);
                ph = &(*pg_h)(ik, ib_local, 0);
                px = &(*pg_x)(ik, ib_local, 0);
            }

            for (int mu = 0; mu < nbs_local; ++mu)
            {
                TK acc = TK(0);
                TK vone = TK(0);
                TK vh = TK(0);
                TK vx = TK(0);
                if (use_one_hart)
                {
                    vone = TK(n) * hone_ptr[mu];
                    vh = TK(n) * hh_ptr[mu];
                    acc += vone + vh;
                }
                if (use_exx)
                {
                    if (is_mixed_exx)
                    {
                        vx = hx_ptr[mu];
                        acc += vx;
                    }
                    else
                    {
                        vx = TK(gn) * hx_ptr[mu];
                        acc += vx;
                    }
                }
                grad_ptr[mu] = TK(2.0) * wk * acc;
                if (print_orb_grad_decomp_)
                {
                    po[mu] = TK(2.0) * wk * vone;
                    ph[mu] = TK(2.0) * wk * vh;
                    px[mu] = TK(2.0) * wk * vx;
                }
            }
        }
    }

    if (occ_entropy_gamma_ > 0.0 && xc_func_.type() == XCFunctionalType::HF)
    {
        for (int ik = 0; ik < nk_; ++ik)
        {
            const double wk = kv_->wk[ik];
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const double n = occ_flat[ik * nbands_ + ib];
                E_entropy_ += occ_entropy_gamma_ * wk * binary_entropy_f(n);
                grad_occ[ik * nbands_ + ib] += occ_entropy_gamma_ * wk * binary_entropy_dfdn(n);
            }
        }
    }

    // h_one_diag / vh_diag / vx_diag are already globally replicated by
    // compute_diagonal(); the per-rank accumulations above all produced the
    // same value, so NO MPI reduction of E_* is needed (that would double-
    // count). The same is true for grad_occ.

    E_ewald_ = pelec_->f_en.ewald_energy;
    E_total_ = E_one_ + E_hartree_ + E_xc_ + E_entropy_ + E_ewald_;

    if (print_orb_grad_decomp_)
    {
        const double res2 = orb_grad_decomp_residual_frob2(grad_wfc, *pg_one, *pg_h, *pg_x);
        const double nfull_c = ambient_frobenius_norm2_psi(grad_wfc);
        const double rel = std::sqrt(res2 / std::max(nfull_c, 1e-300));
        const double n1c = ambient_frobenius_norm2_psi(*pg_one);
        const double nhc = ambient_frobenius_norm2_psi(*pg_h);
        const double nxc = ambient_frobenius_norm2_psi(*pg_x);
        auto sqrtp = [](double x) { return std::sqrt(std::max(0.0, x)); };

        grad_C_to_X(*pg_one);
        grad_C_to_X(*pg_h);
        grad_C_to_X(*pg_x);

        const double n1x = ambient_frobenius_norm2_psi(*pg_one);
        const double nhx = ambient_frobenius_norm2_psi(*pg_h);
        const double nxx = ambient_frobenius_norm2_psi(*pg_x);

        auto proj_can_norm = [&](const psi::Psi<TK>& Gx) -> double {
            psi::Psi<TK> t(Gx);
            project_orbital_gradient(wfc, t);
            const double g2 = stiefel_canonical_inner_product(wfc, t, t);
            return sqrtp(g2);
        };
        const double p_one = proj_can_norm(*pg_one);
        const double p_hart = proj_can_norm(*pg_h);
        const double p_exx = proj_can_norm(*pg_x);

        // Full ambient gradient in X-space (same as returned grad_wfc after C→X).
        grad_C_to_X(grad_wfc);
        const double nfull_x = ambient_frobenius_norm2_psi(grad_wfc);
        psi::Psi<TK> gfull_tmp(grad_wfc);
        project_orbital_gradient(wfc, gfull_tmp);
        const double pfull = sqrtp(stiefel_canonical_inner_product(wfc, gfull_tmp, gfull_tmp));

        std::ostringstream os;
        os << "  RDMFT orb grad decomp (∂E/∂C → U^{-H} for X; EXX = RDMFT exchange / mixed Fock channel):\n";
        os << "    ‖G_one‖_F (C)=" << sqrtp(n1c) << "  ‖G_H‖_F (C)=" << sqrtp(nhc)
           << "  ‖G_EXX‖_F (C)=" << sqrtp(nxc) << "  ‖G_full‖_F (C)=" << sqrtp(nfull_c)
           << "  recon √(‖G_full-Σ‖²)/‖G_full‖=" << rel << "\n";
        os << "    ‖G_one‖_F (X)=" << sqrtp(n1x) << "  ‖G_H‖_F (X)=" << sqrtp(nhx)
           << "  ‖G_EXX‖_F (X)=" << sqrtp(nxx) << "  ‖G_full‖_F (X)=" << sqrtp(nfull_x) << "\n";
        os << "    ‖P(G_one)‖_can=" << p_one << "  ‖P(G_H)‖_can=" << p_hart << "  ‖P(G_EXX)‖_can=" << p_exx
           << "  ‖P(G_full)‖_can=" << pfull << "  (P = canonical Stiefel projection)\n";
        GlobalV::ofs_running << os.str() << std::flush;
    }
    else
    {
        // Transform the orbital gradient from C-space to X-space: G_X = U^{-H} G_C.
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

    // RDMFT operates in X-space. Convert X_k -> C_k for energy evaluation.
    psi::Psi<TK> wfc_C(wfc);
    wfc_X_to_C(wfc_C);
    const psi::Psi<TK>& wfc_eval = wfc_C;

    // Rebuild charge and Hartree (these always depend on occupations)
    build_charge(occ_flat, wfc_eval);

    HR_hartree_->set_zero();
    if (!op_hartree_)
    {
        op_hartree_ = std::make_unique<Veff_rdmft_local<TK, TR>>(
            hsk_hartree_.get(), kv_->kvec_d, pelec_->pot, HR_hartree_.get(), ucell_,
            orb_->cutoffs(), gd_, nspin_, charge_, rho_basis_, vloc_, sf_, "hartree");
    }
    op_hartree_->contributeHR();
    const bool is_mixed_exx = (xc_func_.type() == XCFunctionalType::GEO
                               || xc_func_.type() == XCFunctionalType::HybOpt
                               || xc_func_.type() == XCFunctionalType::CHF
                               || xc_func_.type() == XCFunctionalType::CGA);
    auto mixed_num_terms = [&]() -> int {
        switch (xc_func_.type())
        {
            case XCFunctionalType::GEO: return XCFunctional::num_geo_terms();
            case XCFunctionalType::HybOpt:
            case XCFunctionalType::CHF:
            case XCFunctionalType::CGA: return 2;
            default: return 0;
        }
    };
    auto mixed_term_coeff = [&](int t) -> double {
        switch (xc_func_.type())
        {
            case XCFunctionalType::GEO: return XCFunctional::geo_coef(t);
            case XCFunctionalType::HybOpt:
                return (t == 0) ? XCFunctional::hybopt_hf_weight() : XCFunctional::hybopt_power_weight();
            case XCFunctionalType::CHF:
                return (t == 0) ? XCFunctional::chf_hf_weight() : XCFunctional::chf_corr_weight();
            case XCFunctionalType::CGA:
                return (t == 0) ? XCFunctional::cga_hf_weight() : XCFunctional::cga_corr_weight();
            default: return 0.0;
        }
    };
    auto mixed_term_weight = [&](int t, double n) -> double {
        switch (xc_func_.type())
        {
            case XCFunctionalType::GEO: return xc_func_.pow_reg(n, XCFunctional::geo_alpha(t));
            case XCFunctionalType::HybOpt:
                return (t == 0) ? n : xc_func_.pow_reg(n, XCFunctional::hybopt_power_exponent());
            case XCFunctionalType::CHF:
                return (t == 0) ? n : xc_func_.chf_corr_term(n);
            case XCFunctionalType::CGA:
                return (t == 0) ? n : xc_func_.cga_corr_term(n);
            default: return 0.0;
        }
    };
    std::vector<std::vector<double>> mix_vx_E_acc(nk_, std::vector<double>(nbands_, 0.0));

    // Exchange must also be rebuilt because the modified DM gamma_xc depends
    // on occupations (and on orbitals, but those are fixed for compute_energy
    // use cases). Without this, line-search / finite-difference calls would
    // use a stale H_exx from the last compute() call.
#ifdef __EXX
    if (exx_enabled_ && !is_mixed_exx)
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

        }
    }
    else if (exx_enabled_ && is_mixed_exx)
    {
        const int n_terms = mixed_num_terms();
        const int nb_local_mix = wfc_eval.get_nbands();
        const int nbs_local_mix = wfc_eval.get_nbasis();
        const std::int64_t hpsi_alloc_mix
            = std::max<std::int64_t>(static_cast<std::int64_t>(nb_local_mix * nbs_local_mix), 1);
        std::vector<TK> Hpsi_x_term(static_cast<size_t>(hpsi_alloc_mix), TK(0));
        std::vector<double> vx_diag_term(nbands_, 0.0);
        std::vector<TK> psi_dummy_mix(1, TK(0));
        std::vector<double> term_weights(nk_ * nbands_, 0.0);
        for (int t = 0; t < n_terms; ++t)
        {
            const double coeff = mixed_term_coeff(t);
            for (int ik = 0; ik < nk_; ++ik)
            {
                for (int ib = 0; ib < nbands_; ++ib)
                {
                    const double n = occ_flat[ik * nbands_ + ib];
                    term_weights[ik * nbands_ + ib] = mixed_term_weight(t, n);
                }
            }

            HR_exx_->set_zero();
            std::vector<std::vector<TK>> DM_XC;
            build_DM_xc(occ_flat, wfc_eval, DM_XC, 0.0, &term_weights);

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
            }

            for (int ik = 0; ik < nk_; ++ik)
            {
                const TK* psi_k = psi_k_ptr_or_dummy(wfc_eval, ik, psi_dummy_mix);
                std::fill(Hpsi_x_term.begin(), Hpsi_x_term.end(), TK(0));
                std::fill(vx_diag_term.begin(), vx_diag_term.end(), 0.0);
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
                apply_Hk(hsk_exx_->get_hk(), psi_k, Hpsi_x_term.data());
                compute_diagonal(psi_k, Hpsi_x_term.data(), vx_diag_term.data(), ik);
                for (int ib = 0; ib < nbands_; ++ib)
                {
                    const double n = occ_flat[ik * nbands_ + ib];
                    const double p = mixed_term_weight(t, n);
                    mix_vx_E_acc[ik][ib] += coeff * p * vx_diag_term[ib];
                }
            }
        }
    }
#endif

    // Populate one-body diagonal cache if not valid
    if (!hone_cache_valid_)
    {
        cached_h_one_diag_.resize(nk_);
        const int nb_local = wfc_eval.get_nbands();
        const int nbs_local = wfc_eval.get_nbasis();
        const std::int64_t hpsi_alloc
            = std::max<std::int64_t>(static_cast<std::int64_t>(nb_local * nbs_local), 1);
        std::vector<TK> Hpsi_buf(static_cast<size_t>(hpsi_alloc), TK(0));
        std::vector<TK> psi_dummy(1, TK(0));

        for (int ik = 0; ik < nk_; ++ik)
        {
            cached_h_one_diag_[ik].assign(nbands_, 0.0);
            const TK* psi_k = psi_k_ptr_or_dummy(wfc_eval, ik, psi_dummy);

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
    const std::int64_t hpsi_alloc
        = std::max<std::int64_t>(static_cast<std::int64_t>(nb_local * nbs_local), 1);
    std::vector<TK> Hpsi_buf(static_cast<size_t>(hpsi_alloc), TK(0));
    std::vector<TK> psi_dummy(1, TK(0));

    E_one_ = 0.0;
    E_hartree_ = 0.0;
    E_xc_ = 0.0;
    E_entropy_ = 0.0;

    for (int ik = 0; ik < nk_; ++ik)
    {
        const TK* psi_k = psi_k_ptr_or_dummy(wfc_eval, ik, psi_dummy);

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
        if (exx_enabled_ && !is_mixed_exx)
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
            if (is_mixed_exx)
            {
                E_xc_ += wk * mix_vx_E_acc[ik][ib] * 0.5;
            }
            else if (!rdmft_skip_occ_weight(gn))
            {
                E_xc_ += wk * gn * vx_diag[ib] * 0.5;
            }
        }
    }

    if (occ_entropy_gamma_ > 0.0 && xc_func_.type() == XCFunctionalType::HF)
    {
        for (int ik = 0; ik < nk_; ++ik)
        {
            const double wk = kv_->wk[ik];
            for (int ib = 0; ib < nbands_; ++ib)
            {
                const double n = occ_flat[ik * nbands_ + ib];
                E_entropy_ += occ_entropy_gamma_ * wk * binary_entropy_f(n);
            }
        }
    }

    // As in compute(), the *_diag arrays are replicated across ranks so no MPI
    // reduction of E_* is needed here.

    E_ewald_ = pelec_->f_en.ewald_energy;
    E_total_ = E_one_ + E_hartree_ + E_xc_ + E_entropy_ + E_ewald_;

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

    if (!op_overlap_ || !hsk_overlap_)
    {
        ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
        throw std::runtime_error("RDMFT requires overlap matrix S_k to build X-space variables.");
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
            ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
            throw std::runtime_error("RDMFT X-space setup failed: overlap matrix S_k unavailable at ik="
                                     + std::to_string(ik));
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
                ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
                throw std::runtime_error("RDMFT X-space setup failed: Cholesky factorisation of S_k failed at ik="
                                         + std::to_string(ik) + " info=" + std::to_string(info));
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
                ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
                throw std::runtime_error("RDMFT X-space setup failed: triangular inverse of U_k failed at ik="
                                         + std::to_string(ik) + " info=" + std::to_string(info));
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
            ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
            throw std::runtime_error("RDMFT X-space setup failed: overlap matrix S_k unavailable at ik="
                                     + std::to_string(ik));
        }

        // Copy S_k into Uk_ for in-place Cholesky
        Uk_[ik].assign(SK, SK + nbasis * nbasis);

        // Upper Cholesky: S = U^H U
        int info = detail::potrf_upper(Uk_[ik].data(), nbasis);
        if (info != 0)
        {
            ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
            throw std::runtime_error("RDMFT X-space setup failed: Cholesky factorisation of S_k failed at ik="
                                     + std::to_string(ik) + " info=" + std::to_string(info));
        }

        // Compute U^{-1}
        Uk_inv_[ik] = Uk_[ik];
        info = detail::trtri_upper(Uk_inv_[ik].data(), nbasis);
        if (info != 0)
        {
            ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
            throw std::runtime_error("RDMFT X-space setup failed: triangular inverse of U_k failed at ik="
                                     + std::to_string(ik) + " info=" + std::to_string(info));
        }
    }
#endif

    cholesky_precomputed_ = true;
    GlobalV::ofs_running << "  Cholesky S = U^H U precomputed for all k-points; "
        << "using X_k = U_k C_k as Stiefel variable." << std::endl;

    ModuleBase::timer::end("RDMFT_EG", "precompute_cholesky_S");
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::wfc_C_to_X(psi::Psi<TK>& wfc)
{
    // X_k = U_k * C_k  for each k-point
    if (!cholesky_precomputed_)
    {
        throw std::runtime_error("wfc_C_to_X called before precompute_cholesky_S().");
    }

    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();

#ifdef __MPI
    const int nbasis_const = ParaV_->desc[2];
    const int nbands_const = ParaV_->desc_wfc[3];
    int nbasis = nbasis_const;
    int nbands = nbands_const;
    std::vector<TK> c_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* C = psi_k_ptr_or_dummy(wfc, ik, c_dummy);

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
    std::vector<TK> c_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* C = psi_k_ptr_or_dummy(wfc, ik, c_dummy);

        // X = U * C (in-place via trmm)
        detail::trmm_left_upper_notr(nbasis, nbands, Uk_[ik].data(), nbasis, C, nbasis);
    }
#endif
}

template <typename TK, typename TR>
void EnergyGradient<TK, TR>::wfc_X_to_C(psi::Psi<TK>& wfc)
{
    // C_k = U_k^{-1} * X_k  for each k-point
    if (!cholesky_precomputed_)
    {
        throw std::runtime_error("wfc_X_to_C called before precompute_cholesky_S().");
    }

    const int nb_local = wfc.get_nbands();
    const int nbs_local = wfc.get_nbasis();

#ifdef __MPI
    int nbasis = ParaV_->desc[2];
    int nbands = ParaV_->desc_wfc[3];
    std::vector<TK> x_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* X = psi_k_ptr_or_dummy(wfc, ik, x_dummy);

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
    std::vector<TK> x_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* X = psi_k_ptr_or_dummy(wfc, ik, x_dummy);

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
    if (!cholesky_precomputed_)
    {
        throw std::runtime_error("grad_C_to_X called before precompute_cholesky_S().");
    }

    const int nb_local = grad_wfc.get_nbands();
    const int nbs_local = grad_wfc.get_nbasis();

#ifdef __MPI
    int nbasis = ParaV_->desc[2];
    int nbands = ParaV_->desc_wfc[3];
    std::vector<TK> g_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* G = psi_k_ptr_or_dummy(grad_wfc, ik, g_dummy);

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
    std::vector<TK> g_dummy(1, TK(0));

    for (int ik = 0; ik < nk_; ++ik)
    {
        if (Uk_[ik].empty()) continue;

        TK* G = psi_k_ptr_or_dummy(grad_wfc, ik, g_dummy);

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
