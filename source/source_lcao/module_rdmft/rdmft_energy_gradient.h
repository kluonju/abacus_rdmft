#ifndef RDMFT_ENERGY_GRADIENT_H
#define RDMFT_ENERGY_GRADIENT_H

#include "rdmft_type.h"
#include "rdmft_xc_functional.h"

#include "source_io/module_parameter/parameter.h"
#include "source_psi/psi.h"
#include "source_base/matrix.h"
#include "source_base/parallel_2d.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_basis/module_ao/ORB_read.h"
#include "source_basis/module_nao/two_center_bundle.h"
#include "source_estate/elecstate.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/hs_matrix_k.hpp"
#include "source_lcao/module_operator_lcao/operator_lcao.h"

#ifdef __EXX
#include "source_lcao/module_ri/Exx_LRI.h"
#include "source_lcao/module_ri/module_exx_symmetry/symmetry_rotation.h"
#endif

#include <vector>
#include <complex>
#include <memory>

namespace rdmft
{

/// Computes the RDMFT total energy and its gradients w.r.t. natural occupation
/// numbers and orbital coefficients.
///
/// Template parameters:
///   TK = double (gamma-only) or complex<double> (multi-k)
///   TR = double (nspin <= 2) or complex<double> (nspin = 4)
template <typename TK, typename TR>
class EnergyGradient
{
  public:
    EnergyGradient() = default;
    ~EnergyGradient();

    /// Initialize with ABACUS infrastructure pointers
    void init(const Parallel_Orbitals* ParaV_in,
              const UnitCell* ucell_in,
              const Grid_Driver* gd_in,
              const K_Vectors* kv_in,
              elecstate::ElecState* pelec_in,
              const LCAO_Orbitals* orb_in,
              const TwoCenterBundle* two_center_bundle_in,
              const XCFunctional& xc_func_in);

    /// HF-only occupation entropy prefactor γ (set from INPUT / RDMFTConfig each solve).
    void set_occ_entropy_gamma(double gamma) { occ_entropy_gamma_ = gamma; }

    /// Update ion-step quantities (HR for kinetic/local/nonlocal, EXX ions)
    void update_ion(const UnitCell& ucell,
                    const ModulePW::PW_Basis& rho_basis,
                    const ModuleBase::matrix& vloc,
                    const ModuleBase::ComplexMatrix& sf);

    /// Compute total energy and gradients for given occupations and orbitals.
    /// occ_flat: occupation numbers [nk * nbands], flattened
    /// wfc: natural orbital coefficients
    /// grad_occ: [out] dE/dn for each (ik, ib)
    /// grad_wfc: [out] Euclidean dE/dC* for each (ik, ib, mu)
    /// Returns total energy
    double compute(const std::vector<double>& occ_flat,
                   const psi::Psi<TK>& wfc,
                   std::vector<double>& grad_occ,
                   psi::Psi<TK>& grad_wfc);

    /// Compute energy only (no gradients), useful for line search.
    /// When orbitals are unchanged, uses cached one-body diagonals.
    double compute_energy(const std::vector<double>& occ_flat,
                          const psi::Psi<TK>& wfc);

    /// Invalidate the one-body diagonal cache (call when orbitals change)
    void invalidate_hone_cache() { hone_cache_valid_ = false; }

    /// Get individual energy components from last compute()
    double E_one_body() const { return E_one_; }
    double E_hartree() const { return E_hartree_; }
    double E_xc() const { return E_xc_; }
    double E_ewald() const { return E_ewald_; }
    double E_total() const { return E_total_; }
    /// Binary-entropy regularization energy γ·Σ w_k f(n) from last compute / compute_energy.
    double E_entropy() const { return E_entropy_; }

    int nk() const { return nk_; }
    int nbands() const { return nbands_; }
    int nbasis_local() const { return nbasis_local_; }

  private:
    /// Build charge density from occupations and orbitals
    void build_charge(const std::vector<double>& occ_flat,
                      const psi::Psi<TK>& wfc);

    /// Build the modified DM for exchange: gamma_xc = sum_i w_k g(n_ik) |phi_i><phi_i|
    void build_DM_xc(const std::vector<double>& occ_flat,
                     const psi::Psi<TK>& wfc,
                     std::vector<std::vector<TK>>& DM_XC);

    /// Compute one-body Hamiltonian * wfc and diagonal elements
    void compute_one_body(const psi::Psi<TK>& wfc,
                          std::vector<double>& h_diag,
                          psi::Psi<TK>& H_wfc);

    /// Compute Hartree potential and its action on wfc
    void compute_hartree(const psi::Psi<TK>& wfc,
                         std::vector<double>& vh_diag,
                         psi::Psi<TK>& VH_wfc);

    /// Compute exchange from modified DM and its action on wfc
    void compute_exchange(const std::vector<double>& occ_flat,
                          const psi::Psi<TK>& wfc,
                          std::vector<double>& vx_diag,
                          psi::Psi<TK>& VX_wfc);

    /// Helper: Hk * psi for a given k-point
    void apply_Hk(const TK* HK, const TK* psi_k, TK* Hpsi_k) const;

    /// Helper: diagonal elements <psi_i|H|psi_i>
    void compute_diagonal(const TK* psi_k, const TK* Hpsi_k,
                          double* diag, int ik) const;

  public:
    /// Project orbital gradient onto the tangent space of the Stiefel manifold
    /// with overlap matrix S (generalised Stiefel: C^H S C = I).
    ///   G_R = G - C * sym(C^H S G)
    /// Derivation (canonical/trace metric): seek G_R = G - C*K such that
    /// C^H S G_R is skew-Hermitian.  Since C^H S C = I,
    ///   C^H S G_R = C^H S G - K,
    /// and the skew-Hermitian condition gives K = sym(C^H S G).
    /// When S = I this reduces to G_R = G - C * sym(C^H G) (standard Stiefel).
    /// At any S-orthonormal critical point (e.g. KS eigenstates), G_R == 0,
    /// allowing the orbital inner loop to detect convergence immediately.
    void project_orbital_gradient(const psi::Psi<TK>& wfc,
                                   psi::Psi<TK>& grad_wfc);

    /// S-orthonormalise wfc along the direction grad_wfc with step -alpha:
    ///   C <- C - alpha * G
    ///   C <- C * M^{-1/2}   where M = C^H S C
    /// Ensures the new orbitals lie on the generalised Stiefel manifold
    /// (C^H S C = I) to machine precision.
    void retract_orbitals(psi::Psi<TK>& wfc,
                          const psi::Psi<TK>& grad_wfc,
                          double alpha);

    /// Return pointer to overlap matrix at k-point ik (column-major).
    /// Rebuilt lazily the first time it is requested per ion step.
    const TK* get_SK(int ik);

    /// Compute the Euclidean inner product Re Tr(X^H Y) summed over k-points.
    /// In the X-variable formulation this is the canonical metric used by the
    /// solver. Performs an MPI Allreduce internally.
    double s_inner_product(const psi::Psi<TK>& X, const psi::Psi<TK>& Y);

    /// Per k-point Stiefel Gram residual ||X_k^H X_k - I||_F.
    /// frob_per_ik is resized to nk_.
    void stiefel_gram_residual_frobenius_per_k(const psi::Psi<TK>& wfc,
                                                std::vector<double>& frob_per_ik);

    // ================================================================
    // Cholesky-based S^{1/2} variable transformation
    //
    // RDMFT operates only in X-space, with X_k = U_k C_k and S_k = U_k^H U_k.
    // After precompute_cholesky_S(), all manifold operations assume the
    // orbital variable is X_k and use the standard Stiefel metric X_k^H X_k = I.
    // ================================================================

    /// Precompute the Cholesky factorisation S_k = U_k^H U_k for each
    /// k-point and store U_k, U_k^{-1}. All subsequent manifold operations
    /// assume wfc contains X_k rather than C_k.
    void precompute_cholesky_S();

    /// Transform wfc in-place: C_k -> X_k = U_k C_k.
    /// Requires precompute_cholesky_S() to have been called.
    void wfc_C_to_X(psi::Psi<TK>& wfc);

    /// Transform wfc in-place: X_k -> C_k = U_k^{-1} X_k.
    /// Requires precompute_cholesky_S() to have been called.
    void wfc_X_to_C(psi::Psi<TK>& wfc);

    /// Transform Euclidean gradient in-place:
    ///   G_X = U_k^{-H} G_C   (chain rule from C = U^{-1} X)
    /// Requires precompute_cholesky_S() to have been called.
    void grad_C_to_X(psi::Psi<TK>& grad_wfc);

    // ABACUS infrastructure (non-owning)
    const Parallel_Orbitals* ParaV_ = nullptr;
    const UnitCell* ucell_ = nullptr;
    const Grid_Driver* gd_ = nullptr;
    const K_Vectors* kv_ = nullptr;
    elecstate::ElecState* pelec_ = nullptr;
    const LCAO_Orbitals* orb_ = nullptr;
    const TwoCenterBundle* two_center_bundle_ = nullptr;
    /// Non-owning alias into pelec_->charge. Lifetime is governed by pelec_.
    Charge* charge_ = nullptr;
    const ModulePW::PW_Basis* rho_basis_ = nullptr;
    const ModuleBase::matrix* vloc_ = nullptr;
    const ModuleBase::ComplexMatrix* sf_ = nullptr;

    XCFunctional xc_func_{XCFunctionalType::HF};
    /// HF occupation entropy prefactor; 0 disables (see set_occ_entropy_gamma).
    double occ_entropy_gamma_{0.0};

    int nk_ = 0;
    int nbands_ = 0;
    int nbasis_local_ = 0;
    int nspin_ = 1;

    // Owned Hamiltonian containers
    std::unique_ptr<hamilt::HContainer<TR>> HR_one_;
    std::unique_ptr<hamilt::HContainer<TR>> HR_hartree_;
    std::unique_ptr<hamilt::HContainer<TR>> HR_exx_;
    std::unique_ptr<hamilt::HContainer<TR>> SR_;

    std::unique_ptr<hamilt::HS_Matrix_K<TK>> hsk_one_;
    std::unique_ptr<hamilt::HS_Matrix_K<TK>> hsk_hartree_;
    std::unique_ptr<hamilt::HS_Matrix_K<TK>> hsk_exx_;
    std::unique_ptr<hamilt::HS_Matrix_K<TK>> hsk_overlap_;

    std::unique_ptr<hamilt::OperatorLCAO<TK, TR>> op_ekinetic_;
    std::unique_ptr<hamilt::OperatorLCAO<TK, TR>> op_nonlocal_;
    std::unique_ptr<hamilt::OperatorLCAO<TK, TR>> op_local_;
    std::unique_ptr<hamilt::OperatorLCAO<TK, TR>> op_hartree_;
    std::unique_ptr<hamilt::OperatorLCAO<TK, TR>> op_exx_;
    std::unique_ptr<hamilt::OperatorLCAO<TK, TR>> op_overlap_;

#ifdef __EXX
    std::unique_ptr<Exx_LRI<double>> exx_lri_d_;
    std::unique_ptr<Exx_LRI<std::complex<double>>> exx_lri_c_;
    ModuleSymmetry::Symmetry_rotation symrot_exx_;
    bool exx_spacegroup_symmetry_ = false;
#endif

    Parallel_2D para_Eij_;

    // Energy components
    double E_one_ = 0.0;
    double E_hartree_ = 0.0;
    double E_xc_ = 0.0;
    double E_ewald_ = 0.0;
    double E_entropy_ = 0.0;
    double E_total_ = 0.0;
    double etxc_ = 0.0;
    double vtxc_ = 0.0;

    bool ion_initialized_ = false;
    bool exx_enabled_ = false; // true when RDMFT has its own EXX infrastructure

    // Cache for one-body diagonals (valid when orbitals don't change)
    bool hone_cache_valid_ = false;
    std::vector<std::vector<double>> cached_h_one_diag_; // [nk][nbands]

    // Cholesky S^{1/2} variable transformation data used by the mandatory
    // X-space formulation.
    bool cholesky_precomputed_ = false;
    /// U_k (upper triangular Cholesky factor of S_k) for each k-point.
    /// Stored in column-major format. Only the upper triangle is meaningful;
    /// the lower triangle may contain arbitrary values after the factorisation.
    /// Size is nbasis_local * nbasis_local (non-MPI) or ParaV_->nloc (MPI).
    std::vector<std::vector<TK>> Uk_;
    /// U_k^{-1} for each k-point, same upper-triangular column-major layout.
    std::vector<std::vector<TK>> Uk_inv_;
};

} // namespace rdmft

#endif // RDMFT_ENERGY_GRADIENT_H
