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

    // ABACUS infrastructure (non-owning)
    const Parallel_Orbitals* ParaV_ = nullptr;
    const UnitCell* ucell_ = nullptr;
    const Grid_Driver* gd_ = nullptr;
    const K_Vectors* kv_ = nullptr;
    elecstate::ElecState* pelec_ = nullptr;
    const LCAO_Orbitals* orb_ = nullptr;
    const TwoCenterBundle* two_center_bundle_ = nullptr;
    Charge* charge_ = nullptr;
    const ModulePW::PW_Basis* rho_basis_ = nullptr;
    const ModuleBase::matrix* vloc_ = nullptr;
    const ModuleBase::ComplexMatrix* sf_ = nullptr;

    XCFunctional xc_func_{XCFunctionalType::HF};

    int nk_ = 0;
    int nbands_ = 0;
    int nbasis_local_ = 0;
    int nspin_ = 1;

    // Owned Hamiltonian containers
    hamilt::HContainer<TR>* HR_one_ = nullptr;
    hamilt::HContainer<TR>* HR_hartree_ = nullptr;
    hamilt::HContainer<TR>* HR_exx_ = nullptr;

    hamilt::HS_Matrix_K<TK>* hsk_one_ = nullptr;
    hamilt::HS_Matrix_K<TK>* hsk_hartree_ = nullptr;
    hamilt::HS_Matrix_K<TK>* hsk_exx_ = nullptr;

    hamilt::OperatorLCAO<TK, TR>* op_ekinetic_ = nullptr;
    hamilt::OperatorLCAO<TK, TR>* op_nonlocal_ = nullptr;
    hamilt::OperatorLCAO<TK, TR>* op_local_ = nullptr;
    hamilt::OperatorLCAO<TK, TR>* op_hartree_ = nullptr;
    hamilt::OperatorLCAO<TK, TR>* op_exx_ = nullptr;

#ifdef __EXX
    Exx_LRI<double>* exx_lri_d_ = nullptr;
    Exx_LRI<std::complex<double>>* exx_lri_c_ = nullptr;
    ModuleSymmetry::Symmetry_rotation symrot_exx_;
    bool exx_spacegroup_symmetry_ = false;
#endif

    Parallel_2D para_Eij_;

    // Energy components
    double E_one_ = 0.0;
    double E_hartree_ = 0.0;
    double E_xc_ = 0.0;
    double E_ewald_ = 0.0;
    double E_total_ = 0.0;
    double etxc_ = 0.0;
    double vtxc_ = 0.0;

    bool ion_initialized_ = false;
    bool exx_enabled_ = false; // true when RDMFT has its own EXX infrastructure

    // Cache for one-body diagonals (valid when orbitals don't change)
    bool hone_cache_valid_ = false;
    std::vector<std::vector<double>> cached_h_one_diag_; // [nk][nbands]
};

} // namespace rdmft

#endif // RDMFT_ENERGY_GRADIENT_H
