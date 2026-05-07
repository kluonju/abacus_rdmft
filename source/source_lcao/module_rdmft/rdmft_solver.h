//==========================================================
// Author: Kai Luo
// Email: kluo@njust.edu.cn
// DATE : April 2026
//==========================================================
#ifndef RDMFT_SOLVER_H
#define RDMFT_SOLVER_H

#include "rdmft_type.h"
#include "rdmft_xc_functional.h"
#include "rdmft_occupation.h"
#include "rdmft_stiefel.h"
#include "rdmft_optimizer.h"
#include "rdmft_energy_gradient.h"

#include "source_psi/psi.h"
#include "source_base/matrix.h"
#include "source_cell/klist.h"

#include <vector>
#include <functional>
#include <cmath>
#include <string>
#include <memory>

namespace rdmft
{

/// High-level RDMFT solver.
/// Supports alternating optimization (occupation then orbitals) and
/// joint (product-manifold) optimization (orbitals and occupations are
/// packed together and stepped simultaneously).
template <typename TK, typename TR>
class RDMFTSolver
{
  public:
    RDMFTSolver() = default;
    ~RDMFTSolver() = default;

    /// Initialize the solver with configuration
    void init(const RDMFTConfig& config,
              EnergyGradient<TK, TR>& energy_grad,
              const K_Vectors* kv,
              int nbands,
              double n_electrons,
              const RDMFTNelectronTargetMeta& nelec_meta = RDMFTNelectronTargetMeta());

    /// Run the full RDMFT optimization.
    /// On input: occ_flat contains initial occupations, wfc contains initial orbitals.
    /// On output: occ_flat and wfc contain optimized values.
    /// Returns the final total energy.
    double solve(std::vector<double>& occ_flat, psi::Psi<TK>& wfc);

    /// Run occupation optimization only (orbitals fixed).
    OptResult optimize_occupations(std::vector<double>& occ_flat,
                                    const psi::Psi<TK>& wfc);

    /// Run orbital optimization only (occupations fixed).
    OptResult optimize_orbitals(const std::vector<double>& occ_flat,
                                 psi::Psi<TK>& wfc);

    /// Verify gradient consistency by finite differences
    bool check_gradient_consistency(const std::vector<double>& occ_flat,
                                     const psi::Psi<TK>& wfc,
                                     double epsilon = 1e-5,
                                     double tolerance = 1e-4);

    /// ELK `rdmeval` analogue: log ε_ik = ∂E/∂n at probe n=0.5 per state (expensive).
    void print_elk_style_evalsv(const std::vector<double>& occ_flat, const psi::Psi<TK>& wfc);

    /// Get the config
    const RDMFTConfig& config() const { return config_; }

    /// Get last optimization result
    OptResult last_result() const { return last_result_; }

  private:
    /// Alternating optimization strategy
    double solve_alternating(std::vector<double>& occ_flat, psi::Psi<TK>& wfc);

    /// Joint (product-manifold) optimisation strategy. Orbitals and
    /// occupation parameters are packed into one descent direction and
    /// updated simultaneously at every outer iteration.
    double solve_joint(std::vector<double>& occ_flat, psi::Psi<TK>& wfc);

    /// Single step of orbital optimization on Stiefel manifold
    void orbital_step(const std::vector<double>& occ_flat, psi::Psi<TK>& wfc);

    /// Compute norm of the occupation gradient
    double occ_grad_norm(const std::vector<double>& grad) const;

    /// ‖G‖_can = sqrt(⟨G,G⟩_can) at Stiefel point `wfc_X` (see
    /// `EnergyGradient::stiefel_canonical_inner_product`). Used for
    /// rdmft_orb_grad_tol and logging.
    double orb_grad_norm(const psi::Psi<TK>& wfc_X, const psi::Psi<TK>& rgrad) const;

    RDMFTConfig config_;
    EnergyGradient<TK, TR>* energy_grad_ = nullptr;
    const K_Vectors* kv_ = nullptr;
    int nk_ = 0;
    int nbands_ = 0;
    double n_electrons_ = 0.0;
    RDMFTNelectronTargetMeta nelec_meta_{};

    // Occupation sub-components
    std::unique_ptr<OccupationParam> occ_param_;
    std::unique_ptr<OccupationConstraint> occ_constraint_;
    std::unique_ptr<EuclideanOptimizer> occ_optimizer_;
    /// Sigma-shift parameterization (active when config_.occ_param == SigmaShift).
    std::unique_ptr<SigmaShiftOccParam> sigma_shift_param_;

    // Orbital sub-components
    std::unique_ptr<EuclideanOptimizer> orb_optimizer_;

    // Stiefel manifold per k-point
    std::vector<StiefelManifold<TK>> stiefel_;

    OptResult last_result_;
};

} // namespace rdmft

#endif // RDMFT_SOLVER_H
