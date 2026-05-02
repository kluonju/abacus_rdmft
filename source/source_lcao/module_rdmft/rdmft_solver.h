#ifndef RDMFT_SOLVER_H
#define RDMFT_SOLVER_H

#include "rdmft_type.h"
#include "rdmft_xc_functional.h"
#include "rdmft_occupation.h"
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
    /// Dispatches to optimize_orbitals_riemannian_bb (default) or
    /// optimize_orbitals_simple based on config_.orb_strategy.
    OptResult optimize_orbitals(const std::vector<double>& occ_flat,
                                 psi::Psi<TK>& wfc);

    /// Verify gradient consistency by finite differences
    bool check_gradient_consistency(const std::vector<double>& occ_flat,
                                     const psi::Psi<TK>& wfc,
                                     double epsilon = 1e-5,
                                     double tolerance = 1e-4);

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

    /// RiemannianBB strategy (default): the existing implementation that
    /// combines Riemannian gradient + BB1 spectral step + non-monotone
    /// Armijo + Wen-Yin trust radius + direction blending. Named after
    /// "Riemannian SPG" in the source comments (Iannazzo-Porcelli, IMA
    /// JNA 38 (2018) 495).
    OptResult optimize_orbitals_riemannian_bb(const std::vector<double>& occ_flat,
                                              psi::Psi<TK>& wfc);

    /// Simple strategy: textbook Riemannian SD/CG + monotone Armijo
    /// (Absil-Mahony-Sepulchre 2008, §4.2). No BB step, no non-monotone
    /// history, no trust radius, no descent guard. Only sd / cg are
    /// honoured (lbfgs / adam fall back to cg with a warning).
    OptResult optimize_orbitals_simple(const std::vector<double>& occ_flat,
                                       psi::Psi<TK>& wfc);

    /// Single step of orbital optimization on Stiefel manifold
    void orbital_step(const std::vector<double>& occ_flat, psi::Psi<TK>& wfc);

    /// Compute norm of the occupation gradient
    double occ_grad_norm(const std::vector<double>& grad) const;

    /// Compute norm of the orbital Riemannian gradient
    double orb_grad_norm(const psi::Psi<TK>& rgrad) const;

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

    OptResult last_result_;

    /// Cleared in solve(); set after logging once so the first occ-gradient evaluation point is printed.
    bool logged_initial_occ_for_first_occ_gradient_ = false;
};

} // namespace rdmft

#endif // RDMFT_SOLVER_H
