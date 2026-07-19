#ifndef RDMFT_OCC_OPTIMIZER_H
#define RDMFT_OCC_OPTIMIZER_H

#include "rdmft_backend.h"
#include "rdmft_params.h"

#include <vector>

//! \file rdmft_occ_optimizer.h
//! \brief Occupation-number inner optimiser (SPG2 and EBI).
//!
//! Faithful C++ port of the two occupation blocks in the Quantum ESPRESSO
//! RDMFT reference:
//!
//!   * SPG2 (qe-rdmft/rdmft_spg.f90): the canonical Birgin-Martinez-Raydan
//!     spectral projected gradient with the Euclidean L2 proximal projector
//!     P_w, the BB spectral steplength and a straight-chord Armijo/Wolfe line
//!     search.  This replaces the previous (incorrect) augmented-Lagrangian /
//!     projected-gradient occupation code.
//!
//!   * EBI (qe-rdmft/rdmft_ebi.f90): the explicit-by-implicit erf
//!     parameterisation n = (erf(x+mu)+1)/2 with an implicit chemical
//!     potential mu enforcing the electron count, optimised by steepest
//!     descent in x-space (EBI@GD).

namespace rdmft
{

//! Result of one occupation block.
struct OccBlockResult
{
    bool converged = false;
    double energy = 0.0;
    int iterations = 0;
    double grad_norm = 0.0; //!< ||g1||_inf (SPG2) or ||grad_x|| (EBI)
};

class OccOptimizer
{
  public:
    //! Optimise the occupations in place (orbitals held fixed inside backend).
    //! etot is updated to the final block energy.
    OccBlockResult run(RdmftBackend& backend, const RdmftParams& params,
                       std::vector<double>& occ, double& etot);

  private:
    OccBlockResult run_spg2(RdmftBackend& backend, const RdmftParams& params,
                            std::vector<double>& occ, double& etot);
    OccBlockResult run_ebi(RdmftBackend& backend, const RdmftParams& params,
                           std::vector<double>& occ, double& etot);
    OccBlockResult run_bgd(RdmftBackend& backend, const RdmftParams& params,
                           std::vector<double>& occ, double& etot);
};

} // namespace rdmft

#endif // RDMFT_OCC_OPTIMIZER_H
