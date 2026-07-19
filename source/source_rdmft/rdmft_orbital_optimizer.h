#ifndef RDMFT_ORBITAL_OPTIMIZER_H
#define RDMFT_ORBITAL_OPTIMIZER_H

#include "rdmft_backend.h"
#include "rdmft_params.h"

#include <vector>

//! \file rdmft_orbital_optimizer.h
//! \brief Natural-orbital inner optimiser on the Stiefel manifold.
//!
//! Port of the orbital block of the Quantum ESPRESSO RDMFT reference
//! (qe-rdmft/rdmft_solver.f90::rdmft_orbital_step_joint together with the
//! Stiefel primitives in rdmft_stiefel.f90).  The orbitals live on the
//! generalised Stiefel manifold C^H S C = I; all geometry (tangent
//! projection, retraction via a Cholesky/polar factor, and the metric inner
//! product, including the AO overlap S for LCAO or S = I for plane waves) is
//! delegated to the backend, so the same SD / CG / L-BFGS driver serves both
//! bases and every spin case.

namespace rdmft_core
{

struct OrbBlockResult
{
    bool converged = false;
    double energy = 0.0;
    int iterations = 0;
    double grad_norm = 0.0; //!< ||G_R|| Riemannian gradient norm
};

class OrbitalOptimizer
{
  public:
    //! Optimise the backend's natural orbitals in place (occupations fixed).
    OrbBlockResult run(RdmftBackend& backend, const RdmftParams& params,
                       const std::vector<double>& occ, double& etot);
};

} // namespace rdmft_core

#endif // RDMFT_ORBITAL_OPTIMIZER_H
