#ifndef RDMFT_DRIVER_H
#define RDMFT_DRIVER_H

#include "rdmft_backend.h"
#include "rdmft_params.h"

#include <vector>

//! \file rdmft_driver.h
//! \brief Outer RDMFT driver alternating occupation and orbital blocks.
//!
//! Mirrors qe-rdmft/rdmft_solver.f90::rdmft_run_alternating.  The driver is
//! basis independent: it only talks to the abstract RdmftBackend oracle and
//! the two inner optimisers, so a single implementation drives the
//! plane-wave and LCAO ESolver_RDMFT backends for gamma-only, collinear and
//! non-collinear spin.

namespace rdmft
{

struct DriverResult
{
    bool converged = false;
    double energy = 0.0;
    int outer_iterations = 0;
};

class RdmftDriver
{
  public:
    //! Minimise the RDMFT energy in place over occ (occupations) and the
    //! backend's natural orbitals.  Returns the final energy in etot.
    DriverResult solve(RdmftBackend& backend, const RdmftParams& params,
                       std::vector<double>& occ, double& etot);
};

} // namespace rdmft

#endif // RDMFT_DRIVER_H
