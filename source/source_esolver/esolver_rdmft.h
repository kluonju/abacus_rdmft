#ifndef SOURCE_ESOLVER_ESOLVER_RDMFT_H
#define SOURCE_ESOLVER_ESOLVER_RDMFT_H

#include "esolver_ks_lcao.h"

#ifdef __RDMFT
#include "source_lcao/module_rdmft/rdmft_energy_gradient.h"
#endif

//! \file esolver_rdmft.h
//! \brief New, modular RDMFT energy solver (selected by esolver_type = "rdmft").
//!
//! ESolver_RDMFT is a dedicated ESolver for Reduced Density Matrix Functional
//! Theory built on top of the basis-independent optimisation core in
//! source_rdmft/.  A Kohn-Sham SCF run first provides the starting natural
//! orbitals and occupations; the solver then minimises the RDMFT energy over
//! occupation numbers (SPG2 / EBI) and natural orbitals (Stiefel) through the
//! abstract rdmft::RdmftBackend oracle.
//!
//! The LCAO specialisation reuses the validated rdmft::EnergyGradient physics
//! engine via ModuleESolver::RdmftBackendLCAO.  A plane-wave backend can be
//! plugged into the same driver by implementing rdmft::RdmftBackend for the
//! PW Hamiltonian; the factory dispatch below selects the appropriate class
//! for gamma-only (double), collinear (nspin=2) and non-collinear (nspin=4).

namespace ModuleESolver
{

template <typename TK, typename TR>
class ESolver_RDMFT : public ESolver_KS_LCAO<TK, TR>
{
  public:
    ESolver_RDMFT();
    ~ESolver_RDMFT();

  protected:
    //! After the Kohn-Sham SCF finishes, run the modular RDMFT optimisation.
    void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;

#ifdef __RDMFT
    //! Energy/gradient engine owned by this esolver (independent of the legacy
    //! embedded RDMFT path in ESolver_KS_LCAO).
    rdmft::EnergyGradient<TK, TR> eg_;
    bool eg_ready_ = false;
#endif
};

} // namespace ModuleESolver

#endif // SOURCE_ESOLVER_ESOLVER_RDMFT_H
