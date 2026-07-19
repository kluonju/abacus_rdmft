#ifndef ESOLVER_RDMFT_PW_H
#define ESOLVER_RDMFT_PW_H

#include "esolver_ks_pw.h"

//! \file esolver_rdmft_pw.h
//! \brief Plane-wave RDMFT energy solver (esolver_type = rdmft, basis_type = pw).
//!
//! ESolver_RDMFT_PW performs a Kohn-Sham (hybrid) SCF to obtain the starting
//! natural orbitals and occupations, then minimises the RDMFT energy over the
//! occupation numbers through the basis-independent optimisation core
//! (source_rdmft) using a plane-wave RdmftBackend whose exact-exchange term is
//! evaluated with the existing ACE (Adaptively Compressed Exchange) machinery,
//! preserving its k-pool / G-vector MPI parallelism.
//!
//! Orbitals are held at the converged hybrid-KS solution (a valid RDMFT mode
//! that reproduces Hartree-Fock for the HF functional, where the KS-hybrid
//! orbitals are already stationary).

namespace ModuleESolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class ESolver_RDMFT_PW : public ESolver_KS_PW<T, Device>
{
  public:
    ESolver_RDMFT_PW();
    ~ESolver_RDMFT_PW();

  protected:
    //! Run the modular RDMFT occupation optimisation after the KS SCF.
    void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;
};

} // namespace ModuleESolver

#endif // ESOLVER_RDMFT_PW_H
