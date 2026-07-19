#ifndef RDMFT_BACKEND_H
#define RDMFT_BACKEND_H

#include "rdmft_occ_constraints.h"

#include <vector>

//! \file rdmft_backend.h
//! \brief Basis-independent energy/gradient oracle for the RDMFT optimisers.
//!
//! This is the single seam between the optimisation core (occupation and
//! orbital blocks) and the basis-specific Hamiltonian machinery.  A concrete
//! backend (plane-wave or LCAO) owns the natural orbitals and knows how to
//! build one-body, Hartree and exchange terms; it exposes to the optimisers
//! only what they need:
//!
//!   * the occupation index space and electron-number targets;
//!   * the RDMFT total energy for a given occupation vector (using the
//!     backend's current orbitals);
//!   * the physical occupation gradient dE/dn (already weighted by w_k);
//!   * Stiefel-manifold primitives on the packed orbital tangent space.
//!
//! Occupation vectors are laid out as occ[ib + ik*nbnd] (band fastest), the
//! same convention as OccConstraints.  Orbital tangent vectors are opaque
//! flat double buffers whose interpretation (real/complex, k-blocking,
//! gamma-trick, AO overlap S) is entirely internal to the backend.  This
//! keeps the SPG2 / EBI / Stiefel algorithms identical for gamma-only
//! (double), collinear (nspin=2) and non-collinear (nspin=4) runs and for
//! both the plane-wave and LCAO bases.

namespace rdmft
{

class RdmftBackend
{
  public:
    virtual ~RdmftBackend() = default;

    //! Occupation index space + electron-number constraints.
    virtual const OccConstraints& occ_constraints() const = 0;

    // --- occupation oracle (orbitals held fixed) ----------------------------

    //! RDMFT total energy for the given occupations and the backend's
    //! current natural orbitals.
    virtual double total_energy(const std::vector<double>& occ) = 0;

    //! Physical occupation gradient g_ik = dE/dn_ik (includes the BZ weight
    //! w_k, matching the QE rdmft_grad_n convention).
    virtual void grad_occ(const std::vector<double>& occ, std::vector<double>& grad) = 0;

    // --- orbital (Stiefel) oracle (occupations held fixed) ------------------

    //! Whether this backend supports orbital optimisation.  A pure
    //! occupation-only backend (used in unit tests or frozen-orbital runs)
    //! may return false.
    virtual bool has_orbital_optimization() const { return false; }

    //! Length of the packed orbital tangent representation.
    virtual int orb_dim() const { return 0; }

    //! Riemannian (tangent) gradient of the energy at the current orbitals
    //! for the given occupations.  Fills gR (length orb_dim()) and returns
    //! its squared Stiefel norm sum_k <gR_k, gR_k>_S.
    virtual double riemannian_gradient(const std::vector<double>& occ,
                                       std::vector<double>& gR)
    {
        (void)occ;
        (void)gR;
        return 0.0;
    }

    //! Stiefel inner product of two packed tangent vectors (summed over k).
    virtual double orb_inner(const std::vector<double>& a, const std::vector<double>& b)
    {
        (void)a;
        (void)b;
        return 0.0;
    }

    //! Project a packed vector onto the tangent space at the current orbitals.
    virtual void orb_project_tangent(std::vector<double>& v) { (void)v; }

    //! Retract the current orbitals along the packed tangent direction dir by
    //! step alpha (C <- R_C(alpha*dir)), preserving C^H S C = I.
    virtual void orb_retract(const std::vector<double>& dir, double alpha)
    {
        (void)dir;
        (void)alpha;
    }

    //! Save / restore the current orbitals (used by orbital line searches).
    virtual void orb_save() {}
    virtual void orb_restore() {}
};

} // namespace rdmft

#endif // RDMFT_BACKEND_H
