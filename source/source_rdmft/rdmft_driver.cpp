#include "rdmft_driver.h"

#include "rdmft_occ_optimizer.h"
#include "rdmft_orbital_optimizer.h"

#include <cmath>

namespace rdmft
{

DriverResult RdmftDriver::solve(RdmftBackend& backend, const RdmftParams& params,
                                std::vector<double>& occ, double& etot)
{
    DriverResult res;
    const OccConstraints& con = backend.occ_constraints();

    // Apply the occupation-init mode to the KS-seeded vector, then project it
    // onto the feasible set (box + electron-number constraint).
    con.initialize_occupations(occ, static_cast<int>(params.occ_init_mode), params.occ_init_perturb,
                               params.occ_init_nbands_top);
    etot = backend.total_energy(occ);

    OccOptimizer occ_opt;
    OrbitalOptimizer orb_opt;

    const bool do_occ = (params.strategy == SolverStrategy::Alternating
                         || params.strategy == SolverStrategy::OccOnly);
    const bool do_orb = (params.strategy == SolverStrategy::Alternating
                         || params.strategy == SolverStrategy::OrbOnly)
                        && backend.has_orbital_optimization();

    double e_prev = etot;
    for (int outer = 0; outer < params.outer_maxiter; ++outer)
    {
        res.outer_iterations = outer + 1;
        bool occ_conv = true;
        bool orb_conv = true;

        if (do_occ)
        {
            OccBlockResult r = occ_opt.run(backend, params, occ, etot);
            occ_conv = r.converged;
        }
        if (do_orb)
        {
            OrbBlockResult r = orb_opt.run(backend, params, occ, etot);
            orb_conv = r.converged;
        }

        const double dE = etot - e_prev;
        e_prev = etot;

        // Single-block strategies converge as soon as their block does.
        if (params.strategy == SolverStrategy::OccOnly)
        {
            res.converged = occ_conv;
            if (occ_conv)
            {
                break;
            }
            continue;
        }
        if (params.strategy == SolverStrategy::OrbOnly)
        {
            res.converged = orb_conv;
            if (orb_conv)
            {
                break;
            }
            continue;
        }

        if (outer > 0 && std::fabs(dE) < params.energy_tol && occ_conv && orb_conv)
        {
            res.converged = true;
            break;
        }
        if (occ_conv && orb_conv && !do_orb)
        {
            res.converged = true;
            break;
        }
    }

    res.energy = etot;
    return res;
}

} // namespace rdmft
