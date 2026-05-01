// =============================================================================
// Unit tests for the Spectral Projected Gradient (SPG) occupation optimisation
// path in RDMFT (used by both `projected_gradient` and `active_set` constraint
// methods, which now route to the same SPG core).
//
// The RDMFT solver's full `optimize_occupations()` cannot be exercised in a
// unit test because it depends on EnergyGradient<> (which requires LCAO
// infrastructure). These tests therefore reproduce the same SPG algorithmic
// flow (Birgin-Martinez-Raydan, SIOPT 2000; Grippo-Lampariello-Lucidi
// non-monotone Armijo, SINUM 1986) on a toy problem:
//
//     E(n) = sum_i n_i * eps_i,
//     subject to  sum_i w_i n_i = N_e,  0 <= n_i <= 1.
//
// The exact minimum is obtained by greedily filling from lowest energy. We
// verify that the SPG runner
//   (a) converges to this minimum,
//   (b) satisfies the electron-number constraint at convergence,
//   (c) keeps all n_i in [0, 1].
// =============================================================================

#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_occupation.h"
#include "source_lcao/module_rdmft/rdmft_optimizer.h"
#include "source_lcao/module_rdmft/rdmft_type.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <vector>

using namespace rdmft;

namespace
{

// ---------------------------------------------------------------------------
// Toy occupation-number optimization problem.
//
//   E(n) = sum_i eps_i * n_i
//   subject to:  sum_i w_i * n_i = N_e,   0 <= n_i <= 1.
//
// Exact solution: fill from lowest energy until N_e electrons are placed.
// ---------------------------------------------------------------------------
struct ToyOccProblem
{
    int nb = 0;               // number of bands
    double Ne = 0.0;          // target electron number
    std::vector<double> eps;  // band energies, size nb
    std::vector<double> wk;   // k-point weights, size 1 (single k-point)

    double energy(const std::vector<double>& occ) const
    {
        double E = 0.0;
        for (int i = 0; i < nb; ++i)
            E += eps[i] * occ[i];
        return E;
    }

    std::vector<double> gradient(const std::vector<double>& occ) const
    {
        (void)occ;
        // dE/dn_i = eps_i
        return eps;
    }

    // Greedy analytic minimum: fill from lowest eps.
    double ref_energy() const
    {
        std::vector<int> order(nb);
        for (int i = 0; i < nb; ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return eps[a] < eps[b]; });

        double E = 0.0;
        double rem = Ne;
        for (int i : order)
        {
            double ni = std::min(1.0, rem);
            E += ni * eps[i];
            rem -= ni;
            if (rem <= 0.0) break;
        }
        return E;
    }

    std::vector<double> ref_occ() const
    {
        std::vector<int> order(nb);
        for (int i = 0; i < nb; ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return eps[a] < eps[b]; });

        std::vector<double> n(nb, 0.0);
        double rem = Ne;
        for (int i : order)
        {
            n[i] = std::min(1.0, rem);
            rem -= n[i];
            if (rem <= 0.0) break;
        }
        return n;
    }
};

// Make a simple 4-band problem with N_e = 2, wk = {1.0}.
ToyOccProblem make_4band(double Ne = 2.0)
{
    ToyOccProblem p;
    p.nb = 4;
    p.Ne = Ne;
    p.eps = {1.0, 2.0, 3.0, 4.0};
    p.wk = {1.0};
    return p;
}

// Make a 6-band problem that has fractional occupations at the minimum.
ToyOccProblem make_6band_frac(double Ne = 3.5)
{
    ToyOccProblem p;
    p.nb = 6;
    p.Ne = Ne;
    p.eps = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
    p.wk = {1.0};
    return p;
}

// ---------------------------------------------------------------------------
// Run the Spectral Projected Gradient optimisation on the toy problem.
// Mirrors rdmft_solver.cpp's optimize_occupations SPG block (the case that
// handles both ConstraintMethod::ProjectedGradient and ConstraintMethod::ActiveSet).
// ---------------------------------------------------------------------------
struct SPGRunner
{
    ToyOccProblem prob;
    RDMFTConfig   config;
    ConstraintMethod method = ConstraintMethod::ProjectedGradient;
    mutable std::vector<double> energy_history;

    // Returns final energy after at most max_iter SPG iterations.
    double run(std::vector<double>& occ, int max_iter = 200) const
    {
        OccupationConstraint constraint(method, prob.Ne, prob.wk, prob.nb);

        // SPG safeguarding constants (BMR Algorithm 2.1).
        const double alpha_min = 1e-10;
        const double alpha_max = 1e10;
        const int nm_memory = 10;
        const double rho = (config.line_search_rho > 0.0 && config.line_search_rho < 1.0)
                               ? config.line_search_rho : 0.5;
        const double c1 = config.line_search_c1;

        std::deque<double> f_history;
        std::vector<double> occ_prev;
        std::vector<double> grad_prev;
        bool have_prev = false;

        double E = 0.0;

        for (int iter = 0; iter < max_iter; ++iter)
        {
            E = prob.energy(occ);
            energy_history.push_back(E);
            auto grad = prob.gradient(occ);

            // Spectral (BB1) step.
            double alpha_bb;
            if (have_prev)
            {
                double ss = 0.0, sy = 0.0;
                for (int i = 0; i < prob.nb; ++i)
                {
                    const double s = occ[i] - occ_prev[i];
                    const double y = grad[i] - grad_prev[i];
                    ss += s * s;
                    sy += s * y;
                }
                alpha_bb = (sy > 1e-30 && std::isfinite(ss)) ? ss / sy : 1.0;
                if (!std::isfinite(alpha_bb) || alpha_bb <= 0.0) alpha_bb = 1.0;
            }
            else
            {
                double g_inf = 0.0;
                for (double g : grad) g_inf = std::max(g_inf, std::abs(g));
                alpha_bb = 1.0 / std::max(1.0, g_inf);
            }
            alpha_bb = std::min(alpha_max, std::max(alpha_min, alpha_bb));

            // Spectral projected direction d = P(n - α_BB g) - n.
            std::vector<double> trial_proj(prob.nb);
            for (int i = 0; i < prob.nb; ++i) trial_proj[i] = occ[i] - alpha_bb * grad[i];
            constraint.project(trial_proj);

            std::vector<double> dir(prob.nb);
            double dd = 0.0;
            for (int i = 0; i < prob.nb; ++i)
            {
                dir[i] = trial_proj[i] - occ[i];
                dd += grad[i] * dir[i];
            }

            // Non-monotone Armijo reference.
            if (static_cast<int>(f_history.size()) >= nm_memory) f_history.pop_front();
            f_history.push_back(E);
            double f_max = E;
            for (double f : f_history) if (f > f_max) f_max = f;

            // Non-monotone Armijo backtracking along the convex segment n + λ d.
            // Each n + λ d is a convex combination of two feasible points,
            // hence feasible without re-projection.
            const std::vector<double> occ_old = occ;
            double lambda = 1.0;
            std::vector<double> occ_trial(prob.nb);
            bool success = false;
            for (int ls = 0; ls < config.line_search_max_iter; ++ls)
            {
                for (int i = 0; i < prob.nb; ++i) occ_trial[i] = occ[i] + lambda * dir[i];
                const double E_trial = prob.energy(occ_trial);
                if (E_trial <= f_max + c1 * lambda * dd)
                {
                    success = true;
                    break;
                }
                lambda *= rho;
            }

            occ_prev = occ;
            grad_prev = grad;
            have_prev = true;

            if (success)
            {
                occ = occ_trial;
            }
            // else: zero step; SPG line search in floating point can fail only
            // for degenerate dd ≈ 0 -- treat as a stationary iterate.

            // Convergence on L1 occupation move.
            double sum_dn = 0.0;
            for (int i = 0; i < prob.nb; ++i) sum_dn += std::abs(occ[i] - occ_old[i]);
            if (sum_dn < config.rdmft_occ_tol && sum_dn > 1e-20) break;
        }

        return E;
    }
};

// Initial feasible point: equal distribution clamped to interior.
std::vector<double> initial_occ_uniform(int nb, double Ne)
{
    std::vector<double> occ(nb, Ne / nb);
    // Clamp to interior to avoid gradient vanishing at boundary.
    const double m = 1e-3;
    for (auto& n : occ) n = std::max(m, std::min(1.0 - m, n));
    // Rescale to exactly satisfy constraint.
    double s = 0.0;
    for (auto n : occ) s += n;
    if (s > 1e-12)
        for (auto& n : occ) n *= Ne / s;
    return occ;
}

} // anonymous namespace

// =============================================================================
// Tests
// =============================================================================

class SPGOccupationTest : public ::testing::Test {};

// ---------------------------------------------------------------------------
// SPG converges on the 4-band integer-fill problem.
// ---------------------------------------------------------------------------
TEST_F(SPGOccupationTest, converges_4band)
{
    auto prob = make_4band(2.0);
    RDMFTConfig cfg;
    cfg.rdmft_occ_tol = 1e-10;

    SPGRunner runner{prob, cfg, ConstraintMethod::ProjectedGradient, {}};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    double E = runner.run(occ, 200);

    EXPECT_NEAR(E, prob.ref_energy(), 1e-5);

    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, prob.Ne, 1e-8);

    for (auto n : occ)
    {
        EXPECT_GE(n, -1e-12);
        EXPECT_LE(n, 1.0 + 1e-12);
    }
}

// ---------------------------------------------------------------------------
// SPG: fractional electron count gives fractional minimiser.
// ---------------------------------------------------------------------------
TEST_F(SPGOccupationTest, converges_fractional)
{
    auto prob = make_6band_frac(3.5);
    RDMFTConfig cfg;
    cfg.rdmft_occ_tol = 1e-10;

    SPGRunner runner{prob, cfg, ConstraintMethod::ProjectedGradient, {}};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    double E = runner.run(occ, 300);

    EXPECT_NEAR(E, prob.ref_energy(), 1e-4);

    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, prob.Ne, 1e-8);
}

// ---------------------------------------------------------------------------
// SPG: constraint satisfied at convergence (non-integer Ne).
// ---------------------------------------------------------------------------
TEST_F(SPGOccupationTest, constraint_satisfied)
{
    auto prob = make_4band(1.5);
    RDMFTConfig cfg;
    cfg.rdmft_occ_tol = 1e-6;

    SPGRunner runner{prob, cfg, ConstraintMethod::ProjectedGradient, {}};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    runner.run(occ, 200);

    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, prob.Ne, 1e-8);
    for (auto n : occ)
    {
        EXPECT_GE(n, -1e-10);
        EXPECT_LE(n, 1.0 + 1e-10);
    }
}

// ---------------------------------------------------------------------------
// SPG with `ActiveSet` constraint method routes to the same SPG core and
// reaches the same energy as `ProjectedGradient` (both are aliases).
// ---------------------------------------------------------------------------
TEST_F(SPGOccupationTest, PG_and_AS_constraint_aliases_agree)
{
    auto prob = make_4band(2.5);

    RDMFTConfig cfg;
    cfg.rdmft_occ_tol = 1e-10;

    SPGRunner pg{prob, cfg, ConstraintMethod::ProjectedGradient, {}};
    SPGRunner as{prob, cfg, ConstraintMethod::ActiveSet, {}};

    auto occ_pg = initial_occ_uniform(prob.nb, prob.Ne);
    auto occ_as = occ_pg;

    double E_pg = pg.run(occ_pg, 300);
    double E_as = as.run(occ_as, 300);

    EXPECT_NEAR(E_pg, prob.ref_energy(), 1e-5);
    EXPECT_NEAR(E_as, prob.ref_energy(), 1e-5);
    EXPECT_NEAR(E_pg, E_as, 1e-5);
}

// ---------------------------------------------------------------------------
// Stress test: large initial gradient (toy analogue of n^alpha at n -> 0)
// must not stall.  This is the regression scenario from the NiO logs that
// motivated dropping the trust-region cap and the SD-fallback chain.
// ---------------------------------------------------------------------------
TEST_F(SPGOccupationTest, large_gradient_no_stall)
{
    ToyOccProblem p;
    p.nb = 6;
    p.Ne = 2.5;
    p.eps = {-1.0e4, -5.0e3, 1.0e3, 5.0e3, 1.0e4, 2.0e4};
    p.wk = {1.0};

    RDMFTConfig cfg;
    cfg.rdmft_occ_tol = 1e-8;
    cfg.line_search_max_iter = 20;

    SPGRunner runner{p, cfg, ConstraintMethod::ProjectedGradient, {}};
    auto occ = initial_occ_uniform(p.nb, p.Ne);
    double E = runner.run(occ, 200);

    EXPECT_NEAR(E, p.ref_energy(), std::abs(p.ref_energy()) * 1e-6 + 1e-4);
    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, p.Ne, 1e-8);
}

// ---------------------------------------------------------------------------
// RDMFTConfig: constraint method and occ-optimizer types can be round-tripped.
// ---------------------------------------------------------------------------
TEST(RdmftConfigTest, constraint_method_roundtrip)
{
    RDMFTConfig cfg;
    cfg.constraint_method = ConstraintMethod::ProjectedGradient;
    EXPECT_EQ(cfg.constraint_method, ConstraintMethod::ProjectedGradient);

    cfg.constraint_method = ConstraintMethod::ActiveSet;
    EXPECT_EQ(cfg.constraint_method, ConstraintMethod::ActiveSet);

    cfg.constraint_method = ConstraintMethod::AugmentedLagrangian;
    EXPECT_EQ(cfg.constraint_method, ConstraintMethod::AugmentedLagrangian);
}

TEST(RdmftConfigTest, occ_optimizer_type_roundtrip)
{
    RDMFTConfig cfg;
    for (OptimizerType t : {OptimizerType::SteepestDescent,
                             OptimizerType::ConjugateGradient,
                             OptimizerType::LBFGS,
                             OptimizerType::Adam})
    {
        cfg.occ_optimizer = t;
        EXPECT_EQ(cfg.occ_optimizer, t);
    }
}

// ---------------------------------------------------------------------------
// OccupationConstraint::project satisfies the electron-number constraint.
// ---------------------------------------------------------------------------
TEST(OccupationConstraintProjectTest, project_satisfies_nel)
{
    std::vector<double> wk = {1.0};
    int nb = 5;
    double Ne = 3.0;
    OccupationConstraint c(ConstraintMethod::ProjectedGradient, Ne, wk, nb);

    std::vector<double> occ = {1.2, 0.8, -0.1, 0.5, 0.9};
    c.project(occ);

    // Box constraint.
    for (auto n : occ)
    {
        EXPECT_GE(n, -1e-14);
        EXPECT_LE(n, 1.0 + 1e-14);
    }

    // Equality constraint.
    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, Ne, 1e-12);
}

// ---------------------------------------------------------------------------
// OccupationConstraint::apply_active_set: gradient in null space of constraint
//   sum_k w_k * dg_k = 0 for free variables.
// ---------------------------------------------------------------------------
TEST(OccupationConstraintASTest, reduced_gradient_in_constraint_null_space)
{
    std::vector<double> wk = {1.0};
    int nb = 4;
    double Ne = 2.0;
    OccupationConstraint c(ConstraintMethod::ActiveSet, Ne, wk, nb);

    // All bands are free (no active constraints).
    std::vector<double> occ = {0.3, 0.5, 0.7, 0.5};
    std::vector<double> grad = {-0.5, -0.3, 0.2, 0.8};

    auto info = c.identify_active_set(occ, grad);
    std::vector<double> grad_mod(grad);
    c.apply_active_set(info, grad_mod);

    // For free variables the reduced gradient should satisfy
    //   sum_i w_k * grad_mod_i = 0  (tangency to equality constraint).
    double sum_g = 0.0;
    for (int i = 0; i < nb; ++i)
        if (info.is_free[i]) sum_g += wk[0] * grad_mod[i];
    EXPECT_NEAR(sum_g, 0.0, 1e-12);
}
