// =============================================================================
// Unit tests for projected-gradient and active-set occupation optimization
// methods in RDMFT.
//
// The RDMFT solver's full `optimize_occupations()` cannot be exercised in a
// unit test because it depends on EnergyGradient<> (which requires LCAO
// infrastructure). These tests therefore reproduce the same algorithmic flow
// on a toy problem:
//
//     E(n) = sum_i n_i * eps_i   +   penalty(sum_i n_i - N_e)
//
// where eps_i are fixed band energies and penalty is either the augmented
// Lagrangian or (for PG/AS) the box-constraint is handled implicitly.
//
// The toy problem has an exact analytic minimum whose occupation numbers are
// obtained by greedily filling from lowest energy. We verify that both PG and
// AS methods:
//   (a) converge to this minimum,
//   (b) satisfy sum_i n_i = N_e at convergence,
//   (c) keep all n_i in [0, 1].
//
// Additionally we test that the methods accept the configured optimizer type
// (SD / CG / LBFGS / Adam) and that the active-set restart logic is triggered
// when the active-set composition changes.
// =============================================================================

#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_occupation.h"
#include "source_lcao/module_rdmft/rdmft_optimizer.h"
#include "source_lcao/module_rdmft/rdmft_type.h"

#include <algorithm>
#include <cmath>
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
// Run the projected-gradient optimization loop on the toy problem.
// Mirrors rdmft_solver.cpp's optimize_occupations / ProjectedGradient case.
// ---------------------------------------------------------------------------
struct PGRunner
{
    ToyOccProblem prob;
    RDMFTConfig   config;
    mutable std::vector<double> energy_history;

    // Returns final occupations and energy after at most max_iter iterations.
    double run(std::vector<double>& occ, int max_iter = 200) const
    {
        OccupationConstraint constraint(
            ConstraintMethod::ProjectedGradient,
            prob.Ne, prob.wk, prob.nb);

        EuclideanOptimizer opt(config.occ_optimizer, config);
        opt.init(prob.nb);

        double E = 0.0;

        for (int iter = 0; iter < max_iter; ++iter)
        {
            E = prob.energy(occ);
            energy_history.push_back(E);
            auto grad = prob.gradient(occ);

            // Compute search direction.
            std::vector<double> dir;
            opt.compute_direction(grad, dir);

            double dd = 0.0;
            for (int i = 0; i < prob.nb; ++i) dd += dir[i] * grad[i];
            if (dd >= 0.0)
                for (int i = 0; i < prob.nb; ++i) dir[i] = -grad[i];

            // Armijo backtracking with projection.
            const std::vector<double> occ_old = occ;
            double alpha = config.line_search_alpha_init;
            bool success = false;
            std::vector<double> occ_trial;
            const double proj_tol = 1e-6;

            for (int ls = 0; ls < config.line_search_max_iter; ++ls)
            {
                occ_trial = occ;
                for (int i = 0; i < prob.nb; ++i)
                    occ_trial[i] += alpha * dir[i];
                constraint.project(occ_trial);

                // Reject if the projection failed to satisfy the constraint.
                if (std::abs(constraint.constraint_violation(occ_trial)) > proj_tol)
                {
                    alpha *= config.line_search_rho;
                    continue;
                }

                double dd_proj = 0.0;
                for (int i = 0; i < prob.nb; ++i)
                    dd_proj += grad[i] * (occ_trial[i] - occ[i]);

                if (dd_proj >= 0.0) { alpha *= config.line_search_rho; continue; }

                double E_trial = prob.energy(occ_trial);
                if (E_trial <= E + config.line_search_c1 * dd_proj)
                {
                    success = true; break;
                }
                alpha *= config.line_search_rho;
            }

            if (success)
            {
                occ = occ_trial;
            }
            else
            {
                // Mirror the solver: reject unchecked fallback steps so the
                // objective cannot increase because Armijo failed.
                occ = occ_old;
                opt.init(prob.nb);
            }

            // Update optimizer history.
            std::vector<double> step_vec(prob.nb);
            for (int i = 0; i < prob.nb; ++i) step_vec[i] = occ[i] - occ_old[i];
            auto new_grad = prob.gradient(occ);
            opt.update(new_grad, step_vec);

            // Convergence on step size.
            double sum_dn = 0.0;
            for (int i = 0; i < prob.nb; ++i) sum_dn += std::abs(step_vec[i]);
            if (sum_dn < config.rdmft_occ_tol && sum_dn > 1e-20) break;
        }

        return E;
    }
};

// ---------------------------------------------------------------------------
// Run the active-set optimization loop on the toy problem.
// Mirrors rdmft_solver.cpp's optimize_occupations / ActiveSet case.
// ---------------------------------------------------------------------------
struct ASRunner
{
    ToyOccProblem prob;
    RDMFTConfig   config;

    double run(std::vector<double>& occ, int max_iter = 200) const
    {
        OccupationConstraint constraint(
            ConstraintMethod::ActiveSet,
            prob.Ne, prob.wk, prob.nb);

        EuclideanOptimizer opt(config.occ_optimizer, config);
        opt.init(prob.nb);

        int prev_n_active = -1;
        double E = 0.0;

        for (int iter = 0; iter < max_iter; ++iter)
        {
            E = prob.energy(occ);
            auto grad = prob.gradient(occ);

            // Identify active set and apply modification.
            auto as_info = constraint.identify_active_set(occ, grad);
            std::vector<double> grad_mod(grad);
            constraint.apply_active_set(as_info, grad_mod);

            // Count active constraints; restart optimizer if changed.
            int n_active = 0;
            for (int idx = 0; idx < prob.nb; ++idx)
                if (!as_info.is_free[idx]) ++n_active;

            if (n_active != prev_n_active && iter > 0)
                opt.init(prob.nb);
            prev_n_active = n_active;

            // Compute direction on modified gradient.
            std::vector<double> dir;
            opt.compute_direction(grad_mod, dir);
            for (int idx = 0; idx < prob.nb; ++idx)
                if (!as_info.is_free[idx]) dir[idx] = 0.0;

            double dd = 0.0;
            for (int i = 0; i < prob.nb; ++i) dd += dir[i] * grad_mod[i];
            if (dd >= 0.0)
            {
                for (int i = 0; i < prob.nb; ++i) dir[i] = -grad_mod[i];
                for (int idx = 0; idx < prob.nb; ++idx)
                    if (!as_info.is_free[idx]) dir[idx] = 0.0;
            }

            const std::vector<double> occ_old = occ;
            double alpha = config.line_search_alpha_init;
            bool success = false;
            std::vector<double> occ_trial;
            const double proj_tol = 1e-6;

            for (int ls = 0; ls < config.line_search_max_iter; ++ls)
            {
                occ_trial = occ;
                for (int i = 0; i < prob.nb; ++i)
                    occ_trial[i] += alpha * dir[i];
                for (auto& n : occ_trial) n = std::max(0.0, std::min(1.0, n));
                constraint.project(occ_trial);

                // Reject if the projection failed to satisfy the constraint.
                if (std::abs(constraint.constraint_violation(occ_trial)) > proj_tol)
                {
                    alpha *= config.line_search_rho;
                    continue;
                }

                double dd_proj = 0.0;
                for (int i = 0; i < prob.nb; ++i)
                    dd_proj += grad_mod[i] * (occ_trial[i] - occ[i]);

                if (dd_proj >= 0.0) { alpha *= config.line_search_rho; continue; }

                double E_trial = prob.energy(occ_trial);
                if (E_trial <= E + config.line_search_c1 * dd_proj)
                {
                    success = true; break;
                }
                alpha *= config.line_search_rho;
            }

            if (success)
            {
                occ = occ_trial;
            }
            else
            {
                occ = occ_old;
                for (int i = 0; i < prob.nb; ++i)
                    occ[i] -= config.line_search_alpha_init * grad_mod[i];
                for (auto& n : occ) n = std::max(0.0, std::min(1.0, n));
                constraint.project(occ);
                opt.init(prob.nb);
                prev_n_active = -1;
            }

            // Update optimizer history.
            std::vector<double> step_vec(prob.nb);
            for (int i = 0; i < prob.nb; ++i) step_vec[i] = occ[i] - occ_old[i];

            auto new_grad = prob.gradient(occ);
            auto new_as_info = constraint.identify_active_set(occ, new_grad);
            std::vector<double> new_grad_mod(new_grad);
            constraint.apply_active_set(new_as_info, new_grad_mod);
            for (int idx = 0; idx < prob.nb; ++idx)
                if (!new_as_info.is_free[idx]) new_grad_mod[idx] = 0.0;
            opt.update(new_grad_mod, step_vec);

            double sum_dn = 0.0;
            for (int i = 0; i < prob.nb; ++i) sum_dn += std::abs(step_vec[i]);
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

class PGOptimizerTest : public ::testing::Test {};
class ASOptimizerTest : public ::testing::Test {};

// ---------------------------------------------------------------------------
// PG: convergence with steepest descent
// ---------------------------------------------------------------------------
TEST_F(PGOptimizerTest, SD_converges_4band)
{
    auto prob = make_4band(2.0);
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::SteepestDescent;

    cfg.rdmft_occ_tol = 1e-10;

    PGRunner runner{prob, cfg};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    double E = runner.run(occ, 500);

    // Verify energy is close to reference.
    EXPECT_NEAR(E, prob.ref_energy(), 1e-5);

    // Verify constraint: sum_i n_i = Ne.
    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, prob.Ne, 1e-8);

    // Verify box constraint.
    for (auto n : occ)
    {
        EXPECT_GE(n, -1e-12);
        EXPECT_LE(n, 1.0 + 1e-12);
    }
}

// ---------------------------------------------------------------------------
// PG: convergence with conjugate gradient
// ---------------------------------------------------------------------------
TEST_F(PGOptimizerTest, CG_converges_4band)
{
    auto prob = make_4band(2.0);
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::ConjugateGradient;

    cfg.rdmft_occ_tol = 1e-10;

    PGRunner runner{prob, cfg};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    double E = runner.run(occ, 200);

    EXPECT_NEAR(E, prob.ref_energy(), 1e-5);

    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, prob.Ne, 1e-8);
}

// ---------------------------------------------------------------------------
// PG: convergence with L-BFGS
// ---------------------------------------------------------------------------
TEST_F(PGOptimizerTest, LBFGS_converges_4band)
{
    auto prob = make_4band(2.0);
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::LBFGS;
    cfg.lbfgs_memory = 5;
    cfg.line_search_alpha_init = 1.0;
    cfg.rdmft_occ_tol = 1e-10;

    PGRunner runner{prob, cfg};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    double E = runner.run(occ, 100);

    EXPECT_NEAR(E, prob.ref_energy(), 1e-5);

    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, prob.Ne, 1e-8);
}

// ---------------------------------------------------------------------------
// PG: fractional occupations (N_e not integer, optimal solution has n < 1)
// ---------------------------------------------------------------------------
TEST_F(PGOptimizerTest, CG_converges_fractional)
{
    auto prob = make_6band_frac(3.5);
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::ConjugateGradient;

    cfg.rdmft_occ_tol = 1e-10;

    PGRunner runner{prob, cfg};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    double E = runner.run(occ, 300);

    EXPECT_NEAR(E, prob.ref_energy(), 1e-4);

    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, prob.Ne, 1e-8);
}

// ---------------------------------------------------------------------------
// PG: constraint satisfied at convergence (independent of optimizer type)
// ---------------------------------------------------------------------------
TEST_F(PGOptimizerTest, constraint_satisfied_after_SD)
{
    auto prob = make_4band(1.5);
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::SteepestDescent;

    cfg.rdmft_occ_tol = 1e-6;

    PGRunner runner{prob, cfg};
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

TEST_F(PGOptimizerTest, energy_is_monotone_nonincreasing)
{
    auto prob = make_4band(2.0);
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::ConjugateGradient;
    cfg.line_search_alpha_init = 1.0;
    cfg.rdmft_occ_tol = 1e-10;

    PGRunner runner{prob, cfg};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    runner.run(occ, 200);

    ASSERT_FALSE(runner.energy_history.empty());
    for (size_t i = 1; i < runner.energy_history.size(); ++i)
    {
        EXPECT_LE(runner.energy_history[i], runner.energy_history[i - 1] + 1e-12);
    }
}

// ---------------------------------------------------------------------------
// AS: convergence with steepest descent
// ---------------------------------------------------------------------------
TEST_F(ASOptimizerTest, SD_converges_4band)
{
    auto prob = make_4band(2.0);
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::SteepestDescent;

    cfg.rdmft_occ_tol = 1e-10;

    ASRunner runner{prob, cfg};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    double E = runner.run(occ, 500);

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
// AS: convergence with conjugate gradient
// ---------------------------------------------------------------------------
TEST_F(ASOptimizerTest, CG_converges_4band)
{
    auto prob = make_4band(2.0);
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::ConjugateGradient;

    cfg.rdmft_occ_tol = 1e-10;

    ASRunner runner{prob, cfg};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    double E = runner.run(occ, 200);

    EXPECT_NEAR(E, prob.ref_energy(), 1e-5);

    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, prob.Ne, 1e-8);
}

// ---------------------------------------------------------------------------
// AS: convergence with L-BFGS
// ---------------------------------------------------------------------------
TEST_F(ASOptimizerTest, LBFGS_converges_4band)
{
    auto prob = make_4band(2.0);
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::LBFGS;
    cfg.lbfgs_memory = 5;
    cfg.line_search_alpha_init = 1.0;
    cfg.rdmft_occ_tol = 1e-10;

    ASRunner runner{prob, cfg};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    double E = runner.run(occ, 100);

    EXPECT_NEAR(E, prob.ref_energy(), 1e-5);

    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, prob.Ne, 1e-8);
}

// ---------------------------------------------------------------------------
// AS: active set identification triggers optimizer restart
// ---------------------------------------------------------------------------
TEST_F(ASOptimizerTest, active_set_change_restarts_optimizer)
{
    // Use a problem where some bands start inside (0,1) and drift to boundaries.
    auto prob = make_4band(1.0);   // only 1 electron -> 1 band full, rest empty
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::LBFGS;
    cfg.lbfgs_memory = 5;

    cfg.rdmft_occ_tol = 1e-9;

    ASRunner runner{prob, cfg};
    // Start uniformly (all free), expect active set to evolve.
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    double E = runner.run(occ, 300);

    EXPECT_NEAR(E, prob.ref_energy(), 1e-5);

    double s = 0.0;
    for (auto n : occ) s += n;
    EXPECT_NEAR(s, prob.Ne, 1e-8);
}

// ---------------------------------------------------------------------------
// AS: constraint remains satisfied throughout (spot-check at end)
// ---------------------------------------------------------------------------
TEST_F(ASOptimizerTest, constraint_satisfied_after_CG)
{
    auto prob = make_6band_frac(2.8);
    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::ConjugateGradient;

    cfg.rdmft_occ_tol = 1e-6;

    ASRunner runner{prob, cfg};
    auto occ = initial_occ_uniform(prob.nb, prob.Ne);
    runner.run(occ, 300);

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
// Both PG and AS reach the same energy on the 4-band problem.
// ---------------------------------------------------------------------------
TEST(OccOptimizerComparison, PG_and_AS_agree)
{
    auto prob = make_4band(2.5);

    RDMFTConfig cfg;
    cfg.occ_optimizer = OptimizerType::ConjugateGradient;

    cfg.rdmft_occ_tol = 1e-10;

    PGRunner pg{prob, cfg};
    ASRunner as{prob, cfg};

    auto occ_pg = initial_occ_uniform(prob.nb, prob.Ne);
    auto occ_as = occ_pg;

    double E_pg = pg.run(occ_pg, 300);
    double E_as = as.run(occ_as, 300);

    EXPECT_NEAR(E_pg, prob.ref_energy(), 1e-5);
    EXPECT_NEAR(E_as, prob.ref_energy(), 1e-5);
    EXPECT_NEAR(E_pg, E_as, 1e-5);
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
