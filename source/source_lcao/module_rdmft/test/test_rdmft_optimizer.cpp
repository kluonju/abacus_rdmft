#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_optimizer.h"
#include "source_lcao/module_rdmft/rdmft_stiefel.h"
#include <cmath>
#include <vector>
#include <functional>
#include <random>

using namespace rdmft;

class OptimizerTest : public ::testing::Test {};

// Simple quadratic: f(x) = 0.5 * sum_i a_i * (x_i - b_i)^2
struct Quadratic
{
    std::vector<double> a;
    std::vector<double> b;

    double eval(const std::vector<double>& x) const
    {
        double f = 0.0;
        for (size_t i = 0; i < x.size(); ++i)
            f += 0.5 * a[i] * (x[i] - b[i]) * (x[i] - b[i]);
        return f;
    }

    std::vector<double> grad(const std::vector<double>& x) const
    {
        std::vector<double> g(x.size());
        for (size_t i = 0; i < x.size(); ++i)
            g[i] = a[i] * (x[i] - b[i]);
        return g;
    }
};

TEST_F(OptimizerTest, SteepestDescent_converges_quadratic)
{
    Quadratic q;
    q.a = {1.0, 2.0, 3.0};
    q.b = {1.0, -1.0, 0.5};

    RDMFTConfig config;
    EuclideanOptimizer opt(OptimizerType::SteepestDescent, config);
    opt.init(3);

    std::vector<double> x = {0.0, 0.0, 0.0};
    for (int iter = 0; iter < 200; ++iter)
    {
        auto g = q.grad(x);
        std::vector<double> dir;
        opt.compute_direction(g, dir);

        double step = 0.5;
        for (size_t i = 0; i < x.size(); ++i)
            x[i] += step * dir[i];
    }

    for (size_t i = 0; i < x.size(); ++i)
        EXPECT_NEAR(x[i], q.b[i], 1e-4);
}

TEST_F(OptimizerTest, CG_converges_quadratic)
{
    Quadratic q;
    q.a = {1.0, 2.0, 3.0};
    q.b = {1.0, -1.0, 0.5};

    RDMFTConfig config;
    EuclideanOptimizer opt(OptimizerType::ConjugateGradient, config);
    opt.init(3);

    std::vector<double> x = {0.0, 0.0, 0.0};
    for (int iter = 0; iter < 100; ++iter)
    {
        auto g = q.grad(x);
        std::vector<double> dir;
        opt.compute_direction(g, dir);

        double dd = 0.0;
        for (size_t i = 0; i < g.size(); ++i) dd += dir[i] * g[i];

        auto f_at = [&](double step) -> double {
            std::vector<double> xt(x);
            for (size_t i = 0; i < xt.size(); ++i) xt[i] += step * dir[i];
            return q.eval(xt);
        };

        auto ls = armijo_line_search(f_at, q.eval(x), dd, 1.0);
        std::vector<double> step_vec(x.size());
        for (size_t i = 0; i < x.size(); ++i)
        {
            step_vec[i] = ls.step * dir[i];
            x[i] += step_vec[i];
        }

        auto new_g = q.grad(x);
        opt.update(new_g, step_vec);
    }

    for (size_t i = 0; i < x.size(); ++i)
        EXPECT_NEAR(x[i], q.b[i], 1e-6);
}

TEST_F(OptimizerTest, LBFGS_converges_quadratic)
{
    Quadratic q;
    q.a = {1.0, 2.0, 3.0, 4.0, 5.0};
    q.b = {1.0, -1.0, 0.5, 2.0, -0.5};

    RDMFTConfig config;
    config.lbfgs_memory = 5;
    EuclideanOptimizer opt(OptimizerType::LBFGS, config);
    opt.init(5);

    std::vector<double> x(5, 0.0);
    for (int iter = 0; iter < 50; ++iter)
    {
        auto g = q.grad(x);
        std::vector<double> dir;
        opt.compute_direction(g, dir);

        double dd = 0.0;
        for (size_t i = 0; i < g.size(); ++i) dd += dir[i] * g[i];

        auto f_at = [&](double step) -> double {
            std::vector<double> xt(x);
            for (size_t i = 0; i < xt.size(); ++i) xt[i] += step * dir[i];
            return q.eval(xt);
        };

        auto ls = armijo_line_search(f_at, q.eval(x), dd, 1.0);
        std::vector<double> step_vec(x.size());
        for (size_t i = 0; i < x.size(); ++i)
        {
            step_vec[i] = ls.step * dir[i];
            x[i] += step_vec[i];
        }

        auto new_g = q.grad(x);
        opt.update(new_g, step_vec);
    }

    for (size_t i = 0; i < x.size(); ++i)
        EXPECT_NEAR(x[i], q.b[i], 1e-6);
}

TEST_F(OptimizerTest, Adam_converges_quadratic)
{
    Quadratic q;
    q.a = {1.0, 2.0, 3.0};
    q.b = {1.0, -1.0, 0.5};

    RDMFTConfig config;
    config.adam_lr = 0.1;
    EuclideanOptimizer opt(OptimizerType::Adam, config);
    opt.init(3);

    std::vector<double> x = {0.0, 0.0, 0.0};
    for (int iter = 0; iter < 500; ++iter)
    {
        auto g = q.grad(x);
        std::vector<double> dir;
        opt.compute_direction(g, dir);

        // Adam already includes learning rate in the direction
        for (size_t i = 0; i < x.size(); ++i)
            x[i] += dir[i];
    }

    for (size_t i = 0; i < x.size(); ++i)
        EXPECT_NEAR(x[i], q.b[i], 0.05);
}

TEST_F(OptimizerTest, armijo_line_search_works)
{
    auto f = [](double step) -> double { return (step - 1.0) * (step - 1.0); };
    double f0 = f(0.0);
    double deriv = -2.0; // f'(0) = 2*(0-1) = -2

    auto result = armijo_line_search(f, f0, deriv, 2.0);
    EXPECT_TRUE(result.success);
    EXPECT_GT(result.step, 0.0);
    EXPECT_LT(result.f_new, f0);
}

// ============================================================================
// Stiefel-manifold orbital-optimisation tests.
//
// These tests mirror the solver's optimize_orbitals() flow: flatten the
// Riemannian gradient, feed it to the EuclideanOptimizer, unflatten the
// direction, project to the tangent space, and retract.  The test problem
// is the Rayleigh-quotient minimisation
//     min_{C^T C = I}  Tr(C^T A C)
// whose minimum is attained when C spans the invariant subspace associated
// with the p smallest eigenvalues of A.  We verify that L-BFGS and Adam
// both converge to this minimum, matching the SD / CG baseline.
// ============================================================================

namespace
{
std::vector<double> random_stiefel_point(int n, int p, std::mt19937& rng)
{
    std::vector<double> C(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : C) v = dist(rng);

    // Gram-Schmidt orthonormalisation (n > p assumed)
    for (int j = 0; j < p; ++j)
    {
        for (int k = 0; k < j; ++k)
        {
            double dot = 0.0;
            for (int i = 0; i < n; ++i)
                dot += C[i + k * n] * C[i + j * n];
            for (int i = 0; i < n; ++i)
                C[i + j * n] -= dot * C[i + k * n];
        }
        double norm = 0.0;
        for (int i = 0; i < n; ++i)
            norm += C[i + j * n] * C[i + j * n];
        norm = std::sqrt(norm);
        for (int i = 0; i < n; ++i)
            C[i + j * n] /= norm;
    }
    return C;
}

struct Rayleigh
{
    int n, p;
    std::vector<double> A; // symmetric n x n, column-major

    double eval(const std::vector<double>& C) const
    {
        double f = 0.0;
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < n; ++i)
            {
                double ac = 0.0;
                for (int k = 0; k < n; ++k)
                    ac += A[i + k * n] * C[k + j * n];
                f += C[i + j * n] * ac;
            }
        return f;
    }

    // Euclidean gradient  dE/dC = 2 * A * C
    std::vector<double> grad(const std::vector<double>& C) const
    {
        std::vector<double> G(n * p, 0.0);
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < n; ++i)
            {
                double s = 0.0;
                for (int k = 0; k < n; ++k)
                    s += A[i + k * n] * C[k + j * n];
                G[i + j * n] = 2.0 * s;
            }
        return G;
    }
};

// Build a diagonal SPD matrix (analytically known minimum eigenvalues).
Rayleigh make_rayleigh_diag(int n, int p)
{
    Rayleigh r;
    r.n = n;
    r.p = p;
    r.A.assign(n * n, 0.0);
    // eigenvalues 1, 2, ..., n
    for (int i = 0; i < n; ++i) r.A[i + i * n] = static_cast<double>(i + 1);
    return r;
}

// Reference minimum:  sum of the p smallest eigenvalues of A.
double rayleigh_min_diag(int n, int p)
{
    (void)n;
    double s = 0.0;
    for (int i = 0; i < p; ++i) s += static_cast<double>(i + 1);
    return s;
}

// Run one Riemannian-optimiser step using the EuclideanOptimizer on the
// flattened vector.  Returns the new energy.
//
// This helper encodes exactly the same flow that optimize_orbitals() uses
// in the RDMFT solver, so the tests are a faithful unit test of the
// integration (minus the LCAO infrastructure).
template <typename OPT>
double run_stiefel_optimizer(Rayleigh& problem, OPT& opt, int iters,
                             double ls_alpha_init,
                             const std::vector<double>& C0,
                             std::vector<double>& C_out,
                             bool adam_unit_step)
{
    StiefelManifold<double> manifold(problem.n, problem.p);
    C_out = C0;
    std::vector<double> G, G_R, dir_flat, dir_proj, step_vec, C_trial;

    for (int iter = 0; iter < iters; ++iter)
    {
        G = problem.grad(C_out);
        G_R.assign(problem.n * problem.p, 0.0);
        manifold.project_tangent(C_out.data(), G.data(), G_R.data(),
                                  problem.n, problem.p);

        double gnorm2 = 0.0;
        for (auto v : G_R) gnorm2 += v * v;
        if (std::sqrt(gnorm2) < 1e-8) break;

        opt.compute_direction(G_R, dir_flat);

        // Project direction back onto the tangent space (vector transport
        // by projection), matching the solver's optimize_orbitals() flow.
        dir_proj.assign(problem.n * problem.p, 0.0);
        manifold.project_tangent(C_out.data(), dir_flat.data(),
                                  dir_proj.data(), problem.n, problem.p);

        // Descent safeguard.
        double dd = 0.0;
        for (int i = 0; i < problem.n * problem.p; ++i)
            dd += dir_proj[i] * G_R[i];
        if (dd >= 0.0)
        {
            for (int i = 0; i < problem.n * problem.p; ++i)
                dir_proj[i] = -G_R[i];
            dd = 0.0;
            for (int i = 0; i < problem.n * problem.p; ++i)
                dd += dir_proj[i] * G_R[i];
        }

        double E_cur = problem.eval(C_out);
        auto f_at = [&](double alpha) -> double {
            C_trial.assign(problem.n * problem.p, 0.0);
            manifold.retract(C_out.data(), dir_proj.data(), alpha,
                              C_trial.data(), problem.n, problem.p);
            return problem.eval(C_trial);
        };

        double alpha_init = adam_unit_step ? 1.0 : ls_alpha_init;
        auto ls = armijo_line_search(f_at, E_cur, dd, alpha_init,
                                      1e-4, 0.5, 30);

        C_trial.assign(problem.n * problem.p, 0.0);
        manifold.retract(C_out.data(), dir_proj.data(), ls.step,
                          C_trial.data(), problem.n, problem.p);
        C_out = C_trial;

        // Update optimiser state.  Adam already fully updated its moments
        // inside compute_direction(); calling update() again would double-
        // increment its step counter and corrupt the bias correction.
        if (!adam_unit_step)
        {
            std::vector<double> G_new = problem.grad(C_out);
            std::vector<double> G_R_new(problem.n * problem.p, 0.0);
            manifold.project_tangent(C_out.data(), G_new.data(),
                                      G_R_new.data(), problem.n, problem.p);
            step_vec.resize(dir_proj.size());
            for (size_t i = 0; i < step_vec.size(); ++i)
                step_vec[i] = ls.step * dir_proj[i];
            opt.update(G_R_new, step_vec);
        }
    }
    return problem.eval(C_out);
}
} // namespace

TEST_F(OptimizerTest, SD_converges_on_Stiefel_Rayleigh)
{
    const int n = 8, p = 2;
    auto problem = make_rayleigh_diag(n, p);
    double E_ref = rayleigh_min_diag(n, p);

    std::mt19937 rng(42);
    auto C0 = random_stiefel_point(n, p, rng);

    RDMFTConfig config;
    EuclideanOptimizer opt(OptimizerType::SteepestDescent, config);
    opt.init(n * p);

    std::vector<double> C;
    double E = run_stiefel_optimizer(problem, opt, 200, 0.1, C0, C, false);
    EXPECT_NEAR(E, E_ref, 1e-4);
}

TEST_F(OptimizerTest, LBFGS_converges_on_Stiefel_Rayleigh)
{
    const int n = 8, p = 2;
    auto problem = make_rayleigh_diag(n, p);
    double E_ref = rayleigh_min_diag(n, p);

    std::mt19937 rng(42);
    auto C0 = random_stiefel_point(n, p, rng);

    RDMFTConfig config;
    config.lbfgs_memory = 8;
    EuclideanOptimizer opt(OptimizerType::LBFGS, config);
    opt.init(n * p);

    std::vector<double> C;
    double E = run_stiefel_optimizer(problem, opt, 200, 0.1, C0, C, false);
    EXPECT_NEAR(E, E_ref, 1e-6);
}

TEST_F(OptimizerTest, Adam_converges_on_Stiefel_Rayleigh)
{
    const int n = 8, p = 2;
    auto problem = make_rayleigh_diag(n, p);
    double E_ref = rayleigh_min_diag(n, p);

    std::mt19937 rng(42);
    auto C0 = random_stiefel_point(n, p, rng);

    RDMFTConfig config;
    config.adam_lr = 0.05;
    EuclideanOptimizer opt(OptimizerType::Adam, config);
    opt.init(n * p);

    std::vector<double> C;
    double E = run_stiefel_optimizer(problem, opt, 2000, 1.0, C0, C, true);
    // Adam is a first-order stochastic-style optimiser; accept a looser
    // tolerance than L-BFGS.
    EXPECT_NEAR(E, E_ref, 1e-2);
}

TEST_F(OptimizerTest, LBFGS_and_SD_both_reach_Stiefel_minimum)
{
    // Both SD and L-BFGS should converge to the invariant-subspace minimum
    // of the Rayleigh quotient on the Stiefel manifold.  At the reported
    // iteration budget L-BFGS typically reaches the minimum in far fewer
    // evaluations, but the Riemannian retraction introduces O(alpha^2)
    // noise that makes a strict "L-BFGS <= SD" comparison unreliable for
    // small problems.  We test instead that both converge to within the
    // same tolerance.
    const int n = 10, p = 3;
    auto problem = make_rayleigh_diag(n, p);
    double E_ref = rayleigh_min_diag(n, p);

    std::mt19937 rng(7);
    auto C0 = random_stiefel_point(n, p, rng);

    RDMFTConfig config;
    EuclideanOptimizer sd(OptimizerType::SteepestDescent, config);
    sd.init(n * p);
    EuclideanOptimizer lbfgs(OptimizerType::LBFGS, config);
    lbfgs.init(n * p);

    std::vector<double> C_sd, C_lb;
    double E_sd = run_stiefel_optimizer(problem, sd, 400, 0.1, C0, C_sd, false);
    double E_lb = run_stiefel_optimizer(problem, lbfgs, 400, 0.1, C0, C_lb, false);

    EXPECT_NEAR(E_sd, E_ref, 1e-5);
    EXPECT_NEAR(E_lb, E_ref, 1e-5);
}

TEST_F(OptimizerTest, BBStep_BB1_and_BB2_match_formula)
{
    // x_prev -> x, g_prev -> g
    // s = [1, 0], y = [2, 0]
    // BB1 = (sTs)/(sTy) = 1/2
    // BB2 = (sTy)/(yTy) = 2/4 = 1/2
    std::vector<double> x_prev = {0.0, 0.0};
    std::vector<double> g_prev = {0.0, 0.0};
    std::vector<double> x = {1.0, 0.0};
    std::vector<double> g = {2.0, 0.0};

    BarzilaiBorweinStep bb;
    bb.set_bounds(1e-8, 10.0);
    bb.record_state(x_prev, g_prev);

    bb.set_mode(BBStepMode::BB1);
    const double a1 = bb.suggest(x, g, 0.1);
    EXPECT_NEAR(a1, 0.5, 1e-14);

    bb.set_mode(BBStepMode::BB2);
    const double a2 = bb.suggest(x, g, 0.1);
    EXPECT_NEAR(a2, 0.5, 1e-14);
}

TEST_F(OptimizerTest, BBStep_Alternate_toggles_BB1_BB2)
{
    BarzilaiBorweinStep bb;
    bb.set_mode(BBStepMode::Alternate);
    bb.set_bounds(1e-8, 10.0);

    // State 0
    bb.record_state({0.0, 0.0}, {0.0, 0.0});

    // First query (BB1): s=[1,0], y=[2,0] => 0.5
    const double a_first = bb.suggest({1.0, 0.0}, {2.0, 0.0}, 0.1);
    EXPECT_NEAR(a_first, 0.5, 1e-14);

    // Update recorded state to that accepted point.
    bb.record_state({1.0, 0.0}, {2.0, 0.0});

    // Second query (BB2): s=[1,0], y=[1,0] => 1.0
    const double a_second = bb.suggest({2.0, 0.0}, {3.0, 0.0}, 0.1);
    EXPECT_NEAR(a_second, 1.0, 1e-14);
}
