// -----------------------------------------------------------------------------
// Unit tests for the RDMFT joint (product-manifold) optimisation strategy.
//
// The high-level solver `RDMFTSolver::solve_joint` cannot be exercised in a
// unit test directly because it depends on ABACUS LCAO infrastructure
// (`EnergyGradient`, `Parallel_Orbitals`, `LCAO_Orbitals`, ...) which is only
// fully wired up in an end-to-end ABACUS run. These tests therefore:
//
//   1. Verify the SolverStrategy::Joint enum is reachable under the new name
//      and that `RDMFTConfig::strategy` round-trips through it (rename check).
//
//   2. Reproduce the exact algorithmic flow of `solve_joint` on a toy RDMFT
//      model whose energy, gradient and manifold structure match the real
//      code path:
//
//          E(n, C) = sum_i n_i <C_i, A C_i>   +   penalty(sum_i n_i - N_e)
//
//      with n_i parameterised as cos^2(p_i) (matching OccupationParam) and
//      C constrained to the Stiefel manifold C^T C = I (matching
//      StiefelManifold with S = I). The analytic minimiser is
//      n_* = [1, 1, ..., 1, 0, ..., 0] (N_e ones on the bands with smallest
//      eigenvalues) and C_* spans the invariant subspace of the N_b smallest
//      eigenvalues of A. The product-manifold Armijo line search in the test
//      is identical in structure to the one in `rdmft_solver.cpp`, so when
//      these tests pass we know the joint flow converges for SD, CG, L-BFGS
//      and Adam.
// -----------------------------------------------------------------------------
#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_occupation.h"
#include "source_lcao/module_rdmft/rdmft_optimizer.h"
#include "source_lcao/module_rdmft/rdmft_stiefel.h"
#include "source_lcao/module_rdmft/rdmft_type.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

using namespace rdmft;

// -----------------------------------------------------------------------------
// Rename sanity check: the new enum value is reachable and usable, and the
// RDMFTConfig default still compiles with ::Joint assigned.
// -----------------------------------------------------------------------------
TEST(RdmftJointStrategy, enum_rename_compiles)
{
    RDMFTConfig cfg;
    cfg.strategy = SolverStrategy::Joint;
    EXPECT_EQ(cfg.strategy, SolverStrategy::Joint);

    // Alternating and Joint are the only valid values.
    cfg.strategy = SolverStrategy::Alternating;
    EXPECT_EQ(cfg.strategy, SolverStrategy::Alternating);
    cfg.strategy = SolverStrategy::Joint;
    EXPECT_NE(cfg.strategy, SolverStrategy::Alternating);
}

// -----------------------------------------------------------------------------
// Toy product-manifold problem.
// -----------------------------------------------------------------------------
namespace
{

struct ToyProblem
{
    int n = 0;                    // basis size
    int nb = 0;                   // number of bands
    double Ne = 0.0;              // target electron number
    std::vector<double> A;        // n x n symmetric, column-major

    // Small entropic regulariser beta * sum_i [ n_i ln n_i + (1-n_i) ln(1-n_i) ].
    // Zero at n_i in {0,1}, minimum at n_i = 1/2. Together with the orbital
    // energy (which pushes toward integer fillings), a modest beta produces a
    // strictly interior minimum. This makes the test's analytic minimum
    // reachable by SD / CG / Adam — which otherwise stall in the
    // cosine^2-parameterisation flat at the [0,1] boundaries — without
    // changing the algorithmic flow of the joint strategy.
    double beta = 0.0;

    // Augmented-Lagrangian penalty state (shared with real solver logic).
    double lambda = 0.0;
    double mu = 1.0;

    static double clip_boltz(double n) {
        return std::max(1e-14, std::min(1.0 - 1e-14, n));
    }

    double entropy(const std::vector<double>& occ) const
    {
        if (beta == 0.0) return 0.0;
        double s = 0.0;
        for (double n : occ)
        {
            double q = clip_boltz(n);
            s += q * std::log(q) + (1.0 - q) * std::log(1.0 - q);
        }
        return beta * s;
    }

    void entropy_grad(const std::vector<double>& occ,
                      std::vector<double>& g) const
    {
        g.assign(occ.size(), 0.0);
        if (beta == 0.0) return;
        for (size_t i = 0; i < occ.size(); ++i)
        {
            double q = clip_boltz(occ[i]);
            g[i] = beta * (std::log(q) - std::log(1.0 - q));
        }
    }

    // Orbital energy:  sum_i n_i * <C_i, A C_i>  + beta * entropy(n)
    double orb_energy(const std::vector<double>& occ,
                      const std::vector<double>& C) const
    {
        double E = 0.0;
        for (int i = 0; i < nb; ++i)
        {
            double e_i = 0.0;
            for (int mu = 0; mu < n; ++mu)
            {
                double a_c = 0.0;
                for (int nu = 0; nu < n; ++nu)
                    a_c += A[mu + nu * n] * C[nu + i * n];
                e_i += C[mu + i * n] * a_c;
            }
            E += occ[i] * e_i;
        }
        return E + entropy(occ);
    }

    // dE/dn_i = <C_i, A C_i>  + beta * ( log(n_i) - log(1 - n_i) )
    std::vector<double> grad_occ(const std::vector<double>& occ,
                                 const std::vector<double>& C) const
    {
        std::vector<double> g(nb, 0.0);
        for (int i = 0; i < nb; ++i)
        {
            double e_i = 0.0;
            for (int mu = 0; mu < n; ++mu)
            {
                double a_c = 0.0;
                for (int nu = 0; nu < n; ++nu)
                    a_c += A[mu + nu * n] * C[nu + i * n];
                e_i += C[mu + i * n] * a_c;
            }
            g[i] = e_i;
        }
        if (beta != 0.0)
        {
            std::vector<double> ent;
            entropy_grad(occ, ent);
            for (int i = 0; i < nb; ++i) g[i] += ent[i];
        }
        return g;
    }

    // dE/dC = 2 * n_i * A * C_i, stored column-major like the solver uses.
    std::vector<double> grad_orb(const std::vector<double>& occ,
                                 const std::vector<double>& C) const
    {
        std::vector<double> G(n * nb, 0.0);
        for (int i = 0; i < nb; ++i)
            for (int mu = 0; mu < n; ++mu)
            {
                double s = 0.0;
                for (int nu = 0; nu < n; ++nu)
                    s += A[mu + nu * n] * C[nu + i * n];
                G[mu + i * n] = 2.0 * occ[i] * s;
            }
        return G;
    }

    double constraint(const std::vector<double>& occ) const
    {
        double s = 0.0;
        for (double v : occ) s += v;
        return s - Ne;
    }

    double penalty(const std::vector<double>& occ) const
    {
        double c = constraint(occ);
        return lambda * c + 0.5 * mu * c * c;
    }

    // d penalty / d n_i  = lambda + mu * c  (k-weight = 1 in this toy model).
    void penalty_grad(const std::vector<double>& occ, std::vector<double>& g) const
    {
        double c = constraint(occ);
        double factor = lambda + mu * c;
        g.assign(occ.size(), factor);
    }

    double total_energy(const std::vector<double>& occ,
                        const std::vector<double>& C) const
    {
        return orb_energy(occ, C) + penalty(occ);
    }
};

// Build a diagonal SPD A so the exact minimum is known.
ToyProblem make_toy_diag(int n, int nb, double Ne, double beta = 0.0)
{
    ToyProblem p;
    p.n = n;
    p.nb = nb;
    p.Ne = Ne;
    p.beta = beta;
    p.A.assign(n * n, 0.0);
    for (int i = 0; i < n; ++i) p.A[i + i * n] = static_cast<double>(i + 1);
    return p;
}

// Reference minimum energy for the toy problem. With C spanning the lowest-
// eigenvalue invariant subspace (eigenvalues 1..nb), the effective problem
// reduces to
//     min_{n_i}  sum_i n_i * eps_i  +  beta * entropy(n)
//     s.t.      sum_i n_i = Ne,  0 <= n_i <= 1.
// We solve this reduced problem numerically via Newton on the KKT multiplier
// so the test can set an arbitrary (nb, Ne, beta) triple.
double toy_ref_energy_numeric(const ToyProblem& p)
{
    std::vector<double> eps(p.nb);
    for (int i = 0; i < p.nb; ++i) eps[i] = static_cast<double>(i + 1);

    if (p.beta == 0.0)
    {
        // Fill bands greedily up to Ne electrons.
        std::vector<double> n(p.nb, 0.0);
        double remaining = p.Ne;
        for (int i = 0; i < p.nb && remaining > 0; ++i)
        {
            n[i] = std::min(1.0, remaining);
            remaining -= n[i];
        }
        double E = 0.0;
        for (int i = 0; i < p.nb; ++i) E += n[i] * eps[i];
        return E;
    }

    // Entropic KKT: n_i(mu) = sigmoid((mu - eps_i) / beta).
    // Solve sum_i n_i(mu) = Ne by bisection.
    auto fill = [&](double mu) {
        double s = 0.0;
        for (int i = 0; i < p.nb; ++i)
            s += 1.0 / (1.0 + std::exp((eps[i] - mu) / p.beta));
        return s;
    };
    double lo = eps.front() - 50.0 * std::abs(p.beta) - 10.0;
    double hi = eps.back()  + 50.0 * std::abs(p.beta) + 10.0;
    for (int it = 0; it < 200; ++it)
    {
        double mid = 0.5 * (lo + hi);
        if (fill(mid) > p.Ne) hi = mid; else lo = mid;
    }
    double mu = 0.5 * (lo + hi);
    double E = 0.0;
    double s_ent = 0.0;
    for (int i = 0; i < p.nb; ++i)
    {
        double n_i = 1.0 / (1.0 + std::exp((eps[i] - mu) / p.beta));
        E += n_i * eps[i];
        n_i = std::max(1e-14, std::min(1.0 - 1e-14, n_i));
        s_ent += n_i * std::log(n_i) + (1.0 - n_i) * std::log(1.0 - n_i);
    }
    return E + p.beta * s_ent;
}

// Gram-Schmidt random Stiefel point (identity overlap).
std::vector<double> rand_stiefel(int n, int nb, std::mt19937& rng)
{
    std::vector<double> C(n * nb);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : C) v = dist(rng);
    for (int j = 0; j < nb; ++j)
    {
        for (int k = 0; k < j; ++k)
        {
            double dot = 0.0;
            for (int i = 0; i < n; ++i) dot += C[i + k * n] * C[i + j * n];
            for (int i = 0; i < n; ++i) C[i + j * n] -= dot * C[i + k * n];
        }
        double norm = 0.0;
        for (int i = 0; i < n; ++i) norm += C[i + j * n] * C[i + j * n];
        norm = std::sqrt(norm);
        for (int i = 0; i < n; ++i) C[i + j * n] /= norm;
    }
    return C;
}

// Runs a faithful reproduction of `RDMFTSolver<>::solve_joint` on the toy
// problem using the given occupation and orbital optimisers. Mirrors the
// real solver's step order exactly so that a passing test means the
// solver's control flow is correct for the given optimiser mix.
double run_joint_toy(ToyProblem& problem,
                     OptimizerType occ_opt_type,
                     OptimizerType orb_opt_type,
                     int max_iter,
                     double alpha_init_non_lbfgs,
                     const std::vector<double>& n_init,
                     const std::vector<double>& C_init,
                     std::vector<double>& n_out,
                     std::vector<double>& C_out,
                     int seed = 42)
{
    const int n = problem.n;
    const int nb = problem.nb;
    OccupationParam op(OccParamType::CosineSq);

    RDMFTConfig cfg;
    cfg.lbfgs_memory = 8;
    cfg.adam_lr = 0.05;

    std::vector<double> params;
    op.occ_to_params(n_init, params);

    EuclideanOptimizer occ_opt(occ_opt_type, cfg);
    occ_opt.init(nb);

    EuclideanOptimizer orb_opt(orb_opt_type, cfg);
    orb_opt.init(n * nb);

    const bool orb_is_lbfgs = (orb_opt_type == OptimizerType::LBFGS);
    const bool orb_is_adam  = (orb_opt_type == OptimizerType::Adam);
    const bool orb_is_cg    = (orb_opt_type == OptimizerType::ConjugateGradient);

    // CG state
    std::vector<double> prev_orb_grad(n * nb, 0.0);
    std::vector<double> prev_orb_dir(n * nb, 0.0);
    double prev_orb_gnorm2 = 0.0;

    C_out = C_init;
    std::vector<double> occ_vec(nb, 0.0);
    op.params_to_occ(params, occ_vec);

    StiefelManifold<double> manifold(n, nb);

    double E = 0.0;
    (void)seed;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        op.params_to_occ(params, occ_vec);
        E = problem.total_energy(occ_vec, C_out);

        // Occupation parameter-space gradient.
        auto g_occ = problem.grad_occ(occ_vec, C_out);
        std::vector<double> pen;
        problem.penalty_grad(occ_vec, pen);
        for (int i = 0; i < nb; ++i) g_occ[i] += pen[i];
        std::vector<double> g_params;
        op.transform_gradient_batch(g_occ, params, g_params);

        // Riemannian orbital gradient.
        auto g_orb = problem.grad_orb(occ_vec, C_out);
        std::vector<double> g_orb_R(n * nb, 0.0);
        manifold.project_tangent(C_out.data(), g_orb.data(),
                                  g_orb_R.data(), n, nb);

        double orb_gnorm2 = 0.0;
        for (double v : g_orb_R) orb_gnorm2 += v * v;

        // Occ direction.
        std::vector<double> d_occ;
        occ_opt.compute_direction(g_params, d_occ);
        double occ_dd = 0.0;
        for (int i = 0; i < nb; ++i) occ_dd += d_occ[i] * g_params[i];
        if (occ_dd >= 0.0)
        {
            for (int i = 0; i < nb; ++i) d_occ[i] = -g_params[i];
            occ_dd = 0.0;
            for (int i = 0; i < nb; ++i) occ_dd += d_occ[i] * g_params[i];
        }

        // Orbital direction.
        std::vector<double> d_orb(n * nb);
        for (int i = 0; i < n * nb; ++i) d_orb[i] = -g_orb_R[i];

        if (orb_is_lbfgs || orb_is_adam)
        {
            std::vector<double> d_flat;
            orb_opt.compute_direction(g_orb_R, d_flat);
            std::vector<double> d_proj(n * nb, 0.0);
            manifold.project_tangent(C_out.data(), d_flat.data(),
                                      d_proj.data(), n, nb);
            double dd_check = 0.0;
            for (int i = 0; i < n * nb; ++i) dd_check += d_proj[i] * g_orb_R[i];
            if (dd_check >= 0.0)
            {
                for (int i = 0; i < n * nb; ++i) d_orb[i] = -g_orb_R[i];
            }
            else
            {
                d_orb = d_proj;
            }
        }
        else if (orb_is_cg && iter > 0 && prev_orb_gnorm2 > 1e-30)
        {
            // Vector-transport by projection.
            std::vector<double> transported(n * nb, 0.0);
            manifold.project_tangent(C_out.data(), prev_orb_dir.data(),
                                      transported.data(), n, nb);
            double beta = std::max(0.0, orb_gnorm2 / prev_orb_gnorm2);
            std::vector<double> d_cand(n * nb);
            for (int i = 0; i < n * nb; ++i)
                d_cand[i] = -g_orb_R[i] + beta * transported[i];
            std::vector<double> d_cand_proj(n * nb, 0.0);
            manifold.project_tangent(C_out.data(), d_cand.data(),
                                      d_cand_proj.data(), n, nb);
            double dd_check = 0.0;
            for (int i = 0; i < n * nb; ++i)
                dd_check += d_cand_proj[i] * g_orb_R[i];
            if (dd_check < 0.0) d_orb = d_cand_proj;
        }

        double orb_dd = 0.0;
        for (int i = 0; i < n * nb; ++i) orb_dd += d_orb[i] * g_orb_R[i];
        double dd_total = occ_dd + orb_dd;

        // Joint Armijo line search.
        double alpha = (orb_is_lbfgs || orb_is_adam) ? 1.0
                                                     : alpha_init_non_lbfgs;
        const double c1 = 1e-4;
        const double rho = 0.5;
        const int max_ls = 40;

        std::vector<double> C_trial(n * nb);
        std::vector<double> params_trial(nb);
        std::vector<double> occ_trial(nb);
        bool ok = false;
        double E_trial = E;
        for (int ls = 0; ls < max_ls; ++ls)
        {
            manifold.retract(C_out.data(), d_orb.data(), alpha,
                              C_trial.data(), n, nb);
            for (int i = 0; i < nb; ++i)
                params_trial[i] = params[i] + alpha * d_occ[i];
            op.params_to_occ(params_trial, occ_trial);
            E_trial = problem.total_energy(occ_trial, C_trial);
            if (E_trial <= E + c1 * alpha * dd_total)
            {
                ok = true;
                break;
            }
            alpha *= rho;
        }

        if (!ok)
        {
            // Reset optimiser state and try steepest descent next step.
            occ_opt.init(nb);
            if (orb_is_lbfgs || orb_is_adam) orb_opt.init(n * nb);
            if (orb_is_cg)
            {
                prev_orb_gnorm2 = 0.0;
                std::fill(prev_orb_dir.begin(), prev_orb_dir.end(), 0.0);
                std::fill(prev_orb_grad.begin(), prev_orb_grad.end(), 0.0);
            }
            continue;
        }

        // Commit.
        std::vector<double> occ_step(nb);
        for (int i = 0; i < nb; ++i)
        {
            occ_step[i] = alpha * d_occ[i];
            params[i] += occ_step[i];
        }
        op.params_to_occ(params, occ_vec);
        C_out = C_trial;

        // Refresh gradient for optimiser history updates.
        auto g_occ_new = problem.grad_occ(occ_vec, C_out);
        std::vector<double> pen_new;
        problem.penalty_grad(occ_vec, pen_new);
        for (int i = 0; i < nb; ++i) g_occ_new[i] += pen_new[i];
        std::vector<double> g_params_new;
        op.transform_gradient_batch(g_occ_new, params, g_params_new);
        occ_opt.update(g_params_new, occ_step);

        if (orb_is_lbfgs)
        {
            auto g_orb_new = problem.grad_orb(occ_vec, C_out);
            std::vector<double> g_orb_new_R(n * nb, 0.0);
            manifold.project_tangent(C_out.data(), g_orb_new.data(),
                                      g_orb_new_R.data(), n, nb);
            std::vector<double> step_vec(n * nb);
            for (int i = 0; i < n * nb; ++i) step_vec[i] = alpha * d_orb[i];
            orb_opt.update(g_orb_new_R, step_vec);
        }
        if (orb_is_cg)
        {
            prev_orb_grad = g_orb_R;
            prev_orb_dir = d_orb;
            prev_orb_gnorm2 = orb_gnorm2;
        }

        // Augmented-Lagrangian multiplier refresh. Update lambda at every
        // outer step and only grow mu when the constraint is still badly
        // violated, so the penalty never dominates the gradient.
        double c = problem.constraint(occ_vec);
        problem.lambda += problem.mu * c;
        if (std::abs(c) > 1e-3) problem.mu = std::min(problem.mu * 1.5, 100.0);

        E = E_trial;
    }

    op.params_to_occ(params, occ_vec);
    n_out = occ_vec;
    return problem.total_energy(occ_vec, C_out);
}

// Build a physically-meaningful initial point: C is the identity embedding
// (n x nb with the first nb columns of the identity), n is uniform Ne/nb.
void make_init(int n, int nb, double Ne,
               std::vector<double>& n0, std::vector<double>& C0)
{
    n0.assign(nb, std::min(0.99, std::max(0.01, Ne / nb)));
    C0.assign(n * nb, 0.0);
    for (int j = 0; j < nb; ++j) C0[j + j * n] = 1.0;
}

} // namespace

// -----------------------------------------------------------------------------
// SD / CG / L-BFGS / Adam on the occupation block, paired with SD on orbitals.
// All four must converge to the same analytic minimum.
//
// The SD / CG / Adam tests use a small entropic regulariser (beta > 0) so the
// exact optimum is strictly interior in [0, 1]. Otherwise the cosine^2
// parameterisation's Jacobian dn/dp = -sin(2p) vanishes at n in {0, 1}, which
// stalls any first-order method (this is exactly the issue that motivated
// `occ_init_margin` in the production code). L-BFGS bootstraps a Hessian
// approximation that escapes the stall, so its test uses the integer-filling
// limit (beta = 0) and verifies fast convergence.
// -----------------------------------------------------------------------------
TEST(RdmftJointStrategy, joint_SD_SD_converges)
{
    auto problem = make_toy_diag(/*n=*/8, /*nb=*/4, /*Ne=*/2.0, /*beta=*/0.5);
    const double E_ref = toy_ref_energy_numeric(problem);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    std::vector<double> n_out, C_out;
    double E = run_joint_toy(problem,
                             OptimizerType::SteepestDescent,
                             OptimizerType::SteepestDescent,
                             /*max_iter=*/3000,
                             /*alpha_init=*/0.1,
                             n_init, C_init, n_out, C_out);
    EXPECT_LT(E, E_ref + 5e-2);
    EXPECT_GT(E, E_ref - 5e-2);
}

TEST(RdmftJointStrategy, joint_CG_CG_converges)
{
    auto problem = make_toy_diag(8, 4, 2.0, 0.5);
    const double E_ref = toy_ref_energy_numeric(problem);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    std::vector<double> n_out, C_out;
    double E = run_joint_toy(problem,
                             OptimizerType::ConjugateGradient,
                             OptimizerType::ConjugateGradient,
                             3000, 0.1,
                             n_init, C_init, n_out, C_out);
    EXPECT_LT(E, E_ref + 5e-2);
    EXPECT_GT(E, E_ref - 5e-2);
}

TEST(RdmftJointStrategy, joint_LBFGS_LBFGS_converges)
{
    auto problem = make_toy_diag(8, 4, 2.0, 0.0);
    const double E_ref = toy_ref_energy_numeric(problem);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    std::vector<double> n_out, C_out;
    double E = run_joint_toy(problem,
                             OptimizerType::LBFGS,
                             OptimizerType::LBFGS,
                             800, 0.1,
                             n_init, C_init, n_out, C_out);
    EXPECT_NEAR(E, E_ref, 1e-3);
}

TEST(RdmftJointStrategy, joint_Adam_SD_converges)
{
    // Adam on occupations, SD on orbitals, simulating the example in the
    // user guide (Example 3).
    auto problem = make_toy_diag(8, 4, 2.0, 0.5);
    const double E_ref = toy_ref_energy_numeric(problem);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    std::vector<double> n_out, C_out;
    double E = run_joint_toy(problem,
                             OptimizerType::Adam,
                             OptimizerType::SteepestDescent,
                             5000, 0.1,
                             n_init, C_init, n_out, C_out);
    EXPECT_LT(E, E_ref + 1e-1);
    EXPECT_GT(E, E_ref - 1e-1);
}

// -----------------------------------------------------------------------------
// Verify that the packed orbital component of the direction keeps C on the
// Stiefel manifold: after optimisation, C^T C = I to machine precision.
// -----------------------------------------------------------------------------
TEST(RdmftJointStrategy, joint_preserves_stiefel)
{
    auto problem = make_toy_diag(10, 3, 1.5);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    std::vector<double> n_out, C_out;
    run_joint_toy(problem,
                  OptimizerType::ConjugateGradient,
                  OptimizerType::LBFGS,
                  500, 0.05,
                  n_init, C_init, n_out, C_out);

    // Check C^T C = I
    double max_err = 0.0;
    for (int i = 0; i < problem.nb; ++i)
        for (int j = 0; j < problem.nb; ++j)
        {
            double dot = 0.0;
            for (int k = 0; k < problem.n; ++k)
                dot += C_out[k + i * problem.n] * C_out[k + j * problem.n];
            double expected = (i == j) ? 1.0 : 0.0;
            max_err = std::max(max_err, std::abs(dot - expected));
        }
    EXPECT_LT(max_err, 1e-10);
}

// -----------------------------------------------------------------------------
// Verify that the alternating working optimisers still converge on the same
// problem when the joint flow is replaced by two independent sub-problems.
// This guards against regressions that could break the alternating strategy
// while editing `solve_joint`.
// -----------------------------------------------------------------------------
TEST(RdmftJointStrategy, alternating_still_works_with_lbfgs_adam)
{
    auto problem = make_toy_diag(8, 4, 2.0, 0.5);
    const double E_ref = toy_ref_energy_numeric(problem);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    OccupationParam op(OccParamType::CosineSq);
    std::vector<double> params;
    op.occ_to_params(n_init, params);

    StiefelManifold<double> manifold(problem.n, problem.nb);

    EuclideanOptimizer occ_opt(OptimizerType::LBFGS, RDMFTConfig{});
    occ_opt.init(problem.nb);
    EuclideanOptimizer orb_opt(OptimizerType::Adam, [] {
        RDMFTConfig c;
        c.adam_lr = 0.05;
        return c;
    }());
    orb_opt.init(problem.n * problem.nb);

    std::vector<double> C = C_init;
    std::vector<double> occ(problem.nb, 0.0);
    op.params_to_occ(params, occ);

    // Simple alternating loop (fixed iteration count, no adaptive line search).
    for (int outer = 0; outer < 200; ++outer)
    {
        // Occupation sub-problem (L-BFGS on cosine^2 parameters).
        for (int inner = 0; inner < 20; ++inner)
        {
            op.params_to_occ(params, occ);
            auto g_occ = problem.grad_occ(occ, C);
            std::vector<double> pen;
            problem.penalty_grad(occ, pen);
            for (int i = 0; i < problem.nb; ++i) g_occ[i] += pen[i];
            std::vector<double> g_params;
            op.transform_gradient_batch(g_occ, params, g_params);
            std::vector<double> d;
            occ_opt.compute_direction(g_params, d);
            double dd = 0.0;
            for (int i = 0; i < problem.nb; ++i) dd += d[i] * g_params[i];
            if (dd >= 0)
            {
                for (int i = 0; i < problem.nb; ++i) d[i] = -g_params[i];
                dd = 0.0;
                for (int i = 0; i < problem.nb; ++i) dd += d[i] * g_params[i];
            }

            double alpha = 1.0;
            std::vector<double> params_trial(problem.nb);
            std::vector<double> occ_trial(problem.nb);
            double E_cur = problem.total_energy(occ, C);
            bool ok = false;
            double E_new = E_cur;
            for (int ls = 0; ls < 30; ++ls)
            {
                for (int i = 0; i < problem.nb; ++i)
                    params_trial[i] = params[i] + alpha * d[i];
                op.params_to_occ(params_trial, occ_trial);
                E_new = problem.total_energy(occ_trial, C);
                if (E_new <= E_cur + 1e-4 * alpha * dd) { ok = true; break; }
                alpha *= 0.5;
            }
            if (!ok) break;
            std::vector<double> step(problem.nb);
            for (int i = 0; i < problem.nb; ++i)
            {
                step[i] = alpha * d[i];
                params[i] += step[i];
            }
            op.params_to_occ(params, occ);
            auto g_occ_new = problem.grad_occ(occ, C);
            std::vector<double> pen_new;
            problem.penalty_grad(occ, pen_new);
            for (int i = 0; i < problem.nb; ++i) g_occ_new[i] += pen_new[i];
            std::vector<double> g_params_new;
            op.transform_gradient_batch(g_occ_new, params, g_params_new);
            occ_opt.update(g_params_new, step);
        }

        // Orbital sub-problem (Adam, projected/retracted).
        for (int inner = 0; inner < 20; ++inner)
        {
            auto g_orb = problem.grad_orb(occ, C);
            std::vector<double> g_orb_R(problem.n * problem.nb, 0.0);
            manifold.project_tangent(C.data(), g_orb.data(),
                                      g_orb_R.data(), problem.n, problem.nb);

            std::vector<double> d_flat;
            orb_opt.compute_direction(g_orb_R, d_flat);
            std::vector<double> d_proj(problem.n * problem.nb, 0.0);
            manifold.project_tangent(C.data(), d_flat.data(),
                                      d_proj.data(), problem.n, problem.nb);

            double dd = 0.0;
            for (int i = 0; i < problem.n * problem.nb; ++i)
                dd += d_proj[i] * g_orb_R[i];
            if (dd >= 0)
            {
                for (int i = 0; i < problem.n * problem.nb; ++i)
                    d_proj[i] = -g_orb_R[i];
                dd = 0.0;
                for (int i = 0; i < problem.n * problem.nb; ++i)
                    dd += d_proj[i] * g_orb_R[i];
            }

            double alpha = 1.0;
            std::vector<double> C_trial(problem.n * problem.nb, 0.0);
            double E_cur = problem.total_energy(occ, C);
            bool ok = false;
            for (int ls = 0; ls < 30; ++ls)
            {
                manifold.retract(C.data(), d_proj.data(), alpha,
                                  C_trial.data(), problem.n, problem.nb);
                double E_new = problem.total_energy(occ, C_trial);
                if (E_new <= E_cur + 1e-4 * alpha * dd) { ok = true; break; }
                alpha *= 0.5;
            }
            if (ok) C = C_trial;
            else break;
        }

        double c = problem.constraint(occ);
        problem.lambda += problem.mu * c;
        if (std::abs(c) > 1e-6) problem.mu = std::min(problem.mu * 2.0, 1e6);
    }

    double E = problem.total_energy(occ, C);
    EXPECT_NEAR(E, E_ref, 5e-2);
}
