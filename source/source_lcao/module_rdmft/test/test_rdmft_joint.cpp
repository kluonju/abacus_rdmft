// -----------------------------------------------------------------------------
// Unit tests for the RDMFT joint (product-manifold) optimisation strategy.
//
// The high-level solver `RDMFTSolver::solve_joint` cannot be exercised in a
// unit test directly because it depends on ABACUS LCAO infrastructure
// (`EnergyGradient`, `Parallel_Orbitals`, `LCAO_Orbitals`, ...) which is only
// fully wired up in an end-to-end ABACUS run. These tests therefore:
//
//   1. Verify the SolverStrategy::Joint enum is reachable under the new name
//      and that `RDMFTConfig::strategy` / `RDMFTConfig::joint_optimizer`
//      round-trip through it.
//
//   2. Reproduce the exact algorithmic flow of the refactored `solve_joint`
//      on a toy RDMFT model whose energy, gradient and manifold structure
//      match the real code path:
//
//          E(n, C) = sum_i n_i <C_i, A C_i>   +   penalty(sum_i n_i - N_e)
//
//      with n_i parameterised as cos^2(p_i) (matching OccupationParam) and
//      C constrained to the Stiefel manifold C^T C = I (matching
//      StiefelManifold with S = I). A SINGLE EuclideanOptimizer (sd or cg)
//      is applied to the packed vector z = (p_1, ..., p_nb, flat(C)) using
//      the packed gradient g = (dE/dp, flat(G_R)). The produced direction is
//      split back into occupation and orbital blocks, the orbital block is
//      re-projected onto the Stiefel tangent space, and one joint Armijo line
//      search commits the step (linear update for p, retraction for C).
// -----------------------------------------------------------------------------
#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_occupation.h"
#include "source_lcao/module_rdmft/rdmft_optimizer.h"
#include "source_lcao/module_rdmft/test/test_stiefel_helper.h"
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
// The joint strategy uses a single unified optimiser; RDMFTConfig exposes it.
// -----------------------------------------------------------------------------
TEST(RdmftJointStrategy, joint_optimizer_roundtrips)
{
    RDMFTConfig cfg;
    EXPECT_EQ(cfg.joint_optimizer, OptimizerType::ConjugateGradient);

    for (OptimizerType t : {OptimizerType::SteepestDescent, OptimizerType::ConjugateGradient})
    {
        cfg.joint_optimizer = t;
        EXPECT_EQ(cfg.joint_optimizer, t);
    }
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
    // reachable by SD / CG -- which otherwise stall in the
    // cosine^2-parameterisation flat at the [0,1] boundaries -- without
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

// Runs a faithful reproduction of the refactored `RDMFTSolver<>::solve_joint`
// on the toy problem using a SINGLE unified Euclidean optimiser applied to
// the packed variable z = (p, flat(C)). Mirrors the real solver's step order:
//   1. Evaluate energy + Euclidean (dE/dn, dE/dC) at current z.
//   2. Chain-rule dE/dn -> dE/dp; project dE/dC -> Riemannian gradient G_R.
//   3. Pack (dE/dp, flat(G_R)) -> packed_grad, call optimiser for packed_dir.
//   4. Split packed_dir; re-project orbital block onto tangent space.
//   5. Descent safeguard on joint directional derivative.
//   6. Single Armijo line search along packed direction: params += alpha * d_p,
//      C = retract(C, d_C, alpha).
//   7. Update optimiser state with the new packed gradient (cg history).
//   8. Refresh augmented-Lagrangian multiplier.
double run_joint_toy_single(ToyProblem& problem,
                            OptimizerType joint_opt_type,
                            int max_iter,
                            double alpha_init,
                            const std::vector<double>& n_init,
                            const std::vector<double>& C_init,
                            std::vector<double>& n_out,
                            std::vector<double>& C_out)
{
    const int n = problem.n;
    const int nb = problem.nb;
    OccupationParam op(OccParamType::CosineSq);

    RDMFTConfig cfg;

    std::vector<double> params;
    op.occ_to_params(n_init, params);

    const int n_occ = nb;
    const int n_orb = n * nb;
    const int packed_size = n_occ + n_orb;

    EuclideanOptimizer joint_opt(joint_opt_type, cfg);
    joint_opt.init(packed_size);

    C_out = C_init;
    std::vector<double> occ_vec(nb, 0.0);
    op.params_to_occ(params, occ_vec);

    StiefelManifold<double> manifold(n, nb);

    double E = 0.0;
    for (int iter = 0; iter < max_iter; ++iter)
    {
        op.params_to_occ(params, occ_vec);
        E = problem.total_energy(occ_vec, C_out);

        // Occupation parameter-space gradient (chain-rule dE/dn -> dE/dp).
        auto g_occ = problem.grad_occ(occ_vec, C_out);
        std::vector<double> pen;
        problem.penalty_grad(occ_vec, pen);
        for (int i = 0; i < nb; ++i) g_occ[i] += pen[i];
        std::vector<double> g_params;
        op.transform_gradient_batch(g_occ, params, g_params);

        // Riemannian orbital gradient (project Euclidean dE/dC onto T_C St).
        auto g_orb = problem.grad_orb(occ_vec, C_out);
        std::vector<double> g_orb_R(n_orb, 0.0);
        manifold.project_tangent(C_out.data(), g_orb.data(),
                                  g_orb_R.data(), n, nb);

        double orb_gnorm2 = 0.0;
        for (double v : g_orb_R) orb_gnorm2 += v * v;

        // Pack gradients for the single optimiser.
        std::vector<double> packed_grad(packed_size);
        for (int i = 0; i < n_occ; ++i) packed_grad[i] = g_params[i];
        for (int i = 0; i < n_orb; ++i) packed_grad[n_occ + i] = g_orb_R[i];

        // One call to the unified optimiser.
        std::vector<double> packed_dir;
        joint_opt.compute_direction(packed_grad, packed_dir);

        // Split into occupation and orbital blocks.
        std::vector<double> d_occ(n_occ);
        for (int i = 0; i < n_occ; ++i) d_occ[i] = packed_dir[i];
        std::vector<double> d_orb(n_orb);
        for (int i = 0; i < n_orb; ++i) d_orb[i] = packed_dir[n_occ + i];

        // Re-project orbital block onto tangent space at C.
        std::vector<double> d_orb_proj(n_orb, 0.0);
        manifold.project_tangent(C_out.data(), d_orb.data(),
                                  d_orb_proj.data(), n, nb);
        d_orb = d_orb_proj;

        // Joint directional derivative.
        double dd_total = 0.0;
        for (int i = 0; i < n_occ; ++i) dd_total += d_occ[i] * g_params[i];
        for (int i = 0; i < n_orb; ++i) dd_total += d_orb[i] * g_orb_R[i];

        // Descent safeguard on the packed direction.
        if (dd_total >= 0.0)
        {
            for (int i = 0; i < n_occ; ++i) d_occ[i] = -g_params[i];
            for (int i = 0; i < n_orb; ++i) d_orb[i] = -g_orb_R[i];
            dd_total = 0.0;
            for (int i = 0; i < n_occ; ++i) dd_total += d_occ[i] * g_params[i];
            for (int i = 0; i < n_orb; ++i) dd_total += d_orb[i] * g_orb_R[i];
        }

        // Joint Armijo line search.
        double alpha = alpha_init;
        const double c1 = 1e-4;
        const double rho = 0.5;
        const int max_ls = 40;

        std::vector<double> C_trial(n_orb);
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
            joint_opt.init(packed_size);
            continue;
        }

        // Commit.
        std::vector<double> packed_step(packed_size);
        for (int i = 0; i < nb; ++i)
        {
            packed_step[i] = alpha * d_occ[i];
            params[i] += packed_step[i];
        }
        for (int i = 0; i < n_orb; ++i)
            packed_step[n_occ + i] = alpha * d_orb[i];
        op.params_to_occ(params, occ_vec);
        C_out = C_trial;

        if (joint_opt_type == OptimizerType::ConjugateGradient)
        {
            // Refresh packed gradient for CG state (matches RDMFTSolver joint path).
            auto g_occ_new = problem.grad_occ(occ_vec, C_out);
            std::vector<double> pen_new;
            problem.penalty_grad(occ_vec, pen_new);
            for (int i = 0; i < nb; ++i) g_occ_new[i] += pen_new[i];
            std::vector<double> g_params_new;
            op.transform_gradient_batch(g_occ_new, params, g_params_new);

            auto g_orb_new = problem.grad_orb(occ_vec, C_out);
            std::vector<double> g_orb_new_R(n_orb, 0.0);
            manifold.project_tangent(C_out.data(), g_orb_new.data(),
                                      g_orb_new_R.data(), n, nb);

            std::vector<double> packed_grad_new(packed_size);
            for (int i = 0; i < n_occ; ++i) packed_grad_new[i] = g_params_new[i];
            for (int i = 0; i < n_orb; ++i) packed_grad_new[n_occ + i] = g_orb_new_R[i];
            joint_opt.update(packed_grad_new, packed_step);
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
// A single unified optimiser on the packed vector (p, C_flat) converges to
// the analytic minimum of the toy product-manifold problem for SD and CG.
//
// SD / CG use a small entropic regulariser (beta > 0) so the exact optimum is
// strictly interior in [0, 1]. Otherwise the cosine^2 parameterisation's
// Jacobian dn/dp = -sin(2p) vanishes at n in {0, 1}, which stalls first-order
// methods in the integer-filling limit.
// -----------------------------------------------------------------------------
TEST(RdmftJointStrategy, joint_single_SD_converges)
{
    auto problem = make_toy_diag(/*n=*/8, /*nb=*/4, /*Ne=*/2.0, /*beta=*/0.5);
    const double E_ref = toy_ref_energy_numeric(problem);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    std::vector<double> n_out, C_out;
    double E = run_joint_toy_single(problem,
                                    OptimizerType::SteepestDescent,
                                    /*max_iter=*/3000,
                                    /*alpha_init=*/0.1,
                                    n_init, C_init, n_out, C_out);
    EXPECT_LT(E, E_ref + 5e-2);
    EXPECT_GT(E, E_ref - 5e-2);
}

TEST(RdmftJointStrategy, joint_single_CG_converges)
{
    auto problem = make_toy_diag(8, 4, 2.0, 0.5);
    const double E_ref = toy_ref_energy_numeric(problem);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    std::vector<double> n_out, C_out;
    double E = run_joint_toy_single(problem,
                                    OptimizerType::ConjugateGradient,
                                    3000, 0.1,
                                    n_init, C_init, n_out, C_out);
    EXPECT_LT(E, E_ref + 5e-2);
    EXPECT_GT(E, E_ref - 5e-2);
}

// -----------------------------------------------------------------------------
// Verify that the packed orbital component of the direction keeps C on the
// Stiefel manifold: after optimisation, C^T C = I to machine precision.
// -----------------------------------------------------------------------------
TEST(RdmftJointStrategy, joint_single_preserves_stiefel)
{
    auto problem = make_toy_diag(10, 3, 1.5);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    std::vector<double> n_out, C_out;
    run_joint_toy_single(problem,
                         OptimizerType::ConjugateGradient,
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
// The packed gradient must be computed so the joint directional derivative
// along -grad equals -(||grad_p||^2 + ||grad_C||_R^2). This verifies that
// steepest descent in the packed space yields a proper descent direction on
// the product manifold.
// -----------------------------------------------------------------------------
TEST(RdmftJointStrategy, packed_gradient_is_descent)
{
    auto problem = make_toy_diag(8, 3, 1.5, 0.1);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    OccupationParam op(OccParamType::CosineSq);
    std::vector<double> params;
    op.occ_to_params(n_init, params);
    std::vector<double> occ(problem.nb, 0.0);
    op.params_to_occ(params, occ);

    StiefelManifold<double> manifold(problem.n, problem.nb);

    auto g_occ = problem.grad_occ(occ, C_init);
    std::vector<double> pen;
    problem.penalty_grad(occ, pen);
    for (int i = 0; i < problem.nb; ++i) g_occ[i] += pen[i];
    std::vector<double> g_params;
    op.transform_gradient_batch(g_occ, params, g_params);

    auto g_orb = problem.grad_orb(occ, C_init);
    std::vector<double> g_orb_R(problem.n * problem.nb, 0.0);
    manifold.project_tangent(C_init.data(), g_orb.data(),
                              g_orb_R.data(), problem.n, problem.nb);

    double gp2 = 0.0;
    for (double v : g_params) gp2 += v * v;
    double gc2 = 0.0;
    for (double v : g_orb_R) gc2 += v * v;
    const double expected_dd = -(gp2 + gc2);

    // Finite-difference check of the directional derivative along -grad.
    double alpha = 1e-6;
    std::vector<double> params_plus(problem.nb);
    for (int i = 0; i < problem.nb; ++i)
        params_plus[i] = params[i] - alpha * g_params[i];
    std::vector<double> occ_plus(problem.nb);
    op.params_to_occ(params_plus, occ_plus);
    std::vector<double> C_plus(problem.n * problem.nb);
    for (int i = 0; i < problem.n * problem.nb; ++i)
        C_plus[i] = C_init[i] - alpha * g_orb_R[i];

    const double E0 = problem.total_energy(occ, C_init);
    const double E1 = problem.total_energy(occ_plus, C_plus);
    const double dd_fd = (E1 - E0) / alpha;

    // Allow some tolerance for the retraction vs. Euclidean-step mismatch
    // and the cosine^2 nonlinearity; the sign and magnitude should match.
    EXPECT_LT(dd_fd, 0.0);
    EXPECT_NEAR(dd_fd, expected_dd, 1e-2 * std::abs(expected_dd) + 1e-8);
}

// -----------------------------------------------------------------------------
// Verify that the alternating working optimisers still converge on the same
// problem when the joint flow is replaced by two independent sub-problems.
// This guards against regressions that could break the alternating strategy
// while editing `solve_joint`.
// -----------------------------------------------------------------------------
TEST(RdmftJointStrategy, alternating_still_works_with_cg_sd)
{
    auto problem = make_toy_diag(8, 4, 2.0, 0.5);
    const double E_ref = toy_ref_energy_numeric(problem);
    std::vector<double> n_init, C_init;
    make_init(problem.n, problem.nb, problem.Ne, n_init, C_init);

    OccupationParam op(OccParamType::CosineSq);
    std::vector<double> params;
    op.occ_to_params(n_init, params);

    StiefelManifold<double> manifold(problem.n, problem.nb);

    EuclideanOptimizer occ_opt(OptimizerType::ConjugateGradient, RDMFTConfig{});
    occ_opt.init(problem.nb);
    EuclideanOptimizer orb_opt(OptimizerType::SteepestDescent, RDMFTConfig{});
    orb_opt.init(problem.n * problem.nb);

    std::vector<double> C = C_init;
    std::vector<double> occ(problem.nb, 0.0);
    op.params_to_occ(params, occ);

    // Simple alternating loop (fixed iteration count, no adaptive line search).
    for (int outer = 0; outer < 200; ++outer)
    {
        // Occupation sub-problem (CG on cosine^2 parameters).
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

        // Orbital sub-problem (SD, projected/retracted).
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
