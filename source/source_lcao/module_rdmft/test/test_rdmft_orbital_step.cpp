// -----------------------------------------------------------------------------
// Unit tests for the RDMFT orbital sub-problem:
//
//   1. Three Stiefel retractions exposed by EnergyGradient::retract_orbitals
//      via OrbRetraction = { Polar, QR, Cayley }: orthogonality is preserved
//      to machine precision after one retraction step.
//   2. The Riemannian SD/CG + monotone Armijo orbital optimiser
//      converges on a Rayleigh-quotient minimisation.
//   3. RDMFTConfig defaults, parser semantics and string round-trip.
//
// The retraction tests use the serial reference twins implemented in
// test_stiefel_helper.h (StiefelManifold<TK>::retract / retract_qr /
// retract_cayley); these mirror the math in
// EnergyGradient::retract_polar / retract_qr_serial / retract_cayley_serial
// without requiring the full LCAO infrastructure.
// -----------------------------------------------------------------------------
#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_optimizer.h"
#include "source_lcao/module_rdmft/rdmft_type.h"
#include "source_lcao/module_rdmft/test/test_stiefel_helper.h"

#include <cmath>
#include <complex>
#include <random>
#include <vector>

using namespace rdmft;

namespace
{
// Build a random orthonormal column-major nbasis x norbs matrix.
std::vector<double> random_stiefel_real(int nbasis, int norbs, std::mt19937& rng)
{
    std::vector<double> C(static_cast<size_t>(nbasis) * norbs, 0.0);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : C) v = dist(rng);
    // Modified Gram-Schmidt.
    for (int j = 0; j < norbs; ++j)
    {
        for (int k = 0; k < j; ++k)
        {
            double dot = 0.0;
            for (int i = 0; i < nbasis; ++i)
                dot += C[i + k * nbasis] * C[i + j * nbasis];
            for (int i = 0; i < nbasis; ++i)
                C[i + j * nbasis] -= dot * C[i + k * nbasis];
        }
        double norm = 0.0;
        for (int i = 0; i < nbasis; ++i)
            norm += C[i + j * nbasis] * C[i + j * nbasis];
        norm = std::sqrt(norm);
        for (int i = 0; i < nbasis; ++i)
            C[i + j * nbasis] /= norm;
    }
    return C;
}

double frob_orthonormality_error(const double* C, int nbasis, int norbs)
{
    double err2 = 0.0;
    for (int j = 0; j < norbs; ++j)
        for (int i = 0; i < norbs; ++i)
        {
            double dot = 0.0;
            for (int k = 0; k < nbasis; ++k)
                dot += C[k + i * nbasis] * C[k + j * nbasis];
            const double expected = (i == j) ? 1.0 : 0.0;
            err2 += (dot - expected) * (dot - expected);
        }
    return std::sqrt(err2);
}
} // namespace

// -----------------------------------------------------------------------------
// 1. Retraction orthonormality preservation.
// -----------------------------------------------------------------------------

TEST(RdmftOrbRetraction, polar_preserves_orthonormality)
{
    const int nbasis = 16, norbs = 4;
    StiefelManifold<double> m(nbasis, norbs);

    std::mt19937 rng(123);
    std::normal_distribution<double> dist(0.0, 1.0);

    auto X = random_stiefel_real(nbasis, norbs, rng);
    EXPECT_LT(frob_orthonormality_error(X.data(), nbasis, norbs), 1e-12);

    // Random ambient direction projected to tangent space.
    std::vector<double> G(nbasis * norbs, 0.0);
    for (auto& v : G) v = dist(rng);
    std::vector<double> Gt(nbasis * norbs, 0.0);
    m.project_tangent(X.data(), G.data(), Gt.data(), nbasis, norbs);

    // Use the existing polar retraction (StiefelManifold::retract is
    // Cholesky-QR, matching EnergyGradient::retract_polar in X-space).
    for (double step : {1e-3, 1e-2, 1e-1})
    {
        std::vector<double> Xnew(nbasis * norbs, 0.0);
        m.retract(X.data(), Gt.data(), step, Xnew.data(), nbasis, norbs);
        const double err = frob_orthonormality_error(Xnew.data(), nbasis, norbs);
        EXPECT_LT(err, 1e-12) << "polar step=" << step;
    }
}

TEST(RdmftOrbRetraction, qr_preserves_orthonormality)
{
    const int nbasis = 16, norbs = 4;
    StiefelManifold<double> m(nbasis, norbs);

    std::mt19937 rng(456);
    std::normal_distribution<double> dist(0.0, 1.0);

    auto X = random_stiefel_real(nbasis, norbs, rng);

    std::vector<double> G(nbasis * norbs, 0.0);
    for (auto& v : G) v = dist(rng);
    std::vector<double> Gt(nbasis * norbs, 0.0);
    m.project_tangent(X.data(), G.data(), Gt.data(), nbasis, norbs);

    for (double step : {1e-4, 1e-2, 1.0, 5.0})
    {
        std::vector<double> Xnew(nbasis * norbs, 0.0);
        m.retract_qr(X.data(), Gt.data(), step, Xnew.data(), nbasis, norbs);
        const double err = frob_orthonormality_error(Xnew.data(), nbasis, norbs);
        EXPECT_LT(err, 1e-12) << "qr step=" << step;
    }
}

TEST(RdmftOrbRetraction, cayley_preserves_orthonormality)
{
    const int nbasis = 16, norbs = 4;
    StiefelManifold<double> m(nbasis, norbs);

    std::mt19937 rng(789);
    std::normal_distribution<double> dist(0.0, 1.0);

    auto X = random_stiefel_real(nbasis, norbs, rng);

    std::vector<double> G(nbasis * norbs, 0.0);
    for (auto& v : G) v = dist(rng);
    std::vector<double> Gt(nbasis * norbs, 0.0);
    m.project_tangent(X.data(), G.data(), Gt.data(), nbasis, norbs);

    // Cayley is exactly orthogonality-preserving for any step (sign and
    // magnitude); the only failure mode is a singular 2p x 2p system, which
    // the implementation reports via info != 0. With small-to-moderate
    // gradients we never hit that.
    for (double step : {1e-4, 1e-2, 1e-1, 1.0})
    {
        std::vector<double> Xnew(nbasis * norbs, 0.0);
        m.retract_cayley(X.data(), Gt.data(), step, Xnew.data(), nbasis, norbs);
        const double err = frob_orthonormality_error(Xnew.data(), nbasis, norbs);
        EXPECT_LT(err, 1e-11) << "cayley step=" << step;
    }
}

TEST(RdmftOrbRetraction, qr_first_order_matches_minus_grad)
{
    // For small step the Riemannian gradient flow gives  X_new ≈ X - step * G_R.
    // The sign-fixed QR retraction must agree to first order.
    const int nbasis = 12, norbs = 3;
    StiefelManifold<double> m(nbasis, norbs);
    std::mt19937 rng(2025);
    std::normal_distribution<double> dist(0.0, 1.0);

    auto X = random_stiefel_real(nbasis, norbs, rng);
    std::vector<double> G(nbasis * norbs, 0.0);
    for (auto& v : G) v = dist(rng);
    std::vector<double> Gt(nbasis * norbs, 0.0);
    m.project_tangent(X.data(), G.data(), Gt.data(), nbasis, norbs);

    const double step = 1e-6;
    std::vector<double> Xnew(nbasis * norbs, 0.0);
    // Note: helpers' retract / retract_qr take eta with the convention
    // Y = C + step * eta. Pass the *negative* tangent so the resulting
    // Y = C - step * G_R reproduces a descent step.
    std::vector<double> mGt(Gt.size());
    for (size_t i = 0; i < Gt.size(); ++i) mGt[i] = -Gt[i];

    m.retract_qr(X.data(), mGt.data(), step, Xnew.data(), nbasis, norbs);

    // First-order check: ||X_new - X + step * G_R|| / ||step * G_R|| << 1.
    double num2 = 0.0, den2 = 0.0;
    for (int i = 0; i < nbasis * norbs; ++i)
    {
        const double e = Xnew[i] - X[i] + step * Gt[i];
        num2 += e * e;
        den2 += (step * Gt[i]) * (step * Gt[i]);
    }
    if (den2 < 1e-30) den2 = 1.0;
    EXPECT_LT(std::sqrt(num2 / den2), 1e-3);
}

// -----------------------------------------------------------------------------
// 2. Simple SD / CG + monotone Armijo on the Rayleigh quotient.
// -----------------------------------------------------------------------------

namespace
{
// E(C) = Tr(C^T A C),  C in R^{n x p},  C^T C = I.
// Minimum is the sum of the p smallest eigenvalues of A (Courant-Fischer).
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

Rayleigh make_rayleigh_diag(int n, int p)
{
    Rayleigh r;
    r.n = n;
    r.p = p;
    r.A.assign(static_cast<size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) r.A[i + i * n] = static_cast<double>(i + 1);
    return r;
}
double rayleigh_min_diag(int p)
{
    double s = 0.0;
    for (int i = 0; i < p; ++i) s += static_cast<double>(i + 1);
    return s;
}

// One iteration of Riemannian SD/CG + monotone Armijo, retraction selectable.
// Returns final energy; mirrors RDMFTSolver::optimize_orbitals (sd / cg path)
// but at the small algebraic-test scale.
enum class HelperRetraction
{
    Polar,
    QR,
    Cayley
};

double run_riemannian_sd_cg(Rayleigh& problem,
                            bool use_cg,
                            HelperRetraction retr,
                            int iters,
                            double alpha_init,
                            const std::vector<double>& C0,
                            std::vector<double>& C_out)
{
    const int n = problem.n;
    const int p = problem.p;
    StiefelManifold<double> manifold(n, p);
    C_out = C0;
    std::vector<double> G, Gt, dir, dir_prev, G_prev, C_trial;

    for (int iter = 0; iter < iters; ++iter)
    {
        G = problem.grad(C_out);
        Gt.assign(n * p, 0.0);
        manifold.project_tangent(C_out.data(), G.data(), Gt.data(), n, p);

        double gnorm2 = 0.0;
        for (auto v : Gt) gnorm2 += v * v;
        if (std::sqrt(gnorm2) < 1e-9) break;

        // SD: D = -Gt; CG: Polak-Ribiere+ with vector transport by projection.
        dir.assign(n * p, 0.0);
        for (int i = 0; i < n * p; ++i) dir[i] = -Gt[i];
        if (use_cg && !G_prev.empty())
        {
            std::vector<double> Gp_T(n * p, 0.0);
            manifold.project_tangent(C_out.data(), G_prev.data(), Gp_T.data(), n, p);
            double gn2_prev = 0.0;
            for (auto v : Gp_T) gn2_prev += v * v;
            if (gn2_prev > 1e-30)
            {
                std::vector<double> y(n * p, 0.0);
                for (int i = 0; i < n * p; ++i) y[i] = Gt[i] - Gp_T[i];
                double num = 0.0;
                for (int i = 0; i < n * p; ++i) num += Gt[i] * y[i];
                double beta = std::max(0.0, num / gn2_prev);
                std::vector<double> dp_T(n * p, 0.0);
                manifold.project_tangent(C_out.data(), dir_prev.data(), dp_T.data(), n, p);
                for (int i = 0; i < n * p; ++i) dir[i] = -Gt[i] + beta * dp_T[i];
                std::vector<double> dir_proj(n * p, 0.0);
                manifold.project_tangent(C_out.data(), dir.data(), dir_proj.data(), n, p);
                dir = dir_proj;
            }
        }
        // Descent safeguard.
        double dd = 0.0;
        for (int i = 0; i < n * p; ++i) dd += dir[i] * Gt[i];
        if (dd >= 0.0)
        {
            for (int i = 0; i < n * p; ++i) dir[i] = -Gt[i];
            dd = 0.0;
            for (int i = 0; i < n * p; ++i) dd += dir[i] * Gt[i];
        }

        const double E_cur = problem.eval(C_out);
        auto f_at = [&](double a) -> double {
            C_trial.assign(n * p, 0.0);
            switch (retr)
            {
                case HelperRetraction::Polar:
                    manifold.retract(C_out.data(), dir.data(), a,
                                     C_trial.data(), n, p);
                    break;
                case HelperRetraction::QR:
                    manifold.retract_qr(C_out.data(), dir.data(), a,
                                        C_trial.data(), n, p);
                    break;
                case HelperRetraction::Cayley:
                {
                    // Cayley uses G (gradient) directly; sign flip to map our
                    // descent direction onto the formula's convention.
                    std::vector<double> Gconv(n * p, 0.0);
                    for (int i = 0; i < n * p; ++i) Gconv[i] = -dir[i];
                    manifold.retract_cayley(C_out.data(), Gconv.data(), a,
                                            C_trial.data(), n, p);
                    break;
                }
            }
            return problem.eval(C_trial);
        };
        auto d_at = [&](double a) -> double {
            C_trial.assign(n * p, 0.0);
            switch (retr)
            {
                case HelperRetraction::Polar:
                    manifold.retract(C_out.data(), dir.data(), a,
                                     C_trial.data(), n, p);
                    break;
                case HelperRetraction::QR:
                    manifold.retract_qr(C_out.data(), dir.data(), a,
                                        C_trial.data(), n, p);
                    break;
                case HelperRetraction::Cayley:
                {
                    std::vector<double> Gconv(n * p, 0.0);
                    for (int i = 0; i < n * p; ++i) Gconv[i] = -dir[i];
                    manifold.retract_cayley(C_out.data(), Gconv.data(), a,
                                            C_trial.data(), n, p);
                    break;
                }
            }
            std::vector<double> Gn = problem.grad(C_trial);
            std::vector<double> GRn(n * p, 0.0);
            manifold.project_tangent(C_trial.data(), Gn.data(), GRn.data(), n, p);
            double s = 0.0;
            for (int i = 0; i < n * p; ++i)
                s += GRn[i] * dir[i];
            return s;
        };
        const auto ls = nonmonotone_strong_wolfe_line_search(f_at,
            d_at,
            E_cur,
            dd,
            E_cur,
            alpha_init,
            1e-4,
            0.9,
            40,
            40);
        // Commit to C_trial at the accepted (or last-tried) step.
        C_trial.assign(n * p, 0.0);
        switch (retr)
        {
            case HelperRetraction::Polar:
                manifold.retract(C_out.data(), dir.data(), ls.step,
                                 C_trial.data(), n, p);
                break;
            case HelperRetraction::QR:
                manifold.retract_qr(C_out.data(), dir.data(), ls.step,
                                    C_trial.data(), n, p);
                break;
            case HelperRetraction::Cayley:
            {
                std::vector<double> Gconv(n * p, 0.0);
                for (int i = 0; i < n * p; ++i) Gconv[i] = -dir[i];
                manifold.retract_cayley(C_out.data(), Gconv.data(), ls.step,
                                        C_trial.data(), n, p);
                break;
            }
        }
        G_prev = Gt;
        dir_prev = dir;
        C_out = C_trial;
    }
    return problem.eval(C_out);
}
} // namespace

TEST(RdmftOrbStep, sd_qr_converges_on_rayleigh)
{
    const int n = 8, p = 2;
    auto problem = make_rayleigh_diag(n, p);
    const double E_ref = rayleigh_min_diag(p); // 1 + 2 = 3.

    std::mt19937 rng(31);
    auto C0 = random_stiefel_real(n, p, rng);

    std::vector<double> C;
    const double E = run_riemannian_sd_cg(problem, /*use_cg=*/false,
                                          HelperRetraction::QR, 400, 0.1, C0, C);
    EXPECT_NEAR(E, E_ref, 1e-5);
}

TEST(RdmftOrbStep, cg_cayley_converges_on_rayleigh)
{
    const int n = 8, p = 2;
    auto problem = make_rayleigh_diag(n, p);
    const double E_ref = rayleigh_min_diag(p);

    std::mt19937 rng(73);
    auto C0 = random_stiefel_real(n, p, rng);

    std::vector<double> C;
    const double E = run_riemannian_sd_cg(problem, /*use_cg=*/true,
                                          HelperRetraction::Cayley, 400, 0.1, C0, C);
    EXPECT_NEAR(E, E_ref, 1e-5);
}

TEST(RdmftOrbStep, sd_polar_converges_on_rayleigh)
{
    const int n = 10, p = 3;
    auto problem = make_rayleigh_diag(n, p);
    const double E_ref = rayleigh_min_diag(p); // 1 + 2 + 3 = 6.

    std::mt19937 rng(101);
    auto C0 = random_stiefel_real(n, p, rng);

    std::vector<double> C;
    const double E = run_riemannian_sd_cg(problem, /*use_cg=*/false,
                                          HelperRetraction::Polar, 400, 0.1, C0, C);
    EXPECT_NEAR(E, E_ref, 1e-5);
}

// -----------------------------------------------------------------------------
// 3. RDMFTConfig defaults & parser semantics for the orbital retraction.
// -----------------------------------------------------------------------------

TEST(RdmftOrbConfig, default_orb_retraction_is_polar)
{
    RDMFTConfig cfg;
    EXPECT_EQ(cfg.orb_retraction, OrbRetraction::Polar);
    EXPECT_EQ(orb_retraction_to_string(cfg.orb_retraction), std::string("polar"));
    EXPECT_EQ(cfg.orb_optimizer, OptimizerType::ConjugateGradient);
}

TEST(RdmftOrbConfig, parse_orb_retraction_accepts_documented_values)
{
    EXPECT_EQ(parse_orb_retraction("polar"), OrbRetraction::Polar);
    EXPECT_EQ(parse_orb_retraction("qr"), OrbRetraction::QR);
    EXPECT_EQ(parse_orb_retraction("cayley"), OrbRetraction::Cayley);
    EXPECT_THROW(parse_orb_retraction("nonsense"), std::invalid_argument);
}
