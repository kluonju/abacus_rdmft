#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_stiefel.h"
#include <cmath>
#include <vector>
#include <random>

using namespace rdmft;

class StiefelTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        rng_.seed(42);
    }

    std::mt19937 rng_;

    // Generate a random orthogonal matrix (columns of the Stiefel manifold)
    std::vector<double> random_stiefel_point(int n, int p)
    {
        std::vector<double> C(n * p);
        std::normal_distribution<double> dist(0.0, 1.0);
        for (auto& v : C) v = dist(rng_);

        // Gram-Schmidt
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

    // Check C^T C = I
    double orthogonality_error(const double* C, int n, int p)
    {
        double max_err = 0.0;
        for (int i = 0; i < p; ++i)
            for (int j = 0; j < p; ++j)
            {
                double dot = 0.0;
                for (int k = 0; k < n; ++k)
                    dot += C[k + i * n] * C[k + j * n];
                double expected = (i == j) ? 1.0 : 0.0;
                max_err = std::max(max_err, std::abs(dot - expected));
            }
        return max_err;
    }
};

TEST_F(StiefelTest, tangent_projection_preserves_tangent_condition)
{
    int n = 8, p = 3;
    StiefelManifold<double> st(n, p);

    auto C = random_stiefel_point(n, p);

    // Random ambient vector
    std::vector<double> G(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : G) v = dist(rng_);

    std::vector<double> proj(n * p);
    st.project_tangent(C.data(), G.data(), proj.data(), n, p);

    // Check tangent condition: C^T * proj + proj^T * C = 0
    double max_err = 0.0;
    for (int i = 0; i < p; ++i)
        for (int j = 0; j <= i; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < n; ++k)
                sum += C[k + i * n] * proj[k + j * n] + proj[k + i * n] * C[k + j * n];
            max_err = std::max(max_err, std::abs(sum));
        }
    EXPECT_LT(max_err, 1e-10);
}

TEST_F(StiefelTest, retraction_preserves_orthogonality)
{
    int n = 8, p = 3;
    StiefelManifold<double> st(n, p);

    auto C = random_stiefel_point(n, p);
    EXPECT_LT(orthogonality_error(C.data(), n, p), 1e-12);

    // Create a tangent vector
    std::vector<double> G(n * p);
    std::normal_distribution<double> dist(0.0, 0.1);
    for (auto& v : G) v = dist(rng_);
    std::vector<double> eta(n * p);
    st.project_tangent(C.data(), G.data(), eta.data(), n, p);

    // Retract
    std::vector<double> C_new(n * p);
    st.retract(C.data(), eta.data(), 0.1, C_new.data(), n, p);

    EXPECT_LT(orthogonality_error(C_new.data(), n, p), 1e-10);
}

TEST_F(StiefelTest, inner_product_symmetric)
{
    int n = 6, p = 2;
    StiefelManifold<double> st(n, p);
    auto C = random_stiefel_point(n, p);

    std::vector<double> eta1(n * p), eta2(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : eta1) v = dist(rng_);
    for (auto& v : eta2) v = dist(rng_);

    std::vector<double> proj1(n * p), proj2(n * p);
    st.project_tangent(C.data(), eta1.data(), proj1.data(), n, p);
    st.project_tangent(C.data(), eta2.data(), proj2.data(), n, p);

    double ip12 = st.inner_product(proj1.data(), proj2.data(), n, p);
    double ip21 = st.inner_product(proj2.data(), proj1.data(), n, p);
    EXPECT_NEAR(ip12, ip21, 1e-12);
}

TEST_F(StiefelTest, norm_positive_definite)
{
    int n = 6, p = 2;
    StiefelManifold<double> st(n, p);
    auto C = random_stiefel_point(n, p);

    std::vector<double> eta(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : eta) v = dist(rng_);

    std::vector<double> proj(n * p);
    st.project_tangent(C.data(), eta.data(), proj.data(), n, p);

    double norm = st.norm(proj.data(), n, p);
    EXPECT_GT(norm, 0.0);

    std::vector<double> zero(n * p, 0.0);
    EXPECT_NEAR(st.norm(zero.data(), n, p), 0.0, 1e-15);
}

TEST_F(StiefelTest, reorthogonalize_identity_overlap)
{
    int n = 6, p = 3;
    StiefelManifold<double> st(n, p);

    // Start with a non-orthogonal matrix
    std::vector<double> Y(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : Y) v = dist(rng_);

    std::vector<double> C(n * p);
    st.reorthogonalize(Y.data(), C.data(), n, p);

    EXPECT_LT(orthogonality_error(C.data(), n, p), 1e-10);
}

// ============================================================================
// Generalized Stiefel manifold tests with S != I.
//
// We build a simple diagonal positive-definite overlap matrix S = diag(1,2,...,n)
// and verify that:
//   (a) reorthogonalize produces C^H S C = I
//   (b) project_tangent satisfies the generalised tangent condition
//       C^H S proj + proj^H S C = 0
//   (c) retract preserves C^H S C = I
//   (d) a vector in the generalised normal space (G = C * K, K symmetric)
//       projects to zero
//   (e) the projection is idempotent: proj(proj(G)) = proj(G)
// ============================================================================

namespace
{

// Build a diagonal overlap matrix: S_{ii} = (i+1) for i = 0..n-1
std::vector<double> make_diag_overlap(int n)
{
    std::vector<double> S(n * n, 0.0);
    for (int i = 0; i < n; ++i)
        S[i + i * n] = static_cast<double>(i + 1);
    return S;
}

// Build a tridiagonal SPD overlap: S_{ii} = 2, S_{i,i±1} = -0.5
std::vector<double> make_tridiag_overlap(int n)
{
    std::vector<double> S(n * n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        S[i + i * n] = 2.0;
        if (i > 0)     S[i + (i - 1) * n] = -0.5;
        if (i < n - 1) S[i + (i + 1) * n] = -0.5;
    }
    return S;
}

// Maximum entry of |C^H S C - I|
double gen_orthogonality_error(const double* C, const double* S,
                               int n, int p)
{
    double max_err = 0.0;
    for (int i = 0; i < p; ++i)
        for (int j = 0; j < p; ++j)
        {
            double val = 0.0;
            for (int a = 0; a < n; ++a)
                for (int b = 0; b < n; ++b)
                    val += C[a + i * n] * S[a + b * n] * C[b + j * n];
            double expected = (i == j) ? 1.0 : 0.0;
            max_err = std::max(max_err, std::abs(val - expected));
        }
    return max_err;
}

// Maximum entry of |C^H S G + G^H S C| (generalised tangent condition check)
double gen_tangent_condition_error(const double* C, const double* G,
                                   const double* S, int n, int p)
{
    double max_err = 0.0;
    for (int i = 0; i < p; ++i)
        for (int j = 0; j < p; ++j)
        {
            // (C^H S G)_{ij} = sum_{a,b} C_{ai} S_{ab} G_{bj}
            double chsg = 0.0;
            for (int a = 0; a < n; ++a)
                for (int b = 0; b < n; ++b)
                    chsg += C[a + i * n] * S[a + b * n] * G[b + j * n];
            // The (i,j) entry of C^H S G + G^H S C = chsg + conj(chsg_{ji})
            // For real: this is chsg + (G^H S C)_{ij} = chsg + chsg_transposed
            double ghsc = 0.0;
            for (int a = 0; a < n; ++a)
                for (int b = 0; b < n; ++b)
                    ghsc += G[a + i * n] * S[a + b * n] * C[b + j * n];
            max_err = std::max(max_err, std::abs(chsg + ghsc));
        }
    return max_err;
}

// Generate a starting point on St(p,n;S): begin with random Y, then reorthogonalize
std::vector<double> random_gen_stiefel_point(int n, int p,
                                              const std::vector<double>& S,
                                              std::mt19937& rng)
{
    StiefelManifold<double> st(n, p, S.data());
    std::vector<double> Y(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : Y) v = dist(rng);
    std::vector<double> C(n * p);
    st.reorthogonalize(Y.data(), C.data(), n, p);
    return C;
}

} // namespace

TEST_F(StiefelTest, gen_reorthogonalize_diag_overlap)
{
    // Verify C^H S C = I after reorthogonalizing a random matrix w.r.t.
    // a diagonal overlap S = diag(1,2,...,n).
    const int n = 8, p = 3;
    auto S = make_diag_overlap(n);
    StiefelManifold<double> st(n, p, S.data());

    std::vector<double> Y(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : Y) v = dist(rng_);

    std::vector<double> C(n * p);
    st.reorthogonalize(Y.data(), C.data(), n, p);

    EXPECT_LT(gen_orthogonality_error(C.data(), S.data(), n, p), 1e-10);
}

TEST_F(StiefelTest, gen_reorthogonalize_tridiag_overlap)
{
    // Same check with a tridiagonal SPD overlap.
    const int n = 6, p = 2;
    auto S = make_tridiag_overlap(n);
    StiefelManifold<double> st(n, p, S.data());

    std::vector<double> Y(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : Y) v = dist(rng_);

    std::vector<double> C(n * p);
    st.reorthogonalize(Y.data(), C.data(), n, p);

    EXPECT_LT(gen_orthogonality_error(C.data(), S.data(), n, p), 1e-10);
}

TEST_F(StiefelTest, gen_project_tangent_satisfies_tangent_condition)
{
    // G_R = project_tangent(C, G) must satisfy the generalised tangent condition:
    //   C^H S G_R + G_R^H S C = 0
    const int n = 8, p = 3;
    auto S = make_diag_overlap(n);
    StiefelManifold<double> st(n, p, S.data());

    auto C = random_gen_stiefel_point(n, p, S, rng_);

    std::vector<double> G(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : G) v = dist(rng_);

    std::vector<double> proj(n * p);
    st.project_tangent(C.data(), G.data(), proj.data(), n, p);

    EXPECT_LT(gen_tangent_condition_error(C.data(), proj.data(), S.data(), n, p), 1e-10);
}

TEST_F(StiefelTest, gen_project_tangent_tridiag_overlap)
{
    // Same tangent-condition check for the tridiagonal overlap.
    const int n = 7, p = 2;
    auto S = make_tridiag_overlap(n);
    StiefelManifold<double> st(n, p, S.data());

    auto C = random_gen_stiefel_point(n, p, S, rng_);

    std::vector<double> G(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : G) v = dist(rng_);

    std::vector<double> proj(n * p);
    st.project_tangent(C.data(), G.data(), proj.data(), n, p);

    EXPECT_LT(gen_tangent_condition_error(C.data(), proj.data(), S.data(), n, p), 1e-10);
}

TEST_F(StiefelTest, gen_project_normal_space_gives_zero)
{
    // A vector in the generalised normal space N_C = { C*K : K = K^T }
    // should project exactly to zero.
    //
    // Derivation: G_R = G - C * sym(C^H S G).
    // For G = C * K with K symmetric (real: K = K^T):
    //   C^H S G = C^H S C K = I * K = K  (since C^H S C = I)
    //   sym(C^H S G) = sym(K) = K  (K is already symmetric)
    //   G_R = C*K - C*K = 0
    const int n = 8, p = 3;
    auto S = make_diag_overlap(n);
    StiefelManifold<double> st(n, p, S.data());

    auto C = random_gen_stiefel_point(n, p, S, rng_);

    // Build a random symmetric K (p x p) and form G = C * K
    std::normal_distribution<double> dist(0.0, 1.0);
    std::vector<double> K(p * p, 0.0);
    for (int i = 0; i < p; ++i)
        for (int j = i; j < p; ++j)
        {
            double v = dist(rng_);
            K[i + j * p] = v;
            K[j + i * p] = v;
        }
    // G = C * K  (n x p)
    std::vector<double> G(n * p, 0.0);
    for (int a = 0; a < n; ++a)
        for (int j = 0; j < p; ++j)
            for (int l = 0; l < p; ++l)
                G[a + j * n] += C[a + l * n] * K[l + j * p];

    std::vector<double> proj(n * p, 0.0);
    st.project_tangent(C.data(), G.data(), proj.data(), n, p);

    double max_err = 0.0;
    for (double v : proj) max_err = std::max(max_err, std::abs(v));
    EXPECT_LT(max_err, 1e-10);
}

TEST_F(StiefelTest, gen_project_tangent_is_idempotent)
{
    // Applying project_tangent twice should return the same result:
    //   proj(proj(G)) = proj(G)
    const int n = 8, p = 3;
    auto S = make_diag_overlap(n);
    StiefelManifold<double> st(n, p, S.data());

    auto C = random_gen_stiefel_point(n, p, S, rng_);

    std::vector<double> G(n * p);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (auto& v : G) v = dist(rng_);

    std::vector<double> proj1(n * p), proj2(n * p);
    st.project_tangent(C.data(), G.data(), proj1.data(), n, p);
    st.project_tangent(C.data(), proj1.data(), proj2.data(), n, p);

    double max_err = 0.0;
    for (int i = 0; i < n * p; ++i)
        max_err = std::max(max_err, std::abs(proj1[i] - proj2[i]));
    EXPECT_LT(max_err, 1e-10);
}

TEST_F(StiefelTest, gen_retract_preserves_gen_orthogonality)
{
    // retract(C, eta, step) must return a point on St(p,n;S), i.e. C^H S C = I.
    const int n = 8, p = 3;
    auto S = make_diag_overlap(n);
    StiefelManifold<double> st(n, p, S.data());

    auto C = random_gen_stiefel_point(n, p, S, rng_);
    EXPECT_LT(gen_orthogonality_error(C.data(), S.data(), n, p), 1e-10);

    // Build a tangent vector eta
    std::vector<double> G(n * p);
    std::normal_distribution<double> dist(0.0, 0.1);
    for (auto& v : G) v = dist(rng_);
    std::vector<double> eta(n * p);
    st.project_tangent(C.data(), G.data(), eta.data(), n, p);

    // Retract with a non-trivial step
    std::vector<double> C_new(n * p);
    st.retract(C.data(), eta.data(), 1.0, C_new.data(), n, p);

    EXPECT_LT(gen_orthogonality_error(C_new.data(), S.data(), n, p), 1e-10);
}

TEST_F(StiefelTest, gen_s_fidelity_is_zero_on_manifold)
{
    // After reorthogonalization w.r.t. S, Tr(C^H S C) / p == 1,
    // so the S-fidelity deviation Tr(C^H S C)/p - 1 should be 0.
    // This mirrors the corrected compute_s_fidelity_trace formula in rdmft_solver.cpp.
    const int n = 8, p = 3;
    auto S = make_diag_overlap(n);
    StiefelManifold<double> st(n, p, S.data());

    auto C = random_gen_stiefel_point(n, p, S, rng_);

    // Compute Tr(C^H S C)
    double trace_chsc = 0.0;
    for (int j = 0; j < p; ++j)
    {
        // diagonal entry j: sum_{a,b} C_{aj} S_{ab} C_{bj}
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                trace_chsc += C[a + j * n] * S[a + b * n] * C[b + j * n];
    }
    // s_fidelity_deviation = trace / p - 1 should be 0
    double s_fidelity = trace_chsc / static_cast<double>(p) - 1.0;
    EXPECT_NEAR(s_fidelity, 0.0, 1e-10);
}

