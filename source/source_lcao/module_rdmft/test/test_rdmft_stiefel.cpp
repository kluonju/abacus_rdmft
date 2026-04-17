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
