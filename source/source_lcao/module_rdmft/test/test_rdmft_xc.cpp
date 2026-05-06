#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_xc_functional.h"
#include <cmath>

using namespace rdmft;

namespace {
// Mirrors rdmft_energy_gradient.cpp binary entropy helpers (HF occupation regularization).
constexpr double occ_entropy_n_eps = 1e-12;
double test_binary_entropy_f(double n)
{
    n = std::max(occ_entropy_n_eps, std::min(1.0 - occ_entropy_n_eps, n));
    return n * std::log(n) + (1.0 - n) * std::log(1.0 - n);
}
double test_binary_entropy_dfdn(double n)
{
    n = std::max(occ_entropy_n_eps, std::min(1.0 - occ_entropy_n_eps, n));
    return std::log(n / (1.0 - n));
}
} // namespace

class XCFunctionalTest : public ::testing::Test {};

TEST_F(XCFunctionalTest, HF_g_is_identity)
{
    XCFunctional xc(XCFunctionalType::HF);
    EXPECT_DOUBLE_EQ(xc.g(0.0), 0.0);
    EXPECT_DOUBLE_EQ(xc.g(0.5), 0.5);
    EXPECT_DOUBLE_EQ(xc.g(1.0), 1.0);
}

TEST_F(XCFunctionalTest, HF_dg_is_one)
{
    XCFunctional xc(XCFunctionalType::HF);
    EXPECT_DOUBLE_EQ(xc.dg(0.5), 1.0);
    EXPECT_DOUBLE_EQ(xc.dg(0.9), 1.0);
}

TEST_F(XCFunctionalTest, Muller_g_is_sqrt)
{
    XCFunctional xc(XCFunctionalType::Muller);
    EXPECT_NEAR(xc.g(0.25), 0.5, 1e-12);
    EXPECT_NEAR(xc.g(1.0), 1.0, 1e-12);
    EXPECT_NEAR(xc.g(0.0), 0.0, 1e-12);
}

TEST_F(XCFunctionalTest, Muller_dg_is_half_inv_sqrt)
{
    XCFunctional xc(XCFunctionalType::Muller);
    double n = 0.25;
    double expected = 0.5 / std::sqrt(n);
    EXPECT_NEAR(xc.dg(n), expected, 1e-12);
}

TEST_F(XCFunctionalTest, Power_g_with_alpha)
{
    double alpha = 0.656;
    XCFunctional xc(XCFunctionalType::Power, alpha);
    double n = 0.7;
    EXPECT_NEAR(xc.g(n), std::pow(n, alpha), 1e-12);
}

TEST_F(XCFunctionalTest, Power_dg_with_alpha)
{
    double alpha = 0.656;
    XCFunctional xc(XCFunctionalType::Power, alpha);
    double n = 0.7;
    double expected = alpha * std::pow(n, alpha - 1.0);
    EXPECT_NEAR(xc.dg(n), expected, 1e-10);
}

TEST_F(XCFunctionalTest, GU_separability)
{
    XCFunctional xc(XCFunctionalType::GU);
    EXPECT_FALSE(xc.is_separable());

    XCFunctional hf(XCFunctionalType::HF);
    EXPECT_TRUE(hf.is_separable());
}

TEST_F(XCFunctionalTest, GU_coupling_diagonal)
{
    XCFunctional xc(XCFunctionalType::GU);
    double n = 0.6;
    EXPECT_NEAR(xc.f(n, n, true), n * n, 1e-12);
}

TEST_F(XCFunctionalTest, GU_coupling_offdiag)
{
    XCFunctional xc(XCFunctionalType::GU);
    double ni = 0.4, nj = 0.9;
    EXPECT_NEAR(xc.f(ni, nj, false), std::sqrt(ni) * std::sqrt(nj), 1e-12);
}

TEST_F(XCFunctionalTest, GU_diag_factor)
{
    XCFunctional xc(XCFunctionalType::GU);
    double n = 0.7;
    EXPECT_NEAR(xc.gu_diag_factor(n), n * n - n, 1e-12);
}

TEST_F(XCFunctionalTest, gradient_consistency_numerical)
{
    // Check that dg is consistent with g via finite differences
    for (auto type : {XCFunctionalType::HF, XCFunctionalType::Muller, XCFunctionalType::Power})
    {
        XCFunctional xc(type, 0.656);
        double eps = 1e-7;
        for (double n = 0.1; n <= 0.9; n += 0.1)
        {
            double dg_analytic = xc.dg(n);
            double dg_numerical = (xc.g(n + eps) - xc.g(n - eps)) / (2.0 * eps);
            EXPECT_NEAR(dg_analytic, dg_numerical, 1e-5)
                << "type=" << static_cast<int>(type) << " n=" << n;
        }
    }
}

TEST_F(XCFunctionalTest, BinaryEntropy_interior_matches_formula)
{
    const double n = 0.25;
    const double f = test_binary_entropy_f(n);
    const double f_expected = n * std::log(n) + (1.0 - n) * std::log(1.0 - n);
    EXPECT_NEAR(f, f_expected, 1e-14);
    EXPECT_NEAR(test_binary_entropy_dfdn(n), std::log(n / (1.0 - n)), 1e-14);
}

TEST_F(XCFunctionalTest, BinaryEntropy_clamped_boundaries_are_finite)
{
    EXPECT_FALSE(std::isnan(test_binary_entropy_f(0.0)));
    EXPECT_FALSE(std::isnan(test_binary_entropy_f(1.0)));
    EXPECT_FALSE(std::isnan(test_binary_entropy_dfdn(0.0)));
    EXPECT_FALSE(std::isnan(test_binary_entropy_dfdn(1.0)));
}

TEST_F(XCFunctionalTest, boundary_values)
{
    XCFunctional xc(XCFunctionalType::Power, 0.656);
    EXPECT_NEAR(xc.g(0.0), 0.0, 1e-12);
    EXPECT_NEAR(xc.g(1.0), 1.0, 1e-12);
}
