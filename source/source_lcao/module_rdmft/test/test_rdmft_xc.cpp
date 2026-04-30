#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_type.h"
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

TEST_F(XCFunctionalTest, GU_diag_factor_deriv)
{
    XCFunctional xc(XCFunctionalType::GU);
    double n = 0.35;
    EXPECT_NEAR(xc.gu_diag_factor_deriv(n), 2.0 * n - 1.0, 1e-12);
}

TEST_F(XCFunctionalTest, GU_uses_muller_alpha_for_g)
{
    XCFunctional xc(XCFunctionalType::GU);
    EXPECT_NEAR(xc.alpha(), 0.5, 1e-12);
    EXPECT_NEAR(xc.g(0.25), 0.5, 1e-12);
}

TEST_F(XCFunctionalTest, BBC3_parse_and_string)
{
    EXPECT_EQ(parse_xc_type("bbc3"), XCFunctionalType::BBC3);
    EXPECT_EQ(xc_type_to_string(XCFunctionalType::BBC3), "bbc3");
}

TEST_F(XCFunctionalTest, BBC3_is_nonseparable)
{
    XCFunctional xc(XCFunctionalType::BBC3);
    EXPECT_FALSE(xc.is_separable());
}

TEST_F(XCFunctionalTest, GEO_parse_and_string)
{
    EXPECT_EQ(parse_xc_type("geo"), XCFunctionalType::GEO);
    EXPECT_EQ(xc_type_to_string(XCFunctionalType::GEO), "geo");
}

TEST_F(XCFunctionalTest, GEO_is_nonseparable)
{
    XCFunctional xc(XCFunctionalType::GEO);
    EXPECT_FALSE(xc.is_separable());
}

TEST_F(XCFunctionalTest, GEO_decomposition_constants)
{
    EXPECT_EQ(XCFunctional::num_geo_terms(), 3);
    EXPECT_DOUBLE_EQ(XCFunctional::geo_coef(0), 0.25);
    EXPECT_DOUBLE_EQ(XCFunctional::geo_coef(1), 0.25);
    EXPECT_DOUBLE_EQ(XCFunctional::geo_coef(2), 0.50);
    EXPECT_DOUBLE_EQ(XCFunctional::geo_alpha(0), 1.0);
    EXPECT_DOUBLE_EQ(XCFunctional::geo_alpha(1), 0.5);
    EXPECT_DOUBLE_EQ(XCFunctional::geo_alpha(2), 0.75);
    // Coefficients sum to 1 so f^GEO(n,n) = n at the boundary n=1.
    double sum = 0.0;
    for (int t = 0; t < XCFunctional::num_geo_terms(); ++t)
        sum += XCFunctional::geo_coef(t);
    EXPECT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(XCFunctionalTest, GEO_full_coupling_formula)
{
    XCFunctional xc(XCFunctionalType::GEO);
    // f(n_p, n_q) = [n_p n_q + (n_p n_q)^{1/2} + 2 (n_p n_q)^{3/4}] / 4
    for (double np : {0.1, 0.3, 0.7})
    {
        for (double nq : {0.2, 0.5, 0.9})
        {
            const double prod = np * nq;
            const double expected = 0.25 * prod
                                    + 0.25 * std::sqrt(prod)
                                    + 0.5 * std::pow(prod, 0.75);
            EXPECT_NEAR(xc.f(np, nq, false), expected, 1e-12)
                << "n_p=" << np << " n_q=" << nq;
        }
    }
}

TEST_F(XCFunctionalTest, GEO_boundary_values)
{
    XCFunctional xc(XCFunctionalType::GEO);
    EXPECT_NEAR(xc.f(1.0, 1.0, false), 1.0, 1e-12);
    // f(0, n_q) is bounded by O(eps^{1/2}) because of the linear extrapolation used by
    // pow_reg below the regularisation cutoff (eps=1e-8 by default → ~1e-4·n_q^{1/2}).
    EXPECT_NEAR(xc.f(0.0, 0.5, false), 0.0, 1e-3);
    EXPECT_NEAR(xc.f(0.5, 0.0, false), 0.0, 1e-3);
}

TEST_F(XCFunctionalTest, GEO_symmetric_in_arguments)
{
    XCFunctional xc(XCFunctionalType::GEO);
    EXPECT_NEAR(xc.f(0.4, 0.7, false), xc.f(0.7, 0.4, false), 1e-14);
}

TEST_F(XCFunctionalTest, GEO_df_dni_finite_difference)
{
    XCFunctional xc(XCFunctionalType::GEO);
    const double eps = 1e-7;
    for (double np : {0.2, 0.5, 0.8})
    {
        for (double nq : {0.3, 0.6})
        {
            const double analytic = xc.df_dni(np, nq, false);
            const double numeric = (xc.f(np + eps, nq, false) - xc.f(np - eps, nq, false)) / (2.0 * eps);
            EXPECT_NEAR(analytic, numeric, 1e-5) << "n_p=" << np << " n_q=" << nq;
        }
    }
}

TEST_F(XCFunctionalTest, GEO_pow_reg_matches_pow_above_eps)
{
    XCFunctional xc(XCFunctionalType::GEO, 0.75, 1e-8);
    for (double n : {0.1, 0.5, 0.9})
    {
        EXPECT_NEAR(xc.pow_reg(n, 1.0), n, 1e-14);
        EXPECT_NEAR(xc.pow_reg(n, 0.5), std::sqrt(n), 1e-12);
        EXPECT_NEAR(xc.pow_reg(n, 0.75), std::pow(n, 0.75), 1e-12);
        EXPECT_NEAR(xc.dpow_reg(n, 1.0), 1.0, 1e-14);
        EXPECT_NEAR(xc.dpow_reg(n, 0.5), 0.5 * std::pow(n, -0.5), 1e-10);
        EXPECT_NEAR(xc.dpow_reg(n, 0.75), 0.75 * std::pow(n, -0.25), 1e-10);
    }
}

TEST_F(XCFunctionalTest, GEO_pow_reg_bounded_at_zero)
{
    XCFunctional xc(XCFunctionalType::GEO, 0.75, 1e-8);
    EXPECT_FALSE(std::isnan(xc.pow_reg(0.0, 0.5)));
    EXPECT_FALSE(std::isinf(xc.dpow_reg(0.0, 0.5)));
    EXPECT_FALSE(std::isnan(xc.pow_reg(0.0, 0.75)));
    EXPECT_FALSE(std::isinf(xc.dpow_reg(0.0, 0.75)));
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
