//==========================================================
// Author: Kai Luo
// Email: kluo@njust.edu.cn
// DATE : April 2026
//==========================================================
#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_occupation.h"
#include <cmath>
#include <vector>

using namespace rdmft;

class OccupationParamTest : public ::testing::Test {};

TEST_F(OccupationParamTest, CosineSq_roundtrip)
{
    OccupationParam param(OccParamType::CosineSq);
    for (double n = 0.05; n <= 0.95; n += 0.1)
    {
        double p = param.to_param(n);
        double n_back = param.to_occ(p);
        EXPECT_NEAR(n, n_back, 1e-10) << "n=" << n;
    }
}

TEST_F(OccupationParamTest, Logistic_roundtrip)
{
    OccupationParam param(OccParamType::Logistic);
    for (double n = 0.05; n <= 0.95; n += 0.1)
    {
        double p = param.to_param(n);
        double n_back = param.to_occ(p);
        EXPECT_NEAR(n, n_back, 1e-10) << "n=" << n;
    }
}

TEST_F(OccupationParamTest, CosineSq_range)
{
    OccupationParam param(OccParamType::CosineSq);
    for (double p = -10.0; p <= 10.0; p += 0.5)
    {
        double n = param.to_occ(p);
        EXPECT_GE(n, 0.0);
        EXPECT_LE(n, 1.0);
    }
}

TEST_F(OccupationParamTest, Logistic_range)
{
    OccupationParam param(OccParamType::Logistic);
    for (double p = -10.0; p <= 10.0; p += 0.5)
    {
        double n = param.to_occ(p);
        EXPECT_GE(n, 0.0);
        EXPECT_LE(n, 1.0);
    }
}

TEST_F(OccupationParamTest, CosineSq_jacobian_numerical)
{
    OccupationParam param(OccParamType::CosineSq);
    double eps = 1e-7;
    for (double p = 0.2; p <= 1.4; p += 0.2)
    {
        double jac_analytic = param.jacobian(p);
        double jac_numerical = (param.to_occ(p + eps) - param.to_occ(p - eps)) / (2.0 * eps);
        EXPECT_NEAR(jac_analytic, jac_numerical, 1e-5) << "p=" << p;
    }
}

TEST_F(OccupationParamTest, Logistic_jacobian_numerical)
{
    OccupationParam param(OccParamType::Logistic);
    double eps = 1e-7;
    for (double p = -3.0; p <= 3.0; p += 0.5)
    {
        double jac_analytic = param.jacobian(p);
        double jac_numerical = (param.to_occ(p + eps) - param.to_occ(p - eps)) / (2.0 * eps);
        EXPECT_NEAR(jac_analytic, jac_numerical, 1e-5) << "p=" << p;
    }
}

TEST_F(OccupationParamTest, gradient_transform_chain_rule)
{
    OccupationParam param(OccParamType::CosineSq);
    double dE_dn = 2.5;
    double p = 0.5;
    double expected = dE_dn * param.jacobian(p);
    EXPECT_NEAR(param.transform_gradient(dE_dn, p), expected, 1e-12);
}

TEST_F(OccupationParamTest, Logistic_vector_maps_to_unit_interval)
{
    OccupationParam param(OccParamType::Logistic);
    std::vector<double> params = {-1.2, 0.1, 0.7, -0.4};
    std::vector<double> occ;
    param.params_to_occ(params, occ);
    ASSERT_EQ(occ.size(), params.size());
    for (double n : occ)
    {
        EXPECT_GE(n, 0.0);
        EXPECT_LE(n, 1.0);
    }
}

TEST_F(OccupationParamTest, Logistic_transform_gradient_batch_matches_finite_difference)
{
    OccupationParam param(OccParamType::Logistic);
    const std::vector<double> params = {-0.8, 0.2, 0.5, -0.3};
    const std::vector<double> dE_dn = {0.7, -1.1, 0.4, 1.3};

    std::vector<double> analytic;
    param.transform_gradient_batch(dE_dn, params, analytic);

    auto objective = [&](const std::vector<double>& p) {
        std::vector<double> occ;
        param.params_to_occ(p, occ);
        double value = 0.0;
        for (size_t i = 0; i < occ.size(); ++i)
            value += dE_dn[i] * occ[i];
        return value;
    };

    const double eps = 1e-7;
    for (size_t i = 0; i < params.size(); ++i)
    {
        std::vector<double> p_plus = params;
        std::vector<double> p_minus = params;
        p_plus[i] += eps;
        p_minus[i] -= eps;
        const double fd = (objective(p_plus) - objective(p_minus)) / (2.0 * eps);
        EXPECT_NEAR(analytic[i], fd, 1e-5) << "i=" << i;
    }
}


class OccupationConstraintTest : public ::testing::Test {};

TEST_F(OccupationConstraintTest, constraint_violation)
{
    std::vector<double> wk = {0.5, 0.5};
    int nbands = 3;
    double nel = 4.0;
    OccupationConstraint constraint(ConstraintMethod::AugmentedLagrangian, nel, wk, nbands);

    // All occupations = 1.0: sum = 0.5*3 + 0.5*3 = 3.0, violation = 3-4 = -1
    std::vector<double> occ(6, 1.0);
    EXPECT_NEAR(constraint.constraint_violation(occ), -1.0, 1e-12);

    // Occupations that satisfy: sum_k w_k * sum_i n = nel
    // Need 0.5 * (n1+n2+n3) + 0.5 * (n4+n5+n6) = 4
    // So total sum of n_i = 8, each n_i = 8/6 > 1, can't satisfy with box constraint
    // Use nel = 2.0 instead
    OccupationConstraint c2(ConstraintMethod::AugmentedLagrangian, 2.0, wk, nbands);
    std::vector<double> occ2 = {0.8, 0.7, 0.5, 0.8, 0.7, 0.5};
    double cv = c2.constraint_violation(occ2);
    double expected = 0.5 * (0.8 + 0.7 + 0.5) + 0.5 * (0.8 + 0.7 + 0.5) - 2.0;
    EXPECT_NEAR(cv, expected, 1e-12);
}

TEST_F(OccupationConstraintTest, augmented_lagrangian_gradient)
{
    std::vector<double> wk = {1.0};
    int nbands = 2;
    double nel = 1.5;
    OccupationConstraint constraint(ConstraintMethod::AugmentedLagrangian, nel, wk, nbands);

    std::vector<double> occ = {0.8, 0.9};
    double c = constraint.constraint_violation(occ);

    std::vector<double> grad;
    constraint.augmented_lagrangian_gradient(occ, grad);

    // grad should be (lambda + mu*c) * w_k for each element
    double factor = constraint.lambda() + constraint.mu() * c;
    EXPECT_NEAR(grad[0], factor * wk[0], 1e-12);
    EXPECT_NEAR(grad[1], factor * wk[0], 1e-12);
}

TEST_F(OccupationConstraintTest, projection_clips_and_rescales)
{
    std::vector<double> wk = {1.0};
    int nbands = 3;
    double nel = 2.0;
    OccupationConstraint constraint(ConstraintMethod::ProjectedGradient, nel, wk, nbands);

    std::vector<double> occ = {-0.1, 0.5, 1.5};
    constraint.project(occ);

    for (auto n : occ)
    {
        EXPECT_GE(n, 0.0);
        EXPECT_LE(n, 1.0);
    }
}

TEST_F(OccupationConstraintTest, projected_gradient_stop_metric_uses_fixed_tau)
{
    std::vector<double> wk = {1.0};
    int nbands = 4;
    double nel = 2.0;
    OccupationConstraint constraint(ConstraintMethod::ProjectedGradient, nel, wk, nbands);

    // Near-idempotent but still fractional occupations (sum = nel).
    std::vector<double> occ = {1.0, 0.999, 0.001, 0.0};
    // Toy HF-like gradient (lower bands favored).
    std::vector<double> grad = {-2.0, -1.0, 0.5, 1.0};

    auto map_inf = [&](double tau) {
        std::vector<double> trial(occ.size());
        for (size_t i = 0; i < occ.size(); ++i)
            trial[i] = occ[i] - tau * grad[i];
        constraint.project(trial);
        double linf = 0.0;
        for (size_t i = 0; i < occ.size(); ++i)
            linf = std::max(linf, std::abs(occ[i] - trial[i]));
        return linf;
    };

    const double tau_large = 1e6;
    const double old_scaled_metric = map_inf(tau_large) / tau_large; // old PG stopping metric
    const double fixed_tau_metric = map_inf(1.0); // robust stop metric with τ=1

    // Large τ can artificially shrink ||n-P(n-τg)||/τ and cause false convergence.
    EXPECT_LT(old_scaled_metric, 1e-5);
    // A fixed τ stop metric still reports meaningful non-stationarity here.
    EXPECT_GT(fixed_tau_metric, 1e-4);
}

TEST_F(OccupationConstraintTest, active_set_identification)
{
    std::vector<double> wk = {1.0};
    int nbands = 4;
    double nel = 2.0;
    OccupationConstraint constraint(ConstraintMethod::ActiveSet, nel, wk, nbands);

    std::vector<double> occ = {0.0, 0.5, 0.8, 1.0};
    std::vector<double> grad = {1.0, -0.5, 0.2, -0.3};

    auto info = constraint.identify_active_set(occ, grad);

    EXPECT_TRUE(info.at_lower[0]);    // n=0, grad>0 -> active at lower
    EXPECT_FALSE(info.at_lower[1]);   // n=0.5, free
    EXPECT_FALSE(info.at_upper[2]);   // n=0.8, free
    EXPECT_TRUE(info.at_upper[3]);    // n=1, grad<0 -> active at upper
    EXPECT_TRUE(info.is_free[1]);
    EXPECT_TRUE(info.is_free[2]);
}

class SigmaShiftOccParamTest : public ::testing::Test {};

TEST_F(SigmaShiftOccParamTest, compute_lambda_weighted_occupation_sum)
{
    const std::vector<double> wk = {0.25, 0.25, 0.25, 0.25};
    const int nbands = 3;
    const double nel = 2.0;
    SigmaShiftOccParam param(nel, wk, nbands);

    std::vector<double> z(12);
    for (int i = 0; i < 12; ++i)
    {
        z[static_cast<size_t>(i)] = 0.07 * static_cast<double>(i - 6);
    }

    const double lambda = param.compute_lambda(z);
    std::vector<double> occ;
    param.params_to_occ(z, lambda, occ);

    double sum = 0.0;
    for (int ik = 0; ik < 4; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            sum += wk[static_cast<size_t>(ik)] * occ[static_cast<size_t>(ik * nbands + ib)];
        }
    }
    EXPECT_NEAR(sum, nel, 1e-11) << "lambda=" << lambda;
}

TEST_F(SigmaShiftOccParamTest, transform_gradient_matches_finite_difference_linear_energy)
{
    const std::vector<double> wk = {0.5, 0.5};
    const int nbands = 2;
    const double nel = 1.2;
    SigmaShiftOccParam param(nel, wk, nbands);

    std::vector<double> z = {-0.31, 0.12, 0.37, -0.21};
    const double lambda = param.compute_lambda(z);
    const std::vector<double> h = {0.73, -0.41, 1.07, 0.19};

    std::vector<double> analytic;
    param.transform_gradient_batch(h, z, lambda, analytic);
    ASSERT_EQ(analytic.size(), z.size());

    const double eps = 1e-7;
    for (size_t j = 0; j < z.size(); ++j)
    {
        std::vector<double> z_plus = z;
        std::vector<double> z_minus = z;
        z_plus[j] += eps;
        z_minus[j] -= eps;

        SigmaShiftOccParam p_plus(nel, wk, nbands);
        SigmaShiftOccParam p_minus(nel, wk, nbands);
        const double lam_p = p_plus.compute_lambda(z_plus);
        const double lam_m = p_minus.compute_lambda(z_minus);
        std::vector<double> occ_p;
        std::vector<double> occ_m;
        p_plus.params_to_occ(z_plus, lam_p, occ_p);
        p_minus.params_to_occ(z_minus, lam_m, occ_m);

        double Ep = 0.0;
        double Em = 0.0;
        for (size_t i = 0; i < h.size(); ++i)
        {
            Ep += h[i] * occ_p[i];
            Em += h[i] * occ_m[i];
        }
        const double fd = (Ep - Em) / (2.0 * eps);
        EXPECT_NEAR(analytic[j], fd, 5e-6) << "j=" << j;
    }
}

TEST_F(SigmaShiftOccParamTest, stable_sigmoid_extremes)
{
    EXPECT_NEAR(SigmaShiftOccParam::stable_sigmoid(-80.0), 0.0, 1e-30);
    EXPECT_NEAR(SigmaShiftOccParam::stable_sigmoid(80.0), 1.0, 1e-30);
    EXPECT_NEAR(SigmaShiftOccParam::stable_sigmoid(0.0), 0.5, 1e-15);
}
