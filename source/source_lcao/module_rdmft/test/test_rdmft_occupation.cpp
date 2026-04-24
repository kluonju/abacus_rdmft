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

TEST_F(OccupationParamTest, Logistic_batch_enforces_electron_number)
{
    OccupationParam param(OccParamType::Logistic);
    std::vector<double> wk = {0.4, 0.6};
    OccupationConstraint constraint(ConstraintMethod::AugmentedLagrangian, 1.3, wk, 2);

    std::vector<double> params = {-1.2, 0.1, 0.7, -0.4};
    std::vector<double> occ;
    param.params_to_occ_batch(params, constraint, occ);

    EXPECT_NEAR(constraint.constraint_violation(occ), 0.0, 1e-11);
    for (double n : occ)
    {
        EXPECT_GE(n, 0.0);
        EXPECT_LE(n, 1.0);
    }
}

TEST_F(OccupationParamTest, Logistic_batch_gradient_matches_finite_difference)
{
    OccupationParam param(OccParamType::Logistic);
    std::vector<double> wk = {0.4, 0.6};
    OccupationConstraint constraint(ConstraintMethod::AugmentedLagrangian, 1.1, wk, 2);

    const std::vector<double> params = {-0.8, 0.2, 0.5, -0.3};
    const std::vector<double> dE_dn = {0.7, -1.1, 0.4, 1.3};

    std::vector<double> analytic;
    param.transform_gradient_batch_solver(dE_dn, params, constraint, analytic);

    auto objective = [&](const std::vector<double>& p) {
        std::vector<double> occ;
        param.params_to_occ_batch(p, constraint, occ);
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
