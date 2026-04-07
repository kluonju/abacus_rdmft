#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_optimizer.h"
#include <cmath>
#include <vector>
#include <functional>

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
