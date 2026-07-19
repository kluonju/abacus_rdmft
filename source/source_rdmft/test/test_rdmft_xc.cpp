#include "../rdmft_xc.h"

#include <cmath>

#include "gtest/gtest.h"

using namespace rdmft_core;

namespace
{
// Central finite difference of the coupling derivative.
double fd_dg(const RdmftXC& xc, double n, double h = 1.0e-6)
{
    return (xc.g(n + h) - xc.g(n - h)) / (2.0 * h);
}
} // namespace

TEST(RdmftXC, HartreeFock)
{
    RdmftXC xc(XcType::HF, 1.0, 1.0e-8);
    EXPECT_DOUBLE_EQ(xc.g(0.3), 0.3);
    EXPECT_DOUBLE_EQ(xc.dg(0.7), 1.0);
    EXPECT_TRUE(xc.is_separable());
    EXPECT_EQ(xc.n_channels(), 1);
}

TEST(RdmftXC, MullerSqrt)
{
    RdmftXC xc(XcType::Muller, 0.5, 1.0e-8);
    EXPECT_NEAR(xc.g(0.25), 0.5, 1.0e-12);
    EXPECT_NEAR(xc.dg(0.25), fd_dg(xc, 0.25), 1.0e-5);
    EXPECT_TRUE(xc.is_separable());
}

TEST(RdmftXC, PowerDerivativeMatchesFD)
{
    RdmftXC xc(XcType::Power, 0.656, 1.0e-8);
    for (double n = 0.05; n < 1.0; n += 0.1)
    {
        EXPECT_NEAR(xc.dg(n), fd_dg(xc, n), 1.0e-4) << "n=" << n;
    }
}

TEST(RdmftXC, PowerRegularizationBounded)
{
    RdmftXC xc(XcType::Power, 0.5, 1.0e-8);
    // Derivative near zero is bounded by alpha*eps^(alpha-1).
    const double bound = 0.5 * std::pow(1.0e-8, 0.5 - 1.0);
    EXPECT_LE(xc.dg(0.0), bound + 1.0);
    EXPECT_GT(xc.g(0.0), 0.0); // g(0) = (1-alpha) eps^alpha > 0
}

TEST(RdmftXC, ChannelCounts)
{
    EXPECT_EQ(RdmftXC(XcType::CHF, 1.0, 1e-8).n_channels(), 2);
    EXPECT_EQ(RdmftXC(XcType::CGA, 1.0, 1e-8).n_channels(), 2);
    EXPECT_EQ(RdmftXC(XcType::GEO, 0.75, 1e-8).n_channels(), 3);
    EXPECT_EQ(RdmftXC(XcType::HybOpt, 0.541076, 1e-8).n_channels(), 2);
    EXPECT_EQ(RdmftXC(XcType::BOWMOD, 0.61, 1e-8).n_channels(), 4);
}

TEST(RdmftXC, ChannelDerivativeConsistency)
{
    // For every functional the channel dw must match a finite difference of w.
    struct Case
    {
        XcType type;
        double alpha;
    };
    const Case cases[] = {{XcType::CHF, 1.0},   {XcType::CGA, 1.0}, {XcType::GEO, 0.75},
                          {XcType::HybOpt, 0.541076}, {XcType::Muller, 0.5}};
    const double h = 1.0e-6;
    for (const Case& c : cases)
    {
        RdmftXC xc(c.type, c.alpha, 1.0e-10);
        for (int it = 1; it <= xc.n_channels(); ++it)
        {
            for (double n = 0.2; n < 0.9; n += 0.2)
            {
                const double w_plus = xc.channel(it, n + h).w;
                const double w_minus = xc.channel(it, n - h).w;
                const double fd = (w_plus - w_minus) / (2.0 * h);
                EXPECT_NEAR(xc.channel(it, n).dw, fd, 1.0e-3)
                    << "type=" << (int)c.type << " it=" << it << " n=" << n;
            }
        }
    }
}

TEST(RdmftXC, GUExtraDiagonal)
{
    RdmftXC xc(XcType::GU, 0.5, 1.0e-8);
    EXPECT_TRUE(xc.has_extra_diag());
    EXPECT_NEAR(xc.gu_diag_factor(0.5), 0.25 - 0.5, 1.0e-12);
    EXPECT_NEAR(xc.gu_diag_factor_deriv(0.5), 0.0, 1.0e-12);
}

TEST(RdmftXC, BinaryEntropy)
{
    EXPECT_NEAR(RdmftXC::binary_entropy(0.5), 2.0 * 0.5 * std::log(0.5), 1.0e-12);
    const double h = 1.0e-6;
    const double fd = (RdmftXC::binary_entropy(0.3 + h) - RdmftXC::binary_entropy(0.3 - h)) / (2 * h);
    EXPECT_NEAR(RdmftXC::binary_entropy_deriv(0.3), fd, 1.0e-4);
}
