#include "../rdmft_backend.h"
#include "../rdmft_occ_constraints.h"
#include "../rdmft_occ_optimizer.h"
#include "../rdmft_params.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "gtest/gtest.h"

using namespace rdmft;

namespace
{
//! Separable convex quadratic occupation energy:
//!   E = sum_ik wk_ik * 0.5 * k_i * (n_i - t_i)^2
//! whose physical gradient dE/dn = wk_ik * k_i * (n_i - t_i) exercises the
//! box + weighted-sum constrained optimisers with a known KKT solution.
class QuadraticOccBackend : public RdmftBackend
{
  public:
    QuadraticOccBackend(const OccConstraints& con, std::vector<double> t, std::vector<double> k)
        : con_(con), t_(std::move(t)), k_(std::move(k))
    {
    }

    const OccConstraints& occ_constraints() const override { return con_; }

    double total_energy(const std::vector<double>& occ) override
    {
        double e = 0.0;
        for (int ik = 0; ik < con_.nks; ++ik)
        {
            const double w = con_.wk[ik];
            for (int ib = 0; ib < con_.nbnd; ++ib)
            {
                const int i = ib + ik * con_.nbnd;
                const double d = occ[i] - t_[i];
                e += w * 0.5 * k_[i] * d * d;
            }
        }
        return e;
    }

    void grad_occ(const std::vector<double>& occ, std::vector<double>& grad) override
    {
        grad.assign(con_.size(), 0.0);
        for (int ik = 0; ik < con_.nks; ++ik)
        {
            const double w = con_.wk[ik];
            for (int ib = 0; ib < con_.nbnd; ++ib)
            {
                const int i = ib + ik * con_.nbnd;
                grad[i] = w * k_[i] * (occ[i] - t_[i]);
            }
        }
    }

  private:
    OccConstraints con_;
    std::vector<double> t_;
    std::vector<double> k_;
};

//! Independent KKT reference: n_i = clip(t_i + mu/k_i, 0, 1), mu from the
//! weighted-sum electron constraint (single Fermi level).
std::vector<double> kkt_reference(const OccConstraints& con, const std::vector<double>& t,
                                  const std::vector<double>& k, double target)
{
    auto weighted = [&](double mu) {
        double s = 0.0;
        for (int ik = 0; ik < con.nks; ++ik)
        {
            for (int ib = 0; ib < con.nbnd; ++ib)
            {
                const int i = ib + ik * con.nbnd;
                const double n = std::max(0.0, std::min(1.0, t[i] + mu / k[i]));
                s += con.wk[ik] * n;
            }
        }
        return s;
    };
    double lo = -1.0e3;
    double hi = 1.0e3;
    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (weighted(mid) < target)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }
    const double mu = 0.5 * (lo + hi);
    std::vector<double> n(con.size());
    for (int ik = 0; ik < con.nks; ++ik)
    {
        for (int ib = 0; ib < con.nbnd; ++ib)
        {
            const int i = ib + ik * con.nbnd;
            n[i] = std::max(0.0, std::min(1.0, t[i] + mu / k[i]));
        }
    }
    return n;
}

OccConstraints make_constraints(int nbnd, int nks, double ne)
{
    OccConstraints con;
    con.nbnd = nbnd;
    con.nks = nks;
    con.wk.assign(nks, 1.0 / nks);
    con.isk.assign(nks, 1);
    con.n_target = ne;
    return con;
}
} // namespace

TEST(OccConstraints, ProximalProjectSatisfiesConstraint)
{
    OccConstraints con = make_constraints(4, 2, 3.0);
    std::vector<double> occ = {0.9, 0.2, 1.5, -0.3, 0.4, 0.8, 0.1, 0.95};
    con.proximal_project(occ);
    // Box respected.
    for (double v : occ)
    {
        EXPECT_GE(v, -1.0e-12);
        EXPECT_LE(v, 1.0 + 1.0e-12);
    }
    // Weighted sum equals target.
    EXPECT_NEAR(con.weighted_sum(occ), 3.0, 1.0e-8);
}

TEST(OccConstraints, StoppingMapVanishesAtKKT)
{
    OccConstraints con = make_constraints(4, 2, 2.5);
    std::vector<double> t = {0.6, 0.4, 0.7, 0.3, 0.5, 0.55, 0.45, 0.5};
    std::vector<double> k = {1.0, 2.0, 0.5, 1.5, 1.0, 0.8, 1.2, 1.0};
    std::vector<double> nref = kkt_reference(con, t, k, 2.5);
    QuadraticOccBackend backend(con, t, k);
    std::vector<double> grad;
    backend.grad_occ(nref, grad);
    // At the KKT point the SPG2 projected-gradient map must vanish.
    EXPECT_LT(con.pg_map_grad_inf(nref, grad, 1.0), 1.0e-6);
}

TEST(OccConstraints, EbiErfRoundTrip)
{
    OccConstraints con = make_constraints(3, 1, 1.5);
    std::vector<double> occ = {0.2, 0.5, 0.8};
    std::vector<double> x;
    con.ebi_occ_to_params(occ, x);
    std::vector<double> back(3);
    for (int i = 0; i < 3; ++i)
    {
        back[i] = OccConstraints::ebi_occ_from_arg(x[i]);
    }
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(back[i], occ[i], 1.0e-10);
    }
}

TEST(OccConstraints, EbiSolveMuMatchesTarget)
{
    OccConstraints con = make_constraints(4, 2, 2.0);
    std::vector<double> x = {0.1, -0.2, 0.3, 0.0, -0.1, 0.2, 0.4, -0.3};
    std::vector<double> occ;
    con.ebi_sync_occ_from_x(x, occ);
    EXPECT_NEAR(con.weighted_sum(occ), 2.0, 1.0e-8);
}

TEST(OccOptimizer, SPG2ConvergesToKKT)
{
    OccConstraints con = make_constraints(5, 2, 3.0);
    std::vector<double> t = {0.7, 0.3, 0.9, 0.1, 0.5, 0.6, 0.4, 0.8, 0.2, 0.55};
    std::vector<double> k = {1.0, 1.5, 0.8, 2.0, 1.0, 1.2, 0.9, 1.1, 1.3, 1.0};
    std::vector<double> nref = kkt_reference(con, t, k, 3.0);

    QuadraticOccBackend backend(con, t, k);
    RdmftParams params;
    params.occ_optimizer = OccOptimizerType::SPG2;
    params.occ_maxiter = 200;
    params.occ_grad_tol = 1.0e-6;

    std::vector<double> occ(con.size(), 3.0 / con.weighted_sum(std::vector<double>(con.size(), 1.0)));
    double etot = 0.0;
    OccOptimizer opt;
    OccBlockResult r = opt.run(backend, params, occ, etot);

    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(con.weighted_sum(occ), 3.0, 1.0e-7);
    double max_diff = 0.0;
    for (int i = 0; i < con.size(); ++i)
    {
        max_diff = std::max(max_diff, std::fabs(occ[i] - nref[i]));
    }
    EXPECT_LT(max_diff, 1.0e-4);
    EXPECT_NEAR(etot, backend.total_energy(nref), 1.0e-8);
}

TEST(OccOptimizer, EBIConvergesToInteriorKKT)
{
    // Interior solution so the erf parameterisation can reach it exactly.
    OccConstraints con = make_constraints(4, 2, 2.4);
    std::vector<double> t = {0.55, 0.45, 0.6, 0.4, 0.5, 0.52, 0.48, 0.5};
    std::vector<double> k = {1.0, 1.2, 0.9, 1.1, 1.0, 0.95, 1.05, 1.0};
    std::vector<double> nref = kkt_reference(con, t, k, 2.4);
    // Confirm the reference is interior.
    for (double v : nref)
    {
        ASSERT_GT(v, 1.0e-3);
        ASSERT_LT(v, 1.0 - 1.0e-3);
    }

    QuadraticOccBackend backend(con, t, k);
    RdmftParams params;
    params.occ_optimizer = OccOptimizerType::EBI;
    params.occ_maxiter = 500;
    params.occ_grad_tol = 1.0e-9;
    params.occ_tol = 1.0e-10;

    std::vector<double> occ(con.size(), 0.3);
    double etot = 0.0;
    OccOptimizer opt;
    OccBlockResult r = opt.run(backend, params, occ, etot);

    EXPECT_NEAR(con.weighted_sum(occ), 2.4, 1.0e-6);
    double max_diff = 0.0;
    for (int i = 0; i < con.size(); ++i)
    {
        max_diff = std::max(max_diff, std::fabs(occ[i] - nref[i]));
    }
    EXPECT_LT(max_diff, 1.0e-3);
}

TEST(OccOptimizer, SPG2SpinResolvedConstraint)
{
    OccConstraints con;
    con.nbnd = 3;
    con.nks = 2;
    con.wk = {1.0, 1.0};
    con.isk = {1, 2}; // one up row, one down row
    con.fix_magnetization = true;
    con.n_target_up = 2.0;
    con.n_target_down = 1.0;
    con.n_target = 3.0;

    std::vector<double> t = {0.9, 0.6, 0.4, 0.5, 0.3, 0.2};
    std::vector<double> k(6, 1.0);
    QuadraticOccBackend backend(con, t, k);

    RdmftParams params;
    params.occ_optimizer = OccOptimizerType::SPG2;
    params.occ_maxiter = 300;
    params.occ_grad_tol = 1.0e-8;
    params.fix_magnetization = true;

    std::vector<double> occ(6, 0.5);
    double etot = 0.0;
    OccOptimizer opt;
    opt.run(backend, params, occ, etot);

    EXPECT_NEAR(con.weighted_sum_ispin(occ, 1), 2.0, 1.0e-6);
    EXPECT_NEAR(con.weighted_sum_ispin(occ, 2), 1.0, 1.0e-6);
}
