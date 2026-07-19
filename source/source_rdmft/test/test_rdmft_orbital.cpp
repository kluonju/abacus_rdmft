#include "../rdmft_backend.h"
#include "../rdmft_occ_constraints.h"
#include "../rdmft_orbital_optimizer.h"
#include "../rdmft_params.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "gtest/gtest.h"

using namespace rdmft;

namespace
{
//! Toy natural-orbital energy E = sum_j n_j <c_j|A|c_j> with a fixed real
//! symmetric A and orthonormal columns C (overlap S = I, gamma-only style).
//! The minimiser sends the columns to eigenvectors of A, pairing the largest
//! occupations with the smallest eigenvalues (a trace-minimisation problem).
class RayleighBackend : public RdmftBackend
{
  public:
    RayleighBackend(int nbasis, int nbands, std::vector<double> A, std::vector<double> occ)
        : nbasis_(nbasis), nbands_(nbands), A_(std::move(A)), occ_(std::move(occ))
    {
        con_.nbnd = nbands;
        con_.nks = 1;
        con_.wk = {1.0};
        con_.isk = {1};
        con_.n_target = 0.0;
        // Start from the first nbands unit vectors.
        C_.assign(nbasis_ * nbands_, 0.0);
        for (int j = 0; j < nbands_; ++j)
        {
            C_[j + j * nbasis_] = 1.0;
        }
        orthonormalize(C_);
        saved_ = C_;
    }

    const OccConstraints& occ_constraints() const override { return con_; }

    double total_energy(const std::vector<double>& occ) override
    {
        (void)occ;
        double e = 0.0;
        std::vector<double> Ac(nbasis_);
        for (int j = 0; j < nbands_; ++j)
        {
            matvec(&C_[j * nbasis_], Ac.data());
            double cac = 0.0;
            for (int i = 0; i < nbasis_; ++i)
            {
                cac += C_[i + j * nbasis_] * Ac[i];
            }
            e += occ_[j] * cac;
        }
        return e;
    }

    void grad_occ(const std::vector<double>& occ, std::vector<double>& grad) override
    {
        (void)occ;
        grad.assign(con_.size(), 0.0);
    }

    bool has_orbital_optimization() const override { return true; }
    int orb_dim() const override { return nbasis_ * nbands_; }

    double riemannian_gradient(const std::vector<double>& occ, std::vector<double>& gR) override
    {
        (void)occ;
        // Ambient gradient columns g_j = 2 n_j A c_j.
        std::vector<double> G(nbasis_ * nbands_, 0.0);
        std::vector<double> Ac(nbasis_);
        for (int j = 0; j < nbands_; ++j)
        {
            matvec(&C_[j * nbasis_], Ac.data());
            for (int i = 0; i < nbasis_; ++i)
            {
                G[i + j * nbasis_] = 2.0 * occ_[j] * Ac[i];
            }
        }
        gR = G;
        project_tangent(C_, gR);
        double s = 0.0;
        for (double v : gR)
        {
            s += v * v;
        }
        return s;
    }

    double orb_inner(const std::vector<double>& a, const std::vector<double>& b) override
    {
        double s = 0.0;
        for (size_t i = 0; i < a.size(); ++i)
        {
            s += a[i] * b[i];
        }
        return s;
    }

    void orb_project_tangent(std::vector<double>& v) override { project_tangent(C_, v); }

    void orb_retract(const std::vector<double>& dir, double alpha) override
    {
        for (int i = 0; i < nbasis_ * nbands_; ++i)
        {
            C_[i] += alpha * dir[i];
        }
        orthonormalize(C_);
    }

    void orb_save() override { saved_ = C_; }
    void orb_restore() override { C_ = saved_; }

    double gram_residual() const
    {
        double r = 0.0;
        for (int a = 0; a < nbands_; ++a)
        {
            for (int b = 0; b < nbands_; ++b)
            {
                double s = 0.0;
                for (int i = 0; i < nbasis_; ++i)
                {
                    s += C_[i + a * nbasis_] * C_[i + b * nbasis_];
                }
                r += std::fabs(s - (a == b ? 1.0 : 0.0));
            }
        }
        return r;
    }

  private:
    void matvec(const double* c, double* out) const
    {
        for (int i = 0; i < nbasis_; ++i)
        {
            double s = 0.0;
            for (int j = 0; j < nbasis_; ++j)
            {
                s += A_[i + j * nbasis_] * c[j];
            }
            out[i] = s;
        }
    }

    void project_tangent(const std::vector<double>& C, std::vector<double>& G) const
    {
        // sym = 0.5 (C^T G + G^T C); eta = G - C sym.
        std::vector<double> M(nbands_ * nbands_, 0.0);
        for (int a = 0; a < nbands_; ++a)
        {
            for (int b = 0; b < nbands_; ++b)
            {
                double s = 0.0;
                for (int i = 0; i < nbasis_; ++i)
                {
                    s += C[i + a * nbasis_] * G[i + b * nbasis_];
                }
                M[a + b * nbands_] = s;
            }
        }
        std::vector<double> Msym(nbands_ * nbands_, 0.0);
        for (int a = 0; a < nbands_; ++a)
        {
            for (int b = 0; b < nbands_; ++b)
            {
                Msym[a + b * nbands_] = 0.5 * (M[a + b * nbands_] + M[b + a * nbands_]);
            }
        }
        for (int b = 0; b < nbands_; ++b)
        {
            for (int i = 0; i < nbasis_; ++i)
            {
                double s = 0.0;
                for (int a = 0; a < nbands_; ++a)
                {
                    s += C[i + a * nbasis_] * Msym[a + b * nbands_];
                }
                G[i + b * nbasis_] -= s;
            }
        }
    }

    void orthonormalize(std::vector<double>& C) const
    {
        // Modified Gram-Schmidt on the columns.
        for (int j = 0; j < nbands_; ++j)
        {
            double* cj = &C[j * nbasis_];
            for (int k = 0; k < j; ++k)
            {
                const double* ck = &C[k * nbasis_];
                double proj = 0.0;
                for (int i = 0; i < nbasis_; ++i)
                {
                    proj += ck[i] * cj[i];
                }
                for (int i = 0; i < nbasis_; ++i)
                {
                    cj[i] -= proj * ck[i];
                }
            }
            double nrm = 0.0;
            for (int i = 0; i < nbasis_; ++i)
            {
                nrm += cj[i] * cj[i];
            }
            nrm = std::sqrt(nrm);
            if (nrm > 1.0e-14)
            {
                for (int i = 0; i < nbasis_; ++i)
                {
                    cj[i] /= nrm;
                }
            }
        }
    }

    int nbasis_;
    int nbands_;
    std::vector<double> A_;   // symmetric nbasis x nbasis, column-major
    std::vector<double> occ_; // fixed occupations, length nbands
    std::vector<double> C_;   // nbasis x nbands, column-major
    std::vector<double> saved_;
    OccConstraints con_;
};

//! Jacobi eigenvalue algorithm for a small symmetric matrix (column-major).
std::vector<double> jacobi_eigenvalues(std::vector<double> A, int n)
{
    for (int sweep = 0; sweep < 100; ++sweep)
    {
        double off = 0.0;
        for (int p = 0; p < n; ++p)
        {
            for (int q = p + 1; q < n; ++q)
            {
                off += A[p + q * n] * A[p + q * n];
            }
        }
        if (off < 1.0e-20)
        {
            break;
        }
        for (int p = 0; p < n; ++p)
        {
            for (int q = p + 1; q < n; ++q)
            {
                const double apq = A[p + q * n];
                if (std::fabs(apq) < 1.0e-18)
                {
                    continue;
                }
                const double app = A[p + p * n];
                const double aqq = A[q + q * n];
                const double phi = 0.5 * std::atan2(2.0 * apq, aqq - app);
                const double c = std::cos(phi);
                const double s = std::sin(phi);
                for (int i = 0; i < n; ++i)
                {
                    const double aip = A[i + p * n];
                    const double aiq = A[i + q * n];
                    A[i + p * n] = c * aip - s * aiq;
                    A[i + q * n] = s * aip + c * aiq;
                }
                for (int i = 0; i < n; ++i)
                {
                    const double api = A[p + i * n];
                    const double aqi = A[q + i * n];
                    A[p + i * n] = c * api - s * aqi;
                    A[q + i * n] = s * api + c * aqi;
                }
            }
        }
    }
    std::vector<double> ev(n);
    for (int i = 0; i < n; ++i)
    {
        ev[i] = A[i + i * n];
    }
    std::sort(ev.begin(), ev.end());
    return ev;
}

double reference_min_energy(const std::vector<double>& A, int n, std::vector<double> occ)
{
    std::vector<double> ev = jacobi_eigenvalues(A, n); // ascending
    std::sort(occ.rbegin(), occ.rend());               // descending
    double e = 0.0;
    for (size_t j = 0; j < occ.size(); ++j)
    {
        e += occ[j] * ev[j];
    }
    return e;
}

std::vector<double> make_symmetric(int n, unsigned seed)
{
    std::vector<double> A(n * n, 0.0);
    unsigned s = seed;
    auto rnd = [&]() {
        s = s * 1103515245u + 12345u;
        return ((s >> 16) & 0x7fff) / 32768.0 - 0.5;
    };
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            const double v = rnd();
            A[i + j * n] = v;
            A[j + i * n] = v;
        }
    }
    return A;
}
} // namespace

TEST(OrbitalOptimizer, ConjugateGradientReachesEigenspace)
{
    const int nbasis = 6;
    const int nbands = 3;
    std::vector<double> A = make_symmetric(nbasis, 7u);
    std::vector<double> occ = {1.0, 0.7, 0.4};
    const double e_ref = reference_min_energy(A, nbasis, occ);

    RayleighBackend backend(nbasis, nbands, A, occ);
    RdmftParams params;
    params.orb_optimizer = OrbOptimizerType::CG;
    params.orb_maxiter = 400;
    params.orb_grad_tol = 1.0e-8;

    std::vector<double> occ_vec = occ; // one k-point
    double etot = 0.0;
    OrbitalOptimizer opt;
    OrbBlockResult r = opt.run(backend, params, occ_vec, etot);

    EXPECT_LT(backend.gram_residual(), 1.0e-8);
    EXPECT_NEAR(etot, e_ref, 1.0e-5);
    (void)r;
}

TEST(OrbitalOptimizer, SteepestDescentDecreasesEnergy)
{
    const int nbasis = 5;
    const int nbands = 2;
    std::vector<double> A = make_symmetric(nbasis, 42u);
    std::vector<double> occ = {1.0, 0.5};
    const double e_ref = reference_min_energy(A, nbasis, occ);

    RayleighBackend backend(nbasis, nbands, A, occ);
    RdmftParams params;
    params.orb_optimizer = OrbOptimizerType::SD;
    params.orb_maxiter = 800;
    params.orb_grad_tol = 1.0e-9;

    std::vector<double> occ_vec = occ;
    double etot = 0.0;
    OrbitalOptimizer opt;
    opt.run(backend, params, occ_vec, etot);

    EXPECT_LT(backend.gram_residual(), 1.0e-8);
    EXPECT_NEAR(etot, e_ref, 1.0e-4);
}

TEST(OrbitalOptimizer, LBFGSReachesEigenspace)
{
    const int nbasis = 7;
    const int nbands = 3;
    std::vector<double> A = make_symmetric(nbasis, 123u);
    std::vector<double> occ = {1.0, 0.8, 0.3};
    const double e_ref = reference_min_energy(A, nbasis, occ);

    RayleighBackend backend(nbasis, nbands, A, occ);
    RdmftParams params;
    params.orb_optimizer = OrbOptimizerType::LBFGS;
    params.orb_maxiter = 400;
    params.orb_grad_tol = 1.0e-8;
    params.lbfgs_memory = 8;

    std::vector<double> occ_vec = occ;
    double etot = 0.0;
    OrbitalOptimizer opt;
    opt.run(backend, params, occ_vec, etot);

    EXPECT_LT(backend.gram_residual(), 1.0e-8);
    EXPECT_NEAR(etot, e_ref, 1.0e-5);
}
