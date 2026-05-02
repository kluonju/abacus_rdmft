#ifndef RDMFT_TEST_STIEFEL_HELPER_H
#define RDMFT_TEST_STIEFEL_HELPER_H

// Test-only serial reference implementation of Stiefel-manifold operations
// used by the RDMFT unit tests. The production solver does NOT use this
// file; it provides a minimal, standalone (non-MPI, non-LAPACK) reference
// for the projected-tangent / retract / vector-transport algebra so tests
// can exercise EuclideanOptimizer + Riemannian wrapper logic without
// pulling in the full LCAO infrastructure.

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

extern "C" {
    void dgeqrf_(const int* m, const int* n, double* A, const int* lda,
                 double* tau, double* work, const int* lwork, int* info);
    void zgeqrf_(const int* m, const int* n, std::complex<double>* A, const int* lda,
                 std::complex<double>* tau, std::complex<double>* work,
                 const int* lwork, int* info);
    void dorgqr_(const int* m, const int* n, const int* k, double* A,
                 const int* lda, const double* tau, double* work,
                 const int* lwork, int* info);
    void zungqr_(const int* m, const int* n, const int* k, std::complex<double>* A,
                 const int* lda, const std::complex<double>* tau,
                 std::complex<double>* work, const int* lwork, int* info);
    void dgesv_(const int* n, const int* nrhs, double* A, const int* lda,
                int* ipiv, double* B, const int* ldb, int* info);
    void zgesv_(const int* n, const int* nrhs, std::complex<double>* A,
                const int* lda, int* ipiv, std::complex<double>* B,
                const int* ldb, int* info);
    void dgemm_(const char* ta, const char* tb, const int* m, const int* n,
                const int* k, const double* alpha, const double* A, const int* lda,
                const double* B, const int* ldb, const double* beta,
                double* C, const int* ldc);
    void zgemm_(const char* ta, const char* tb, const int* m, const int* n,
                const int* k, const std::complex<double>* alpha,
                const std::complex<double>* A, const int* lda,
                const std::complex<double>* B, const int* ldb,
                const std::complex<double>* beta,
                std::complex<double>* C, const int* ldc);
}

namespace rdmft
{

template <typename TK>
class StiefelManifold
{
  public:
    StiefelManifold() = default;

    StiefelManifold(int nbasis, int norbs, const TK* S_data = nullptr)
        : nbasis_(nbasis), norbs_(norbs)
    {
        if (S_data)
        {
            S_.assign(S_data, S_data + nbasis * nbasis);
            has_S_ = true;
        }
    }

    int nbasis() const { return nbasis_; }
    int norbs() const { return norbs_; }

    void project_tangent(const TK* C, const TK* G, TK* proj, int nbasis, int norbs) const
    {
        std::vector<TK> SC_or_C;
        const TK* effective_C = C;

        if (has_S_)
        {
            SC_or_C.resize(nbasis * norbs);
            mat_mul(S_.data(), C, SC_or_C.data(), nbasis, nbasis, norbs);
            effective_C = SC_or_C.data();
        }

        std::vector<TK> A(norbs * norbs, TK(0));
        mat_mul_herm(effective_C, G, A.data(), nbasis, norbs, norbs);

        for (int i = 0; i < norbs; ++i)
            for (int j = 0; j <= i; ++j)
            {
                TK sym_val = TK(0.5) * (A[i * norbs + j] + conj_val(A[j * norbs + i]));
                A[i * norbs + j] = sym_val;
                A[j * norbs + i] = conj_val(sym_val);
            }

        for (int i = 0; i < nbasis * norbs; ++i)
            proj[i] = G[i];

        mat_mul_sub(C, A.data(), proj, nbasis, norbs, norbs);
    }

    void retract(const TK* C, const TK* eta, double step, TK* C_new,
                 int nbasis, int norbs) const
    {
        std::vector<TK> Y(nbasis * norbs);
        for (int i = 0; i < nbasis * norbs; ++i)
            Y[i] = C[i] + step * eta[i];

        reorthogonalize(Y.data(), C_new, nbasis, norbs);
    }

    /// Householder QR retraction with sign-fixed R diagonal.
    /// Y = C + step * eta; Y = Q R; C_new = Q * diag(phase(R[i,i])).
    /// Mirrors EnergyGradient::retract_qr_serial. S must be identity (the
    /// helper does not expose a generalised Stiefel QR; callers should use
    /// X-space variables).
    void retract_qr(const TK* C, const TK* eta, double step, TK* C_new,
                    int nbasis, int norbs) const;

    /// Wen-Yin low-rank Cayley retraction (Math. Prog. 142 (2013) 397, Alg. 1).
    /// G interpreted as the Riemannian gradient; C_new = X(step) where
    /// X(t) = (I + (t/2) W)^{-1} (I - (t/2) W) C, W = G C^H - C G^H.
    /// Implemented via the SMW reduction to a 2p x 2p system.
    /// Mirrors EnergyGradient::retract_cayley_serial. S must be identity.
    void retract_cayley(const TK* C, const TK* G, double step, TK* C_new,
                         int nbasis, int norbs) const;

    void reorthogonalize(const TK* Y, TK* C_new, int nbasis, int norbs) const
    {
        std::vector<TK> SY(nbasis * norbs);
        if (has_S_)
            mat_mul(S_.data(), Y, SY.data(), nbasis, nbasis, norbs);
        else
            std::copy(Y, Y + nbasis * norbs, SY.begin());

        std::vector<TK> M(norbs * norbs, TK(0));
        mat_mul_herm(Y, SY.data(), M.data(), nbasis, norbs, norbs);

        std::vector<TK> L = M;
        cholesky_lower(L.data(), norbs);

        std::vector<TK> Linv_H(norbs * norbs, TK(0));
        invert_lower_conj_transpose(L.data(), Linv_H.data(), norbs);

        mat_mul(Y, Linv_H.data(), C_new, nbasis, norbs, norbs);
    }

    void vector_transport(const TK* C_new, const TK* xi, TK* transported,
                          int nbasis, int norbs) const
    {
        project_tangent(C_new, xi, transported, nbasis, norbs);
    }

    double inner_product(const TK* eta1, const TK* eta2, int nbasis, int norbs) const
    {
        double result = 0.0;
        if (has_S_)
        {
            std::vector<TK> S_eta2(nbasis * norbs);
            mat_mul(S_.data(), eta2, S_eta2.data(), nbasis, nbasis, norbs);
            for (int i = 0; i < nbasis * norbs; ++i)
                result += real_part(conj_val(eta1[i]) * S_eta2[i]);
        }
        else
        {
            for (int i = 0; i < nbasis * norbs; ++i)
                result += real_part(conj_val(eta1[i]) * eta2[i]);
        }
        return result;
    }

    double norm(const TK* eta, int nbasis, int norbs) const
    {
        return std::sqrt(inner_product(eta, eta, nbasis, norbs));
    }

  private:
    int nbasis_ = 0;
    int norbs_ = 0;
    bool has_S_ = false;
    std::vector<TK> S_;

    static TK conj_val(TK v);
    static double real_part(TK v);

    static void mat_mul(const TK* A, const TK* B, TK* C, int m, int k, int n);
    static void mat_mul_herm(const TK* A, const TK* B, TK* C, int m, int k, int n);
    static void mat_mul_sub(const TK* A, const TK* B, TK* C, int m, int k, int n);
    static void cholesky_lower(TK* A, int n);
    static void invert_lower_conj_transpose(const TK* L, TK* Linv_H, int n);
};

template <>
inline double StiefelManifold<double>::conj_val(double v) { return v; }

template <>
inline double StiefelManifold<double>::real_part(double v) { return v; }

template <>
inline void StiefelManifold<double>::mat_mul(const double* A, const double* B,
                                              double* C, int m, int k, int n)
{
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
        {
            double sum = 0.0;
            for (int l = 0; l < k; ++l)
                sum += A[i + l * m] * B[l + j * k];
            C[i + j * m] = sum;
        }
}

template <>
inline void StiefelManifold<double>::mat_mul_herm(const double* A, const double* B,
                                                   double* C, int m, int k, int n)
{
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < n; ++j)
        {
            double sum = 0.0;
            for (int l = 0; l < m; ++l)
                sum += A[l + i * m] * B[l + j * m];
            C[i + j * k] = sum;
        }
}

template <>
inline void StiefelManifold<double>::mat_mul_sub(const double* A, const double* B,
                                                  double* C, int m, int k, int n)
{
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
        {
            double sum = 0.0;
            for (int l = 0; l < k; ++l)
                sum += A[i + l * m] * B[l + j * k];
            C[i + j * m] -= sum;
        }
}

template <>
inline void StiefelManifold<double>::cholesky_lower(double* A, int n)
{
    for (int j = 0; j < n; ++j)
    {
        double sum = 0.0;
        for (int k = 0; k < j; ++k)
            sum += A[j + k * n] * A[j + k * n];
        A[j + j * n] = std::sqrt(std::max(1e-30, A[j + j * n] - sum));

        for (int i = j + 1; i < n; ++i)
        {
            sum = 0.0;
            for (int k = 0; k < j; ++k)
                sum += A[i + k * n] * A[j + k * n];
            A[i + j * n] = (A[i + j * n] - sum) / A[j + j * n];
        }
        for (int i = 0; i < j; ++i)
            A[i + j * n] = 0.0;
    }
}

template <>
inline void StiefelManifold<double>::invert_lower_conj_transpose(const double* L,
                                                                   double* Linv_H, int n)
{
    std::vector<double> Linv(n * n, 0.0);
    for (int j = 0; j < n; ++j)
    {
        Linv[j + j * n] = 1.0 / L[j + j * n];
        for (int i = j + 1; i < n; ++i)
        {
            double sum = 0.0;
            for (int k = j; k < i; ++k)
                sum += L[i + k * n] * Linv[k + j * n];
            Linv[i + j * n] = -sum / L[i + i * n];
        }
    }
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            Linv_H[i + j * n] = Linv[j + i * n];
}

template <>
inline std::complex<double> StiefelManifold<std::complex<double>>::conj_val(std::complex<double> v)
{
    return std::conj(v);
}

template <>
inline double StiefelManifold<std::complex<double>>::real_part(std::complex<double> v)
{
    return v.real();
}

template <>
inline void StiefelManifold<std::complex<double>>::mat_mul(
    const std::complex<double>* A, const std::complex<double>* B,
    std::complex<double>* C, int m, int k, int n)
{
    using cd = std::complex<double>;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
        {
            cd sum(0.0, 0.0);
            for (int l = 0; l < k; ++l)
                sum += A[i + l * m] * B[l + j * k];
            C[i + j * m] = sum;
        }
}

template <>
inline void StiefelManifold<std::complex<double>>::mat_mul_herm(
    const std::complex<double>* A, const std::complex<double>* B,
    std::complex<double>* C, int m, int k, int n)
{
    using cd = std::complex<double>;
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < n; ++j)
        {
            cd sum(0.0, 0.0);
            for (int l = 0; l < m; ++l)
                sum += std::conj(A[l + i * m]) * B[l + j * m];
            C[i + j * k] = sum;
        }
}

template <>
inline void StiefelManifold<std::complex<double>>::mat_mul_sub(
    const std::complex<double>* A, const std::complex<double>* B,
    std::complex<double>* C, int m, int k, int n)
{
    using cd = std::complex<double>;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
        {
            cd sum(0.0, 0.0);
            for (int l = 0; l < k; ++l)
                sum += A[i + l * m] * B[l + j * k];
            C[i + j * m] -= sum;
        }
}

template <>
inline void StiefelManifold<std::complex<double>>::cholesky_lower(std::complex<double>* A, int n)
{
    using cd = std::complex<double>;
    for (int j = 0; j < n; ++j)
    {
        double sum = 0.0;
        for (int k = 0; k < j; ++k)
            sum += std::norm(A[j + k * n]);
        double diag = std::sqrt(std::max(1e-30, A[j + j * n].real() - sum));
        A[j + j * n] = cd(diag, 0.0);

        for (int i = j + 1; i < n; ++i)
        {
            cd s(0.0, 0.0);
            for (int k = 0; k < j; ++k)
                s += A[i + k * n] * std::conj(A[j + k * n]);
            A[i + j * n] = (A[i + j * n] - s) / diag;
        }
        for (int i = 0; i < j; ++i)
            A[i + j * n] = cd(0.0, 0.0);
    }
}

template <>
inline void StiefelManifold<std::complex<double>>::invert_lower_conj_transpose(
    const std::complex<double>* L, std::complex<double>* Linv_H, int n)
{
    using cd = std::complex<double>;
    std::vector<cd> Linv(n * n, cd(0.0, 0.0));
    for (int j = 0; j < n; ++j)
    {
        Linv[j + j * n] = cd(1.0, 0.0) / L[j + j * n];
        for (int i = j + 1; i < n; ++i)
        {
            cd sum(0.0, 0.0);
            for (int k = j; k < i; ++k)
                sum += L[i + k * n] * Linv[k + j * n];
            Linv[i + j * n] = -sum / L[i + i * n];
        }
    }
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            Linv_H[i + j * n] = std::conj(Linv[j + i * n]);
}

// --- QR retraction: serial helper (Householder via LAPACK) ---

namespace test_helper_detail
{
inline double qr_diag_phase(double r) { return r >= 0.0 ? 1.0 : -1.0; }
inline std::complex<double> qr_diag_phase(const std::complex<double>& r)
{
    const double m = std::abs(r);
    if (m < 1.0e-300) return std::complex<double>(1.0, 0.0);
    return r / m;
}

inline int geqrf(int m, int n, double* A, int lda, double* tau,
                 double* work, int lwork)
{
    int info = 0;
    dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);
    return info;
}
inline int geqrf(int m, int n, std::complex<double>* A, int lda,
                 std::complex<double>* tau, std::complex<double>* work, int lwork)
{
    int info = 0;
    zgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);
    return info;
}
inline int orgqr(int m, int n, int k, double* A, int lda, const double* tau,
                 double* work, int lwork)
{
    int info = 0;
    dorgqr_(&m, &n, &k, A, &lda, tau, work, &lwork, &info);
    return info;
}
inline int orgqr(int m, int n, int k, std::complex<double>* A, int lda,
                 const std::complex<double>* tau, std::complex<double>* work, int lwork)
{
    int info = 0;
    zungqr_(&m, &n, &k, A, &lda, tau, work, &lwork, &info);
    return info;
}
inline int gesv(int n, int nrhs, double* A, int lda, int* ipiv,
                double* B, int ldb)
{
    int info = 0;
    dgesv_(&n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
    return info;
}
inline int gesv(int n, int nrhs, std::complex<double>* A, int lda, int* ipiv,
                std::complex<double>* B, int ldb)
{
    int info = 0;
    zgesv_(&n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
    return info;
}
inline double real_of(double v) { return v; }
inline double real_of(const std::complex<double>& v) { return v.real(); }
inline void gemm(char ta, char tb, int m, int n, int k,
                 double alpha, const double* A, int lda,
                 const double* B, int ldb,
                 double beta, double* C, int ldc)
{
    dgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}
inline void gemm(char ta, char tb, int m, int n, int k,
                 std::complex<double> alpha, const std::complex<double>* A, int lda,
                 const std::complex<double>* B, int ldb,
                 std::complex<double> beta, std::complex<double>* C, int ldc)
{
    zgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}
inline char trans_char(double) { return 'T'; }
inline char trans_char(std::complex<double>) { return 'C'; }
} // namespace test_helper_detail

template <typename TK>
inline void StiefelManifold<TK>::retract_qr(const TK* C, const TK* eta,
                                             double step, TK* C_new,
                                             int nbasis, int norbs) const
{
    using namespace test_helper_detail;
    // Y = C + step * eta (column-major, nbasis x norbs).
    std::vector<TK> Y(static_cast<size_t>(nbasis) * norbs, TK(0));
    for (int i = 0; i < nbasis * norbs; ++i)
        Y[i] = C[i] + TK(step) * eta[i];

    std::vector<TK> tau(norbs, TK(0));
    int info = 0;

    // Workspace size queries.
    TK qq = TK(0);
    info = geqrf(nbasis, norbs, nullptr, nbasis, tau.data(), &qq, -1);
    int lwork = std::max(static_cast<int>(real_of(qq)), norbs);
    TK qo = TK(0);
    info = orgqr(nbasis, norbs, norbs, nullptr, nbasis, tau.data(), &qo, -1);
    lwork = std::max(lwork, static_cast<int>(real_of(qo)));
    if (lwork < norbs) lwork = norbs;
    std::vector<TK> work(static_cast<size_t>(std::max(lwork, 1)), TK(0));

    info = geqrf(nbasis, norbs, Y.data(), nbasis, tau.data(),
                 work.data(), lwork);
    if (info != 0)
    {
        // Fall back: copy as-is (test will then fail orthogonality).
        std::copy(Y.begin(), Y.end(), C_new);
        return;
    }
    std::vector<TK> r_diag(norbs, TK(0));
    for (int j = 0; j < norbs; ++j) r_diag[j] = Y[j + j * nbasis];

    info = orgqr(nbasis, norbs, norbs, Y.data(), nbasis, tau.data(),
                 work.data(), lwork);
    if (info != 0)
    {
        std::copy(Y.begin(), Y.end(), C_new);
        return;
    }
    for (int j = 0; j < norbs; ++j)
    {
        const TK phase = qr_diag_phase(r_diag[j]);
        if (phase == TK(1.0)) continue;
        for (int i = 0; i < nbasis; ++i)
            Y[i + j * nbasis] *= phase;
    }
    std::copy(Y.begin(), Y.end(), C_new);
}

template <typename TK>
inline void StiefelManifold<TK>::retract_cayley(const TK* C, const TK* G,
                                                 double step, TK* C_new,
                                                 int nbasis, int norbs) const
{
    using namespace test_helper_detail;
    const int p = norbs;
    const int two_p = 2 * p;
    const TK one = TK(1.0);
    const TK zero = TK(0.0);
    const char tc = trans_char(TK());
    const double half_tau = 0.5 * step;

    std::vector<TK> XtX(p * p, TK(0));
    std::vector<TK> XtG(p * p, TK(0));
    std::vector<TK> GtX(p * p, TK(0));
    std::vector<TK> GtG(p * p, TK(0));
    gemm(tc, 'N', p, p, nbasis, one, C, nbasis, C, nbasis, zero, XtX.data(), p);
    gemm(tc, 'N', p, p, nbasis, one, C, nbasis, G, nbasis, zero, XtG.data(), p);
    gemm(tc, 'N', p, p, nbasis, one, G, nbasis, C, nbasis, zero, GtX.data(), p);
    gemm(tc, 'N', p, p, nbasis, one, G, nbasis, G, nbasis, zero, GtG.data(), p);

    std::vector<TK> Tmat(two_p * two_p, TK(0));
    auto Tref = [&](int gi, int gj) -> TK& {
        return Tmat[gi + gj * two_p];
    };
    for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i)
        {
            Tref(i, j)         = TK(half_tau) * XtG[i + j * p];
            Tref(i, p + j)     = TK(half_tau) * XtX[i + j * p];
            Tref(p + i, j)     = TK(-half_tau) * GtG[i + j * p];
            Tref(p + i, p + j) = TK(-half_tau) * GtX[i + j * p];
        }
    for (int d = 0; d < two_p; ++d) Tref(d, d) += TK(1.0);

    std::vector<TK> RHS(two_p * p, TK(0));
    for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i)
        {
            RHS[i + j * two_p]     = XtX[i + j * p];
            RHS[p + i + j * two_p] = -GtX[i + j * p];
        }
    std::vector<int> ipiv(two_p, 0);
    const int info = gesv(two_p, p, Tmat.data(), two_p,
                          ipiv.data(), RHS.data(), two_p);
    if (info != 0)
    {
        // Fall back: copy C unchanged.
        for (int i = 0; i < nbasis * p; ++i) C_new[i] = C[i];
        return;
    }

    std::vector<TK> Mtop(p * p, TK(0));
    std::vector<TK> Mbot(p * p, TK(0));
    for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i)
        {
            Mtop[i + j * p] = RHS[i + j * two_p];
            Mbot[i + j * p] = RHS[p + i + j * two_p];
        }
    std::vector<TK> Xnew(static_cast<size_t>(nbasis) * p, TK(0));
    gemm('N', 'N', nbasis, p, p, one, C, nbasis, Mbot.data(), p,
         zero, Xnew.data(), nbasis);
    gemm('N', 'N', nbasis, p, p, one, G, nbasis, Mtop.data(), p,
         one, Xnew.data(), nbasis);
    for (int i = 0; i < nbasis * p; ++i)
        C_new[i] = C[i] - TK(step) * Xnew[i];
}

} // namespace rdmft

#endif // RDMFT_TEST_STIEFEL_HELPER_H
