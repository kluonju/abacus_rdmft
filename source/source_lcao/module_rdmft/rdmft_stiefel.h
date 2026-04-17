#ifndef RDMFT_STIEFEL_H
#define RDMFT_STIEFEL_H

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

namespace rdmft
{

/// Operations on the Stiefel manifold St(p, n; S) where S is the overlap matrix.
/// C in R^{n x p} (or C^{n x p}) with C^H S C = I.
///
/// For planewave basis, S = I (identity).
/// For LCAO basis, S is the overlap matrix.
///
/// TK = double (gamma-only) or complex<double> (multi-k).
template <typename TK>
class StiefelManifold
{
  public:
    StiefelManifold() = default;

    /// Initialize with basis size and number of orbitals.
    /// S_data: overlap matrix in column-major, size nbasis*nbasis (nullptr for PW/identity).
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

    /// Project Euclidean gradient G onto the tangent space at C.
    /// Riemannian gradient = G - C * sym(C^H * G) for canonical metric on St(p,n).
    /// With overlap S: grad = S^{-1} G - C * sym(C^H G), but typically we work
    /// in the S-weighted space so: grad = G - S*C * sym(C^H * G).
    /// Here we use the simpler projection for the embedded metric:
    ///   proj_C(Z) = Z - C * sym(C^H S Z)
    /// which gives a tangent vector satisfying C^H S Z + Z^H S C = 0.
    void project_tangent(const TK* C, const TK* G, TK* proj, int nbasis, int norbs) const
    {
        // Compute A = C^H * G (if no S) or C^H * S * G
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

        // Symmetrize: A_sym = 0.5*(A + A^H)
        for (int i = 0; i < norbs; ++i)
            for (int j = 0; j <= i; ++j)
            {
                TK sym_val = TK(0.5) * (A[i * norbs + j] + conj_val(A[j * norbs + i]));
                A[i * norbs + j] = sym_val;
                A[j * norbs + i] = conj_val(sym_val);
            }

        // proj = G - C * A_sym
        for (int i = 0; i < nbasis * norbs; ++i)
            proj[i] = G[i];

        mat_mul_sub(C, A.data(), proj, nbasis, norbs, norbs);
    }

    /// Polar retraction: R_C(eta) = (C + eta) * [(C+eta)^H S (C+eta)]^{-1/2}
    /// For simplicity, we use QR-based retraction as a fallback.
    void retract(const TK* C, const TK* eta, double step, TK* C_new,
                 int nbasis, int norbs) const
    {
        std::vector<TK> Y(nbasis * norbs);
        for (int i = 0; i < nbasis * norbs; ++i)
            Y[i] = C[i] + step * eta[i];

        reorthogonalize(Y.data(), C_new, nbasis, norbs);
    }

    /// Reorthogonalize columns of Y w.r.t. overlap S using Cholesky:
    /// Y^H S Y = L L^H, then C_new = Y * L^{-H}
    void reorthogonalize(const TK* Y, TK* C_new, int nbasis, int norbs) const
    {
        // Compute M = Y^H S Y
        std::vector<TK> SY(nbasis * norbs);
        if (has_S_)
            mat_mul(S_.data(), Y, SY.data(), nbasis, nbasis, norbs);
        else
            std::copy(Y, Y + nbasis * norbs, SY.begin());

        std::vector<TK> M(norbs * norbs, TK(0));
        mat_mul_herm(Y, SY.data(), M.data(), nbasis, norbs, norbs);

        // Cholesky: M = L*L^H
        std::vector<TK> L = M;
        cholesky_lower(L.data(), norbs);

        // C_new = Y * L^{-H}  (solve L^H * X^T = Y^T)
        // Equivalently: C_new = Y * inv(L)^H
        std::vector<TK> Linv_H(norbs * norbs, TK(0));
        invert_lower_conj_transpose(L.data(), Linv_H.data(), norbs);

        // C_new = Y * Linv_H
        mat_mul(Y, Linv_H.data(), C_new, nbasis, norbs, norbs);
    }

    /// Vector transport by projection: T_eta(xi) = proj_{C_new}(xi)
    void vector_transport(const TK* C_new, const TK* xi, TK* transported,
                          int nbasis, int norbs) const
    {
        project_tangent(C_new, xi, transported, nbasis, norbs);
    }

    /// Inner product on the tangent space: <eta1, eta2> = Re Tr(eta1^H S eta2)
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

    /// Norm of a tangent vector
    double norm(const TK* eta, int nbasis, int norbs) const
    {
        return std::sqrt(inner_product(eta, eta, nbasis, norbs));
    }

  private:
    int nbasis_ = 0;
    int norbs_ = 0;
    bool has_S_ = false;
    std::vector<TK> S_;

    // Helper: conjugate for complex, identity for double
    static TK conj_val(TK v);
    static double real_part(TK v);

    // Simple dense matrix operations (for serial/small problems)
    // For production, these should use BLAS/LAPACK
    static void mat_mul(const TK* A, const TK* B, TK* C, int m, int k, int n);
    static void mat_mul_herm(const TK* A, const TK* B, TK* C, int m, int k, int n);
    static void mat_mul_sub(const TK* A, const TK* B, TK* C, int m, int k, int n);
    static void cholesky_lower(TK* A, int n);
    static void invert_lower_conj_transpose(const TK* L, TK* Linv_H, int n);
};

// --- Template specializations for double ---

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
                sum += A[i + l * m] * B[l + j * k]; // col-major
            C[i + j * m] = sum;
        }
}

template <>
inline void StiefelManifold<double>::mat_mul_herm(const double* A, const double* B,
                                                   double* C, int m, int k, int n)
{
    // C = A^T * B, A is m x k, B is m x n, C is k x n (col-major)
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
    // C -= A * B, A is m x k, B is k x n, C is m x n (col-major)
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
    // Compute inv(L), then transpose (for real, transpose = conjugate transpose)
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
    // Transpose
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            Linv_H[i + j * n] = Linv[j + i * n];
}

// --- Template specializations for complex<double> ---

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

} // namespace rdmft

#endif // RDMFT_STIEFEL_H
