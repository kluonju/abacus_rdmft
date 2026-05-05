#ifndef RDMFT_PSI_FLAT_H
#define RDMFT_PSI_FLAT_H

#include "source_psi/psi.h"
#include <complex>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace rdmft
{

/// Flatten psi::Psi for binary I/O / Euclidean optimisers: real → N doubles;
/// complex → 2N (Re, Im) per element.
inline void psi_to_flat(const psi::Psi<double>& P, std::vector<double>& flat)
{
    const int nk = P.get_nk();
    const int nb = P.get_nbands();
    const int nbs = P.get_nbasis();
    const int n = nk * nb * nbs;
    flat.resize(static_cast<size_t>(n));
    if (n == 0)
    {
        return;
    }
    const double* p = &P(0, 0, 0);
    std::memcpy(flat.data(), p, static_cast<size_t>(n) * sizeof(double));
}

inline void psi_to_flat(const psi::Psi<std::complex<double>>& P, std::vector<double>& flat)
{
    const int nk = P.get_nk();
    const int nb = P.get_nbands();
    const int nbs = P.get_nbasis();
    const int n = nk * nb * nbs;
    flat.resize(static_cast<size_t>(2 * n));
    if (n == 0)
    {
        return;
    }
    const std::complex<double>* p = &P(0, 0, 0);
    for (int i = 0; i < n; ++i)
    {
        flat[static_cast<size_t>(2 * i)] = p[i].real();
        flat[static_cast<size_t>(2 * i + 1)] = p[i].imag();
    }
}

inline void flat_to_psi(const std::vector<double>& flat, psi::Psi<double>& P)
{
    const int nk = P.get_nk();
    const int nb = P.get_nbands();
    const int nbs = P.get_nbasis();
    const int n = nk * nb * nbs;
    if (n == 0)
    {
        return;
    }
    if (flat.size() != static_cast<size_t>(n))
    {
        throw std::runtime_error("flat_to_psi(double): flat size mismatch");
    }
    double* p = &P(0, 0, 0);
    std::memcpy(p, flat.data(), static_cast<size_t>(n) * sizeof(double));
}

inline void flat_to_psi(const std::vector<double>& flat, psi::Psi<std::complex<double>>& P)
{
    const int nk = P.get_nk();
    const int nb = P.get_nbands();
    const int nbs = P.get_nbasis();
    const int n = nk * nb * nbs;
    if (n == 0)
    {
        return;
    }
    if (flat.size() != static_cast<size_t>(2 * n))
    {
        throw std::runtime_error("flat_to_psi(complex): flat size mismatch");
    }
    std::complex<double>* p = &P(0, 0, 0);
    for (int i = 0; i < n; ++i)
    {
        p[i] = std::complex<double>(flat[static_cast<size_t>(2 * i)], flat[static_cast<size_t>(2 * i + 1)]);
    }
}

} // namespace rdmft

#endif
