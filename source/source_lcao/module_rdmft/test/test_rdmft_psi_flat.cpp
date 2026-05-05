#include "gtest/gtest.h"
#include "source_lcao/module_rdmft/rdmft_psi_flat.h"
#include "source_psi/psi.h"

TEST(RdmftPsiFlat, RealRoundTrip)
{
    psi::Psi<double> p;
    p.resize(2, 3, 4);
    double* pp = p.get_pointer();
    for (int i = 0; i < 2 * 3 * 4; ++i)
    {
        pp[i] = 0.01 * static_cast<double>(i);
    }
    std::vector<double> flat;
    rdmft::psi_to_flat(p, flat);
    ASSERT_EQ(flat.size(), static_cast<size_t>(24));

    psi::Psi<double> q;
    q.resize(2, 3, 4);
    rdmft::flat_to_psi(flat, q);
    const double* qp = q.get_pointer();
    for (int i = 0; i < 24; ++i)
    {
        EXPECT_DOUBLE_EQ(pp[i], qp[i]);
    }
}

TEST(RdmftPsiFlat, ComplexRoundTrip)
{
    psi::Psi<std::complex<double>> p;
    p.resize(2, 2, 3);
    std::complex<double>* pp = p.get_pointer();
    for (int i = 0; i < 2 * 2 * 3; ++i)
    {
        pp[i] = std::complex<double>(0.1 * i, -0.2 * i);
    }
    std::vector<double> flat;
    rdmft::psi_to_flat(p, flat);
    ASSERT_EQ(flat.size(), static_cast<size_t>(2 * 12));

    psi::Psi<std::complex<double>> q;
    q.resize(2, 2, 3);
    rdmft::flat_to_psi(flat, q);
    const std::complex<double>* qp = q.get_pointer();
    for (int i = 0; i < 12; ++i)
    {
        EXPECT_DOUBLE_EQ(pp[i].real(), qp[i].real());
        EXPECT_DOUBLE_EQ(pp[i].imag(), qp[i].imag());
    }
}
