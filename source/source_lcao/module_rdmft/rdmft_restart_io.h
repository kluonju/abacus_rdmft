#ifndef RDMFT_RESTART_IO_H
#define RDMFT_RESTART_IO_H

#include "rdmft_psi_flat.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_io/module_restart/restart.h"
#include "source_psi/psi.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_common.h"
#include "source_base/tool_quit.h"

#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace rdmft
{

/// Text: `Restart_rdmft.txt`; binary: `rdmft_occ_0_0` (rank 0), `rdmft_wfc_RANK_0` per rank.
template <typename TK>
void write_rdmft_restart(const Parallel_Orbitals& pv,
                         const int nk,
                         const int nbands,
                         const int nspin,
                         const bool gamma_only,
                         const double n_electrons,
                         const std::string& rdmft_functional,
                         const std::vector<double>& occ_flat,
                         const psi::Psi<TK>& wfc)
{
    if (GlobalC::restart.folder.empty())
    {
        ModuleBase::WARNING_QUIT("rdmft::write_rdmft_restart", "GlobalC::restart.folder is empty.");
    }
    const int nocc = nk * nbands;
    if (static_cast<int>(occ_flat.size()) != nocc)
    {
        ModuleBase::WARNING_QUIT("rdmft::write_rdmft_restart",
                                 "occ_flat size does not match nk*nbands.");
    }
    const int is_complex = std::is_same<TK, std::complex<double>>::value ? 1 : 0;
    const int gamma_int = gamma_only ? 1 : 0;
    const int nbasis_global = pv.get_wfc_global_nbasis();

    if (GlobalV::MY_RANK == 0)
    {
        const std::string meta_path = GlobalC::restart.folder + "Restart_rdmft.txt";
        std::ofstream ofs(meta_path.c_str(), std::ofstream::out | std::ofstream::trunc);
        if (!ofs)
        {
            ModuleBase::WARNING_QUIT("rdmft::write_rdmft_restart",
                                     "Cannot open " + meta_path + " for write.");
        }
        ofs << "ABACUS_RDMFT_RESTART 1\n";
        ofs << nk << " " << nbands << " " << nspin << " " << gamma_int << " " << is_complex << " "
            << nbasis_global << " " << nocc << "\n";
        ofs << std::scientific << std::setprecision(17) << n_electrons << "\n";
        ofs << rdmft_functional << "\n";
        ofs.close();

        if (!GlobalC::restart.save_disk("rdmft_occ", 0, nocc, const_cast<double*>(occ_flat.data())))
        {
            ModuleBase::WARNING_QUIT("rdmft::write_rdmft_restart", "Failed to write rdmft_occ binary.");
        }
    }

    std::vector<double> flat;
    psi_to_flat(wfc, flat);
    if (!GlobalC::restart.save_disk("rdmft_wfc", 0, static_cast<int>(flat.size()), flat.data()))
    {
        ModuleBase::WARNING_QUIT("rdmft::write_rdmft_restart", "Failed to write rdmft_wfc binary.");
    }
}

namespace detail
{
inline bool read_rdmft_meta_rank0(const std::string& folder,
                                  int& nk,
                                  int& nbands,
                                  int& nspin,
                                  int& gamma_only,
                                  int& is_complex,
                                  int& nbasis_global,
                                  int& nocc,
                                  double& n_electrons,
                                  std::string& rdmft_functional)
{
    const std::string meta_path = folder + "Restart_rdmft.txt";
    std::ifstream ifs(meta_path.c_str());
    if (!ifs)
    {
        return false;
    }
    std::string magic;
    int ver = 0;
    if (!(ifs >> magic >> ver) || magic != "ABACUS_RDMFT_RESTART" || ver != 1)
    {
        return false;
    }
    if (!(ifs >> nk >> nbands >> nspin >> gamma_only >> is_complex >> nbasis_global >> nocc))
    {
        return false;
    }
    if (nocc != nk * nbands)
    {
        return false;
    }
    if (!(ifs >> n_electrons))
    {
        return false;
    }
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if (!std::getline(ifs, rdmft_functional))
    {
        return false;
    }
    return true;
}
} // namespace detail

template <typename TK>
void read_rdmft_restart(const Parallel_Orbitals& pv,
                        const int nk,
                        const int nbands,
                        const int nspin,
                        const bool gamma_only,
                        const double n_electrons_expected,
                        const std::string& rdmft_functional_expected,
                        std::vector<double>& occ_flat,
                        psi::Psi<TK>& wfc)
{
    if (GlobalC::restart.folder.empty())
    {
        ModuleBase::WARNING_QUIT("rdmft::read_rdmft_restart", "GlobalC::restart.folder is empty.");
    }

    int fnk = 0, fnbands = 0, fnspin = 0, fgamma = 0, fis_complex = 0, fnbasis = 0, fnocc = 0;
    double fn_elec = 0.0;
    std::string fxc;

    bool meta_ok = true;
    if (GlobalV::MY_RANK == 0)
    {
        meta_ok = detail::read_rdmft_meta_rank0(GlobalC::restart.folder, fnk, fnbands, fnspin, fgamma,
                                                fis_complex, fnbasis, fnocc, fn_elec, fxc);
        if (meta_ok)
        {
            const int expect_complex = std::is_same<TK, std::complex<double>>::value ? 1 : 0;
            const int gamma_int = gamma_only ? 1 : 0;
            if (fnk != nk || fnbands != nbands || fnspin != nspin || fgamma != gamma_int
                || fis_complex != expect_complex || fnbasis != pv.get_wfc_global_nbasis())
            {
                meta_ok = false;
            }
            if (std::abs(fn_elec - n_electrons_expected) > 1.0e-8)
            {
                meta_ok = false;
            }
            if (fxc != rdmft_functional_expected)
            {
                meta_ok = false;
            }
        }
    }

#ifdef __MPI
    int ok_i = meta_ok ? 1 : 0;
    Parallel_Common::bcast_int(ok_i);
    if (ok_i == 0)
    {
        ModuleBase::WARNING_QUIT(
            "rdmft::read_rdmft_restart",
            "RDMFT restart metadata mismatch or missing Restart_rdmft.txt (check nk, nbands, nspin, "
            "gamma_only, basis, n_electrons, rdmft_functional, and same MPI layout as the saving run).");
    }
#else
    if (!meta_ok)
    {
        ModuleBase::WARNING_QUIT(
            "rdmft::read_rdmft_restart",
            "RDMFT restart metadata mismatch or missing Restart_rdmft.txt (check nk, nbands, nspin, "
            "gamma_only, basis, n_electrons, rdmft_functional).");
    }
#endif

    const int nocc = nk * nbands;
    occ_flat.resize(static_cast<size_t>(nocc));
    if (GlobalV::MY_RANK == 0)
    {
        if (!GlobalC::restart.load_disk("rdmft_occ", 0, nocc, occ_flat.data(), false))
        {
            ModuleBase::WARNING_QUIT("rdmft::read_rdmft_restart", "Failed to read rdmft_occ binary.");
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_double(occ_flat.data(), nocc);
#endif

    std::vector<double> flat;
    psi_to_flat(wfc, flat);
    if (!GlobalC::restart.load_disk("rdmft_wfc", 0, static_cast<int>(flat.size()), flat.data(), false))
    {
        ModuleBase::WARNING_QUIT("rdmft::read_rdmft_restart", "Failed to read rdmft_wfc binary.");
    }
    try
    {
        flat_to_psi(flat, wfc);
    }
    catch (const std::exception& e)
    {
        ModuleBase::WARNING_QUIT("rdmft::read_rdmft_restart", std::string("flat_to_psi: ") + e.what());
    }
}

} // namespace rdmft

#endif
