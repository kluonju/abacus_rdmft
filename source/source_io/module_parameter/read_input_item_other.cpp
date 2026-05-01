#include "source_base/global_function.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iostream>

namespace ModuleIO
{
void ReadInput::item_others()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    // non-collinear spin-constrained
    {
        Input_Item item("sc_mag_switch");
        item.annotation = "switch to control spin-constrained DFT";
        item.category = "Spin-Constrained DFT";
        item.type = "Boolean";
        item.description = "Switch to control spin-constrained DFT calculation";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.sc_mag_switch);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_mag_switch)
            {
//                ModuleBase::WARNING_QUIT("ReadInput",
//                                         "This feature is not stable yet and might lead to "
//                                         "erroneous results.\n"
//                                         " Please wait for the official release version.");
                // if (para.input.nspin != 4 && para.input.nspin != 2)
                // {
                //     ModuleBase::WARNING_QUIT("ReadInput", "nspin must be 2 or
                //     4 when sc_mag_switch > 0");
                // }
                // if (para.input.calculation != "scf")
                // {
                //     ModuleBase::WARNING_QUIT("ReadInput", "calculation must
                //     be scf when sc_mag_switch > 0");
                // }
                // if (para.input.nupdown > 0.0)
                // {
                //     ModuleBase::WARNING_QUIT("ReadInput", "nupdown should not
                //     be set when sc_mag_switch > 0");
                // }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("decay_grad_switch");
        item.annotation = "switch to control gradient break condition";
        item.category = "Spin-Constrained DFT";
        item.type = "Boolean";
        item.description = "Switch to control gradient break condition in spin-constrained DFT";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.decay_grad_switch);
        this->add_item(item);
    }
    {
        Input_Item item("sc_thr");
        item.annotation = "Convergence criterion of spin-constrained iteration (RMS) in uB";
        item.category = "Spin-Constrained DFT";
        item.type = "Real";
        item.description = "Convergence criterion of spin-constrained iteration (RMS) in uB";
        item.default_value = "1.0e-6";
        item.unit = "uB";
        item.availability = "sc_mag_switch is true";
        read_sync_double(input.sc_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_thr < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sc_thr must >= 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nsc");
        item.annotation = "Maximal number of spin-constrained iteration";
        item.category = "Spin-Constrained DFT";
        item.type = "Integer";
        item.description = "Maximal number of spin-constrained iteration";
        item.default_value = "100";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_int(input.nsc);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nsc <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nsc must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nsc_min");
        item.annotation = "Minimum number of spin-constrained iteration";
        item.category = "Spin-Constrained DFT";
        item.type = "Integer";
        item.description = "Minimum number of spin-constrained iteration";
        item.default_value = "2";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_int(input.nsc_min);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nsc_min <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nsc_min must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("sc_scf_nmin");
        item.annotation = "Minimum number of outer scf loop before "
                          "initializing lambda loop";
        item.category = "Spin-Constrained DFT";
        item.type = "Integer";
        item.description = "Minimum number of outer scf loop before initializing lambda loop";
        item.default_value = "2";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_int(input.sc_scf_nmin);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_scf_nmin < 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sc_scf_nmin must >= 2");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("alpha_trial");
        item.annotation = "Initial trial step size for lambda in eV/uB^2";
        item.category = "Spin-Constrained DFT";
        item.type = "Real";
        item.description = "Initial trial step size for lambda in eV/uB^2";
        item.default_value = "0.01";
        item.unit = "eV/uB^2";
        item.availability = "sc_mag_switch is true";
        read_sync_double(input.alpha_trial);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.alpha_trial <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "alpha_trial must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("sccut");
        item.annotation = "Maximal step size for lambda in eV/uB";
        item.category = "Spin-Constrained DFT";
        item.type = "Real";
        item.description = "Maximal step size for lambda in eV/uB";
        item.default_value = "3.0";
        item.unit = "eV/uB";
        item.availability = "sc_mag_switch is true";
        read_sync_double(input.sccut);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sccut <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sccut must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("sc_drop_thr");
        item.annotation = "Convergence criterion ratio of lambda iteration in Spin-constrained DFT";
        item.category = "Spin-Constrained DFT";
        item.type = "Real";
        item.description = "Convergence criterion ratio of lambda iteration in Spin-constrained DFT";
        item.default_value = "1.0e-2";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_double(input.sc_drop_thr);
        this->add_item(item);
    }
    {
        Input_Item item("sc_scf_thr");
        item.annotation = "Density error threshold for inner loop of spin-constrained SCF";
        item.category = "Spin-Constrained DFT";
        item.type = "Real";
        item.description = "Density error threshold for inner loop of spin-constrained SCF";
        item.default_value = "1.0e-4";
        item.unit = "";
        item.availability = "sc_mag_switch is true";
        read_sync_double(input.sc_scf_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.sc_scf_thr <= 0.0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "sc_scf_thr must > 0.0");
            }
        };
        this->add_item(item);
    }

    // Quasiatomic Orbital analysis
    {
        Input_Item item("qo_switch");
        item.annotation = "switch to control quasiatomic orbital analysis";
        item.category = "Quasiatomic Orbital (QO) analysis";
        item.type = "Boolean";
        item.description = "Whether to let ABACUS output QO analysis required files";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.qo_switch);
        this->add_item(item);
    }
    {
        Input_Item item("qo_basis");
        item.annotation = "type of QO basis function: hydrogen: hydrogen-like "
                          "basis, pswfc: read basis from pseudopotential";
        item.category = "Quasiatomic Orbital (QO) analysis";
        item.type = "String";
        item.description = R"(Type of QO basis function:
* hydrogen: hydrogen-like basis
* pswfc: read basis from pseudopotential
* szv: single-zeta valence basis)";
        item.default_value = "szv";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.qo_basis);
        this->add_item(item);
    }
    {
        Input_Item item("qo_strategy");
        item.annotation = "strategy to generate generate radial orbitals";
        item.category = "Quasiatomic Orbital (QO) analysis";
        item.type = "Vector of String (1 or n values where n is the number of atomic types)";
        item.description = "Strategy to generate radial orbitals for QO analysis. For hydrogen: energy-valence, for pswfc and szv: all";
        item.default_value = "for hydrogen: energy-valence, for pswfc and szv: all";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.qo_strategy.push_back(item.str_values[i]);
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.qo_strategy.size() != para.input.ntype)
            {
                if (para.input.qo_strategy.size() == 1)
                {
                    para.input.qo_strategy.resize(para.input.ntype, para.input.qo_strategy[0]);
                }
                else
                {
                    std::string default_strategy;
                    if (para.input.qo_basis == "hydrogen")
                    {
                        default_strategy = "energy-valence";
                    }
                    else if ((para.input.qo_basis == "pswfc") || (para.input.qo_basis == "szv"))
                    {
                        default_strategy = "all";
                    }
                    else
                    {
                        ModuleBase::WARNING_QUIT("ReadInput",
                                                 "When setting default values for qo_strategy, "
                                                 "unexpected/unknown "
                                                 "qo_basis is found. Please check it.");
                    }
                    para.input.qo_strategy.resize(para.input.ntype, default_strategy);
                }
            }
        };
        sync_stringvec(input.qo_strategy, para.input.ntype, "all");
        this->add_item(item);
    }
    {
        Input_Item item("qo_screening_coeff");
        item.annotation = "rescale the shape of radial orbitals";
        item.category = "Quasiatomic Orbital (QO) analysis";
        item.type = "Vector of Real (n values where n is the number of atomic types; 1 value allowed for qo_basis=pswfc)";
        item.description = "The screening coefficient for each atom type to rescale the shape of radial orbitals";
        item.default_value = "0.1";
        item.unit = "Bohr^-1";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.qo_screening_coeff.push_back(std::stod(item.str_values[i]));
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (!item.is_read())
            {
                return;
            }
            if (para.input.qo_screening_coeff.size() != para.input.ntype)
            {
                if (para.input.qo_basis == "pswfc")
                {
                    double default_screening_coeff
                        = (para.input.qo_screening_coeff.size() == 1) ? para.input.qo_screening_coeff[0] : 0.1;
                    para.input.qo_screening_coeff.resize(para.input.ntype, default_screening_coeff);
                }
                else
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "qo_screening_coeff should have the same number of "
                                             "elements as ntype");
                }
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            for (auto screen_coeff: para.input.qo_screening_coeff)
            {
                if (screen_coeff < 0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "screening coefficient must >= 0 "
                                             "to tune the pswfc decay");
                }
                if (std::fabs(screen_coeff) < 1e-6)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "every low screening coefficient might yield very high "
                                             "computational cost");
                }
            }
        };
        sync_doublevec(input.qo_screening_coeff, para.input.ntype, 0.1);
        this->add_item(item);
    }
    {
        Input_Item item("qo_thr");
        item.annotation = "accuracy for evaluating cutoff radius of QO basis function";
        item.category = "Quasiatomic Orbital (QO) analysis";
        item.type = "Real";
        item.description = "The convergence threshold determining the cutoff of generated orbital. Lower threshold will yield orbital with larger cutoff radius.";
        item.default_value = "1.0e-6";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.qo_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.qo_thr > 1e-6)
            {
                ModuleBase::WARNING("ReadInput",
                                    "too high the convergence threshold might "
                                    "yield unacceptable result");
            }
        };
        this->add_item(item);
    }

    // PEXSI
    {
        Input_Item item("pexsi_npole");
        item.annotation = "Number of poles in expansion";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "The number of poles used in the pole expansion method, should be a even number.";
        item.default_value = "40";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_npole);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_inertia");
        item.annotation = "Whether inertia counting is used at the very "
                          "beginning of PEXSI process";
        item.category = "PEXSI";
        item.type = "Boolean";
        item.description = "Whether inertia counting is used at the very beginning.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pexsi_inertia);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nmax");
        item.annotation = "Maximum number of PEXSI iterations after each "
                          "inertia counting procedure";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "Maximum number of PEXSI iterations after each inertia counting procedure.";
        item.default_value = "80";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_comm");
        item.annotation = "Whether to construct PSelInv communication pattern";
        item.category = "PEXSI";
        item.type = "Boolean";
        item.description = "Whether to construct PSelInv communication pattern.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pexsi_comm);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_storage");
        item.annotation = "Storage space used by the Selected Inversion "
                          "algorithm for symmetric matrices";
        item.category = "PEXSI";
        item.type = "Boolean";
        item.description = "Whether to use symmetric storage space used by the Selected Inversion algorithm for symmetric matrices.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pexsi_storage);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_ordering");
        item.annotation = "Ordering strategy for factorization and selected inversion";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "Ordering strategy for factorization and selected inversion. 0: Parallel ordering using ParMETIS, 1: Sequential ordering using METIS, 2: Multiple minimum degree ordering";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_ordering);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_row_ordering");
        item.annotation = "Row permutation strategy for factorization and "
                          "selected inversion, 0: NoRowPerm, 1: LargeDiag";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "Row permutation strategy for factorization and selected inversion, 0: No row permutation, 1: Make the diagonal entry of the matrix larger than the off-diagonal entries.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_row_ordering);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nproc");
        item.annotation = "Number of processors for parmetis";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "Number of processors for PARMETIS. Only used if pexsi_ordering == 0.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_nproc);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_symm");
        item.annotation = "Matrix symmetry";
        item.category = "PEXSI";
        item.type = "Boolean";
        item.description = "Whether the matrix is symmetric.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pexsi_symm);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_trans");
        item.annotation = "Whether to transpose";
        item.category = "PEXSI";
        item.type = "Boolean";
        item.description = "Whether to factorize the transpose of the matrix.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.pexsi_trans);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_method");
        item.annotation = "pole expansion method, 1: Cauchy Contour Integral, "
                          "2: Moussa optimized method";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "The pole expansion method to be used. 1 for Cauchy Contour Integral method, 2 for Moussa optimized method.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_method);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_nproc_pole");
        item.annotation = "Number of processes used by each pole";
        item.category = "PEXSI";
        item.type = "Integer";
        item.description = "The point parallelizaion of PEXSI. Recommend two points parallelization.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.pexsi_nproc_pole);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_temp");
        item.annotation = "Temperature, in the same unit as H";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Temperature in Fermi-Dirac distribution, in Ry, should have the same effect as the smearing sigma when smearing method is set to Fermi-Dirac.";
        item.default_value = "0.015";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_temp);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_gap");
        item.annotation = "Spectral gap";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Spectral gap, this can be set to be 0 in most cases.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_gap);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_delta_e");
        item.annotation = "An upper bound for the spectral radius of S^{-1} H";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Upper bound for the spectral radius of S^{-1}H.";
        item.default_value = "20";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_delta_e);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_lower");
        item.annotation = "Initial guess of lower bound for mu";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Initial guess of lower bound for mu.";
        item.default_value = "-10";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu_lower);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_upper");
        item.annotation = "Initial guess of upper bound for mu";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Initial guess of upper bound for mu.";
        item.default_value = "10";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu_upper);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu");
        item.annotation = "Initial guess for mu (for the solver)";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Initial guess for mu (for the solver).";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_thr");
        item.annotation = "Stopping criterion in terms of the chemical "
                          "potential for the inertia counting procedure";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Stopping criterion in terms of the chemical potential for the inertia counting procedure.";
        item.default_value = "0.05";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu_thr);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_expand");
        item.annotation = "If the chemical potential is not in the initial "
                          "interval, the interval is expanded by "
                          "muInertiaExpansion";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "If the chemical potential is not in the initial interval, the interval is expanded by this value.";
        item.default_value = "0.3";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu_expand);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_mu_guard");
        item.annotation = "Safe guard criterion in terms of the chemical potential to "
                          "reinvoke the inertia counting procedure";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Safe guard criterion in terms of the chemical potential to reinvoke the inertia counting procedure.";
        item.default_value = "0.2";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_mu_guard);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_elec_thr");
        item.annotation = "Stopping criterion of the PEXSI iteration in terms "
                          "of the number of electrons compared to "
                          "numElectronExact";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "Stopping criterion of the PEXSI iteration in terms of the number of electrons compared to numElectronExact.";
        item.default_value = "0.001";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_elec_thr);
        this->add_item(item);
    }
    {
        Input_Item item("pexsi_zero_thr");
        item.annotation = "if the absolute value of matrix element is less "
                          "than ZERO_Limit, it will be considered as 0";
        item.category = "PEXSI";
        item.type = "Real";
        item.description = "if the absolute value of CCS matrix element is less than this value, it will be considered as zero.";
        item.default_value = "1e-10";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.pexsi_zero_thr);
        this->add_item(item);
    }

    // Only for Test
    {
        Input_Item item("out_alllog");
        item.annotation = "output information for each processor, when parallel";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print information into individual logs from all ranks in an MPI run.\n* True: Information from each rank will be written into individual files named OUT.{calculation}_{suffix}/running_${calculation}.log.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_alllog);
        this->add_item(item);
    }
    {
        Input_Item item("nurse");
        item.annotation = "for coders";
        item.category = "Variables useful for debugging";
        item.type = "Integer";
        item.description = "Debugging flag for developers";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.nurse);
        this->add_item(item);
    }
    {
        Input_Item item("t_in_h");
        item.annotation = "calculate the kinetic energy or not";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to include kinetic term in obtaining the Hamiltonian matrix.
* 0: No.
* 1: Yes.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.t_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vl_in_h");
        item.annotation = "calculate the local potential or not";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to include local pseudopotential term in obtaining the Hamiltonian matrix.
* 0: No.
* 1: Yes.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.vl_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vnl_in_h");
        item.annotation = "calculate the nonlocal potential or not";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to include non-local pseudopotential term in obtaining the Hamiltonian matrix.
* 0: No.
* 1: Yes.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.vnl_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vh_in_h");
        item.annotation = "calculate the hartree potential or not";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to include Hartree potential term in obtaining the Hamiltonian matrix.
* 0: No.
* 1: Yes.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.vh_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("vion_in_h");
        item.annotation = "calculate the local ionic potential or not";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to include local ionic potential term in obtaining the Hamiltonian matrix.
* 0: No.
* 1: Yes.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.vion_in_h);
        this->add_item(item);
    }
    {
        Input_Item item("test_force");
        item.annotation = "test the force";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to output the detailed components in forces.
* 0: No.
* 1: Yes.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.test_force);
        this->add_item(item);
    }
    {
        Input_Item item("test_stress");
        item.annotation = "test the stress";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to output the detailed components in stress.
* 0: No.
* 1: Yes.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.test_stress);
        this->add_item(item);
    }
    {
        Input_Item item("test_skip_ewald");
        item.annotation = "whether to skip ewald";
        item.category = "Variables useful for debugging";
        item.type = "Boolean";
        item.description = R"(Specify whether to skip the calculation of the ewald energy.
* 0: No.
* 1: Yes.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.test_skip_ewald);
        this->add_item(item);
    }
    {
        Input_Item item("ri_hartree_benchmark");
        item.annotation = "whether to use the RI approximation for the Hartree term in LR-TDDFT for benchmark (with FHI-aims/ABACUS read-in style)";
        item.category = "Linear Response TDDFT";
        item.type = "String";
        item.description = "Whether to use the RI approximation for the Hartree term in LR-TDDFT for benchmark (with FHI-aims/ABACUS read-in style)";
        item.default_value = "none";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.ri_hartree_benchmark);
        this->add_item(item);
    }
    {
        Input_Item item("aims_nbasis");
        item.annotation = "the number of basis functions for each atom type used in FHI-aims (for benchmark)";
        item.category = "Linear Response TDDFT";
        item.type = "A number(ntype) of Integers";
        item.description = "Atomic basis set size for each atom type (with the same order as in STRU) in FHI-aims.";
        item.default_value = "{} (empty list, where ABACUS use its own basis set size)";
        item.unit = "";
        item.availability = "ri_hartree_benchmark = aims";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.aims_nbasis.push_back(std::stod(item.str_values[i]));
            }
            };
        sync_intvec(input.aims_nbasis, para.input.aims_nbasis.size(), 0);
        this->add_item(item);
    }

    // RDMFT, added by jghan, 2024-10-16
#ifdef __RDMFT
    {
        Input_Item item("rdmft");
        item.annotation = "whether to perform rdmft calculation, default is false";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Boolean";
        item.description = "Whether to perform rdmft calculation (reduced density matrix funcional theory)";
        item.default_value = "false";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.rdmft);
        this->add_item(item);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && para.input.nspin == 4)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "rdmft is not available for nspin = 4");
            }
        };
    }
    {
        Input_Item item("rdmft_power_alpha");
        item.annotation = "the alpha parameter of power-functional, g(occ_number) = occ_number^alpha"
                          " used in exx-type functionals such as muller and power";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "The alpha parameter of power-functional(or other exx-type/hybrid functionals) which used in RDMFT, g(occ_number) = occ_number^alpha";
        item.default_value = "0.656";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.rdmft_power_alpha);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            const std::string& func = para.input.rdmft_functional.empty()
                                      ? para.input.dft_functional
                                      : para.input.rdmft_functional;
            if( func == "hf" || func == "pbe0" )
            {
                para.input.rdmft_power_alpha = 1.0;
            }
            else if( func == "muller" )
            {
                para.input.rdmft_power_alpha = 0.5;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if( (para.input.rdmft_power_alpha < 0) || (para.input.rdmft_power_alpha > 1) )
            {
                ModuleBase::WARNING_QUIT("ReadInput", "rdmft_power_alpha should be greater than 0.0 and less than 1.0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_functional");
        item.annotation = "RDMFT exchange-correlation functional: hf, muller, power, gu, bbc3, geo, optgm";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "String";
        item.description = "RDMFT XC functional used for the RDMFT optimisation stage. "
                           "Supported values: hf (Hartree-Fock), muller (Muller/BBC1), "
                           "power (power functional, uses rdmft_power_alpha), gu (Goedecker-Umrigar), "
                           "bbc3 (BBC3-inspired rank-separated), geo "
                           "(f(n_p, n_q) = [n_p n_q + sqrt(n_p n_q) + 2 (n_p n_q)^{3/4}]/4), "
                           "optgm (same three powers as geo with calibrated mixture weights). "
                           "If empty, the old single-step RDMFT code path is used.";
        item.default_value = "";
        item.unit = "";
        item.availability = "rdmft == true";
        read_sync_string(input.rdmft_functional);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                const std::string& f = para.input.rdmft_functional;
                if (f != "hf" && f != "muller" && f != "power" && f != "gu"
                    && f != "bbc3" && f != "geo" && f != "optgm")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_functional must be one of: hf, muller, power, gu, bbc3, geo, optgm");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_solver_strategy");
        item.annotation = "RDMFT optimisation strategy: alternating or joint";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "String";
        item.description = "Strategy for the RDMFT energy minimisation. "
                           "'alternating': alternate between occupation and orbital sub-problems. "
                           "'joint' (formerly 'product_manifold'): pack occupations and orbitals into a "
                           "single point on the product manifold and optimise them simultaneously. "
                           "The legacy value 'product_manifold' is still accepted as an alias.";
        item.default_value = "alternating";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            std::string v = strvalue;
            if (v == "product_manifold") v = "joint";
            para.input.rdmft_solver_strategy = v;
        };
        sync_string(input.rdmft_solver_strategy);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                const std::string& s = para.input.rdmft_solver_strategy;
                if (s != "alternating" && s != "joint")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_solver_strategy must be 'alternating' or 'joint' "
                        "(legacy alias 'product_manifold' is also accepted)");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_occ_optimizer");
        item.annotation = "Optimiser for occupation numbers in RDMFT: sd, cg, lbfgs, adam";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "String";
        item.description = "Gradient-based optimiser for the occupation sub-problem. "
                           "sd: steepest descent, cg: conjugate gradient, lbfgs: lbfgs, adam: Adam. "
                           "Used by augmented_lagrangian only (Strong Wolfe if lbfgs, Armijo otherwise). "
                           "Ignored by projected_gradient / active_set: these route to the textbook Spectral "
                           "Projected Gradient method, which uses a Barzilai-Borwein spectral step and a "
                           "non-monotone Armijo line search and does not consult this keyword.";
        item.default_value = "cg";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_string(input.rdmft_occ_optimizer);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_orb_optimizer");
        item.annotation = "Optimiser for orbitals in RDMFT: sd, cg, lbfgs, adam";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "String";
        item.description = "Gradient-based optimiser for the orbital (Stiefel manifold) sub-problem. "
                           "sd: steepest descent, cg: conjugate gradient, lbfgs: lbfgs, adam: Adam. "
                           "Alternating orbital inner always uses Armijo line search (no Strong Wolfe).";
        item.default_value = "cg";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_string(input.rdmft_orb_optimizer);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_joint_optimizer");
        item.annotation = "Single unified optimiser for the RDMFT joint strategy: sd, cg, lbfgs, adam";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "String";
        item.description = "For rdmft_solver_strategy = joint, the occupation parameters and "
                           "orbital coefficients are packed into one point on the product manifold "
                           "and stepped simultaneously by a single optimiser of this type. "
                           "Allowed values: sd (steepest descent), cg (conjugate gradient), "
                           "lbfgs, adam. Joint line search: Strong Wolfe if lbfgs, Armijo if sd/cg/adam.";
        item.default_value = "lbfgs";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_string(input.rdmft_joint_optimizer);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                const std::string& s = para.input.rdmft_joint_optimizer;
                if (s.empty())
                {
                    return;
                }
                std::string t = s;
                for (char& c : t)
                {
                    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
                }
                if (t != "sd" && t != "cg" && t != "lbfgs" && t != "adam")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_joint_optimizer must be one of: sd, cg, lbfgs, adam");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_outer_maxiter");
        item.annotation = "Maximum RDMFT outer iterations (alternating / product-manifold)";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Integer";
        item.description = "Maximum number of outer cycles: each cycle is one occupation sub-problem "
                           "plus one orbital sub-problem (alternating), or one joint step (product_manifold).";
        item.default_value = "200";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_int(input.rdmft_outer_maxiter);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_orb_maxiter");
        item.annotation = "Maximum iterations per RDMFT orbital sub-problem";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Integer";
        item.description = "Maximum inner iterations while optimising orbitals with occupations fixed "
                           "(Stiefel manifold step).";
        item.default_value = "50";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_int(input.rdmft_orb_maxiter);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_occ_maxiter");
        item.annotation = "Maximum occupation iterations per RDMFT outer step";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Integer";
        item.description = "Maximum iterations for the occupation sub-problem within one outer step.";
        item.default_value = "50";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_int(input.rdmft_occ_maxiter);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_occ_init_mode");
        item.annotation = "Initial occupation source for RDMFT: ks, perturbed, binary, uniform";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "String";
        item.description = "ks: use KS occupations directly. "
                   "perturbed: in a Fermi window, apply +/-delta based on the KS occupation. "
                   "binary: in the same Fermi window, set to 1-delta if occ>=0.5, or delta if occ<0.5. "
                   "uniform: in the same Fermi window, set all selected occupations "
                   "to the local average N_top / N_selected.";
        item.default_value = "ks";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_string(input.rdmft_occ_init_mode);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                const std::string& s = para.input.rdmft_occ_init_mode;
                if (s != "ks" && s != "perturbed" && s != "binary" && s != "uniform")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                    "rdmft_occ_init_mode must be one of: ks, perturbed, binary, uniform");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_occ_init_perturb");
        item.annotation = "Initial occupation perturbation magnitude for RDMFT";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Perturbation delta used when rdmft_occ_init_mode = perturbed. "
                           "In each k-point, selected bands above Fermi are increased by +delta "
                           "and selected bands below Fermi are reduced by -delta, then projected "
                           "to the feasible set. 0 disables the perturbation.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_occ_init_perturb);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                if (para.input.rdmft_occ_init_perturb < 0.0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_occ_init_perturb must be >= 0.0");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_occ_init_nbands_top");
        item.annotation = "Initial occupation Fermi-window half-width";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Integer";
        item.description = "Fermi-window half-width K used by rdmft_occ_init_mode = perturbed/uniform. "
                           "For each k-point, select K bands above and K bands below the Fermi boundary. "
                           "K <= 0 disables these windowed initialisation modes and falls back to KS occupations.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_int(input.rdmft_occ_init_nbands_top);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                if (para.input.rdmft_occ_init_nbands_top < 0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_occ_init_nbands_top must be >= 0");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_energy_tol");
        item.annotation = "RDMFT convergence threshold on energy change (Ry)";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Alternating and joint RDMFT: after the first outer iteration, the outer loop declares "
                           "convergence if |E - E_prev| (Ry) < rdmft_energy_tol (used only when > 0) **or** both "
                           "occupation- and orbital-inner criteria are met (alternating: inner `converged` flags; "
                           "joint: occ/orb stationarity flags). If <= 0, the energy test is omitted and outer "
                           "convergence uses the inner criteria only.";
        item.default_value = "1e-8";
        item.unit = "Ry";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_energy_tol);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_orb_grad_tol");
        item.annotation = "RDMFT convergence threshold on gradient norm";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Orbital inner loop: Riemannian gradient norm ||G_R|| convergence threshold.";
        item.default_value = "1e-6";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_orb_grad_tol);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_orb_energy_tol");
        item.annotation = "RDMFT orbital inner: |dE| convergence (Ry)";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Alternating strategy, orbital inner loop: when this value is > 0, "
                           "convergence requires both ||G_R|| < rdmft_orb_grad_tol and energy stability: "
                           "|E_k - E_{k-1}| before the step and |E_new - E| after an accepted line-search "
                           "step must be below this (Ry). Use <= 0 to disable the energy criterion "
                           "(gradient-only stopping).";
        item.default_value = "1e-8";
        item.unit = "Ry";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_orb_energy_tol);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_occ_tol");
        item.annotation = "RDMFT occupation inner: sum of |Δn| convergence";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Occupation inner: stop a step when sum_i |n_i^{new}-n_i^{old}| in "
                           "one inner iteration is below this. Used by augmented_lagrangian, "
                           "projected_gradient, and active_set paths (PG / AS additionally exit "
                           "the inner loop on rdmft_occ_grad_tol vs ||g_proj||_inf).";
        item.default_value = "1e-6";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_occ_tol);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_occ_grad_tol");
        item.annotation = "RDMFT occupation inner: gradient norm threshold (non-ALM / joint non-ALM)";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Projected-gradient (SPG) occupation inner loop converges when "
                           "||r||_inf = ||n - P(n - ∇E)||_inf is below this (Bertsekas projected-gradient "
                           "residual at unit step, the textbook SPG stationarity measure). Also used for ALM "
                           "first-inner ||dL/dp||, joint, and other gradient checks as in the solver.";
        item.default_value = "1e-6";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_occ_grad_tol);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_occ_energy_tol");
        item.annotation = "RDMFT PG: unused (inner stop is ||g_proj|| vs rdmft_occ_grad_tol)";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Not used when rdmft_constraint = projected_gradient: PG occupation inner loop "
                           "always stops on ||g_proj||_inf < rdmft_occ_grad_tol (g_proj = (n - P(n - τ∇E))/τ with "
                           "τ from the occupation line search as for rdmft_occ_grad_tol). Kept for backward-compatible "
                           "INPUT files.";
        item.default_value = "1e-8";
        item.unit = "Ry";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_occ_energy_tol);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_occ_entropy_gamma");
        item.annotation = "RDMFT HF: binary-entropy regularization γ for occupation optimization";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "When rdmft_functional is hf and this is > 0, adds γ·Σ_k w_k Σ_i (n ln n + (1-n)ln(1-n)) "
                           "to the objective and matching ∂E/∂n so CG/lbfgs has curvature (typical 1e-6). "
                           "Ignored for non-HF functionals. Use 0 for muller (default).";
        item.default_value = "0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_occ_entropy_gamma);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_occ_param");
        item.annotation = "Occupation parameterisation for RDMFT: cosine_sq, logistic, sigma_shift";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "String";
        item.description = "Occupation map n in [0,1] from unconstrained parameters. "
                           "cosine_sq: n = cos^2(theta). "
                           "logistic: n = sigma(x) = 1/(1+exp(-x)) per state; with augmented_lagrangian, "
                           "the electron-number constraint is enforced by the ALM term (not a global mu solve). "
                           "sigma_shift: n_{ik} = sigma(z_{ik}+lambda) where lambda is solved each step by "
                           "bisection to satisfy sum_k w_k sum_i n_{ik} = N_e exactly; no ALM penalty needed.";
        item.default_value = "cosine_sq";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_string(input.rdmft_occ_param);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                const std::string& s = para.input.rdmft_occ_param;
                if (s != "cosine_sq" && s != "logistic" && s != "sigma_shift")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_occ_param must be 'cosine_sq', 'logistic', or 'sigma_shift'");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_constraint");
        item.annotation = "Electron-number constraint method for RDMFT";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "String";
        item.description = "Method to enforce the electron-number constraint sum_k w_k sum_i n_ik = N_e. "
                           "augmented_lagrangian: ALM on the occupation parameters (rdmft_occ_param cosine_sq or logistic) with Armijo line search. "
                           "projected_gradient: Spectral Projected Gradient (Birgin-Martinez-Raydan, SIOPT 2000) "
                           "in occupation space with non-monotone Armijo line search (Grippo-Lampariello-Lucidi). "
                           "active_set: alias for projected_gradient in the current implementation -- kept for backward INPUT compatibility.";
        item.default_value = "augmented_lagrangian";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_string(input.rdmft_constraint);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                const std::string& s = para.input.rdmft_constraint;
                if (s != "augmented_lagrangian" && s != "projected_gradient" && s != "active_set")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_constraint must be 'augmented_lagrangian', 'projected_gradient', or 'active_set'");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_alm_lambda_init");
        item.annotation = "Initial ALM Lagrange multiplier lambda for RDMFT occupations";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Initial lambda value for the augmented Lagrangian method "
                           "(used only when rdmft_constraint = augmented_lagrangian).";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_alm_lambda_init);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_alm_mu_init");
        item.annotation = "Initial ALM penalty parameter mu for RDMFT occupations";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Initial mu value for the augmented Lagrangian method "
                           "(used only when rdmft_constraint = augmented_lagrangian).";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_alm_mu_init);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                if (para.input.rdmft_alm_mu_init <= 0.0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_alm_mu_init must be > 0.0");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_alm_mu_factor");
        item.annotation = "ALM mu growth factor per update for RDMFT occupations";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Multiplicative factor for updating mu in augmented Lagrangian mode. "
                           "Effective update is mu <- min(mu * factor, mu_max).";
        item.default_value = "2.0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_alm_mu_factor);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                if (para.input.rdmft_alm_mu_factor < 1.0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_alm_mu_factor must be >= 1.0");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_alpha_step");
        item.annotation = "Initial line-search step length for RDMFT";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Initial trial step length for the Armijo backtracking line search in RDMFT. Used by "
                           "ALM occupation Armijo and the orbital Armijo line search; the SPG (projected_gradient / "
                           "active_set) occupation block manages its own spectral (Barzilai-Borwein) step length and "
                           "does not consult this value.";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_alpha_step);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_line_search_polynomial");
        item.annotation = "Use polynomial (quadratic/cubic) step in RDMFT Armijo line search";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Boolean";
        item.description = "If true, after a failed Armijo trial the next step is suggested by a "
                           "quadratic model (first failure) and a cubic (later failures) along the "
                           "line; if false, use geometric reduction (multiply by rdmft_line_search_rho).";
        item.default_value = "true";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_bool(input.rdmft_line_search_polynomial);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_line_search_c1");
        item.annotation = "Armijo sufficient-decrease c1 for RDMFT line searches";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Armijo constant in (0, 1): E_trial <= E + c1 * alpha * (directional derivative). "
                           "Used for alternating orbital Armijo, ALM/PG/AS occupation line searches, joint Armijo, "
                           "and the Armijo condition inside Strong Wolfe. Smaller values demand a steeper energy cut.";
        item.default_value = "1e-4";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_line_search_c1);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                if (para.input.rdmft_line_search_c1 <= 0.0 || para.input.rdmft_line_search_c1 >= 1.0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "rdmft_line_search_c1 must be in (0, 1)");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_line_search_c2");
        item.annotation = "Strong Wolfe curvature parameter c2 for RDMFT (lbfgs line search)";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Used with strong Wolfe when rdmft_occ_optimizer = lbfgs (ALM) or "
                           "joint + lbfgs: require |g · d| <= c2 |g0 · d| at the trial point.";
        item.default_value = "0.9";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_line_search_c2);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                if (para.input.rdmft_line_search_c2 < 0.0 || para.input.rdmft_line_search_c2 > 1.0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "rdmft_line_search_c2 must be in [0, 1]");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_line_search_max_zoom");
        item.annotation = "Maximum zoom iterations in RDMFT Strong Wolfe line search";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Int";
        item.description = "Zoom sub-iteration cap (Nocedal & Wright) for the strong Wolfe search.";
        item.default_value = "20";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_int(input.rdmft_line_search_max_zoom);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                if (para.input.rdmft_line_search_max_zoom < 1)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "rdmft_line_search_max_zoom must be >= 1");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_alm_bb_enabled");
        item.annotation = "Enable Barzilai-Borwein step seed for ALM occupations";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Boolean";
        item.description = "If true, ALM occupation optimization seeds Armijo backtracking with a "
                           "Barzilai-Borwein step estimate in occupation-parameter space.";
        item.default_value = "true";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_bool(input.rdmft_alm_bb_enabled);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_alm_bb_mode");
        item.annotation = "Barzilai-Borwein mode for ALM occupations: bb1, bb2, alternate";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "String";
        item.description = "BB seed policy for ALM occupation Armijo initial step. "
                           "bb1: alpha = (s^T s)/(s^T y); bb2: alpha = (s^T y)/(y^T y); "
                           "alternate: alternate bb1 and bb2 each inner iteration.";
        item.default_value = "alternate";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_string(input.rdmft_alm_bb_mode);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                const std::string& s = para.input.rdmft_alm_bb_mode;
                if (s != "bb1" && s != "bb2" && s != "alternate")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_alm_bb_mode must be one of: bb1, bb2, alternate");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_alm_bb_alpha_min");
        item.annotation = "Lower bound for ALM BB step seed";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Lower clamp for the Barzilai-Borwein Armijo initial step in ALM occupation optimization.";
        item.default_value = "1e-8";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_alm_bb_alpha_min);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                if (para.input.rdmft_alm_bb_alpha_min <= 0.0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_alm_bb_alpha_min must be > 0.0");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_alm_bb_alpha_max");
        item.annotation = "Upper bound for ALM BB step seed";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Upper clamp for the Barzilai-Borwein Armijo initial step in ALM occupation optimization.";
        item.default_value = "10.0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_alm_bb_alpha_max);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                if (para.input.rdmft_alm_bb_alpha_max <= 0.0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_alm_bb_alpha_max must be > 0.0");
                }
                if (para.input.rdmft_alm_bb_alpha_max < para.input.rdmft_alm_bb_alpha_min)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                        "rdmft_alm_bb_alpha_max must be >= rdmft_alm_bb_alpha_min");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_lbfgs_memory");
        item.annotation = "lbfgs history vectors for RDMFT";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Integer";
        item.description = "Number of past gradient/step pairs stored by lbfgs optimiser.";
        item.default_value = "10";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_int(input.rdmft_lbfgs_memory);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_adam_lr");
        item.annotation = "Adam learning rate for RDMFT";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Learning rate for the Adam optimiser in RDMFT.";
        item.default_value = "0.001";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_adam_lr);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_grad_check");
        item.annotation = "Finite-difference gradient check before RDMFT optimisation";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Boolean";
        item.description = "If true, run a finite-difference gradient verification before the "
                           "RDMFT optimisation begins. Results are written to the running log.";
        item.default_value = "false";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_bool(input.rdmft_grad_check);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_print_stiefel_gram");
        item.annotation = "Print per-k-point Stiefel Gram residual ||G_k - I||_F (alternating strategy)";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Boolean";
        item.description = "If true, each RDMFT alternating outer iteration logs the Frobenius norm of "
                           "G_k - I (G_k = X_k^H X_k or C_k^H S_k C_k). This runs an extra distributed "
                           "Gram-matrix multiply per k-point per outer step and increases memory traffic; "
                           "leave false for production runs.";
        item.default_value = "false";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_bool(input.rdmft_print_stiefel_gram);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_joint_orb_scale");
        item.annotation = "Joint-strategy scaling factor between orbital and occupation blocks";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "Diagonal pre-conditioner on the packed (occupation parameter, "
                           "orbital coefficient) vector of the joint RDMFT optimiser. "
                           "The orbital gradient is multiplied by this factor before the "
                           "optimiser compute_direction, and the resulting orbital "
                           "direction is multiplied by it again before applying the step, "
                           "so the physical orbital step is attenuated/amplified by "
                           "joint_orb_scale^2 relative to the occupation-parameter step. "
                           "Default 1.0 (no rescaling); values less than one damp the "
                           "orbital block if its gradient is much larger than the "
                           "occupation gradient; values greater than one boost the orbital "
                           "block otherwise.";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\" && rdmft_solver_strategy == \"joint\"";
        read_sync_double(input.rdmft_joint_orb_scale);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_nelec_use_input");
        item.annotation = "Use PARAM.inp.nelec as base for RDMFT electron equality constraint";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Boolean";
        item.description = "If false (default), the constraint target is N_e = sum_{ik} wg(ik,ib) + "
                           "rdmft_nelec_delta (KS occupation weights before RDMFT). If true, "
                           "N_e = PARAM.inp.nelec + rdmft_nelec_delta. Use true to align the RDMFT "
                           "constraint with the global charge target when it differs from the "
                           "current wg sum.";
        item.default_value = "false";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_bool(input.rdmft_nelec_use_input);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_nelec_delta");
        item.annotation = "Additive offset for RDMFT electron equality target N_e";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "N_e = (rdmft_nelec_use_input ? PARAM.inp.nelec : sum wg) + rdmft_nelec_delta. "
                           "Must yield a positive N_e at runtime.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_nelec_delta);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty())
            {
                const double d = para.input.rdmft_nelec_delta;
                if (!std::isfinite(d))
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "rdmft_nelec_delta must be finite");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_hybrid_dft_xc");
        item.annotation = "Add semilocal DFT XC on RDMFT density to RDMFT total energy";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Boolean";
        item.description = "When true, after building rho from natural orbitals and occupations, "
                           "evaluate semilocal XC via the same PotXC path as LCAO-DFT (dft_functional) "
                           "and add λ·E_xc[ρ] to the energy with V_xc in the one-body Hamiltonian; "
                           "the RI RDMFT exchange contribution is scaled by (1−λ). Requires LIBXC when "
                           "the chosen dft_functional needs it.";
        item.default_value = "false";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_bool(input.rdmft_hybrid_dft_xc);
        this->add_item(item);
    }
    {
        Input_Item item("rdmft_hybrid_dft_xc_lambda");
        item.annotation = "Weight λ for semilocal DFT XC vs RDMFT RI exchange (0 to 1)";
        item.category = "Reduced Density Matrix Functional Theory";
        item.type = "Real";
        item.description = "E += λ·E_xc^DFT[ρ]; RDMFT exchange energy and occupation/orbital exchange "
                           "gradients use coupling scaled by (1−λ). Ignored when rdmft_hybrid_dft_xc is false.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "rdmft == true && rdmft_functional != \"\"";
        read_sync_double(input.rdmft_hybrid_dft_xc_lambda);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.rdmft && !para.input.rdmft_functional.empty() && para.input.rdmft_hybrid_dft_xc)
            {
                const double x = para.input.rdmft_hybrid_dft_xc_lambda;
                if (x < 0.0 || x > 1.0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "rdmft_hybrid_dft_xc_lambda must be in [0, 1]");
                }
            }
        };
        this->add_item(item);
    }
#endif

    // EXX PW by rhx0820, 2025-03-10
    {
        Input_Item item("exxace");
        item.annotation = "whether to perform ace calculation in exxpw";
        item.category = "Exact Exchange (PW)";
        item.type = "Boolean";
        item.description = R"(Whether to use the ACE method (https://doi.org/10.1021/acs.jctc.6b00092) to accelerate the calculation the Fock exchange matrix. Should be set to true most of the time.
* True: Use the ACE method to calculate the Fock exchange operator.
* False: Use the traditional method to calculate the Fock exchange operator.)";
        item.default_value = "True";
        item.unit = "";
        item.availability = "exx_separate_loop==True.";
        read_sync_bool(input.exxace);
        this->add_item(item);
    }
    {
        Input_Item item("exx_gamma_extrapolation");
        item.annotation = "whether to perform gamma extrapolation in exxpw";
        item.category = "Exact Exchange (PW)";
        item.type = "Boolean";
        item.description = "Whether to use the gamma point extrapolation method to calculate the Fock exchange operator. See https://doi.org/10.1103/PhysRevB.79.205114 for details. Should be set to true most of the time.";
        item.default_value = "True";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.exx_gamma_extrapolation);
        this->add_item(item);
    }
    {
        Input_Item item("ecutexx");
        item.annotation = "energy cutoff for exx calculation, Ry";
        item.category = "Exact Exchange (PW)";
        item.type = "Real";
        item.description = "The energy cutoff for EXX (Fock) exchange operator in plane wave basis calculations. Reducing ecutexx below ecutrho may significantly accelerate EXX computations. This speed improvement comes with a reduced numerical accuracy in the exchange energy calculation.";
        item.default_value = "same as ecutrho";
        item.unit = "Ry";
        item.availability = "";
        read_sync_double(input.ecutexx);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.ecutexx < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ecutexx must >= 0");
            }
        };
        this->add_item(item);
    }

    {
        Input_Item item("exx_thr_type");
        item.annotation = "threshold type for exx outer loop, energy or density";
        item.category = "Exact Exchange (PW)";
        item.type = "String";
        item.description = R"(The type of threshold used to judge whether the outer loop has converged in the separate loop EXX calculation.
* energy: use the change of exact exchange energy to judge convergence.
* density: if the change of charge density difference between two successive outer loop iterations is seen as converged according to scf_thr, then the outer loop is seen as converged.)";
        item.default_value = "density";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.exx_thr_type);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            std::string thr_type = para.input.exx_thr_type;
            std::transform(thr_type.begin(), thr_type.end(), thr_type.begin(), ::tolower);
            if (thr_type != "energy" && thr_type != "density")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exx_thr_type should be energy or density");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exx_ene_thr");
        item.annotation = "threshold for exx outer loop when exx_thr_type = energy";
        item.category = "Exact Exchange (PW)";
        item.type = "Real";
        item.description = "The threshold for the change of exact exchange energy to judge convergence of the outer loop in the separate loop EXX calculation.";
        item.default_value = "1e-5";
        item.unit = "Ry";
        item.availability = "exx_thr_type==energy";
        read_sync_double(input.exx_ene_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.exx_ene_thr <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exx_ene_thr must > 0");
            }
        };
        this->add_item(item);
    }
}
} // namespace ModuleIO
