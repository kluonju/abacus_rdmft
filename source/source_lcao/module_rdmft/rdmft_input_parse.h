//==========================================================
// Author: Kai Luo
// Email: kluo@njust.edu.cn
// DATE : April 2026
//==========================================================
// Helpers to parse RDMFT INPUT strings and quit with a running-log message on error.

#ifndef RDMFT_INPUT_PARSE_H
#define RDMFT_INPUT_PARSE_H

#include "source_base/tool_quit.h"
#include "source_lcao/module_rdmft/rdmft_type.h"

#include <cctype>
#include <stdexcept>
#include <string>

namespace rdmft
{

/// Parse `rdmft_functional` after INPUT checks; converts `std::invalid_argument` from
/// `parse_xc_type` into `WARNING_QUIT` so the message appears in the log and the run exits cleanly.
inline XCFunctionalType parse_xc_type_or_quit(const std::string& name, const char* file_tag)
{
    try
    {
        return parse_xc_type(name);
    }
    catch (const std::invalid_argument& e)
    {
        ModuleBase::WARNING_QUIT(file_tag, e.what());
    }
}

/// Parse occupation / orbital / joint optimiser keywords (trim, ASCII-lowercase, common aliases).
/// All-whitespace input selects conjugate gradient (INPUT default behaviour). Any other unrecognised
/// string triggers `WARNING_QUIT` with `file_tag` (e.g. `"ReadInput"` or `"ESolver_KS_LCAO::after_scf"`).
inline OptimizerType parse_optimizer_input_or_quit(const std::string& s_in,
                                                   const char* input_keyword,
                                                   const char* file_tag)
{
    const char* const ws = " \t\n\r\f\v";
    const auto first = s_in.find_first_not_of(ws);
    if (first == std::string::npos)
    {
        return OptimizerType::ConjugateGradient;
    }
    const auto last = s_in.find_last_not_of(ws);
    std::string s = s_in.substr(first, last - first + 1);
    for (char& c : s)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (s == "sd" || s == "steepest" || s == "steepest_descent" || s == "gd")
    {
        return OptimizerType::SteepestDescent;
    }
    if (s == "lbfgs" || s == "l-bfgs" || s == "l_bfgs" || s == "bfgs")
    {
        return OptimizerType::LBFGS;
    }
    if (s == "adam")
    {
        return OptimizerType::Adam;
    }
    if (s == "cg" || s == "conjugate_gradient" || s == "conjugate" || s == "pr" || s == "fr"
        || s == "polak" || s == "fletcher_reeves")
    {
        return OptimizerType::ConjugateGradient;
    }

    const std::string msg = std::string(input_keyword) + " has invalid value \"" + s_in
                            + "\". Allowed: sd, cg, lbfgs, adam (see INPUT documentation for aliases).";
    ModuleBase::WARNING_QUIT(file_tag, msg);
}

/// Parse `rdmft_occ_ls_type` / `rdmft_orb_ls_type` keywords.
inline RdmftLineSearchPreset parse_line_search_preset_or_quit(const std::string& s_in,
                                                              const char* input_keyword,
                                                              const char* file_tag)
{
    const char* const ws = " \t\n\r\f\v";
    const auto first = s_in.find_first_not_of(ws);
    if (first == std::string::npos)
    {
        return RdmftLineSearchPreset::Auto;
    }
    const auto last = s_in.find_last_not_of(ws);
    std::string s = s_in.substr(first, last - first + 1);
    for (char& c : s)
    {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (s == "auto" || s == "default")
    {
        return RdmftLineSearchPreset::Auto;
    }
    if (s == "armijo" || s == "backtrack" || s == "bt")
    {
        return RdmftLineSearchPreset::Armijo;
    }
    if (s == "sw" || s == "strong" || s == "strong_wolfe" || s == "wolfe_strong")
    {
        return RdmftLineSearchPreset::StrongWolfe;
    }
    if (s == "wolfe" || s == "weak" || s == "weak_wolfe" || s == "ww" || s == "wolfe_weak")
    {
        return RdmftLineSearchPreset::WeakWolfe;
    }

    const std::string msg = std::string(input_keyword) + " has invalid value \"" + s_in
                            + "\". Allowed: auto, armijo, sw (strong Wolfe), wolfe (weak Wolfe).";
    ModuleBase::WARNING_QUIT(file_tag, msg);
}

} // namespace rdmft

#endif // RDMFT_INPUT_PARSE_H
