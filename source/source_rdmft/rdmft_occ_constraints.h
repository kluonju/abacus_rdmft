#ifndef SOURCE_RDMFT_RDMFT_OCC_CONSTRAINTS_H
#define SOURCE_RDMFT_RDMFT_OCC_CONSTRAINTS_H

#include <vector>

//! \file rdmft_occ_constraints.h
//! \brief Box + electron-number constraints for the occupation sub-problem.
//!
//! Ported from qe-rdmft/rdmft_occupation.f90 and the EBI parameterisation in
//! qe-rdmft/rdmft_ebi.f90.  Occupations are stored as a flat vector of length
//! nbnd*nks with band index fastest (occ[ib + ik*nbnd]), matching the QE
//! RESHAPE convention.  Per-spin electron-number constraints are supported
//! through the isk() spin label of each (k,spin) row.
//!
//! All reductions here are process-local: in the QE reference they carry
//! mp_sum over the k-point pool, but the redesigned backend gathers the full
//! (k,spin) list on every rank before calling into this layer, so a plain
//! serial sum is both correct and deadlock-free.

namespace rdmft
{

//! Immutable description of the occupation index space and constraint targets.
struct OccConstraints
{
    int nbnd = 0;
    int nks = 0;                 //!< number of (k,spin) rows
    std::vector<double> wk;      //!< BZ weight per row, length nks
    std::vector<int> isk;        //!< spin label (1 or 2) per row, length nks
    bool fix_magnetization = false;
    double n_target = 0.0;       //!< total electron target (weighted sum)
    double n_target_up = 0.0;
    double n_target_down = 0.0;

    int size() const { return nbnd * nks; }

    // --- weighted sums / constraint residual --------------------------------
    double weighted_sum(const std::vector<double>& occ) const;
    double weighted_sum_ispin(const std::vector<double>& occ, int ispin) const;
    double magnetization_electrons(const std::vector<double>& occ) const;
    double constraint_violation(const std::vector<double>& occ) const;

    // --- Euclidean L2 proximal projector P_w --------------------------------
    //! In-place projection onto {0<=n<=1, sum_k w_k sum_i n_ik = target}.
    void proximal_project(std::vector<double>& occ) const;
    //! Uniform-shift proximal projector P_u (auxiliary).
    void proximal_project_uniform(std::vector<double>& occ) const;

    // --- SPG2 stopping map and KKT residual ---------------------------------
    //! ||n - P_w(n - alpha*grad)||_inf.
    double pg_map_grad_inf(const std::vector<double>& occ,
                           const std::vector<double>& grad,
                           double alpha) const;
    //! Box KKT residual max(0, V-W).
    double pg_kkt_residual(const std::vector<double>& occ,
                           const std::vector<double>& grad,
                           double boundary_tol = 1.0e-10) const;

    // --- EBI erf parameterisation -------------------------------------------
    static double ebi_occ_from_arg(double arg);
    static double ebi_occ_prime(double arg);
    static double ebi_erfinv(double y);

    //! Solve the scalar shift mu so the (filtered) weighted EBI sum matches target.
    double ebi_solve_mu(const std::vector<double>& x, double target, int ispin,
                        bool use_filter) const;
    //! Map EBI params x to occupations using per-spin shifts.
    void ebi_params_to_occ(const std::vector<double>& x, double mu_up, double mu_dw,
                           std::vector<double>& occ) const;
    //! Map occupations to EBI params x = erfinv(2n-1).
    void ebi_occ_to_params(const std::vector<double>& occ, std::vector<double>& x) const;
    //! Fold the implicit shift into x (x <- x + mu) so the next mu solve starts near 0.
    void ebi_fold_mu(std::vector<double>& x) const;
    //! Transform dE/dn into dE/dx via the implicit-function chain rule.
    void ebi_transform_gradient(const std::vector<double>& grad_n,
                                const std::vector<double>& x, double mu_up, double mu_dw,
                                std::vector<double>& grad_x) const;
    //! Solve mu(s) for x and write occupations n(x) into occ.
    void ebi_sync_occ_from_x(const std::vector<double>& x, std::vector<double>& occ) const;

  private:
    double eval_weighted_shift_sum(const std::vector<double>& x, double shift, int ispin,
                                   bool use_filter) const;
    double eval_weighted_uniform_sum(const std::vector<double>& x, double mu, int ispin,
                                     bool use_filter) const;
    void proximal_project_core(std::vector<double>& occ, double target_ne, int ispin,
                               bool use_filter) const;
    void proximal_project_uniform_core(std::vector<double>& occ, double target_ne, int ispin,
                                       bool use_filter) const;
    double ebi_weighted_occ_sum(const std::vector<double>& x, double mu, int ispin,
                                bool use_filter) const;
};

} // namespace rdmft

#endif // SOURCE_RDMFT_RDMFT_OCC_CONSTRAINTS_H
