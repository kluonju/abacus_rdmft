#ifndef RDMFT_LINE_SEARCH_H
#define RDMFT_LINE_SEARCH_H

#include <functional>
#include <vector>

//! \file rdmft_line_search.h
//! \brief Line searches and the Barzilai-Borwein spectral steplength.
//!
//! Ported from qe-rdmft/rdmft_linesearch.f90 and rdmft_spg.f90.  The line
//! search is expressed through an abstract evaluator so the same driver serves
//! the occupation (SPG2 straight chord, EBI x-space) and orbital (Stiefel
//! retraction) inner blocks.

namespace rdmft
{

//! Evaluator signature: given a step alpha it must set the trial energy f,
//! the directional derivative g = phi'(alpha), and an error flag (nonzero to
//! reject the trial, e.g. constraint violation).
using LineSearchEval = std::function<void(double alpha, double& f, double& g, int& ierr)>;

struct LineSearchResult
{
    bool success = false;
    double step = 0.0;
    double f_new = 0.0;
    double g_new = 0.0;
    int n_eval = 0;
};

//! Monotone/nonmonotone Armijo backtracking on a straight segment.
//! f0 is the energy at alpha=0, phi0 = phi'(0) < 0 the initial slope.
//! f_ref allows a nonmonotone reference (GLL / Zhang-Hager); pass f0 for
//! the standard monotone rule.
LineSearchResult armijo_line_search(const LineSearchEval& eval, double f0, double phi0,
                                    double alpha0, double c1, double rho, int max_iter,
                                    double f_ref, double alpha_floor);

//! Strong-Wolfe line search with bracketing + zoom (Nocedal & Wright Alg 3.5/3.6).
LineSearchResult strong_wolfe_line_search(const LineSearchEval& eval, double f0, double phi0,
                                          double alpha0, double c1, double c2, int max_iter,
                                          int max_zoom, double f_ref);

//! Barzilai-Borwein state: remembers the last accepted (x, g) pair.
class BarzilaiBorwein
{
  public:
    void init(double alpha_min, double alpha_max);
    void reset();
    void record(const std::vector<double>& x, const std::vector<double>& g);
    bool has_prev() const { return have_prev_; }

    //! Spectral steplength alpha = clamp((s.s)/(s.y), [amin, amax]); falls
    //! back to alpha_init when s.y <= 0 (BMR SPG2 Algorithm 2.2, Step 3).
    double spectral_step(const std::vector<double>& x, const std::vector<double>& g,
                         double alpha_init) const;

  private:
    bool have_prev_ = false;
    double alpha_min_ = 1.0e-8;
    double alpha_max_ = 1.0e2;
    std::vector<double> x_prev_;
    std::vector<double> g_prev_;
};

} // namespace rdmft

#endif // RDMFT_LINE_SEARCH_H
