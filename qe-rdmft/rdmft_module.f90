!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE rdmft_module
  !------------------------------------------------------------------
  !! Top-level data and control variables for the RDMFT (Reduced
  !! Density Matrix Functional Theory) solver in PWscf.
  !!
  !! The naming conventions for the input keywords and the algorithmic
  !! choices follow the ABACUS reference implementation in
  !! \texttt{abacus\_source/source\_lcao/module\_rdmft}; see
  !! \texttt{PW/src/rdmft/README.md} for a side-by-side mapping.
  !!
  !! The solver runs after the standard KS-SCF (or KS-EXX) loop has
  !! converged, takes the converged Kohn-Sham orbitals and occupations
  !! as the initial guess for the natural orbitals \(C^k\) and natural
  !! occupation numbers \(n_{ik}\in[0,1]\), and minimises the RDMFT
  !! total energy
  !!
  !!  E[n,C] = E_one + E_H + E_xc[gamma] + E_ewald
  !!
  !! where \(E_{xc}\) uses one of the available RDMFT exchange-correlation
  !! functionals listed in \texttt{rdmft\_xc}.  The exchange-like part
  !! of the energy is built using the existing PWscf EXX-via-ACE
  !! machinery (\texttt{aceinit}, \texttt{vexxace\_*}) by replacing the
  !! KS occupation weights with \(g(n_{ik})\,w_k\) for the duration of
  !! one ACE build.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ----------------------------------------------------------------
  ! Master switch (read from the &RDMFT namelist).
  ! ----------------------------------------------------------------
  !
  LOGICAL :: do_rdmft = .FALSE.
  !! Master switch: run the RDMFT alternating optimisation after the
  !! KS-SCF converges.  Set to TRUE in the &RDMFT namelist.
  !
  ! ----------------------------------------------------------------
  ! XC functional and Power-functional exponent.
  ! ----------------------------------------------------------------
  !
  CHARACTER(LEN=16) :: rdmft_functional = 'muller'
  !! 'hf', 'muller', 'power', 'gu', 'chf', 'cga', 'geo', 'hybopt', 'bow', 'bowmod'.
  !! Same accepted spellings as in the ABACUS reader.
  !
  INTEGER :: rdmft_xc_id = 2
  !! Integer code corresponding to \texttt{rdmft\_functional} (set by
  !! \texttt{rdmft\_xc\_set\_type}).
  !
  REAL(DP) :: rdmft_power_alpha = 0.656_DP
  !! Power-functional exponent \(\alpha\); ignored for HF / Muller / GU /
  !! CHF / CGA / GEO / HybOpt (their \(\alpha\) is fixed in
  !! \texttt{rdmft\_xc\_set\_type}).
  !
  REAL(DP) :: rdmft_reg_eps = 1.0e-8_DP
  !! Cutoff \(\varepsilon\) for ABACUS-style piecewise power regularisation
  !! (see \texttt{rdmft\_xc.f90} and ABACUS \texttt{rdmft\_xc\_functional.h}):
  !! \(g(n)=n^\alpha\) for \(n\ge\varepsilon\), linear Taylor below
  !! \(\varepsilon\); \(g'(n)=\alpha\max(n,\varepsilon)^{\alpha-1}\).
  !! Used by Muller / Power / GU / GEO / HybOpt.  No effect for HF
  !! (\(\alpha=1\)).
  !!
  !! Default matches ABACUS (\(10^{-8}\)).  For Muller (\(\alpha=0.5\)):
  !! \(|g'(0)|=0.5/\sqrt{\texttt{rdmft\_reg\_eps}}\).  Larger values (e.g.
  !! \(10^{-3}\)) soften the occupation block; combine with Fermi-window
  !! perturbed init for Müller (see \texttt{PW/src/rdmft/README.md}).
  !
  ! ----------------------------------------------------------------
  ! Solver strategy and constraint method.
  ! ----------------------------------------------------------------
  !
  CHARACTER(LEN=16) :: rdmft_solver_strategy = 'alternating'
  !! 'alternating' (default) optimises occupations with orbitals fixed
  !! and orbitals with occupations fixed in alternation.  'joint'
  !! packs the occupation parameters and orbital coefficients into a
  !! single point on the product manifold and steps both blocks
  !! simultaneously with a single Armijo line search.  Mirror of
  !! ABACUS' \texttt{rdmft\_solver\_strategy}.
  !
  CHARACTER(LEN=16) :: rdmft_block_order = 'occ_orb'
  !! Block order within each outer alternating cycle (ignored for
  !! \texttt{joint}).  \texttt{'occ\_orb'} (default): occupation block
  !! then orbital block; \texttt{'orb\_occ'}: orbitals first, then
  !! occupations.
  !
  CHARACTER(LEN=16) :: rdmft_constraint = 'projected_gradient'
  !! 'augmented_lagrangian', 'projected_gradient'.  Default mirrors
  !! the user request to use projected gradient on occupations.
  !
  CHARACTER(LEN=8) :: rdmft_occ_optimizer = 'spg2'
  !! Occupation optimiser inside the alternating block.  Recognised
  !! values:
  !!
  !! * ``'spg2'`` (default, canonical Birgin--Mart\'inez--Raydan
  !!   Algorithm 2.2): spectral steepest descent in the box-constrained
  !!   occupation space.  At each inner step the direction is
  !!   \(\mathbf{d}_k = P_w(\mathbf{n}_k-\alpha_k^{\mathrm{BB}}\mathbf{g}_k)-\mathbf{n}_k\)
  !!   with \(\alpha_k\) the BMR SPG2 spectral steplength (inverse
  !!   Rayleigh quotient from consecutive accepted iterates; Birgin--
  !!   Mart\'inez--Raydan 1999/2000 Algorithm~2.2 Step~3) and
  !!   \(P_w\) the Euclidean L\(_2\) proximal map onto \(\mathcal{F}\)
  !!   (not \(P_u\); see doc §5.4.1).
  !!   The line search walks the straight chord
  !!   \(\mathbf{n}_k+\lambda\mathbf{d}_k\), \(\lambda\in[0,1]\).
  !!   ``'sd'`` is kept as a synonym so old input files still work.
  !! * ``'cg'``: SPG2 framework with a Polak--Ribi\`ere CG direction
  !!   substituted for the raw gradient.  More robust than canonical
  !!   SPG2 on wide-Hessian metals but adds momentum bookkeeping; not
  !!   the BMR canonical choice.
  !! * ``'lbfgs'``: SPG2 framework with an L-BFGS direction
  !!   (history depth ``rdmft_lbfgs_memory``).
  !! * ``'bgd'``: ELK-style box-aware occupation gradient descent
  !!   (port of \texttt{rdmvaryn.f90}; see :subroutine:`rdmft_bgd_occ_block`).
  !!   Uses the charge-conserving, box-aware reduced gradient
  !!   \(\gamma_{ik}\) with chemical-potential shift \(\kappa\) chosen
  !!   so \(\sum_k w_k\sum_i\gamma_{ik}=0\), followed by an energy
  !!   line search along the (straight, feasible) segment
  !!   \(n_0+\tau\gamma\).
  !! * ``'ebi'``: explicit-by-implicit (EBI) erf parameterisation with
  !!   gradient descent on the unconstrained \(x_{ik}\) parameters
  !!   (EBI@GD per Yao \textit{et al.}, J. Phys. Chem. A
  !!   \textbf{2022}, \textbf{126}, 5654--5662;
  !!   :subroutine:`rdmft_ebi_occ_block`).
  !!
  !! All SPG2 directions (``spg2`` / ``sd`` / ``cg`` / ``lbfgs``)
  !! share the same straight-chord projected line search.
  !
  REAL(DP) :: rdmft_bgd_tau = 1.0_DP
  !! Base occupation step for the ``bgd`` and ``ebi`` occupation optimisers.
  !! For ``bgd``: seeds \(\tau_0=\min(\texttt{rdmft\_bgd\_tau},\tau_{\max})\)
  !! on the ELK feasible segment (ELK ``taurdmn``, default ``1.0``).
  !! For ``ebi``: initial line-search step in parameter space when
  !! ``rdmft_occ_ls_init_step = fixed``.
  !
  REAL(DP) :: rdmft_bgd_backtrack = 0.75_DP
  !! Geometric backtracking factor for the ELK box-aware gradient-descent
  !! occupation block (``rdmft_occ_optimizer = 'bgd'``).  Default
  !! ``0.75`` matches ELK's hard-coded factor in \texttt{rdmvaryn.f90}.
  !! Lives separately from the SPG/EBI ``rdmft_line_search_rho`` so
  !! changing the global backtracking factor does not silently
  !! perturb the ELK-mirror behaviour of the ``bgd`` block.
  !
  CHARACTER(LEN=8) :: rdmft_orb_optimizer = 'cg'
  !! Orbital block optimiser in the alternating strategy.  Values:
  !! ``'sd'`` (steepest descent on the Stiefel tangent),
  !! ``'cg'`` (Polak--Ribi\`ere; default), ``'lbfgs'`` (L-BFGS).
  !! Unrelated to the SPG2 occupation-block default
  !! (``rdmft_occ_optimizer = 'spg2'``).
  !
  INTEGER :: rdmft_lbfgs_memory = 10
  !! L-BFGS history length (number of stored ``(s, y)`` pairs).
  !! Used when \texttt{rdmft\_occ\_optimizer}, \texttt{rdmft\_orb\_optimizer}
  !! or \texttt{rdmft\_joint\_optimizer} is set to ``lbfgs``.
  !
  CHARACTER(LEN=8) :: rdmft_joint_optimizer = 'cg'
  !! 'sd', 'cg' (default), 'lbfgs', 'spg'.  Used by
  !! \texttt{rdmft\_solver\_strategy = 'joint'}.
  !!
  !! For 'sd'/'cg'/'lbfgs' the joint optimiser drives the packed
  !! ``(params, evc\_flat)`` vector (cosine-squared occupation
  !! reparameterisation) with a single search direction and a single
  !! Armijo line search (:subroutine:`rdmft_run_joint`).
  !!
  !! For 'spg' the occupations are optimised in native n-space with the
  !! canonical SPG2 projected-gradient step (Euclidean proximal
  !! projection ``P_w`` + BMR spectral steplength + straight feasible
  !! chord) and the orbitals with a Stiefel Polak--Ribiere CG step,
  !! coupled through one joint Armijo line search
  !! (:subroutine:`rdmft_run_joint_spg`).
  !
  ! ----------------------------------------------------------------
  ! Iteration limits and tolerances.
  ! ----------------------------------------------------------------
  !
  INTEGER :: rdmft_outer_maxiter = 50
  !! Outer alternating cycles.
  !
  INTEGER :: rdmft_occ_maxiter = 20
  !! Inner occupation iterations per outer cycle.
  !
  INTEGER :: rdmft_orb_maxiter = 20
  !! Inner orbital iterations per outer cycle.
  !
  REAL(DP) :: rdmft_energy_tol = 1.0e-7_DP
  !! Outer convergence: \(|E_{out}-E_{out-1}| < \) this (alternating
  !! mode: sufficient on its own once \texttt{outer > 1}).
  !
  REAL(DP) :: rdmft_orb_grad_tol = 1.0e-5_DP
  !! Stop the orbital block when \(\|G_R\|<\) this.
  !
  REAL(DP) :: rdmft_occ_grad_tol = 1.0e-5_DP
  !! Stop the occupation SPG2 block when
  !! \(\|g_1\|_\infty = \|P_w(n-\alpha\,\partial E/\partial n)-n\|_\infty\)
  !! with \(\alpha=1\) (BMR) falls below this.  No chord step or line
  !! search is attempted once the map norm is this small.  Set \(\le 0\)
  !! to disable this criterion.
  !
  REAL(DP) :: rdmft_occ_tol = 1.0e-6_DP
  !! Inner occupation stopping tolerance: exit the OCC block when an
  !! accepted step has \(\sum_{ik}|n_{ik}^{\mathrm{new}}-n_{ik}^{\mathrm{old}}|
  !! <\) this.  Set \(\le 0\) to disable this criterion (the block
  !! then runs until \texttt{rdmft\_occ\_maxiter} or line-search
  !! failure).
  !
  REAL(DP) :: rdmft_line_search_c1 = 1.0e-4_DP
  !! Armijo sufficient-decrease constant (used by all backends).
  !
  REAL(DP) :: rdmft_line_search_c2 = 0.9_DP
  !! Wolfe curvature condition constant: \(|\langle g(\alpha),d\rangle|
  !! \le c_2 |\langle g_0,d\rangle|\).  Only used by the
  !! \texttt{'wolfe'} line-search backend.
  !
  REAL(DP) :: rdmft_line_search_rho = 0.5_DP
  !! Geometric backtracking factor.
  !
  LOGICAL :: rdmft_line_search_polynomial = .TRUE.
  !! Use **quadratic** (\(f_0,\,\phi'(0),\,f(\alpha)\)) interpolation
  !! to pick the next trial step in Armijo backtracking, with
  !! safeguarded clamps to \([\rho_{\min}\alpha,\rho_{\max}\alpha]\)
  !! (\(\rho_{\min}=0.1,\,\rho_{\max}=0.5\)).  When the quadratic fit
  !! is degenerate (small denominator, non-descent slope) the helper
  !! falls back to plain geometric shrinking by
  !! :var:`rdmft_line_search_rho`.  Set to \texttt{.FALSE.} to
  !! recover the legacy behaviour
  !! \(\alpha \leftarrow \alpha \cdot \texttt{rdmft\_line\_search\_rho}\)
  !! at every reject.  Mirrors ABACUS's ``rdmft_line_search_polynomial =
  !! true`` default.
  !
  INTEGER :: rdmft_line_search_max_iter = 20
  !! Maximum line-search trials.
  !
  CHARACTER(LEN=16) :: rdmft_line_search = 'wolfe'
  !! Line-search backend:
  !!
  !! * ``'armijo'`` -- classical monotone Armijo backtracking;
  !! * ``'zhang_hager'`` (default) -- non-monotone Zhang-Hager
  !!   running-average Armijo, more tolerant of curved energy
  !!   landscapes that the orbital block of RDMFT routinely
  !!   produces;
  !! * ``'wolfe'`` -- bracket-and-zoom strong Wolfe conditions
  !!   with cubic interpolation, recommended for L-BFGS.
  !
  REAL(DP) :: rdmft_zhang_hager_eta = 0.85_DP
  !! Memory parameter for the Zhang-Hager moving reference
  !! \(C_k = (\eta Q_{k-1} C_{k-1} + f_k)/Q_k\).  ``eta = 0``
  !! recovers monotone Armijo; ``eta = 1`` keeps the maximum
  !! reference over the entire history.
  !
  REAL(DP) :: rdmft_orb_ls_stepsize = 1.0_DP
  !! Initial Armijo trial step for the orbital block.
  !
  CHARACTER(LEN=16) :: rdmft_orb_strategy = 'joint'
  !! Multi-k orbital descent strategy.  Two values:
  !!
  !! * ``'joint'`` (default) -- Riemannian descent on the **product**
  !!   Stiefel manifold ``\prod_k \mathrm{St}(\text{nbnd}, \text{npwx})``
  !!   with one global Armijo line search per inner iteration: a
  !!   single scalar ``alpha`` retracts every ``C^k`` simultaneously
  !!   and is accepted only when the total energy decreases
  !!   sufficiently.  Most accurate but the per-trial cost scales as
  !!   nks.
  !!
  !! * ``'block_k'`` -- block-coordinate Stiefel descent: sweep
  !!   through k-points (or per-pool k-groups) and at each k take
  !!   ONE Riemannian step, holding ``\{C^{k'\ne k}\}`` fixed.  The
  !!   line search at every k still uses the **global** energy
  !!   (``rdmft_total_energy`` over all k), but only the current
  !!   k's evc is varied per trial.  Each k carries its own
  !!   optimiser state across outer iterations (per-k SD / CG /
  !!   L-BFGS history), giving Gauss-Seidel-style descent.  This is
  !!   the natural mode for QE's k-point pool parallelism: each MPI
  !!   pool optimises its own k-points and only the final
  !!   global-energy reduction needs an inter-pool ``mp_sum``.
  !
  LOGICAL :: rdmft_active = .FALSE.
  !! INTERNAL flag.  Set to ``TRUE`` for the entire duration of
  !! :subroutine:`rdmft_run` and reset to ``FALSE`` on exit.  Used by
  !! :subroutine:`sum_band` to skip the standard ``weights()`` call,
  !! which would otherwise recompute ``wg`` from the band eigenvalues
  !! and the Fermi level (smearing, fixed_occ, tetrahedra, ...).  The
  !! RDMFT solver controls ``wg`` itself via
  !! :subroutine:`rdmft_set_wg_from_n` / :subroutine:`rdmft_set_wg_for_channel`,
  !! so allowing ``weights()`` to run would silently overwrite those
  !! weights and decouple ``rho`` / ``ehart`` from the natural
  !! occupations.  Mirrors the DMFT ``dmft`` / ``dmft_updated`` guard
  !! already present in :subroutine:`sum_band`.
  !
  LOGICAL :: rdmft_exxbuff_stale = .TRUE.
  INTEGER :: rdmft_exxbuff_nbnd_occ = 0
  !! Upper band index ``exxbuff`` was sized for at the last
  !! :subroutine:`rdmft_reinit_exxbuff` / ``exxinit`` call.
  !! When ``x_nbnd_occ`` grows (e.g. TSM DOS probe at \(n=0.5\) on a
  !! nearly-empty band), ``exxbuff`` must be rebuilt before ``vexx``.
  !! INTERNAL flag.  ``exxbuff`` (the real-space EXX wavefunction
  !! buffer used by ACE / vexx) is set up once at the start of
  !! :subroutine:`rdmft_run`, and again every time the orbitals
  !! change.  This flag is set to ``TRUE`` after every accepted
  !! orbital step (so the next energy evaluation will refresh
  !! ``exxbuff`` via :subroutine:`rdmft_refresh_exxbuff`) and to
  !! ``FALSE`` after such a refresh.  The occupation block leaves it
  !! ``FALSE`` (the orbitals do not change there), so the per-channel
  !! cost in the occupation block stays cheap.
  !
  LOGICAL :: rdmft_orb_no_ls = .FALSE.
  !! Skip the Armijo line search in the orbital block and apply a
  !! single fixed-step Stiefel SD update (`alpha = rdmft_orb_ls_stepsize
  !! / max(|G_R|, 1)`).  Cuts the h_psi call count per outer cycle by
  !! ~10x and avoids a heap-corruption pattern in QE's vloc_psi_k_acc
  !! that surfaces after several hundred invocations on multi-k
  !! complex pseudopotential paths.  The solver flips this flag on
  !! automatically when `nks > 1`.
  !
  REAL(DP) :: rdmft_occ_ls_stepsize = 1.0_DP
  !! Initial line-search trial \(\alpha_0\) for the SPG occupation
  !! block when \texttt{rdmft\_occ\_ls\_init\_step = 'fixed'}.
  !
  CHARACTER(LEN=24) :: rdmft_occ_ls_type = 'auto'
  !! SPG occupation line search: ``auto`` (BMR SPG2 default: monotone
  !! Armijo backtracking with ``rdmft_spg_gamma``), ``armijo`` (same
  !! monotone Armijo), ``gll``/``nm_armijo`` (GLL nonmonotone Armijo
  !! with ``rdmft_spg_ls_memory``),
  !! ``strong_wolfe``/``sw``, ``nonmonotone_strong_wolfe``/``nm_sw``
  !! (Zhang--Hager reference), ``weak_wolfe``/``wolfe``.
  !
  CHARACTER(LEN=24) :: rdmft_orb_ls_type = 'auto'
  !! Alternating orbital line search: ``auto`` (monotone strong Wolfe),
  !! ``armijo``, ``strong_wolfe``/``sw``, ``weak_wolfe``/``wolfe``.
  !
  CHARACTER(LEN=16) :: rdmft_occ_ls_init_step = 'bb'
  !! Initial \(\alpha_0\) policy: ``fixed``, ``barzilai_borwein``/``bb``,
  !! ``quadratic``.
  !
  CHARACTER(LEN=16) :: rdmft_orb_ls_init_step = 'bb'
  !! Initial \(\alpha_0\) policy for the orbital line search.  Same
  !! values as :var:`rdmft_occ_ls_init_step`: ``fixed`` (use the
  !! literal :var:`rdmft_orb_ls_stepsize`), ``barzilai_borwein``/``bb``
  !! (suggest :math:`\alpha_0 = (s^T s)/(s^T y)` from the previous
  !! orbital iterate, fall back to :var:`rdmft_orb_ls_stepsize` on the
  !! very first inner step where no history is available),
  !! ``quadratic`` (quadratic-fit suggestion from the previous
  !! line-search history).  The BB / quadratic suggestions are formed
  !! on the **Euclidean** packed (Re/Im, per-k flattened) view of the
  !! orbital iterate -- standard practice for Stiefel optimisers
  !! coupled with Stiefel retraction; see ABACUS \texttt{joint\_bb} for
  !! the same pattern.
  !
  INTEGER :: rdmft_line_search_max_zoom = 30
  !! Zoom iteration cap in Wolfe line searches.
  !
  REAL(DP) :: rdmft_bb_alpha_min = 1.0e-8_DP
  REAL(DP) :: rdmft_bb_alpha_max = 1.0e2_DP
  !! Barzilai--Borwein / quadratic \(\alpha_0\) clamp for BGD / EBI /
  !! orbital line searches (not used by the SPG occupation block; see
  !! :var:`rdmft_spg_alpha_min` / :var:`rdmft_spg_alpha_max`).
  !
  ! ----------------------------------------------------------------
  ! BMR SPG2 defaults (Birgin--Mart\'inez--Raydan 1999/2000).
  ! ----------------------------------------------------------------
  !
  REAL(DP) :: rdmft_spg_gamma = 1.0e-4_DP
  !! Sufficient-decrease constant \(\gamma\) in the SPG Armijo test
  !! (BMR numerical experiments).
  !
  REAL(DP) :: rdmft_spg_alpha_min = 1.0e-30_DP
  REAL(DP) :: rdmft_spg_alpha_max = 1.0e30_DP
  !! Spectral-steplength safeguard ``[\alpha_{\min}, \alpha_{\max}]``
  !! for BMR SPG2 Step~3.
  !
  REAL(DP) :: rdmft_spg_sigma1 = 0.1_DP
  REAL(DP) :: rdmft_spg_sigma2 = 0.9_DP
  !! Chord backtracking safeguards ``\sigma_1, \sigma_2`` for the
  !! one-dimensional quadratic fit on \(\lambda\) (BMR eq.~(2)).
  !
  INTEGER :: rdmft_spg_ls_memory = 10
  !! GLL nonmonotone memory \(M\) for the SPG2 chord line search
  !! (BMR Algorithm~2.2, Grippo--Lampariello--Lucidi reference).
  !! Set \(M=1\) to recover monotone Armijo against \(f(x_k)\).
  !
  ! ----------------------------------------------------------------
  ! Preconditioning controls.
  ! ----------------------------------------------------------------
  !
  LOGICAL :: rdmft_occ_precond = .FALSE.
  !! Enable ELK-style diagonal scaling of the SPG occupation
  !! direction (\texttt{rdmft\_compute\_elk\_occ\_scale}).  Default
  !! \texttt{.FALSE.}: use raw \(\partial E/\partial n_{ik}\) from
  !! \texttt{rdmft\_grad\_n}.
  !
  REAL(DP) :: rdmft_occ_precond_shift = 1.0_DP
  !! Level shift (Ry) for legacy \texttt{rdmft\_apply\_occ\_precond}
  !! (band-diagonal); unused by the SPG occupation block.
  !
  LOGICAL :: rdmft_orb_precond = .FALSE.
  !! Apply a level-shift rotation preconditioner to the orbital
  !! gradient before building the Stiefel CG / L-BFGS direction.
  !! The preconditioner is
  !!
  !!     \(\tilde G_i(\mathbf g) = G_i(\mathbf g) /
  !!         (\tfrac12 (k+\mathbf g)^2 + V_0 - \epsilon_i + \mathrm{shift})\),
  !!
  !! mirroring QE's ``g_psi`` band-by-band preconditioner.  Without
  !! it, plane-wave components with very large kinetic energy
  !! dominate the Stiefel descent direction and force the orbital
  !! block to take very small steps.
  !
  REAL(DP) :: rdmft_orb_precond_shift = 1.0_DP
  !! Level shift (in Ry) added to the orbital preconditioner
  !! denominator to keep it strictly positive.  Increasing this
  !! value softens the preconditioner toward the identity.
  !
  ! ----------------------------------------------------------------
  ! Initial occupation seed (matches ABACUS rdmft_occ_init_*).
  ! ----------------------------------------------------------------
  !
  CHARACTER(LEN=16) :: rdmft_occ_init_mode = 'ks'
  !! 'ks', 'perturbed', 'binary', or 'uniform' (ABACUS
  !! \texttt{OccInitMode}).
  !
  REAL(DP) :: rdmft_occ_init_perturb = 1.0e-2_DP
  !! Perturbation magnitude ``delta`` for \texttt{perturbed} and
  !! \texttt{binary} modes.
  !
  INTEGER :: rdmft_occ_init_nbands_top = 0
  !! Fermi-window half-width $K$ for \texttt{perturbed},
  !! \texttt{binary}, and \texttt{uniform} (ABACUS convention):
  !! when $K>0$, only the $K$ bands above and $K$ bands below the
  !! approximate Fermi boundary are modified; \texttt{binary} and
  !! \texttt{uniform} leave KS occupations unchanged when $K\le 0$
  !! (and \texttt{binary} also requires $\delta>0$).  For
  !! \texttt{perturbed} with $K\le 0$ the legacy QE recipe perturbs
  !! **every** band instead.
  !
  ! ----------------------------------------------------------------
  ! Verbosity / debugging.
  ! ----------------------------------------------------------------
  !
  INTEGER :: rdmft_verbose = 1
  !! 0 = silent; 1 = one line per outer/inner iteration plus one Wolfe
  !! line-search summary per call; 2 = adds OCC/ORB block wall times,
  !! per-macro-iteration wall time, optimisation total, and per-trial
  !! Wolfe/Armijo bracket/zoom lines.
  !
  LOGICAL :: rdmft_grad_check = .FALSE.
  !! If TRUE, run ABACUS-style finite-difference gradient checks on
  !! occupations and orbitals before the main optimisation loop.
  !!
  ! NOTE: the occupation coupling fast path (formerly
  ! ``rdmft_occ_reuse_orbitals`` + ``rdmft_frozen_fock``) has been
  ! removed.  It cached per-band exchange couplings and reused them
  ! across the occupation inner loop, but was not k-point-pool
  ! reproducible on metallic / fractional-occupation systems (it
  ! introduced an ``-nk``-dependent energy spread of order 1e-4 Ry that
  ! the SPG loop amplified to the mRy level).  The solver now always
  ! rebuilds the exchange via the per-trial ACE path, which is
  ! pool-invariant.
  ! ----------------------------------------------------------------
  ! Post-solve DOS / magnetization (rdmft_dos.x only; not pw.x).
  ! ----------------------------------------------------------------
  !
  LOGICAL :: rdmft_compute_dos = .FALSE.
  !! Internal flag set by \texttt{rdmft\_dos.x} only.  Must remain
  !! \texttt{.FALSE.} during \texttt{pw.x}; use the dedicated
  !! \texttt{rdmft\_dos.x} post-processor for DOS and Mulliken
  !! magnetization after the RDMFT save is written.
  !!
  REAL(DP) :: rdmft_dos_emin = -1.0e6_DP
  REAL(DP) :: rdmft_dos_emax =  1.0e6_DP
  !! Energy window for DOS plots in **Ry**; leave at the sentinel
  !! $\pm 10^6$ to auto-detect from band extrema.
  !!
  REAL(DP) :: rdmft_dos_deltae = 1.0e-3_DP
  !! Energy grid step (Ry).
  !!
  REAL(DP) :: rdmft_dos_degauss = 1.0e-3_DP
  !! Gaussian broadening width (Ry).
  !!
  LOGICAL :: rdmft_dos_pdos = .TRUE.
  !! Write $l$-resolved partial DOS per atom site.
  !!
  LOGICAL :: rdmft_dos_occ_weighted = .FALSE.
  !! In \texttt{rdmft\_dos\_spectral='elk'} (default):
  !! \texttt{.FALSE.} weights each state by
  !! \texttt{sc}$\times$\texttt{occmax}; \texttt{.TRUE.} uses
  !! \texttt{sc}$\times$\texttt{occsv}=\texttt{sc}$\times$\(n\times\)
  !! \texttt{occmax}.  In \texttt{'sharma'} mode this flag selects
  !! occupied-only vs.\ two-branch PRL Eq.~(7) spectral sum instead.
  !! The per-atom Mulliken magnetisation always uses the occupation,
  !! independently of this flag.
  !!
  CHARACTER(LEN=8) :: rdmft_dos_spectral = 'elk'
  !! Spectral assembly model.  \texttt{'elk'} (default) mirrors ELK
  !! \texttt{dos.f90}: one \(\delta(\omega-\varepsilon)\) per state
  !! with \texttt{sc}$\times$\texttt{occmax} (or \texttt{sc}$\times$\(n\)
  !! \texttt{occmax} when \texttt{rdmft\_dos\_occ\_weighted}).\ 
  !! \texttt{'sharma'} uses PRL 110, 116403 Eq.~(7) two-branch sum.
  !!
  CHARACTER(LEN=8) :: rdmft_dos_eref = 'auto'
  !! Energy reference for DOS plots (Ry).  \texttt{'auto'} (default):
  !! \texttt{RDM\_DEDN.OUT} chemical potential
  !! \(\mu=\texttt{dedn}(i_\mu)\) for both ELK and Sharma spectral modes.
  !! \texttt{'efermi'} and \texttt{'mu'} are aliases for the same
  !! \texttt{dedn}-table recipe; \texttt{'efermig'} uses legacy
  !! \texttt{efermig} bisection on TSM energies instead.
  !!
  CHARACTER(LEN=256) :: rdmft_dos_prefix = ' '
  !! Output file prefix; blank means \texttt{trim(prefix)//'.rdmft'}.
  !!
  CHARACTER(LEN=8) :: rdmft_dos_integration = 'brzint'
  !! BZ integration method.  \texttt{'brzint'} (default) uses ELK's
  !! trilinear BZ interpolation + histogram integration
  !! (\texttt{brzint.f90}).  \texttt{'gauss'} uses Gaussian broadening
  !! on the SCF k-mesh (QE \texttt{dos.x} style; \texttt{'sharma'}
  !! mode only).
  !!
  INTEGER :: rdmft_dos_ngrkf = 100
  !! BZ subdivision factor for \texttt{brzint}
  !! (\texttt{nsk(i)=max(ngrkf/ngridk(i),1)}).
  !!
  INTEGER :: rdmft_dos_nswplot = 0
  !! Number of 3-point running-average smoothing passes applied after
  !! \texttt{brzint} (\texttt{fsmooth}).
  !!
  LOGICAL :: rdmft_dos_msum = .TRUE.
  !! When \texttt{.TRUE.} (default) partial DOS is summed over magnetic
  !! quantum number $m$ and written per $l$.
  !!
  LOGICAL :: rdmft_dos_ssum = .FALSE.
  !! When \texttt{.TRUE.} partial and total DOS are summed over spin
  !! and written as a single column.
  !!
  REAL(DP) :: rdmft_dos_sqaxis(3) = (/ 0.0_DP, 0.0_DP, 1.0_DP /)
  !! Spin quantisation axis for noncollinear DOS weights; used to rotate
  !! the spin-density matrix before taking the diagonal weights
  !! \texttt{sc}.
  !!
  INTEGER :: rdmft_dos_nwplot = 0
  !! Energy-mesh point count.  When $\ge 2$, the DOS grid uses this
  !! many points between \texttt{rdmft\_dos\_emin} and
  !! \texttt{rdmft\_dos\_emax} instead of \texttt{rdmft\_dos\_deltae}.
  !! $0$ = use \texttt{rdmft\_dos\_deltae}.
  !!
  ! ----------------------------------------------------------------
  ! Post-solve EMD / Compton (rdmft_emd.x only; not pw.x).
  ! ----------------------------------------------------------------
  !
  LOGICAL :: rdmft_compute_emd = .FALSE.
  !! Internal flag set by \texttt{rdmft\_emd.x} only.
  !!
  LOGICAL :: rdmft_emd_plot1d = .FALSE.
  !! When \texttt{.TRUE.}, compute a twice-integrated 1D Compton profile
  !! after writing the momentum density (ELK task 171).
  !!
  INTEGER :: rdmft_emd_line_npt = 150
  !! Number of points along the 1D plotting line.
  !!
  REAL(DP) :: rdmft_emd_line_start(3) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
  REAL(DP) :: rdmft_emd_line_end(3)   = (/ 1.0_DP, 0.0_DP, 0.0_DP /)
  !! Plotting-line endpoints in crystal reciprocal coordinates (ELK
  !! \texttt{plot1d}).
  !!
  CHARACTER(LEN=256) :: rdmft_emd_prefix = ' '
  !! Output file prefix; blank means \texttt{trim(prefix)//'.rdmft'}.
  !!
  LOGICAL :: rdmft_emd_write_elk_fmt = .FALSE.
  !! Also write unformatted ELK-style \texttt{EMD.OUT}.
  !!
  LOGICAL :: rdmft_emd_spin_sum = .TRUE.
  !! When \texttt{.TRUE.} (default), sum spin channels into a single
  !! $n(p)$; when \texttt{.FALSE.}, write/store spin-resolved columns.
  !!
  INTEGER :: rdmft_emd_band_min = 1
  INTEGER :: rdmft_emd_band_max = 0
  !! Natural-orbital band window included in the summed EMD/Compton
  !! output.  ``band_max <= 0`` means through \texttt{nbnd}.  Use this
  !! to exclude semicore bands from pseudopotentials with many valence
  !! electrons, e.g. SG15 Na where band 5 carries the conduction electron.
  !!
  LOGICAL :: rdmft_emd_write_band_resolved = .FALSE.
  !! Also write a sparse band-resolved ``*.emd_band`` file with rows
  !! ``ik ib ig p1 p2 p3 |p| n_ib(p)`` (plus spin columns when requested).
  !!
  REAL(DP) :: rdmft_emd_sqaxis(3) = (/ 0.0_DP, 0.0_DP, 1.0_DP /)
  !! Spin quantisation axis for noncollinear spin-resolved EMD (same
  !! convention as \texttt{rdmft\_dos\_sqaxis}).
  !!
  REAL(DP) :: rdmft_emd_pmax_kf_factor = 5.0_DP
  !! Default momentum range for EMD / Compton integration:
  !! $|p|_{\max} = \texttt{pmax\_kf\_factor}\,k_F$ using the
  !! selected-band equivalent Fermi momentum.  Set $\le 0$ to use the
  !! plane-wave cutoff from \texttt{ecutwfc} instead.
  !!
  LOGICAL :: rdmft_emd_do_compton = .FALSE.
  !! When \texttt{.TRUE.}, compute directional Compton profile $J(q)$
  !! on a uniform $q$-grid (Bohr$^{-1}$) and write \texttt{*.emd\_compton}
  !! plus per-band profiles \texttt{*.emd\_compton\_band}.
  !!
  INTEGER :: rdmft_emd_cp_nq = 256
  !! Number of points on the uniform Compton $q$-grid ($\ge 64$).
  !!
  REAL(DP) :: rdmft_emd_cp_qmax = 0.0_DP
  !! Upper $q$ bound (Bohr$^{-1}$); $\le 0$ uses $|p|_{\max}$ from EMD.
  !!
  REAL(DP) :: rdmft_emd_cp_qmin = 0.0_DP
  !! Lower $q$ bound (Bohr$^{-1}$); default is $0$.  Only $q \ge 0$ is used.
  !!
  LOGICAL :: rdmft_emd_do_autocorr = .FALSE.
  !! When \texttt{.TRUE.}, compute directional autocorrelation $B(r)$
  !! from the Compton profile and write \texttt{*.emd\_af}.
  !!
  INTEGER :: rdmft_emd_af_npt = 150
  !! Number of $r$ points for the autocorrelation output ($\ge 2$).
  !!
  REAL(DP) :: rdmft_emd_af_r_start = 0.0_DP
  REAL(DP) :: rdmft_emd_af_r_end = 6.0_DP
  !! Real-space range (Bohr) along the direction conjugate to
  !! \texttt{line\_start}$\to$\texttt{line\_end}.
  !!
  ! ----------------------------------------------------------------
  ! Denser-k restart refinement (pw.x only).
  ! ----------------------------------------------------------------
  !!
  CHARACTER(LEN=256) :: rdmft_source_prefix = ' '
  !! Coarse RDMFT save prefix.  Together with
  !! \texttt{rdmft\_source\_outdir}, triggers automatic k-grid
  !! interpolation when the current run uses a denser Monkhorst--Pack
  !! mesh (e.g.\ 2$\times$2$\times$2 $\to$ 3$\times$3$\times$3).
  !! Must differ from the current run
  !! \texttt{prefix}.
  !!
  CHARACTER(LEN=256) :: rdmft_source_outdir = ' '
  !! Coarse RDMFT save directory (paired with
  !! \texttt{rdmft\_source\_prefix}).  Must differ from the current
  !! run \texttt{outdir} so the fine job does not overwrite the coarse
  !! \texttt{.save/} tree and \texttt{.rdmft.save} file.
  !!
  LOGICAL :: rdmft_do_krefine = .FALSE.
  !! Set by \texttt{rdmft\_resolve\_krefine\_from\_source} when source
  !! paths are given and the fine k-mesh is denser than the coarse mesh
  !! read from the source save (each \texttt{nk} not smaller; integer
  !! commensurability not required).
  !!
  LOGICAL :: rdmft_krefine_resolved = .FALSE.
  !! Guards one-shot k-mesh resolution per \texttt{pw.x} invocation.
  !!
  INTEGER :: rdmft_krefine_orb_maxiter = 5
  !! Fixed-occupation orbital polish iterations after k-interpolation.
  !!
  REAL(DP) :: rdmft_krefine_orb_tol = 1.0e-5_DP
  !! Stop orbital polish when the max Stiefel gradient norm falls below
  !! this threshold.
  !!
  LOGICAL :: rdmft_krefine_done = .FALSE.
  !! Set by \texttt{rdmft\_krefine\_for\_restart} when the in-memory
  !! state already holds the fine-grid \((n, evc)\).
  !!
  ! ----------------------------------------------------------------
  ! pp.x density / 1-RDM coherency post-processing (plot_num 126-128).
  ! ----------------------------------------------------------------
  !
  REAL(DP) :: rdmft_xhole_origin(3) = (/ -1.0_DP, -1.0_DP, -1.0_DP /)
  !! Fractional-crystal coordinates of the reference point
  !! \(\mathbf r_0\) for the 1-RDM coherency slice \(h(\mathbf r_0,\mathbf r)\)
  !! \(g(\mathbf r_0,\mathbf r)\) (\texttt{plot\_num=127} in
  !! \texttt{pp.x}).  All components $\le -1$ select the grid point
  !! of maximum \(\rho(\mathbf r)\) automatically.
  !!
  LOGICAL :: rdmft_update_save = .FALSE.
  !! When \texttt{.TRUE.} and \texttt{pp.x} rebuilds the RDMFT density
  !! (\texttt{plot\_num=126--128}), refresh the QE save-directory
  !! charge density via \texttt{write\_scf}.
  !!
  ! ----------------------------------------------------------------
  ! Restart / checkpointing.
  ! ----------------------------------------------------------------
  !
  LOGICAL :: rdmft_restart = .FALSE.
  !! Read the natural occupations (and rely on QE's wavefunction
  !! restart mechanism for the natural orbitals) from a previous
  !! RDMFT run.  When \texttt{.TRUE.}, \texttt{rdmft\_initial\_n\_from\_ks}
  !! reads \texttt{rdmft\_n} from the file
  !! ``{outdir}/{prefix}.rdmft.save`` instead of seeding from the
  !! Kohn-Sham weights.  The orbitals are loaded automatically by
  !! PWscf if the user sets \texttt{startingwfc = 'file'} in the
  !! \texttt{\&electrons} namelist and reuses the same ``outdir`` (and
  !! \texttt{prefix}); the wavefunction buffer is then populated with
  !! the natural orbitals saved by the previous run.
  !
  INTEGER :: rdmft_save_every = 0
  !! Save the RDMFT state (natural occupations and natural orbitals)
  !! every N outer alternating / joint cycles.  Default 1 (save after
  !! every macro iteration).  Set to a larger value (e.g. 5) to reduce
  !! the I/O overhead on long runs; set $\le 0$ to disable periodic
  !! saving entirely (a final save still happens at the end of
  !! :subroutine:`rdmft_run`).
  !
  CHARACTER(LEN=256) :: rdmft_restart_file = ' '
  !! Pathless basename for the restart file (default
  !! ``{prefix}.rdmft.save`` written under ``outdir``).  When non-blank,
  !! the literal filename is used as-is (no ``outdir`` prefix added).
  !
  !
  ! ----------------------------------------------------------------
  ! Internal state.
  ! ----------------------------------------------------------------
  !
  REAL(DP), ALLOCATABLE :: rdmft_n(:,:)
  !! Natural occupation numbers \(n_{ik}\), shape (nbnd, nks).
  !! Each entry is in [0, 1] with the closed-shell convention
  !! \(\sum_{ik} w_k n_{ik} = N_e/spin\_factor\).
  !
  REAL(DP) :: rdmft_etot = 0.0_DP
  REAL(DP) :: rdmft_e_one = 0.0_DP
  REAL(DP) :: rdmft_e_har = 0.0_DP
  REAL(DP) :: rdmft_e_xc = 0.0_DP
  REAL(DP) :: rdmft_e_const = 0.0_DP
  REAL(DP) :: rdmft_e_entropy = 0.0_DP
  !! Energy decomposition reported by the solver.
  !
  REAL(DP) :: rdmft_temp = 0.0_DP
  !! HF-only binary-entropy regularisation temperature in **Kelvin**.
  !! \(E_{\mathrm{ent}}=k_B T\sum_k w_k\sum_i f_{\mathrm{bin}}(n_{ik})\)
  !! with \(k_B=\) \texttt{K\_BOLTZMANN\_RY} (Ry/K).  Default \(0\) disables.
  !
  REAL(DP) :: rdmft_n_target = 0.0_DP
  !! Total-electron equality target
  !! \(\sum_{ik} w_k n_{ik} = \texttt{nelec}\) when magnetisation is
  !! not fixed.
  !
  LOGICAL :: rdmft_fix_magnetization = .FALSE.
  !! If \texttt{.TRUE.} (LSDA and \texttt{tot\_magnetization} set in
  !! the PW input so \texttt{two\_fermi\_energies}=\texttt{.TRUE.}),
  !! occupations are projected separately onto
  !! \(\sum_{k,\uparrow} w_k\sum_i n_{ik}=\texttt{nelup}\) and
  !! \(\sum_{k,\downarrow} w_k\sum_i n_{ik}=\texttt{neldw}\).
  !!
  REAL(DP) :: rdmft_n_target_up = 0.0_DP
  REAL(DP) :: rdmft_n_target_down = 0.0_DP
  !! Spin-resolved targets when \texttt{rdmft\_fix\_magnetization}.
  !
  REAL(DP) :: rdmft_spin_factor = 2.0_DP
  !! 2 for \texttt{nspin}=1 closed-shell, 1 for spin-polarised
  !! per-channel.
  !!
  !! ---------------------------------------------------------------
  !! ACE-projector ``xi`` cache (single-slot).  Set when
  !! :subroutine:`rdmft_compute_xc_channel` finishes a full all-k
  !! :subroutine:`aceinit` call; the next ``rdmft_compute_xc_channel``
  !! call can then skip the (very expensive) all-k aceinit if the
  !! channel ID is the same and the natural occupations have not
  !! changed since.  Invalidated whenever the orbitals change (via
  !! :subroutine:`rdmft_mark_orbitals_changed`) or whenever the
  !! single-k ``aceinit_k`` path runs (which leaves the all-k ``xi``
  !! buffer inconsistent across k-points).
  !!
  !! The freshness check on the natural-occupation side is done by
  !! comparing ``rdmft_n`` against ``rdmft_xi_cache_n`` element-wise;
  !! this is robust against the many sites that update
  !! ``rdmft_n`` (joint optimiser, line search trial, occupation
  !! solvers, gradient checks, SPG sweeps, ...) without having to
  !! sprinkle explicit ``rdmft_invalidate_xi_cache`` calls in all of
  !! them.
  LOGICAL :: rdmft_xi_cache_valid = .FALSE.
  INTEGER :: rdmft_xi_cache_channel = 0
  REAL(DP), ALLOCATABLE :: rdmft_xi_cache_n(:,:)
  !!
  PUBLIC :: rdmft_mark_orbitals_changed, rdmft_invalidate_xi_cache
  !!
CONTAINS
  !
  SUBROUTINE rdmft_mark_orbitals_changed()
    !! Orbitals moved: refresh EXX buffers and drop the xi cache.
    rdmft_exxbuff_stale = .TRUE.
    rdmft_exxbuff_nbnd_occ = 0
    ! Orbital change invalidates exxbuff -> xi (which is built from
    ! exxbuff) is also stale.
    rdmft_xi_cache_valid = .FALSE.
    rdmft_xi_cache_channel = 0
  END SUBROUTINE rdmft_mark_orbitals_changed
  !
  SUBROUTINE rdmft_invalidate_xi_cache()
    !! Drop the ACE-projector cache.  Call from any path that builds
    !! ``xi(:, :, only_k)`` for a single k or otherwise leaves
    !! ``xi`` in a state where it is no longer the all-k aceinit
    !! result for the current channel.  The ``rdmft_n`` snapshot is
    !! kept around to avoid an extra alloc/dealloc cycle; only the
    !! validity flag and channel tag are cleared.
    rdmft_xi_cache_valid = .FALSE.
    rdmft_xi_cache_channel = 0
  END SUBROUTINE rdmft_invalidate_xi_cache
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_allocate(nbnd_in, nks_in)
    !---------------------------------------------------------------
    !! Allocate the per-(band, k-point) state arrays.
    !
    INTEGER, INTENT(IN) :: nbnd_in, nks_in
    !
    IF (.NOT. ALLOCATED(rdmft_n))  ALLOCATE(rdmft_n (nbnd_in, nks_in))
    rdmft_n  = 0.0_DP
    !
  END SUBROUTINE rdmft_allocate
  !
  !-----------------------------------------------------------------
  SUBROUTINE rdmft_deallocate()
    !---------------------------------------------------------------
    !! Release the per-(band, k-point) state arrays.
    !
    IF (ALLOCATED(rdmft_n))  DEALLOCATE(rdmft_n)
    IF (ALLOCATED(rdmft_xi_cache_n)) DEALLOCATE(rdmft_xi_cache_n)
    rdmft_xi_cache_valid = .FALSE.
    rdmft_xi_cache_channel = 0
    !
  END SUBROUTINE rdmft_deallocate
  !
END MODULE rdmft_module
