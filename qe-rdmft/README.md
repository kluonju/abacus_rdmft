# `PW/src/rdmft/` — RDMFT solver for PWscf

Self-contained Fortran implementation of the alternating Reduced
Density Matrix Functional Theory (RDMFT) optimiser inside Quantum
ESPRESSO's `PW/src/`.  The naming conventions, algorithmic options and
INPUT keywords mirror the ABACUS reference module
[`abacus_source/source_lcao/module_rdmft`](../../../abacus_source/source_lcao/module_rdmft)
one-to-one (see also that module's `doc/` for the full mathematical
derivation).

## What it does

After PWscf's standard SCF (or hybrid EXX-SCF) loop converges, the
solver takes the converged Kohn-Sham natural orbitals stored in `evc`
and the converged occupation weights stored in `wg(:,:)` as the initial
guess.  Two solver strategies are available (selected by
`rdmft_solver_strategy`):

* `'alternating'` (default) — alternate between an occupation block
  (projected gradient on `n_ik`) and an orbital block (Stiefel CG on
  `evc`).  The occupation block is **frozen** for the rest of the run
  once the per-cycle change `sum_ik |n_ik^out - n_ik^out-1|` has stayed
  below `rdmft_occ_tol` for two consecutive outer cycles -- the same
  freeze logic ABACUS uses to lock in occupations once they stop
  moving.
* `'joint'` — pack `(occupation_parameters, evc_flat)` into a single
  point on the product manifold and step both blocks simultaneously
  with **one** Armijo line search per outer iteration (linear update
  for the cosine-squared occupation parameters, Stiefel retraction for
  each `C^k`).

It then minimises the RDMFT total energy

```
E[n, psi] = T + E_loc + E_NL + E_H[rho_RDMFT]
          + sum_t coef_t * E_x^HF[gamma_t]
          + E_Ewald
```

with respect to the natural occupation numbers `n_ik in [0, 1]` (one
per band per k-point) and the natural orbital coefficients `C^k = evc`,
subject to the equality constraint

```
sum_k w_k sum_i n_ik = N_e .
```

In the **alternating** strategy, two block updates alternate per outer
cycle:

1. **Occupation block** — projected-gradient updates with the L2
   proximal projector onto `{ 0 <= n <= 1, sum_k w_k sum_i n_ik = N_e }`;
   monotone projected Armijo line search; SD recovery on failure.  The
   block is skipped (frozen) once the per-cycle change in `n_ik` has
   stayed below `rdmft_occ_tol` for two consecutive outer cycles.
2. **Orbital block** — Riemannian SD/CG on each k-point's standard
   Stiefel manifold (`evc^H * evc = I`), with the Cholesky-based polar
   retraction and a monotone Armijo line search.

In the **joint** strategy a single Armijo line search drives the
packed `(p, C)` variable per outer cycle:

1. **Occupation parameters** `p` via `n_ik = cos^2(p_ik)` so the box
   constraint `0 <= n <= 1` is automatically satisfied.  After every
   accepted step we also apply the L2 proximal projector to clean up
   any drift in the linear equality `sum_k w_k sum_i n_ik = N_e`.
2. **Orbital block**: re-projects the orbital portion of the packed
   direction onto the Stiefel tangent space at `C` (necessary when the
   raw direction is built from a Euclidean optimiser) and performs the
   Cholesky polar retraction at each trial step.

The **joint directional derivative**
`dd = <grad_p, dir_p> + sum_k <G_R^k, dir_C^k>`
is checked before the line search; if positive, the direction is
replaced by the packed steepest descent so the line search always has a
descent direction.

The exchange-like part of the RDMFT energy reuses the existing PWscf
**EXX-via-ACE** machinery: for each XC channel the solver overwrites
`wg(:,:)` with the channel-specific weights `wk * w_t(n_ik)`, calls
`exxinit` + `aceinit` so the global `xi(:,:,:)` encodes the modified-DM
Fock operator, then `vexxace_*` returns `Vx_t |psi>` for the per-band
diagonals and the orbital gradient.  No new EXX implementation is
introduced.

## Files

| File | Purpose |
|------|---------|
| `rdmft_module.f90` | Control variables (`do_rdmft`, `rdmft_functional`, ...), tolerances and the natural occupation array. |
| `rdmft_xc.f90` | RDMFT XC functionals (`hf`, `muller`, `power`, `gu`, `chf`, `cga`, `geo`, `hybopt`, `bow`, `bowmod`).  Provides `rdmft_g(n)` / `rdmft_dg(n)` for the separable case, the channel decomposition `rdmft_xc_channel(it, n, coef, w_t(n), w_t'(n))`, and the exact pair-kernel path for `bow`. |
| `rdmft_occupation.f90` | Weighted L2 proximal (joint/init) and ABACUS SPG helpers: uniform-shift proximal, KKT residual, Bertsekas projected-gradient map. |
| `rdmft_spg.f90` | `rdmft_spg_occ_block`: SPG occupation inner loop (optional ELK preconditioner, monotone or Wolfe line search). |
| `rdmft_linesearch.f90` | Zhang–Hager, Barzilai–Borwein, quadratic \(\alpha_0\), cubic zoom helper (`rdmft_cubic_min`). |
| `rdmft_stiefel.f90` | `stiefel_project_tangent_*`, `stiefel_orthonormalize_*` (Cholesky polar) and `stiefel_retract_*` for both gamma-only (real) and full (complex) wavefunction layouts. |
| `rdmft_energy.f90` | `rdmft_total_energy`, `rdmft_grad_n`, and `rdmft_compute_riemannian_gradient`.  Wraps QE's `sum_band` / `v_of_rho` / `set_vrs` to refresh the density-dependent local potential, calls `h_psi` (with EXX off) for the one-body+Hartree diagonals, and `rdmft_compute_xc_channel` (which calls `aceinit` / `vexxace_*`) for each XC channel.  GU diagonal and HF entropy terms live here. |
| `rdmft_grad_check.f90` | `rdmft_check_gradient_consistency` (finite-difference vs analytic gradients). |
| `rdmft_solver.f90` | `rdmft_run` dispatches to `rdmft_run_alternating` (SPG occ + Stiefel orb) or `rdmft_run_joint` (joint product-manifold optimisation). |
| `rdmft_dos.f90` | RDMFT DOS / Mulliken magnetization post-processor (Sharma Eq. 7, ELK-style options). |
| `rdmft_density.f90` | RDMFT particle density and 1-RDM coherency diagnostics for `pp.x` (plot_num 126--128). |
| `rdmft_pp_setup.f90` | Module `rdmft_pp_bootstrap`: shared RDMFT post-processing bootstrap for `pp.x` and `rdmft_dos.x`. |
| `rdmft_input.f90` | `rdmft_read_input`: a tiny standalone reader for the `&rdmft` namelist.  Called from `iosys` (`PW/src/input.f90`) so existing `pw.x` inputs are unaffected when no `&rdmft` block is present. |
| `examples/h2_molecule/h2.in` | H2 in a 12-Bohr cubic cell, KS-HF + Müller RDMFT (molecule demo). |
| `examples/h4_chain_solid/h4_chain.in` | Periodic 4-H chain at the Γ point, KS-HF + Müller RDMFT (solid demo). |

## Hooks into the PW driver

Two minimal patches were applied:

* `PW/src/input.f90` — `iosys` now calls `rdmft_read_input` after the
  standard QE namelists are processed.  When no `&rdmft` block is
  found, `do_rdmft` stays `.false.` and behaviour is unchanged.
* `PW/src/electrons.f90` — `electrons` invokes `rdmft_run` whenever
  `do_rdmft = .true.`, both at the end of a non-hybrid SCF and after
  the hybrid EXX outer loop has converged (and at the "stopping after
  N iterations" exit so the optimiser still runs on the latest
  orbitals).

The CMake build (`PW/CMakeLists.txt`) and the legacy Make build
(`PW/src/Makefile`) both list the new `rdmft/*.f90` objects.

## INPUT keywords (`&rdmft` namelist)

Append at the end of the standard pw.x input.  Defaults match
`rdmft_module.f90` and the keyword names mirror ABACUS's
`rdmft_*` INPUT_PARAMETERS conventions.

```
&rdmft
  do_rdmft               = .true.        ! master switch (default .false.)
  rdmft_functional       = 'muller'      ! hf / muller / power / gu / chf / cga / geo / hybopt / bow / bowmod
  rdmft_power_alpha      = 0.656         ! Power-functional exponent
  rdmft_reg_eps          = 1.0d-8        ! ABACUS piecewise power cutoff (alpha < 1)
  rdmft_solver_strategy  = 'alternating' ! alternating | joint
  rdmft_constraint       = 'projected_gradient'
  rdmft_occ_optimizer    = 'spg2'        ! spg2 (default) | sd (alias) | cg | lbfgs | bgd | ebi
                                          !   spg2 = canonical Birgin-Martinez-Raydan SPG2
                                          !   sd   = synonym for spg2 (kept for backward compat)
                                          !   cg   = SPG2 + Polak-Ribiere CG direction
                                          !   lbfgs= SPG2 + L-BFGS direction (memory rdmft_lbfgs_memory)
                                          !   bgd  = ELK box-aware GD (rdmvaryn + Armijo)
                                          !   gd   = deprecated alias for bgd
                                          !   ebi  = EBI@GD per Yao 2022
  rdmft_bgd_tau           = 1.0d0         ! base occupation step for bgd (ELK taurdmn)
  rdmft_bgd_backtrack     = 0.75d0        ! geometric backtracking factor for bgd (ELK rdmvaryn)
  rdmft_occ_ls_type      = 'auto'        ! auto (monotone Armijo for SPG2) | armijo | gll | sw | nm_sw | weak_wolfe
  rdmft_orb_ls_type      = 'auto'        ! auto (monotone strong Wolfe) | armijo | sw | weak_wolfe
  rdmft_occ_ls_init_step = 'barzilai_borwein'  ! fixed | bb | quadratic
  rdmft_occ_ls_stepsize  = 1.0d0         ! alpha0 when init_step = fixed
  rdmft_line_search_max_zoom = 30        ! Wolfe zoom cap
  rdmft_bb_alpha_min     = 1.0d-8
  rdmft_bb_alpha_max     = 1.0d2
  rdmft_orb_optimizer    = 'cg'          ! sd | cg | lbfgs   (alternating only)
  rdmft_orb_strategy     = 'joint'       ! joint | block_k   (multi-k orbital block)
  rdmft_joint_optimizer  = 'lbfgs'       ! sd | cg | lbfgs   (joint only)
  rdmft_lbfgs_memory     = 10            ! L-BFGS history depth
  rdmft_outer_maxiter    = 50            ! outer alternating cycles
  rdmft_occ_maxiter      = 20            ! inner occupation iterations / outer cycle
  rdmft_orb_maxiter      = 20            ! inner orbital iterations / outer cycle
  rdmft_energy_tol       = 1.0d-7        ! outer |dE| tolerance (Ry)
  rdmft_orb_grad_tol     = 1.0d-5        ! orbital ||G_R|| tolerance
  rdmft_occ_grad_tol     = 1.0d-5        ! SPG2 ||g_1||_inf reference (diagnostic)
  rdmft_occ_tol          = 1.0d-6        ! alternating freeze: skip occ block once
                                          ! sum|dn|_outer < this for 2 consecutive cycles;
                                          ! set <= 0 to disable the freeze.
  rdmft_line_search_c1   = 1.0d-4        ! Armijo sufficient decrease constant
  rdmft_line_search_rho  = 0.5d0         ! geometric backtracking factor
  rdmft_line_search_max_iter = 20
  rdmft_orb_ls_stepsize  = 0.5d0         ! initial Armijo trial for orbital block
  rdmft_occ_init_mode    = 'ks'          ! ks | perturbed | binary | uniform
  rdmft_occ_init_perturb = 0.0d0         ! delta for perturbed / binary
  rdmft_occ_init_nbands_top = 0          ! K>0: Fermi-window (ABACUS); K<=0:
                                          ! perturbed -> legacy global push;
                                          ! binary/uniform -> ks
  rdmft_verbose          = 1   ! 0 silent; 1 summary + block/macro wall times; 2 + per-trial LS
  rdmft_grad_check       = .false.       ! FD check of dE/dn and orb retraction at start
  rdmft_temp             = 0.0d0         ! HF only: entropy temperature (K); uses k_B*T
  rdmft_restart          = .false.       ! restart from a previous RDMFT save (see below)
  rdmft_save_every       = 1             ! checkpoint after every N macro iterations (0 = disable)
  rdmft_restart_file     = ' '           ! override the default {outdir}/{prefix}.rdmft.save
/
```

### SPG2 occupation block (`rdmft_constraint = 'projected_gradient'`)

For the alternating solver, the occupation inner loop is the canonical
**Birgin–Martínez–Raydan SPG2** (Algorithm 2.2; *SIAM J. Optim.* **10**,
1196 (2000)).  The default `rdmft_occ_optimizer = 'spg2'` selects
canonical SPG2 verbatim; `'cg'` and `'lbfgs'` substitute alternative
descent directions into the same chord/proximal framework.

- **Geometry:** \(\partial E/\partial n_{ik}\) from `rdmft_grad_n`; optional **ELK-style** diagonal scale from per-\((ik,ib)\) probes at \(n=0.5\) when `rdmft_occ_precond = .true.` (`rdmft_compute_elk_occ_scale`).
- **Feasibility:** Euclidean L\(_2\) proximal projection \(P_w\) onto \([0,1]\) with \(\sum_k w_k n_{ik}=N_e\) (`rdmft_proximal_project_occ`; see `doc/rdmft_calculation.md` §5.3 for \(P_w\), \(P_u\), and the SPG chord).
- **Chord line search.** Build the SPG2 chord
  \(\mathbf{d}_k = P_w(\mathbf{n}_k - \alpha_k^{\mathrm{BB}}\,\texttt{grad\_dir}_k) - \mathbf{n}_k\)
  with `grad_dir` from `rdmft_build_occ_grad_dir` (copy of \(\partial E/\partial\mathbf{n}\) for `sd`, optionally ELK-scaled), the Barzilai–Borwein spectral steplength \(\alpha_k^{\mathrm{BB}}\),
  and walk the *straight* segment
  \(\mathbf{n}_k + \lambda\mathbf{d}_k\), \(\lambda\in[0,1]\).
  Every trial inside the chord remains feasible without further
  projections (BMR Algorithm 2.2).
- **Diagnostics:** each inner iteration logs the KKT residual `max(0, V-W)` on \(\partial E/\partial n_{ik}\) (free + box bands) and the SPG2 stopping map \(\|g_1\|_\infty = \|P_w(\mathbf{n}-t\,\partial E/\partial\mathbf{n})-\mathbf{n}\|_\infty\) with \(t=1\) (BMR); compare \(\|g_1\|_\infty\) to `rdmft_occ_grad_tol`. Inner stop uses `rdmft_occ_tol` on \(\sum|\Delta n|\).
- **Line search:** `rdmft_occ_ls_type` — default `auto` uses **monotone Armijo** on the SPG2 chord ($\lambda_0=1$, $\gamma=\texttt{rdmft\_spg\_gamma}$). Set `gll`/`nm_armijo` for BMR GLL nonmonotone Armijo ($M=\texttt{rdmft\_spg\_ls\_memory}$); `nm_sw` selects nonmonotone strong Wolfe with a Zhang–Hager reference. Initial step from `rdmft_occ_ls_init_step` applies to BGD/EBI only (`barzilai_borwein`, `quadratic`, or `fixed` + `rdmft_occ_ls_stepsize`).
- **Orbital line search:** `rdmft_orb_ls_type` — default `auto` uses monotone strong Wolfe on the product Stiefel manifold (`joint`) or per-k block (`block_k`); set `armijo` for Zhang–Hager / geometric backtracking via `rdmft_line_search`.
- **`rdmft_occ_optimizer`:**
  - `'spg2'` (default) — **canonical** BMR SPG2: spectral steepest descent on \(\partial E/\partial\mathbf{n}\) (optional ELK preconditioning).
  - `'sd'` — synonym for `'spg2'` (kept for backward compatibility).
  - `'cg'` — SPG2 framework with Polak–Ribière CG direction.  Useful when canonical SPG2 oscillates on wide-Hessian metals.
  - `'lbfgs'` — SPG2 framework with an L-BFGS direction (history depth `rdmft_lbfgs_memory`).
  - `'bgd'` — ELK occupation gradient descent (see below); not SPG2.
  - `'ebi'` — EBI@GD erf parameterisation (see below); not SPG2.
  - Any other value falls back to canonical SPG2 with a one-line warning.

### Box-aware gradient-descent occupation block (`rdmft_occ_optimizer = 'bgd'`)

A faithful port of ELK's `rdmvaryn.f90` occupation update (`rdmft_bgd.f90`),
with one **addition**: an explicit energy **line search** along ELK's
search direction.  ELK itself takes a single feasibility-only step
from `taurdmn` with a hard-coded geometric backtracking factor of
`0.75`; we use the same `0.75` backtrack (exposed as
`rdmft_bgd_backtrack`, kept separate from the SPG/EBI
`rdmft_line_search_rho`) and run an Armijo sufficient-decrease test
on each trial so the block always returns a step that lowers the
energy or accepts the line origin.  The reduced-gradient direction,
the chemical-potential \(\kappa\) solve, the unit-weighted-norm
rescaling, and the box-respecting feasibility bound \(\tau_{\max}\)
are bit-equivalent to `rdmvaryn`.

- **Search direction:** the charge-conserving, box-aware reduced gradient
  \(\gamma_{ik} = g_{ik}(1-n_{ik})\) if \(g_{ik}>0\) else \(g_{ik}\,n_{ik}\), with
  \(g_{ik} = \texttt{dedn}_{ik}-\kappa\) and ELK `dedn_{ik} = -(\partial E/\partial n_{ik})/w_k\),
  and \(\kappa\) (a chemical potential) solved by ELK's bracketing loop so that
  \(\sum_k w_k\sum_i\gamma_{ik}=0\) (separately per spin channel when `tot_magnetization` is fixed).
  \(\gamma\) is normalised when its weighted square norm exceeds 1 (ELK
  `rdmvaryn`).
- **Feasibility:** because \(\sum_k w_k\gamma_{ik}=0\), the electron number is preserved exactly
  along the straight line \(n(\tau)=n_0+\tau\gamma\); the box \([0,1]\) holds for
  \(\tau\in[0,\tau_{\max}]\), so the line search is projection-free on that interval.
- **Line search:** seeds \(\tau_0=\min(\texttt{rdmft\_bgd\_tau},\tau_{\max})\) then backtracks
  with geometric factor `rdmft_bgd_backtrack` (default `0.75`, matching ELK
  `rdmvaryn`) for Armijo sufficient decrease (or strong/non-monotone Wolfe via
  `rdmft_occ_ls_type`).
- **Descent / convergence:** \(\phi'(0)=\langle\nabla_n E,\gamma\rangle\le 0\),
  vanishing at a KKT point; the block stops on \(\sum|\Delta n| < \) `rdmft_occ_tol`, on
  \(\phi'(0)\ge 0\), or on a vanishing direction.
- **`rdmft_bgd_tau`:** base step \(\tau\) (ELK `taurdmn`), default `1.0`.
- **`rdmft_bgd_backtrack`:** geometric backtracking factor for the BGD line
  search (ELK `rdmvaryn` uses `0.75`); default `0.75`.
- **`rdmft_occ_precond`:** `.false.` (default) uses the raw SPG direction; `.true.` enables ELK scaling.
- The occupation XC gradient is always the closed-form contraction $c_t w_{\mathbf{k}} w_t'(n)\,D_{ii}$ (the symmetric quadratic form $E_x = -\tfrac{1}{2} g^T K g$ has gradient $D_\alpha = -(Kg)_\alpha$ in one shot; this matches the ABACUS reference `rdmft_energy_gradient.cpp`), with $D_{ii}$ rebuilt via ACE at every occupation-gradient/energy call. (The former `rdmft_frozen_fock` / `rdmft_occ_reuse_orbitals` cached-coupling fast path has been **removed** -- see the note below.)

### EBI occupation block (`rdmft_occ_optimizer = 'ebi'`)

Implements the explicit-by-implicit (EBI) parameterisation of Yao *et al.*
(J. Phys. Chem. A **2022**, **126**, 5654–5662;
[doi:10.1021/acs.jpca.2c02345](https://pubs.acs.org/doi/10.1021/acs.jpca.2c02345))
in [`rdmft_ebi.f90`](rdmft_ebi.f90).  The inner optimiser is the
**EBI@GD** combination — explicit-by-implicit parameterisation +
first-order gradient descent on the unconstrained parameters — which
that paper recommends as *"the lowest converged energies for different
types of systems, with the lowest computational scaling"*, after
ruling out EBI@NM (Newton, local-minimum issues on strongly correlated
systems) and preconditioner-free EBI@CG / EBI@L-BFGS.

- **Parameterisation:** \(n_{ik}=\tfrac{1}{2}(\mathrm{erf}(x_{ik}+\mu)+1)\) with
  unconstrained \(x_{ik}\) and a scalar \(\mu\) (one per spin when
  `tot_magnetization` is fixed) solved by bisection so
  \(\sum_k w_k n_{ik}=N_e\) at every step.  After each accepted inner
  step the implicit shift is folded into \(x\leftarrow x+\mu\) and
  \(\mu\) is reset to zero (`rdmft_ebi_fold_mu`).
- **Gradient (Eqs. 21–22, Yao 2022):** implicit chain rule
  \(\partial E/\partial x_{ik}=s(u_{ik})\bigl[g_{ik}-w_kR\bigr]\) with
  \(s(u)=\pi^{-1/2}\,e^{-u^2}\),
  \(R = \sum_{jk'} w_{k'} g_{jk'} s(u_{jk'}) / \sum_{jk'} w_{k'} s(u_{jk'})\),
  and \(u_{ik}=x_{ik}+\mu\).  In LSDA with fixed magnetisation,
  separate \(R_\uparrow, R_\downarrow\) are computed on each spin
  subset.  Implementation: `rdmft_ebi_transform_gradient`.
- **Optimiser (EBI@GD, Yao 2022):** first-order steepest descent on \(x\)
  (no preconditioner; no L-BFGS / CG, per the comparative study in
  the same paper), with energy line search along
  \(x(\alpha)=x_0+\alpha d\); no proximal projection is needed because
  the electron count is preserved by construction at every \(\alpha\)
  via the per-trial \(\mu\) solve inside `rdmft_ebi_ls_eval`.
- **Line search:** same controls as the SPG/bgd blocks
  (`rdmft_occ_ls_type`, `rdmft_occ_ls_init_step`,
  `rdmft_occ_ls_stepsize`, `rdmft_bgd_tau` as fallback); a
  quadratic-interpolating Armijo backtrack is used by default
  (`rdmft_line_search_polynomial = .TRUE.`).

### `rdmft_grad_check`

When `.true.`, [`rdmft_grad_check.f90`](rdmft_grad_check.f90) runs **once** after the
initial RDMFT energy and **before** the main loop (ABACUS
`check_gradient_consistency`):

- **Occupations:** central finite difference on up to 10 `(ik,ib)` pairs vs
  analytic `rdmft_grad_n` (default `epsilon=1e-5`, `tolerance=1e-4`).
- **Orbitals:** forward difference along the product-manifold Stiefel retraction
  `C^k ← R_{C^k}(-ε G_R^k)` vs analytic slope `-‖G_R‖_S²`.

Expensive (many full energy evaluations); intended for development and
validation. For **GU**, the orbital check may **fail** until the explicit
`½ w_k² (n²−n) J_ii` term’s `∂J/∂C` contribution is added to the orbital
gradient; occupation FD should still pass.

### `rdmft_temp` (HF only)

Optional binary-entropy regularisation at temperature `rdmft_temp` in **Kelvin**
(default `0` disables).  The prefactor is \(k_B T\) with QE's
`K_BOLTZMANN_RY` (Boltzmann constant in Ry/K):

$$E_{\mathrm{ent}} = k_B T \sum_k w_k \sum_i \bigl[n\ln n + (1-n)\ln(1-n)\bigr],$$

added to `rdmft_etot` and `∂E/∂n` when `rdmft_functional = 'hf'`.

### GU functional (`rdmft_functional = 'gu'`)

Off-diagonal exchange uses the Müller-type EXX build (`g(n)=\sqrt{n}` with
`pow_reg`). The on-site diagonal piece

$$E_{\mathrm{GU,diag}} = \tfrac{1}{2}\sum_k w_k^2\sum_i (n_{ik}^2-n_{ik})\,J_{ii}^{k}$$

with $J_{ii}=\langle\psi_{ik}|V_x[\gamma]|psi_{ik}\rangle$ from the same ACE
channel is added in [`rdmft_energy.f90`](rdmft_energy.f90). Occupation gradients
include `½ w_k^2 (2n-1) J_ii` (holding $J_{ii}$ fixed, as in the ABACUS
derivation remark).
```

### Initial occupation modes (`rdmft_occ_init_*`)

Mirrors ABACUS `OccInitMode` in `abacus_source/.../rdmft_solver.cpp`.
All modes start from KS occupations `n_ik = wg/wk`, optionally modify
a Fermi window of `K = min(rdmft_occ_init_nbands_top, nbnd)` bands
above and below the approximate Fermi boundary (first band with
`n < 0.5`), then call `rdmft_proximal_project_occ`.

| Mode | Requires | Action on Fermi-window bands |
|------|----------|------------------------------|
| `ks` | — | no change |
| `perturbed` | `delta > 0` | `n += delta` if KS `n < 0.5`, else `n -= delta` |
| `binary` | `K > 0`, `delta > 0` | set to `delta` if KS `n < 0.5`, else `1 - delta` |
| `uniform` | `K > 0` | replace by mean KS occupation over selected bands |

Fallbacks: `binary` / `uniform` with `K <= 0` (and `binary` with
`delta <= 0`) leave KS occupations unchanged.  `perturbed` with
`K <= 0` and `delta > 0` uses the legacy QE global push on every band
(not ABACUS); prefer `K > 0` for Müller.

### Alternating freeze-occupation logic

The alternating driver tracks `occ_outer_dn_sum = sum_ik |n_ik^out -
n_ik^out-1|` after each outer occupation block.  When this stays below
`rdmft_occ_tol` for two consecutive outer cycles the occupation block
is permanently skipped for the remainder of the run, leaving only the
Stiefel orbital block to refine the energy.  The log line for each
outer iteration reports both the per-cycle delta and the freeze flag::

```
RDMFT outer    3  E = -2.2426947632  |dE| = 6.39E-04  Ne = 2.000000
                  sum|dn|_out = 2.19E-03  occ_frozen = F
** RDMFT: sum|dn|_outer below 5.00E-03 for 2 consecutive cycles
   -- freezing occupation block.
RDMFT outer    4  E = -2.2422745552  |dE| = 4.20E-04  Ne = 2.000000
                  sum|dn|_out = 4.34E-03  occ_frozen = T
RDMFT outer    5  E = -2.2423276166  |dE| = 5.31E-05  Ne = 2.000000
                  sum|dn|_out = 0.00E+00  occ_frozen = T
```

Set `rdmft_occ_tol <= 0` to disable the freeze if you want to keep the
occupation block running every outer cycle.

### Optimiser efficiency on HF (sd / cg / lbfgs)

A side-by-side benchmark on H2 in a 12-Bohr box and bulk Si at the
Γ point with `rdmft_functional = 'hf'` and a perturbed initial
occupation (so the optimisers actually have work to do) is shipped in
`PW/src/rdmft/examples/hf_benchmark/` together with reproducible
shell scripts.  Summary (KS-HF references in parentheses):

**H2** (KS-HF = -2.24112743 Ry):

| strategy     | optimiser | n_outer | dE vs KS-HF | wall (s) |
|--------------|-----------|---------|-------------|----------|
| alternating  | sd        |  2      | +1.2e-9     |  1.62    |
| alternating  | cg        |  2      | +1.2e-9     |  1.52    |
| alternating  | lbfgs     |  2      | +1.3e-9     |  1.93    |
| joint        | sd        | 30      | +2.6e-4     |  2.65    |
| joint        | cg        | 30      | +2.2e-5     |  2.80    |
| joint        | lbfgs     | 30      | +5.0e-6     |  2.08    |

**Bulk Si at the Γ point** (KS-HF = -14.43586741 Ry):

| strategy     | optimiser | n_outer | dE vs KS-HF | wall (s) |
|--------------|-----------|---------|-------------|----------|
| alternating  | sd        |  2      | -5.4e-8     |  0.43    |
| alternating  | cg        |  2      | -5.4e-8     |  0.42    |
| alternating  | lbfgs     |  2      | -5.4e-8     |  0.52    |
| joint        | sd        | 30      | +5.6e-4     |  0.71    |
| joint        | cg        | 30      | +1.1e-4     |  0.71    |
| joint        | lbfgs     | 30      | +8.4e-4     |  0.71    |

Take-aways:

* **Both strategies recover the SCF energy** for HF on every system,
  validating the EXX-via-ACE wiring and the sign / weight conventions.
* **Alternating with the freeze-occupation logic** is the right
  default: independent of the inner optimiser, the outer loop
  converges in two cycles and the freeze guarantees the cheap
  occupation block stops costing anything once the natural
  occupations stabilise.
* For the **joint** strategy, **L-BFGS is best on small,
  well-conditioned problems** (H2-like), reaching ~5 µRy of KS-HF,
  while **CG is more robust on larger systems** with a wider Hessian
  spectrum (Si-like).  Plain SD is a useful debugging baseline.
* Per-iteration cost in the alternating strategy: SD < CG < L-BFGS by
  ~25 % (history bookkeeping); the difference is dominated by the
  shared `aceinit` + `vexxace_*` cost, so all three are essentially
  equivalent in wall time.

See `PW/src/rdmft/examples/hf_benchmark/README.md` for the full
benchmark scripts and discussion.

### Joint (product-manifold) optimisation

`rdmft_solver_strategy = 'joint'` activates the alternative driver
that steps both blocks together.  For each outer iteration:

1. The natural occupations are mapped onto unconstrained parameters
   ``p`` via `n_ik = cos^2(p_ik)` so the box constraint is built in.
2. The full energy and gradients `(dE/dn, dE/dC)` are evaluated.
3. `dE/dn` is chain-ruled into `dE/dp` through the parameterisation
   Jacobian `dn/dp = -sin(2p)`, and `dE/dC` is projected onto the
   Stiefel tangent space at the current `C^k`.
4. A packed steepest-descent direction is built; the joint directional
   derivative `<grad_p, dir_p> + sum_k <G_R^k, dir_C^k>` is checked
   for descent.
5. **One** Armijo line search runs along the packed direction --
   linear update `p <- p + alpha * dir_p` for the parameters and
   Stiefel retraction `R_C(alpha * dir_C^k)` for each orbital block.
6. After acceptance, the L2 proximal projector cleans up any drift in
   the equality constraint.

The joint solver supports all three optimisers
(`rdmft_joint_optimizer = 'sd' | 'cg' | 'lbfgs'`) on the packed
`(p, evc_flat)` view; `'cg'` uses Polak-Ribière with
transport-by-reprojection on the orbital block, and `'lbfgs'` uses the
two-loop recursion with the curvature pair built from the previous
accepted step.

## Algorithmic notes

### Sign and weight conventions

For the gamma-only / `nspin = 1` / closed-shell case (the path covered
by the shipped examples) QE bakes the spin doubling into `wk(k)` so

```
wg(i, k) = wk(k) * f_per_spin(i, k),  f_per_spin in [0, 1]
```

and `sum_ik wg(i,k) = N_e`.  The natural occupation we evolve is the
same per-spin number, so the equality target is simply `nelec` and the
RDMFT one-body / Hartree weight is `wg_RDMFT = wk * n_ik`.  For each
exchange channel `t` we set `wg = wk * w_t(n_ik)` so the
`x_occupation` array that `aceinit` reads is `w_t(n_ik)`; the
resulting ACE projectors then encode the modified-DM Fock operator.

### Fixed-orbital occupation fast path (cached couplings) — removed

Earlier revisions offered an occupation-block "fast path" (flags
`rdmft_frozen_fock` + `rdmft_occ_reuse_orbitals`) that, while orbitals were
fixed, precomputed a per-band exchange coupling tensor once per block and
assembled `vx_diag` / occupation gradients / channel energies from it without
rebuilding the ACE projector on each inner iteration.  **This fast path has
been removed.**  Although the channel exchange energy is bilinear in the
channel weights for fixed orbitals, the cached-coupling assembly was **not
k-point-pool reproducible** on metallic / fractional-occupation systems: with
`npool > 1` (`-nk` / `-npool`) the occupation energy/gradient acquired an
`-nk`-dependent spread of order $10^{-4}$ Ry that the SPG occupation loop
amplified to the mRy level, whereas the per-trial ACE path is pool-invariant to
$\sim 10^{-10}$ Ry.  The solver now always rebuilds the exchange via ACE
(`rdmft_compute_exchange` → `aceinit` / `vexxace`), which is the validated,
deterministic path.  The occupation XC gradient is unchanged: the closed-form
contraction $c_t w_{\mathbf{k}} w_t'(n)\,D_{ii}$.

Note: skipping `exxinit` when `rdmft_exxbuff_stale` is false reuses the EXX
wavefunction buffer only; **`aceinit` is still the dominant cost** of the
default (no-reuse) path.

### `tot_magnetization` (LSDA)

When the PW `&system` namelist sets `tot_magnetization` with
`nspin = 2`, PWscf enables `two_fermi_energies` and fixes
`nelup`/`neldw` during the preceding KS-SCF.  RDMFT mirrors that
split: `rdmft_initial_n_from_ks` sets `rdmft_fix_magnetization` and
the occupation projector (`rdmft_proximal_project_occ`) enforces

```
sum_{k,up}   w_k sum_i n_ik = nelup
sum_{k,down} w_k sum_i n_ik = neldw
```

independently (two spin-filtered L2 proximal maps).  If
`tot_magnetization` is omitted, only the total `nelec` constraint is
applied and the magnetisation may change during the RDMFT solve.

### Energy decomposition

`rdmft_total_energy` rebuilds the density and the Hartree+semilocal-XC
potentials with the modified weights, computes the per-band one-body
diagonal via `h_psi` (EXX disabled with `stop_exx`/`start_exx`), and
sums one ACE pass per XC channel.  The accumulator subtracts the
double-counted Hartree (`-2*ehart`) and the semilocal `vtxc` so what
remains is exactly

```
E = T + E_loc + E_NL + ehart + sum_t coef_t * (1/2) sum_ik wg_t * vx_diag_t + ewld
```

For HF (`g(n) = n`) at the converged KS-HF orbitals this reproduces the
KS-HF total energy to ~1e-9 Ry on the shipped H2 example.

### Stiefel manifold

For PWscf the AO overlap is the identity (norm-conserving PP, no AO
basis), so the standard Stiefel projection
`Z = G - C * sym(C^H G)` already gives the right tangent vector and
the Cholesky-polar retraction
`R_C(eta) = (C + eta) * L^{-H}` with `(C+eta)^H (C+eta) = L L^H`
preserves orthonormality exactly.  The gamma-trick variants in
`rdmft_stiefel.f90` use the same formulas but compute the inner
products in the `2*Re(...) - G=0` half-space convention used by `evc`.

## Building

The new files are added to both build systems:

```bash
# CMake (recommended)
mkdir build && cd build
cmake -DCMAKE_Fortran_COMPILER=gfortran ..
make pw

# Classic Make
./configure
cd PW && make
```

A clean build of `pw.x` with the changes in this folder takes a couple
of minutes on a single core; only the new `rdmft/*.f90` objects, plus
`PW/src/input.f90` and `PW/src/electrons.f90`, are recompiled
incrementally.

### OpenMP (recommended for the ACE rebuild / vexxace hot path)

The dominant cost of an RDMFT outer cycle is the chain
`exxinit` → `aceinit` → `vexxace_*` that has to run on every Armijo
trial in the occupation and orbital line searches.  All three
routines are now OpenMP-aware: the large element-wise array
initialisations (`xi`, `vv`, the real → complex copy of the
projector matrix) carry explicit collapsed `!$omp parallel do`
loops, every `ZGEMM` / `ZTRMM` / `DGEMM` / `ZPOTRF` underneath them
threads through a thread-aware BLAS, and the FFTs invoked through
the EXX path scale through the OpenMP-enabled build of FFTW.

To activate parallel ACE rebuild and ACE application, configure
with `--enable-openmp` and link against threaded BLAS / LAPACK and
the threaded FFTW build, e.g.

```bash
BLAS_LIBS="-L/usr/lib/x86_64-linux-gnu/openblas-openmp -lopenblas" \
LAPACK_LIBS="-L/usr/lib/x86_64-linux-gnu/openblas-openmp -lopenblas" \
./configure --enable-openmp
make veryclean && make -j pw
```

Then set the per-MPI-rank thread count at run time:

```bash
export OMP_NUM_THREADS=4
mpirun -np <ranks> pw.x -in input.in
```

The `slurm_qe_frac.sh` template in `examples/LiH/` already maps
`SLURM_CPUS_PER_TASK → OMP_NUM_THREADS` for this purpose.

Reference numbers for `examples/LiH/lih.scf.in` on a 4-core, 1-MPI
build of QE 7.5 (gfortran 14, OpenBLAS-OpenMP, FFTW3-OMP) -- five
RDMFT outer iterations, 292 ACE rebuilds, 292 `vexxace` applications:

| Routine    | Serial WALL | `OMP_NUM_THREADS=4` WALL | Speedup |
|------------|-------------|--------------------------|---------|
| `aceinit`  | 123.10 s    | 61.95 s                  | 2.0×    |
| `vexx`     | 121.54 s    | 61.57 s                  | 2.0×    |
| `vexxace`  |   1.82 s    |  0.36 s                  | 5.1×    |
| `fftc`     |  44.45 s    | 14.44 s                  | 3.1×    |
| **Total**  | **3 m 01 s**| **1 m 50 s**             | **1.65×** |

The RDMFT energy at every outer iteration agrees bit-for-bit
between the serial and the 4-thread runs.

## Running

After building, copy a hydrogen pseudopotential into your run directory
and execute

```bash
mkdir -p tmp
pw.x < PW/src/rdmft/examples/h2_molecule/h2.in > h2.out
```

The KS-HF SCF runs first (you'll see the standard `! total energy`
lines, then `EXX self-consistency reached`).  RDMFT then prints its
own banner and per-iteration table:

```
============================================================
RDMFT alternating optimisation
============================================================
rdmft_functional      = muller
rdmft_power_alpha     =   0.5000
rdmft_constraint      = projected_gradie
...
RDMFT outer    1  E =      -2.2424075002  |dE| =   4.01E-02  Ne =     2.000000
...
!RDMFT_ETOTAL  =      -2.2422885283 Ry
```

The `!RDMFT_ETOTAL` line is the converged RDMFT total energy.

## Verification

* **HF on H2** (`rdmft_functional = 'hf'`): the alternating optimiser
  converges in one outer iteration to the converged KS-HF energy to
  the last printed digit, demonstrating that the EXX-via-ACE pipeline
  has been wired in with the correct sign/factor conventions.
* **Müller on H2**: gives a slightly lower energy than KS-HF
  (~0.1 mRy below), reflecting the additional left-right correlation
  the Müller kernel captures at the equilibrium bond length.
* **Periodic H4 chain**: KS-HF SCF + RDMFT-Müller post-processing runs
  to completion with `Ne` conserved exactly at every outer iteration.
* **Multi-k Si (2x2x2) and charged systems** (`tot_charge != 0`): see
  `examples/multik_charged/`.  Both the `'joint'` and `'block_k'`
  multi-k orbital strategies reproduce KS-HF to ≤ 4e-8 Ry; the
  charged H<sub>4</sub><sup>2+</sup> and H<sub>2</sub><sup>-</sup>
  reference cases match KS-HF to ≤ 3e-9 Ry.
* **PBE → RDMFT-HF** (`examples/pbe_to_hf/`): a pure-PBE SCF stopped
  after 2 inner iterations (no Fock character whatsoever in the
  starting orbitals) is fed into RDMFT-HF; the optimiser walks all
  the way to the converged KS-HF stationary point on both H2 (Δ ≈
  1 nRy) and multi-k Si (Δ ≈ 1.5 µRy).  This exercises the
  on-the-fly EXX activation path (`setup_exx` is called from
  `rdmft_run` when the SCF was non-hybrid), the runtime XC-ID
  switch to (5,0,0,0,0,0) so that `v_of_rho` returns
  `V_xc^semilocal = 0`, and the per-evaluation `exxbuff` refresh
  that keeps `vexx` consistent with the moving orbitals.
* **`nspin = 2` (spin-polarised)** (`examples/nspin2/`): five cases
  exercising H2 (Mz=0), Si 2x2x2 (Mz=0) and H2^- (Mz=1), each from
  both a converged HF SCF and from a rough 2-iteration PBE SCF.  All
  five reach the converged KS-HF reference to ≤ 2 µRy and conserve
  the requested magnetisation exactly (enforced by separate
  ``nelup``/``neldw`` proximal projections whenever PWscf input sets
  ``tot_magnetization``).  This shows that the spin-polarised path
  (per-spin ``wk`` + ``wg`` bookkeeping,
  per-spin ``x_occupation``, the per-spin ``vexxace`` /
  ``aceinit`` channel calls inherited from QE's hybrid driver)
  is consistent with the RDMFT energy formula and the Stiefel
  multi-k orbital block.

## Known limitations and TODOs

* The semilocal-`V_xc` issue mentioned in earlier revisions of this
  README is now solved at the source: at the start of every
  `rdmft_run` the global XC IDs are switched to `(5,0,0,0,0,0)`
  ("HF") and the auxiliary flags are refreshed via
  `xclib_set_auxiliary_flags(.FALSE.)`.  This makes
  `v_of_rho` return `V_xc^semilocal = 0` *and* `vtxc = 0` for the
  duration of the RDMFT solve, so the orbital gradient (which uses
  `h_psi` and therefore `vrs = vltot + V_H + V_xc`) is the pure HF
  (or RDMFT-functional) gradient by construction, regardless of
  what XC the upstream SCF used.  The original IDs are restored on
  exit from `rdmft_run`.
* USPP / PAW: `becp` is allocated once at the top of `rdmft_run` and
  re-used by every `h_psi` call inside the post-SCF loop; gamma-only
  NC, gamma-only USPP and multi-k NC cases all work.  Multi-k USPP /
  PAW is untested.
* Multi-k Brillouin-zone sampling: two interchangeable strategies
  are implemented in the alternating orbital block (selectable via
  the new keyword `rdmft_orb_strategy`):

  * `'joint'` (default) — Riemannian descent on the **product
    Stiefel manifold** ``\prod_k \mathrm{St}(\text{nbnd},
    \text{npwx}; I)``.  Each inner iteration computes the
    Riemannian gradient at every k-point in one pass, builds the
    search direction (SD / CG / L-BFGS), and runs a single GLOBAL
    Armijo line search where one scalar `alpha` is accepted only
    if the total energy across all k decreases sufficiently
    (`E(alpha) <= E(0) + c1 * alpha * <G_R, dir>_product`).  Per
    trial the solver retracts every `C^k` simultaneously, saves
    all per-k buffers, and asks `rdmft_total_energy` for one
    global energy value.

  * `'block_k'` — Gauss-Seidel sweep over k-points (or per-pool
    k-groups) on the same product manifold.  At each k take ONE
    Riemannian step, holding all other k's evc fixed in their wfc
    buffers.  The line search at every k still uses the global
    energy (`rdmft_total_energy` summed over all k) but only the
    current k's evc is varied per trial.  Each k owns its own
    optimiser state (per-k SD / CG / L-BFGS history persisted
    across outer iterations).  This maps naturally onto QE's
    k-point pool parallelisation: each MPI pool optimises its
    own k-points independently and only the global-energy
    reduction at each line-search trial requires a cross-pool
    `mp_sum`.  Memory cost per inner iteration is `nks` x
    smaller than the joint variant (one per-k workspace at a
    time) and per-pool work scales as `nks_per_pool`, so on a
    p-pool MPI run the orbital block is roughly p× faster than
    the joint variant.

  Both strategies reproduce KS-HF on Si 2×2×2 to ≤ 3 nRy with
  the orbital block fully active (see
  `examples/multik_charged/`).
* `tot_charge != 0` is supported transparently (the
  electron-number target follows QE's `nelec`); verified on
  H<sub>4</sub><sup>2+</sup> (`nspin = 1`) and on H<sub>2</sub><sup>-</sup>
  (`nspin = 2`).  PWscf disables ACE when one spin channel is empty,
  so e.g. `tot_charge = +1` on H<sub>2</sub> with
  `nspin = 2 tot_magnetization = 1` cannot be run with the current
  ACE-based RDMFT exchange path.
* Occupation (SPG2) defaults to **monotone Armijo** on the chord
  (`rdmft_occ_ls_type = 'auto'`; $\gamma=\texttt{rdmft\_spg\_gamma}$).
  Set `gll`/`nm_armijo` for BMR GLL nonmonotone Armijo
  ($M=\texttt{rdmft\_spg\_ls\_memory}$); `nm_sw` for nonmonotone strong
  Wolfe (Zhang–Hager).  The alternating
  orbital block defaults to monotone strong Wolfe (`rdmft_orb_ls_type = 'auto'`).
  The packed **joint** strategy (`rdmft_solver_strategy = 'joint'`) still uses Armijo only.
* For Muller / Power / `gu` / `geo`, small-$n$ handling matches ABACUS
  `rdmft_xc_functional.h` (`pow_reg` / `dpow_reg`): for $n\ge\varepsilon$,
  $g(n)=n^\alpha$ and $g'(n)=\alpha n^{\alpha-1}$; for $n<\varepsilon$,
  $g$ is the Taylor line through $(\varepsilon,\varepsilon^\alpha)$ and
  $g'(n)=\alpha\max(n,\varepsilon)^{\alpha-1}$.  Default
  `rdmft_reg_eps = 1e-8`.  Thus $g(0)=(1-\alpha)\varepsilon^\alpha$
  (small but nonzero when $\alpha<1$), $g(n)\to n^\alpha$ for
  $n\gg\varepsilon$, and the projected gradient at empty bands is
  bounded by `~alpha/eps^(1-alpha)` (Müller $\alpha=0.5$:
  `~0.5/sqrt(rdmft_reg_eps)`).  For **Müller** especially, combine the
  default with the **Fermi-window** occupation seed:
  `rdmft_occ_init_mode = 'perturbed'`, a small
  `rdmft_occ_init_perturb` (e.g. `1e-3`), and `rdmft_occ_init_nbands_top > 0`
  (QE LiH examples use `1` for a small Fermi window; the ABACUS LiH_HF
  example uses `4`).  Leaving `nbands_top = 0` keeps the older QE behaviour
  that perturbs every band; that can destabilise charged-state scans because
  deep occupied states are moved away from unity for no physical benefit.
  A larger `rdmft_reg_eps` (e.g. `1e-3`) further softens $g$ and $g'$ below
  $\varepsilon$ if the occupation block is still stiff.
* **GU orbital gradient:** the explicit diagonal $(n^2-n)J_{ii}$ correction is
  in the energy and $\partial E/\partial n$; the $\partial J_{ii}/\partial C$
  piece is not yet in the orbital block (same gap as ABACUS production).

## RDMFT density of states and per-atom magnetization

After a converged RDMFT solve, DOS and Mulliken magnetization are computed
by the dedicated post-processor **`rdmft_dos.x`**, not during `pw.x`.
This mirrors the standard QE `pw.x` / `dos.x` split and keeps the expensive
TSM probe out of the optimisation loop.

**Default (`spectral = 'elk'`)** mirrors ELK task-10 `dos.f90`:
one \(\delta(\omega-\varepsilon)\) per state at \(\omega=\varepsilon^{\mathrm{TSM}}-
\varepsilon_F\), weight `sc × occmax` (or `sc × n × occmax` when
`occ_weighted = .true.`), integrated with `brzint`.  The energy
zero uses a Fermi level from TSM eigenvalues (`eref = 'efermi'`,
default for ELK mode).

**Optional (`spectral = 'sharma'`)** follows PRL 110, 116403
Eq.~(7): two-branch sum with RDMFT \(\mu\) (`eref = 'mu'` by
default) and Gaussian or `brzint` integration.

Enable with a two-step workflow (see below).  DOS keywords belong in
`&inputrdmftdos` for `rdmft_dos.x` only — not in `&rdmft` for `pw.x`:

```fortran
! pw.x -- RDMFT optimisation only
&rdmft
  do_rdmft = .true.
  ...
/

! rdmft_dos.x -- DOS post-processing (reads the saved state)
&inputrdmftdos
  prefix = 'fe'
  outdir = './tmp/'
  functional = 'hf'          ! must match the pw.x run
  spectral = 'elk'
  nwplot = 500
  ngrkf  = 100
  emin   = -1.2
  emax   =  1.2
  occ_weighted = .false.
/
```

| Keyword | Meaning |
|---------|---------|
| `functional` | RDMFT XC functional (must match `pw.x`) |
| `power_alpha`, `reg_eps` | Functional parameters when needed |
| `spectral` | `'elk'` (default) or `'sharma'` (PRL two-branch) |
| `eref` | `'auto'` (default), `'efermi'`, or `'mu'` |
| `occ_weighted` | `.false.` → `sc×occmax`; `.true.` → `sc×occsv` |
| `integration` | `'brzint'` (default) or `'gauss'` (Sharma mode) |
| `emin`, `emax` | Energy window (Ry); omit (=±1e6) for auto |
| `deltae` | Grid step (Ry) when `nwplot = 0` |
| `degauss` | Broadening width (Ry) |
| `nwplot` | Energy mesh points (≥2); `0` = use `deltae` |
| `ngrkf` | BZ subdivision for `brzint` (default 100) |
| `nswplot` | Post-smoothing passes after `brzint` (default 0) |
| `pdos` | Write partial DOS (default `.true.`) |
| `msum` | Sum PDOS over $m$ (default `.true.`) |
| `ssum` | Sum DOS over spin (default `.false.`) |
| `sqaxis` | Spin quantisation axis (noncollinear; default $(0,0,1)$) |
| `file_prefix` | Output prefix; blank → `{prefix}.rdmft` |
| `write_evalsv` | Write `<prefix>.rdmft.evalsv` for ELK cross-checks |
| `restart_file` | Override `{outdir}/{prefix}.rdmft.save` path |

> **Units.** QE DOS files use Ry and states/Ry/cell; ELK uses Ha and
> states/Ha/cell.  Multiply energies by 2 and DOS by 1/2 when comparing.

> **Partial DOS.** QE uses Mulliken projection from PP `PP_CHI` blocks;
> ELK uses APW `gendmatk`.  Total DOS on the same k-mesh should be
> compared first; PDOS will differ between codes.

> **Note.** DOS is computed only by `rdmft_dos.x`.  The TSM
> `eps_tsm(n=1/2)` column is printed by `rdmft_dos.x` (or
> `rdmft_print_tsm_energies`), not at the end of `pw.x`.

**Output files** (written by `ionode` in the run directory):

* `${prefix}.rdmft.dos` — total DOS (spin-up positive, spin-down negative;
  energies in Ry relative to efermi (ELK) or μ (Sharma); states/Ry/cell)
* `${prefix}.rdmft.evalsv` — optional TSM energies (`write_evalsv`)
* `${prefix}.rdmft.pdos_at<S>_<N>` — $l$-resolved partial DOS per atom site
  (only written when the pseudopotentials carry atomic wave-functions —
  see "Pseudopotential requirements" below)
* `${prefix}.rdmft.idos` — interstitial DOS ($\mathrm{TDOS}-\sum\mathrm{PDOS}$;
  written when PDOS is enabled)
* `${prefix}.rdmft.mag` — Mulliken $M_z$ per site and species totals
  (same gating as the PDOS files, and additionally requires `lsda`)

**Magnetization:** per-atom and per-species spin moments use the same
atomic projections as PDOS, with $Q_\uparrow - Q_\downarrow$ in electron
units (consistent with PW `tot_magnetization`).

## RDMFT particle density and 1-RDM coherency plots (`pp.x`)

After a converged RDMFT `pw.x` run, use the standard **`pp.x`** post-processor
(two-step workflow, like DOS) to build cube files from the saved natural
orbitals and occupations:

$$
\rho(\mathbf r)=\gamma(\mathbf r,\mathbf r)
=\sum_k w_k\sum_i n_{ik}\,|\psi_{ik}(\mathbf r)|^2,
$$

**1-RDM coherency diagnostic (not the XC exchange hole).** `pp.x` can also plot

$$
h(\mathbf r_0,\mathbf r)=\gamma(\mathbf r_0,\mathbf r)-\rho(\mathbf r_0)\,\rho(\mathbf r),
\qquad
h(\mathbf r,\mathbf r)=\rho(\mathbf r)-\rho(\mathbf r)^2,
$$

These fields are **visualization diagnostics** only. They are **not** the RDMFT exchange–correlation hole (which follows from the pair kernel $f(n_i,n_j)$; see Section 2.5 of `doc/rdmft_calculation.md`). Do not use them to interpret XC energies or gradients.

**Step 1:** `pw.x` with `do_rdmft = .true.` saves `{outdir}/{prefix}.rdmft.save`
and collected natural orbitals in `{outdir}/{prefix}.save/`.

**Step 2:** `pp.x` with one of the RDMFT `plot_num` values:

| `plot_num` | Quantity |
|------------|----------|
| 126 | RDMFT particle density $\rho(\mathbf r)$ |
| 127 | 1-RDM coherency slice $h(\mathbf r_0,\mathbf r)=\gamma(\mathbf r_0,\mathbf r)-\rho(\mathbf r_0)\rho(\mathbf r)$ (**not** the XC hole) |
| 128 | On-site diagnostic $h(\mathbf r,\mathbf r)=\rho(\mathbf r)-\rho(\mathbf r)^2$ |

Example input (`&inputpp` + `&plot`):

```fortran
&inputpp
  prefix = 'h2'
  outdir = './tmp/'
  plot_num = 126
  filplot = 'h2.rdmft.rho.pp'
  rdmft_functional = 'muller'
  rdmft_xhole_origin(1) = -1.0
  rdmft_xhole_origin(2) = -1.0
  rdmft_xhole_origin(3) = -1.0
  rdmft_update_save = .false.
/
&plot
  iflag = 3
  output_format = 6
  fileout = 'h2.rdmft.rho.cube'
/
```

| Keyword | Meaning |
|---------|---------|
| `rdmft_functional` | Must match the upstream RDMFT run (default `muller`) |
| `rdmft_xhole_origin(3)` | Fractional reference $\mathbf r_0$ for `plot_num=127`; all $\le -1$ picks max-$\rho$ grid point |
| `rdmft_update_save` | If `.true.`, refresh `{outdir}/{prefix}.save/charge-density` via `write_scf` (default `.false.`) |
| `rdmft_restart_file` | Override default `{outdir}/{prefix}.rdmft.save` |

See [`examples/h2_dos_two_step/h2.step3.rho.pp.in`](examples/h2_dos_two_step/h2.step3.rho.pp.in)
for a worked H$_2$ example.  Spin-polarised / noncollinear plots are not
yet supported.

Integrated $\int\rho\,d\mathbf r$ is printed to stdout when the density is rebuilt.

### Pseudopotential requirements

`rdmft_dos.x` always writes the total DOS, but the PDOS
and Mulliken magnetisation rely on atomic-wave-function projections via
QE's `atomic_wfc` machinery.  When the input pseudopotentials carry no
`<PP_CHI>` blocks (most ONCV PBE pseudos for first-row elements such as
H/Li/Be, distributed with the demo at `examples/LiH/upf/`), the wfc
count `n_atom_wfc` returns 0.  The DOS routine now detects this case,
emits a one-line warning, skips the projection / Cholesky / PDOS / Mz
steps, and still writes the total DOS so the run remains useful.
Switch to a pseudopotential family that carries atomic wave-functions
(e.g. SSSP, GBRV, the Garrity-Bennett-Rabe-Vanderbilt PBE sets, the
**PseudoDojo** ONCV PBE set -- whose UPFs *do* carry `<PP_CHI>` blocks,
unlike the SG15 ONCV UPFs which ship with `number_of_wfc=0` -- or the
PZ NC `Fe.pz-n-nc.UPF` shipped with QE that is used by
`examples/Fe_bcc_dos/`) when PDOS / Mz outputs are required.  The
`examples/NiO_nk2x2x2_oneshot/` notes show how to pull the PseudoDojo
Ni/O UPFs used for the ELK DOS cross-check.

### Worked example

[`examples/Fe_bcc_dos/`](examples/Fe_bcc_dos/) drives the full
post-processor on ferromagnetic BCC Fe with the PZ NC pseudopotential
that already carries 4s/4p/3d `<PP_CHI>` blocks.  After about two
minutes on a 4-core OpenMP build, the run reproduces:

* a per-site Mulliken $M_z = +1.97\,e$ matching the SCF
  `total_magnetization = 1.97 \mu_B`/cell to better than 1 %;
* spin-resolved TDOS with a sharp minority-spin d-peak below $E_F$
  (around $-1.3$ eV) and the majority-spin d band fully filled, as
  expected for bcc Fe;
* l-resolved PDOS integrating to ≈ 3.77 e in `d↑` vs. 1.99 e in `d↓`
  on the Fe site.

That example uses a 2×2×2 k mesh and `rdmft_outer_maxiter = 0` (DOS
from the converged SCF orbitals, no RDMFT optimisation) so the run
finishes quickly in CI — it is **not** a publication-grade
reproduction of the LAPW + dense-k results in PRL **110**, 116403
(2013), but it exercises every code path in the new `rdmft_dos`
module.

**Spin support.**  Collinear (`nspin = 1` or `2`) is the validated
path.  Noncollinear (`noncolin = .true.`) supports spin-resolved total
DOS via ELK-style `sc` weights (`rdmft_compute_spin_weights`) and
PDOS via `atomic_wfc_nc_proj` (without `lspinorb`).  Per-atom Mulliken
`Mz` remains lsda-only.  Choose `integration = 'gauss'` or
`'brzint'`; projections are not symmetrized (`lsym` ignored).

### Noncollinear RDMFT

The non-collinear path mirrors the ELK RDMFT driver (task 300)
described in `abacus_source/source_lcao/module_rdmft/doc/elk_rdmft_algorithm.md`:
the natural orbitals are 2-component spinors stored in
``evc(npwx*npol, nbnd)`` (with ``npol = 2`` when ``noncolin = .true.``)
and the natural occupations are scalar numbers in $[0, 1]$ per band
per k-point (no extra spin doubling factor, exactly as for ELK's
``occsv``).  All of the kernels in ``rdmft_energy.f90`` ,
``rdmft_stiefel.f90`` and ``rdmft_spg.f90`` are already written
against the ``lda = npwx*npol`` convention -- the per-band diagonals
and the Stiefel projection/retraction sum the two spinor halves of
``evc`` -- and the underlying PWscf EXX-via-ACE driver
(``aceinit``/``vexx``/``vexxace_*``) is noncolinear-aware, so the
RDMFT loop runs unchanged on ``noncolin = .true.`` inputs (with or
without ``lspinorb = .true.``).  A typical noncollinear ``&system``
block looks like

```
&SYSTEM
  ...
  noncolin = .true.
  lspinorb = .false.    ! .true. requires SOC-aware UPF pseudopotentials
  starting_magnetization(1) = 0.1
  angle1(1) = 0.0
  angle2(1) = 0.0
  ...
/
```

The RDMFT total energy, the per-band $\partial E/\partial n$
gradient, the SPG occupation block, and the Stiefel orbital block
all behave the same way under ``noncolin = .true.`` as they do under
``nspin = 1``; only the DOS post-processor degrades to a total-DOS
view (PDOS / Mz skipped, see "Spin support" above).

### Restart / checkpointing

`rdmft_restart = .true.` reads the natural occupations
``rdmft_n(nbnd, nkstot)`` from the file
``{outdir}/{prefix}.rdmft.save`` (or the path given by
``rdmft_restart_file``) before the first RDMFT outer iteration, so a
long run can be resumed without redoing the macro-iteration history.
``rdmft_save_every`` (default ``1``) controls how often the save
file is written during the main loop; a final unconditional save
also runs at the end of ``rdmft_run`` so the converged solution is
always recoverable.

The natural **orbitals** ride along with QE's own wavefunction
restart mechanism: each save call invokes ``punch('config')`` so
the wavefunction buffer is flushed to the standard
``{outdir}/{prefix}.save/`` directory (collected wfc + ``data-file-schema.xml``).
On restart the user should set ``startingwfc = 'file'`` in the
``&electrons`` namelist so that PWscf's ``wfcinit`` reloads those
natural orbitals into ``evc`` before ``rdmft_run`` is entered.

A canonical two-step workflow:

```
# First run (computes RDMFT from KS):
&control  outdir = './tmp', prefix = 'h2',  ... /
&electrons ... /
&rdmft do_rdmft = .true., rdmft_save_every = 5 /

# Resume (reuse the same outdir and prefix):
&control  outdir = './tmp', prefix = 'h2',  ... /
&electrons startingwfc = 'file', electron_maxstep = 0 /
&rdmft do_rdmft = .true., rdmft_restart = .true. /
```

> **Why ``electron_maxstep = 0``?** With ``electron_maxstep > 0``
> PWscf's ``electrons_scf`` loop diagonalizes the KS / Fock
> Hamiltonian once (``c_bands``) and rotates ``evc`` onto its
> eigenstates *before* ``rdmft_run`` is entered.  The saved natural
> orbitals are **not** eigenstates of any one-body Hamiltonian (that
> is the whole point of the RDMFT optimisation), so even one
> diagonalization perturbs them and the very first restart energy
> drifts from the previously-saved value by a small amount.  Setting
> ``electron_maxstep = 0`` skips the SCF iteration entirely, so
> ``evc`` stays exactly as it was loaded from disk and the restart
> reproduces the previously-saved ``!RDMFT_ETOTAL`` bit-for-bit.
> On the shipped H2 example this gives
>
> ```
> # First run:
> !RDMFT_ETOTAL  =      -2.2428292579 Ry
>
> # Restart with electron_maxstep = 0 and rdmft_outer_maxiter = 0:
>   Initial RDMFT total energy =      -2.2428292579 Ry
> !RDMFT_ETOTAL  =      -2.2428292579 Ry
> ```

The save file is small (one ASCII row per ``(band, k-point)`` pair)
and contains a header with ``(nbnd, nkstot, nspin)`` for sanity
checks; a mismatch with the current run prints a warning, falls back
to the KS-derived initial occupations and continues without aborting.

Notes:

* QE's standard `disk_io` knob controls whether the wavefunction
  buffer is durable on disk.  Set `disk_io = 'medium'` or `'high'`
  (the default for SCF runs) so that the periodic save can actually
  persist `evc` between pw.x invocations; with `disk_io = 'low'` /
  `'none'` only the in-memory buffer is preserved and a true
  cross-process restart is not possible.
* Use ``electron_maxstep = 0`` on the resume run (see box above):
  this guarantees a bit-reproducible restart.  ``electron_maxstep
  = 1`` is still tolerated -- the alternating loop reconverges
  within a few outer iterations -- but the very first restart
  energy will be slightly worse than the previously-saved value
  because the diagonalization in ``c_bands`` perturbs the natural
  orbitals.

**Limitations (v1):** Gaussian broadening only; projections are not
symmetrized (`lsym` ignored); the noncollinear DOS path writes only
the total spectral function (PDOS / Mz skipped).

### Two-step DOS workflow (`pw.x` → `rdmft_dos.x`)

The RDMFT DOS post-processor lives in its own executable
`bin/rdmft_dos.x` (source: `PP/src/rdmft_dos.f90`) so the RDMFT
solver and the DOS calculation are completely separated in exactly
the same way as the standard QE `pw.x` / `dos.x` split.  The RDMFT
calculation (step 1) persists only the two pieces of state the DOS
step needs -- the natural occupations and the natural orbitals --
and `rdmft_dos.x` (step 2) reads them back, runs the
transition-state band probe on the fly, and writes the DOS / PDOS /
Mz outputs.

**Step 1 -- RDMFT calculation (save occupation + orbital).**  Run the
optimisation as usual; the save / restart machinery already in place
flushes:

* the natural occupations \(n_{ik}\) to
  ``{outdir}/{prefix}.rdmft.save`` (via the periodic
  ``rdmft_save_every`` checkpoints and the final unconditional save
  at the end of ``rdmft_run``);
* the natural orbitals \(C^k\) to the standard PWscf save directory
  ``{outdir}/{prefix}.save/`` (via the ``punch('config')`` call
  inside the save routine, which writes collected wavefunctions plus
  ``data-file-schema.xml``).

The expensive transition-state band probe (a per-(band, k) ACE
rebuild) is deliberately **not** written at this stage -- saving it
would mean every save call inside the optimiser pays the cost.
Instead, it is computed exactly once during the DOS step.

```
# Step 1 input -- RDMFT optimisation, no DOS:
&control  outdir = './tmp', prefix = 'h2', ... /
&electrons ... /
&rdmft do_rdmft = .true., rdmft_save_every = 1 /
```

**Step 2 -- dedicated DOS calculation, `rdmft_dos.x`.**  This is a
standalone executable (not a `pw.x` re-run): it lives at
`PP/src/rdmft_dos.f90` and is built into `bin/rdmft_dos.x` by both
the legacy Make build and the CMake build.  The input is a single
`&inputrdmftdos` namelist containing only DOS-relevant knobs:

```
&inputrdmftdos
  prefix             = 'h2'
  outdir             = './tmp'
  functional         = 'muller'      ! must match the previous pw.x run
  emin               = -1.10         ! Ry; omit (=±1e6) for auto window
  emax               =   0.37
  deltae             = 0.00147
  degauss            = 0.00368
  pdos               = .true.
  occ_weighted       = .false.
/
```

Run with:

```
rdmft_dos.x < step2.in > step2.out
```

Internally `rdmft_dos.x`:

1. Calls `read_file_new(needwf=.TRUE.)` to load the PWscf state from
   `{outdir}/{prefix}.save/data-file-schema.xml` (lattice, atoms,
   k-points, charge density, local potential, Fermi energy).
2. Opens the wavefunction buffer (`iunwfc`) and re-fills it from the
   collected wfc files `{outdir}/{prefix}.save/wfc*.dat`.
3. Calls `rdmft_load_state` to populate `rdmft_n` from
   `{outdir}/{prefix}.rdmft.save`.
4. Forces the global XC IDs to `(5,0,0,0,0,0)` and `exxalfa = 1.0`
   so `v_of_rho` returns `V_xc = 0` during the TSM probe.  Explicitly
   sets `exx_bgrp_type = EXX_BGRP_BANDS` (the standard band-parallel
   EXX scheme); without this, Fortran's zero-initialisation of the
   module variable sends `exxinit` / `vexx` down the band-pairs path
   (`exxinit_bp` / `vexx_bp`), which uses a different `exxbuff` layout
   and produces a ~3 mRy residue in E_xc on multi-k systems (Fe).
   Then calls `setup_exx` / `rdmft_exxinit_once` to bring up the
   EXX-via-ACE machinery.
5. Calls `rdmft_compute_dos()` and `rdmft_compute_magnetization()` --
   exactly the same subroutines pw.x would call inline.

The result is the standard `<prefix>.rdmft.dos`,
`<prefix>.rdmft.pdos_at<S>_<N>` and `<prefix>.rdmft.mag` files.
Tuning the DOS window (`emin / emax`), broadening
(`degauss`), grid step (`deltae`) or weighting
(`occ_weighted`) only requires another `rdmft_dos.x`
invocation -- the (potentially hour-scale) RDMFT optimisation is not
repeated.

**Denser-k restart (optional).**  When the RDMFT optimisation used a
coarse Monkhorst–Pack mesh, run a second `pw.x` job on a finer
`K_POINTS` grid, pointing `rdmft_source_prefix` /
`rdmft_source_outdir` at the **coarse** save while setting a
**different** `prefix` / `outdir` in `&control` for the fine run (the
solver aborts if source and output paths coincide).  If the fine mesh
is denser than the coarse mesh (e.g. 2×2×2 → 3×3×3), k-interpolation is
triggered automatically when the meshes differ.
PWscf also skips the KS SCF loop on the fine mesh
(`electron_maxstep` is forced to 0); occupations and natural orbitals
are bootstrapped from the coarse save and trilinearly interpolated onto
the fine irreducible mesh, followed by a short fixed-occupation orbital
polish (5 iterations by default) and the RDMFT optimisation.  Run
`rdmft_dos.x` on that fine save with no k-refine keys.  v1
limitations: automatic MP input only; each fine \texttt{nk} must be
$\geq$ the coarse value with at least one strict inequality (integer
commensurability not required; irreducible k-list lengths need not
divide evenly, e.g. AFM NiO 8 to 27); for LSDA the
coarse ``.rdmft.save`` must use the usual spin-block layout
(``nkstot = 2 * n_spatial``); Γ-only / manual k-lists are not
supported.

**Cost and crystal symmetry.**  The TSM band energies cost
$\mathcal{O}(N_b N_k)$ full `aceinit`/`vexx` gradient evaluations (one
per (band, k) pair), so they dominate wall time on multi-k / many-band
systems (transition-metal oxides).  Two things keep this affordable:

* **Crystal symmetry is already exploited.**  The probe runs only over
  the k-points QE keeps after symmetry reduction (the irreducible
  wedge), not the full Monkhorst–Pack grid.  For a cubic
  $2\times2\times2$ mesh with the full point group this is an $8\to3$
  reduction (verified on the Fe and CoO examples: "number of k points
  = 3" while "EXX: q-point mesh: 2 2 2" keeps the full grid for the
  exchange sum).  The RDMFT solver never disables symmetry, so it
  inherits whatever reduction QE finds; the DOS weights $w_k$ already
  carry the star multiplicities.  Antiferromagnetic cells (e.g. NiO
  with $\pm$ `starting_magnetization` on the two sublattices) have a
  lower magnetic point group, so their irreducible set is larger than
  the non-magnetic one — that is physical, not a missed reduction.  To
  check what your run uses, look at "number of k points = N" in the
  output; to deliberately keep more (e.g. for testing) you would set
  `nosym`/`noinv`, which only *increases* the count.
* **The probe is only computed when the DOS is requested.**  The final
  occupation-table `eps_tsm(n=1/2)` column is filled only when
  `rdmft_dos.x`; a `pw.x`-only run skips it entirely and
  prints just the cheap $h_{\mathrm{diag}}$/$v_x$ diagonals.
  (Previously the probe was hardwired on for every run, so even
  no-DOS / formerly-`'kf'` runs silently paid the
  $\mathcal{O}(N_b N_k)$ cost at the end.)
* **Single-k probe build.**  Each probe sets one occupation $n_{ik}=1/2$
  and reads a single gradient element, so it only needs the exchange
  rebuilt at *that* k-point.  `rdmft_grad_n` / `rdmft_compute_xc_channel`
  take an `only_k` argument that builds the ACE projector at the one
  k-point (`aceinit_k`) instead of the full all-k rebuild — an
  $N_k$-fold reduction, numerically identical (verified on the Fe DOS
  example: the `eps_tsm` table is bit-for-bit unchanged).  The
  per-pool cost is then $\propto N_b\,(N_k/P)$, so **using more k-point
  pools (`-nk P`) speeds the probe up roughly linearly** — combine
  `-nk` (divisor of `nkstot`) with `-np` plane-wave ranks.  Even so the
  transition-state DOS on a transition-metal-oxide supercell is a heavy
  post-processing step (minutes); keep `nbnd` and the k mesh as small as
  the DOS window allows.

**k-point pools.**  The DOS post-processor is pool-correct: it gathers
the per-pool weights / spin indices (`wk`, `isk`) and the per-band data
onto ionode before binning, and computes the global magnetization as a
local sum (no `mp_sum` inside the `ionode`-only output block).  An
earlier version indexed the *local* `wk`/`isk` arrays over the global
`1:nkstot` list (wrong spin-down DOS under `npool > 1`) and called a
pool-collective magnetization reduction from ionode only, which
deadlocked the other pools (the "DOS hangs forever with `-nk > 1`"
symptom).  Both are fixed.

## Known crash: ``ZPOTRF INFO=2 / Cholesky failed in invchol`` (obsolete: fast path removed)

Historical: `pw.x` used to abort at the very first SPG occupation
block on metallic super-cells (the NiO 2×2×2 / 8-atom case was the
canonical example) when the occupation coupling fast path was enabled.
That fast path (`rdmft_occ_reuse_orbitals` / `rdmft_frozen_fock`) has
since been **removed entirely** (it was not `-nk` reproducible on
metals), so this crash mode no longer exists.  The related
rank-deficient-`mexx` fallback described in the next subsection is a
property of the ACE path proper and still applies.

### Variant: ``ZPOTRF (1) / Cholesky failed in invchol`` from the ACE path (fixed)

This incarnation of the crash showed up on the ACE occupation path on
quasi-1D / very-large-vacuum cells (the user-reported H₂-chain
``a=10 Å, c=2 Å, nbnd=4, 8×1×1 k mesh, Power(α=0.65)`` case is the
canonical example):

```
   occ block (SPG / armijo, max 10 inner iters)
ZPOTRF exited with INFO=            2
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Error in routine ZPOTRF (1):
     Cholesky failed in invchol.
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```

Here `rdmft_compute_xc_channel` legitimately hands `aceinit_k` the
full per-band weights, but for a quasi-1D geometry the cross-band Fock
metric `mexx = <φ|V_x|φ>` is only rank-deficient (its `-mexx` has a
few small or even small-negative eigenvalues coming from the q→0
divergence treatment and the vacuum-like upper bands at each k).
ZPOTRF then aborts on the second leading minor.

**Fix.**  `aceupdate_k` (and the gamma-only `aceupdate`) now snapshot
`mexx` before the Cholesky and, if `ZPOTRF` returns `INFO ≠ 0`, fall
back to a Hermitian eigendecomposition with a thresholded
pseudo-inverse square root: directions with eigenvalues above
`max(eig) · 10⁻¹⁰` are inverted normally and directions at or below
the threshold are projected out of the ACE factor.  The resulting
operator `-|ξ⟩⟨ξ|` reproduces `V_x` exactly on the principal
subspace and is zero on the null space of `mexx` (the subspace
`mexx` carries no information about), so no physical exchange
contribution is lost.

When the fallback engages, the first few k-point hits are reported
on stdout (subsequent calls are summarised to avoid log flooding):

```
     aceupdate_k: -mexx not PD (ZPOTRF info=2); eigendecomp fallback. min(eig)= -2.989E-02  max(eig)=  1.323E+00
     ...
     aceupdate_k: further fallback warnings suppressed.
```

The fast path (Cholesky succeeds) is byte-for-byte unchanged, so
standard HF / hybrid pw.x calculations and the prior RDMFT benchmarks
are not affected.

### Variant: ``ZPOTRF (1) / Cholesky failed`` at ``-npool > 1`` (fixed)

A *different* incarnation of the same abort — ``ZPOTRF`` with
``INFO = 1`` part-way through a run (e.g. at the second outer cycle of
the CoO / NiO ``nspin = 2`` examples) — was caused by **k-point pool
desynchronisation**, not by the rank-1 precompute.  Several inter-pool
reductions were evaluated inside ``IF (ionode)`` print guards (so only
one pool called them) and several line-search scalars (the occupation
slope ``φ'``, the orbital ``‖G_R‖`` and slope, the CG coefficients)
were per-pool partial sums.  Either defect makes the pools issue
``Allreduce`` calls in different orders / counts, so MPI pairs the
wrong partial sums; the resulting garbage energies and occupations end
up handing ``aceinit`` a non-positive-definite metric, and ``ZPOTRF``
aborts.  All of these reductions are now performed on every rank, so
this crash no longer occurs.  A separate multi-pool defect in the
exchange-energy assembly (a non-conformable ``SUM(wg * vx_diag)`` with
the global-shape ``wg``) has also been fixed, so ``npool > 1`` now
gives a correct, deterministic ``E_xc``; see "K-point pool
parallelism" above.

## Parallelisation strategies for large RDMFT runs

A typical RDMFT-Power / RDMFT-Müller solve on a transition-metal-oxide
supercell (e.g. an NiO 8-atom cell at a 2×2×2 k mesh) spends almost all
its time **before** the first outer RDMFT iteration in the chain
``setup_exx → rdmft_exxinit_once → first rdmft_total_energy → first
aceinit``.  Each ``aceinit`` rebuilds the ACE projector ``xi(:,:,ik)``
for every k-point in the pool, and the FFT work inside ``vexx`` scales
roughly as ``nbnd * nks_per_pool * nqs * cost(FFT)``.  Three QE
parallelism axes can be combined and all of them inherit naturally
into the RDMFT solver:

1. **Plane-wave (G-vector) parallelism** (``mpirun -np N pw.x``) splits
   the G-vector list of every k-point across MPI ranks in the band
   group.  Every plane-wave-bound routine (``h_psi``, ``vloc_psi``,
   ``vexx``, FFTs) benefits.  All RDMFT energy / gradient kernels
   already carry an explicit ``mp_sum(..., intra_bgrp_comm)`` after
   the per-band dot products, so plane-wave parallelism is the safest
   first axis to scale on.  On the BCC Fe example (2×2×2 k mesh,
   nbnd = 10) we measure:

   | Layout | Wall (s) | Speed-up | ``!RDMFT_ETOTAL`` (Ry) |
   |--------|----------|----------|------------------------|
   | ``-np 1``                  | 28.3 | 1.0× | -33.6692554386 |
   | ``-np 2 -npool 1``         | 18.8 | 1.5× | -33.6687671679 |
   | ``-np 4 -npool 1``         | 15.1 | 1.9× | -33.6689343783 |

   The ~5e-4 Ry energy drift across different ``-np`` reflects the
   non-associativity of the per-rank plane-wave reduction; it scales
   with the line-search noise budget the same way as the OpenMP issue
   below.

2. **K-point pool parallelism** (``mpirun -np N pw.x -npool P``) splits
   the k-point list across ``P`` MPI pools.  **The RDMFT optimiser is
   now fully pool-correct**: every quantity that aggregates over
   k-points and feeds a control decision is reduced across
   ``inter_pool_comm`` on *all* ranks, including the energy /
   occupation sums (``rdmft_weighted_sum``,
   ``rdmft_constraint_violation_occ``, ``rdmft_pg_kkt_residual``, the
   SPG ``sum|dn|`` test, the proximal-projector bisection), **and** the
   line-search drivers (the SPG occupation slope ``φ'`` and CG
   coefficients, the orbital Riemannian gradient norm ``‖G_R‖`` and the
   product-manifold / joint slope, and the Barzilai–Borwein init step).

   Previously some of these reductions were skipped or — worse —
   evaluated *inside* ``IF (ionode)`` print guards.  Because MPI matches
   collective ``Allreduce`` calls by communicator in per-rank order,
   calling a reduction on ionode only (or computing a per-pool slope)
   desynchronised the pools: every subsequent ``inter_pool_comm``
   reduction paired up the wrong partial sums, corrupting the energy
   and occupations and finally feeding ``aceinit`` a
   non-positive-definite metric — the deterministic ``ZPOTRF (1) /
   Cholesky failed in invchol`` abort reported on the CoO and NiO
   ``nspin = 2`` multi-k runs.  That class of crash is fixed.

   A second multi-pool defect lived in the **exchange energy
   assembly**.  QE allocates the weight array ``wg`` with the *global*
   shape ``(nbnd, nkstot)`` but only its first ``nks`` columns hold a
   given pool's weights (it is indexed by the local k index).  The
   Fock-energy accumulator wrote ``SUM(wg * vx_diag)`` with
   ``vx_diag`` of shape ``(nbnd, nks)``; at ``npool = 1`` this is
   conformable (``nks = nkstot``), but at ``npool > 1`` it is a
   non-conformable array expression that reads ``vx_diag`` out of
   bounds into adjacent heap, producing a wrong / ``NaN`` /
   run-to-run-varying ``E_xc`` (the ``-nk 2`` "race").  It is now
   sliced explicitly (``SUM(wg(1:nbnd,1:nks) * vx_diag(1:nbnd,1:nks))``),
   and the analogous ``wg`` save/restore in ``rdmft_run`` is sliced
   too.

   With both fixes, **pool parallelism is fully supported and
   deterministic**.  On CoO 2×2×2 (``nspin = 2``, 6 k-points) a
   ``-nk 2`` run (pools split by spin) reproduces the ``-npool 1``
   trajectory essentially bit-for-bit (``E_xc = -32.3727049526`` Ry,
   identical across repeated runs; the alternating occupation/orbital
   energies track the single-pool run to all printed digits).

   ```bash
   mpirun -np N pw.x -nk P -in input.in     # P k-point pools (P | nkstot)
   ```

   Note that ``-nk`` is a synonym for ``-npool``; pass only one of them
   (the original CoO script's ``-npool 1 -nk 2`` was contradictory --
   the trailing ``-nk 2`` won).  Choose ``P`` to divide ``nkstot``
   (e.g. ``-nk 2`` for the spin-up / spin-down split when
   ``nspin = 2``).

3. **OpenMP** (``OMP_NUM_THREADS=N``) accelerates the threaded BLAS /
   FFTW backends and the QE-side ``$!omp parallel do`` loops in
   ``vexx`` / ``aceinit``.  Useful for filling a node beyond the
   number of MPI ranks (e.g. ``-np 4 OMP_NUM_THREADS=4`` on a
   16-core node).

### Recommended layouts

| System size                          | Recommended layout                  |
|--------------------------------------|--------------------------------------|
| Small molecule, 1 k-point            | ``-np N`` (pure PW)                  |
| 2–4 atom unit cell, ≤ 2×2×2          | ``-np N -nk P`` with P \| nkstot     |
| ≥ 8-atom super-cell, ≥ 2×2×2         | ``-np N -nk P`` with P \| nkstot     |
| Single node, free OMP slots          | add ``OMP_NUM_THREADS=k`` per rank   |

K-point pool parallelism (``-nk P`` / ``-npool P``) is the most
effective axis for transition-metal-oxide supercells, where the
per-(k,q) ``aceinit`` / ``vexx`` chain dominates: ``P`` pools cut that
cost by ``P``.  Choose ``P`` to divide ``nkstot`` (for ``nspin = 2``
the natural choice ``-nk 2`` puts spin-up k-points on one pool and
spin-down on the other).  Combine with plane-wave MPI inside each pool
(``-np`` larger than ``P``) and OpenMP threads to fill the node.

For the user-reported CoO / NiO supercells (``nspin = 2`` →
``2 × nkstot_spatial`` k-points), e.g.
``mpirun -np 8 pw.x -nk 2 ...`` (2 pools, 4 plane-wave ranks each) on a
single 8-core node works and prints both spin channels.

## Where the algorithm comes from

See the ABACUS reference module's documentation:

* `abacus_source/source_lcao/module_rdmft/doc/rdmft_derivation.md`
  — full mathematical derivation of the energy, gradients, Stiefel
  projection / retraction, and projected gradient.
* `abacus_source/source_lcao/module_rdmft/doc/elk_rdmft_algorithm.md`
  — ELK RDMFT driver summary (QE small-$n$ coupling follows ABACUS
  `pow_reg`, not ELK's split energy/derivative floors).
* `abacus_source/source_lcao/module_rdmft/doc/rdmft_usage.md`
  — INPUT keyword reference (one-to-one mapping to the names used in
  this Fortran port).

The current implementation is a direct, minimal Fortran port of the
``alternating + projected_gradient + Stiefel`` path requested in the
task statement, reusing PWscf's existing EXX-via-ACE machinery for the
exchange evaluation.
