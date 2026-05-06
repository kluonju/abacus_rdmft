# RDMFT User Guide

This document explains how to set up and run a Reduced Density Matrix Functional
Theory (RDMFT) calculation in ABACUS.

---

## Quick-start example

A minimal `INPUT` file for an RDMFT calculation on top of a Kohn–Sham starting
point with the Müller functional:

```
INPUT_PARAMETERS
calculation         scf
basis_type          lcao
dft_functional      pbe
rdmft_functional     muller
rdmft               1
```

ABACUS will first run a standard KS-SCF loop using `dft_functional` (e.g. PBE
for a good starting guess), then hand the converged orbitals and occupation
numbers to the RDMFT optimiser.  The RDMFT stage uses `rdmft_functional`.  All
RDMFT-specific keywords below are optional and have sensible defaults.

### k-point sampling

Both `gamma_only = 1` (gamma-only LCAO) and `gamma_only = 0` (multi-k LCAO,
optionally with a single Γ point) are supported and produce RDMFT energies
that match the KS reference to machine precision at the converged HF orbital.
On the shipped `H2_HF` and `LiH_HF` examples we verify
`E_RDMFT_HF == E_KS_HF` byte-for-byte at the converged HF point.

---

## Complete keyword reference

### Core switches

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft` | bool | `false` | Master switch. Set to `1` (or `true`) to enable RDMFT. |
| `rdmft_functional` | string | `""` (empty) | RDMFT XC functional. Supported: `hf`, `muller`, `power`, `gu`, `chf`, `cga`, `geo`, `optgm`. **When left empty only the legacy single-step RDMFT evaluator is run — the RDMFT energy is NOT minimised. Set this keyword to activate the full RDMFT optimisation engine.** |
| `rdmft_power_alpha` | real | `0.656` | Exponent α in the coupling function g(n) = n^α. Only used when `rdmft_functional` is `power`. For HF and Müller, α is fixed automatically (1.0 and 0.5). Valid range: (0, 1). |

The **initial KS-SCF** uses `dft_functional` (LDA, PBE, SCAN, etc.).  The **RDMFT optimisation** uses
`rdmft_functional`.  Supported values:

| `rdmft_functional` | g(n) | Notes |
|--------------------|------|-------|
| `hf` | n | Hartree–Fock, α forced to 1.0 |
| `muller` | √n | Müller / BBC1, α forced to 0.5 |
| `power` | n^α | Power functional, α from `rdmft_power_alpha` |
| `gu` | √(n_i)√(n_j) off-diag, n_i² diagonal | Goedecker–Umrigar with self-interaction correction |
| `chf` | 1/2 n_i n_j + 1/2 √(n_i(1−n_i))√(n_j(1−n_j)) | Corrected Hartree–Fock (Csanyi–Arias) |
| `cga` | 1/4 n_i n_j + 1/4 √(n_i(2−n_i))√(n_j(2−n_j)) | Csanyi–Goedecker–Arias (PRA 65, 032510) |
| `geo` | [n_i n_j + (n_i n_j)^(1/2) + 2(n_i n_j)^(3/4)] / 4 | GEO mixture functional |
| `optgm` | (1-λ) n_i n_j + λ n_i^α n_j^α | OptGM fixed-parameter convex mixture |

### Solver strategy

| Keyword | Type | Default | Allowed values |
|---------|------|---------|----------------|
| `rdmft_solver_strategy` | string | `alternating` | `alternating`, `joint` |

- **`alternating`** — Alternate between optimising occupation numbers (orbitals
  fixed) and orbitals (occupations fixed).  This is the standard approach and
  the most robust.
- **`joint`** — Pack occupations (Euclidean space, after parameterisation)
  and orbitals (Stiefel manifold) into a single point on the product manifold
  and optimise everything simultaneously.  The packed variable
  `z = (p, C^1, ..., C^{Nk})` is optimised by **one** unified optimiser,
  configured via `rdmft_joint_optimizer`.  The optimiser consumes the packed
  gradient `g = (dE/dp, G_R^1, ..., G_R^{Nk})`, where `dE/dp` is the
  chain-rule transformation of `dE/dn` through the occupation parameterisation
  and each `G_R^k` is the Riemannian (tangent-space) projection of the
  Euclidean orbital gradient at `C^k`.  The produced search direction is
  re-projected onto each tangent space and a single Armijo line search runs
  along the packed direction (linear update for the parameters, Stiefel
  retraction for each `C^k`).  `rdmft_occ_optimizer` and `rdmft_orb_optimizer`
  are ignored by this strategy.  The legacy value `product_manifold` is still
  accepted as an alias for `joint`.

### Occupation parameterisation

| Keyword | Type | Default | Allowed values |
|---------|------|---------|----------------|
| `rdmft_occ_param` | string | `cosine_sq` | `cosine_sq`, `logistic`, `sigma_shift` |

The first two options map unconstrained real parameters to the interval [0, 1]:

- **`cosine_sq`** — n = cos²(θ).  The gradient transforms as
  dE/dθ = −sin(2θ) · dE/dn.
- **`logistic`** — per-state sigmoid n = σ(x) = 1/(1+e^{−x}) with independent
  parameters x. With `rdmft_constraint = augmented_lagrangian`, the
  electron-number constraint is handled by the augmented Lagrangian (not a
  global μ solve). Projected gradient and active set work in occupation space
  and do not use this map.
- **`sigma_shift`** — \(n_{ik}=\sigma(z_{ik}+\lambda)\) with a **scalar**
  \(\lambda\) solved each evaluation so that \(\sum_k w_k\sum_i n_{ik}=N_e\)
  exactly (bisection).  This is the unconstrained occupation block of the
  **product-manifold** setup and is wired only in **`rdmft_solver_strategy joint`**.
  Using `sigma_shift` with `alternating` is rejected at solver initialisation.
  The joint driver skips the augmented-Lagrangian occupation penalty when
  `sigma_shift` is active; `rdmft_constraint` still selects the method for
  **alternating** runs (not applicable together with `sigma_shift`).

### Electron-number constraint

| Keyword | Type | Default | Allowed values |
|---------|------|---------|----------------|
| `rdmft_constraint` | string | `augmented_lagrangian` | `augmented_lagrangian`, `projected_gradient`, `active_set` |

The total electron number must satisfy Σ_k w_k Σ_i n_{ik} = N_e.

- **`augmented_lagrangian`** — Augmented Lagrangian on occupation parameters
  (`rdmft_occ_param` = `cosine_sq` or `logistic`), Armijo backtracking, and
  automatic updates of λ and μ.
- **`projected_gradient`** — Occupations are **physical** band weights
  \(n_{ik}\in[0,1]\) (not the `cosine_sq` / `logistic` maps).  Each inner
  iteration uses the Euclidean `rdmft_occ_optimizer` direction in \(n\)-space,
  a **projected Armijo** line search onto the feasible set (clip to \([0,1]\)
  then dual rescaling so \(\sum_k w_k\sum_i n_{ik}=N_e\)), and a **Bertsekas
  projected-gradient map** for convergence; see
  [Projected gradient details](#projected-gradient-details) below.
- **`active_set`** — Legacy reduced equality-constrained solver on the free
  occupations. It remains available, but it is not one of the main recommended
  paths above.

<a id="projected-gradient-details"></a>
### Projected gradient details (`rdmft_constraint = projected_gradient`)

PG works **directly in occupation space** on the convex set
\(0\le n_{ik}\le 1\) with the linear equality \(\sum_k w_k\sum_i n_{ik}=N_e\).
The `rdmft_occ_param` keyword is ignored for this path (ALM uses `cosine_sq` /
`logistic`; `sigma_shift` is only used with `joint`).

**Feasible set and projection.**  After a trial move \(n+\alpha d\), the code
clips each \(n_{ik}\) to \([0,1]\) and applies a **one-dimensional dual solve**
(rescaling with a per–\(k\)-weight shift \(\lambda\)) so the weighted electron
sum matches \(N_e\).  If the bracket for \(\lambda\) fails (too many states
pinned at 0 or 1), the trial is rejected and the line search shrinks \(\alpha\).

**Search direction.**  The gradient \(\partial E/\partial n\) is taken at fixed
orbitals.  The same Euclidean optimisers as for ALM (`sd`, `cg`, `lbfgs`, `adam`)
build a descent direction \(d\); if \(d^\top \nabla_n E\ge 0\), the code falls
back to **steepest descent**, \(d = -\nabla_n E\).

**Line search (projected Armijo).**  For trial \(n' = P(n+\alpha d)\) with
feasible projection \(P\):

1. If the equality residual after projection is too large, reduce \(\alpha\)
   (no energy evaluation).
2. Let \(\Delta n = n'-n\).  If \(\nabla_n E\cdot \Delta n \ge 0\), reduce
   \(\alpha\) (the projected segment is not a first-order descent at linear
   order).
3. Otherwise evaluate \(E(n')\) and accept the step if
   \(E(n') \le E(n) + c_1\,(\nabla_n E\cdot \Delta n)\) with Armijo constant
   \(c_1 =\) `rdmft_line_search_c1` (same keyword as for orbital Armijo and AS
   occupations).
4. On rejection, multiply \(\alpha\) by a fixed factor (default \(0.5\) in the
   solver) until the iteration cap is hit.

If **no** trial in that backtracking loop is accepted, the solver applies a
**single projected steepest step**
\(n \leftarrow P\bigl(n - \tau\,\nabla_n E\bigr)\) with \(\tau =\)
`rdmft_alpha_step`, then **resets** the occupation optimiser state (CG / L-BFGS
/ Adam history).  This mirrors the recovery used in the `active_set` path when
its line search fails.

**Convergence criterion (inner).**  Define the **projected-gradient** (Bertsekas)
\(\mathbf{g}_{\mathrm{proj}} = \tfrac{1}{\tau}\bigl(\mathbf{n} - P(\mathbf{n} -
\tau \nabla_{\mathbf{n}}E)\bigr)\) with **\(\tau\) tied to the occupation line
search**: before stepping, \(\tau\) equals the first Armijo trial \(\alpha_0\)
from `rdmft_occ_ls_init_step` (with fallback to `rdmft_alpha_step` if \(\alpha_0\)
is invalid); after a successful line search, \(\tau\) is the **accepted** Armijo
step \(\alpha\); after an SD fallback (failed line search), \(\tau =\)
`rdmft_alpha_step`, matching the recovery step.  The PG inner loop stops only
when \(\|\mathbf{g}_{\mathrm{proj}}\|_\infty <\) `rdmft_occ_grad_tol` (evaluated
**post-step**; the log reports the unscaled map \(\|\mathbf{n} - P(\cdot)\|_\infty\)
and the \(\tau\) used).  The first inner may exit early if the **pre-step**
\(\|\mathbf{g}_{\mathrm{proj}}\|_\infty\) is already below the tolerance.  INPUT
`rdmft_occ_energy_tol` is **not** used for PG stopping.

**Outer (alternating / joint).**  A run is not considered **globally** converged
on `rdmft_energy_tol` alone: the outer loop requires (after the first cycle)
**both** inner sub-problem convergence flags and, when `rdmft_energy_tol` \(>0\), a
small \(|E - E_{\mathrm{prev}}|\) between outer iterations.  If
`rdmft_energy_tol` \(\le 0\), the outer energy test is off and the loop stops
when the inner flags alone indicate convergence.

**When to use PG.**  Useful for experiments or when you want updates purely in
\(n\)-space without ALM penalties.  For difficult functionals (e.g. Müller at
stretched geometries), **ALM (`augmented_lagrangian`) is usually more robust**.
If the log shows repeated `PG line search failed` / `applied SD fallback`, try
smaller `rdmft_alpha_step`, looser `rdmft_line_search_c1`, more inner iterations
`rdmft_occ_maxiter`, or switch constraint method.

### Optimiser selection

| Keyword | Type | Default | Allowed values |
|---------|------|---------|----------------|
| `rdmft_occ_optimizer` | string | `cg` | `sd`, `cg`, `lbfgs`, `adam` |
| `rdmft_orb_optimizer` | string | `cg` | `sd`, `cg`, `lbfgs`, `adam` |
| `rdmft_joint_optimizer` | string | `lbfgs` | `sd`, `cg`, `lbfgs`, `adam` |

`rdmft_occ_optimizer` and `rdmft_orb_optimizer` control the two sub-problem
optimisers used by the `alternating` strategy.  `rdmft_joint_optimizer`
selects the **single** unified optimiser used by the `joint` strategy on the
packed `(p, C)` variable.  Only the keyword matching the active strategy is
consulted; the others are ignored.

You can choose different optimisers for the occupation and orbital
sub-problems:

| Value | Algorithm | Notes |
|-------|-----------|-------|
| `sd` | Steepest descent | Most robust, slowest convergence |
| `cg` | Conjugate gradient (Polak–Ribière / Fletcher–Reeves with Powell restart) | Good balance of speed and reliability. For **orbitals**, if an Armijo line search fails, the next inner iteration restarts along steepest descent (-G_R). |
| `lbfgs` | limited-memory lbfgs | Fast for smooth landscapes, uses `rdmft_lbfgs_memory` history vectors. For the orbital sub-problem, the quasi-Newton direction is projected back onto the Stiefel tangent space (Riemannian lbfgs by projection). |
| `adam` | Adam | Adaptive learning rate, useful for noisy or ill-conditioned problems. For the orbital sub-problem Adam's Euclidean update is projected onto the tangent space and retracted onto the Stiefel manifold at each step. |

All four optimisers are available for both `rdmft_occ_optimizer` and
`rdmft_orb_optimizer` and may be mixed freely (e.g. `cg` on occupations and
`lbfgs` on orbitals).

### Convergence control

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_outer_maxiter` | int | `200` | Maximum **outer** RDMFT cycles: each cycle is one occupation optimisation plus one orbital optimisation (`alternating`), or one joint product-manifold step (`joint`). |
| `rdmft_occ_maxiter` | int | `50` | Maximum **inner** iterations for the occupation sub-problem (orbitals fixed) within one outer cycle. |
| `rdmft_orb_maxiter` | int | `50` | Maximum **inner** iterations for the orbital sub-problem (occupations fixed) within one outer cycle. |
| `rdmft_energy_tol` | real | `1e-8` | Outer: require **smaller** than this for \(\|E - E_{\mathrm{prev}}\|\) (Ry) **and** both occupation- and orbital-inner `converged` flags, after the first outer step.  Set `<=0` to disable the energy part (stopping uses inner flags only). |
| `rdmft_orb_grad_tol` | real | `1e-6` | Alternating orbital inner loop: stop when Riemannian gradient norm `||G_R||` is below this. |
| `rdmft_orb_energy_tol` | real | `1e-8` | Alternating orbital inner loop: when `> 0`, stopping requires **both** `||G_R|| <` `rdmft_orb_grad_tol` **and** small energy moves: `|E_k - E_{k-1}|` before the step and `|E_{\mathrm{new}} - E|` after an accepted line search must stay below this (Ry). Set `<= 0` for gradient-only stopping. |
| `rdmft_occ_tol` | real | `1e-7` | **ALM / active_set:** tolerance on occupation inner updates (e.g. sum \(|\Delta n|\) and KKT-style checks for AS). **Not** the PG stopping rule. |
| `rdmft_occ_grad_tol` | real | `1e-6` | **PG:** stop when \(\|\frac{1}{\tau}(n-P(n-\tau\nabla_n E))\|_\infty <\) this, with \(\tau\) from the occupation line search (initial \(\alpha_0\), accepted \(\alpha\), or `rdmft_alpha_step` after SD fallback—see projected-gradient details).  Also used for other occupation-gradient checks (e.g. ALM diagnostics) as in the code. |
| `rdmft_occ_energy_tol` | real | `1e-8` | **Not used** for `projected_gradient` (PG inner stop is `rdmft_occ_grad_tol` vs \(\|g_{\mathrm{proj}}\|_\infty\) only). Kept for backward-compatible INPUT. |

### Initial occupation setup

Before RDMFT optimisation starts, occupations are initialised uniformly as
one of the modes below.

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_occ_init_mode` | string | `ks` | Initial occupation mode: `ks`, `perturbed`, `uniform`. |
| `rdmft_occ_init_perturb` | real | `0.0` | Perturbation magnitude `delta` used by `rdmft_occ_init_mode = perturbed`. |
| `rdmft_occ_init_nbands_top` | int | `0` | Fermi-window half-width `K` used by `perturbed`/`uniform`: select `K` bands above and `K` bands below the Fermi boundary (per k-point). |

- `ks`: use KS occupations directly as RDMFT initial occupations.
- `perturbed`: for each k-point and selected Fermi-window bands, apply `+delta` to bands above Fermi and `-delta` to bands below Fermi, then project back to the feasible set.
- `uniform`: for each k-point, compute `N_top` over the selected `2K` Fermi-window bands and set each selected occupation to `N_top / (2K)` (or `N_top / N_selected` near band edges), then project.
- If `K <= 0`, `perturbed` and `uniform` fall back to `ks`.

### Line search and optimiser tuning

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_alpha_step` | real | `1.0` | Default / fallback trial step for RDMFT line searches. Orbitals (Armijo) use it for non–quasi-Newton optimisers. For **PG** occupations it is the fallback if the line-search initial \(\alpha_0\) is non-positive / non-finite, the step length in the **SD fallback** after a failed line search, and (matching that move) the post-step Bertsekas \(\tau\) when fallback runs; the usual post-step \(\tau\) is the **accepted** Armijo \(\alpha\). For ALM it is also the fallback when `rdmft_occ_ls_init_step` does not supply \(\alpha_0\). |
| `rdmft_occ_ls_init_step` | string | `bb` | Occupation line-search **first trial** step: `fixed` (always 1.0), `bb` (Barzilai–Borwein estimate; uses `rdmft_alm_bb_*` clamps), or `quad` (quadratic interpolation from the previous inner iteration’s first energy trial). Used for `projected_gradient` and `active_set`. |
| `rdmft_line_search_c1` | real | `1e-4` | Armijo sufficient-decrease \(c_1\) in \((0,1)\): smaller is stricter. Used for alternating **orbital** Armijo, ALM / PG / AS **occupation** line searches, and joint Armijo; also the Armijo side of **Strong Wolfe** when \(c_2\) is used. |
| `rdmft_line_search_c2` | real | `0.9` | Strong Wolfe curvature \(c_2\): used when **Strong Wolfe** runs (`rdmft_occ_optimizer = lbfgs` with ALM, or `rdmft_joint_optimizer = lbfgs` in joint mode). Require \(\lvert g^\top d\rvert \le c_2 \lvert g_0^\top d\rvert\). |
| `rdmft_line_search_max_zoom` | int | `20` | Maximum **zoom** iterations in Strong Wolfe (Nocedal & Wright). |
| `rdmft_lbfgs_memory` | int | `10` | Number of past gradient/step pairs stored by lbfgs. |
| `rdmft_adam_lr` | real | `0.001` | Learning rate for the Adam optimiser. |
| `rdmft_alm_lambda_init` | real | `0.0` | Initial ALM Lagrange multiplier `lambda` (only for `rdmft_constraint = augmented_lagrangian`). |
| `rdmft_alm_mu_init` | real | `1.0` | Initial ALM penalty parameter `mu` (only for `rdmft_constraint = augmented_lagrangian`). |
| `rdmft_alm_mu_factor` | real | `2.0` | Multiplicative ALM penalty update factor: `mu <- min(mu * factor, mu_max)`. |

`rdmft_alm_bb_enabled`, `rdmft_alm_bb_mode`, `rdmft_alm_bb_alpha_min`, and `rdmft_alm_bb_alpha_max` control Barzilai–Borwein step estimates in occupation **parameter** space for the augmented-Lagrangian path (when `rdmft_alm_bb_enabled` is true) and bound the BB branch of `rdmft_occ_ls_init_step = bb`. **Strong Wolfe** (bracket + zoom) is used for lbfgs in ALM occupations and for the **joint** strategy when `rdmft_joint_optimizer = lbfgs`; otherwise **Armijo** applies with `rdmft_line_search_c1`.  **ALM** occupation Armijo may use optional **polynomial** suggestions when `rdmft_line_search_polynomial` is true; **PG** and **active_set** occupation line searches use **geometric** backtracking only (same \(c_1\), no polynomial branch). Alternating **orbital** optimisation uses Armijo only (no polynomial there). `rdmft_occ_ls_init_step` is not read for orbital optimisation.

### Debugging

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_grad_check` | bool | `false` | If `true`, run a finite-difference gradient verification **at the start of the RDMFT solve**, immediately after the LCAO coefficients are converted to internal **X-space** (`precompute_cholesky_S` / `wfc_C_to_X`) and before the main alternating or joint loop. Occupations use a central difference in `n`; the **orbital** check uses the projected Riemannian gradient `G_R` and a **forward** difference along the same polar retraction as the line search `(E(t)-E_0)/t`, which matches the Armijo slope `-||G_R||^2` more reliably than a symmetric `±t` probe. This is expensive but essential for development and validation. Results are written to the running log. |

---

## Worked examples

### Example 1: Müller functional for a semiconductor

```
INPUT_PARAMETERS
calculation         scf
basis_type          lcao
dft_functional      pbe
rdmft_functional     muller
rdmft               1
rdmft_solver_strategy   alternating
rdmft_occ_optimizer     cg
rdmft_orb_optimizer     cg
rdmft_outer_maxiter        100
rdmft_energy_tol        1e-7
```

This runs a KS-SCF with PBE to convergence, then optimises natural occupations
(CG with augmented-Lagrangian constraint) and natural orbitals (CG on the Stiefel
manifold) in an alternating fashion using the Müller functional.

### Example 2: Power functional with lbfgs and gradient check

```
INPUT_PARAMETERS
calculation         scf
basis_type          lcao
dft_functional      pbe
rdmft_functional     power
rdmft               1
rdmft_power_alpha       0.55
rdmft_solver_strategy   alternating
rdmft_occ_optimizer     lbfgs
rdmft_orb_optimizer     lbfgs
rdmft_lbfgs_memory      20
rdmft_outer_maxiter        200
rdmft_energy_tol        1e-8
rdmft_grad_check        1
```

The `rdmft_grad_check 1` line triggers a finite-difference comparison of the
analytic gradients at the KS starting point.  Look for the "Gradient
Consistency Check" section in the running log to verify that all relative errors
are small (typically < 1e-4).

### Example 3: Joint (product-manifold) optimisation with logistic parameterisation

```
INPUT_PARAMETERS
calculation         scf
basis_type          lcao
dft_functional      pbe
rdmft_functional     hf
rdmft               1
rdmft_solver_strategy   joint
rdmft_occ_param         logistic
rdmft_constraint        augmented_lagrangian
rdmft_joint_optimizer   lbfgs
rdmft_alpha_step         0.01
rdmft_outer_maxiter        300
```

The explicit `rdmft_alpha_step 0.01` is smaller than the default `1.0` and
can stabilise joint HF steps; omit it to use the default.

Here occupations and orbitals are optimised simultaneously on the product
manifold (new `joint` keyword; `product_manifold` is still accepted as a
deprecated alias).  The packed variable `(p, C)` is handled by a single
lbfgs optimiser (`rdmft_joint_optimizer lbfgs`): every outer iteration the
solver evaluates the packed gradient `(dE/dp, G_R)`, asks lbfgs for a
descent direction, re-projects the orbital block onto the Stiefel tangent
space, and runs one Armijo line search along the packed direction (linear
update for the occupation parameters, Stiefel retraction for the orbitals).

### Example 4: Gamma-only solid with projected gradient

```
INPUT_PARAMETERS
calculation         scf
basis_type          lcao
gamma_only          1
dft_functional      pbe
rdmft_functional     power
rdmft               1
rdmft_power_alpha       0.656
rdmft_constraint        projected_gradient
rdmft_occ_optimizer     sd
rdmft_orb_optimizer     cg
rdmft_outer_maxiter        150
```

For gamma-only calculations all matrices are real (`double`), which is
automatically detected.  With `projected_gradient`, occupations stay in
\([0,1]\) with \(\sum_k w_k n = N_e\) after each accepted PG step or SD
fallback (clip + dual rescaling); see
[Projected gradient details](#projected-gradient-details).

---

## How the calculation proceeds

1. ABACUS runs a standard KS-DFT self-consistent field (SCF) calculation using
   `dft_functional`.  The number of SCF steps is controlled by `scf_nmax` as usual.

2. After SCF convergence (or after `scf_nmax` steps), the converged KS
   orbitals and occupation numbers are passed to the RDMFT module as the
   initial guess for the natural orbitals and natural occupation numbers.

3. The RDMFT module then minimises the RDMFT total energy functional (using
   `rdmft_functional`)
   E[{n_i}, {φ_i}] with respect to the natural occupation numbers {n_i}
   and natural orbital coefficients {C^k}, subject to:
   - 0 ≤ n_i ≤ 1  (N-representability, enforced by parameterisation)
   - Σ_k w_k Σ_i n_{ik} = N_e  (electron number, enforced by constraint method)
   - (C^k)† S^k C^k = I  (orbital orthonormality, enforced by Stiefel manifold)

4. The optimised RDMFT total energy is written back into the ABACUS energy
   record.  Energy components are printed in the running log with labels
   `E_TV_RDMFT`, `E_hartree_RDMFT`, `Exc_*_RDMFT`, and `Etotal_RDMFT`.

---

## Tips

- **Three iteration limits.** `rdmft_outer_maxiter` limits the **outer** RDMFT loop (each cycle: occupation step + orbital step in `alternating`, or one joint step in `joint`/`product_manifold`). `rdmft_occ_maxiter` and `rdmft_orb_maxiter` limit the **inner** optimisations of occupations (fixed orbitals) and orbitals (fixed occupations) within each outer cycle of the alternating strategy; they are unused for the `joint` strategy (which does a single joint step per outer iteration). In older versions, `rdmft_orb_maxiter` incorrectly doubled as the outer-loop limit; use `rdmft_outer_maxiter` for that now.

- **Always set `rdmft_functional`.** Without it only the legacy single-step
  RDMFT evaluator runs, which simply reports the RDMFT energy at the KS
  solution; the occupations and orbitals are **not** optimised.
  This is the most common cause of "I enabled RDMFT but nothing
  happens" reports.

- **Start with a good KS guess.** Use `dft_functional` (e.g. PBE, LDA) for the
  initial KS-SCF.  A well-converged KS-DFT starting point greatly accelerates
  RDMFT convergence.

- **Do NOT set `dft_functional` to a hybrid functional just to run RDMFT.**
  The RDMFT XC functional is controlled by `rdmft_functional`, not
  `dft_functional`.  Setting `dft_functional = muller` (or `hf`, `pbe0`, ...)
  forces the KS-SCF itself to do exact exchange, which is significantly
  more expensive (≈1–2 min of LibRI initialisation plus EXX in every SCF
  step for a small cell) and brings no benefit over a cheap LDA/PBE
  starting guess followed by a muller RDMFT optimisation.

- **Check gradients during development.** Set `rdmft_grad_check 1` when
  implementing or testing new functionals.  If the analytic and finite-difference
  gradients disagree, the optimisation will not converge.

- **Augmented Lagrangian is usually best** for the electron-number constraint,
  especially combined with cosine-squared parameterisation.

- **Projected gradient** ignores `rdmft_occ_param` and optimises raw \(n_{ik}\).
  If PG line searches fail often, reduce `rdmft_alpha_step`, relax
  `rdmft_line_search_c1`, or prefer ALM.

- **CG or lbfgs** are the recommended optimisers for production.  Steepest
  descent is useful for debugging.  Adam can help with difficult convergence
  but may require tuning `rdmft_adam_lr`.

- **For periodic systems with multiple k-points**, the constraint
  Σ_k w_k Σ_i n_{ik} = N_e is automatically handled using the k-point weights.

---

## Troubleshooting

### PG occupation line search always fails or resets the optimiser

The running log may show `occ line search (PG Armijo): ... fail` and
`applied SD fallback and reset optimizer` when no trial satisfies projected
Armijo within the backtracking cap.  Typical mitigations: smaller
`rdmft_alpha_step`, slightly **larger** `rdmft_line_search_c1` (e.g. `1e-3`–`1e-2`),
more `rdmft_occ_maxiter`, or `rdmft_occ_optimizer sd`.  For stubborn cases switch
to `rdmft_constraint augmented_lagrangian`.  If failures persist with plausible
inputs, run `rdmft_grad_check 1` to verify \(\partial E/\partial n\) against
finite differences.

### ABACUS runs for a long time with no RDMFT output

Typical causes (in roughly decreasing likelihood):

1. **`rdmft_functional` is not set.** Only the legacy single-step evaluator
   runs and no `===== RDMFT Alternating Optimization =====` banner is printed.
   Set `rdmft_functional muller` (or `hf`/`power`/`gu`/`chf`/`cga`/`geo`/`optgm`) to activate the
   optimisation engine.
2. **`dft_functional` is set to a hybrid (e.g. `muller`, `hf`, `pbe0`).**
   The KS-SCF itself then runs with exact exchange, and LibRI initialisation
   alone takes ~1–2 minutes of CPU time before the first SCF step even
   starts. Remove the `dft_functional` line (or set it to `pbe`/`lda`) to
   fall back to a cheap semi-local KS starting guess.
3. **ABACUS was built without RDMFT support.** Re-configure with
   `-DENABLE_RDMFT=ON` and rebuild; otherwise `INPUT: rdmft=true` triggers a
   `WARNING_QUIT`.

### The log shows `XC_fun: default` and a second `Etotal_RDMFT` at the end

This happened in older builds where the legacy single-step
`RDMFT::run()` was always invoked after the new solver, overwriting its
output.  Current builds guard the legacy path on
`rdmft_functional == ""`, so when the new engine is active only its
`FINAL_ETOT_IS` is reported.  If you still see this, make sure your
build includes this guard in `source/source_io/module_ctrl/ctrl_scf_lcao.cpp`.

### The example in `source/source_lcao/module_rdmft/example/ZnO_RDMFT`

The shipped `INPUT` uses the pseudopotential's native LDA for the KS-SCF
and `rdmft_functional muller` for the RDMFT stage.  On a single MPI rank
the first ≈2 min are spent initialising LibRI for the RDMFT Fock operator,
then the KS-SCF converges in ~25 iterations and the RDMFT alternating
CG optimisation runs for `rdmft_outer_maxiter` outer steps, printing
intermediate occupation and orbital sub-problem progress before writing
the final `Etotal_RDMFT` and `!FINAL_ETOT_IS` lines.
