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

---

## Complete keyword reference

### Core switches

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft` | bool | `false` | Master switch. Set to `1` (or `true`) to enable RDMFT. |
| `rdmft_functional` | string | `""` (empty) | RDMFT XC functional. Supported: `hf`, `muller`, `power`, `gu`. **When left empty only the legacy single-step RDMFT evaluator is run — the RDMFT energy is NOT minimised. Set this keyword to activate the full RDMFT optimisation engine.** |
| `rdmft_power_alpha` | real | `0.656` | Exponent α in the coupling function g(n) = n^α. Only used when `rdmft_functional` is `power`. For HF and Müller, α is fixed automatically (1.0 and 0.5). Valid range: (0, 1). |

The **initial KS-SCF** uses `dft_functional` (LDA, PBE, SCAN, etc.).  The **RDMFT optimisation** uses
`rdmft_functional`.  Supported values:

| `rdmft_functional` | g(n) | Notes |
|--------------------|------|-------|
| `hf` | n | Hartree–Fock, α forced to 1.0 |
| `muller` | √n | Müller / BBC1, α forced to 0.5 |
| `power` | n^α | Power functional, α from `rdmft_power_alpha` |
| `gu` | √(n_i)√(n_j) off-diag, n_i² diagonal | Goedecker–Umrigar with self-interaction correction |

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
  iteration is a **Spectral Projected Gradient** (Birgin–Martínez–Raydan,
  SIAM J. Optim. 10 (2000) 1196) step with a Barzilai–Borwein spectral step
  length and a **monotone** Armijo line search along the spectral projected
  direction; see
  [Projected gradient details](#projected-gradient-details) below.
- **`active_set`** — Alias for `projected_gradient` in the current
  implementation: both route to the same SPG inner loop.  Kept for
  backward INPUT compatibility.

<a id="projected-gradient-details"></a>
### Projected gradient details (`rdmft_constraint = projected_gradient`)

PG (and its alias `active_set`) is the **Spectral Projected Gradient (SPG)**
method of Birgin–Martínez–Raydan (SIAM J. Optim. 10 (2000) 1196), with a
**monotone** Armijo line search along the spectral projected direction
(sufficient decrease vs the current energy \(E(\mathbf{n}_k)\)).  It works **directly in occupation space** on the
convex set \(0\le n_{ik}\le 1\) with the linear equality
\(\sum_k w_k\sum_i n_{ik}=N_e\); the `rdmft_occ_param` keyword is ignored
on this path (ALM uses `cosine_sq` / `logistic`; `sigma_shift` is only used
with `joint`).

**Per inner iteration** \(k\):

1. Compute \(\mathbf{g}_k = \partial E/\partial \mathbf{n}\) at \(\mathbf{n}_k\).
2. **Spectral (BB1) step length.**  With \(\mathbf{s}_k = \mathbf{n}_k - \mathbf{n}_{k-1}\),
   \(\mathbf{y}_k = \mathbf{g}_k - \mathbf{g}_{k-1}\),
   \(\alpha_k^{\mathrm{BB}} = (\mathbf{s}_k\cdot\mathbf{s}_k)/(\mathbf{s}_k\cdot\mathbf{y}_k)\)
   safeguarded into \([10^{-10}, 10^{10}]\).  At the first iteration the fallback
   is \(\alpha_0 = 1/\max(1, \|\mathbf{g}\|_\infty)\).
3. **Spectral projected direction**
   \(\mathbf{d}_k = P_{\mathcal{C}}(\mathbf{n}_k - \alpha_k^{\mathrm{BB}}\, \mathbf{g}_k)
   - \mathbf{n}_k\).  This is provably a descent direction for any closed
   convex \(\mathcal{C}\) (BMR Lemma 2.1), with
   \(\mathbf{g}_k\cdot\mathbf{d}_k \le -\|\mathbf{d}_k\|^2 / \alpha_k^{\mathrm{BB}} \le 0\).
4. **Monotone Armijo** along \(\mathbf{n}_k(\lambda) = \mathbf{n}_k + \lambda \mathbf{d}_k\),
   \(\lambda = 1, \rho, \rho^2, \ldots\): accept the smallest \(\lambda\) with
   \(E(\mathbf{n}_k(\lambda)) \le E(\mathbf{n}_k) + c_1 \lambda\, \mathbf{g}_k\cdot\mathbf{d}_k\).
   Each trial is the convex combination
   \((1-\lambda)\mathbf{n}_k + \lambda P_{\mathcal{C}}(\cdots)\) of two feasible
   points, hence automatically feasible — **no per-trial re-projection**.

`rdmft_occ_optimizer` is **not** consulted on this path: SPG already has a
proven globally convergent step, so no SD/CG direction or curvature
history is used.

**Closed form of \(P_{\mathcal{C}}\).**  The Euclidean foot point on
\(\mathcal{C}\) has the standard shifted-clipping form
\(y_{ik}(\lambda) = \min(1, \max(0, x_{ik} - \lambda\, w_k))\) with a single
scalar \(\lambda\) determined by bisection from
\(\sum_k w_k \sum_i y_{ik}(\lambda) = N_e\).  We do **not** pre-clip \(x\)
to \([0,1]\) before solving for \(\lambda\) — the box clip is part of the
dual map.

**Convergence criterion (inner).**  The textbook SPG stationarity measure is
the Bertsekas projected-gradient residual at unit step,
\(\mathbf{r}(\mathbf{n}) = \mathbf{n} - P_{\mathcal{C}}(\mathbf{n} - \nabla_{\mathbf{n}}E)\)
(BMR Eq. 2.6).  At a KKT point \(\mathbf{r} = 0\).  The inner loop stops when
\(\|\mathbf{r}\|_\infty \le\) `rdmft_occ_proj_tol` (no energy-difference test).
The iteration cap is `rdmft_occ_maxiter`.

**Outer (alternating / joint).**  A run is not considered **globally** converged
on `rdmft_energy_tol` alone: the outer loop requires (after the first cycle)
**both** inner sub-problem convergence flags and, when `rdmft_energy_tol` \(>0\), a
small \(|E - E_{\mathrm{prev}}|\) between outer iterations.  If
`rdmft_energy_tol` \(\le 0\), the outer energy test is off and the loop stops
when the inner flags alone indicate convergence.

### Optimiser selection

| Keyword | Type | Default | Allowed values |
|---------|------|---------|----------------|
| `rdmft_occ_optimizer` | string | `cg` | `sd`, `cg` |
| `rdmft_orb_optimizer` | string | `cg` | `sd`, `cg` |
| `rdmft_orb_retraction` | string | `polar` | `polar`, `qr`, `cayley` (`qr` / `cayley` serial-only; MPI builds fall back to `polar` with a one-time warning) |
| `rdmft_joint_optimizer` | string | `cg` | `sd`, `cg` |

`rdmft_occ_optimizer` and `rdmft_orb_optimizer` control the two sub-problem
optimisers used by the `alternating` strategy.  `rdmft_joint_optimizer`
selects the **single** unified optimiser used by the `joint` strategy on the
packed `(p, C)` variable.  Only the keyword matching the active strategy is
consulted; the others are ignored.

`rdmft_orb_optimizer` chooses the Riemannian optimiser used by the
`alternating` orbital sub-problem (Stiefel SD or CG; see *Optimiser semantics*
below).  Both are first-order Riemannian methods on `St(N_b, N_basis)` with
non-monotone Strong Wolfe line search along the retraction.

`rdmft_orb_retraction` selects the Stiefel retraction used by every orbital
step (both `alternating` and `joint`):

| Value | Algorithm | Notes |
|-------|-----------|-------|
| `polar` (default) | $R_X(\eta) = (X+\eta)\bigl[(X+\eta)^\dagger (X+\eta)\bigr]^{-1/2}$ via Cholesky-QR. | Cheap (one Cholesky + one triangular solve) and fully MPI-parallelised (ScaLAPACK `pdpotrf` + `pdtrsm`). Loses about half the working precision when $(X+\eta)^\dagger (X+\eta)$ is ill-conditioned. |
| `qr` | Householder QR of $X+\eta = QR$ with sign-fixed $R$ diagonal. | More numerically robust than `polar` when $(X+\eta)^\dagger(X+\eta)$ is near-singular. **Serial-only**; MPI builds emit a one-time warning and fall back to `polar`. |
| `cayley` | Wen–Yin low-rank Cayley retraction (Wen & Yin, *Math. Prog.* **142** (2013) 397, Algorithm 1). Solves a $2p \times 2p$ system via Sherman–Morrison–Woodbury; preserves orthogonality exactly without a triangular factor. | Attractive when $N_{\mathrm{basis}} \gg N_{\mathrm{bands}}$. **Serial-only**; MPI builds emit a one-time warning and fall back to `polar`. |

You can choose different optimisers for the occupation and orbital
sub-problems:

| Value | Algorithm | Notes |
|-------|-----------|-------|
| `sd` | Steepest descent | Most robust, slowest convergence |
| `cg` | Conjugate gradient (Polak–Ribière / Fletcher–Reeves with Powell restart) | Good balance of speed and reliability. For **orbitals**, if the line search fails, the next inner iteration restarts along steepest descent (-G_R). |

Both optimisers are available for `rdmft_occ_optimizer` and `rdmft_orb_optimizer`
and may be mixed (e.g. `cg` on occupations and `sd` on orbitals).

### Convergence control

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_outer_maxiter` | int | `200` | Maximum **outer** RDMFT cycles: each cycle is one occupation optimisation plus one orbital optimisation (`alternating`), or one joint product-manifold step (`joint`). |
| `rdmft_occ_maxiter` | int | `50` | Maximum **inner** iterations for the occupation sub-problem (orbitals fixed) within one outer cycle. |
| `rdmft_orb_maxiter` | int | `50` | Maximum **inner** iterations for the orbital sub-problem (occupations fixed) within one outer cycle. |
| `rdmft_energy_tol` | real | `1e-8` | Outer: require **smaller** than this for \(\|E - E_{\mathrm{prev}}\|\) (Ry) **and** both occupation- and orbital-inner `converged` flags, after the first outer step.  Set `<=0` to disable the energy part (stopping uses inner flags only). |
| `rdmft_orb_grad_tol` | real | `1e-5` | Alternating orbitals and **joint** orbital block: relative factor \(\varepsilon_g\) — stop when \(\|G_R\|_F \le \varepsilon_g \max(1, \|G_R(x_0)\|_F)\), with \(x_0\) the iterate at the start of the orbital inner loop (joint: reference norms at outer iteration 0). |
| `rdmft_occ_grad_tol` | real | `1e-5` | **ALM** occupations and **joint** occupation block: \(\varepsilon_g\) for \(\|\nabla_p L\| \le \varepsilon_g \max(1, \|\nabla_p L(x_0)\|)\) with \(x_0\) at the first ALM inner iteration (respectively joint outer iter 0). Not used on the SPG path. |
| `rdmft_occ_proj_tol` | real | `1e-5` | **SPG** (`projected_gradient` / `active_set`): stop when \(\|n - P_\Omega(n - \nabla_n E)\|_\infty \le \varepsilon_{\mathrm{proj}}\) (Bertsekas residual; BMR SIOPT 2000). |

### Initial occupation setup

Before RDMFT optimisation starts, occupations are initialised uniformly as
one of the modes below.

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_occ_init_mode` | string | `ks` | Initial occupation mode: `ks`, `perturbed`, `binary`, `uniform`. |
| `rdmft_occ_init_perturb` | real | `0.0` | Perturbation magnitude `delta` used by `rdmft_occ_init_mode = perturbed`. |
| `rdmft_occ_init_nbands_top` | int | `0` | Fermi-window half-width `K` (per k-point): `K` bands above and `K` bands below the Fermi boundary. Used by `perturbed`, `binary`, and `uniform`. For `perturbed`, `K > 0` sets the window; `K <= 0` uses an automatic small window (default half-width 3, capped by `nbands`). For `binary`/`uniform`, `K <= 0` falls back to `ks`. |

- `ks`: use KS occupations directly as RDMFT initial occupations.
- `perturbed` (ELK/Exciting-style): for each k-point, apply `+delta` to selected empty-side bands and `-delta` to selected occupied-side bands **only within the Fermi window** around the boundary (same band set as `binary`/`uniform` for a given `K`), then project back to the feasible set. This avoids shifting deep valence and high empty bands, which is closer to only perturbing fractional `OCCSV` near the Fermi level in codes such as ELK.
- `binary`: same Fermi window as above; set occupations to `1-delta` if the KS value is `>= 0.5`, otherwise `delta`, then project.
- `uniform`: for each k-point, compute `N_top` over the selected `2K` Fermi-window bands and set each selected occupation to `N_top / (2K)` (or `N_top / N_selected` near band edges), then project.
- If `K <= 0`, `binary` and `uniform` fall back to `ks`. If `delta <= 0`, `perturbed` falls back to `ks`. When `K <= 0` and `perturbed` is active, a default Fermi window is used (see table).

### Line search and optimiser tuning

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_alpha_step` | real | `1.0` | Initial trial step \(\alpha_0\) for **Strong Wolfe** on alternating **orbitals**, **ALM** occupations, and **joint** steps (when ALM BB seed is off, orbitals/joint use this directly). The SPG occupation block (`projected_gradient` / `active_set`) uses Barzilai–Borwein plus **monotone** Armijo and does not consult this keyword. |
| `rdmft_line_search_c1` | real | `1e-4` | Sufficient-decrease \(c_1\in(0,1)\): \(\varphi(\alpha)\le f_{\mathrm{ref}}+c_1\alpha\varphi'(0)\) in Strong Wolfe (joint / ALM / orbitals), and the same constant in **SPG** monotone Armijo (here \(f_{\mathrm{ref}}=E(\mathbf{n}_k)\) at the current iterate). Smaller is stricter. |
| `rdmft_line_search_c2` | real | `0.9` | Strong Wolfe curvature: \(|\varphi'(\alpha)|\le c_2|\varphi'(0)|\). Must satisfy `rdmft_line_search_c1` \(< c_2 < 1\). Not used by SPG. |
| `rdmft_line_search_max_iter` | int | `30` | Max **bracket expansion** steps in Strong Wolfe (joint, ALM, orbitals). |
| `rdmft_line_search_max_zoom` | int | `30` | Max **zoom** (interval refinement) iterations inside Strong Wolfe. |
| `rdmft_line_search_nm_memory` | int | `10` | Non-monotone memory \(M\) for Strong Wolfe: \(f_{\mathrm{ref}}=\max\) of the last \(M\) energies at the line-search iterate. Use `1` for monotone decrease vs the current energy only. Not used by SPG occupations (monotone Armijo vs the current energy only). |
| `rdmft_alm_lambda_init` | real | `0.0` | Initial ALM Lagrange multiplier `lambda` (only for `rdmft_constraint = augmented_lagrangian`). |
| `rdmft_alm_mu_init` | real | `1.0` | Initial ALM penalty parameter `mu` (only for `rdmft_constraint = augmented_lagrangian`). |
| `rdmft_alm_mu_factor` | real | `2.0` | Multiplicative ALM penalty update factor: `mu <- min(mu * factor, mu_max)`. |

`rdmft_alm_bb_enabled`, `rdmft_alm_bb_mode`, `rdmft_alm_bb_alpha_min`, and `rdmft_alm_bb_alpha_max` seed the **Strong Wolfe** initial step for **ALM** occupations in occupation-parameter space. **ALM**, alternating **orbitals**, and **joint** use **non-monotone Strong Wolfe** with `rdmft_line_search_c1`, `rdmft_line_search_c2`, and the bracket/zoom iteration caps. The **SPG** path (`projected_gradient` / `active_set`) keeps the standard **BB1 spectral step + monotone Armijo** along the projected direction; it does **not** use Strong Wolfe or `rdmft_occ_optimizer`.

The **alternating orbital** sub-problem is Riemannian SD or Polak–Ribière⁺ CG on the Stiefel manifold (Absil–Mahony–Sepulchre, *Optimization Algorithms on Matrix Manifolds*, Princeton 2008), with **non-monotone Strong Wolfe** along `rdmft_orb_retraction` (default `polar`). For ill-conditioned regularised functionals (Müller / Power / GEO), reduce `rdmft_alpha_step` or switch to `rdmft_solver_strategy = joint` if the alternating map diverges (see `rdmft_derivation.md` §7.1).

### Debugging

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_grad_check` | bool | `false` | If `true`, run a finite-difference gradient verification **at the start of the RDMFT solve**, immediately after the LCAO coefficients are converted to internal **X-space** (`precompute_cholesky_S` / `wfc_C_to_X`) and before the main alternating or joint loop. Occupations use a central difference in `n`; the **orbital** check uses the projected Riemannian gradient `G_R` and a **forward** difference along the same polar retraction as the line search `(E(t)-E_0)/t`, which matches the directional derivative \(\langle G_R, D\rangle\) at \(t=0\) more reliably than a symmetric `±t` probe. This is expensive but essential for development and validation. Results are written to the running log. |

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

### Example 2: Power functional with CG and gradient check

```
INPUT_PARAMETERS
calculation         scf
basis_type          lcao
dft_functional      pbe
rdmft_functional     power
rdmft               1
rdmft_power_alpha       0.55
rdmft_solver_strategy   alternating
rdmft_occ_optimizer     cg
rdmft_orb_optimizer     cg
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
rdmft_joint_optimizer   cg
rdmft_alpha_step         0.01
rdmft_outer_maxiter        300
```

The explicit `rdmft_alpha_step 0.01` is smaller than the default `1.0` and
can stabilise joint HF steps; omit it to use the default.

Here occupations and orbitals are optimised simultaneously on the product
manifold (new `joint` keyword; `product_manifold` is still accepted as a
deprecated alias).  The packed variable `(p, C)` is handled by a single
optimiser (`rdmft_joint_optimizer` = `sd` or `cg`): every outer iteration the
solver evaluates the packed gradient `(dE/dp, G_R)`, computes a descent
direction, re-projects the orbital block onto the Stiefel tangent space, and
runs one Armijo line search along the packed direction (linear update for the
occupation parameters, Stiefel retraction for the orbitals).

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

- **Projected gradient (SPG)** ignores `rdmft_occ_param` and optimises raw
  \(n_{ik}\) with a Barzilai–Borwein spectral step and a **monotone** Armijo line
  search along the spectral projected direction; no tuning is normally required.

- **CG** is the recommended optimiser for production in most cases.  Steepest
  descent is useful for debugging.

- **For periodic systems with multiple k-points**, the constraint
  Σ_k w_k Σ_i n_{ik} = N_e is automatically handled using the k-point weights.

---

## Troubleshooting

### SPG occupation line search fails

The SPG inner loop logs lines of the form
`occ line search (SPG monotone Armijo): ... ok|fail`.  Because the
spectral projected direction is *guaranteed* to be a descent direction
(BMR Lemma 2.1), the line search can only fail in floating-point
arithmetic when `g · d ≈ 0` (the iterate is essentially stationary).  When
this happens a zero step is accepted and the outer loop's convergence
criteria still apply.  If failures occur at non-stationary points, run
`rdmft_grad_check 1` to verify \(\partial E/\partial n\) against finite
differences.

### ABACUS runs for a long time with no RDMFT output

Typical causes (in roughly decreasing likelihood):

1. **`rdmft_functional` is not set.** Only the legacy single-step evaluator
   runs and no `===== RDMFT Alternating Optimization =====` banner is printed.
   Set `rdmft_functional muller` (or `hf`/`power`/`gu`) to activate the
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
