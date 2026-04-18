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
| `rdmft_solver_strategy` | string | `alternating` | `alternating`, `product_manifold` |

- **`alternating`** — Alternate between optimising occupation numbers (orbitals
  fixed) and orbitals (occupations fixed).  This is the standard approach and
  the most robust.
- **`product_manifold`** — Treat occupations (Euclidean space, after
  parameterisation) and orbitals (Stiefel manifold) as a single product
  manifold and optimise everything simultaneously.

### Occupation parameterisation

| Keyword | Type | Default | Allowed values |
|---------|------|---------|----------------|
| `rdmft_occ_param` | string | `cosine_sq` | `cosine_sq`, `logistic` |

Both options map an unconstrained real parameter to the interval [0, 1]:

- **`cosine_sq`** — n = cos²(θ).  The gradient transforms as
  dE/dθ = −sin(2θ) · dE/dn.
- **`logistic`** — n = σ(x) = 1/(1+e^{−x}).  The gradient transforms as
  dE/dx = n(1−n) · dE/dn.

### Electron-number constraint

| Keyword | Type | Default | Allowed values |
|---------|------|---------|----------------|
| `rdmft_constraint` | string | `augmented_lagrangian` | `augmented_lagrangian`, `projected_gradient`, `active_set` |

The total electron number must satisfy Σ_k w_k Σ_i n_{ik} = N_e.

- **`augmented_lagrangian`** — Adds a penalty λ·c + (μ/2)·c² to the energy.
  The multiplier λ and penalty μ are updated automatically.  Works well with
  the cosine_sq or logistic parameterisation where the box constraint [0,1]
  is already built in.
- **`projected_gradient`** — After each gradient step, clip occupations to
  [0,1] and rescale to satisfy the electron-number constraint.  Simple but can
  be slow to converge.
- **`active_set`** — Track which occupations are pinned at 0 or 1, solve the
  reduced equality-constrained problem on the free variables.

### Optimiser selection

| Keyword | Type | Default | Allowed values |
|---------|------|---------|----------------|
| `rdmft_occ_optimizer` | string | `cg` | `sd`, `cg`, `lbfgs`, `adam` |
| `rdmft_orb_optimizer` | string | `cg` | `sd`, `cg`, `lbfgs`, `adam` |

You can choose different optimisers for the occupation and orbital
sub-problems:

| Value | Algorithm | Notes |
|-------|-----------|-------|
| `sd` | Steepest descent | Most robust, slowest convergence |
| `cg` | Conjugate gradient (Polak–Ribière / Fletcher–Reeves with Powell restart) | Good balance of speed and reliability |
| `lbfgs` | Limited-memory BFGS | Fast for smooth landscapes, uses `rdmft_lbfgs_memory` history vectors. For the orbital sub-problem, the quasi-Newton direction is projected back onto the Stiefel tangent space (Riemannian L-BFGS by projection). |
| `adam` | Adam | Adaptive learning rate, useful for noisy or ill-conditioned problems. For the orbital sub-problem Adam's Euclidean update is projected onto the tangent space and retracted onto the Stiefel manifold at each step. |

All four optimisers are available for both `rdmft_occ_optimizer` and
`rdmft_orb_optimizer` and may be mixed freely (e.g. `cg` on occupations and
`lbfgs` on orbitals).

### Convergence control

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_outer_maxiter` | int | `200` | Maximum **outer** RDMFT cycles: each cycle is one occupation optimisation plus one orbital optimisation (`alternating`), or one joint product-manifold step (`product_manifold`). |
| `rdmft_occ_maxiter` | int | `50` | Maximum **inner** iterations for the occupation sub-problem (orbitals fixed) within one outer cycle. |
| `rdmft_orb_maxiter` | int | `50` | Maximum **inner** iterations for the orbital sub-problem (occupations fixed) within one outer cycle. |
| `rdmft_energy_tol` | real | `1e-8` | Convergence threshold on the change in total energy (Ry) between outer steps. |
| `rdmft_grad_tol` | real | `1e-6` | Convergence threshold on the norm of the gradient. |

### Initial KS occupation adjustment

Before RDMFT optimisation starts, occupations may be pushed slightly away from 0 and 1 (margin `m` in `RDMFTConfig`, default `1e-3`) so the cosine-squared / logistic parameterisation has a non-vanishing Jacobian at the KS seed. By default only the **highest** bands are perturbed; lower bands keep the KS values.

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_occ_init_nbands_top` | int | `5` | Number of **highest** bands (per k-point) that receive the `[m, 1−m]` clamp and participate in the electron-count rescaling. Bands with index `ib < nbands − K` are left unchanged. Set to `0` to apply the adjustment to **all** bands (legacy behaviour). |

### Line search and optimiser tuning

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_alpha_step` | real | `0.1` | Initial trial step length for the Armijo backtracking line search. |
| `rdmft_lbfgs_memory` | int | `10` | Number of past gradient/step pairs stored by L-BFGS. |
| `rdmft_adam_lr` | real | `0.001` | Learning rate for the Adam optimiser. |

### Debugging

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `rdmft_grad_check` | bool | `false` | If `true`, run a finite-difference gradient verification before the optimisation begins. This is expensive (one extra energy evaluation per checked variable) but essential for development and validation. Results are written to the running log. |

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

### Example 2: Power functional with L-BFGS and gradient check

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

### Example 3: Product-manifold with logistic parameterisation

```
INPUT_PARAMETERS
calculation         scf
basis_type          lcao
dft_functional      pbe
rdmft_functional     hf
rdmft               1
rdmft_solver_strategy   product_manifold
rdmft_occ_param         logistic
rdmft_constraint        augmented_lagrangian
rdmft_occ_optimizer     adam
rdmft_orb_optimizer     sd
rdmft_adam_lr            0.005
rdmft_alpha_step         0.01
rdmft_outer_maxiter        300
```

Here occupations and orbitals are optimised simultaneously on the product
manifold.  Occupations use the logistic parameterisation with Adam, while
orbitals use steepest descent.

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
automatically detected.  The projected-gradient method clips-and-rescales
occupations after every step.

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

- **Three iteration limits.** `rdmft_outer_maxiter` limits the **outer** RDMFT loop (each cycle: occupation step + orbital step in `alternating`, or one joint step in `product_manifold`). `rdmft_occ_maxiter` and `rdmft_orb_maxiter` limit the **inner** optimisations of occupations (fixed orbitals) and orbitals (fixed occupations) within each outer cycle. In older versions, `rdmft_orb_maxiter` incorrectly doubled as the outer-loop limit; use `rdmft_outer_maxiter` for that now.

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

- **CG or L-BFGS** are the recommended optimisers for production.  Steepest
  descent is useful for debugging.  Adam can help with difficult convergence
  but may require tuning `rdmft_adam_lr`.

- **For periodic systems with multiple k-points**, the constraint
  Σ_k w_k Σ_i n_{ik} = N_e is automatically handled using the k-point weights.

---

## Troubleshooting

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
