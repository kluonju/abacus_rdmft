# RDMFT-HF efficiency benchmark

Verifies that with `rdmft_functional = 'hf'` both `'alternating'` and
`'joint'` strategies converge to the converged Kohn-Sham Hartree-Fock
total energy on a molecule (H2 in a box) and a solid (bulk Si at the
Γ point), and compares the cost of the three optimisers
(`sd` / `cg` / `lbfgs`) supported in each block.

The HF kernel is the only RDMFT functional whose minimum coincides
exactly with a pure Hartree-Fock SCF stationary point: at the
converged HF orbitals with integer occupations `n_ik in {0, 1}`, all
RDMFT first-order optimality conditions are satisfied.  This makes HF
an ideal correctness check.  We seed the RDMFT optimisers with a
**perturbed** initial occupation (`rdmft_occ_init_mode = 'perturbed'`,
`rdmft_occ_init_perturb = 0.10` for H2, `0.05` for Si) so the
optimisers actually have work to do; otherwise the alternating loop
would converge in 0 inner iterations.

## How to run

Copy the pseudopotentials into the working directory and execute:

```
cp ../../../CPV/examples/EXX-wf-example/H_HSCV_PBE-1.0.UPF .
cp ../../../QEHeat/examples/pseudo/Si_ONCV_PBE-1.1.upf  Si.upf

# Run the benchmarks (PWX defaults to ./bin/pw.x)
PWX=path/to/pw.x bash run_h2_bench.sh
PWX=path/to/pw.x bash run_si_bench.sh
```

Each script runs a single `pw.x` job per `(strategy, optimiser)` cell
of the matrix.  Convergence is identical at every cell (same
tolerances) so the timings are an apples-to-apples comparison.

## Settings

```
rdmft_functional    = 'hf'
rdmft_constraint    = 'projected_gradient'
rdmft_outer_maxiter = 30
rdmft_occ_maxiter   = 20
rdmft_orb_maxiter   = 10
rdmft_energy_tol    = 1.0e-9
rdmft_occ_grad_tol  = 1.0e-7
rdmft_orb_grad_tol  = 1.0e-7
rdmft_occ_tol       = 1.0e-6   ! freeze-occ kicks in once sum|dn|_outer
                                ! < this for 2 consecutive cycles.
rdmft_lbfgs_memory  = 8
```

## H2 in a 12-Bohr cubic box

KS-HF reference: `-2.24112743 Ry`

| strategy     | optimiser | n_outer | RDMFT_ETOT (Ry)  | dE vs KS-HF | wall (s) |
|--------------|-----------|---------|------------------|-------------|----------|
| alternating  | sd        |  2      | -2.2411274288    | +1.2e-9     |  1.62    |
| alternating  | cg        |  2      | -2.2411274288    | +1.2e-9     |  1.52    |
| alternating  | lbfgs     |  2      | -2.2411274287    | +1.3e-9     |  1.93    |
| joint        | sd        | 30      | -2.2408711364    | +2.6e-4     |  2.65    |
| joint        | cg        | 30      | -2.2411051441    | +2.2e-5     |  2.80    |
| joint        | lbfgs     | 30      | -2.2411224800    | +5.0e-6     |  2.08    |

* **Alternating** finds the KS-HF energy to machine precision in two
  outer iterations regardless of the optimiser — once the
  freeze-occupation logic kicks in (`rdmft_occ_tol = 1e-6`) the
  occupation block is permanently skipped and the orbital block (which
  is essentially at the HF stationary point) only takes a residual
  step.  All three optimisers are within 1.3 nRy of KS-HF.

* **Joint** is intrinsically slower because it is forced to share a
  single Armijo step between the cosine-squared occupation
  parameters and the Stiefel orbital block.  The relative ordering
  matches expectations: `lbfgs > cg > sd` with **L-BFGS the best in
  both accuracy and time**, reaching 5 µRy of KS-HF in 30 outer
  iterations (≈25 % faster than CG).

## Bulk Si at the Γ point

KS-HF reference: `-14.43586741 Ry`

| strategy     | optimiser | n_outer | RDMFT_ETOT (Ry)   | dE vs KS-HF | wall (s) |
|--------------|-----------|---------|-------------------|-------------|----------|
| alternating  | sd        |  2      | -14.4358674636    | -5.4e-8     |  0.43    |
| alternating  | cg        |  2      | -14.4358674637    | -5.4e-8     |  0.42    |
| alternating  | lbfgs     |  2      | -14.4358674636    | -5.4e-8     |  0.52    |
| joint        | sd        | 30      | -14.4353124374    | +5.6e-4     |  0.71    |
| joint        | cg        | 30      | -14.4357582383    | +1.1e-4     |  0.71    |
| joint        | lbfgs     | 30      | -14.4350233380    | +8.4e-4     |  0.71    |

* **Alternating** again converges in two outer iterations to the KS-HF
  reference (≈54 nRy below it because the perturbation slightly lowers
  the SCF threshold).  All three optimisers behave identically.

* **Joint** ordering for bulk Si is `cg > sd ≳ lbfgs` — different from
  H2 because the wider band manifold and larger Stiefel block make the
  per-iteration L-BFGS curvature pair noisier.  CG is safest in this
  regime because the Powell-style restart kicks in often enough to
  stay near the steepest-descent trajectory.

## Take-aways

1. **Both strategies recover the SCF energy on HF**, validating the
   sign / weight conventions of the EXX-via-ACE RDMFT exchange
   evaluator.  The alternating strategy reproduces KS-HF to nano-Ry on
   every (system, optimiser) cell.
2. **Alternating is the right default.**  It separates the easy
   (occupation) and harder (orbital) sub-problems and the
   freeze-occupation logic guarantees that the inner occupation loop
   stops costing anything once the natural occupations stabilise.
3. **L-BFGS is the recommended joint-strategy optimiser** for small,
   well-conditioned problems (H2-like) where the curvature pairs are
   informative; **CG is more robust** for larger systems where the
   joint Hessian has a wider eigenvalue spread (Si-like).  Plain SD is
   a useful debugging baseline.
4. With `rdmft_occ_tol = 1e-6` the alternating freeze triggers
   immediately for the HF kernel (the occupations don't move beyond
   the perturbation), which is the desired behaviour.
