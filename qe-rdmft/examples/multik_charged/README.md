# Multi-k and tot_charge verification

Verifies that the alternating RDMFT-HF solver gives the correct
Kohn-Sham Hartree-Fock total energy in two regimes that are absent
from the small `examples/h2_molecule` and `examples/h4_chain_solid`
demos:

1. **Multi-k Si** — bulk silicon with a 2x2x2 Monkhorst-Pack mesh
   (3 irreducible k-points), reduced from the gamma-only
   `examples/hf_benchmark` setup.
2. **Charged systems** (`tot_charge != 0`):
   - Closed-shell H<sub>4</sub><sup>2+</sup> (4 H atoms, +2 charge,
     2 electrons, `nspin = 1`).
   - Open-shell H<sub>2</sub><sup>-</sup> (2 H atoms, -1 charge,
     3 electrons, `nspin = 2`, `tot_magnetization = 1`).

## How to run

```
PWX=path/to/pw.x bash run.sh
```

The script builds the `pw.x` inputs, runs each case, and prints a
single comparison table.

## Reference results

```
case                                                       KS-HF         RDMFT-HF         Ne          dE_KS
---------------------------------------- ---------------- ---------------- ---------- --------------
Multi-k Si (2x2x2), alt, occ-only                    -15.13128573   -15.1312857322   8.000000      -2.20e-09
H4^2+ (tot_charge=+2), nspin=1, alt                   -1.80795595    -1.8079559531   2.000000      -3.10e-09
H2^- (tot_charge=-1), nspin=2, alt                    -2.23755180    -2.2375518001   3.000000      -1.00e-10
Multi-k Si (2x2x2), alt, orb_strategy=joint          -15.13128573   -15.1312857322   8.000000      -2.20e-09
Multi-k Si (2x2x2), alt, orb_strategy=block_k        -15.13128573   -15.1312857322   8.000000      -2.20e-09
```

`dE_KS` is the difference between the converged RDMFT-HF total energy
and the converged KS-HF total energy reported by the same `pw.x`
binary on the same input.  All three cases agree to ≤ 3 nRy.

## Settings used

```
rdmft_functional       = 'hf'
rdmft_solver_strategy  = 'alternating'
rdmft_constraint       = 'projected_gradient'
rdmft_outer_maxiter    = 4
rdmft_occ_maxiter      = 10
rdmft_orb_maxiter      = 0          ! occupation-only (see Notes below)
rdmft_energy_tol       = 1.0e-8
rdmft_occ_grad_tol     = 1.0e-6
rdmft_occ_tol          = 1.0e-6     ! freeze-occ if sum|dn|_outer below this
                                     ! for two consecutive cycles
```

## What the fixes do

Multi-k support required two PWscf-specific patches that are now part
of `PW/src/rdmft/`:

* `g2_kin(ik)` is called before every `h_psi` invocation that targets
  a k-point different from the current one (the kinetic operator
  `(k+G)^2` lives in module `wvfct` and is updated only by
  `g2_kin`).  Without this the per-band one-body diagonal at every
  k-point except the SCF's last one was wrong by O(eV/electron).
* The orbital-block Armijo line search saves the trial `evc` to the
  wfc buffer **before** `rdmft_total_energy` rebuilds the density;
  otherwise `sum_band` reads the OLD `evc` for the current k from the
  buffer and the trial energy is unchanged regardless of how big the
  step is.

For `tot_charge != 0` no special code path was needed — the
electron-number target follows QE's `nelec` directly, so positive and
negative charge states map to `Ne = nelec` automatically.

## Notes

The alternating orbital block on multi-k runs supports two
strategies, both of which reproduce KS-HF on Si 2×2×2 to ≤ 3 nRy
(see the table above):

* `rdmft_orb_strategy = 'joint'` (default) — Riemannian descent on
  the **product** Stiefel manifold with one GLOBAL Armijo line
  search per inner iteration.  A single scalar `alpha` retracts
  every `C^k` simultaneously and is accepted only if the total
  energy (summed over all k) decreases sufficiently.

* `rdmft_orb_strategy = 'block_k'` — block-coordinate
  Gauss-Seidel sweep over k-points: at each k take ONE Riemannian
  step, holding all other k's evc fixed (in their wfc buffers).
  Each k owns its own optimiser state (per-k SD / CG / L-BFGS
  history persisted across outer iterations).  The line search
  still uses the global energy through `rdmft_total_energy`, but
  only the current k's evc is varied per trial.

The `block_k` strategy maps naturally onto QE's k-point pool
parallelisation: each MPI pool owns a subset of k-points and can
optimise them independently.  Only the global-energy reduction at
each line-search trial requires a cross-pool `mp_sum`, which is
already performed inside `rdmft_total_energy` through `sum_band` /
`v_of_rho`.

See `PW/src/rdmft/rdmft_solver.f90` :: `rdmft_orbital_step_joint`
and `rdmft_orbital_step_block_k` for the implementations and
`PW/src/rdmft/README.md` for a textual description of both.

* H<sub>2</sub><sup>+</sup> with `nspin=2 tot_magnetization=1` is
  **not** a working test case because PWscf disables ACE when one
  spin channel has zero electrons (the spin-down channel for
  H<sub>2</sub><sup>+</sup> has no electrons).  Use a closed-shell or
  symmetric magnetic charged system instead, as in this folder.
