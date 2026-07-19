# RDMFT example inputs

Two minimal `pw.x` inputs that exercise the alternating RDMFT solver.
Each ships an ABACUS-style ``&rdmft`` namelist appended at the end of
the standard PWscf input file.  Both rely on a hydrogen pseudopotential
shipped with QE itself; copy `H_HSCV_PBE-1.0.UPF` from
`CPV/examples/EXX-wf-example/` into the directory you run from (or
adjust `pseudo_dir`).

## Molecule case: H2 (`h2_molecule/h2.in`)

```
mkdir tmp && pw.x < h2.in > h2.out
```

* PWscf first runs a Hartree-Fock SCF (`input_dft = 'hf'`,
  `ace = .true.`) so the EXX-via-ACE machinery is initialised.
* The post-SCF RDMFT solver is then invoked from `electrons.f90` (via
  the `do_rdmft = .true.` switch in the `&rdmft` block).
* The chosen `rdmft_functional = 'muller'` reuses the already-built
  ACE projectors with band weights `wk * sqrt(n_ik)` to compute the
  RDMFT exchange Hamiltonian.

The expected RDMFT energy is slightly below the converged HF energy
(Müller captures some left-right correlation), e.g.

```
!RDMFT_ETOTAL  =      -2.2422885283 Ry
```

vs the KS-HF energy reported just above the RDMFT block:

```
!    total energy              =      -2.24112743 Ry
```

## Solid case: 4-atom H chain (`h4_chain_solid/h4_chain.in`)

A periodic 4-H chain in a tetragonal box with vacuum padding (Γ-point
sampling).  Same KS-HF + RDMFT-Müller workflow as the molecule case.

```
mkdir tmp && pw.x < h4_chain.in > h4_chain.out
```

Note that for very stretched / strongly-correlated occupation patterns,
projected-gradient + Müller can oscillate when many bands sit at
`n -> 0` (stiff `g'(n)`).  The recommended robust setup mirrors ELK /
ABACUS practice: a **Fermi-window** occupation seed
(`rdmft_occ_init_mode = 'perturbed'`, `rdmft_occ_init_perturb = 1e-3`,
`rdmft_occ_init_nbands_top = 1` for a minimal window, or larger `K`
when needed) so deep occupied bands are not artificially shifted.
Set `rdmft_grad_check = .true.` in `&RDMFT` to verify analytic gradients before
the optimisation loop (development only).  For HF with occupation entropy, use
`rdmft_temp` in Kelvin for HF occupation entropy (HF functional only; uses `k_B*T`).

The default `rdmft_reg_eps = 1e-8` uses ABACUS piecewise power regularisation
(`g=n^alpha` for `n >= eps`, linear Taylor below `eps`); raise
`rdmft_reg_eps` (e.g. `1e-3`) if the occupation block remains stiff.  See
`PW/src/rdmft/README.md` for the full keyword reference.

## Multi-k and tot_charge verification

`multik_charged/` confirms that the alternating RDMFT-HF solver
reproduces KS-HF to nano-Ry on a multi-k bulk Si calculation
(2x2x2 mesh, 3 irreducible k-points) and on charged molecular
systems (`tot_charge = +2` closed-shell H<sub>4</sub><sup>2+</sup>
and `tot_charge = -1` open-shell H<sub>2</sub><sup>-</sup>).  See
`multik_charged/README.md` for the result table and the multi-k /
charge-handling notes (in particular the multi-k g2_kin /
buffer-save fixes added to `PW/src/rdmft/`).

## HF efficiency benchmark

`hf_benchmark/` holds a side-by-side timing of the three optimisers
(`sd` / `cg` / `lbfgs`) for both `'alternating'` and `'joint'`
strategies on H2 and bulk Si with `rdmft_functional = 'hf'`.  See
`hf_benchmark/README.md` for the result tables and conclusions.

## Trying other XC functionals

Change `rdmft_functional` to one of:

| value     | g(n) / kernel form                                          |
|-----------|--------------------------------------------------------------|
| `hf`      | identity (g(n) = n) — exact at the SCF orbital              |
| `muller`  | `sqrt(n)` (default in the example)                          |
| `power`   | `n^alpha`, set `rdmft_power_alpha = 0.656` for solids       |
| `gu`      | Goedecker-Umrigar (sqrt branch in the production code)      |
| `chf`     | Corrected Hartree-Fock (Csányi-Arias)                       |
| `cga`     | Csányi-Goedecker-Arias                                      |
| `geo`     | Three-channel mixture                                        |
| `hybopt`  | HF + Power(α=0.541) convex mixture                          |
| `bow`     | Baldsiefen BOW kernel (exact pair path; expensive)          |
| `bowmod`  | Separable BOW approximation (channel path)                  |

For mixed kernels (`gu`, `chf`, `cga`, `geo`, `hybopt`, `bowmod`) the energy
evaluator runs one ACE build per channel and accumulates the
contributions with the fixed coefficients listed in
`rdmft_xc_channel`.
