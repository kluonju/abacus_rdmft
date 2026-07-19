# RDMFT on bcc Fe — three-step DOS workflow (`pw.x` → `pw.x` fine → `rdmft_dos.x`)

End-to-end demonstration of QE's standard SCF / DOS split, applied
to RDMFT on ferromagnetic bcc Fe (1-atom unit cell, 2×2×2 coarse k
mesh refined to 4×4×4, nspin = 2, PZ-LDA SCF + RDMFT-HF spectral DOS).
This is the multi-step companion to the existing single-binary
`PW/src/rdmft/examples/Fe_bcc_dos/` example.

| Step | Executable | Role |
|---|---|---|
| 1 | `pw.x`        | PZ-LDA SCF + RDMFT setup on 2×2×2; saves natural occupations + orbitals.  No DOS computed here. |
| 2 | `pw.x`        | Fine 4×4×4 mesh: k-interpolate `(n, C^k)` from coarse save, optional orbital polish, RDMFT save. |
| 3 | `rdmft_dos.x` | dedicated DOS post-processor on the **fine** save: TSM probe → `fe_fine.rdmft.dos` / PDOS / Mz. |

The point of this example is to exercise `rdmft_dos.x` on a non-trivial
spin-polarised, multi-k, transition-metal d-band system -- the
characteristic majority / minority Fe d-band splitting must reproduce
in the saved DOS exactly the same way the original `Fe_bcc_dos/`
example does inline.

## Files in this directory

| File | Purpose |
|------|---------|
| `fe.step1.in` | `pw.x` input: PZ SCF + RDMFT post-loop save on 2×2×2.  No DOS at this stage. |
| `fe.step2.in` | `rdmft_dos.x` input: DOS on the **coarse** save (same k mesh as step 1). |
| `fe.step2.fine.in` | `pw.x` input: 4×4×4 mesh + auto k-refine from step-1 coarse save. |
| `fe.step3.in` | `rdmft_dos.x` input: DOS on the **fine** save from step 2. |
| `Fe.pz-n-nc.UPF` | PZ norm-conserving Fe pseudopotential (with `<PP_CHI>` blocks needed for PDOS / Mz). |
| `README.md` | This file. |

## How to run

```bash
# 1) Coarse RDMFT calculation (pw.x, 2x2x2)
mkdir -p tmp
pw.x < fe.step1.in > fe.step1.out

# 2a) Coarse DOS (optional sanity check on the 2x2x2 save)
rdmft_dos.x < fe.step2.in > fe.step2.out

# 2b) Fine-grid pw.x restart with k-refinement from step 1
mkdir -p tmp_fine
pw.x < fe.step2.fine.in > fe.step2.fine.out

# 3) Fine-grid DOS on the converged fine save
rdmft_dos.x < fe.step3.in > fe.step3.out
```

Step 1 finishes in a couple of seconds on a 4-core machine (10 bands,
6 k-points, no RDMFT iterations).  Step 2b adds k-interpolation and
optional orbital polish; step 3 runs the TSM probe on the finer
irreducible k-set.

## What each input does

**`fe.step1.in` — RDMFT calculation.**  Plain PZ-LDA SCF
(`input_dft` is **not** set so QE picks up the PZ XC from the
pseudopotential), with `do_rdmft = .true.` and
`rdmft_outer_maxiter = 0` so the optimiser is skipped.  At the end of
`rdmft_run` the save routines write:

* `tmp/fe.rdmft.save` -- natural occupations seeded from the
  KS-smeared SCF weights;
* `tmp/fe.save/wfc*.dat` + `data-file-schema.xml` -- the KS
  wavefunctions (which double as natural orbitals when
  `rdmft_outer_maxiter = 0`).

`rdmft_dos.x` is run separately in step 2, so no DOS files are produced
during step 1 and the expensive transition-state band probe is **not** run.

**`fe.step2.in` — DOS calculation.**  Single `&inputrdmftdos` namelist
consumed by `rdmft_dos.x`:

```
&inputrdmftdos
  prefix             = 'fe'
  outdir             = './tmp/'
  functional         = 'hf'        ! must match the previous pw.x run
  emin               = -0.882
  emax               =   0.588
  deltae             = 0.00368
  degauss            = 0.00735
  pdos               = .true.
/
```

The post-processor internally:

1. Calls `read_file_new(needwf=.TRUE.)` to load PWscf state.
2. Refills the wfc buffer from `tmp/fe.save/wfc*.dat`.
3. Calls `sym_rho_init` so `sum_band → sym_rho` symmetrises the
   density correctly on the bcc Fe point group.
4. Forces XC to `(5,0,0,0,0,0)` (HF), `exxalfa = 1.0`,
   `use_ace = .TRUE.`, sets up the EXX grid via `setup_exx`, runs
   `stop_exx` so the next `exxinit` call computes the Coulomb
   divergence treatment from scratch, then `rdmft_exxinit_once`.
5. Loads `rdmft_n` from `tmp/fe.rdmft.save`.
6. Computes the Ewald energy.
7. Calls `rdmft_total_energy` to refresh the density / vrs against
   the loaded `(n, evc)` state.
8. Runs `rdmft_compute_dos()` and `rdmft_compute_magnetization()` --
   the TSM probe is `O(N_b * N_k)` single-k ACE builds, executed
   only here, never during step 1.

**`fe.step2.fine.in` — fine-grid pw.x restart.**  Uses a **new**
`prefix = 'fe_fine'` and `outdir = './tmp_fine/'` so the coarse
step-1 save in `./tmp/fe` is never overwritten.  Finer
`K_POINTS automatic 4 4 4` and

```
&rdmft
  ...
  rdmft_source_prefix   = 'fe'        ! coarse save (read-only)
  rdmft_source_outdir   = './tmp/'    ! coarse save (read-only)
/
```

Because the fine `K_POINTS` mesh is denser than the coarse step-1 mesh,
k-interpolation is triggered automatically (no
extra flag needed).  KS SCF on the fine mesh is
skipped (`electron_maxstep = 0`); `pw.x` trilinearly interpolates
`(n, C^k)` from the coarse save, runs a short fixed-occupation orbital
polish (default 5 iterations), and writes `tmp_fine/fe_fine.rdmft.save`
plus `tmp_fine/fe_fine.save/`.
The fine Monkhorst–Pack grid comes **only** from `K_POINTS`; coarse
mesh metadata is read from the source save.

**`fe.step3.in` — fine-grid DOS.**  Same namelist as `fe.step2.in` but
points at `prefix = 'fe_fine'`, `outdir = './tmp_fine/'`.  No k-refine
keys — `rdmft_dos.x` loads the converged fine save directly.

### Interpreting energies during k-refinement

When `rdmft_krefine_orb_maxiter > 0`, step 2b prints an orbital-polish
block after k-interpolation and EXX init.  Important points:

* **`RDMFT orb inner ... global E =`** reports the **global** RDMFT
  total energy (summed over all k-points, including Ewald) after each
  per-k line search — **not** a per-k contribution and not a TSM band
  energy.
* Compare the **final** `polish iter N: global E =` line or the
  **`Initial RDMFT total energy`** line in step 2b to the coarse step-1
  `!RDMFT_ETOTAL` — not the mid-sweep inner lines from the first
  k-sweep, which can look very different while `||G_R^k||` is still
  large.
* Interpolated fine-grid orbitals start far from the variational
  minimum; increase `rdmft_krefine_orb_maxiter` if gradient norms
  remain O(0.1) or larger after the last polish iteration.

## Expected output

```bash
grep -E '!RDMFT_ETOTAL|RDMFT total energy at loaded' fe.step1.out fe.step2.out
# fe.step1.out:     !RDMFT_ETOTAL  =     -33.6692554395 Ry
# fe.step2.out:     RDMFT total energy at loaded (n, evc) =     -33.6692554395 Ry

grep -A 3 'Per-atom magnet' fe.step1.out fe.step2.out
# both files report:
#       atom  species      Mz        Q_up       Q_dn
#         1  Fe      1.970213    4.977756    3.007543

ls fe.rdmft.*
# fe.rdmft.dos   fe.rdmft.mag   fe.rdmft.pdos_at001_001
```

The per-atom Mulliken Mz reproduces the inline `Fe_bcc_dos/`
reference value of 1.97 µB exactly, and matches the SCF
`total_magnetization` from `fe.step1.out` to ~ 1 %.  The total
energy and all energy components (including E_xc) are bit-identical
between the two runs after the `exx_bgrp_type` fix described below.

## Why this matters

* This is the dedicated `pw.x` / `dos.x`-style two-step workflow:
  the heavy SCF + RDMFT setup runs once, then `rdmft_dos.x` can be
  re-run as many times as needed to scan
  `degauss` / `emin/max` / `deltae` /
  `occ_weighted` / `pdos` without redoing the
  optimiser.
* The transition-state-model band probe (the dominant cost of an
  RDMFT-DOS calculation on a transition-metal cell) is computed
  **only in step 2**, on demand, never during step 1.
* The wavefunction buffer (`tmp/fe.save/wfc*.dat`) used by step 2 is
  the same file QE writes for any SCF run, so the saved orbitals can
  also be consumed by `dos.x`, `projwfc.x`, `bands.x`, `pp.x` etc.
  for the underlying KS spectrum / projections.
