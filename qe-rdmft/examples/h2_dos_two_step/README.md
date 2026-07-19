# RDMFT on H2 — two-step DOS workflow (`pw.x` → `rdmft_dos.x`)

End-to-end demonstration of QE's standard SCF / DOS split, applied
to RDMFT on H2 in a 12-Bohr cubic cell.  This is the same physics as
the `h2_molecule/` example; only the DOS plumbing is different.

| Step | Executable | Role |
|---|---|---|
| 1   | `pw.x`        | RDMFT calculation: runs the alternating Müller optimisation, saves natural occupations + natural orbitals.  No DOS computed here. |
| 2   | `rdmft_dos.x` | dedicated DOS post-processor: reads the saved `(occupation, orbital)` pair, computes the transition-state DOS, writes `*.rdmft.dos` / `*.rdmft.pdos_at*` / `*.rdmft.mag`. |

The two executables are completely separate, exactly the same way
`pw.x` and `dos.x` are separate in the standard QE workflow:
`rdmft_dos.x` is a new program living in `PP/src/rdmft_dos.f90`,
linked against `qe_pw` (so it can re-use the RDMFT modules for the
TSM probe and the DOS / Mz routines) and built into `bin/rdmft_dos.x`
by both the legacy Make build and the CMake build.

## Files in this directory

| File | Purpose |
|------|---------|
| `h2.step1.in` | `pw.x` input: RDMFT optimisation, saves only natural occupations + natural orbitals.  No DOS computed here. |
| `h2.step2.in` | `rdmft_dos.x` input: dedicated DOS post-processor, reads the saved state and writes `h2.rdmft.dos`. |
| `README.md`   | This file. |

The pseudopotential `H_HSCV_PBE-1.0.UPF` is **not** shipped in this
directory; copy it from `CPV/examples/EXX-wf-example/` (or from
`PW/src/rdmft/examples/h2_molecule/`) before running.

## How to run

```bash
# 0) Get the pseudo
cp ../../../../CPV/examples/EXX-wf-example/H_HSCV_PBE-1.0.UPF .

# 1) RDMFT calculation (pw.x, runs the alternating Müller optimisation)
mkdir -p tmp
pw.x < h2.step1.in > h2.step1.out

# 2) DOS calculation (rdmft_dos.x, the dedicated post-processor) --
#    reads the (n, evc) state saved by step 1; no pw.x, no SCF, no
#    optimisation.
rdmft_dos.x < h2.step2.in > h2.step2.out
```

Inspect the outputs:

```bash
grep -E '!RDMFT_ETOTAL|RDMFT total energy at loaded' h2.step1.out h2.step2.out
# h2.step1.out:     !RDMFT_ETOTAL  =      -2.2429169157 Ry
# h2.step2.out:     RDMFT total energy at loaded (n, evc) =      -2.2429169157 Ry
#                                                                   ^^^^^^^^^^ bit identical

ls h2.rdmft.dos     # produced by step 2 only
head -3 h2.rdmft.dos
# # E(Ry)  DOS  (states/Ry/cell; non-spin-polarized / noncollinear)
#  -1.50000000E+01  0.00000000E+00
#  -1.49800000E+01  0.00000000E+00
```

## What each input does

**`h2.step1.in` — RDMFT calculation.**  Standard alternating
optimisation, `rdmft_outer_maxiter = 8`.  At the end of every outer
cycle (and once unconditionally on exit) the natural occupations are
written to `tmp/h2.rdmft.save`; the natural orbitals are flushed to
`tmp/h2.save/wfc*.dat` (collected wavefunctions + XML schema) via the
automatic `punch('config')` call inside `rdmft_save_state`.  DOS is
computed only in step 2 via `rdmft_dos.x`, so no DOS files are written
during step 1 and the expensive transition-state probe is **not** run.

**`h2.step2.in` — DOS calculation.**  This input is read by
`rdmft_dos.x`, not `pw.x`.  It is a single `&inputrdmftdos` namelist
containing only DOS-relevant knobs:

```
&inputrdmftdos
  prefix             = 'h2'
  outdir             = './tmp'
  functional         = 'muller'   ! must match what step 1 used
  emin               = -1.10      ! Ry; omit (= ±1e6) for auto window
  emax               =   0.37
  deltae             = 0.00147    ! Ry grid step
  degauss            = 0.00368    ! Ry Gaussian broadening
  pdos               = .false.
  occ_weighted       = .false.
/
```

`rdmft_dos.x` internally:

1. Calls `read_file_new(needwf=.TRUE.)` to load the PWscf state
   (lattice, atoms, k-points, charge density, local potential, Fermi
   energy) from `tmp/h2.save/data-file-schema.xml`.
2. Opens the wavefunction buffer (`iunwfc`) and re-fills it from the
   collected wavefunctions `tmp/h2.save/wfc*.dat` written by step 1.
3. Calls `rdmft_load_state` to populate `rdmft_n` from
   `tmp/h2.rdmft.save`.
4. Forces the global XC IDs to `(5,0,0,0,0,0)` and `exxalfa = 1.0`
   so `v_of_rho` returns `V_xc = 0` during the transition-state
   probe, then calls `setup_exx` / `rdmft_exxinit_once` to bring up
   the EXX-via-ACE machinery.
5. Calls `rdmft_compute_dos()` and `rdmft_compute_magnetization()` --
   the TSM probe runs only in this post-processing step.

No SCF, no optimiser; the only RDMFT work that actually runs is the
on-the-fly TSM probe (`O(N_b N_k)` single-k ACE builds) and the
binning into `<prefix>.rdmft.dos`.

## Iterating on DOS parameters

Once `tmp/h2.rdmft.save` and `tmp/h2.save/wfc*.dat` exist, you can
re-run step 2 as many times as you like to scan
`degauss`, `emin / emax`, `deltae`,
`pdos` or `occ_weighted`.  Each step-2 invocation
takes ~1 s on this H2 setup; the (potentially hour-scale) RDMFT
optimisation is **not** repeated.

## What is NOT saved at step 1, by design

The transition-state-model band energies $\varepsilon^{\mathrm{TSM}}_{ik} =
(\partial E/\partial n_{ik})|_{n_{ik}=1/2}$ are the dominant DOS
cost: one ACE rebuild per `(band, k-point)` pair.  Saving them at
every optimiser checkpoint would mean every save call pays an
`O(N_b N_k)` cost.  Instead the optimiser saves only the cheap
`(occupation, orbital)` pair and the DOS step computes the TSM probe
exactly once -- when the user actually asks for the DOS.
