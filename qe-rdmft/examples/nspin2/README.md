# RDMFT-HF on spin-polarised inputs (`nspin = 2`)

## What this example tests

Verifies that the RDMFT-HF post-SCF stationary point coincides with
the converged KS-HF total energy when the SCF was run with
`nspin = 2`, both for closed-shell (`tot_magnetization = 0`) and
spin-polarised (`tot_magnetization /= 0`) inputs, and from rough
non-HF starting orbitals as well.

The script `run.sh` runs five cases:

  1. H2 in a 12-Bohr cubic Martyna-Tuckerman cell, `nspin = 2`,
     `tot_magnetization = 0`, `input_dft = 'hf'` SCF
     (`conv_thr = 1e-8`) followed by RDMFT-HF.
  2. Diamond-Si 2x2x2 multi-k crystal (3 irreducible k-points),
     `nspin = 2`, `tot_magnetization = 0`, `input_dft = 'hf'`
     SCF (`conv_thr = 1e-7`) followed by RDMFT-HF.
  3. H2^- (`tot_charge = -1`, 3 electrons), `nspin = 2`,
     `tot_magnetization = 1`, `input_dft = 'hf'` SCF followed by
     RDMFT-HF.  This is the non-trivial spin-polarised case where
     spin up has 2 occupied bands and spin down has 1.
  4. H2 closed shell `nspin = 2`, `tot_magnetization = 0`, but
     started from a *pure-PBE* SCF stopped after only 2 inner
     iterations (no Fock character whatsoever), then RDMFT-HF.
     Tests the on-the-fly EXX bring-up path together with
     spin-polarisation bookkeeping.
  5. Si 2x2x2 closed shell `nspin = 2`, `tot_magnetization = 0`,
     same as (4): rough pure-PBE SCF then RDMFT-HF.

## Expected result

```
case                                           KS-HF / KS-PBE         RDMFT-HF       Mz          dE_KS
-------------------------------------------- ---------------- ---------------- -------- --------------
H2 nspin=2 Mz=0, HF -> RDMFT-HF                   -2.24112743    -2.2411274288    -0.00      +1.20e-09
Si 2x2x2 nspin=2 Mz=0, HF -> RDMFT-HF            -15.13128573   -15.1312857709    -0.00      -4.09e-08
H2^- nspin=2 Mz=1, HF -> RDMFT-HF                 -2.24197631    -2.2419763535     1.00      -4.35e-08
H2 nspin=2 Mz=0, PBE(2it) -> RDMFT-HF             -2.24112743    -2.2411274288    -0.00      +1.20e-09
Si 2x2x2 nspin=2 Mz=0, PBE(2it) -> RDMFT-HF      -15.13128573   -15.1312875306    -0.00      -1.80e-06
```

All five cases reach the converged KS-HF total energy to better
than ~2 µRy (~ 25 nK) and conserve the requested magnetization
exactly.

## Notes / known limitation

* PWscf disables ACE when one spin channel is *empty*
  (`x_nbnd_occ = 0` for that spin), because the
  `MatChol(-mexx)` step in `aceupdate` breaks down on a zero
  Fock matrix.  **RDMFT Müller/channel exchange** can likewise
  build a **numerically zero** `⟨φ|Vx|φ⟩` block when that spin’s
  channel weights (`wg ∝ √(n)`, etc.) are all below `eps_occ`, so no
  EXX pairs contribute; PW now skips Cholesky in that case (`Vx ≡ 0`
  on that slice). This means cases like H2 with
  `tot_magnetization = 2` (two same-spin electrons, empty
  opposite-spin channel) cannot be run with the current
  ACE-based RDMFT exchange path; the same input also fails in
  pure KS-HF SCF for the same reason.  All cases here keep both
  spin channels populated.

## Running

```bash
cd PW/src/rdmft/examples/nspin2
ln -sf <path-to-pseudos>/H_HSCV_PBE-1.0.UPF .
ln -sf <path-to-pseudos>/Si.upf .
ln -sf <path-to-build>/bin/pw.x .
PWX=./pw.x ./run.sh
```
