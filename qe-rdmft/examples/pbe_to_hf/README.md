# PBE → RDMFT-HF: recovering the Hartree-Fock minimum from
# non-Hartree-Fock starting orbitals

## What this example tests

Many production workflows already have a converged DFT (PBE) charge
density and orbitals on disk, and the user would like to know
whether the RDMFT-HF solver can use those orbitals as a *starting
point* and walk all the way to the proper Hartree-Fock minimum -- in
other words, whether the RDMFT-HF Stiefel/Riemannian descent is
robust against being seeded with orbitals that have **no Fock
character** at all.

The script `run.sh` runs the following test for two systems
(molecular H2 and a 2x2x2 multi-k Si crystal):

  1. Run PBE SCF for **only 2 inner iterations** (so it stops well
     before convergence).  This gives a very rough, fully PBE-
     flavoured set of starting orbitals.  Pure PBE has *no* exact
     exchange in it, so the EXX module is dormant on entry.
  2. Hand control to `rdmft_run`, which:
        * forces the global `exx_fraction` to 1.0 (full Fock),
        * switches the XC IDs to `(5,0,0,0,0,0)` (HF) so that
          `v_of_rho` returns `V_xc^semilocal = 0` for the rest of
          the RDMFT solve,
        * if necessary, calls `setup_exx` so that the EXX
          symmetry / FFT machinery is set up on the fly (this is
          the work that QE's hybrid-loop driver normally does on
          the first hybrid SCF iteration), and
        * runs the alternating RDMFT-HF optimiser starting from
          those PBE orbitals.
  3. Compares the converged RDMFT-HF energy with a tightly-
     converged KS-HF reference (`input_dft='hf'`, `conv_thr=1e-8`).

## Expected result

```
KS-HF reference:  H2 = -2.24112743 Ry,  Si = -15.13128573 Ry

case (PBE 2 SCF iters -> RDMFT-HF)         RDMFT-HF              dE_KSHF
H2 in 12-Bohr cubic, MT, gamma          -2.2411274288         +1.20e-09
Multi-k Si (2x2x2)                     -15.1312871797         -1.45e-06
```

Both systems converge back to the proper KS-HF stationary point to
better than ~1 µRy, despite being seeded with orbitals that were
optimised for an entirely different functional.  This demonstrates
that:

  * the RDMFT-HF Riemannian gradient + Stiefel retraction descent
    is a true minimiser of the HF energy on the orbital manifold,
    not just a small local correction around the SCF point;
  * the implementation correctly handles the
    `exxbuff` / `x_occupation` refresh per energy + gradient
    evaluation, so the optimiser is using the *current* orbitals
    everywhere (an earlier "single-shot exxinit" caching strategy
    silently produced energies below the variational HF minimum
    because `vexx`'s exxbuff was stale -- now fixed);
  * the multi-k orbital block correctly restores `evc` /
    `current_k` after the per-channel walk over k inside
    `rdmft_compute_xc_channel`, so the Stiefel tangent projection
    acts at the right `C^k` (an earlier silent bug used the last
    k-point's orbitals for the projection at every k).

## Running

```bash
cd PW/src/rdmft/examples/pbe_to_hf
ln -sf <path-to-pseudos>/H_HSCV_PBE-1.0.UPF .
ln -sf <path-to-pseudos>/Si.upf .
ln -sf <path-to-build>/bin/pw.x .
PWX=./pw.x ./run.sh
```
