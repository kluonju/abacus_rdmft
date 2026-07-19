# `examples/Fe_bcc_dos/` — BCC Fe RDMFT density of states smoke test

Two-step exercise of the ELK-style RDMFT DOS / Mulliken magnetization
post-processor via **`rdmft_dos.x`**.  Spin-polarised norm-conserving BCC Fe
(`Fe.pz-n-nc.UPF`, ferromagnetic, 2×2×2 k mesh, `nbnd = 10`).

## Run

```bash
pw.x -in fe.scf.in > fe.out
rdmft_dos.x < fe.dos.in > fe.dos.out
```

The first step runs SCF (+ optional zero-iteration RDMFT save).  The second
step reads the saved state and writes:

| File | Contents |
|------|----------|
| `fe.rdmft.dos`            | total DOS (spin-up positive, spin-down negative, Ry) |
| `fe.rdmft.pdos_at001_001` | l-resolved partial DOS at atom 1 (`Fe`) |
| `fe.rdmft.mag`            | Mulliken Mz per site and per species |

## What this demonstrates

* The Mulliken Mz of the single Fe site (≈ +1.97 e⁻) matches the SCF
  `total_magnetization` of 1.97 μB / cell.  Site Mz - global Mz < 1 %.
* The RDMFT TDOS shows the characteristic exchange-split Fe d band:
  the minority-spin (down) peak sits below the Fermi level (~ −1.3 eV),
  the majority-spin (up) peak is fully below it.  See `fe.rdmft.dos`.
* PDOS integrated over occupied energies gives ≈ 3.77 e⁻ in d↑ vs.
  ≈ 1.99 e⁻ in d↓ (the s and p channels supply the small remainder
  needed to reach the Bader spin moment).

## Relation to the literature

The DOS post-processor is patterned after the ELK RDMFT driver used in
Sharma, Dewhurst, Lathiotakis, Gross, _Phys. Rev. B_ **78**, 201103(R)
(2008) and the spectral-density analysis in Sharma, Dewhurst, Shallcross,
Gross, _Phys. Rev. Lett._ **110**, 116403 (2013) (arXiv:0912.1118).  Those
calculations used the LAPW basis of ELK plus dense k meshes and the
Power/BB1 functionals; this smoke test uses a small NC pseudopotential,
a 2×2×2 k mesh and `rdmft_outer_maxiter = 0` (DOS from the SCF orbitals,
no RDMFT update), so the spectra here cannot be compared one-to-one with
the PRL data.  The goal is to verify the QE-side machinery (atomic-wfc
projections, brzint TDOS / PDOS, per-site Mulliken Mz) and to provide a
starting point for full RDMFT-DOS production runs.

## Notes on cost

The transition-state DOS probe costs `O(nbnd * nks)` gradient evaluations
and runs only in **`rdmft_dos.x`**, not during `pw.x`.  Tune DOS windows and
options in `fe.dos.in` / `&inputrdmftdos`.
