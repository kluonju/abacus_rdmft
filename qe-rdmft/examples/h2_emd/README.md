# RDMFT on H2 — EMD / Compton / autocorrelation workflow (`pw.x` → `rdmft_emd.x`)

Two-step workflow mirroring the DOS example, but for the electron momentum
density (EMD), directional Compton profiles \(J(q)\), and autocorrelation
functions \(B(r)\).

| Step | Program | Purpose |
|------|---------|---------|
| 1 | `pw.x` | RDMFT optimisation; writes `tmp/h2.rdmft.save` and `tmp/h2.save/` |
| 2 | `rdmft_emd.x` | Momentum density + optional Compton profile and \(B(r)\) |

## Requirements

- A **uniform** Monkhorst–Pack mesh with **zero offset** (`0 0 0` in
  `K_POINTS automatic`).  Symmetry-reduced meshes are unfolded internally,
  but k-pool parallelization is not supported.
- The momentum grid is fixed by **`ecutwfc`** from the `pw.x` run (read from
  the saved state; no separate cutoff in `&inputrdmftemd`).
- Pseudopotential orbitals give **pseudo** momentum densities, not
  all-electron Compton profiles.

## Run

```bash
# copy H_HSCV_PBE-1.0.UPF into this directory (see h2_dos_two_step example)
pw.x < h2.step1.in > h2.step1.out
rdmft_emd.x < h2.step2.in > h2.step2.out
```

## Outputs

| File | Contents |
|------|----------|
| `h2.rdmft.emd` | Per-k \(n(p)\) on the \(G+k\) sphere up to `ecutwfc` |
| `h2.rdmft.emd_compton` | Directional Compton profile \(J(q)\) vs \(q\) (Bohr\(^{-1}\)) |
| `h2.rdmft.emd_af` | Directional autocorrelation \(B(r)\) vs \(r\) (Bohr) |
| `h2.rdmft.emd1d` | Legacy parametric-line Compton profile (`do_emd_plot1d = .true.`) |
| `h2.rdmft.emd_band` | Optional band-resolved \(n_i(p)\) (`write_band_resolved = .true.`) |
| `h2.rdmft.emd_compton_band` | Per-band Compton \(J_i(q)\) (always with `do_emd_compton = .true.`) |
| `EMD.OUT` / `EMD1D.OUT` / `EMDCOMPTON.OUT` | Optional ELK-compatible text (`write_elk_fmt = .true.`) |

The stdout line `sum_k w_k sum_G n(p)` should match the electron count (2 for H2).
The Compton header reports `integral(J) dq` and autocorr reports `B(0)`; both
should be close to the selected-band electron count.

Set `band_min` / `band_max` in `&inputrdmftemd` to restrict outputs to a
natural-orbital band window (e.g. valence-only Compton by setting
`band_min = 1`, `band_max = N_valence`).

### Band-resolved Compton (valence / core separation)

Per-natural-orbital Compton profiles \(J_i(q)\) are written to
`*.emd_compton_band` whenever `do_emd_compton = .true.`.  Each band block is
preceded by `# band ib  n_ib = ...` with the k-weighted occupation.  Data rows:
`ib  q  J_ib(q)  [spin]`.

The optional sparse momentum-density file `*.emd_band` is controlled separately
by `write_band_resolved = .true.`.

To obtain a **valence-only** summed profile, set `band_min` / `band_max`.
To decompose within the window, sum selected bands from `*.emd_compton_band`
in post-processing.

```fortran
  do_emd_compton = .true.
  band_min       = 1
  band_max       = 2    ! H2: both occupied valence bands
  ! optional sparse n_ib(p) on G+k sphere:
  ! write_band_resolved = .true.
```

## Compton and autocorrelation keywords

Recommended input for physical-\(q\) Compton plus \(B(r)\) along a direction:

```fortran
&inputrdmftemd
  do_emd_compton  = .true.
  do_emd_autocorr = .true.
  line_start      = 0.0 0.0 0.0
  line_end        = 1.0 0.0 0.0
  cp_nq           = 256
  cp_qmax         = 0.0          ! <= 0 uses |p|_max from EMD
  af_r_start      = 0.0
  af_r_end        = 6.0
  af_npt          = 150
/
```

- **`line_start` / `line_end`**: crystallographic direction in reciprocal
  coordinates (e.g. `(1,1,1)` for `<111>`).
- **`cp_qmax`**: upper bound on \(q\) in Bohr\(^{-1}\); grid runs
  \(q \in [\texttt{cp\_qmin}, q_\mathrm{max}]\) with default `cp_qmin = 0`.
- **`af_r_start` / `af_r_end`**: real-space range in Bohr for \(B(r)\).

Relationship: \(B(r) = \int J(q)\, \cos(qr)\, dq\) with even \(J(q)\); the
\(q<0\) contribution is included via a factor of 2 when the grid starts at
\(q=0\).

## Plotting (gnuplot)

```gnuplot
# Compton profile
plot "h2.rdmft.emd_compton" u 1:2 w l title "J(q)"

# Autocorrelation function
plot "h2.rdmft.emd_af" u 1:2 w l title "B(r)"
```

Legacy `.emd1d` uses a parametric coordinate along `line_start`→`line_end`
(not physical Bohr\(^{-1}\)); prefer `.emd_compton` for new workflows.
