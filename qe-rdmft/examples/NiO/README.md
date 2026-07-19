# NiO RDMFT DOS — ELK 10.9.5 cross-check

Matched to [`elk-10.9.5/examples/RDMFT/NiO/elk.in`](../../../../elk-10.9.5/examples/RDMFT/NiO/elk.in):

| Setting | ELK | QE (`&inputrdmftdos`) |
|---------|-----|----------------------|
| Structure | 2-atom primitive NiO | same (`rdmft.in`) |
| k-mesh | 4×4×4, vkloff 0.25/0.5/0.625 | `K_POINTS automatic 4 4 4 2 4 5` |
| RDMFT XC | power, α=0.65 | `functional='power'`, `power_alpha=0.65` |
| DOS window | −0.6…0.6 Ha | `emin/emax = -1.2/1.2` Ry |
| Grid | nwplot=500, ngrkf=100 | same |
| dosocc | `.false.` | `occ_weighted=.false.` |
| Spectral | ELK `dos.f90` | `spectral='elk'` (default) |
| Integration | brzint | `integration='brzint'` (default) |

## Workflow

```bash
# 1) RDMFT optimisation (no DOS inside pw.x)
pw.x < rdmft.in > rdmft.out

# 2) ELK-mode DOS post-processor
rdmft_dos.x < dos.in > dos.out

# 3) Compare against ELK TDOS.OUT (run ELK separately in elk-10.9.5/examples/RDMFT/NiO/)
python3 compare_elk_dos.py --elk TDOS.OUT --qe nio.rdmft.dos
```

Point `--elk` at the `TDOS.OUT` produced by ELK tasks 0+300+10 on the same
k-mesh and functional.  Exact agreement is not expected (LAPW vs PP basis), but
peak positions and gross shape should align once the ELK spectral model is used.

Optional: set `write_evalsv = .true.` to emit `nio.rdmft.evalsv` for
direct TSM-energy comparison with ELK `EIGVAL.OUT`.

## Pseudopotentials

Provide Ni and O UPF files with `PP_CHI` blocks if partial DOS is needed.
Symlink from another example, e.g.:

```bash
ln -sf ../NiO_nk1x1x1/Ni.upf .
ln -sf ../NiO_nk1x1x1/O.upf .
```

(Adjust paths if those files live elsewhere on your system.)
