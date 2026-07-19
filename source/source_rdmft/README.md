# `source_rdmft` — modular, basis-independent RDMFT core

This module implements Reduced Density Matrix Functional Theory (RDMFT) as a
self-contained optimisation core that is **independent of the plane-wave / LCAO
Hamiltonian machinery**.  It is the shared engine behind the new
`ESolver_RDMFT` (selected with `esolver_type = rdmft`) and is designed so the
same occupation and orbital optimisers drive both the plane-wave and the LCAO
backends, for gamma-only (real, `double`), collinear (`nspin = 2`) and
non-collinear (`nspin = 4`) spin.

## Design

Everything the optimisers need is expressed through a single abstract oracle,
`rdmft::RdmftBackend` (`rdmft_backend.h`):

```
              +-----------------------------+
              |        RdmftDriver          |   alternating outer loop
              +--------------+--------------+
                             |
         +-------------------+-------------------+
         |                                       |
+--------v--------+                     +--------v---------+
|  OccOptimizer   |  SPG2 / EBI         | OrbitalOptimizer |  Stiefel SD/CG/LBFGS
+--------+--------+                     +--------+---------+
         |                                       |
         +-------------------+-------------------+
                             |
                    +--------v---------+
                    |  RdmftBackend    |  <-- abstract energy/gradient oracle
                    +--------+---------+
                             |
            +----------------+-----------------+
            |                                  |
   +--------v---------+              +---------v----------+
   | RdmftBackendLCAO |              |  (PW backend:      |
   | (wraps           |              |   pluggable via    |
   |  EnergyGradient) |              |   RdmftBackend)    |
   +------------------+              +--------------------+
```

Occupation vectors are flat arrays `occ[ib + ik*nbnd]` (band fastest); orbital
tangent vectors are opaque flat `double` buffers whose interpretation
(real/complex, k-blocking, gamma trick, AO overlap `S`) is entirely internal to
the backend.  This keeps the algorithms identical across bases and spin cases.

## Files

| File | Contents |
|------|----------|
| `rdmft_params.h` | Configuration struct and enum/string parsing. |
| `rdmft_xc.{h,cpp}` | Natural-orbital XC channel couplings (HF, Müller, Power, GU, CHF, CGA, GEO, HybOpt, BOWMOD). |
| `rdmft_occ_constraints.{h,cpp}` | Euclidean L2 proximal projector `P_w`, SPG2 stopping map, KKT residual, EBI erf parameterisation + implicit `mu`. |
| `rdmft_line_search.{h,cpp}` | Armijo, strong-Wolfe, Barzilai–Borwein spectral step. |
| `rdmft_occ_optimizer.{h,cpp}` | SPG2 and EBI occupation blocks. |
| `rdmft_orbital_optimizer.{h,cpp}` | Stiefel SD/CG/L-BFGS orbital block. |
| `rdmft_backend.h` | Abstract energy/gradient oracle. |
| `rdmft_driver.{h,cpp}` | Alternating outer loop. |

## Provenance

The algorithms are ported from the Quantum ESPRESSO RDMFT reference bundled in
`qe-rdmft/`:

- SPG2 — `rdmft_spg.f90` (Birgin–Martínez–Raydan spectral projected gradient).
- EBI  — `rdmft_ebi.f90` (Yao et al. 2022 explicit-by-implicit erf map).
- Projector / stopping map — `rdmft_occupation.f90`.
- Stiefel geometry — `rdmft_stiefel.f90`, `rdmft_solver.f90`.
- XC channels — `rdmft_xc.f90`.

## Tests

`test/` contains gtest unit tests with toy analytic backends that verify:

- XC channel derivatives against finite differences;
- the projector / constraint machinery and EBI erf/`mu` round trips;
- SPG2 and EBI convergence to an independently computed KKT point;
- Stiefel SD/CG/L-BFGS convergence to the trace-minimising eigenspace with
  orthonormality preserved.
