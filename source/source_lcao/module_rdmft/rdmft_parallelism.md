# RDMFT Parallelism

Status summary of the MPI parallelism used inside `module_rdmft`, what is
shared with the rest of ABACUS, and what would be required to extend to true
k-point or band parallelism.

## What is parallel today

RDMFT inherits the standard ABACUS LCAO parallel layout:

* **Basis / matrix (row-column) parallelism** via `Parallel_Orbitals` + BLACS.
  All k-dependent matrices (`H_TV(k)`, `H_H(k)`, `H_xc(k)`, `S(k)`, `DM(k)`,
  the per-band blocks of `psi(k)` and `H * psi(k)`, and the `<psi|H|psi>_ij`
  blocks in `para_Eij`) are stored in the 2D block-cyclic layout described by
  `ParaV->desc` and `para_Eij.desc`. All per-k GEMMs go through ScaLAPACK
  (`pdgemm` / `pzgemm`); Cholesky-QR retraction uses `pdpotrf`/`pzpotrf` +
  `pdtrsm`/`pztrsm` across the same grid.
* **Real-space Hamiltonian parallelism** through `hamilt::HContainer<TR>`
  which owns the set of `<I,J,R>` atom pairs assigned to the rank via the
  same `Parallel_Orbitals`. Grid integrals (`ModuleGint::cal_gint_vl`,
  `cal_gint_rho`) and the EXX density/operator builders reuse the LCAO
  infrastructure.
* **Charge / symmetry / Ewald reduction** is performed by the shared ABACUS
  code paths invoked from `update_charge` / `EnergyGradient::build_charge`.

## What is *not* parallel today

* **k-point parallelism.** ABACUS LCAO disables `kpar > 1` globally (see
  `source_io/module_parameter/read_input_item_system.cpp` line 326), so every
  rank holds the full k-point list. Inside RDMFT, every loop over `ik` is
  therefore serial across ranks even though the BLAS kernels inside each `ik`
  are distributed.
* **Explicit band parallelism.** Bands are already distributed through
  `ParaV->ncol_bands` (column-cyclic block layout on the same 2D grid used
  for the basis), so there is no extra outer-level band MPI split beyond what
  ScaLAPACK provides.

## Upgrade path

A clean extension to (k, band) parallelism consists of two steps that can be
done independently.

### 1. k-point pools

Drop the LCAO `kpar == 1` restriction for the RDMFT code path and iterate
only over the pool-local k-points (`kv->get_nks()`), then reduce energies
and gradients over `MPI_COMM_WORLD` (not the intra-pool `POOL_WORLD`) at
the end of each `compute()` / `cal_Hk_Hpsi` call.

Touch points for the legacy path (`rdmft.cpp`):
* `init()` — remove the `nk_total *= nspin` over-count and keep an explicit
  `nks_local = kv->get_nks()` and `nks_global`.
* `cal_Hk_Hpsi()` — loop is already per-ik, just needs the outer range to be
  `nks_local` and a final `Parallel_Reduce::reduce_all(wfcHwfc_*)` across
  `MPI_COMM_WORLD` before the energy accumulator.
* `cal_Energy()` — energies are already reduced with `reduce_all`, but the
  call happens inside `POOL_WORLD`; move to `MPI_COMM_WORLD` for kpar > 1.

Touch points for the new engine (`rdmft_energy_gradient.cpp`):
* `compute()` — `grad_occ` is built per `ik` locally; the orbital gradient
  `grad_wfc` is already distributed band × basis but replicated over k. With
  k pools, each rank's `grad_wfc(ik, ib, mu)` only contains the k slice it
  owns; reduce `E_one_`, `E_hartree_`, `E_xc_` across `MPI_COMM_WORLD` (the
  `compute_diagonal` reduction inside a pool is unchanged).
* `build_charge()` and `build_DM_xc()` — these already use the ABACUS
  DensityMatrix infrastructure, which is designed to work inside a k pool.
* `get_SK()` cache is per `ik` — it remains correct because only the pool
  that owns `ik` will ask for its `SK`.

Once this is in place, the LCAO constraint (`kpar > 1 has not been
supported for lcao calculation`) needs to be relaxed for the RDMFT branch
specifically, otherwise `GlobalV::KPAR` is forced to 1 at input parsing.

### 2. Band parallelism beyond the 2D block cycle

The current 2D layout already distributes bands along the column axis
(`ParaV->ncol_bands`). For the *RDMFT-specific* loops (e.g. the per-band
`for (int ib_local = 0; ib_local < nb_local; ++ib_local)` loop inside
`compute()` that assembles `grad_wfc`) the work is already proportional to
`nbands / ncol_procs`. What is still serial is the `ik` loop and the
per-band accumulation of `E_one_ / E_hartree_ / E_xc_` where the diagonal
comes from `compute_diagonal`'s Allreduce (small, `O(nbands)`).

A band-parallel extension would mean splitting the 2D grid into a (row,
col, band) tensor parallelism. The cleanest way to do so without touching
ScaLAPACK is to run multiple BLACS contexts simultaneously, one per band
subset, and reduce the energy contributions at the end. The existing
`para_Eij` context makes this straightforward because band diagonals are
already computed independently per BLACS column.

## Redundancy cleanup in this change

* Removed never-read buffers: `HK_XC`, `DM_XC_pass`, `H_wfc_XC`,
  `wfcHwfc_XC`, and the legacy `HK_RDMFT_pass` / `HK_XC_pass` scaffolding.
* Removed `add_wfcHwfc()` (no callers remain after the refactor to
  `add_occNum`).
* Unified the `gamma_only` and multi-k branches in
  `rdmft_pot.cpp::cal_V_TV / cal_V_hartree / cal_V_XC`,
  `update_state_rdmft.cpp::update_charge`, and the new engine's
  `build_charge` / `build_DM_xc` (the two branches were textually
  identical apart from the `DensityMatrix` constructor signature).
* Factored the triple-identical `contributeHR()` specialisations of
  `Veff_rdmft` and `Veff_rdmft_local` into a single `build_HR_gint`
  helper shared by all template instantiations.
* Replaced the per-term `HkPsi / cal_bra_op_ket / _diagonal_in_serial`
  sequence in `cal_Hk_Hpsi()` with a single `diag_action` lambda applied
  to each of TV, Hartree, DFT-XC, and EXX. Same semantics, ~4x less
  code.
