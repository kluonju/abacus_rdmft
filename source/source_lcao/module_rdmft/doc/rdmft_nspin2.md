# RDMFT for collinear spin-polarized calculations (`nspin = 2`)

This note summarises the structural fixes and the spin-resolved equality
constraint added to the ABACUS LCAO RDMFT module so that
`nspin = 2` runs reproduce the KS-LCAO HF energy and gradient identities,
and so the solver can independently constrain `N_↑` and `N_↓` (the same
"split-Fermi" treatment KS uses when `nupdown` is set, see
`docs/advanced/scf/spin.md`).

The companion derivation note `rdmft_derivation.md` § 11 already states
the physical formulas; this file is the implementation companion.

---

## 1. Why the old `nspin = 2` path was wrong

### 1.1 Hartree (and grid-XC) potential was added twice

`Veff_rdmft_local::contributeHR` (in
`rdmft_energy_gradient.cpp`) and the legacy `Veff_rdmft` helper
(`rdmft_tools.cpp`) accumulate the local potential into a **single-spin**
`HContainer<TR>` (`HR_hartree_`, `HR_one_`, `HR_exx_`).  The old code
contained an inner loop

```cpp
for (int is = 0; is < nspin; ++is) {
    vr_eff = &v(is, 0);
    ModuleGint::cal_gint_vl(vr_eff, hR);
}
```

For `nspin = 2`, `PotHartree::cal_v_eff` returns `v(0) = v(1) = V_H`
(the Hartree potential of the total density, identical for both spin
components in the collinear case), so the loop adds `V_H` to the single
shared `hR` **twice**.  Every band matrix element becomes `<ψ|2 V_H|ψ>`,
which doubles `E_H`, `∂E/∂n_iks`, and the orbital block of `∂E/∂C`.
Because both the analytic gradient and finite-difference probes of the
running energy walk on the same buggy potential, internal FD checks
could still pass while the absolute energies disagreed with KS-LCAO HF.

The KS-LCAO `Veff` operator does **not** have this issue because it
holds an `HContainer` with two spin slots (`hRS2`) and the operator
chain calls `contributeHR` once per spin with the data pointer
swapped between slots.  RDMFT's `EnergyGradient` is single-slot by
design, so it must accumulate `V_H` exactly once per call.

**Fix.**  In `rdmft_energy_gradient.cpp::local_build_HR_gint` and
`rdmft_tools.cpp::build_HR_gint`, drop the `is` loop and call
`cal_gint_vl(v(0), hR)` exactly once for both the `"hartree"` and
`"xc"` branches.

### 1.2 Multi-k `DensityMatrix` constructor was off by `nspin`

In `EnergyGradient::build_charge` and `EnergyGradient::build_DM_xc`,
RDMFT created the multi-k density matrix with

```cpp
new elecstate::DensityMatrix<TK, double>(ParaV_, nspin_, kv_->kvec_d, nk_);
```

`nk_` (= `pelec_->wg.nr`) is the **flattened** `(k, spin)` row count, so
for `nspin = 2` it is already `2 * nks_orig`.  The `DensityMatrix`
multi-k constructor expects the per-spin k count (`_DMK` is then sized
`_nk * _nspin`); passing the flattened count over-allocates `_DMK` and
mis-orders `cal_DMR`, so the spin-up DMR receives both spin DMKs while
the spin-down DMR stays empty.  The total density `ρ_↑ + ρ_↓` is still
correct (the spin-up DMR contains the sum) but `m_z = ρ_↑ - ρ_↓` is not.

**Fix.**  Mirror the canonical KS-LCAO pattern from
`source/source_lcao/setup_dm.cpp:19`:

```cpp
const int nk_per_spin = nk_ / nspin_;
new elecstate::DensityMatrix<TK, double>(ParaV_, nspin_, kv_->kvec_d, nk_per_spin);
```

The same correction is applied to the legacy
`update_state_rdmft.cpp::update_charge` path (`nk_total / nspin`).

### 1.3 Single equality constraint cannot pin `(N_↑, N_↓)`

A single chemical-potential multiplier μ enforces

  Σ_k w_k Σ_i n_iks = N_e

but for `nspin = 2` the spin populations are not coupled by the
constraint, so the projection (and the augmented-Lagrangian gradient)
can drift to any feasible split, e.g. `(N_e, 0)` instead of the
KS-converged `(N_↑, N_↓)`.  This is exactly the same issue that
KS-LCAO solves with split Fermi levels when `nupdown` is set
(`PARAM.globalv.two_fermi`): you have *two* equality constraints, you
need *two* multipliers.

---

## 2. Spin-resolved `OccupationConstraint`

`OccupationConstraint` now has a second constructor

```cpp
OccupationConstraint(method, kweights, nbands, nspin, isk, n_electrons_per_spin);
```

When `nspin > 1`:

* `isk[ik] ∈ [0, nspin)` is the spin index for `(k, σ)` row `ik` (read
  from `K_Vectors::isk`).
* Two equality constraints are imposed,

      Σ_{k(s)} w_k Σ_i n_iks = N_s   for s = 0, 1.

  The "scalar" `n_electrons` is just the sum (Σ N_s); it is kept as a
  back-compat hook for the per-band `(λ + μ c) w_k` gradient terms.
* `lambda_per_spin_[s]`, `mu_per_spin_[s]` carry the two multipliers and
  penalty parameters (a single `μ` `set_mu(μ)` broadcasts to both
  spins, mirroring the original API).

### 2.1 Augmented Lagrangian penalty / gradient

The penalty becomes

```
P(n) = Σ_s [ λ_s c_s(n) + ½ μ_s c_s(n)² ],  c_s = Σ_{k(s)} w_k Σ_i n_iks - N_s
```

so the per-band gradient is

```
∂P/∂n_ik = (λ_{σ(ik)} + μ_{σ(ik)} c_{σ(ik)}) · w_k
```

with `σ(ik) = isk[ik]`.  `update_multiplier(occ)` advances each
`λ_s ← λ_s + μ_s c_s(n)` independently, and `increase_penalty` scales
all `μ_s` by the same factor (same convergence schedule per spin).

### 2.2 Projection (clip + bisection) and proximal projection

Both `project` (pre-clip + bisection) and `proximal_project` (no
pre-clip) now run an **independent** 1D dual bisection per spin block:

```
y_iks(α_s) = clip(x_iks - α_s · w_k, 0, 1)        (project / proximal)
y_iks(μ_s) = clip(x_iks - μ_s,        0, 1)       (project_uniform / proximal_uniform)

solve Σ_{k(s)} w_k Σ_i y_iks(α_s) = N_s   independently per s
```

This is the L2 nearest feasible point because the two equality
constraints touch disjoint coordinate blocks, so the joint projection
factorises.

### 2.3 Active set

`identify_active_set` tags lower / upper / free coordinates the same
way as before, then computes a separate Lagrange-multiplier estimate
per spin from the free indices in that spin block:

```
λ_s = - (Σ_{i ∈ Free, σ(i) = s} (∇E)_i w_{k(i)})
       /  (Σ_{i ∈ Free, σ(i) = s} w_{k(i)}²)
```

`apply_active_set` adds `λ_{σ(i)} w_{k(i)}` to the gradient on free
indices, so every spin block ends up orthogonal (in the `w`-weighted
sense) to its own electron-sum constraint row.  The scalar
`info.lagrange_mult` is populated with `max_s |λ_s|` for diagnostic
logging; downstream consumers that need the per-spin values read
`info.lagrange_mult_per_spin`.

---

## 3. Spin-resolved `SigmaShiftOccParam`

The sigma-shift parameterisation (joint solver only) now solves a
**separate** scalar shift λ_s per spin so that

```
Σ_{k(s)} w_k Σ_i σ(z_iks + λ_s) = N_s        for s = 0, 1
```

independently.  The map `z → n` uses the matching λ_{σ(ik)}, and the
implicit-derivative chain rule

```
∂E/∂z_iks = σ'(z_iks + λ_s) · ( ∂E/∂n_iks - w_k · (Σ h σ' / Σ w σ')_s )
```

is computed per spin block.

When `nspin = 1` the class falls back to the original single-shift
constructor / single ratio.

---

## 4. Driver wiring

* `RDMFTNelectronTargetMeta` now carries `two_fermi_active`,
  `n_electrons_per_spin`, and `nupdown` so the LCAO driver can pass the
  per-spin targets through to the solver.
* `ESolver_KS_LCAO::after_scf` populates per-spin targets when
  `PARAM.inp.nspin == 2`:
  * Prefer the KS occupation seed `N_s = Σ_{ik(s), ib} pelec->wg(ik, ib)`
    (rescaled to the RDMFT total `N_e` via `rdmft_nelec_delta`) so the
    solver locks onto the same magnetisation KS converged to.
  * If no KS seed is available (e.g. `scf_nmax = 0`) or
    `rdmft_nelec_use_input` is `true`, fall back to the symmetric
    `nupdown` split `N_↑ = (N_e + nupdown) / 2`,
    `N_↓ = (N_e - nupdown) / 2`, mirroring KS's `init_nelec_spin`.
* `RDMFTSolver::init`:
  * Detects `nspin = 2` from the meta block (or as a safety fallback,
    from `PARAM.inp.nspin == 2` + a consistent `kv_->isk`).
  * Builds the per-spin `isk` from `kv_->isk` and constructs the
    spin-resolved `OccupationConstraint` (and `SigmaShiftOccParam` when
    `rdmft_occ_param sigma_shift` is selected).
  * Logs a single `RDMFT spin constraint:` line summarising whether the
    solver is in single- or two-Fermi mode.

No new INPUT keywords are required: `nspin` and `nupdown` reuse the
existing KS controls, and the per-spin targets are derived from the KS
seed (or the symmetric split) automatically.

---

## 5. Verification

Built with

```
cmake -B build -DENABLE_RDMFT=ON -DUSE_ELPA=OFF -DBUILD_TESTING=OFF \
                -DCMAKE_CXX_COMPILER=/usr/bin/g++
cmake --build build -j4
```

### 5.1 Closed-shell `nspin = 1` ↔ `nspin = 2` invariance

For a closed-shell system with integer occupations 1 / 0 (HF Aufbau) and
the same KS reference, the RDMFT objective must give the **same** total
energy independent of `nspin`. Pre-fix, `nspin = 2` was systematically
higher than `nspin = 1` because the `V_H` double-add inflated `E_H` by
exactly 2× for `nspin = 2` (Section 1.1). With the fix in place the
identity is recovered to within optimisation tolerance.

LiH (gamma-only) with `dft_functional = pbe` for the KS reference, then
RDMFT-HF on top (`rdmft_outer_maxiter 50`, `rdmft_occ_maxiter 5`,
`rdmft_orb_maxiter 5`, `rdmft_energy_tol 1e-7`,
`rdmft_orb_grad_tol 1e-5`, `rdmft_occ_grad_tol 1e-5`,
`rdmft_occ_init_mode = ks` so occupations stay at 1 / 0):

| Quantity                | `nspin = 1`              | `nspin = 2`              |
| ----------------------- | ------------------------ | ------------------------ |
| `E_one_elec` (Ry)       | -21.32323095             | -21.32315340             |
| `E_Hartree` (Ry)        |   9.18592165             |   9.18579375             |
| `E_xc` (Ry)             |  -3.84473043             |  -3.84468006             |
| `E_total` (Ry)          | **-15.19029180**         | **-15.19029177**         |
| `E_total` (eV)          | -206.6745227             | -206.6745224             |

Both stop on `|dE| < rdmft_energy_tol = 1e-7 Ry`; the residual gap is
3 × 10⁻⁷ eV (≈ 2 × 10⁻⁸ Ry), far below the convergence tolerance.

H₂ (gamma-only) with the same recipe and tighter tolerances
(`rdmft_outer_maxiter 20`, `rdmft_orb_maxiter 30`,
`rdmft_orb_grad_tol 1e-7`, `rdmft_energy_tol 1e-9`):

| Quantity                | `nspin = 1`              | `nspin = 2`              |
| ----------------------- | ------------------------ | ------------------------ |
| `E_total` (Ry)          | -2.0945324138            | -2.0945324138            |
| `E_total` (eV)          | -28.4975754736           | -28.4975754737           |

Match to ~10⁻¹⁰ Ry (machine precision).

Pre-fix the same H₂ closed-shell run gave `E_total(nspin=2)` ≈
`E_total(nspin=1) + E_H` (a ~12 eV gap on H₂, ~125 eV on LiH); seeing
that "`nspin = 2` has higher energy" with otherwise correct integer 1 / 0
occupations is the canonical fingerprint of the V_H double-add bug.

### 5.2 Open-shell H₂ (KS-HF reference)

A gamma-only `nspin = 2` H₂ test with KS-HF SCF followed by an
RDMFT-HF one-shot evaluation (4 OMP threads, 1 MPI rank):

| Quantity                               | KS-HF (`dft_functional=hf`) | RDMFT-HF (`rdmft_functional=hf`) one-shot |
| -------------------------------------- | --------------------------- | ------------------------------------------ |
| `E_one_elec` (eV)                      | -22.392                     | -22.392                                    |
| `E_Hartree` (eV)                       | 12.055                      | 12.055                                     |
| `E_xc` (eV) / RDMFT exchange (eV)      | -0.577 (LDA-corr residue) + -13.599 (E_exx) | -13.599 (single E_xc bucket)               |
| `E_Ewald` (eV)                         | -4.551                      | -4.551                                     |
| `E_total` (eV)                         | -29.064                     | -28.487 (= KS - LDA residue)               |

The 0.577 eV gap is the LDA-correlation residue ABACUS still
accumulates from the LDA pseudopotential's local XC even when
`dft_functional=hf`; RDMFT only contributes the EXX-style exchange and
hence reports the energy with that residue removed.  When `nspin = 1`
runs under the same setup the same identity holds (RDMFT one-shot E
matches `E_KohnSham - E_xc(LDA-residue)` to within SCF tolerance), so
the two values now scale identically across `nspin = 1` / `nspin = 2`,
recovering the spin-equivalence the buggy path violated.

`rdmft_grad_check` reports `PASS` for every probed occupation
(rel_err = 0 to floating-point precision against central differences)
and for the directional orbital derivative (rel_err ≈ 1.5e-5, set by
the FD truncation; the canonical-metric Riemannian gradient is
recomputed and compared against the same retraction the line search
uses).

For `nspin = 1` the same gamma-only run still passes the gradient
check (rel_err ≈ 3e-5 in orbitals, 0 in occupations) and the
single-multiplier path is untouched.

---

## 6. Where to look in the source

| Concern                                    | File                                                               |
| ------------------------------------------ | ------------------------------------------------------------------ |
| `V_H` single-add fix                       | `rdmft_energy_gradient.cpp::local_build_HR_gint`, `rdmft_tools.cpp::build_HR_gint` |
| `DensityMatrix` multi-k argument           | `rdmft_energy_gradient.cpp::build_charge`, `EnergyGradient::build_DM_xc`, `update_state_rdmft.cpp::update_charge` |
| Spin-resolved constraint type              | `rdmft_occupation.h::OccupationConstraint`, `SigmaShiftOccParam`   |
| Per-spin meta plumbing                     | `rdmft_type.h::RDMFTNelectronTargetMeta`                           |
| Driver wiring (`nspin = 2` detection)      | `esolver_ks_lcao.cpp::after_scf`, `rdmft_solver.cpp::RDMFTSolver::init` |

---

## 7. Reproducer (for regressions)

```text
INPUT_PARAMETERS
calculation       scf
gamma_only        1
nspin             2
basis_type        lcao
dft_functional    hf

rdmft             1
rdmft_functional  hf
rdmft_outer_maxiter 1
rdmft_occ_maxiter   0
rdmft_orb_maxiter   0
rdmft_grad_check    1
```

Run with `OMP_NUM_THREADS=4 mpirun -np 1 abacus_std_para` from
`source/source_lcao/module_rdmft/example/H2_HF/` (any STRU with
two electrons works; use `dft_functional=hf` to land on a KS-HF
reference).  Look for

```
RDMFT spin constraint: nspin=2, two equality constraints active. ...
RDMFT init occupations: ... constraint_after_project ≈ 0
Gradient check PASSED
```

Pre-fix this would either fail the gradient check or report `E_H` /
`E_total` ~ 2× the KS reference and a `(N_↑, N_↓)` projection that drifts
away from KS.
