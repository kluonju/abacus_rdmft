# Reduced Density Matrix Functional Theory (RDMFT): Mathematical Derivation

This note is written to match the **ABACUS LCAO implementation** in
`source/source_lcao/module_rdmft` (energy/gradients in `rdmft_energy_gradient.*`,
XC kernels in `rdmft_xc_functional.h`, manifold ops in `rdmft_stiefel.h`,
constraints in `rdmft_occupation.h`, outer loops in `rdmft_solver.cpp`). Where
the physics literature uses a compact integral notation, we spell out the same
weighting and signs as in the code (including optional HF occupation entropy and
mixed EXX channels).

## 1. One-Body Reduced Density Matrix

The one-body reduced density matrix (1-RDM) in spectral representation is

$$
\gamma(\mathbf{r}, \mathbf{r}') = \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_{i=1}^{N_b} n_{i\mathbf{k}}  \phi_{i\mathbf{k}}(\mathbf{r})  \phi_{i\mathbf{k}}^*(\mathbf{r}')
$$

where $n_{i\mathbf{k}} \in [0,1]$ are **natural occupation numbers**, $\phi_{i\mathbf{k}}$ are
**natural orbitals**, and $w_{\mathbf{k}}$ are k-point weights normalized so that
$\sum_{\mathbf{k}} w_{\mathbf{k}} = 1$ (or 2 accounting for spin in the restricted case).

### LCAO Expansion

Each natural orbital is expanded in localized atomic orbitals (LCAO):

$$
\phi_{i\mathbf{k}}(\mathbf{r}) = \sum_{\mu=1}^{N_{\mathrm{basis}}} C_{\mu i}^{\mathbf{k}}  \chi_{\mu \mathbf{k}}(\mathbf{r})
$$

The overlap matrix is $S_{\mu\nu}^{\mathbf{k}} = \langle \chi_{\mu\mathbf{k}} | \chi_{\nu\mathbf{k}} \rangle$.
Orthonormality of natural orbitals requires:

$$
(C^{\mathbf{k}})^\dagger S^{\mathbf{k}} C^{\mathbf{k}} = I_{N_b}
$$

Hence $C^{\mathbf{k}}$ lies on the **generalized Stiefel manifold**
$\mathrm{St}(N_b, N_{\mathrm{basis}}; S^{\mathbf{k}})$.

### Density Matrix in LCAO Basis

$$
\gamma_{\mu\nu}^{\mathbf{k}} = \sum_i n_{i\mathbf{k}}  C_{\mu i}^{\mathbf{k}} (C_{\nu i}^{\mathbf{k}})^*
= C^{\mathbf{k}}  \mathrm{diag}(\mathbf{n}_{\mathbf{k}})  (C^{\mathbf{k}})^\dagger
$$

The electron density is:

$$
\rho(\mathbf{r}) = \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_{\mu\nu} \gamma_{\mu\nu}^{\mathbf{k}}  \chi_{\mu\mathbf{k}}(\mathbf{r})  \chi_{\nu\mathbf{k}}^*(\mathbf{r})
$$

### Role of the planewave grid (ABACUS)

The **natural orbitals optimised in RDMFT are expanded in the LCAO basis** as in
the previous subsection. The electronic **density** is accumulated on the
real-space grid and the **Hartree** potential uses the same PW-based Hartree
machinery as the rest of LCAO (`H_Hartree_pw`, `rho_basis` passed into
`EnergyGradient::update_ion`). The helper `StiefelManifold` can be constructed
with `S = I` (identity overlap) as a *mathematical* special case (see comments
in `rdmft_stiefel.h`); the shipped RDMFT driver, however, always works in LCAO
with overlap $S^{\mathbf{k}}$ and the $S$-weighted Stiefel constraint
$(C^{\mathbf{k}})^\dagger S^{\mathbf{k}} C^{\mathbf{k}} = I$.

---

## 2. Total Energy Functional

$$
E[n_{i\mathbf{k}}, C^{\mathbf{k}}] = E_{\mathrm{one}} + E_H[\rho] + E_{xc}[\gamma]
+ E_{\mathrm{ent}}(n) + E_{\mathrm{Ewald}}
$$

where $E_{\mathrm{ent}} = \gamma \sum_{\mathbf{k},i} w_{\mathbf{k}} f_{\mathrm{bin}}(n_{i\mathbf{k}})$
is **optional HF-only** entropy (`rdmft_occ_entropy_gamma`, default $\gamma=0$).
The running total returned by `EnergyGradient::compute` is
`E_one_ + E_hartree_ + E_xc_ + E_entropy_ + E_ewald_`.

### 2.1 One-Body Energy

$$
E_{\mathrm{one}} = \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i n_{i\mathbf{k}}  h_{ii}^{\mathrm{one}}(\mathbf{k})
$$

where

$$
h_{ii}^{\mathrm{one}}(\mathbf{k}) = (C^{\mathbf{k}})^\dagger h^{\mathbf{k}} C^{\mathbf{k}} \big|_{ii}
$$

and $h^{\mathbf{k}} = T^{\mathbf{k}} + V_{\mathrm{local}}^{\mathbf{k}} + V_{\mathrm{NL}}^{\mathbf{k}}$
is the one-body Hamiltonian matrix (kinetic + local pseudopotential + nonlocal pseudopotential).

### 2.2 Hartree Energy

$$
E_H[\rho] = \frac{1}{2} \iint \frac{\rho(\mathbf{r})\rho(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|} d\mathbf{r}d\mathbf{r}'
$$

This is a functional of $\rho$, which in turn depends on $n_{i\mathbf{k}}$, $C^{\mathbf{k}}$.

### 2.3 Exchange-Correlation Functionals

In RDMFT, the xc energy is typically an explicit functional of the 1-RDM.
The code implements several kernels (`rdmft_functional` / `XCFunctionalType`),
many of which are defined through a coupling function $f(n_i, n_j)$ in the
exchange-like integral notation:

$$
E_{xc}[\gamma] = -\frac{1}{2} \sum_{\mathbf{k}\mathbf{k}'} w_{\mathbf{k}} w_{\mathbf{k}'}
\sum_{ij} f(n_{i\mathbf{k}}, n_{j\mathbf{k}'})  K_{ij}^{\mathbf{k}\mathbf{k}'}
$$

where $K_{ij}^{\mathbf{k}\mathbf{k}'} = \langle \phi_{i\mathbf{k}} \phi_{j\mathbf{k}'} | \hat{v}_{c} | \phi_{j\mathbf{k}'} \phi_{i\mathbf{k}} \rangle$
are two-electron exchange integrals.

The sign and $\tfrac{1}{2}$ prefactor in this **continuum** notation are conventional;
the **code** evaluates exchange through the modified density matrix and
`Exx_LRI` as in §2.4 (no explicit dense $K_{ij}$ tensor).

#### 2.3.1 Hartree-Fock (HF)

$$
f^{\mathrm{HF}}(n_i, n_j) = n_i  n_j, \qquad g^{\mathrm{HF}}(n) = n
$$

#### 2.3.2 Müller Functional

$$
f^{\mathrm{M}}(n_i, n_j) = \sqrt{n_i} \sqrt{n_j}, \qquad g^{\mathrm{M}}(n) = \sqrt{n} = n^{1/2}
$$

#### 2.3.3 Power Functional

$$
f^{\mathrm{P}}(n_i, n_j) = n_i^\alpha  n_j^\alpha, \qquad g^{\mathrm{P}}(n) = n^\alpha
$$

with $\alpha \in (0.5, 1)$. Standard choices: $\alpha \approx 0.656$ for solids,
$\alpha \approx 0.525$ for molecular dissociation.

#### 2.3.4 Goedecker-Umrigar (GU) Functional

$$
f^{\mathrm{GU}}(n_i, n_j) = \begin{cases}
\sqrt{n_i}\sqrt{n_j} & i \neq j \\
n_i^2 & i = j
\end{cases}
$$

The diagonal modification ensures correct self-interaction cancellation.
The off-diagonal part uses the Müller form.

#### 2.3.5 Corrected Hartree-Fock (CHF, INPUT `chf`)

$$
f^{\mathrm{CHF}}(n_i, n_j)
= \tfrac{1}{2} n_i n_j
+ \tfrac{1}{2} \sqrt{n_i(1-n_i)} \, \sqrt{n_j(1-n_j)}
$$

(as implemented in `XCFunctional::f`; derivatives use `chf_corr_term_deriv`).

#### 2.3.6 Csanyi–Goedecker–Arias (CGA, INPUT `cga`)

$$
f^{\mathrm{CGA}}(n_i, n_j)
= \tfrac{1}{4} n_i n_j
+ \tfrac{1}{4} \sqrt{n_i(2-n_i)} \, \sqrt{n_j(2-n_j)}
$$

#### 2.3.7 GEO functional (INPUT `geo`)

$$
f^{\mathrm{GEO}}(n_i, n_j)
= \frac{1}{4} n_i n_j
+ \frac{1}{4} n_i^{1/2} n_j^{1/2}
+ \frac{1}{2} n_i^{3/4} n_j^{3/4}
$$

implemented as a **sum of three separable power pieces** with coefficients
$(\tfrac{1}{4}, \tfrac{1}{4}, \tfrac{1}{2})$ and exponents $(1, \tfrac{1}{2}, \tfrac{3}{4})$,
each regularised at small $n$ like the Power functional (`pow_reg` / `dpow_reg`).

#### 2.3.8 HybOpt (INPUT `hybopt`)

Convex combination of HF and a fixed Power-like channel:

$$
f^{\mathrm{HybOpt}}(n_i, n_j)
= w_{\mathrm{HF}} \, n_i n_j
+ w_{\mathrm{P}} \, n_i^{\alpha_{\mathrm{HybOpt}}} n_j^{\alpha_{\mathrm{HybOpt}}},
\qquad
\alpha_{\mathrm{HybOpt}} = 0.541076,\quad
w_{\mathrm{P}} = 0.938328
$$

(with $w_{\mathrm{HF}} = 1 - w_{\mathrm{P}}$).

### 2.4 Exchange energy in the ABACUS implementation

#### Modified density matrix (single EXX build)

In `EnergyGradient::build_DM_xc`, the exchange density matrix at each
$\mathbf{k}$ is assembled with **band weights** $w_{\mathbf{k}}\,g(n_{i\mathbf{k}})$
fed into the same `cal_dm_psi` machinery as the KS density matrix. In LCAO
components,

$$
\gamma_{\mathrm{xc},\mu\nu}^{\mathbf{k}}
= \sum_i w_{\mathbf{k}} \, g(n_{i\mathbf{k}}) \,
C_{\mu i}^{\mathbf{k}} (C_{\nu i}^{\mathbf{k}})^*
$$

Bands with $|g(n)|$ below a tiny cutoff are skipped (`rdmft_occ_weight_eps`).

The **Fock exchange** Hamiltonian $H_{\mathrm{exx}}$ is then built from
$\gamma_{\mathrm{xc}}$ via the **range-separated / density-fitting EXX stack**
(`Exx_LRI`, `RI_2D_Comm::add_Hexx`, `cal_exx_elec`), not a literal dense
four-index $K_{ij}$ build.

For **HF** in code, $g(n)=\max(0,n)$ (occupations are still clipped to $[0,1]$
for other functionals; HF keeps the raw weight for the KS seed convention).

#### Energy accumulation (separable HF / Müller / Power / GU)

Let $\varepsilon^{\mathrm{exx}}_{i\mathbf{k}} = \langle \phi_{i\mathbf{k}} |
H_{\mathrm{exx}}[\gamma_{\mathrm{xc}}] | \phi_{i\mathbf{k}} \rangle$ be the
diagonal returned as `vx_diag` in `EnergyGradient::compute`. The exchange part
of the RDMFT energy is accumulated **per band** as

$$
E_{xc} = \frac{1}{2} \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i
g(n_{i\mathbf{k}}) \, \varepsilon^{\mathrm{exx}}_{i\mathbf{k}}
$$

(the same $\tfrac{1}{2}$ prefactor as the Hartree term uses for the explicit
$n$-weighted diagonal).

**GU functional (implementation detail).** The analytic GU kernel uses
$f^{\mathrm{GU}}(n_i,n_j)=\sqrt{n_i}\sqrt{n_j}$ off-diagonal and $n_i^2$ on the
diagonal; `XCFunctional` also exposes `gu_diag_factor(n)=n^2-n$ and its
derivative for that diagonal correction. In the **current** `EnergyGradient`
path, GU is assigned the same Müller exponent $\alpha=\tfrac{1}{2}$ for $g$ and
$\mathrm{d}g$ as Müller, and **only** the single modified DM above is passed to
EXX: the explicit extra energy $\tfrac{1}{2}\sum_{\mathbf{k}} w_{\mathbf{k}}^2
\sum_i (n_{i\mathbf{k}}^2-n_{i\mathbf{k}})\,J_{ii}^{\mathbf{k}\mathbf{k}}$ is
**not** added in `compute` / `compute_energy`. Thus GU in production runs
coincides with the Müller-type EXX treatment until a future patch wires in the
diagonal correction.

#### Mixed channels (GEO, CHF, CGA, HybOpt)

For these types, `EnergyGradient::compute` performs **several** EXX builds in
sequence: for each channel $t$ it forms a modified DM with per-band weights
$w_{\mathbf{k}}$ times the channel weight $p_t(n_{i\mathbf{k}})$ (linear $n$ for
the HF-like piece, $\sqrt{n(1-n)}$, $\sqrt{n(2-n)}$, or regularised
$n^{\alpha_t}$ as coded), accumulates $H_{\mathrm{exx}}\psi$ and diagonal
contributions, and sums energy and gradients with the fixed coefficients
(`geo_coef` / `hybopt_*_weight` / CHF / CGA weights).

---

## 3. Constraints

### 3.1 Electron Number Conservation

$$
\sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i n_{i\mathbf{k}} = N_e
$$

where $N_e$ is the number of electrons per unit cell.

The solver target $N_e$ may follow KS-weighted sums or `PARAM.inp.nelec` plus an
optional `rdmft_nelec_delta` offset (`RDMFTNelectronTargetMeta` in `rdmft_type.h`).

### 3.2 Ensemble N-representability

$$
0 \leq n_{i\mathbf{k}} \leq 1 \quad \forall i, \mathbf{k}
$$

### 3.3 Orthonormality of Natural Orbitals

$$
(C^{\mathbf{k}})^\dagger S^{\mathbf{k}} C^{\mathbf{k}} = I \quad \forall \mathbf{k}
$$

---

## 4. Occupation Number Parameterizations

To handle the box constraint $n \in [0,1]$ and electron number constraint via
unconstrained optimization, we introduce parameterizations.

### 4.1 Cosine-Squared Parameterization

$$
n_{i\mathbf{k}} = \cos^2(\theta_{i\mathbf{k}})
$$

**Jacobian:**

$$
\frac{\partial n}{\partial \theta} = -2\cos\theta\sin\theta = -\sin(2\theta)
$$

**Inverse:** $\theta = \arccos(\sqrt{n})$ 

### 4.2 Logistic (Sigmoid) Parameterization

$$
n_{i\mathbf{k}} = \sigma(x_{i\mathbf{k}}) = \frac{1}{1 + e^{-x_{i\mathbf{k}}}}
$$

**Jacobian:**

$$
\frac{\partial n}{\partial x} = n(1-n)
$$

**Inverse:** $x = \ln(n/(1-n))$ 

### 4.3 Sigma-shift parameterization (INPUT `sigma_shift`)

Used with the **joint** solver only (`RDMFTSolver::init` rejects `sigma_shift` when
`rdmft_solver_strategy` is alternating). Unconstrained parameters
$z_{i\mathbf{k}}\in\mathbb{R}$ are mapped through a **shared** scalar shift
$\lambda$ and the stable sigmoid $\sigma$ (`SigmaShiftOccParam::stable_sigmoid`):

$$
n_{i\mathbf{k}} = \sigma(z_{i\mathbf{k}} + \lambda),
\qquad
\sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i n_{i\mathbf{k}} = N_e .
$$

Each outer step solves for $\lambda$ by **bisection** on the monotone map
$\lambda \mapsto \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i \sigma(z_{i\mathbf{k}}+\lambda)$.
The chain rule for $\partial E/\partial z$ subtracts a single electron-count
term proportional to $w_{\mathbf{k}}$; see `SigmaShiftOccParam` in
`rdmft_occupation.h`.

### 4.4 Transformed gradient

Given the energy gradient $\partial E / \partial n_{i\mathbf{k}}$, the gradient in the
unconstrained parameter is:

$$
\frac{\partial E}{\partial \theta_{i\mathbf{k}}} = \frac{\partial E}{\partial n_{i\mathbf{k}}} \cdot \frac{\partial n}{\partial \theta} = -\sin(2\theta_{i\mathbf{k}})  \frac{\partial E}{\partial n_{i\mathbf{k}}}
$$

$$
\frac{\partial E}{\partial x_{i\mathbf{k}}} = n_{i\mathbf{k}}(1 - n_{i\mathbf{k}})  \frac{\partial E}{\partial n_{i\mathbf{k}}}
$$

---

## 5. Gradients of the Total Energy

### 5.1 Gradient w.r.t. natural occupation numbers (as in `EnergyGradient::compute`)

Let $h_{ii}^{\mathrm{one}}(\mathbf{k})$ and $v_{H,ii}(\mathbf{k})$ denote the
**diagonal** matrix elements (code: `h_one_diag`, `vh_diag`) of the one-body
Hamiltonian and the Hartree potential in the natural-orbital basis at
$(\mathbf{k},i)$, and let $\varepsilon^{\mathrm{exx}}_{i\mathbf{k}}$ be the
corresponding diagonal of $H_{\mathrm{exx}}[\gamma_{\mathrm{xc}}]$ (`vx_diag`).

The **accumulated** partial derivatives use a **plus** Hartree diagonal (this is
what the implementation differentiates consistently with the $n$-weighted
Hartree energy term and the internal potential rebuild):

$$
\frac{\partial E}{\partial n_{i\mathbf{k}}}
= w_{\mathbf{k}} \Bigl(
h_{ii}^{\mathrm{one}}(\mathbf{k}) + v_{H,ii}(\mathbf{k})
\Bigr)
+ \frac{\partial E_{xc}}{\partial n_{i\mathbf{k}}}
$$

For **separable** HF / Müller / Power / GU (single EXX build),

$$
\frac{\partial E_{xc}}{\partial n_{i\mathbf{k}}}
= w_{\mathbf{k}} \, g'(n_{i\mathbf{k}}) \,
\varepsilon^{\mathrm{exx}}_{i\mathbf{k}} .
$$

For **mixed** GEO / CHF / CGA / HybOpt, the code adds the weighted derivatives of
each channel’s occupation weights to the same diagonal exchange response
(`mix_vx_G_acc` in `rdmft_energy_gradient.cpp`).

**Optional HF entropy** (INPUT `rdmft_occ_entropy_gamma` $=\gamma>0$, HF only):

$$
E_{\mathrm{ent}} = \gamma \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i f_{\mathrm{bin}}(n_{i\mathbf{k}}),
\qquad
\frac{\partial E_{\mathrm{ent}}}{\partial n_{i\mathbf{k}}}
= \gamma \, w_{\mathbf{k}} \, f'_{\mathrm{bin}}(n_{i\mathbf{k}}),
$$

with $f_{\mathrm{bin}}(n)=n\ln n+(1-n)\ln(1-n)$ (`binary_entropy_f`).

**GU remark.** The analytic extra term from the diagonal GU kernel is **not**
added on top of the Müller-type $\mathrm{d}g$ contribution in `compute`; the
occupation derivative follows the same formula as Müller/Power with GU’s
regularised $g,\mathrm{d}g$ at $\alpha=\tfrac{1}{2}$.

### Derivatives of $g(n)$ (closed form on $n\ge\varepsilon$; see code for $n<\varepsilon$)


| Functional | $g(n)$ (code, $n\in[0,1]$) | $g'(n)$ on $n\ge\varepsilon$ |
| ---------- | ------------------------- | ------------------------------ |
| HF         | $\max(0,n)$               | $1$                            |
| Müller     | $n^{1/2}$ (regularised)   | $\tfrac{1}{2} n^{-1/2}$        |
| Power      | $n^{\alpha}$ (regularised) | $\alpha n^{\alpha-1}$        |
| GU         | same Müller branch        | same Müller branch             |

For $\alpha<1$, `XCFunctional` replaces $g$ and $g'$ below a cutoff
$\varepsilon_{\mathrm{reg}}$ (default $10^{-8}$) by a linear Taylor extrapolation
from $n=\varepsilon_{\mathrm{reg}}$ so $g'(0)$ stays finite (`rdmft_xc_functional.h`).


### 5.2 Gradient w.r.t. orbital coefficients (C-space, before Stiefel projection)

The **returned** LCAO gradient in `EnergyGradient::compute` includes the usual
factor **2** from the real pairing / Wirtinger convention documented in the
source (gamma-only real case matches $\partial E/\partial C = 2 n H C$ for a
single quadratic orbital energy):

$$
G_{\mu i}^{\mathbf{k}}
= 2 \, w_{\mathbf{k}} \Bigl[
n_{i\mathbf{k}} \bigl((h^{\mathbf{k}} + V_H^{\mathbf{k}}) C^{\mathbf{k}}\bigr)_{\mu i}
+ \bigl(\text{exchange column}\bigr)_{\mu i}
\Bigr],
$$

where the exchange column is $g(n_{i\mathbf{k}})\,(H_{\mathrm{exx}}^{\mathbf{k}} C^{\mathbf{k}})_{\mu i}$
for a single EXX build, and for mixed functionals it is the pre-summed
$H_{\mathrm{exx}}\psi$ contribution per band (`mix_Hpsi_x_acc`). Bands with
zero occupation (or zero $g(n)$ in the exchange-only column) are skipped.

After assembly, $G$ is converted to **X-space** (`grad_C_to_X`) and projected on
the Stiefel tangent space (`project_orbital_gradient`).


### 5.3 Riemannian gradient on the Stiefel manifold (code: `StiefelManifold::project_tangent`)

With overlap $S^{\mathbf{k}}$, the tangent space at $C$ is
$T_C \mathrm{St} = \{ Z : C^\dagger S Z + Z^\dagger S C = 0 \}$.
The implementation projects the ambient (Euclidean) gradient $G$ as

$$
\mathrm{proj}_C(G) = G - C \,\mathrm{sym}\!\bigl(C^\dagger S G\bigr),
\qquad
\mathrm{sym}(A)=\tfrac{1}{2}(A+A^\dagger),
$$

i.e. the first Hermitian factor uses $S C$ when $S$ is present (`rdmft_stiefel.h`).

The **inner product** on tangent vectors is
$\langle \eta_1, \eta_2 \rangle = \mathrm{Re}\,\mathrm{Tr}(\eta_1^\dagger S \eta_2)$.

### 5.4 Retraction (code: `StiefelManifold::retract` / `reorthogonalize`)

Given a tangent step $\eta$ and line-search parameter $\alpha$, set
$Y = C + \alpha\eta$, form the Gram matrix $M = Y^\dagger S Y$, Cholesky factor
$M = L L^\dagger$, and retract to

$$
C_{\mathrm{new}} = Y \, L^{-\dagger}.
$$

(Comments in the header also mention a polar form; the **default** retraction
path used by the solver is this Cholesky-based $S$-orthogonalisation.)

---

## 6. Constraint Handling for Occupations

### 6.1 Augmented Lagrangian Method

The constrained problem

$$
\min_{\mathbf{n}} E(\mathbf{n}) \quad \text{s.t.} \quad c(\mathbf{n}) = \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i n_{i\mathbf{k}} - N_e = 0, \quad 0 \leq n_{i\mathbf{k}} \leq 1
$$

is solved via the augmented Lagrangian:

$$
\mathcal{L}_A(\mathbf{n}, \lambda, \mu) = E(\mathbf{n}) + \lambda  c(\mathbf{n}) + \frac{\mu}{2} c(\mathbf{n})^2
$$

With the parameterization $n = \cos^2\theta$ (or logistic), the box constraints are
automatically satisfied, and we iteratively:

1. Minimize $\mathcal{L}_A$ w.r.t. $\theta$ (unconstrained)
2. Update $\lambda \leftarrow \lambda + \mu  c(\mathbf{n})$
3. Optionally increase $\mu$

The gradient of $\mathcal{L}_A$ w.r.t. $n_{i\mathbf{k}}$ is:

$$
\frac{\partial \mathcal{L}_A}{\partial n_{i\mathbf{k}}} = \frac{\partial E}{\partial n_{i\mathbf{k}}} + (\lambda + \mu  c(\mathbf{n}))  w_{\mathbf{k}}
$$

### 6.2 Projected Gradient Method

Optimize $n_{i\mathbf{k}}$ on the feasible set using a **projected search**
along a descent direction, with monotone Armijo backtracking.

**Classical PG (reference).**  The textbook iterate is a gradient step followed
by projection:
$$
\mathbf{x} = \mathbf{n}^{(t)} - \alpha_t \, \nabla_{\mathbf{n}} E(\mathbf{n}^{(t)}),\qquad
\mathbf{n}^{(t+1)} = P_{\mathcal{C}}(\mathbf{x}),
$$
i.e. the map $\mathbf{n} \mapsto P_{\mathcal{C}}(\mathbf{n} - \alpha \nabla E)$.

**ABACUS implementation (`rdmft_constraint = projected_gradient`).**  Each inner
iteration uses a search direction $\mathbf{d}$ in occupation space from a
configurable Euclidean optimizer (steepest descent, nonlinear conjugate
gradient, L-BFGS, or Adam), controlled by INPUT `rdmft_occ_optimizer`.  If
$\mathbf{d}^\top \mathbf{g} \ge 0$ with $\mathbf{g} = \nabla_{\mathbf{n}} E$,
the code **replaces** $\mathbf{d}$ by $-\mathbf{g}$ (descent safeguard).  The
trial is the **curvilinear** projected point
$$
\mathbf{n}'(\alpha) = P_{\mathcal{C}}\bigl(\mathbf{n} + \alpha \mathbf{d}\bigr),
$$
not only $P_{\mathcal{C}}(\mathbf{n} - \alpha \mathbf{g})$.  Classical PG is
the special case $\mathbf{d} = -\mathbf{g}$.

The projector $P_{\mathcal{C}}$ is the same for all variants; it maps onto
$$\mathcal{C} = \left\{\mathbf{n}:\; 0\le n_{i\mathbf{k}}\le 1,\; \sum_{\mathbf{k}} w_{\mathbf{k}}\sum_i n_{i\mathbf{k}} = N_e\right\}.$$

**Closed form of $P_{\mathcal{C}}$.**  The projection is the unique minimizer
$$
P_{\mathcal{C}}(\mathbf{x})
= \arg\min_{\mathbf{y}} \; \tfrac{1}{2} \lVert \mathbf{y} - \mathbf{x} \rVert_2^2
\quad \text{s.t.} \quad
\sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i y_{i\mathbf{k}} = N_e,\; 0 \le y_{i\mathbf{k}} \le 1.
$$
The equality is linear and the normal direction is the same for all bands at a
given $\mathbf{k}$: $(\nabla_{\mathbf{n}} c)_i = w_{\mathbf{k}}$.  Writing KKT
conditions (one scalar Lagrange multiplier $\lambda$ for the equality) yields the
**shifted clipping** form
$$
y_{i\mathbf{k}}(\lambda)
= \min\bigl(1, \max(0, x_{i\mathbf{k}} - \lambda  w_{\mathbf{k}}) \bigr),
$$
and the scalar $\lambda$ is fixed by the scalar equation
$$
\sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i y_{i\mathbf{k}}(\lambda) = N_e.
$$
The left-hand side is a **non-increasing, piecewise-linear** function of
$\lambda$, so a bracket followed by **bisection** (or a monotone 1D root find)
gives $\lambda$ efficiently.  (This is the correct PG projection; it is **not**
the same as “clip then multiplicatively rescale” free components: that rescaling
is a different feasible repair, not the Euclidean foot point.)  The feasible
range of the weighted sum is
$\bigl[0,\; \sum_{\mathbf{k}} w_{\mathbf{k}} N_b\bigr]$
with $N_b$ the number of bands per $\mathbf{k}$; if $N_e$ lies outside, the
constraint set is empty and no exact projection exists.

**Optional approximations (not the Euclidean projector):**

- **Box only:**  $y_{i\mathbf{k}} = \min(1,\max(0, x_{i\mathbf{k}}))$ — does not
  enforce the electron sum.
- **Sort-and-shift on a vector:**  standard $O(m\log m)$ simplex algorithms apply
  to $\sum_i z_i = \text{const}$ with $z_i \ge 0$ (one multiplier per
  *component*).  Here the constraint ties **all** $(i,\mathbf{k})$ through
  $w_{\mathbf{k}}$ in the *same* way for every band at $\mathbf{k}$, so the
  correct projector is the single-$\lambda$ clipping form above, not a
  per-dimension water-filling sort unless the problem is reparameterized
  accordingly.

**Step size and line search (ABACUS).**  The first trial step $\alpha_0$ for
backtracking is chosen by INPUT `rdmft_occ_ls_init_step`: fixed (`rdmft_occ_ls_stepsize`), Barzilai–Borwein
(uses `rdmft_alm_bb_mode` and clamps `rdmft_alm_bb_alpha_min` /
`rdmft_alm_bb_alpha_max`, with fallback to the fixed default step `1`; unlike the ALM
path this BB seed is **not** gated on `rdmft_alm_bb_enabled`), or a quadratic
model from the previous inner iteration’s first energy trial (see
`rdmft_usage.md`).  Monotone **Armijo**
accepts $\alpha$ when
$$
E\bigl(\mathbf{n}'(\alpha)\bigr) \le E(\mathbf{n}) + c_1\, \mathbf{g}^\top \bigl(\mathbf{n}'(\alpha) - \mathbf{n}\bigr),
$$
with $c_1$ from INPUT `rdmft_line_search_c1` and geometric backtracking (multiply $\alpha$ by
`line_search_rho`, default $0.5$ in `RDMFTConfig`; PG uses pure geometric shrink, not the polynomial
Armijo refinements controlled by `rdmft_line_search_polynomial`).
Trials with $\mathbf{g}^\top(\mathbf{n}'(\alpha)-\mathbf{n}) \ge 0$ are rejected
(step not first-order descent along the projected segment).  If after
projection the weighted electron sum differs from $N_e$ by more than $10^{-6}$,
the trial is rejected (pathological / near-infeasible bracket).  **Special case
(steepest descent):** $\mathbf{d}=-\mathbf{g}$ gives
the familiar map $\phi(\alpha)=E(P_{\mathcal{C}}(\mathbf{n}-\alpha\mathbf{g}))$
with the same Armijo inequality in terms of $\mathbf{n}'(\alpha)-\mathbf{n}$.

**Line search failure recovery.**  If no Armijo step is found, the solver applies
one **projected steepest** move $\mathbf{n} \leftarrow P_{\mathcal{C}}(\mathbf{n} - \tau_{\mathrm{ls}}\,\mathbf{g})$
with $\tau_{\mathrm{ls}} = 1$, resets the occupation optimizer
curvature state (CG / L-BFGS history), and continues.

**Stopping criteria (ABACUS inner loop).**  Let $\mathbf{g}_{\mathrm{proj}} =
\bigl(\mathbf{n} - P_{\mathcal{C}}(\mathbf{n} - \tau \mathbf{g})\bigr)/\tau$
(Bertsekas projected-gradient vector).  **$\tau$ is taken from the occupation
line search** so the stationarity measure uses the same scale as the step:
**pre-step**, $\tau = \alpha_0$ (the first Armijo trial from
`rdmft_occ_ls_init_step`: fixed (`rdmft_occ_ls_stepsize`), Barzilai–Borwein, or quad, with fallback to
$1$ if the estimate is non-positive or non-finite); **post-step**
after a successful line search, $\tau = \alpha_{\mathrm{acc}}$ (accepted
Armijo step along $\mathbf{d}$); after **SD fallback** (failed line search),
$\tau = 1$ (same step length as the recovery move).  The inner
loop stops when $\|\mathbf{g}_{\mathrm{proj}}\|_\infty <$ `rdmft_occ_grad_tol`,
evaluated **after** an accepted or fallback step (post-step $\mathbf{n}$ and
gradient).  On the first inner iteration, if the **pre-step**
$\|\mathbf{g}_{\mathrm{proj}}\|_\infty$ is already below tolerance, the inner
loop exits immediately without a line search.
**`rdmft_occ_energy_tol` is not used** for the PG occupation inner (kept for INPUT
compatibility).  The iteration cap is `rdmft_occ_maxiter`.

Numerical tips and implementation notes:

- Use analytic gradients $\partial E/\partial n$; finite differences are costly.
- **INPUT keywords (PG, ABACUS):** `rdmft_occ_optimizer`, `rdmft_occ_ls_init_step`
  (sets Bertsekas pre-step $\tau=\alpha_0$ together with BB/quad clamps),
  fixed $\alpha=1$ (fallback for invalid $\alpha_0$ and SD-fallback step length),
  `rdmft_occ_grad_tol`, `rdmft_occ_maxiter`, `rdmft_line_search_c1`, and for the
  `bb` branch of `rdmft_occ_ls_init_step` also `rdmft_alm_bb_mode`,
  `rdmft_alm_bb_alpha_min`, `rdmft_alm_bb_alpha_max`.  Defaults and semantics
  are tabulated in [`rdmft_usage.md`](rdmft_usage.md).
- In practice, $P_{\mathcal{C}}$ is evaluated by solving the scalar dual
  $\lambda$ in the shift–clip form above (bisection on a monotone
  piecewise-linear map).  PG and active-set occupation updates in the code apply
  this projection to trial occupation vectors after each line-search step.
- If many occupations are interior (not at bounds), projected gradient is efficient.
- If orthonormal orbitals are optimized concurrently, use separate step sizes for
   orbitals and occupations or alternate updates (block coordinate style).
- For orbital constraints (Stiefel), use the Cholesky-based $S$-orthogonalisation
  retraction in production (`StiefelManifold::retract`); see §5.4.

Pseudocode (ABACUS projected-gradient inner loop; `project` denotes
$P_{\mathcal{C}}$; neglect spin labels):

```text
initialize n feasible; configure occ_optimizer, tol = rdmft_occ_grad_tol
for inner = 0 .. occ_maxiter-1:
   g = grad_n(E, n)
   alpha0 = initial_step_from(rdmft_occ_ls_init_step, BB/quad, rdmft_occ_ls_stepsize, clamps)
   tau_pre = alpha0 if alpha0 > 0 else rdmft_occ_ls_stepsize
   g_proj_pre = (n - project(n - tau_pre * g)) / tau_pre
   if inner == 0 and ||g_proj_pre||_inf < tol: break     // early exit
   d = direction_from_optimizer(g, history)               // SD / CG / lbfgs / adam
   if dot(d, g) >= 0: d = -g                              // descent safeguard
   alpha = alpha0
   n_save = n
   success = false
   alpha_acc = 0
   for trial = 1 .. line_search_max_iter:
      n_try = project(n + alpha * d)
      if |sum(w*n_try) - N_e| > 1e-6: alpha *= rho; continue
      delta = n_try - n
      if dot(g, delta) >= 0: alpha *= rho; continue
      if E(n_try) <= E(n) + c1 * dot(g, delta): success = true; n = n_try; alpha_acc = alpha; break
      alpha *= rho
   if not success:
      n = project(n_save - 1 * g)         // SD fallback
      alpha_acc = 1
      reset_optimizer_curvature_state()
   recompute g at new n
   tau_post = alpha_acc                                   // accepted α, or 1 after fallback
   g_proj_post = (n - project(n - tau_post * g)) / tau_post
   if ||g_proj_post||_inf < tol: break
end
```

When not to use:

- If the active-set (set of variables at bounds) is small and changes infrequently,
   an active-set method or a second-order reduced Newton solve on the free set
   may converge much faster near the solution.

### 6.3 Active set (`rdmft_constraint = active_set`)

This path is **not** a full reduced-space Newton KKT solve on the free variables;
it is a **gradient-projection style** loop that shares much machinery with the PG
case (`rdmft_solver.cpp`, `ConstraintMethod::ActiveSet`).

**Active identification** (`OccupationConstraint::identify_active_set`): for each
flattened index $I=(i,\mathbf{k})$ with tolerance `tol` (default $10^{-8}$),

- **lower active** if $n_I \le \texttt{tol}$ and $(\nabla E)_I > 0$ (gradient points
  into the interior of $[0,1]$ from the $n=0$ face);
- **upper active** if $n_I \ge 1-\texttt{tol}$ and $(\nabla E)_I < 0$;
- otherwise **free**.

**Reduced gradient** (`apply_active_set`): components on active indices are set to
zero; on free indices the code adds $\lambda\, w_{\mathbf{k}}$ with
$\lambda = -(\sum_{I\in\mathcal{F}} w_{\mathbf{k}(I)} (\nabla E)_I) /
\sum_{I\in\mathcal{F}} w_{\mathbf{k}(I)}^2$ so the modified vector is orthogonal
(in the $w$-weighted sense) to the electron-sum constraint row.

**Search direction** uses the same `EuclideanOptimizer` family as PG (SD / CG /
L-BFGS / Adam) on this modified gradient; the direction is then **zeroed on active
indices** so bound-pinned occupations do not move.

**Line search** follows the same `rdmft_occ_ls_type` policy as PG (Armijo by
default; optional strong / weak Wolfe for CG / L-BFGS in alternating mode, etc.),
with **trial points clipped and re-projected** onto
$\mathcal{C}=\{0\le n\le 1,\ \sum w n = N_e\}$ via the same shift–clip dual solve
as in §6.2.

**Stopping / diagnostics**: early exit on the first inner iteration if the
constraint residual, free-set stationarity norm of the modified gradient, and a
**dual complementarity** surrogate (`active_set_dual_complementarity_violation` in
`rdmft_solver.cpp`) are all below `rdmft_occ_tol`. Otherwise the loop continues up
to `rdmft_occ_maxiter`. When the **number of active constraints changes**, the
occupation optimiser state is **reset** (same pattern as PG line-search failure).

For a textbook reduced active-set Newton method on the free face, see standard
optimisation references; the paragraph above describes what is actually coded.


---

## 7. Optimization Strategies

### 7.1 Alternating Optimization

Alternate between:

- **Orbital step**: Fix $n_{i\mathbf{k}}$, optimize $C^{\mathbf{k}}$ on Stiefel manifold
- **Occupation step**: Fix $C^{\mathbf{k}}$, optimize $n_{i\mathbf{k}}$ with constraints

### 7.2 Joint (Product Manifold) Optimization

Define the product manifold:

$$
\mathcal{M} = \underbrace{\mathbb{R}^{N_k \times N_b}}_{\text{occupation (Euclidean)}} \times \prod_{\mathbf{k}} \underbrace{\mathrm{St}(N_b, N; S^{\mathbf{k}})}_{\text{orbitals (Stiefel)}}
$$

With cosine-squared or logistic parameters, the occupation block is unconstrained
Euclidean **together with** the augmented Lagrangian for $\sum w n = N_e$ when
using `rdmft_constraint = augmented_lagrangian`. With **`rdmft_occ_param sigma_shift`**
and the **joint** strategy, the electron sum is enforced by the implicit shift
(§4.3) and **no** augmented-Lagrangian penalty on $N_e$ is used for occupations.

The tangent vector at a point $(\boldsymbol{p}, C^{\mathbf{k}})$ is
$(\delta\boldsymbol{p}, \eta^{\mathbf{k}})$ where
$\delta\boldsymbol{p} \in \mathbb{R}^{N_k \times N_b}$ and
$\eta^{\mathbf{k}} \in T_{C^{\mathbf{k}}} \mathrm{St}$.

Retraction: update occupation parameters in Euclidean fashion (linear in $p$ for
cosine/logistic; implicit $\lambda$ solve each step for sigma-shift), and apply
the Stiefel retraction of §5.4 to each $C^{\mathbf{k}}$.

In the implementation (`rdmft_solver.cpp::solve_joint`), the occupation
parameters and orbital coefficients are packed into **one** vector
$z = (\,p\,,\,\mathrm{flat}(C^1), \ldots, \mathrm{flat}(C^{N_k})\,)$ and a
**single** Euclidean optimiser, selected via the `rdmft_joint_optimizer`
keyword (default lbfgs), drives its evolution. Each outer iteration:

1. evaluates the energy $E$ and Euclidean gradients $(\nabla_n E, \nabla_C E)$;
2. chain-rules $\nabla_n E \to \nabla_p E$ through the occupation
   parameterisation, and projects $\nabla_C E$ onto the Stiefel tangent space
   at $C^{\mathbf{k}}$ to obtain $G_R^{\mathbf{k}}$;
3. packs the gradient as $g = (\nabla_p E,\,\mathrm{flat}(G_R^1), \ldots)$
   and asks the single unified optimiser for a packed descent direction
   $d = (d_p, \mathrm{flat}(d_{C^1}), \ldots)$;
4. re-projects each orbital block $d_{C^{\mathbf{k}}}$ onto the tangent space
   at the current $C^{\mathbf{k}}$ (necessary because Euclidean preconditioners
   used by lbfgs / Adam generally leave the tangent space) and falls back to
   the packed steepest-descent direction $d = -g$ if the joint directional
   derivative
   $dd_{\text{total}} = \langle \nabla_p E, d_p\rangle + \sum_{\mathbf{k}} \langle G_R^{\mathbf{k}}, d_{C^{\mathbf{k}}}\rangle_{S^{\mathbf{k}}}$
   is not negative;
5. performs **geometric** Armijo backtracking along the packed direction (trial
   $\alpha$, shrink by `line_search_rho` until Armijo holds or `line_search_max_iter`
   is hit—**no** polynomial interpolation branch, unlike `armijo_line_search` in
   `rdmft_optimizer.h`): occupation parameters update linearly
   $p \leftarrow p + \alpha\, d_p$, while each $C^{\mathbf{k}}$ is retracted onto
   the generalised Stiefel manifold via
   $C^{\mathbf{k}} \leftarrow R_{C^{\mathbf{k}}}(\alpha\, d_{C^{\mathbf{k}}})$;
6. builds the new packed gradient at the step's end-point and feeds
   $(g_\text{new}, \alpha\, d)$ to the unified optimiser's `update()`, so its
   history (lbfgs $(s,y)$ pairs, Adam moments, CG previous gradient) is
   updated once with a coherent product-manifold view; then refreshes the
   augmented-Lagrangian multiplier.

This is precisely Riemannian optimisation on the product manifold
$\mathcal{M}$ with a single shared optimiser and line search, rather than
two independent block updates.

---

## 8. Optimisation algorithms and line searches

The **Euclidean** directions for occupations (after parameterisation) and for the
**packed** joint vector are produced by `EuclideanOptimizer` in `rdmft_optimizer.h`.
Orbital-only steps in the **alternating** strategy additionally **project** the
ambient gradient to the Stiefel tangent space and use **Riemannian** inner
products in the line search (`stiefel_canonical_inner_product`).

### 8.1 Steepest descent (SD)

Search direction $d_k = -g_k$ (negative Euclidean gradient in the current
working variables: unconstrained occupation parameters, packed joint vector, or
ambient orbital gradient before projection—depending on the loop).

### 8.2 Nonlinear conjugate gradient (CG)

Implemented as **Polak–Ribière** on the current gradient $g_k$ and the previous
gradient $g_{k-1}$ (`compute_cg_direction`):

$$
\beta_k^{\mathrm{PR}}
= \frac{g_k^\top (g_k - g_{k-1})}{\|g_{k-1}\|^2},
\qquad
d_k = -g_k + \beta_k^{\mathrm{PR}} \, d_{k-1}.
$$

Details in code:

- **First step** ($k=0$): $d_0 = -g_0$.
- **Restart:** $\beta_k^{\mathrm{PR}} \leftarrow \max(0, \beta_k^{\mathrm{PR}})$ (Hager–style
  safeguard so PR never reverses the search direction).
- **Descent test:** if $d_k^\top g_k \ge 0$, the code **replaces** $d_k$ by $-g_k$
  (same safeguard as the projected occupation loop).

**Alternating orbital CG** (`optimize_orbitals` in `rdmft_solver.cpp`): after forming
$d_k$ from projected Riemannian gradients, a linear combination
$d_k \leftarrow -g_k + \beta_k \,\mathcal{T}(d_{k-1})$ uses the **previous** tangent
direction as a stand-in for transport, then **re-projects** $d_k$ onto the tangent
space at the current $C$ to remove drift from the CG combination.

### 8.3 Limited-memory BFGS (L-BFGS)

`EuclideanOptimizer::compute_lbfgs_direction` implements the **two-loop recursion**
(Nocedal & Wright) in the **standard Euclidean** metric of the working vector:

- After an accepted step, `update()` stores $s_k = \Delta x$ (the `step_vec` passed
  in by the solver) and $y_k = g_{k+1} - g_k$.
- Pairs are kept only if $s_k^\top y_k > 10^{-12}$; at most **`rdmft_lbfgs_memory`**
  pairs are retained (oldest dropped from the front of the deque).
- Initial inverse-Hessian diagonal scaling uses
  $\gamma_k = (s_{k-1}^\top y_{k-1}) / \|y_{k-1}\|^2$ from the **most recent** pair.

The returned direction is **minus** the two-loop result so that, for positive-definite
curvature information, $d_k$ is a descent direction in the Euclidean sense. For
**joint** optimisation this is applied to the packed $(p,\mathrm{flat}(C))$ vector;
for **alternating orbitals**, the line search still evaluates the energy along a
**retraction**, but the lbfgs history is built from **ambient** orbital increments
after the step—consistent with a pragmatic “projected lbfgs” pattern rather than a
full manifold-aware lbfgs transport.

### 8.4 Adam

`compute_adam_direction` maintains bias-corrected moments $(\hat m_k,\hat v_k)$ and sets

$$
d_k = -\texttt{adam\_lr}\; \hat m_k / (\sqrt{\hat v_k} + \texttt{adam\_eps})
$$

with $\beta_1,\beta_2$ from `RDMFTConfig` (`rdmft_adam_*` INPUT). The joint solver
uses $\alpha_0 = 1$ for the outer line search when `joint_optimizer` is Adam.

---

### 8.5 Line search routines (`armijo_line_search`, Wolfe family)

All are free functions in `rdmft_optimizer.h`. They share INPUT **`line_search_c1`**
(Armijo / Wolfe sufficient-decrease constant $c_1$, default $10^{-4}$),
**`line_search_c2`** (curvature parameter, default $0.9$),
**`line_search_max_iter`**, and for Wolfe searches **`line_search_max_zoom`**
(zoom-phase iteration cap).

#### 8.5.1 Armijo backtracking

Accepts the first $\alpha$ such that

$$
\phi(\alpha) \le \phi(0) + c_1 \,\alpha\, \phi'(0),
$$

where $\phi(\alpha)$ is the 1D objective along the search ray and $\phi'(0)$ is the
directional derivative $g^\top d$ supplied by the caller.

- If **`rdmft_line_search_polynomial`** is `true` (default), a failed trial uses a
  **quadratic** model through $(0,\phi(0),\phi'(0))$ and $(\alpha,\phi(\alpha))$ for the
  first failure, then **cubic** models using the last two failed points; the suggested
  next $\alpha$ is **clamped** to $[0.1,\,0.5]$ times the last failed step (fixed
  factors in code). If the model is unusable, fall back to $\alpha \leftarrow \rho\alpha$.
- If `false`, use **pure geometric** shrinking $\alpha \leftarrow \texttt{line\_search\_rho}\,\alpha$
  (default $\rho=\tfrac{1}{2}$).

Used for: **ALM** occupation steps when `rdmft_occ_ls_type` is `auto` or `armijo`;
**PG / AS** occupation monotone mode; **orbital SD / Adam**; and anywhere else the
solver calls `armijo_line_search` explicitly.

#### 8.5.2 Strong Wolfe

Implements Nocedal & Wright **Algorithm 3.5 / 3.6** (`strong_wolfe_line_search`):

1. **Sufficient decrease:** $\phi(\alpha) \le \phi(0) + c_1 \alpha \phi'(0)$.
2. **Strong curvature:** $|\phi'(\alpha)| \le c_2 |\phi'(0)|$.

The bracketing phase **doubles** $\alpha$ until violation, then **zooms** with a
cubic interpolant between bracket endpoints. Requires $\phi'(0) < 0$.

#### 8.5.3 Non-monotone strong Wolfe (Zhang–Hager-style reference)

`strong_wolfe_nm_line_search` is the same bracket–zoom machinery, but the sufficient
decrease test uses a **fixed** reference value $f_{\mathrm{ref}}$ instead of $\phi(0)$:

$$
\phi(\alpha) \le f_{\mathrm{ref}} + c_1 \alpha \phi'(0),
$$

with curvature conditions unchanged. In `rdmft_solver.cpp` this is wired to **PG /
AS occupation** line search when `rdmft_occ_ls_type = sw` **and**
`rdmft_occ_optimizer = cg`; $f_{\mathrm{ref}}$ is updated from a Zhang–Hager-like
non-monotone rule using **`rdmft_occ_cg_nonmonotone_eta`**.

#### 8.5.4 Weak Wolfe

`weak_wolfe_line_search` enforces

$$
\phi(\alpha) \le \phi(0) + c_1 \alpha \phi'(0),
\qquad
\phi'(\alpha) \ge c_2 \,\phi'(0)
$$

(with $\phi'(0) < 0$, the second inequality allows the slope to become **less negative**,
not necessarily small in absolute value). Used for **L-BFGS orbital** line search when
`rdmft_orb_ls_type` selects weak Wolfe.

---

### 8.6 INPUT mapping (`rdmft_occ_ls_type`, `rdmft_orb_ls_type`)

Parsed as `RdmftLineSearchPreset` (`auto`, `armijo`, `sw`, `wolfe`—see
`rdmft_input_parse.h`). Effective policies are computed in `rdmft_solver.cpp`.

**Occupations, augmented Lagrangian** (`effective_rdmft_ls_policy` — optimiser type
is **ignored**):

| `rdmft_occ_ls_type` | Line search on augmented Lagrangian $L$ |
| ------------------- | ---------------------------------------- |
| `auto` or `armijo`  | `armijo_line_search` (polynomial per `rdmft_line_search_polynomial`) |
| `sw`                | `strong_wolfe_line_search` |
| `wolfe`             | `weak_wolfe_line_search` |

**Occupations, projected gradient / active set** (`occ_proj_ls_mode_from_preset`):

| `rdmft_occ_ls_type` | `rdmft_occ_optimizer` | Mode |
| ------------------- | ---------------------- | ---- |
| `auto` or `armijo` | any | **Monotone** projected Armijo (`\alpha \leftarrow \rho\alpha` only—no polynomial branch) |
| `sw` | **cg** | **Non-monotone strong Wolfe** on the projected arc (`strong_wolfe_nm_line_search`) |
| `sw` | **sd / lbfgs / adam** | **Standard strong Wolfe** (`strong_wolfe_line_search`) |
| `wolfe` | any | **Weak Wolfe** (`weak_wolfe_line_search`) |

Additional PG knobs: **`rdmft_pg_occ_cg_ls_alpha_cap`** caps the initial Wolfe trial
when using CG; **`rdmft_pg_occ_ls_recovery_alpha`** controls a small Armijo recovery
after failed Wolfe (when enabled).

**Orbitals, alternating inner** (`effective_orbital_ls_policy`):

| `rdmft_orb_ls_type` | Optimiser | Line search |
| ------------------- | --------- | ----------- |
| `auto` | **cg** | Strong Wolfe |
| `auto` | **lbfgs** | Weak Wolfe |
| `auto` | **sd / adam** | Armijo (`armijo_line_search` with polynomial flag) |
| explicit `armijo` / `sw` / `wolfe` | — | Forces that family regardless of optimiser |

**Joint strategy** (`solve_joint`): a **single geometric Armijo** loop on the packed
step (§7.2); `rdmft_occ_ls_type` / `rdmft_orb_ls_type` do **not** switch the joint line
search to Wolfe. Initial $\alpha_0$ is **`rdmft_alm_bb_enabled`** Barzilai–Borwein on
the packed iterate when `joint_optimizer` is **sd/cg**; it is fixed to **1** for
**lbfgs/adam** joint runs.

**Barzilai–Borwein** (`BarzilaiBorweinStep` in `rdmft_optimizer.h`): cold start uses
$\alpha \approx \|x\|/\|g\|$; thereafter **BB1** $\|s\|^2/(s^\top y)$, **BB2**
$(s^\top y)/\|y\|^2$, or **alternate** between them (`rdmft_alm_bb_mode`), clamped to
[`rdmft_alm_bb_alpha_min`, `rdmft_alm_bb_alpha_max`]. Used for ALM occupation
`rdmft_occ_ls_init_step = bb`, joint SD/CG initial step when enabled, and PG/AS
**initial** Armijo trial via `rdmft_occ_ls_init_step` (see `rdmft_usage.md`).

---

## 9. Gradient consistency check (`rdmft_grad_check`)

When `RDMFTConfig::grad_check` is true, `RDMFTSolver::solve` runs finite-difference
checks **after** `precompute_cholesky_S()` and conversion of orbitals to internal
**X-space** (`wfc_C_to_X`), before the main optimisation loop.

- **Occupations:** central difference in each $n_{i\mathbf{k}}$ vs analytic
  `grad_occ` from `EnergyGradient::compute`.
- **Orbitals:** compares the **projected** Riemannian gradient $G_R$ to a **forward**
  difference along the same retraction used in line search,
  $(E(R_C(t\eta)) - E_0)/t$, which matches the Armijo slope $-\lVert G_R\rVert^2$
  more reliably than a symmetric $\pm t$ probe (`rdmft_usage.md`).

These checks are expensive but mirror the validation path used in development.

---

## 10. Gamma-Only Simplification

For $\Gamma$-point-only calculations ($\mathbf{k} = 0$):

- All matrices are real
- $C \in \mathbb{R}^{N \times N_b}$, use `double` instead of `complex<double>`
- The Stiefel manifold becomes $\mathrm{St}(N_b, N; S)$ over the reals
- The electron number constraint simplifies to $\sum_i n_i = N_e$ 

---

## 11. Spin-Polarized Case (nspin=2)

The derivation above uses a spin-restricted notation in which spin is either
suppressed or absorbed into an overall factor. For a collinear spin-polarized
calculation with `nspin=2`, the natural orbitals and occupation numbers carry
an additional spin label $\sigma \in \{\uparrow, \downarrow\}$, and the 1-RDM is
block-diagonal in spin.

### 11.1 Spin-Resolved 1-RDM

Introduce the combined coordinate $x = (\mathbf{r}, \sigma)$. Then

$$
\gamma(x, x') = \gamma(\mathbf{r}\sigma, \mathbf{r}'\sigma')
= \delta_{\sigma\sigma'} \, \gamma_{\sigma}(\mathbf{r}, \mathbf{r}')
$$

with

$$
\gamma_{\sigma}(\mathbf{r}, \mathbf{r}')
= \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_{i=1}^{N_b}
n_{i\mathbf{k}\sigma} \, \phi_{i\mathbf{k}\sigma}(\mathbf{r}) \, \phi_{i\mathbf{k}\sigma}^*(\mathbf{r}')
$$

The spin densities and magnetization density are

$$
\rho_{\sigma}(\mathbf{r}) = \gamma_{\sigma}(\mathbf{r}, \mathbf{r}),
\qquad
\rho(\mathbf{r}) = \rho_{\uparrow}(\mathbf{r}) + \rho_{\downarrow}(\mathbf{r}),
\qquad
m_z(\mathbf{r}) = \rho_{\uparrow}(\mathbf{r}) - \rho_{\downarrow}(\mathbf{r})
$$

In an LCAO basis with spin-independent spatial orbitals $\chi_{\mu\mathbf{k}}$
and orthonormal spinors $\zeta_{\sigma}$, we write

$$
\phi_{i\mathbf{k}\sigma}(\mathbf{r}, \sigma')
= \delta_{\sigma\sigma'} \sum_{\mu=1}^{N_{\mathrm{basis}}}
C_{\mu i}^{\mathbf{k}\sigma} \, \chi_{\mu\mathbf{k}}(\mathbf{r})
$$

so each spin channel satisfies its own orthonormality constraint:

$$
(C^{\mathbf{k}\sigma})^\dagger S^{\mathbf{k}} C^{\mathbf{k}\sigma} = I_{N_b}
$$

Hence, for `nspin=2`, the orbital manifold is the product

$$
\prod_{\sigma \in \{\uparrow,\downarrow\}} \prod_{\mathbf{k}}
\mathrm{St}(N_b, N_{\mathrm{basis}}; S^{\mathbf{k}})
$$

and the spin-resolved LCAO density matrix is

$$
\gamma_{\mu\nu}^{\mathbf{k}\sigma}
= C^{\mathbf{k}\sigma}
\operatorname{diag}(\mathbf{n}_{\mathbf{k}\sigma})
(C^{\mathbf{k}\sigma})^\dagger
$$

### 11.2 Total Energy Functional

The total energy becomes

$$
E[\{n_{i\mathbf{k}\sigma}, C^{\mathbf{k}\sigma}\}]
= E_{\mathrm{one}} + E_H[\rho] + E_{xc}[\gamma_{\uparrow}, \gamma_{\downarrow}] + E_{\mathrm{Ewald}}
$$

with $\rho = \rho_{\uparrow} + \rho_{\downarrow}$.

The one-body term is

$$
E_{\mathrm{one}}
= \sum_{\sigma} \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i
n_{i\mathbf{k}\sigma} \, h_{ii}^{\mathrm{one},\sigma}(\mathbf{k})
$$

In the usual collinear case without spin-orbit coupling, the kinetic and ionic
one-body operators are spin-independent, so the explicit $\sigma$ label on
$h^{\mathrm{one},\sigma}$ is often only a bookkeeping label.

The Hartree term depends on the total density:

$$
E_H[\rho]
= \frac{1}{2} \iint
\frac{[\rho_{\uparrow}(\mathbf{r}) + \rho_{\downarrow}(\mathbf{r})]
[\rho_{\uparrow}(\mathbf{r}') + \rho_{\downarrow}(\mathbf{r}')]}{|\mathbf{r} - \mathbf{r}'|}
d\mathbf{r} \, d\mathbf{r}'
$$

For the exchange-like RDMFT contribution, the **implementation** evaluates
Fock exchange from the spin-resolved modified density matrices via `Exx_LRI`
(§2.4) rather than assembling dense $K_{ij}$ four-centre integrals. In collinear
`nspin=2`, `EnergyGradient` passes `nspin` through to the density-matrix and EXX
pipelines so $\rho_\uparrow$, $\rho_\downarrow$, and the exchange Hamiltonian are
built with the same spin bookkeeping as KS-LCAO.

For each separable / single-channel build, the exchange energy is accumulated as
in §2.4 **per spin channel** (equal-spin exchange only):

$$
E_{xc} = \frac{1}{2} \sum_{\sigma} \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i
g(n_{i\mathbf{k}\sigma}) \,
\varepsilon^{\mathrm{exx}}_{i\mathbf{k}\sigma},
$$

with $\gamma_{\mathrm{xc}}^{\mathbf{k}\sigma}$ built using weights
$w_{\mathbf{k}}\, g(n_{i\mathbf{k}\sigma})$ as in the spinless formula.

**GU (implementation).** As in the spinless case, the explicit diagonal
$(n^2-n)\,J_{ii}$ GU correction is **not** added on top of the Müller-type EXX
build in `EnergyGradient::compute`; the same caveat applies spin-by-spin.

### 11.3 Constraints

The total electron-number constraint becomes

$$
\sum_{\sigma} \sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i n_{i\mathbf{k}\sigma} = N_e
$$

The ensemble $N$-representability bounds remain

$$
0 \le n_{i\mathbf{k}\sigma} \le 1
$$

and each spin channel has its own Stiefel constraint:

$$
(C^{\mathbf{k}\sigma})^\dagger S^{\mathbf{k}} C^{\mathbf{k}\sigma} = I
$$

If one additionally imposes a fixed collinear magnetization, then one adds the
constraint

$$
\sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i
\bigl(n_{i\mathbf{k}\uparrow} - n_{i\mathbf{k}\downarrow}\bigr) = M_z
$$

with a second Lagrange multiplier (or augmented-Lagrangian penalty) associated
with $M_z$.

### 11.4 Gradients

The occupation-number gradient matches §5.1 **per spin** (again with a **plus**
Hartree diagonal as implemented):

$$
\frac{\partial E}{\partial n_{i\mathbf{k}\sigma}}
= w_{\mathbf{k}} \Bigl(
h_{ii}^{\mathrm{one},\sigma}(\mathbf{k}) + v_{H,ii}^{\sigma}(\mathbf{k})
\Bigr)
+ \frac{\partial E_{xc}}{\partial n_{i\mathbf{k}\sigma}} .
$$

For separable / single-channel exchange,

$$
\frac{\partial E_{xc}}{\partial n_{i\mathbf{k}\sigma}}
= w_{\mathbf{k}} \, g'(n_{i\mathbf{k}\sigma}) \,
\varepsilon^{\mathrm{exx}}_{i\mathbf{k}\sigma} .
$$

Mixed functionals (GEO / CHF / CGA / HybOpt) use the same channel-summing scheme
as §5.1, independently per $\sigma$.

The Euclidean orbital gradient for each spin block matches §5.2, including the
explicit factor **2** from the Wirtinger / real-pairing convention:

$$
G^{\mathbf{k}\sigma}
= 2 \, w_{\mathbf{k}} \Bigl[
\bigl(h^{\mathbf{k}\sigma} + V_H^{\mathbf{k}}\bigr)
C^{\mathbf{k}\sigma} \operatorname{diag}(\mathbf{n}_{\mathbf{k}\sigma})
+ \text{(mixed or } g\cdot H_{\mathrm{exx}} \text{ exchange column)}
\Bigr].
$$

The Riemannian projection and retraction are then applied independently to each
$C^{\mathbf{k}\sigma}$:

$$
\operatorname{grad} E^{\mathbf{k}\sigma}
= \Pi_{T_{C^{\mathbf{k}\sigma}}\mathrm{St}}\bigl(G^{\mathbf{k}\sigma}\bigr),
\qquad
C_{\mathrm{new}}^{\mathbf{k}\sigma}
= R_{C^{\mathbf{k}\sigma}}\bigl(\eta^{\mathbf{k}\sigma}\bigr)
$$

The occupation parameterizations from §4 (including sigma-shift, §4.3) extend
by carrying the spin index on $z$, $\theta$, or $x$ when applicable.

### 11.5 Reduction to the current implementation

`EnergyGradient::init` sets `nk_` from `pelec->wg.nr`, i.e. the same **flattened
$(\mathbf{k},\sigma)$ row count** used for KS occupation weights in the running
calculation. The occupation vector `occ_flat` therefore has length
`nk_ * nbands` with the ABACUS LCAO convention for collinear spin (rather than a
separate explicit list of $\sigma$ indices in the solver API).

The spin-resolved formulas in §11.1–§11.4 are still the correct physics; they map
onto the code’s layout after identifying how `ik` encodes spin in the current
`K_Vectors` / `wg` convention.

---

## 12. Summary of Implementation Requirements


| Component                         | Variables | Manifold / discretisation        | Where in code / this note        |
| --------------------------------- | --------- | -------------------------------- | -------------------------------- |
| One-body + Hartree                | n, C      | —                                | §5.1–§5.2, `EnergyGradient::compute` |
| XC (HF / Müller / Power)          | n, C      | Modified DM + `Exx_LRI`          | §2.4, §5.1–§5.2                  |
| XC (GU)                           | n, C      | Same EXX path as Müller ($g=\sqrt{n}$ reg.); analytic $(n^2-n)J_{ii}$ **not** in `compute` | §2.4 GU remark |
| XC (GEO / CHF / CGA / HybOpt)      | n, C      | Multiple EXX builds, summed      | §2.4 mixed channels              |
| HF occupation entropy (optional)  | n         | —                                | §5.1, `occ_entropy_gamma_`       |
| Occupation constraints            | n         | ALM / PG / active-set + §4.3 shift | §4, §6                         |
| Orbital orthogonality             | C (X)     | Stiefel: §5.3–§5.4               | `StiefelManifold`, `project_orbital_gradient` |


