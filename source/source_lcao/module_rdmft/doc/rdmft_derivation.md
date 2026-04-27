# Reduced Density Matrix Functional Theory (RDMFT): Mathematical Derivation

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

### Planewave Expansion

For planewave (PW) basis, each natural orbital is:

$$
\phi_{i\mathbf{k}}(\mathbf{r}) = \frac{1}{\sqrt{\Omega}} \sum_{\mathbf{G}} c_{i\mathbf{k}}(\mathbf{G})  e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}
$$

where $\Omega$ is the unit cell volume and $\mathbf{G}$ are reciprocal lattice vectors within
the kinetic energy cutoff. Orthonormality is simply
$\sum_{\mathbf{G}} c_{i\mathbf{k}}^*(\mathbf{G}) c_{j\mathbf{k}}(\mathbf{G}) = \delta_{ij}$,
so the coefficients lie on the standard Stiefel manifold $\mathrm{St}(N_b, N_{\mathrm{pw}}^{\mathbf{k}})$.

---

## 2. Total Energy Functional

$$
E[n_{i\mathbf{k}}, C^{\mathbf{k}}] = E_{\mathrm{one}} + E_H[\rho] + E_{xc}[\gamma] + E_{\mathrm{Ewald}}
$$

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
We consider four functionals, all defined through the coupling function
$f(n_i, n_j)$ appearing in the exchange-like integral:

$$
E_{xc}[\gamma] = -\frac{1}{2} \sum_{\mathbf{k}\mathbf{k}'} w_{\mathbf{k}} w_{\mathbf{k}'}
\sum_{ij} f(n_{i\mathbf{k}}, n_{j\mathbf{k}'})  K_{ij}^{\mathbf{k}\mathbf{k}'}
$$

where $K_{ij}^{\mathbf{k}\mathbf{k}'} = \langle \phi_{i\mathbf{k}} \phi_{j\mathbf{k}'} | \hat{v}_{c} | \phi_{j\mathbf{k}'} \phi_{i\mathbf{k}} \rangle$
are two-electron exchange integrals.

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

### 2.4 Unified Exchange via Modified Density Matrix

For the separable functionals (HF, Müller, Power), we can write
$f(n_i, n_j) = g(n_i) g(n_j)$, and define the **modified density matrix**:

$$
\gamma_{\mathrm{xc},\mu\nu}^{\mathbf{k}} = \sum_i w_{\mathbf{k}}  g(n_{i\mathbf{k}})  C_{\mu i}^{\mathbf{k}} (C_{\nu i}^{\mathbf{k}})^*
$$

The exchange energy is then

$$
E_{xc} = -\frac{1}{2} \mathrm{Tr}[\gamma_{\mathrm{xc}}  K  \gamma_{\mathrm{xc}}]
= \frac{1}{2} \mathrm{Tr}[\gamma_{\mathrm{xc}}  H_{\mathrm{exx}}[\gamma_{\mathrm{xc}}]]
$$

where $H_{\mathrm{exx}}[\gamma_{\mathrm{xc}}]$ is the exchange Hamiltonian computed from
$\gamma_{\mathrm{xc}}$ using LibRI.

For the GU functional, the energy has an additional correction:

$$
E_{xc}^{\mathrm{GU}} = E_{xc}^{\mathrm{Müller}} + \frac{1}{2} \sum_{\mathbf{k}} w_{\mathbf{k}}^2 \sum_i (n_{i\mathbf{k}}^2 - n_{i\mathbf{k}})  J_{ii}^{\mathbf{k}\mathbf{k}}
$$

where $J_{ii}^{\mathbf{k}\mathbf{k}} = K_{ii}^{\mathbf{k}\mathbf{k}}$ is the self-exchange integral.

---

## 3. Constraints

### 3.1 Electron Number Conservation

$$
\sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i n_{i\mathbf{k}} = N_e
$$

where $N_e$ is the number of electrons per unit cell.

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

### 4.3 Transformed Gradient

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

### 5.1 Gradient w.r.t. Natural Occupation Numbers

$$
\frac{\partial E}{\partial n_{i\mathbf{k}}} = w_{\mathbf{k}}  h_{ii}^{\mathrm{one}}(\mathbf{k})
- w_{\mathbf{k}}  v_{H,ii}(\mathbf{k})
- \frac{\partial E_{xc}}{\partial n_{i\mathbf{k}}}
$$

where $v_{H,ii}(\mathbf{k}) = (C^{\mathbf{k}})^\dagger V_H^{\mathbf{k}} C^{\mathbf{k}} |_{ii}$.

For the separable functionals (HF, Müller, Power):

$$
\frac{\partial E_{xc}}{\partial n_{i\mathbf{k}}} = w_{\mathbf{k}}  g'(n_{i\mathbf{k}})  \langle \phi_{i\mathbf{k}} | H_{\mathrm{exx}}[\gamma_{\mathrm{xc}}] | \phi_{i\mathbf{k}} \rangle
$$

For the GU functional, there is an additional self-interaction correction term:

$$
\frac{\partial E_{xc}^{\mathrm{GU}}}{\partial n_{i\mathbf{k}}} = w_{\mathbf{k}} \Bigl[
g'_{\mathrm{M}}(n_{i\mathbf{k}})  h_{\mathrm{exx},ii}^{\mathbf{k}}
- w_{\mathbf{k}} (2n_{i\mathbf{k}} - 1)  J_{ii}^{\mathbf{k}\mathbf{k}}
\Bigr]
$$

### Derivatives of g(n):


| Functional | $g(n)$     | $g'(n)$                |
| ---------- | -------- | -------------------- |
| HF         | $n$        | $1$                    |
| Müller     | $n^{1/2}$  | $\frac{1}{2} n^{-1/2}$ |
| Power      | $n^\alpha$ | $\alpha n^{\alpha-1}$  |


### 5.2 Gradient w.r.t. Orbital Coefficients (Euclidean)

The Euclidean gradient in the ambient space is:

$$
G_{\mu i}^{\mathbf{k}} = \frac{\partial E}{\partial (C_{\mu i}^{\mathbf{k}})^*}
= w_{\mathbf{k}} \Bigl[
n_{i\mathbf{k}} \bigl(h^{\mathbf{k}} + V_H^{\mathbf{k}}\bigr) C^{\mathbf{k}} \big|_{\mu i}

- g(n_{i\mathbf{k}})  H_{\mathrm{exx}}^{\mathbf{k}} C^{\mathbf{k}} \big|_{\mu i}
\Bigr]
$$

In compact notation:

$$
G^{\mathbf{k}} = w_{\mathbf{k}} \Bigl[
\bigl(h^{\mathbf{k}} + V_H^{\mathbf{k}}\bigr) C^{\mathbf{k}}  \mathrm{diag}(\mathbf{n}_{\mathbf{k}})

- H_{\mathrm{exx}}^{\mathbf{k}} C^{\mathbf{k}}  \mathrm{diag}(g(\mathbf{n}_{\mathbf{k}}))
\Bigr]
$$

For GU, the exchange Hamiltonian $H_{\mathrm{exx}}$ is from the Müller DM,
plus a diagonal self-interaction correction:

$$
G^{\mathbf{k}}_{\mathrm{GU}} = G^{\mathbf{k}}_{\mathrm{Müller}}
- w_{\mathbf{k}}^2 \sum_i (n_{i\mathbf{k}}^2 - n_{i\mathbf{k}}) \frac{\partial J_{ii}}{\partial (C^{\mathbf{k}})^*}
$$

### 5.3 Riemannian Gradient on Stiefel Manifold

The Stiefel manifold $\mathrm{St}(N_b, N; S)$ with metric induced by $S$ has the
tangent space at $C$:

$$
T_C \mathrm{St} = \{ Z \in \mathbb{C}^{N \times N_b} : C^\dagger S Z + Z^\dagger S C = 0 \}
$$

The Riemannian gradient is the projection of the Euclidean gradient onto the tangent space:

$$
\mathrm{grad} E = G - S^{-1} C  \mathrm{sym}(C^\dagger G)
$$

where $\mathrm{sym}(A) = \frac{1}{2}(A + A^\dagger)$.

For the **canonical metric** on Stiefel:

$$
\mathrm{grad} E = G - C (C^\dagger S G)_{\mathrm{sym}}
$$

Wait, more precisely, with the metric $\langle Z_1, Z_2 \rangle = \mathrm{Re}\,\mathrm{Tr}(Z_1^\dagger S Z_2)$:

$$
\mathrm{grad} E = S^{-1} G - C  \mathrm{sym}(C^\dagger G)
$$

### 5.4 Retraction on Stiefel Manifold

Given a tangent vector $\eta \in T_C \mathrm{St}$, the QR-based retraction is:

$$
R_C(\eta) = \mathrm{qf}(C + \eta)
$$

where $\mathrm{qf}$ denotes the Q-factor of the QR decomposition (with positive diagonal in R),
followed by S-orthogonalization: solve $S^{1/2} (C+\eta) = QR$, return $S^{-1/2} Q$.

For the **polar retraction**:

$$
R_C(\eta) = (C + \eta) \bigl[ (C + \eta)^\dagger S (C + \eta) \bigr]^{-1/2}
$$

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
backtracking is chosen by INPUT `rdmft_occ_ls_init_step`: fixed $1$, Barzilai–Borwein
(uses `rdmft_alm_bb_mode` and clamps `rdmft_alm_bb_alpha_min` /
`rdmft_alm_bb_alpha_max`, with fallback to `rdmft_alpha_step`; unlike the ALM
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
with $\tau_{\mathrm{ls}} =$ `rdmft_alpha_step`, resets the occupation optimizer
curvature state (CG / L-BFGS history), and continues.

**Stopping criteria (ABACUS inner loop).**  Let $\mathbf{g}_{\mathrm{proj}} =
\bigl(\mathbf{n} - P_{\mathcal{C}}(\mathbf{n} - \tau \mathbf{g})\bigr)/\tau$
(Bertsekas projected-gradient vector).  **$\tau$ is taken from the occupation
line search** so the stationarity measure uses the same scale as the step:
**pre-step**, $\tau = \alpha_0$ (the first Armijo trial from
`rdmft_occ_ls_init_step`: fixed $1$, Barzilai–Borwein, or quad, with fallback to
`rdmft_alpha_step` if the estimate is non-positive or non-finite); **post-step**
after a successful line search, $\tau = \alpha_{\mathrm{acc}}$ (accepted
Armijo step along $\mathbf{d}$); after **SD fallback** (failed line search),
$\tau =$ `rdmft_alpha_step` (same step length as the recovery move).  The inner
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
  `rdmft_alpha_step` (fallback for invalid $\alpha_0$ and SD-fallback step length),
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
- For orbital constraints (Stiefel), prefer Riemannian retraction (QR or polar)
   rather than Euclidean clipping; see §5.4.

Pseudocode (ABACUS projected-gradient inner loop; `project` denotes
$P_{\mathcal{C}}$; neglect spin labels):

```text
initialize n feasible; configure occ_optimizer, tol = rdmft_occ_grad_tol
for inner = 0 .. occ_maxiter-1:
   g = grad_n(E, n)
   alpha0 = initial_step_from(rdmft_occ_ls_init_step, BB/quad, rdmft_alpha_step, clamps)
   tau_pre = alpha0 if alpha0 > 0 else rdmft_alpha_step
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
      n = project(n_save - rdmft_alpha_step * g)         // SD fallback
      alpha_acc = rdmft_alpha_step
      reset_optimizer_curvature_state()
   recompute g at new n
   tau_post = alpha_acc                                   // accepted α, or rdmft_alpha_step after fallback
   g_proj_post = (n - project(n - tau_post * g)) / tau_post
   if ||g_proj_post||_inf < tol: break
end
```

When not to use:

- If the active-set (set of variables at bounds) is small and changes infrequently,
   an active-set method or a second-order reduced Newton solve on the free set
   may converge much faster near the solution.

### 6.3 Active Set Method

Maintain active sets
$$\mathcal{A}_0 = \{I:\; n_I = 0\},\qquad \mathcal{A}_1 = \{I:\; n_I = 1\}$$
and the free set $\mathcal{F} = \{I:\; 0 < n_I < 1\}$, where again $I=(i,\mathbf{k})$.

Overview:

- The active-set method iteratively guesses which bounds are active (occupied at
   0 or 1) and solves a reduced equality-constrained optimization on the free set.
   It then updates Lagrange multipliers for the active constraints and adjusts the
   active set until KKT conditions are satisfied.

Core algorithm (bound/simplex case):

1. Choose an initial active set (for example from the current projected-gradient iterate).
2. Solve the reduced problem on free variables: minimize $E(n)$ subject to
    $\sum_{I\in\mathcal{F}} w_I n_I = N_e - \sum_{I\in\mathcal{A}_1} w_I$ and
    $n_I$ fixed at 0 or 1 on active indices. This can be done via a Newton step
    on free variables or by solving the KKT linear system for a quadratic model.
3. If the step violates a bound for some free index, move along the step until
    the first bound is hit; add that index to the corresponding active set and go to 2.
4. Compute multipliers $\lambda_I$ for active constraints. If any multiplier
    violates complementarity (wrong sign), remove its constraint from the active
    set and go to 2.
5. Stop when primal feasibility, complementary slackness and dual feasibility
    (KKT residuals) are below tolerances.

Pseudocode (sketch):

```text
initialize n, form A0, A1, F
while not converged:
   solve reduced Newton system on F (or perform CG on Hessian-free model)
   compute candidate step and max step length before hitting bounds
   if bound hit:
      step to bound, add index to A0 or A1
      continue
   accept full step
   compute multipliers for active constraints
   if any multiplier violates sign condition:
      remove violating index from active set
      continue
   check KKT residuals -> break if small
end
```

Computing multipliers and KKT system:

- If the reduced problem is solved by Newton, form the KKT linear system
   (H_F  A^T; A 0) for Hessian on free set $H_F$ and equality constraint matrix
   $A$ (the weighted-sum row). Solve for primal step and multiplier update.
- For large systems use iterative solvers (CG, MINRES) preconditioned by a
   diagonal or limited-memory factor.

Numerical tips:

- Warm-start linear solves: cache factorizations of the reduced Hessian and update
   incrementally when the active set changes.
- Use limited-memory quasi-Newton (lbfgs) on the free set if exact Hessians are
   expensive; form a small KKT system for the equality constraint.
- Add trust-region safeguards or fallback to projected-gradient when the
   reduced-step increases the objective (nonconvexity caution).

When to prefer active-set:

- When only a small fraction of occupations are at the bounds (sparse active set),
   the reduced Newton/QUASI-NEWTON solves can converge in very few outer iterations
   and achieve fast (superlinear) local convergence.

Hybrid strategies:

- A practical pattern is to run projected-gradient iterations to approach
   a neighborhood of the solution, then switch to an active-set solver to
   enforce exact complementary slackness and remove the residual projected gradient.
- For RDMFT: use projected gradient for several outer iterations, detect when
   many occupations settle near 0 or 1, then invoke active-set on the remaining
   free occupations while holding orbitals fixed (or solved together in a
   reduced joint solve).

Stopping and tolerances:

- KKT residual tolerances for active-set: primal feasibility ~1e-8–1e-6,
   dual complementarity ~1e-6–1e-4 depending on problem scale.
- Use looser tolerances during early iterations and tighten near convergence.

Examples and diagnostics:

- Log active-set entries and multiplier signs each iteration to diagnose
   oscillations (add/remove cycles). If oscillations occur, increase damping
   or use a small trust-region.
- Compare final active-set with projected-gradient saturations to validate
   the hybrid strategy.


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

With cosine-squared parameterization, the occupation part becomes unconstrained Euclidean
(plus the augmented Lagrangian for electron number). The product manifold approach optimizes
all variables simultaneously using Riemannian optimization.

The tangent vector at a point $(\boldsymbol{\theta}, C^{\mathbf{k}})$ is
$(\delta\boldsymbol{\theta}, \eta^{\mathbf{k}})$ where
$\delta\boldsymbol{\theta} \in \mathbb{R}^{N_k \times N_b}$ and
$\eta^{\mathbf{k}} \in T_{C^{\mathbf{k}}} \mathrm{St}$.

Retraction: apply Euclidean update to $\boldsymbol{\theta}$ and Stiefel retraction to each $C^{\mathbf{k}}$.

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
5. performs a single Armijo backtracking line search along the packed
   direction: the occupation parameters are updated linearly
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

## 8. Optimization Algorithms

All algorithms below work generically on the manifold $\mathcal{M}$ (or its submanifolds
in the alternating case).

### 8.1 Steepest Descent (SD)

$$
x_{t+1} = R_{x_t}(-\alpha_t  \mathrm{grad} f(x_t))
$$

with step size $\alpha_t$ chosen by line search (Armijo backtracking).

### 8.2 Conjugate Gradient (CG)

Riemannian CG with vector transport $\mathcal{T}$:

$$
\eta_t = -\mathrm{grad} f(x_t) + \beta_t  \mathcal{T}_{x_{t-1} \to x_t}(\eta_{t-1})
$$

with Fletcher-Reeves or Polak-Ribière $\beta_t$.

### 8.3 lbfgs

limited-memory lbfgs adapted to Riemannian setting using vector transport
to move previous gradients and steps to the current tangent space.

### 8.4 Adam

Riemannian Adam with bias-corrected first and second moment estimates,
transported to the current tangent space.

---

## 9. Gradient Consistency Check

For numerical verification, use finite differences:

$$
\frac{\partial E}{\partial n_{i\mathbf{k}}} \approx \frac{E(n_{i\mathbf{k}} + \epsilon) - E(n_{i\mathbf{k}} - \epsilon)}{2\epsilon}
$$

For orbital gradients, perturb along a tangent direction $\eta$:

$$
\langle \mathrm{grad} E, \eta \rangle \approx \frac{E(R_C(t\eta)) - E(R_C(-t\eta))}{2t}
$$

These checks are critical for validating the implementation before running production calculations.

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

For the exchange-like RDMFT contribution, only equal-spin pairs contribute:

$$
E_{xc}[\gamma_{\uparrow}, \gamma_{\downarrow}]
= -\frac{1}{2} \sum_{\sigma}
\sum_{\mathbf{k}\mathbf{k}'} w_{\mathbf{k}} w_{\mathbf{k}'}
\sum_{ij} f(n_{i\mathbf{k}\sigma}, n_{j\mathbf{k}'\sigma})
K_{ij}^{\mathbf{k}\mathbf{k}',\sigma}
$$

because the spin functions are orthogonal and there is therefore no exchange
between $\uparrow$ and $\downarrow$ blocks.

For the separable functionals (HF, Müller, Power), define the spin-resolved
modified density matrix

$$
\gamma_{\mathrm{xc},\mu\nu}^{\mathbf{k}\sigma}
= \sum_i w_{\mathbf{k}} \, g(n_{i\mathbf{k}\sigma})
C_{\mu i}^{\mathbf{k}\sigma} (C_{\nu i}^{\mathbf{k}\sigma})^*
$$

Then the exchange energy is the sum of the two spin-channel contributions:

$$
E_{xc} = \frac{1}{2} \sum_{\sigma}
\operatorname{Tr}\bigl[\gamma_{\mathrm{xc},\sigma}
H_{\mathrm{exx}}[\gamma_{\mathrm{xc},\sigma}]\bigr]
$$

For the GU functional, the self-interaction correction is also spin-resolved:

$$
E_{xc}^{\mathrm{GU}}
= E_{xc}^{\mathrm{Müller}}
+ \frac{1}{2} \sum_{\sigma} \sum_{\mathbf{k}} w_{\mathbf{k}}^2 \sum_i
(n_{i\mathbf{k}\sigma}^2 - n_{i\mathbf{k}\sigma}) J_{ii}^{\mathbf{k}\mathbf{k},\sigma}
$$

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

The occupation-number gradient is the direct spin-resolved analogue of §5.1:

$$
\frac{\partial E}{\partial n_{i\mathbf{k}\sigma}}
= w_{\mathbf{k}} \, h_{ii}^{\mathrm{one},\sigma}(\mathbf{k})
- w_{\mathbf{k}} \, v_{H,ii}^{\sigma}(\mathbf{k})
- \frac{\partial E_{xc}}{\partial n_{i\mathbf{k}\sigma}}
$$

For the separable functionals,

$$
\frac{\partial E_{xc}}{\partial n_{i\mathbf{k}\sigma}}
= w_{\mathbf{k}} \, g'(n_{i\mathbf{k}\sigma})
\langle \phi_{i\mathbf{k}\sigma}
| H_{\mathrm{exx}}[\gamma_{\mathrm{xc},\sigma}] |
\phi_{i\mathbf{k}\sigma} \rangle
$$

The Euclidean orbital gradient for each spin block is likewise

$$
G^{\mathbf{k}\sigma}
= w_{\mathbf{k}} \Bigl[
\bigl(h^{\mathbf{k}\sigma} + V_H^{\mathbf{k}}\bigr)
C^{\mathbf{k}\sigma} \operatorname{diag}(\mathbf{n}_{\mathbf{k}\sigma})
- H_{\mathrm{exx}}^{\mathbf{k}\sigma}
C^{\mathbf{k}\sigma} \operatorname{diag}(g(\mathbf{n}_{\mathbf{k}\sigma}))
\Bigr]
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

The occupation parameterizations from §4 also extend trivially by carrying the
spin index:

$$
n_{i\mathbf{k}\sigma} = \cos^2(\theta_{i\mathbf{k}\sigma})
\qquad \text{or} \qquad
n_{i\mathbf{k}\sigma} = \frac{1}{1 + e^{-x_{i\mathbf{k}\sigma}}}
$$

### 11.5 Reduction to the Current Implementation

For `nspin=2`, the implementation can be viewed as folding the spin label into
a composite index $I = (\sigma, \mathbf{k})$. Equivalently,

$$
\sum_I \equiv \sum_{\sigma} \sum_{\mathbf{k}}
$$

Then the spin-polarized formulas reduce to the same algebraic form as the
spin-restricted ones after the replacements

$$
n_{i\mathbf{k}} \to n_{iI},
\qquad
C^{\mathbf{k}} \to C^I,
\qquad
w_{\mathbf{k}} \to w_I = w_{\mathbf{k}}
$$

The only physical caveat is that the exchange operator remains block-diagonal in
spin, so $H_{\mathrm{exx}}[\gamma_{\mathrm{xc}}]$ is built independently for the
$\uparrow$ and $\downarrow$ channels before their contributions are summed.

---

## 12. Summary of Implementation Requirements


| Component                   | Variables | Manifold               | Gradient                |
| --------------------------- | --------- | ---------------------- | ----------------------- |
| One-body energy             | n, C      | —                      | Eqs. in §5.1, §5.2, §11 |
| Hartree energy              | n, C      | —                      | Eqs. in §5.1, §5.2, §11 |
| XC energy (HF/Müller/Power) | n, C      | —                      | Eqs. in §5.1, §5.2, §11 |
| XC energy (GU)              | n, C      | —                      | Additional SIC terms, §11 |
| Occupation constraint       | n         | $[0,1]$, $\sum w n = N_e$ | §4, §6, §11            |
| Orbital orthogonality       | C         | Stiefel                | §5.3, §5.4, §11        |


