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

#### 2.3.5 GEO Functional

$$
f^{\mathrm{GEO}}(n_p, n_q) = \frac{1}{4}\Bigl[ n_p n_q + (n_p n_q)^{1/2} + 2 (n_p n_q)^{3/4} \Bigr]
$$

This is **non-separable in the single-$g$ sense** but decomposes as a weighted sum
of three separable Power-like terms:

$$
f^{\mathrm{GEO}}(n_p, n_q) = \frac{1}{4} g_1(n_p) g_1(n_q)
                          + \frac{1}{4} g_{1/2}(n_p) g_{1/2}(n_q)
                          + \frac{1}{2} g_{3/4}(n_p) g_{3/4}(n_q),
\qquad g_\alpha(n) = n^{\alpha}.
$$

The triplet $\{(c_t, \alpha_t)\} = \{(1/4, 1), (1/4, 1/2), (1/2, 3/4)\}$ is the
constructive recipe used by the implementation: at each energy / gradient
evaluation, the three modified density matrices
$\gamma^{(t)}_{\mathrm{xc}} = \sum_i w_{\mathbf{k}} n_{i\mathbf{k}}^{\alpha_t}
C_{\mu i}^{\mathbf{k}}(C_{\nu i}^{\mathbf{k}})^*$ are built and three EXX
evaluations are performed; energies, occupation gradients, and orbital
gradients are accumulated weighted by $c_t$ and the appropriate
$g_{\alpha_t}(n)$ or $g'_{\alpha_t}(n)$.

By construction $f^{\mathrm{GEO}}(1,1) = 1$ (the coefficients
$1/4 + 1/4 + 1/2 = 1$) and $f^{\mathrm{GEO}}(0, n_q) = 0$ analytically,
matching the HF/Müller/Power behaviour at the $n=0$ and $n=1$ extremes.

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

After the Cholesky variable change $X = U C$ (with $S = U^\dagger U$, see §3.3),
the constraint becomes $X^\dagger X = I$ and the orbital sub-problem lives on
the standard Stiefel manifold $\mathrm{St}(N_b, N)$.  Three retractions are
implemented; the choice is selected at run time by the `rdmft_orb_retraction`
INPUT keyword (default `polar`).  All preserve $X_{\text{new}}^\dagger X_{\text{new}} = I$
to machine precision when applicable.

| `rdmft_orb_retraction` | Formula (one Stiefel step with tangent $\eta = -\alpha\, G_R$) | Cost / robustness |
|---|---|---|
| **`polar`** (default) | $R_X(\eta) = (X + \eta)\bigl[(X+\eta)^\dagger (X+\eta)\bigr]^{-1/2}$, computed as Cholesky-QR: $M = (X+\eta)^\dagger (X+\eta) = L L^\dagger$, then $X_{\text{new}} = (X+\eta)\, L^{-\dagger}$. | One Cholesky + one triangular solve per k-point. Fully MPI-parallelised (ScaLAPACK `pdpotrf` + `pdtrsm`). Loses about half the working precision when $M$ is ill-conditioned, e.g. when $\alpha \|\eta\|$ is too large; the outer line search keeps $\alpha$ small enough that this rarely matters. |
| **`qr`** | $R_X(\eta) = Q\, \mathrm{diag}\bigl(\mathrm{phase}(R_{ii})\bigr)$, where $X+\eta = Q R$ is the Householder QR (LAPACK `?geqrf` + `?orgqr`). The diagonal sign-fix makes $R$ have positive real diagonal so the retraction is uniquely defined. | More numerically robust than `polar` when $(X+\eta)^\dagger(X+\eta)$ is near-singular. Serial-only; MPI builds emit a one-time warning and fall back to `polar`. |
| **`cayley`** | Wen-Yin low-rank Cayley retraction (Wen & Yin, *Math. Prog.* 142 (2013) 397, Algorithm 1). With $W = G_R X^\dagger - X G_R^\dagger$ rank-$2p$, $W = U V^\dagger$ for $U = [G_R \mid X]$, $V = [X \mid -G_R]$ (each $N \times 2p$): $$X_{\text{new}} = X - \alpha\, U \bigl(I_{2p} + \tfrac{\alpha}{2} V^\dagger U\bigr)^{-1}\, V^\dagger X.$$ | Solves only one $2p \times 2p$ system; preserves orthogonality exactly without a triangular factor. Attractive when $N \gg p$. Serial-only; MPI builds emit a one-time warning and fall back to `polar`. |

The retraction choice affects every orbital step in both the alternating
strategy (§7.1) and the joint strategy (§7.2); a single
`EnergyGradient::set_orb_retraction` call during solver initialisation forwards
the INPUT value to `EnergyGradient::retract_orbitals`, which dispatches to
`retract_polar`, `retract_qr_serial`, or `retract_cayley_serial`.

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

### 6.2 Spectral Projected Gradient (SPG)

The `rdmft_constraint = projected_gradient` (and the `active_set` alias, which
routes here as well) optimisation of $n_{i\mathbf{k}}$ on the feasible set
$\mathcal{C}$ uses the **Spectral Projected Gradient** method of
Birgin–Martínez–Raydan (SIAM J. Optim. 10 (2000) 1196), with the non-monotone
Armijo line search of Grippo–Lampariello–Lucidi (SIAM J. Numer. Anal. 23 (1986)
707).  This is the textbook standard for box-and-equality-constrained smooth
optimisation; it has none of the trust-region caps, multi-stage SD-fallback
chains, multiple `α₀` policies, or per-trial re-projection rejections of the
earlier ABACUS PG/AS implementation.

**Per inner iteration $k$:**

1. Compute $\mathbf{g}_k = \partial E/\partial \mathbf{n}$ at $\mathbf{n}_k$.
2. Spectral (BB1) step length
   $$
   \alpha_k^{\mathrm{BB}} = \frac{\langle \mathbf{s}_k, \mathbf{s}_k\rangle}{\langle \mathbf{s}_k, \mathbf{y}_k\rangle},\qquad
   \mathbf{s}_k = \mathbf{n}_k - \mathbf{n}_{k-1},\;\; \mathbf{y}_k = \mathbf{g}_k - \mathbf{g}_{k-1},
   $$
   safeguarded into $[\alpha_{\min}, \alpha_{\max}]$.  At the first inner step the
   fallback is $\alpha_0 = 1/\max(1, \|\mathbf{g}\|_\infty)$.
3. **Spectral projected direction**
   $$
   \mathbf{d}_k = P_{\mathcal{C}}\bigl(\mathbf{n}_k - \alpha_k^{\mathrm{BB}}\, \mathbf{g}_k\bigr) - \mathbf{n}_k.
   $$
   For any closed convex $\mathcal{C}$ this is provably a descent direction with
   $$
   \mathbf{g}_k^\top \mathbf{d}_k \le -\frac{\|\mathbf{d}_k\|^2}{\alpha_k^{\mathrm{BB}}} \le 0
   $$
   (BMR Lemma 2.1).  No descent-direction safeguard or trust-region cap is
   needed.
4. **Non-monotone Armijo** along the convex segment
   $\mathbf{n}_k(\lambda) = \mathbf{n}_k + \lambda\, \mathbf{d}_k$,
   $\lambda \in (0, 1]$: find the smallest $j \ge 0$ such that
   $$
   E\bigl(\mathbf{n}_k(\rho^j)\bigr)
   \le f_{\max} + c_1 \rho^j\, \mathbf{g}_k^\top \mathbf{d}_k,
   \qquad
   f_{\max} = \max_{0 \le i \le \min(k,\, M-1)} E(\mathbf{n}_{k-i}),
   $$
   with $c_1$ from `rdmft_line_search_c1`, $\rho$ from `rdmft_line_search_rho`,
   and a fixed history length $M = 10$ (the value used in the original SPG
   paper).  Each trial point $\mathbf{n}_k(\lambda)$ is the convex combination
   $(1-\lambda)\mathbf{n}_k + \lambda P_{\mathcal{C}}(\cdots)$ of two feasible
   points; convexity of $\mathcal{C}$ guarantees $\mathbf{n}_k(\lambda) \in
   \mathcal{C}$ without re-projection.

`rdmft_occ_optimizer` (sd / cg) is **not** consulted on this
path: SPG already has a proven globally convergent step length.

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

**Stopping criterion.**  Let
$\mathbf{r}(\mathbf{n}) = \mathbf{n} - P_{\mathcal{C}}(\mathbf{n} - \nabla_{\mathbf{n}}E)$
be the **τ-free Bertsekas projected-gradient residual** at unit step (the
textbook SPG stationarity measure; BMR Eq. (2.6)).  At a KKT point of the
constrained problem $\mathbf{r} = 0$.  The inner loop stops when
$\|\mathbf{r}\|_\infty \le$ `rdmft_occ_proj_tol`.  The iteration cap is `rdmft_occ_maxiter`.

**INPUT keywords (SPG path):** `rdmft_occ_proj_tol`,
`rdmft_occ_maxiter`, `rdmft_line_search_c1`.  Backtracking uses a fixed ratio
from `RDMFTConfig::line_search_rho` (default `0.5`; not a separate INPUT keyword).

Pseudocode (SPG occupation inner loop; `P` denotes $P_{\mathcal{C}}$; no spin
labels):

```text
initialize n feasible; tol = rdmft_occ_proj_tol
have_prev = false
f_history = empty deque (length M = 10)
for inner = 0 .. occ_maxiter-1:
   g = grad_n(E, n)
   r = n - P(n - g)
   if ||r||_inf <= tol: break
   if have_prev:
      s = n - n_prev; y = g - g_prev
      alpha_BB = clamp((s.s)/(s.y), alpha_min, alpha_max)   // BB1 + safeguard
   else:
      alpha_BB = 1 / max(1, ||g||_inf)
   d = P(n - alpha_BB * g) - n         // spectral projected direction (descent)
   dd = g . d                          // dd <= -||d||^2 / alpha_BB <= 0
   push E(n) to f_history (drop oldest if size > M)
   f_max = max f_history
   lambda = 1
   for ls = 1 .. line_search_max_iter:
      if E(n + lambda*d) <= f_max + c1 * lambda * dd: accept; break
      lambda *= rho
   n_prev = n; g_prev = g; have_prev = true
   if accepted: n <- n + lambda * d
end
```

When not to use:

- For convex inner sub-problems with a small active set near the solution,
  an active-set / reduced-Newton solver (not currently implemented) would
  achieve faster (superlinear) local convergence.

### 6.3 Active Set Method (alias)

The `active_set` value of `rdmft_constraint` is now an **alias for
`projected_gradient`**: both route to the SPG implementation in §6.2.  The
projection $P_{\mathcal{C}}$ already handles bound activation implicitly; the
bespoke active-set bookkeeping (with its own descent safeguard, optimiser
restart on active-set change, two-stage SD-fallback line search and KKT
complementarity stop) was removed in favour of the simpler SPG formulation,
which has the same theoretical guarantees and is more robust on regularised
separable functionals (Müller / Power / GEO) where $\partial E/\partial n$
is bounded but very large near the regularisation cutoff.


---

## 7. Optimization Strategies

### 7.1 Alternating Optimization

Alternate between:

- **Occupation step**: Fix $C^{\mathbf{k}}$, optimise $n_{i\mathbf{k}}$ on the
  feasible set $\mathcal{C} = \{0 \le n_{ik} \le 1,\; \sum_k w_k\sum_i n_{ik} = N_e\}$.
  Implementation: the Spectral Projected Gradient method described in §6.2.

- **Orbital step**: Fix $n_{i\mathbf{k}}$, optimise $X^{\mathbf{k}}$ on the
  Stiefel manifold $\mathrm{St}(N_b, N_{\mathrm{basis}})$ in X-space
  (after the Cholesky $S = U^H U$ change of variable, §3.4 + §5.4).

  The orbital sub-problem uses a textbook Riemannian gradient method
  (Absil–Mahony–Sepulchre, *Optimization Algorithms on Matrix
  Manifolds*, Princeton 2008, §4.2). Per inner iteration:

  1. Evaluate $E$ and the Euclidean gradient $G$, then project to the
     Riemannian gradient
     $G_R = G - X\,\mathrm{sym}(X^H G)$ at the current iterate $X_k$.
  2. Build the search direction $D \in T_{X_k}\mathrm{St}$ according to
     `rdmft_orb_optimizer`:
       - `sd`: $D = -G_R$.
       - `cg`: Polak–Ribière⁺ with vector transport by tangent-space
         projection of $(G_R^{\mathrm{prev}}, D_{\mathrm{prev}})$ at
         $X_k$.
     Descent safeguard: if $\langle G_R, D\rangle \ge 0$ (or non-finite),
     reset $D = -G_R$.
  3. Line search along the retracted curve
     $X(\alpha) = R_{X_k}(\alpha\,D)$ where $R$ is one of three Stiefel
     retractions selected by `rdmft_orb_retraction` (§5.4; default
     `polar`): **non-monotone Strong Wolfe** (`nonmonotone_strong_wolfe_line_search`
     in `rdmft_optimizer.h`), i.e. sufficient decrease vs a reference
     $f_{\mathrm{ref}}$ from the last $M$ energies and a curvature bound on
     $|\varphi'(\alpha)|$. Initial trial uses INPUT `rdmft_alpha_step` (mapped to
     `line_search_alpha_init` in code; with CG scaling from the previous accepted
     step when available).
  4. Commit $X_{k+1} = R_{X_k}(\alpha\,D)$ with the accepted $\alpha$.

  No Barzilai–Borwein spectral step on this path (unlike SPG occupations), no
  Wen–Yin trust radius, no suspicious-descent guard. Convergence
  criterion: $\lVert G_R\rVert_F \le \varepsilon_g \max(1, \lVert G_R(X_0)\rVert_F)$
  with `rdmft_orb_grad_tol` as $\varepsilon_g$ and $X_0$ the iterate at the start
  of the orbital inner loop.

  **Caveat for regularised functionals (Müller / Power / GEO).**  Because
  the regularised functionals at fractional occupations have no global
  lower bound as a function of $X$ (the exchange Coulomb integral
  $K_{ij} = \langle \phi_i\phi_j|r_{12}^{-1}|\phi_i\phi_j\rangle$ can be
  made arbitrarily large by orbital concentration), an aggressive line search can
  in principle accept huge "descents" $E_{\mathrm{trial}} - E \sim -10^4$
  Ry into a spurious unphysical basin. Mitigations: (i) tighten
  `rdmft_orb_grad_tol` only as far as the physical basin's descent map
  warrants; (ii) reduce `rdmft_alpha_step` so the first trial
  $\alpha_0\,D$ is small relative to $\lVert X\rVert_F$; (iii) when the
  alternating outer loop diverges in this way, switch to
  `rdmft_solver_strategy = joint`, which scales the orbital block by
  `joint_orb_scale` and avoids the alternating amplification.

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
keyword (default `cg`), drives its evolution. Each outer iteration:

1. evaluates the energy $E$ and Euclidean gradients $(\nabla_n E, \nabla_C E)$;
2. chain-rules $\nabla_n E \to \nabla_p E$ through the occupation
   parameterisation, and projects $\nabla_C E$ onto the Stiefel tangent space
   at $C^{\mathbf{k}}$ to obtain $G_R^{\mathbf{k}}$;
3. packs the gradient as $g = (\nabla_p E,\,\mathrm{flat}(G_R^1), \ldots)$
   and asks the single unified optimiser for a packed descent direction
   $d = (d_p, \mathrm{flat}(d_{C^1}), \ldots)$;
4. re-projects each orbital block $d_{C^{\mathbf{k}}}$ onto the tangent space
   at the current $C^{\mathbf{k}}$ (the packed direction may drift slightly off
   the tangent space numerically) and falls back to
   the packed steepest-descent direction $d = -g$ if the joint directional
   derivative
   $dd_{\text{total}} = \langle \nabla_p E, d_p\rangle + \sum_{\mathbf{k}} \langle G_R^{\mathbf{k}}, d_{C^{\mathbf{k}}}\rangle_{S^{\mathbf{k}}}$
   is not negative;
5. performs a single **non-monotone Strong Wolfe** line search along the packed
   direction: the occupation parameters are updated linearly
   $p \leftarrow p + \alpha\, d_p$, while each $C^{\mathbf{k}}$ is retracted onto
   the generalised Stiefel manifold via
   $C^{\mathbf{k}} \leftarrow R_{C^{\mathbf{k}}}(\alpha\, d_{C^{\mathbf{k}}})$;
6. when `rdmft_joint_optimizer` is `cg`, builds the new packed gradient at the
   step's end-point and calls `update()` so the conjugate-gradient state is
   advanced; then refreshes the augmented-Lagrangian multiplier.

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

with step size $\alpha_t$ chosen by line search (here: non-monotone Strong Wolfe
along the retraction / product update, except SPG occupations which use BB + non-monotone Armijo).

### 8.2 Conjugate Gradient (CG)

Riemannian CG with vector transport $\mathcal{T}$:

$$
\eta_t = -\mathrm{grad} f(x_t) + \beta_t  \mathcal{T}_{x_{t-1} \to x_t}(\eta_{t-1})
$$

with Fletcher-Reeves or Polak-Ribière $\beta_t$.

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


