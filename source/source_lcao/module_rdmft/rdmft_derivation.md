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

Perform gradient descent on $n_{i\mathbf{k}}$ and project back onto the feasible set:

$$
n_{i\mathbf{k}}^{(t+1)} = \mathrm{Proj}_{[0,1]}\bigl(n_{i\mathbf{k}}^{(t)} - \alpha_t \frac{\partial E}{\partial n_{i\mathbf{k}}}\bigr)
$$

followed by rescaling to enforce the electron number constraint:

$$
n_{i\mathbf{k}}^{(t+1)} \leftarrow n_{i\mathbf{k}}^{(t+1)} \cdot \frac{N_e}{\sum_{\mathbf{k}} w_{\mathbf{k}} \sum_i n_{i\mathbf{k}}^{(t+1)}}
$$

(with re-clipping if rescaling violates box constraints).

### 6.3 Active Set Method

Maintain active sets $\mathcal{A}_0 = \{(i,\mathbf{k}): n_{i\mathbf{k}} = 0\}$ and
$\mathcal{A}_1 = \{(i,\mathbf{k}): n_{i\mathbf{k}} = 1\}$.

At each step:

1. On the free set $\mathcal{F} = \{(i,\mathbf{k}): 0 < n_{i\mathbf{k}} < 1\}$, solve the
  reduced problem with equality constraint.
2. Check KKT multipliers for active constraints; release violated ones.
3. Check feasibility for the updated free variables; add newly violated constraints.

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
keyword (default L-BFGS), drives its evolution. Each outer iteration:

1. evaluates the energy $E$ and Euclidean gradients $(\nabla_n E, \nabla_C E)$;
2. chain-rules $\nabla_n E \to \nabla_p E$ through the occupation
   parameterisation, and projects $\nabla_C E$ onto the Stiefel tangent space
   at $C^{\mathbf{k}}$ to obtain $G_R^{\mathbf{k}}$;
3. packs the gradient as $g = (\nabla_p E,\,\mathrm{flat}(G_R^1), \ldots)$
   and asks the single unified optimiser for a packed descent direction
   $d = (d_p, \mathrm{flat}(d_{C^1}), \ldots)$;
4. re-projects each orbital block $d_{C^{\mathbf{k}}}$ onto the tangent space
   at the current $C^{\mathbf{k}}$ (necessary because Euclidean preconditioners
   used by L-BFGS / Adam generally leave the tangent space) and falls back to
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
   history (L-BFGS $(s,y)$ pairs, Adam moments, CG previous gradient) is
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

### 8.3 L-BFGS

Limited-memory BFGS adapted to Riemannian setting using vector transport
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

## 11. Summary of Implementation Requirements


| Component                   | Variables | Manifold               | Gradient             |
| --------------------------- | --------- | ---------------------- | -------------------- |
| One-body energy             | n, C      | —                      | Eqs. in §5.1, §5.2   |
| Hartree energy              | n, C      | —                      | Eqs. in §5.1, §5.2   |
| XC energy (HF/Müller/Power) | n, C      | —                      | Eqs. in §5.1, §5.2   |
| XC energy (GU)              | n, C      | —                      | Additional SIC terms |
| Occupation constraint       | n         | $[0,1]$, $\sum w n = N_e$ | §4, §6               |
| Orbital orthogonality       | C         | Stiefel                | §5.3, §5.4           |


