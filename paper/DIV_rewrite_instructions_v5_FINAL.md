# D-IV Paper Rewrite Instructions: Final Version (v5)

## Overview

Rewrite the D-IV paper with:
- **Pointwise-then-OLS** as the main estimator (closest to Fréchet regression)
- **Endpoint projection** as a complement (tighter finite-sample bound, global feasibility)
- Both asymptotically equivalent under the strictness assumption

Key references:
- `fivr_v4.tex` (current draft) — contains ALL ℓ∞-Hadamard differentiability machinery, covariance continuity, projection properties, and existing proofs. REUSE directly.
- R3D paper (Van Dijcke 2025, arxiv:2504.03992) — asymptotic theory structure
- `hadamard_proof_v5.tex` — coefficient projection CLT (supplementary appendix only)

---

## Part 1: The Narrative Arc

**Population (Section 3.1):** The Fréchet IV problem $\argmin_{m \in \mathcal{P}} \mathbb{E}[s(Z,x) W_2^2(Y, m)]$ has solution $\Pi_\mathcal{Q}(\psi_x)$. Under correct specification, $\psi_x \in \mathcal{Q}$ so the projection is inactive. The structural quantile function is the Fréchet IV barycenter, a valid QF, and linear in $x$ simultaneously.

**Finite-sample tension:** In finite samples, $\hat\psi_x$ may violate monotonicity. PAVA is nonlinear, so one cannot simultaneously maintain the barycenter interpretation at every $x$, linearity in $x$, and global monotonicity.

**The D-IV estimator (Section 3.3):** We solve the Fréchet IV problem at every observed covariate value, obtaining the IV-weighted Wasserstein barycenters $\Pi_\mathcal{Q}(\hat\psi_{X_j})$, and recover the linear parametrization by OLS. This preserves the Fréchet interpretation at every evaluation point. However, the recovered coefficients need not produce monotone QFs at all $x$ values in the support.

**Endpoint reparametrization (Section 4.1):** For finite-sample improvement guarantees, we also consider the endpoint reparametrization: project at the two support extremes and recover coefficients. This produces a linear model with valid QFs everywhere. The parallelogram identity gives a weighted improvement guarantee that is strictly tighter than the pointwise oracle inequality (by Popoviciu's inequality).

---

## Part 2: What Changes and What Stays

### Unchanged:
- Sections 1-2 (Introduction, Setup, Model)
- Section 3.1 (Population problem, Fréchet IV identification, Lemma 1)
- Section 3.2 (Identification interpretation, decomposition, CLP/MD comparison)
- All appendix results in fivr_v4.tex

### Major rewrites:
- Section 3.3 (Estimator definition) → pointwise-then-OLS as main, endpoint as complement
- Section 4.1 (Finite-sample properties) → both oracle inequalities + parallelogram
- Section 4.2 (Asymptotic distribution) → restructured
- Section 4.3 (Inference) → simplified

---

## Part 3: The Main Estimator (Section 3.3)

### 3.1 The D-IV estimator (pointwise-then-OLS)

Present as the primary estimator in four steps.

**Step 1: Estimate the IV weights.** For each observation $j$ and target $x$:
```
ŝⱼ(x) = 1 + (x - μ̂_X)ᵀ (Σ̂_{ZX}ᵀ Σ̂_{ZZ}⁻¹ Σ̂_{ZX})⁻¹ Σ̂_{ZX}ᵀ Σ̂_{ZZ}⁻¹ (Zⱼ - μ̂_Z)
```

**Step 2: Compute the IV-weighted average quantile curve at each observation:**
```
ψ̂_{Xⱼ}(u) = (1/n) Σᵢ ŝᵢ(Xⱼ) Q_{Yᵢ}(u)
```

**Step 3: Solve the Fréchet IV problem at each observation.** Project onto the space of valid quantile functions:
```
Q̂(Xⱼ, u) = Π_Q(ψ̂_{Xⱼ})(u)
```
Each $\hat Q(X_j, \cdot)$ is the closest valid probability distribution (in $W_2$) to the IV-weighted average at covariate value $X_j$ — the sample Fréchet IV barycenter.

**Step 4: Recover coefficient functions by OLS:**
```
β̂₁(u) = Σ̂_{XX}⁻¹ · (1/n) Σⱼ (Xⱼ - μ̂_X) Q̂(Xⱼ, u)
β̂₀(u) = (1/n) Σⱼ Q̂(Xⱼ, u)
```

**Key points to emphasize:**
- Steps 1-3 are the sample analogue of the population Fréchet IV problem (Eq. [projection-fivr_intercept_main]). At each $X_j$, the estimator finds the closest valid distribution to the IV-weighted average.
- Step 4 recovers the linear parametrization from the projected evaluations.
- The projection step guarantees each $\hat Q(X_j, \cdot)$ is a valid QF and closer to its target than $\hat\psi_{X_j}$ in all $L^p$ norms (Lemma [improvement]).
- The recovered coefficients $(\hat\beta_0, \hat\beta_1)$ need not lie in $\mathcal{C}$: the implied quantile function $\hat q(x, \cdot) = \hat\beta_0 + \hat\beta_1(x - \mu_X)$ may violate monotonicity at covariate values other than $\{X_1, \ldots, X_n\}$. For applications requiring global feasibility, see the endpoint reparametrization (Section 4.1).

### 3.2 Multivariate X via reference distributions

For $X = (D, W)$ with scalar endogenous treatment $D$ and additional covariates $W \in \mathbb{R}^{p-1}$.

Under the linear model:
```
Q_{Yⱼ}(u) = β₀(u) + β_D(u)(Dⱼ - μ_D) + γ(u)ᵀ(Wⱼ - μ_W) + ηⱼ(u)
```

Choose a common reference distribution $H$ for $W$ (natural: empirical marginal of $W$). Define endpoint-averaged quantile curves:
```
q₊ᴴ(u) = ∫ q(d₊, w, u) dH(w) = β₀(u) + β_D(u)d₊ + γ(u)ᵀμ_H
q₋ᴴ(u) = ∫ q(d₋, w, u) dH(w) = β₀(u) + β_D(u)d₋ + γ(u)ᵀμ_H
```

Key identity (since $H$ is the SAME at both endpoints):
```
β_D(u) = [q₊ᴴ(u) - q₋ᴴ(u)] / (d₊ - d₋)
```

Estimation: form $\hat\gamma_\pm^H = \tilde\beta_0 + \tilde\beta_D d_\pm + \tilde\gamma^\top \bar W$, project each, recover $\hat\beta_D = [\hat\gamma_+^{H*} - \hat\gamma_-^{H*}] / (d_+ - d_-)$.

---

## Part 4: Finite-Sample Properties (Section 4.1)

### 4.1 Pointwise improvement (existing Lemma)

Keep: $\|\Pi_\mathcal{Q}(\hat Q) - Q_0\|_{L^p} \le \|\hat Q - Q_0\|_{L^p}$ for any monotone target $Q_0$, all $p$.

### 4.2 Pointwise-then-OLS oracle inequality

**Theorem (Distributional improvement).** For any $b \in \mathcal{C}_n$ and any realization:
```
‖β̂₀ᵖʷ - b₀‖²_{L²} + σ̂²_X ‖β̂₁ᵖʷ - b₁‖²_{L²} ≤ ‖β̃₀ - b₀‖²_{L²} + σ̂²_X ‖β̃₁ - b₁‖²_{L²}
```
Under correct specification, $\beta \in \mathcal{C}$ a.s., so taking $b = \beta$ gives improvement toward the truth. Under misspecification, taking $b = \Pi_{\mathcal{C}_n}(\beta^{\text{unc}})$ gives improvement toward the best constrained approximation.

**Proof:** Three steps.
1. At each $j$, Lemma (improvement) with target $q_b(X_j, \cdot) \in \mathcal{Q}$ (valid QF since $b \in \mathcal{C}_n$):
   $\int |\Pi_\mathcal{Q}(\hat\psi_{X_j})(u) - q_b(X_j, u)|^2 du \le \int |\hat\psi_{X_j}(u) - q_b(X_j, u)|^2 du$
2. OLS Pythagorean decomposition: regressing projected curves on $(1, \tilde X_j)$ decomposes the sample-averaged squared error into fitted + residual. Setting the comparison point to $b$ and dropping the non-negative residual term gives the inequality after integrating over $u$.
3. Identify the LHS as $\|\hat\beta_0^{pw} - b_0\|^2 + \hat\sigma_X^2\|\hat\beta_1^{pw} - b_1\|^2$ and the RHS similarly.

**Remark (Fréchet interpretation).** By the identity connecting the weighted coefficient norm to expected Wasserstein distance:
```
(1/n) Σⱼ W₂²(Q̂(Xⱼ, ·), q_b(Xⱼ, ·)) ≤ (1/n) Σⱼ W₂²(ψ̂_{Xⱼ}, q_b(Xⱼ, ·))
```
The projected Fréchet barycenters are on average closer to the target distributions in Wasserstein distance.

### 4.3 Endpoint reparametrization and parallelogram identity

**Definition (Endpoint estimator).** Define the support endpoints $x̲, x̄$ (population or sample; see Remark below). Form endpoint curves $\hat\gamma_\pm = \tilde\beta_0 \pm \delta\tilde\beta_1$ where $\delta = (x̄ - x̲)/2$. Project each: $\hat\gamma_\pm^* = \Pi_\mathcal{Q}(\hat\gamma_\pm)$. Recover: $\hat\beta_1^{\text{end}} = (\hat\gamma_+^* - \hat\gamma_-^*) / (2\delta)$.

**Properties:**
- The implied QF $\hat q^{\text{end}}(x, \cdot)$ is a valid QF for all $x \in [x̲, x̄]$ (by vertex characterization).
- $(\hat\beta_0^{\text{end}}, \hat\beta_1^{\text{end}}) \in \mathcal{C}$: global feasibility.

**Theorem (Endpoint improvement / parallelogram identity).** Suppose the population endpoint curves $\gamma_+, \gamma_- \in \mathcal{Q}$. (Holds under correct specification; more generally, whenever the pseudo-true IV-weighted QFs at the endpoints are non-decreasing.) Then for any realization:
```
‖β̂₀ᵉⁿᵈ - β₀‖²_{L²} + δ² ‖β̂₁ᵉⁿᵈ - β₁‖²_{L²} ≤ ‖β̃₀ - β₀‖²_{L²} + δ² ‖β̃₁ - β₁‖²_{L²}
```
where $\delta = (x̄ - x̲)/2$.

**Proof:** Parallelogram identity applied to endpoint errors $e_\pm = \varepsilon_0 \pm \delta\varepsilon_1$. [Same proof as in v4 instructions.]

**Remark (Popoviciu comparison).** By Popoviciu's variance inequality, $\sigma_X^2 \le \delta^2 = ((x̄ - x̲)/2)^2$ with equality iff $X$ is binary. The endpoint bound has weight $\delta^2 \ge \sigma_X^2$ on the slope, so it is always at least as tight as the pointwise oracle inequality on the slope component. For uniform $X$ on $[x̲, x̄]$, $\sigma_X^2 = \delta^2/3$: the pointwise leakage is three times worse.

**Remark (What we can and cannot prove).** Both the pointwise and endpoint improvement results guarantee a weighted joint improvement in $(\beta_0, \beta_1)$. Neither guarantees improvement of $\beta_1$ alone in unweighted $L^2$. The leakage from the intercept into the slope bound is dampened by larger endpoint separation (or larger $\sigma_X^2$) but is of the same order $O_P(n^{-1})$ as the slope MSE. In all simulations, both $\beta_0$ and $\beta_1$ improve individually.

**Remark (Two estimators, one asymptotic distribution).** The pointwise-then-OLS and endpoint estimators are distinct in finite samples but asymptotically equivalent under Assumption 5 (Section 4.2). In simulations for scalar $X$, their MSE differs by less than 0.5 percentage points. The pointwise estimator preserves the Fréchet barycenter interpretation at every evaluation point. The endpoint estimator provides global feasibility and a tighter finite-sample bound. Both are valid and complementary.

---

## Part 5: Asymptotic Distribution (Section 4.2)

### 5.1 Existing machinery (ALL from fivr_v4.tex)

- **Theorem `thm:convergence_unprojected`**: CLT for $\hat\psi_x$ at fixed $x$
- **Theorem `thm:convergence_projected`**: CLT for $\Pi_\mathcal{Q}(\hat\psi_x)$ at fixed $x$
- Supporting lemmas: `lem:projection_properties`, `lem:Gamma-cont-strong`, `cor:cont-version`, `lem:iso-hadamard`

### 5.2 Uniform convergence of the pointwise-then-OLS estimator

**Theorem (Asymptotic distribution of the pointwise D-IV estimator).** Under Assumptions 1-4 and Average Strictness:
```
√n(β̂ᵖʷ(·) - β*(·)) ⟹ G_β(·)   in ℓ∞([a,b])^{p+1}
```
the same Gaussian process as the unconstrained 2SLS estimator.

**Proof:** This proof is already essentially in fivr_v4.tex (the proof of Theorem `thm:convergence_beta`). The key structure:

**Step 1 (Decomposition).** Write $\hat\beta_1^{pw}(u) = \tilde\beta_1(u) + \Delta_{1,n}(u)$ where $\Delta_{1,n}(u) = \hat\Sigma_{XX}^{-1} (1/n) \sum_j \tilde X_j D_{X_j}(u)$ and $D_x(u) = \Pi_\mathcal{Q}(\hat\psi_x)(u) - \hat\psi_x(u)$ is the PAVA correction.

**Step 2 (The correction is negligible in $\ell^\infty$).** Under Assumption 5, $\psi_x$ is strictly increasing for each $x \in [x̲, x̄]$. By Lemma `lem:iso-hadamard`, $\Pi_\mathcal{Q}$ is Hadamard differentiable at $\psi_x$ with derivative equal to the identity. Therefore the Hadamard remainder satisfies:
```
‖D_x‖_∞ = o(‖ψ̂_x - ψ_x‖_∞) = o_P(n^{-1/2}) · (1 + ‖x - μ_X‖)
```

**Critical point about uniformity over $x$:** This remainder bound uses the Hadamard differentiability at each fixed $x$, not uniformity over $x$. However, it does NOT require a uniform-in-$x$ Hadamard differentiability result. The reason: the OLS step AVERAGES the corrections $D_{X_j}$ weighted by $\tilde X_j$. We need:
```
√n ‖Δ_{1,n}‖_∞ = √n · ‖Σ̂_{XX}⁻¹ (1/n) Σⱼ X̃ⱼ D_{X_j}‖_∞ = o_P(1)
```

This follows from:
- $\|D_{X_j}\|_\infty = o_P(n^{-1/2}) \cdot (1 + \|X_j - \mu_X\|)$ at each $j$ (Hadamard remainder at fixed $X_j$)
- $\hat\Sigma_{XX}^{-1} = O_P(1)$
- $(1/n) \sum_j \|X_j - \hat\mu_X\| \cdot (1 + \|X_j - \mu_X\|) = O_P(1)$ (by finite 4th moments, Assumption 4)

Combining: $\sqrt{n}\|\Delta_{1,n}\|_\infty \le O_P(1) \cdot o_P(1) \cdot O_P(1) = o_P(1)$.

**The $o_P(n^{-1/2})$ rate of $D_{X_j}$ at each $X_j$:** This is the key step that might seem to require uniformity but does not. For each fixed $x$ in the interior of the support, the Hadamard differentiability at $\psi_x$ gives $\|D_x\|_\infty = o(\|\hat\psi_x - \psi_x\|_\infty)$. The "little-o" rate may depend on $x$ through $\psi_x$ — but since $x \mapsto \psi_x$ is continuous and $[x̲, x̄]$ is compact, the strict increase constant $\kappa$ is uniform over $x$ (it's a minimum of $\kappa(x)$ over a compact set), and the Hadamard remainder rate is controlled by $\kappa$ and the modulus of continuity of $\psi_x'$, both of which are bounded on the compact support. So the "little-o" is in fact uniform over $x \in [x̲, x̄]$ under the maintained assumptions — but this follows from compactness and continuity, not from a separate uniform-Hadamard theorem.

**Step 3 (Slutsky).** $\sqrt{n}(\hat\beta^{pw} - \beta^*) = \sqrt{n}(\tilde\beta - \beta^*) + o_P(1)$ in $\ell^\infty([a,b])^{p+1}$. The unconstrained CLT (Theorem `thm:convergence_beta`) gives the limit.

### 5.3 Uniform convergence of the endpoint estimator

**Theorem (Asymptotic distribution of the endpoint D-IV estimator).** Under the same assumptions:
```
√n(β̂ᵉⁿᵈ(·) - β*(·)) ⟹ G_β(·)   in ℓ∞([a,b])^{p+1}
```
the same Gaussian process.

**Proof (with population endpoints $x̲, x̄$):**

**Step 1 (Joint convergence of unprojected endpoint processes).** Both $\hat\psi_{x̲}$ and $\hat\psi_{x̄}$ are continuous linear functionals of the same empirical process. Joint convergence:
```
(√n(ψ̂_{x̲} - ψ_{x̲}), √n(ψ̂_{x̄} - ψ_{x̄})) ⟹ (G₋, G₊)   in ℓ∞([a,b])²
```
follows from the Donsker theorem applied to the union of function classes (Example 2.10.7 of vdVW).

**Step 2 (Joint convergence of projected processes).** $\Pi_\mathcal{Q}$ is Hadamard differentiable at $\psi_{x̲}$ and $\psi_{x̄}$ tangentially to $C([a,b])$ with derivative the identity (Lemma `lem:iso-hadamard`). The joint limit $(G_-, G_+)$ has a.s. continuous paths (Corollary `cor:cont-version`). Extended functional delta method:
```
(√n(Π_Q(ψ̂_{x̲}) - ψ_{x̲}), √n(Π_Q(ψ̂_{x̄}) - ψ_{x̄})) ⟹ (G₋, G₊)
```

**Step 3 (Continuous mapping).** Coefficient recovery $(f_-, f_+) \mapsto ((a_+ f_- - a_- f_+)/(a_+-a_-), (f_+ - f_-)/(a_+-a_-))$ is continuous linear $\ell^\infty \times \ell^\infty \to \ell^\infty \times \ell^\infty$. CMT gives $\sqrt{n}(\hat\beta^{end} - \beta^*) \Rightarrow G_\beta$.

**Step 4 (Identification with unconstrained limit).** The same linear map applied to the UNPROJECTED processes gives $\sqrt{n}(\tilde\beta - \beta^*)$. Since the Hadamard derivative is the identity at both endpoints, the image of $(G_-, G_+)$ under the recovery map is the same in both cases. SHOW this explicitly — write out that the recovery map is the same, the inputs have the same limit, so the outputs have the same limit.

**Key advantage: no uniformity over $x$ needed.** The endpoint CLT requires only the existing fixed-$x$ CLT at two points, joint convergence (union of Donsker classes), and a continuous linear map. No uniformity over a continuum of $x$ values is needed.

### 5.4 Sample endpoints perturbation argument

**Proposition (Sample endpoints are asymptotically equivalent).** If sample endpoints $x̲_n = \min_i X_i$ and $x̄_n = \max_i X_i$ replace population endpoints, the asymptotic distribution is unchanged.

**Proof:** Write:
```
√n[Π_Q(ψ̂_{x̄_n}) - ψ_{x̄_n}] = √n[Π_Q(ψ̂_{x̄}) - ψ_{x̄}] + Rₙ
```
where the remainder is:
```
Rₙ = √n[Π_Q(ψ̂_{x̄_n}) - Π_Q(ψ̂_{x̄})] - √n[ψ_{x̄_n} - ψ_{x̄}]
```

**Projection term:** By the $\ell^\infty$-Lipschitz property of $\Pi_\mathcal{Q}$:
```
‖Π_Q(ψ̂_{x̄_n}) - Π_Q(ψ̂_{x̄})‖_∞ ≤ ‖ψ̂_{x̄_n} - ψ̂_{x̄}‖_∞ = |x̄_n - x̄| · ‖β̃₁‖_∞
```

**Population term:** $\|\psi_{x̄_n} - \psi_{x̄}\|_\infty = |x̄_n - x̄| \cdot \|\beta_1\|_\infty$.

**Rate:** For compact support with $f_X(x̄) > 0$: $|x̄_n - x̄| = O_P(n^{-1})$. And $\|\tilde\beta_1\|_\infty = O_P(1)$. Therefore:
```
‖Rₙ‖_∞ ≤ √n · O_P(n⁻¹) · O_P(1) = O_P(n⁻¹/²) = o_P(1)
```

Same for $x̲_n$. The denominator $x̄_n - x̲_n \to x̄ - x̲$ a.s. By Slutsky in $\ell^\infty([a,b])$:
```
√n(β̂₁ᵉⁿᵈ(·; x̲_n, x̄_n) - β₁(·)) ⟹ (G₊ - G₋) / (x̄ - x̲)
```

**Assumption required:** Compact support with $f_X > 0$ near $x̲$ and $x̄$ (so sample extrema converge at rate $n^{-1}$). Standard in grouped-data IV settings. For unbounded support, treat endpoints as known or add a convergence rate assumption.

NOTE: This argument uses only the $\ell^\infty$-Lipschitz property of $\Pi_\mathcal{Q}$ and the $O_P(n^{-1})$ rate of sample extrema. No uniformity over $x$ is needed — we are just saying "moving the evaluation point by $O(n^{-1})$ moves the projected curve by $O(n^{-1})$ in sup-norm."

### 5.5 Discussion of Assumption 5

Keep existing discussion. Add: for scalar $X$, $\psi_x$ is a convex combination of $\psi_{x̲}$ and $\psi_{x̄}$ for $x \in [x̲, x̄]$, so strict increase at endpoints implies strict increase everywhere.

---

## Part 6: Inference (Section 4.3)

### 6.1 Multiplier bootstrap

Unchanged. Operates on unconstrained influence functions. No recomputation of the projection.

### 6.2 Bootstrap validity

Explicit proof:
1. Bootstrap validity for the unconstrained process.
2. $\sqrt{n}\|\hat\beta - \tilde\beta\|_\infty = o_P(1)$ (from Section 5.2 or 5.3, depending on which estimator).
3. Bounded-Lipschitz argument.
4. Bootstrap critical values consistent for all three estimators (unconstrained, pointwise, endpoint).

---

## Part 7: What NOT to Include

1. **Do NOT present the QP-based coefficient projection as a separate main estimator.** Mention briefly in a remark: it enforces the same constraint set $\mathcal{C}$ with a different objective ($\Sigma_{XX}$-weighted $L^2$), and its continuum CLT requires substantially stronger assumptions (see Supplementary Appendix).

2. **Do NOT claim equivalence between endpoint and coefficient projection.** They enforce the same constraints but minimize different objectives. They coincide for binary treatment only.

3. **Do NOT claim slope-only improvement.** Both the pointwise and endpoint bounds give weighted joint improvement.

4. **Do NOT include shrinkage-toward-zero or monotone-$\beta_1$ results.** Tangential.

5. **Do NOT include tent-bound/modulus/anchor machinery in main text.** Supplementary only.

---

## Part 8: Appendix Structure

**A.1:** Identification proofs (keep from fivr_v4.tex)

**A.2:** Projection properties + ℓ∞-Hadamard differentiability (keep from fivr_v4.tex)

**A.3:** Covariance continuity + continuous modification (keep from fivr_v4.tex)

**A.4:** Unprojected CLT proof (keep from fivr_v4.tex)

**A.5:** Projected CLT at fixed $x$ (keep from fivr_v4.tex)

**A.6:** NEW: Asymptotic equivalence of the pointwise-then-OLS estimator (the $\Delta_{1,n} = o_P(n^{-1/2})$ argument from Section 5.2, with the compactness/uniformity discussion)

**A.7:** NEW: Joint convergence of endpoint processes + continuous mapping to coefficients (Section 5.3)

**A.8:** NEW: Sample endpoints perturbation (Section 5.4)

**A.9:** NEW: Parallelogram identity proof + pointwise oracle inequality proof

**A.10:** Bootstrap validity proof (updated)

**A.11:** Covariance derivation (keep from fivr_v4.tex)

**A.12 (Supplementary):** Coefficient projection CLT in the continuum (from hadamard_proof_v5.tex)

---

## Part 9: Simulation Updates

Scalar code `estimate_coef_proj_scalar` = endpoint estimator. `estimate_pw_fivr_scalar` = pointwise-then-OLS.

Updates:
1. Present pointwise-then-OLS as "D-IV" (main estimator) in tables.
2. Present endpoint as "D-IV (endpoint)" as additional comparison.
3. Add parallelogram-weighted norm as additional metric.
4. Existing results show both are nearly identical for scalar $X$ — report this.
5. For multivariate case, update pointwise code to reference distribution approach.

---

## Part 10: Checklist

- [ ] Theorem `thm:convergence_projected` applies at fixed $x$ including $x̲, x̄$ — verify
- [ ] Joint convergence of $(ψ̂_{x̲}, ψ̂_{x̄})$ via union of Donsker classes — write explicitly
- [ ] The "little-o is uniform over compact support" argument (Section 5.2 Step 2) — verify that $\kappa$ uniform + compact support suffices
- [ ] Parallelogram identity algebra — verify
- [ ] Sample endpoints perturbation — verify $O_P(n^{-1})$ rate and $o_P(1)$ in ℓ∞
- [ ] Step 4 of endpoint CLT — SHOW identification with unconstrained limit explicitly
- [ ] Pointwise oracle inequality proof (3 steps) — verify OLS Pythagorean decomposition
- [ ] Multivariate reference distribution — verify consistency with 2SLS
- [ ] Bootstrap proof uses bounded-Lipschitz explicitly
- [ ] Honest remark about individual components included
- [ ] "No uniformity over $x$ needed" advantage of endpoint approach stated
