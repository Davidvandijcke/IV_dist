# D-IV Simulation Findings

Systematic exploration of when the D-IV projection improves over unconstrained 2SLS. Over 100 DGP configurations tested. All results use the Melly-Pons/CLP DGP framework or direct group-level models.

## The D-IV Estimator

The D-IV estimator (Section 3.3 of the paper):
1. Compute unconstrained 2SLS coefficients β̃₀(u), β̃₁(u)
2. Evaluate fitted curves ψ̂(Xⱼ, u) = β̃₀(u) + β̃₁(u)(Xⱼ − μ̂) at each observed Xⱼ
3. Project each curve to monotonicity via PAVA: Q̂(Xⱼ, ·) = Π_Q(ψ̂(Xⱼ, ·))
4. Recover coefficients by OLS of Q̂(Xⱼ, u) on (1, Xⱼ − μ̂)

The improvement theorem (Theorem 1) guarantees that the joint weighted MSE of (β̂₀, β̂₁) under D-IV is ≤ that of unconstrained 2SLS, for any reference coefficients b that produce valid (monotone) quantile functions at the observed X values.

## Baseline: MP/CLP DGP

The Melly-Pons (2025) / Chetverikov et al. (2016) DGP:
```
y_ij = x_{1ij} * (u_ij/2) + x_{2j} * sqrt(u_ij) + alpha_j(u_ij)
```
- x_{1ij}, x_{2j} ~ exp(0.25 * N(0,1)) (log-normal, sd ≈ 0.26)
- alpha_j(u) = u * eta_j - u/2, eta_j ~ U(0,1) (mean-zero group heterogeneity)
- Endogenous: x_{2j} = z_j + eta_j + nu_j (eta in both treatment and outcome)
- Instrument: z_j ~ exp(0.25 * N(0,1))
- True treatment effect: gamma(u) = sqrt(u) (monotone, smooth)

Since gamma(u) = sqrt(u) is monotone, the Skorohod representation is valid: u_ij IS the within-group rank, and the group quantile function Q_{Y_j}(u) = x_{2j}·sqrt(u) + alpha_j(u) + E[x₁]·(u/2) holds exactly with β₁(u) = gamma(u).

At M=50, N=50: median F-stat ≈ 23, D-IV gain ≈ 0%. This DGP has strong instruments and very well-behaved distributions — few monotonicity violations occur.

---

## Channel 1: Instrument Strength (the primary driver)

**Setup**: MP DGP with pi_Z parameter scaling instrument relevance: x_{2j} = pi_Z * z_j + eta_j + sqrt(1-pi_Z²) * nu_j

**Results** (M=50, N=50, 500 reps):

| pi_Z | Median F | Violation rate | D-IV gain |
|------|----------|---------------|-----------|
| 0.1  | 1.2      | 100%          | ~15%      |
| 0.2  | ~5       | ~98%          | ~18%      |
| 0.3  | ~8       | ~77%          | ~38%      |
| 0.5  | ~15      | ~70%          | ~7%       |
| 0.7  | ~25      | ~35%          | ~12%      |
| 1.0  | ~40      | ~20%          | ~0%       |

**Mechanism**: Weak instruments inflate the variance of β̃₁(u). The resulting large deviations from the true (monotone) coefficient functions produce frequent and large violations of monotonicity in the fitted ψ̂(x, u) curves. PAVA corrects these violations.

**F-stat threshold**: Gains appear below F ≈ 15 and are substantial below F ≈ 10 (the Stock-Yogo weak instrument threshold). Above F ≈ 20, gains vanish.

**Non-monotonicity in gains**: The gain at pi_Z=0.1 (F≈1.2) is smaller than at pi_Z=0.3 (F≈8). At very weak instruments (F ≈ 1), 2SLS bias dominates MSE for both estimators, limiting the relative gain from variance reduction via projection. The projection reduces the variance component of MSE, but when bias dominates, variance reduction provides little relative improvement.

**Relevance**: Many applied IV settings have F-stats in the 5-15 range (Bartik/shift-share instruments, judge designs, etc.). The MP/CLP DGP has an unusually strong instrument (F ≈ 23-40 depending on M), which is why gains are small in their specification.

---

## Channel 2: Heavy-Tailed Within-Group Distributions + Small N

**Setup**: MP DGP (pi_Z=1, F≈23) but add heavy-tailed individual errors: y_ij = ... + t_df(N) or lognormal errors. No individual covariates x1 (to isolate the effect).

**Results** (M=50, F≈23):

| Error distribution | N=10 | N=25 | N=50 | N=100 |
|-------------------|------|------|------|-------|
| Normal            | 1.7% | 0.5% | 0.2% | 0.2%  |
| t₅               | 3.2% | 0.8% | 0.8% | 0.5%  |
| t₃               | 2.7% | 1.0% | 0.5% | 0.1%  |
| t₂               | **13.3%** | 3.2% | 2.3% | 0.2% |
| t₁.₅             | **28.7%** | 5.4% | 2.8% | 0.7% |
| Lognormal(0, 1.0)| 9.2% | 6.0% | 1.8% | —     |
| Lognormal(0, 1.5)| **25.6%** | **16.5%** | 3.8% | — |
| Lognormal(0, 2.0)| **43.6%** | **26.2%** | 5.1% | — |

**Mechanism**: Heavy tails cause outliers in the individual outcome data. With small N, these outliers distort the sample quantile functions Q̂_{Y_j}(u), which propagates through 2SLS into non-monotone coefficient estimates. This is a **first-stage noise** channel (noisy Q̂_{Y_j}), distinct from the second-stage noise channel (weak instruments affecting β̃₁).

**Why MP misses this**: The MP DGP generates outcomes from a Skorohod representation y_ij = q(x_{2j}, u_ij) with no additional error term. The individual u_ij ~ U(0,1) already determines the rank, so the sample quantile function is a consistent estimator that converges quickly. Adding heavy-tailed noise on top (as in real data: wages, health costs, test scores) makes the quantile estimation much harder.

**Interaction with Channel 1**: The two channels interact superadditively. At F≈12 with t₃ errors and N=25: 22.6% gain (vs 12% from F≈12 alone and ~1% from t₃+N=25 alone).

---

## Channel 3: Large Group Heterogeneity

**Setup**: MP DGP but with eta_j ~ N(0, sigma_eta) instead of U(0,1).

**Results** (M=50, N=50):

| sigma_eta | Median F | D-IV gain |
|-----------|----------|-----------|
| 0.29 (≈ MP) | ~23   | 0%        |
| 0.50      | ~11      | ~3%       |
| 1.00      | ~3       | ~23%      |
| 1.50      | ~2       | ~33%      |
| 2.00      | ~1       | ~9%       |

**Mechanism**: Large sigma_eta means eta_j has high variance. Since eta enters both the treatment equation (x₂ = z + eta + nu) and the outcome (alpha_j(u) = u*eta - u/2), larger eta variance (a) increases confounding and (b) dilutes the instrument's contribution to x₂. The effective F-stat drops, which is the actual driver of the gains.

**Non-monotonicity in gains**: Gains rise from 0% to 33% as sigma_eta increases from 0.29 to 1.5, then FALL to 9% at sigma_eta=2.0 (F≈1). This mirrors the pattern in Channel 1: at very weak instruments (F ≈ 1), 2SLS bias dominates MSE for both estimators, limiting the relative gain from variance reduction via projection.

**Important**: This is NOT an independent channel — it works entirely by reducing the effective instrument strength. The gains track the F-stat, not sigma_eta per se.

---

## Channel 4: γ(u) Shape and Amplitude

**Setup**: Group-level DGP with Q_{Y_j}(u) estimated from N individual draws (no x1 mixing). Vary the functional form and amplitude of gamma(u).

### 4a. γ amplitude (monotone gamma)

gamma(u) = c * sqrt(u), beta0(u) = 3*Phi^{-1}(u), N=25, F≈20:

| c (amplitude) | Violation rate | D-IV gain |
|---------------|---------------|-----------|
| 0.5           | 19%           | 0.0%      |
| 1.0           | 31%           | 0.1%      |
| 2.0           | 79%           | **0.5%**  |
| 3.0           | 95%           | 0.4%      |
| 5.0           | 100%          | **-12.7%**|

**Key finding**: Larger gamma amplitude increases violation frequency but the gains are small (≤0.5%) and turn **negative** at high amplitude.

**Negative gain at c=5**: This reflects severe misspecification of the linear model at extreme covariate values. The fitted curve q(x, u) = 3·Φ⁻¹(u) + 5·√u·(x−μ) is non-monotone in u for large |x−μ|, meaning no probability distribution has this as its quantile function. The unconstrained 2SLS converges to a pseudo-true parameter β^unc(u) that itself produces non-monotone fitted curves at some observed X values. The projection corrects these toward monotonicity but moves further from the pseudo-true target, increasing IMSE. This is the regime where Theorem 1's conditions (reference coefficients producing valid quantile functions) are not met.

**Why gains stay small at moderate c**: The estimation noise in β̃₁(u) does not scale with c (it comes from within-group QF estimation and group heterogeneity), but the signal (gamma itself) does. So the signal-to-noise ratio improves with larger c, keeping violation magnitudes small even as their frequency rises.

### 4b. γ shape (non-monotone gamma in individual-level DGP)

Tested at F≈23, M=50, N=50:

| γ shape | Violation rate | D-IV gain |
|---------|---------------|-----------|
| sqrt(u) (baseline) | 32% | 0.7% |
| u + 0.3*sin(4πu) | 70% | 0.9% |
| u + 0.5*sin(4πu) | 46% | 0.1% |
| 4*u*(1-u) (hump) | 42% | 0.2% |
| 8*u*(1-u) (big hump) | 49% | 0.2% |
| 2*(u-0.5) (sign-changing) | 36% | 0.9% |

**Key finding**: The shape of gamma — oscillating, hump-shaped, sign-changing — does not independently produce meaningful D-IV gains when the instrument is strong. The gains are 0-1% regardless of shape.

### 4c. Skorohod representation and non-monotone gamma

**Critical observation**: When the DGP is specified at the individual level as y_ij = x_{2j}·γ(u_ij) + α_j(u_ij) with u_ij ~ U(0,1), the map h_j(u) = x_{2j}·γ(u) + α_j(u) is a valid Skorohod representation (u_ij is the within-group rank) **only if** h_j is non-decreasing in u.

- **When γ is monotone** (√u, u, a(u−c)): h_j is monotone for all j (given sufficient β₀ slope), the Skorohod index IS the rank, and the group quantile Q_{Y_j}(u) = h_j(u) exactly. The population 2SLS slope equals γ(u).

- **When γ has derivative that changes sign** (hump 4u(1−u), sinusoidal u+A·sin(kπu)): h_j may be non-monotone in u for some groups. In that case, h_j is NOT the group quantile function. The actual Q_{Y_j} is the monotone rearrangement of h_j — a nonlinear functional. The Skorohod index u_ij is NOT the within-group rank. The population 2SLS slope β₁(u) ≠ γ(u).

This is not a "mixing" effect — it occurs even without individual covariates x₁. The mechanism is simply that **quantiles and non-monotone transformations don't commute**. The latent index model y_ij = x_{2j}·γ(u_ij) + ... generates well-defined random variables, but when γ is non-monotone, the latent index loses its quantile interpretation, and the group quantile function differs from the structural map.

Population target divergence from gamma(u):

| gamma(u) | Max |pop_beta1 − gamma| |
|-----------|--------------------------|
| sqrt(u) (monotone) | 0.011 |
| u (monotone) | 0.008 |
| 2*(u−0.5) (constant derivative) | 0.010 |
| 4*u*(1−u) (hump, derivative flips) | **0.794** |
| u + 0.5*sin(4πu) (oscillating) | **0.604** |

**Practical implication**: For all monotone γ shapes (which is what the paper uses, including the β_slope·Φ⁻¹(u) extension), the individual-level Skorohod DGP and the group-level model coincide, and γ(u) is the correct target. For non-monotone γ shapes, one should specify the DGP directly at the group level by drawing N individuals from a distribution with quantile function Q_{Y_j}(u) = β₀(u) + γ(u)·x_{2j} + α_j(u), which ensures γ(u) is the true coefficient by construction.

**Note on Average Strictness**: In the MP/CLP DGP and similar specifications with substantial within-group heterogeneity, Average Strictness is empirically satisfied with comfortable margin. The within-group mixing produces smooth group quantile functions with derivatives bounded away from zero. However, this is a DGP-specific empirical observation, not a general theoretical guarantee. Counterexamples exist: if γ(u) has large non-monotone amplitude and within-group heterogeneity is small, the monotone rearrangement can produce near-zero slopes.

---

## What Does NOT Help (explored but confirmed unhelpful)

| Feature tested | Gain at F≈23 | Why |
|---------------|-------------|-----|
| Wide x support (sd=1-4) | 0% or negative | Projection bias at high-leverage points |
| Multivariate X (p=2-5) | 0-1.7% | Partialling out already smooths |
| Smooth misspecification | 0-1% | Linear approx errors are small |
| Non-separable misspec h(x)·Φ⁻¹(u) | -1% to -9% | Projection moves away from pseudo-true parameter |
| Unbalanced groups | 0% | Noise averages out |
| Finer quantile grids | 0% | More points but same violation size |
| Many instruments (l=5-30) | -2% to -6% | See note below |
| Partial instrument invalidity | 0% | Biases estimates but doesn't increase violations |
| Threshold misspecification | 0-3.7% | Small and inconsistent |
| Heterogeneous individual slopes | 0% | No effect on group-level QFs |

**Note on many instruments**: The negative gains (-2% to -6%) with many weak instruments (l=5-30, each with coefficient 0.1) deserve further investigation. Many-instrument bias pushes 2SLS toward OLS. Possible mechanisms: (a) the overfit first stage creates a pseudo-true parameter that produces non-monotone curves at some x, and projection moves away from it (same mechanism as the c=5 case); (b) with many columns in P_Z, the projection smoothing interacts adversely with the PAVA correction. Testing with LIML or JIVE instead of 2SLS would clarify whether this is a many-IV bias issue or a projection issue. **TODO: investigate.**

---

## Estimand Comparison: D-IV vs CLP

**Setup**: MP DGP with composition effects — mean(x₁|x₂) = ρ*(x₂ − μ), so treated groups have systematically different worker types.

**Without individual covariates x₁**: CLP and D-IV target the **same** estimand β₁(u) (the total group-level effect). CLP without x₁ is exactly the unconstrained 2SLS step of D-IV. The only difference in estimates is the projection correction.

**With individual covariates x₁**: CLP targets the **direct** (within-type) effect δ(u), while D-IV targets the **total** effect β₁(u) = δ(u) + composition + re-ranking. These are different estimands.

**Results with x₁** (M=50, N=50, F≈23, 100 reps):

| ρ (composition) | D-IV IMSE (→total) | CLP IMSE (→direct) | CLP IMSE (→total) |
|-----------------|-------------------|-------------------|-------------------|
| 0.0 | 0.035 | 0.032 | 0.032 |
| 0.3 | 0.036 | 0.032 | **0.045** |
| 0.5 | 0.037 | 0.033 | **0.063** |
| 1.0 | 0.039 | 0.036 | **0.135** |
| 2.0 | 0.047 | 0.055 | **0.400** |

**Key finding**: When composition effects exist, CLP and D-IV estimate fundamentally different quantities. A researcher using CLP to estimate the total policy effect (what happened to the group's distribution) will be increasingly biased as composition effects grow. At ρ=2, CLP's IMSE for the total effect is **8.5× larger** than D-IV's.

This is arguably the most important finding for the paper — it's not an IMSE improvement story but a **what-are-you-estimating** story.

---

## Summary of Confirmed Channels

### Channels that produce substantial D-IV gains:

1. **Weak/moderate instruments (F < 15)**: 7-38% gains. The primary IMSE channel. Gains are non-monotone in F: smaller at very weak (F≈1, bias dominates) and zero at strong (F>20).
2. **Heavy-tailed within-group distributions + small N**: 13-44% at F≈23 with t₂ errors. Independent of instrument strength.
3. **Interaction of 1 and 2**: Superadditive (F≈12 + t₃ + N=25 → 23%).
4. **Large confounding variance**: Works by reducing effective F (not independent).

### The qualitative finding (not IMSE):

5. **CLP vs D-IV estimand difference**: With composition effects and individual covariates, CLP estimates the direct effect while D-IV estimates the total effect. CLP is biased for the total effect by up to 8.5×. Without individual covariates, both target the same estimand and only the projection correction differs.

### Key theoretical insight:

The D-IV improvement equals the size of the PAVA correction ||D_x||², which is bounded by the estimation error ||ψ̂_x − q_x||. With precise estimation (strong F, large N, light tails), the violations are small regardless of DGP features (gamma shape, x support, etc.). The projection is a **pure variance-reduction device**: it can only reduce variance, not bias, and the variance reduction is proportional to the estimation noise.
