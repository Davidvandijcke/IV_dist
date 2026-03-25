# D-IV Simulation Findings

Systematic exploration of when the D-IV projection improves over unconstrained 2SLS. Over 100 DGP configurations tested. All results use the Melly-Pons/CLP DGP framework.

## The D-IV Estimator

The D-IV estimator (Section 3.3 of the paper):
1. Compute unconstrained 2SLS coefficients β̃₀(u), β̃₁(u)
2. Evaluate fitted curves ψ̂(Xⱼ, u) = β̃₀(u) + β̃₁(u)(Xⱼ − μ̂) at each observed Xⱼ
3. Project each curve to monotonicity via PAVA: Q̂(Xⱼ, ·) = Π_Q(ψ̂(Xⱼ, ·))
4. Recover coefficients by OLS of Q̂(Xⱼ, u) on (1, Xⱼ − μ̂)

The improvement theorem (Theorem 1) guarantees that the joint weighted MSE of (β̂₀, β̂₁) under D-IV is ≤ that of unconstrained 2SLS, for any reference coefficients b that produce valid (monotone) quantile functions at the observed X values.

## The DGP

All experiments use the MP (2025) / CLP (2016) DGP (Section 4.2), implemented as `dgp_mp()` with parameters controlling the design:

```
y_ij = beta0(u_ij) + x_{1ij} * (u_ij/2) + x_{2j} * gamma(u_ij) + alpha_j(u_ij)
```

- u_ij ~ U(0,1) is the within-group rank (Skorohod representation)
- x_{1ij} ~ exp(0.25 * N(0,1)), x_{2j} group-level treatment
- gamma(u) = sqrt(u) + beta_slope * Phi^{-1}(u) (treatment effect)
- alpha_j(u) = u * eta_j - u/2, eta_j ~ U(0,1) (group heterogeneity)
- Endogenous: x_{2j} = pi_Z * z_j + eta_j + sqrt(1-pi_Z²) * nu_j
- z_j, nu_j ~ exp(0.25 * N(0,1))

The `base_dist` parameter controls beta0(u):
- `"normal"`: beta0(u) = 0 (original MP specification)
- `"t2"`, `"t3"`: beta0(u) = qt(u, df) (heavy-tailed base distribution)
- `"lognormal"`: beta0(u) = qlnorm(u, 0, 1.5) - mean (skewed, heavy right tail)

Since gamma(u) = sqrt(u) is monotone, the Skorohod map is always valid: u_ij IS the within-group rank, and gamma(u) IS the true population 2SLS coefficient. This was verified numerically: max|pop_β₁ − sqrt(u)| < 0.03 for all base distributions (M=2000, N=500, 10 reps).

---

## Channel 1: Instrument Strength (the primary driver)

**Setup**: MP DGP with pi_Z scaling instrument relevance.

**Results** (M=50, N=50, 500 reps):

| pi_Z | Median F | IMSE(2SLS) | IMSE(D-IV) | D-IV gain |
|------|----------|-----------|-----------|-----------|
| 0.10 | ~1       | 1,073     | 909       | 15.3%     |
| 0.20 | ~5       | 59.5      | 48.7      | 18.3%     |
| 0.30 | ~8       | 102       | 63.2      | 37.9%     |
| 0.50 | ~15      | 0.79      | 0.74      | 6.5%      |
| 0.70 | ~25      | 0.09      | 0.08      | 11.9%     |
| 1.00 | ~40      | 0.02      | 0.02      | 0.0%      |

**Mechanism**: Weak instruments inflate the variance of β̃₁(u). The resulting large deviations from the true (monotone) coefficient functions produce frequent and large violations of monotonicity in the fitted ψ̂(x, u) curves. PAVA corrects these violations.

**F-stat threshold**: Gains appear below F ≈ 15 and are substantial below F ≈ 10. Above F ≈ 20, gains vanish.

**Non-monotonicity in gains**: The gain at pi_Z=0.1 (F≈1) is smaller than at pi_Z=0.3 (F≈8). At very weak instruments (F ≈ 1), 2SLS bias dominates MSE for both estimators, limiting the relative gain from variance reduction via projection.

**Relevance**: Many applied IV settings have F-stats in the 5-15 range (Bartik/shift-share instruments, judge designs, etc.). The original MP/CLP DGP has F ≈ 23-40 depending on M, which is stronger than many applied settings.

---

## Channel 2: Heavy-Tailed Base Distribution + Small N

**Setup**: MP DGP with pi_Z=1 (strong instrument, F≈23) but heavy-tailed base distribution beta0(u). The heavy tails are built into the structural quantile function, so the population target gamma(u) = sqrt(u) is unchanged. The heavier tails make within-group distributions harder to estimate from N individual observations.

**Results** (M=50, pi_Z=1, 500 reps):

| Base distribution | N=10 | N=25 | N=50 |
|-------------------|------|------|------|
| Normal (original MP) | 0.1% | 0.1% | 0.0% |
| t₃               | 1.1% | 0.1% | 0.0% |
| t₂               | **20.1%** | 0.4% | 0.0% |
| Lognormal(0, 1.5) | **7.1%** | 0.9% | 0.1% |

**Mechanism**: Heavy tails in the base distribution cause outliers in the individual outcome data. With small N, these outliers distort the sample quantile functions Q̂_{Y_j}(u), which propagates through 2SLS into non-monotone coefficient estimates. This is a **first-stage noise channel** (noisy Q̂_{Y_j}), distinct from the second-stage noise channel (weak instruments affecting β̃₁).

**Why MP misses this**: The original MP DGP has beta0(u) = 0, so the within-group distribution comes entirely from x₁·(u/2) + x₂·gamma(u) + alpha_j(u), which has light tails (lognormal × uniform, bounded group effects). The sample quantile function converges quickly. Real economic outcomes (wages, health costs, test scores) often have much heavier tails.

**Key results**: t₂ base + N=10 gives **20.1% gain even with a strong instrument (F≈23)**. This is a genuinely independent channel from instrument strength.

---

## Channel 3: Treatment Effect Heterogeneity (amplification)

**Setup**: MP DGP with gamma(u) = sqrt(u) + beta_slope · Phi^{-1}(u). Tested at both pi_Z=0.5 (moderate IV) and pi_Z=1.0 (strong IV).

**Results** (M=50, N=50, 500 reps):

| pi_Z | beta_slope=0 | beta_slope=0.1 | beta_slope=0.2 | beta_slope=0.3 |
|------|-------------|----------------|----------------|----------------|
| 0.5  | 6.5%        | 6.7%           | 7.5%           | **8.4%**       |
| 1.0  | 0.0%        | 0.0%           | 0.0%           | 0.0%           |

**Mechanism**: More variation in gamma(u) across quantiles means β̃₁(u) has more "room" for noise-induced violations. At moderate IV (pi_Z=0.5), increasing beta_slope from 0 to 0.3 raises the gain from 6.5% to 8.4% — a 29% amplification. At strong IV (pi_Z=1.0), the effect is zero regardless of beta_slope.

**Important**: Heterogeneity does NOT independently produce gains. It **amplifies** the weak-IV channel. The noise magnitude is fixed by F; heterogeneity makes the signal more fragile, so the same noise causes more violations. But when the noise is small (strong F), heterogeneity doesn't matter.

---

## Channel 4: Large Group Heterogeneity

**Setup**: MP DGP but with eta_j ~ N(0, sigma_eta) instead of U(0,1).

**Results** (M=50, N=50, 500 reps):

| sigma_eta | Median F | D-IV gain |
|-----------|----------|-----------|
| 0.29 (≈ MP) | ~23   | 0%        |
| 0.50      | ~11      | ~3%       |
| 1.00      | ~3       | ~23%      |
| 1.50      | ~2       | ~33%      |
| 2.00      | ~1       | ~9%       |

**Mechanism**: NOT an independent channel. Large sigma_eta dilutes the instrument's contribution to x₂ (eta variance swamps z variance), reducing the effective F-stat. The gains track F, not sigma_eta.

**Non-monotonicity**: Gains fall at sigma_eta=2.0 (F≈1) because 2SLS bias dominates MSE for both estimators, limiting relative gains from variance reduction.

---

## Interaction of Channels 1 and 2

The weak-IV and heavy-tail channels interact **superadditively** (tested earlier with inline code):

| Setting | Gain |
|---------|------|
| F≈12 alone (normal base, N=50) | 12% |
| Heavy tails alone (t₃ base, N=25, F≈23) | 1% |
| F≈12 + t₃ + N=25 | **23%** |
| F≈12 + t₂ + N=15 | **32%** |

---

## What Does NOT Help (explored but confirmed unhelpful at F≈23)

| Feature tested | Gain | Why |
|---------------|------|-----|
| γ(u) shape: sinusoidal, hump, sign-changing | 0-1% | Violations exist but tiny when β̃₁ is precise |
| Wide x support (sd=1-4) | 0% or negative | Projection bias at high-leverage points |
| Multivariate X (p=2-5) | 0-1.7% | Partialling out already smooths |
| Smooth misspecification | 0-1% | Linear approx errors are small |
| Non-separable misspec h(x)·Φ⁻¹(u) | -1% to -9% | Projection moves away from pseudo-true parameter |
| Unbalanced groups | 0% | Noise averages out |
| Finer quantile grids | 0% | More points but same violation size |
| Many instruments (l=5-30) | -2% to -6% | Needs investigation (possible many-IV bias interaction) |
| Partial instrument invalidity | 0% | Biases but doesn't create violations |

### Note on γ shape at high amplitude

When γ has very large amplitude (e.g., gamma = 5·sqrt(u)) relative to the beta0 slope, the fitted curve q(x, u) = beta0(u) + 5·sqrt(u)·(x−μ) can be non-monotone in u for large |x−μ|. This means no probability distribution has that quantile function — the linear model is **misspecified** at extreme covariate values. The unconstrained 2SLS converges to a pseudo-true parameter that itself produces non-monotone fitted curves at some observed X. The projection moves toward monotonicity but away from the pseudo-true target, producing **negative** IMSE gains (-13% at c=5). This is the regime where Theorem 1's conditions (reference coefficients producing valid QFs) are not met.

---

## Skorohod Representation and Non-Monotone γ

When the DGP is specified at the individual level as y_ij = x_{2j}·γ(u_ij) + α_j(u_ij), the map h_j(u) = x_{2j}·γ(u) + α_j(u) is a valid Skorohod representation (u_ij is the within-group rank) **only if** h_j is non-decreasing in u.

- **Monotone γ** (√u, u, a(u−c), √u + β_slope·Φ⁻¹(u) for moderate β_slope): h_j is monotone, the Skorohod index IS the rank, Q_{Y_j}(u) = h_j(u), and pop_β₁(u) = γ(u).

- **Non-monotone γ** (4u(1−u), u+0.5·sin(4πu)): h_j can be non-monotone. Then h_j is NOT the group QF. The actual Q_{Y_j} is the monotone rearrangement of h_j, and pop_β₁(u) ≠ γ(u). Measured divergence: max|pop_β₁ − γ| up to 0.8 for hump-shaped γ.

This is not a "mixing" effect — it occurs even without individual covariates x₁. The mechanism is that **quantiles and non-monotone transformations don't commute**. The latent index model generates well-defined random variables, but when γ is non-monotone, the latent index loses its quantile interpretation.

**For all monotone γ shapes used in the paper** (including the β_slope·Φ⁻¹(u) extension), the Skorohod representation is valid and γ(u) is the correct target.

### Average Strictness

In the MP/CLP DGP and similar specifications with substantial within-group heterogeneity, Average Strictness is empirically satisfied with comfortable margin. The within-group mixing produces smooth group quantile functions with derivatives bounded away from zero. However, this is a DGP-specific empirical observation, not a general theoretical guarantee.

---

## Estimand Comparison: D-IV vs CLP

**Without individual covariates x₁**: CLP and D-IV target the **same** estimand β₁(u). CLP without x₁ is exactly the unconstrained 2SLS step of D-IV. The only difference in estimates is the projection correction.

**With individual covariates x₁ and composition effects**: CLP targets the **direct** (within-type) effect δ(u), while D-IV targets the **total** effect β₁(u) = δ(u) + composition + re-ranking.

**Setup**: MP DGP with mean(x₁|x₂) = ρ·(x₂ − μ).

**Results** (M=50, N=50, F≈23, 100 reps):

| ρ | D-IV IMSE (→total) | CLP IMSE (→direct) | CLP IMSE (→total) |
|---|-------------------|-------------------|-------------------|
| 0.0 | 0.035 | 0.032 | 0.032 |
| 0.3 | 0.036 | 0.032 | **0.045** |
| 0.5 | 0.037 | 0.033 | **0.063** |
| 1.0 | 0.039 | 0.036 | **0.135** |
| 2.0 | 0.047 | 0.055 | **0.400** |

At ρ=2, CLP's IMSE for the total effect is **8.5× larger** than D-IV's. This is not an IMSE improvement story but a **what-are-you-estimating** story.

---

## Summary

### Channels that produce substantial D-IV gains:

1. **Weak/moderate instruments (F < 15)**: 6-38% gains. The primary IMSE channel.
2. **Heavy-tailed base distribution + small N**: 7-20% at F≈23 with t₂ or lognormal base. Independent of instrument strength.
3. **Interaction of 1 and 2**: Superadditive (F≈12 + t₃ + N=25 → 23%).
4. **Treatment effect heterogeneity at moderate IV**: Amplifies weak-IV gains by ~29% (6.5% → 8.4%). Zero effect at strong IV.

### The qualitative finding:

5. **CLP vs D-IV estimand difference**: With composition effects and individual covariates, CLP estimates the direct effect while D-IV estimates the total effect. CLP is biased for the total effect by up to 8.5×.

### Key theoretical insight:

The D-IV improvement equals the size of the PAVA correction ||D_x||², which is bounded by the estimation error ||ψ̂_x − q_x||. With precise estimation (strong F, large N, light tails), violations are small regardless of DGP features. The projection is a **variance-reduction device** proportional to estimation imprecision. Two sources of imprecision drive gains: (i) weak instruments (noisy β̃₁), and (ii) heavy-tailed within-group distributions with small N (noisy Q̂_{Y_j}).
