# D-IV Simulation Findings

Systematic exploration of when the D-IV projection improves over unconstrained 2SLS. Over 150 DGP configurations tested across two DGP families: the original Melly-Pons/CLP DGP (lognormal, all-positive x₂) and a centered-x₂ DGP (normal instruments, x₂ spans negative values).

## The D-IV Estimator

The D-IV estimator (Section 3.3 of the paper):
1. Compute unconstrained 2SLS coefficients β̃₀(u), β̃₁(u)
2. Evaluate fitted curves ψ̂(Xⱼ, u) = β̃₀(u) + β̃₁(u)(Xⱼ − μ̂) at each observed Xⱼ
3. Project each curve to monotonicity via PAVA: Q̂(Xⱼ, ·) = Π_Q(ψ̂(Xⱼ, ·))
4. Recover coefficients by OLS of Q̂(Xⱼ, u) on (1, Xⱼ − μ̂)

Theorem 1 guarantees: joint Σ̂_XX-weighted MSE of D-IV ≤ that of 2SLS, for any reference coefficients producing valid QFs. Lemma 2 guarantees: W₂²(Q̂(Xⱼ,·), Q_true) ≤ W₂²(ψ̂_Xⱼ, Q_true) in every replication.

## Two DGP Families

### dgp_mp(): Original Melly-Pons/CLP DGP
- x₂ always positive (lognormal components), support ≈ [0.5, 4]
- gamma(u) = √u (unbounded derivative near 0)
- Strong monotonicity margin: steep intercept from lognormal mixture
- At M=50, pi_Z=1: F≈23, D-IV gain ≈ 0%

### dgp_centered(): Centered-x₂ DGP (more realistic for applied IV)
- x₂ centered (normal instruments), spans negative values
- gamma(u) = u + beta_slope·Φ⁻¹(u) (bounded derivative)
- β₀(u) = σ·Φ⁻¹(u) (or t-distributed for heavy tails)
- Heterogeneous first stage: x₂ = (π_Z + δ(η−0.5))·z + η + noise
- Optional controls: W_k ~ N(0,1) with effects u/(k+2)

The centered DGP is more representative of applied settings (import shocks, policy doses, relative prices) where the treatment variable can be positive or negative.

---

## Why the MP DGP Suppresses D-IV Gains

Three features of the MP DGP create an unusually large monotonicity margin:

1. **All-positive x₂**: With x₂ > 0 always, (Xⱼ − μ) is bounded below by ≈ −2.5. The violation condition β₀'(u) + β̃₁'(u)·(Xⱼ−μ) < 0 can only trigger at groups far below the mean, where the margin is thin.

2. **Steep intercept**: β₀(u) ≈ μ_X·√u has derivative ≈ 1.25/√u, providing a large buffer that absorbs estimation noise.

3. **Strong instrument**: At M=50, pi_Z=1, the median F≈23-40. The 2SLS is very precise.

With centered x₂, groups with negative Xⱼ have the second term *subtract* from the margin, making violations much more likely even at moderate F.

---

## Results: Centered-x₂ DGP

All results below use dgp_centered() with 500 replications unless noted.

### Channel 1: Instrument Strength (M=50, N=50)

| π_Z | F≈ | β₁ gain | W₂² gain | Invalid groups |
|-----|-----|---------|----------|----------------|
| 0.1 | ~1 | 14.2% | 11.1% | 42.9% |
| 0.2 | ~4 | **39.1%** | **31.4%** | 26.4% |
| 0.3 | ~8 | **27.7%** | **20.6%** | 14.1% |
| 0.5 | ~15 | 3.7% | 1.7% | 2.9% |
| 0.7 | ~25 | 0.1% | 0.0% | 0.7% |
| 1.0 | ~40 | 0.0% | 0.0% | 0.1% |

Non-monotonicity in gains at very weak IV (π_Z=0.1): 2SLS bias dominates MSE, limiting relative gains from variance reduction.

### Channel 2: Sample Size (π_Z=0.5, F≈12)

| (M, N) | β₁ gain | W₂² gain | Invalid groups |
|---------|---------|----------|----------------|
| (25, 25) | **38.4%** | **25.7%** | 14.7% |
| (25, 50) | **39.5%** | **22.6%** | 9.4% |
| (50, 50) | 3.7% | 1.7% | 2.9% |
| (100, 50) | 0.1% | 0.0% | 0.7% |

Gains concentrated at small M. At M=25, ~10-15% of groups have invalid fitted distributions.

### Channel 3: Controls + Heterogeneous First Stage (M=50, N=50, π_Z=0.5)

| p | δ (hetero FS) | β₁ gain | W₂² gain |
|---|---------------|---------|----------|
| 1 | 0 | 3.7% | 1.7% |
| 2 | 0 | 8.7% | 3.2% |
| 5 | 0 | **35.9%** | **19.8%** |
| 1 | 0.5 | 6.5% | 3.7% |
| 2 | 0.5 | **16.1%** | **6.9%** |
| 5 | 0.5 | **56.8%** | **44.3%** |
| 1 | 1.0 | 12.6% | 8.9% |
| 2 | 1.0 | **32.3%** | **20.1%** |
| 5 | 1.0 | **42.5%** | **25.7%** |

**Controls amplify gains dramatically.** More regressors add dimensions to ψ̂(x,u), making the monotonicity constraint harder to satisfy. Combined with first-stage heterogeneity, gains reach 57%.

**Heterogeneous first stage** (δ > 0) is a powerful independent channel. With δ=1, the first-stage coefficient varies by group: x₂ = (π_Z + δ(η−0.5))·z + η + noise. The linear 2SLS uses a single first-stage slope, creating group-specific misfit. This is realistic: Bartik/shift-share instruments have heterogeneous first stages across regions.

### Channel 4: Heavy Tails (M=50, π_Z=0.5)

| Base distribution | N=25 | N=50 |
|-------------------|------|------|
| Normal | 14.3% / 5.1% | 3.7% / 1.7% |
| t₅ | 15.5% / 6.1% | 4.5% / 2.3% |
| t₃ | **17.9%** / **7.9%** | 5.4% / 2.9% |

(Format: β₁ gain / W₂² gain)

Heavy tails in the base distribution make within-group quantile estimation harder, creating additional noise in Q̂_{Y_j}(u). Effect is modest (~3pp above normal) but compounds with other channels.

### Realistic Combination

p=3, δ=0.5, t₅ base, β_slope=0.2, F≈12, M=50, N=50:

| Metric | 2SLS | D-IV | Gain |
|--------|------|------|------|
| β₁ IMSE | 0.062 | 0.054 | **13.9%** |
| W₂² | 0.100 | 0.094 | **5.6%** |
| Invalid groups | 4.8% | 0% | — |

This mimics a labor/health economics IV study with a Bartik-style instrument, centered treatment, moderate IV strength, and a few controls.

---

## Results: MP Baseline (pi_Z=1, lognormal x₂)

| (M, N) | β₁ gain | Invalid groups |
|---------|---------|----------------|
| (25, 25) | 7.3% | 6.6% |
| (25, 50) | 5.0% | 3.3% |
| (50, 50) | 0.0% | 0.6% |
| (100, 50) | 0.0% | 0.1% |

Gains only at M=25, where MP note F<10 in 40% of draws.

---

## W₂² Improvement: The Headline Metric

Lemma 2 guarantees W₂²(Q̂, Q_true) ≤ W₂²(ψ̂, Q_true) in **every replication**. This was confirmed empirically: `always_better=YES` across all DGPs tested.

W₂² gains are always smaller than β₁ IMSE gains because W₂² includes both the improvement from projection AND the OLS residuals (which don't change). But W₂² is never negative, which makes it the cleanest result for the paper.

---

## Fraction of Groups with Invalid Fitted Distributions

This is the most interpretable diagnostic for applied researchers. At moderate IV (F≈12):
- p=1, no hetero FS: 2.9% invalid (1-2 groups out of 50)
- p=2, δ=0.5: 3.2% invalid
- p=5, δ=0.5: 3.7% invalid
- Weak IV (π_Z=0.2): 26.4% invalid (13 out of 50 groups)
- Small M=25: 10-15% invalid

The unconstrained 2SLS produces invalid probability distributions for these groups. D-IV fixes all of them.

---

## What Does NOT Help (confirmed at F≈23)

| Feature | Gain | Why |
|---------|------|-----|
| γ(u) shape (oscillating, hump, sign-changing) | 0-1% | Violations tiny when β̃₁ precise |
| Wide x support alone (all-positive) | 0% | High-leverage projection bias |
| Nonlinear treatment effect (x² interaction) | 0% | IV still valid, linear coef consistent |
| Nonlinear control effects (w² misspec) | 0% | Doesn't affect endogenous coef |
| Many instruments (l=5-30) | -2% to -6% | Many-IV bias interacts with projection |

---

## Skorohod Representation Notes

For monotone γ shapes (√u, u, a(u−c), u + β_slope·Φ⁻¹(u)):
- The Skorohod map h_j(u) is monotone → u_ij IS the within-group rank
- Q_{Y_j}(u) = h_j(u) exactly → γ(u) IS the population 2SLS coefficient

For non-monotone γ (e.g., 4u(1−u)):
- h_j can be non-monotone → u_ij is NOT the rank
- Q_{Y_j} ≠ h_j → pop_β₁ ≠ γ (divergence up to 0.8)
- Must specify DGP at the group level to use non-monotone γ

**Centered x₂ + √u**: The combination x₂<0 with γ=√u (unbounded derivative near 0) can cause h_j'(u) < 0 near u=0 for negative-x₂ groups. Use bounded γ' (like γ=u) with centered x₂.

---

## Estimand Comparison: D-IV vs CLP

With composition effects (mean(x₁|x₂) = ρ·(x₂−μ)):
- D-IV targets **total effect**: γ(u) + ρ·(u/2)
- CLP targets **direct effect**: γ(u)

At ρ=2: CLP's IMSE for the total effect is **8.5× larger** than D-IV's. Without x₁ controls, both target the same estimand.

---

## Summary of Channels

| Channel | β₁ gain range | W₂² gain | Independent? |
|---------|--------------|----------|-------------|
| Weak IV (F<15) | 4-39% | 2-31% | Primary driver |
| Small M (M=25) | 5-40% | 9-26% | Via effective F |
| Controls (p>1) | 6-36% | 2-20% | Yes (more ψ dimensions) |
| Hetero first stage (δ>0) | 7-13% | 4-9% | Yes (group-specific misfit) |
| Heavy tails (t₃) | +3pp | +2pp | Amplifies other channels |
| β_slope (heterogeneity) | +2pp at moderate F | — | Amplifies weak IV |
| Controls × Hetero FS | Up to 57% | Up to 44% | Superadditive |

### TODO
- Coverage experiments for confidence bands (code ready in R/inference.R and simulations/run_coverage.R)
- Discrete within-group distributions
- Over-identification (2-3 instruments for 1 endogenous)
