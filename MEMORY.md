# MEMORY.md - FIVR Project Memory

## Core Principle

### ⚠️ ALWAYS Use Heterogeneous β(u) in Simulations
**There is no point estimating quantile-specific treatment effects if β is constant across quantiles.** If someone wants a constant treatment effect, they'd just run standard 2SLS.

All simulation DGPs must use heterogeneous β(u):
```r
# CORRECT: β(u) = β_base + β_slope * Φ⁻¹(u)
dgp_heterogeneous(M, N, q_grid, beta_base = 1.0, beta_slope = 0.3)

# WRONG: Constant β — pointless for quantile IV
dgp_homogeneous(...)  # Never use this as main simulation
```

Default heterogeneity: `beta_slope = 0.3` (moderate, realistic)

## Technical Insights

### DGP Bug (v3 → v4, 2026-01-25)
v3 DGP was broken — generated non-monotone QFs. Fixed in v4 by bounding heterogeneity.

### When Projection Helps (with heterogeneous β)
- **Weak IV** (π_Z=0.1): **22% MSE improvement**
- **Moderate IV** (π_Z=0.2-0.25): **4-8% improvement**
- **Strong IV** (π_Z≥0.5): <1% improvement (estimates already precise)
- **Small M** (15-25 groups): Larger gains than large M

### CoefProj vs PW-OLS Equivalence
**Equivalent when**: X ≥ μ (or ≤ μ) w.p.1
**CoefProj wins when**: X spans both sides of μ

Reason: PAVA adjustments are opposite on different sides of μ. When X is centered, they cancel in PW-OLS regression.

### Multivariate Implementation
- QP efficient up to p=8
- CoefProj beats PW-OLS by 0.5-3% consistently
- Improvement doesn't clearly increase with p (open question)

### Simulation Performance (2026-01-25)
**Bottleneck**: PW-FIVR calls PAVA M times per estimation (91% of runtime)
**Fixes applied**:
1. Vectorized PW-FIVR: matrix outer product for weights (1.4x)
2. Parallel reps via `mclapply` across N-1 cores (~7x on 8-core machine)
**Result**: ~10x total speedup (3.1s for 7 DGPs × 300 reps)

## Key Files
- **`simulations/fivr_simulations.R`** — consolidated simulation framework
- `simulations/README.md` — documentation
- `code/archive/` — old simulation versions (v1-v6)
- `memory/2026-01-25.md` — detailed session notes

### CoefProj vs PW-OLS Divergence (2026-01-25)
**With heterogeneous β(u)**, CoefProj beats PW-OLS in multivariate case:
- p=3, moderate hetero (β_slope=0.3): **+2.2%** over PW
- p=3, high hetero (β_slope=0.5): **+1.7%** over PW
- p=2, very weak IV (π_Z=0.1): **+0.7%** over PW

**Scalar X**: CoefProj ≈ PW (essentially equivalent, <0.5% difference)

**Divergence grows with:**
1. Higher p (more vertices for joint projection)
2. Weaker IV (more violations)
3. Stronger β heterogeneity (steeper β(u))

### Benchmark Simulation Results (2026-01-25, all heterogeneous β)

**IV Strength Study** (M=30, N=25, β_slope=0.3):
| π_Z | 2SLS MSE | Coef Gain |
|-----|----------|-----------|
| 0.10 | 162.2 | **22.5%** |
| 0.20 | 6.0 | 5.4% |
| 0.30 | 1.5 | 5.1% |
| 0.50 | 0.04 | 0.8% |

**Sample Size Study** (π_Z=0.25, β_slope=0.3):
- M=25, N=10: **10.2%** improvement
- M=25, N=25: **6.7%** improvement
- M=100, N=25: **1.1%** improvement

## Open Questions
- Formal equivalence proof for PW-OLS ≈ CoefProj (scalar X case)
- Why does CoefProj advantage not monotonically increase with p?

## Paper Framing (2026-01-XX)

### Unified Target Framework
The key insight: define a single target that works regardless of specification:
```
β* = Π_C(β^unc) = best constrained linear IV approximation
```

- **Correct specification:** β^unc ∈ C ⟹ β* = β^unc = β (true parameter)
- **Misspecification:** β^unc ∉ C ⟹ β* = projection of β^unc onto C (best we can do while maintaining valid distributions)

### L² Improvement Theorem (General)
**Theorem:** Under support containment (X ⊆ X_n):
```
‖β̂* - β*‖ ≤ ‖β̃ - β*‖  with prob 1
```

**Proof:** Uses contraction property of projections onto closed convex sets:
- β* ∈ C ⊆ C_n (by definition of β* as projection onto C)
- Therefore Π_{C_n}(β*) = β*
- By contraction: ‖Π_{C_n}(β̃) - Π_{C_n}(β*)‖ ≤ ‖β̃ - β*‖

**Corollary (Correct Specification):** When β^unc = β ∈ C, we have β* = β and the improvement is relative to truth.

**Corollary (Misspecification):** Even when linear model is wrong, we improve in estimating the best constrained approximation.

### Fréchet Regression Connection
- **Σ_XX-weighted L² norm = expected W₂² distance** (integrated over X distribution)
- **Coefficient projection = constrained Fréchet IV regression**: minimize expected W₂² to IV estimates, subject to output being valid probability measure
- The projection target β* is the Wasserstein barycenter path closest to unconstrained solution, constrained to valid distributions

### First Stage is IV, Not OLS
The quantile function is the *outcome*, not the target:
```
Q_{Y_j}(u) = β₀(u) + β₁(u)ᵀ(X_j - μ_X) + ε_j(u)
E[ε_j(u) | Z_j] = 0  (IV exclusion)
```
This is 2SLS applied to quantile functions as functional outcomes.
