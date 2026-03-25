# MEMORY.md - D-IV Project Memory

## Core Principle

### ALWAYS Use Heterogeneous β(u) in Simulations
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

### Simulation Performance (2026-01-25)
**Bottleneck**: D-IV calls PAVA M times per estimation (91% of runtime)
**Fixes applied**:
1. Vectorized D-IV: matrix outer product for weights (1.4x)
2. Parallel reps via `mclapply` across N-1 cores (~7x on 8-core machine)
**Result**: ~10x total speedup (3.1s for 7 DGPs × 300 reps)

### Benchmark Simulation Results (2026-01-25, all heterogeneous β)

**IV Strength Study** (M=30, N=25, β_slope=0.3):
| π_Z | 2SLS MSE | D-IV Gain |
|-----|----------|-----------|
| 0.10 | 162.2 | **22.5%** |
| 0.20 | 6.0 | 5.4% |
| 0.30 | 1.5 | 5.1% |
| 0.50 | 0.04 | 0.8% |

**Sample Size Study** (π_Z=0.25, β_slope=0.3):
- M=25, N=10: **10.2%** improvement
- M=25, N=25: **6.7%** improvement
- M=100, N=25: **1.1%** improvement

## CLP Empirical Application (2026-03-25)

### Data
- `CLP_regression_data.dta`: 1444 obs per class (CZ × decade), pre-computed quantile changes
- Classes: all, m, f, c, nc, mfg, nmfg, etc.
- Key vars: `d_tradeusch_pw` (endogenous), `d_tradeotch_pw_lag` (instrument), `timepwt48` (weights), `statefip` (cluster)
- Quantile outcomes: `d_p5_lnwkwage` through `d_p95_lnwkwage` (19 quantiles)
- ADH mean: `adh_d_avg_lnwkwage_` (note trailing underscore in raw data)

### Results
- CLP replication exact: our `estimate_2sls` with controls + weights matches `ivreg` to machine precision
- ADH mean effects: -0.759 (all), -0.614 (females), -0.892 (males)
- D-IV projection active at 100% of observations; max coefficient adjustment ~0.15-0.19 log points
- D-IV smooths wiggles in 2SLS estimates, especially at upper quantiles

## Paper Framing

### Unified Target Framework
The key insight: define a single target that works regardless of specification:
```
β* = Π_C(β^unc) = best constrained linear IV approximation
```

- **Correct specification:** β^unc ∈ C ⟹ β* = β^unc = β (true parameter)
- **Misspecification:** β^unc ∉ C ⟹ β* = projection of β^unc onto C (best we can do while maintaining valid distributions)

### Fréchet Regression Connection
- **Σ_XX-weighted L² norm = expected W₂² distance** (integrated over X distribution)
- **Coefficient projection = constrained Fréchet IV regression**: minimize expected W₂² to IV estimates, subject to output being valid probability measure

### First Stage is IV, Not OLS
The quantile function is the *outcome*, not the target:
```
Q_{Y_j}(u) = β₀(u) + β₁(u)ᵀ(X_j - μ_X) + ε_j(u)
E[ε_j(u) | Z_j] = 0  (IV exclusion)
```
This is 2SLS applied to quantile functions as functional outcomes.

## Open Questions
- Why does D-IV adjustment magnitude not monotonically increase with sample size?
