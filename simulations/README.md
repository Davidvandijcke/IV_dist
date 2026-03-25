# FIVR Simulations

Comprehensive Monte Carlo study for the Functional Instrumental Variable Regression paper.

## Design Principle

> **All simulations use heterogeneous β(u).** There is no point estimating quantile-specific treatment effects if β is constant across quantiles.

Default DGP: `β(u) = β_base + β_slope × Φ⁻¹(u)` with `β_slope = 0.3`

## Main Script

**`fivr_simulations.R`** — consolidated simulation framework

### Structure

| Section | Description |
|---------|-------------|
| PART 0 | Setup and Configuration |
| PART 1 | Estimator Functions (2SLS, CoefProj, PW-FIVR) |
| PART 2 | Data Generating Processes |
| PART 3 | Simulation Framework (parallel execution) |
| PART 4 | Main Simulations |
| PART 5 | Results Analysis and Plots |

### Simulations

All use **heterogeneous β(u)** by default:

1. **Sample Size Study** — Varying M (groups) and N (individuals per group)
2. **IV Strength Study** — Varying π_Z (instrument relevance)
3. **Heterogeneity Study** — Varying β_slope (degree of effect heterogeneity)
4. **Multivariate X Study** — Testing with p > 1 covariates
5. **Extreme Cases** — Small samples, weak IV, high heterogeneity

### Usage

```r
# Full run (takes ~30 mins with 500 reps)
source("simulations/fivr_simulations.R")

# Or load functions and run custom simulations
source("simulations/fivr_simulations.R")
results <- run_simulation(
  config = list(name = "my_test"),
  dgp_fn = dgp_heterogeneous,
  dgp_args = list(M = 50, N = 50, pi_Z = 0.3, beta_base = 1.0, beta_slope = 0.4),
  n_reps = 500
)
```

### Outputs

- `results/simulation_results.rds` — Full results object
- `results/plot_*.png` — Publication-ready figures

## Key Findings (2026-01-25)

### CoefProj/PW vs Unconstrained 2SLS

| IV Strength (π_Z) | MSE Improvement |
|-------------------|-----------------|
| 0.10 (weak) | **22%** |
| 0.20 (moderate) | **5-8%** |
| 0.25 | **4-7%** |
| 0.50 (strong) | **<1%** |

### CoefProj vs PW-FIVR

| Condition | CoefProj Advantage |
|-----------|-------------------|
| Scalar X | ~0% (equivalent) |
| Multivariate X (p=2) | 0.5-0.7% |
| Multivariate X (p=3) | **1.5-2.2%** |

**CoefProj advantage grows with:**
- Higher p (more vertices for joint projection)
- Weaker IV (more monotonicity violations)
- Higher β_slope (steeper β(u) → larger violations)

**Bottom line:** 
- Use projection methods when IV is weak-moderate (π_Z < 0.3)
- CoefProj > PW-FIVR when p > 1
- Strong IV (π_Z ≥ 0.5): projection barely helps (estimates already precise)
