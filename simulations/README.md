# D-IV Simulations

Comprehensive Monte Carlo study for the D-IV paper.

## Design Principle

> **All simulations use heterogeneous β(u).** There is no point estimating quantile-specific treatment effects if β is constant across quantiles.

Default DGP: `β(u) = β_base + β_slope × Φ⁻¹(u)` with `β_slope = 0.3`

## Structure

| File | Description |
|------|-------------|
| `dgp.R` | Unified Melly-Pons DGP |
| `config.R` | Experiment configurations |
| `run_simulations.R` | MC runner (parallel via mclapply) |
| `analysis.R` | Paper output: tables + figures |
| `results/` | .rds output (gitignored) |

### Estimators

1. **2SLS** — unconstrained IV estimator (baseline)
2. **D-IV** — pointwise projection then OLS

### Simulations

All use **heterogeneous β(u)** by default:

1. **Baseline** — Validate implementation on MP DGP variants
2. **IV Strength** — Varying π_Z (instrument relevance)
3. **Sample Size** — Varying M (groups) and N (individuals per group)
4. **Heterogeneity** — Varying β_slope (degree of effect heterogeneity)

### Usage

```r
# Run all experiments
cd IV_dist && Rscript simulations/run_simulations.R

# Run one experiment
Rscript simulations/run_simulations.R iv_strength

# Interactive
source("simulations/run_simulations.R")
results <- run_experiment("iv_strength")
```

### Outputs

- `results/<experiment>.rds` — Per-experiment results
- `results/all_results.rds` — Combined results

## Key Findings

### D-IV vs Unconstrained 2SLS

| IV Strength (π_Z) | MSE Improvement |
|-------------------|-----------------|
| 0.10 (weak) | **22%** |
| 0.20 (moderate) | **5-8%** |
| 0.25 | **4-7%** |
| 0.50 (strong) | **<1%** |

**Bottom line:**
- Use D-IV when IV is weak-moderate (π_Z < 0.3)
- Strong IV (π_Z ≥ 0.5): projection barely helps (estimates already precise)
