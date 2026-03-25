# PROJECT.md - D-IV Status

## Overview

**D-IV** — Distribution-valued Instrumental Variable Regression

Developing econometric methods for IV estimation when outcomes are distributional/functional.

## Current Status (2026-03-25)

### Estimators Implemented
1. **2SLS (unconstrained)** — baseline IV estimator
2. **D-IV** — pointwise projection then OLS (the paper's primary estimator)

Both support multivariate X/Z (controls) and analytic weights.

### Inference Implemented (`R/inference.R`)
- **Pointwise CIs** — sandwich (HC0 + cluster-robust), verified against `ivreg` + `sandwich`
- **Uniform confidence bands** — multiplier bootstrap (Gaussian/Rademacher)
- **Projected bootstrap** — PAVA inside each bootstrap draw:
  - Uniform bands (for D-IV only)
  - Pointwise CIs (normal-based + percentile): ~6-9% tighter than sandwich
- **Cluster-robust** — clustered sandwich + cluster wild bootstrap

### Key Findings (Simulations)
- **All simulations use heterogeneous β(u)** — no point to quantile IV with constant effects
- **D-IV beats 2SLS** by:
  - 22% with weak IV (π_Z = 0.1)
  - 5-8% with moderate IV (π_Z = 0.2)
  - <1% with strong IV (π_Z ≥ 0.5)
- **Coverage** (run_coverage.R, 50 reps): pointwise ~95%, uniform ~90-94%
- **Projected bootstrap** pointwise CIs: 6% width reduction at π_Z=0.3, M=50

### Empirical Application
Replicated CLP (2016, Econometrica) — distributional effects of Chinese import competition on U.S. wages. Applied D-IV to same data.
- CLP replication is exact (max |2SLS - ivreg| = 0)
- D-IV smooths quantile-by-quantile estimates, max adjustment ~0.15-0.19 log points
- Figures include three-layer D-IV confidence ribbons (projected boot PW, sandwich PW, uniform)
- All inference clustered by state, weighted by CZ population
- Projected bootstrap SE ratio: 0.91 (full sample), 0.95 (F), 0.92 (M)
- See `empirical/clp_replication.R`

### File Structure
```
R/                        # Clean estimator code
├── estimators.R          # 2SLS, D-IV (weights, return_internals)
├── inference.R           # Pointwise CIs, uniform bands, bootstrap
├── utils.R               # PAVA, helpers
└── first_stage.R         # Group quantile computation
simulations/
├── dgp.R                 # Unified MP DGP (with base_dist)
├── config.R              # Experiment configurations
├── run_simulations.R     # MC runner
├── run_coverage.R        # Coverage simulation for CIs
├── analysis.R            # Paper output: tables + figures
└── results/              # .rds output
empirical/
├── clp_replication.R     # CLP (2016) replication + D-IV + confidence bands
└── figures/              # Output figures (3 PDFs with CI ribbons)
```

## Next Steps

1. Run full simulation study (500+ reps) with updated DGP configs
2. Generate publication-ready tables (including coverage)
3. Write up empirical application section for paper
4. Add inference section to paper (Section 4.3 + coverage table)

## Key Decisions

- **Bound β heterogeneity** in DGP to ensure valid quantile functions
- **Parallel processing** via mclapply for speed
- **Controls via multivariate X/Z**, weights via `weights` parameter
- **Weighted residuals** in sandwich: WLS influence function uses `sqrt(w) * xi`, not `xi`
- **Cluster bootstrap**: one omega per cluster (not per observation)

## References

(Papers, resources, etc.)
