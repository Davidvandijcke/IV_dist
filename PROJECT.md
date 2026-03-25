# PROJECT.md - FIVR Status

## Overview

**FIVR** — Functional Instrumental Variable Regression

Developing econometric methods for IV estimation when outcomes are distributional/functional.

## Current Status (2026-01-25)

### Consolidated Simulation Framework
Replaced multiple scattered scripts with a single comprehensive simulation file:
- **`simulations/fivr_simulations.R`** — all estimators, DGPs, and simulation studies

### Estimators Implemented
1. **2SLS (unconstrained)** — baseline IV estimator
2. **CoefProj** — joint monotonicity projection via QP at vertices
3. **PW-FIVR** — pointwise projection then OLS

### Key Findings
- **All simulations use heterogeneous β(u)** — no point to quantile IV with constant effects
- **CoefProj/PW beats 2SLS** by:
  - 22% with weak IV (π_Z = 0.1)
  - 5-8% with moderate IV (π_Z = 0.2)
  - <1% with strong IV (π_Z ≥ 0.5)
- **CoefProj beats PW-FIVR** with multivariate X:
  - p=3: +2% advantage
  - Scalar X: essentially equivalent

### File Structure
```
simulations/
├── fivr_simulations.R    # Main consolidated script
├── results/              # Output directory
└── README.md             # Documentation

code/
├── archive/              # Old simulation versions (v1-v6)
└── [other project code]
```

## Next Steps

1. Run full simulation study (500+ reps)
2. Generate publication-ready tables
3. Theoretical analysis of CoefProj vs PW-FIVR equivalence conditions

## Key Decisions

- **Joint QP projection** preferred over ADMM (faster for p ≤ 8)
- **Bound β heterogeneity** in DGP to ensure valid quantile functions
- **Parallel processing** via mclapply for speed

## References

(Papers, resources, etc.)
