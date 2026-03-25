# CLAUDE.md - FIVR (Functional IV Regression) Project

## Overview

Developing econometric methods for IV estimation when outcomes are distributional/functional. This is David's job market paper.

## User Context

- **Name:** David Van Dijcke
- **Role:** PhD Economics, University of Michigan (6th year)
- **Status:** On the job market 2025-26
- **Timezone:** America/Detroit (EST/EDT)
- **Research Focus:** Econometrics, statistics, causal inference
- **Languages:** R, Python, Stata, LaTeX

## Working Mode

Research mode on FIVR project:
- Technical, efficient, no fluff
- Deep on econometrics/statistics
- Proactive about simulation design and code quality
- Remember the math — this is rigorous research

## Session Startup

Before doing anything else each session:
1. Read `PROJECT.md` — current project status & goals
2. Read `memory/YYYY-MM-DD.md` (today + yesterday) for recent context
3. Read `MEMORY.md` for project-specific long-term memory

## Memory Architecture

**Dual memory system:**
- `MEMORY.md` — FIVR-specific learnings, decisions, technical insights
- `memory/*.md` — daily logs of FIVR work
- `/Users/davidvandijcke/clawd/MEMORY.md` — cross-project learnings (global)

**Rule:** Always update local MEMORY.md. Update global MEMORY.md only for insights that apply beyond this project.

## Project Structure

```
IV_dist/
├── simulations/
│   ├── fivr_simulations.R    # Main consolidated script
│   ├── results/              # Output directory
│   └── README.md             # Documentation
├── code/
│   └── archive/              # Old simulation versions (v1-v6)
├── paper/                    # LaTeX/manuscript files
├── data/                     # Datasets
└── memory/                   # Daily session logs
```

## Core Principle

### ⚠️ ALWAYS Use Heterogeneous β(u) in Simulations
There is no point estimating quantile-specific treatment effects if β is constant across quantiles. All simulation DGPs must use heterogeneous β(u):
```r
# CORRECT: β(u) = β_base + β_slope * Φ⁻¹(u)
dgp_heterogeneous(M, N, q_grid, beta_base = 1.0, beta_slope = 0.3)
```

## Key Files

- **Paper:** `current_version_v4.tex`
- **New coefficient projection doc:** `coefficient_projection_v1.tex`
- **Main simulation script:** `simulations/fivr_simulations.R`
- **Code:** `code/` directory

## Estimators

1. **2SLS (unconstrained)** — baseline IV estimator
2. **CoefProj** — joint monotonicity projection via QP at vertices
3. **PW-FIVR** — pointwise projection then OLS

## Key Technical Insights

### Coefficient Projection Approach
For scalar X ∈ [x̲, x̄], only need to check monotonicity at endpoints:
- γ₋(u) = β₀(u) + (x̲ - μ)β₁(u) must be monotone
- γ₊(u) = β₀(u) + (x̄ - μ)β₁(u) must be monotone

### When Projection Helps (with heterogeneous β)
- **Weak IV** (π_Z=0.1): **22% MSE improvement**
- **Moderate IV** (π_Z=0.2-0.25): **4-8% improvement**
- **Strong IV** (π_Z≥0.5): <1% improvement
- **Small M** (15-25 groups): Larger gains than large M

### Benchmark Results (2026-01-25, all heterogeneous β)

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

### CoefProj vs PW-OLS

**Scalar X**: Essentially equivalent (<0.5% difference)

**Multivariate X**: CoefProj wins by 0.5-3%
- p=3, moderate hetero (β_slope=0.3): **+2.2%** over PW
- p=3, high hetero (β_slope=0.5): **+1.7%** over PW
- p=2, very weak IV (π_Z=0.1): **+0.7%** over PW

**Equivalence conditions:**
- **Equivalent when**: X ≥ μ (or ≤ μ) w.p.1
- **CoefProj wins when**: X spans both sides of μ
- Reason: PAVA adjustments are opposite on different sides of μ; they cancel in PW-OLS regression

**Divergence grows with:**
1. Higher p (more vertices for joint projection)
2. Weaker IV (more violations)
3. Stronger β heterogeneity (steeper β(u))

## Working Style

**Be direct.** Skip filler phrases.

**Be resourceful.** Check existing code and memory before asking questions.

**For simulations:** Always verify the DGP uses heterogeneous β(u). Check `MEMORY.md` for past simulation bugs.

**Write things down.** Update this file with important findings and decisions.

## Development Notes

### Implementation
- Joint QP projection preferred over ADMM (faster for p ≤ 8)
- Parallel processing via mclapply for simulation speed (N-1 cores)
- Results in `simulations/results/`

### Performance (2026-01-25)
**Bottleneck**: PW-FIVR calls PAVA M times per estimation (91% of runtime)

**Fixes applied:**
1. Vectorized PW-FIVR: matrix outer product for weights (1.4x speedup)
2. Parallel reps via `mclapply` across N-1 cores (~7x on 8-core machine)

**Result**: ~10x total speedup (3.1s for 7 DGPs × 300 reps)

### DGP Evolution
- **v3 bug (fixed in v4)**: Generated non-monotone quantile functions
- **v4 fix**: Bounded heterogeneity to ensure valid quantile functions
- Default heterogeneity: `beta_slope = 0.3` (moderate, realistic)

## Paper Framing

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

**Corollaries:**
- **Correct specification:** When β^unc = β ∈ C, we have β* = β and the improvement is relative to truth
- **Misspecification:** Even when linear model is wrong, we improve in estimating the best constrained approximation

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

## Current Status (2026-01-25)

### Next Steps
1. Run full simulation study (500+ reps)
2. Generate publication-ready tables
3. Theoretical analysis of CoefProj vs PW-FIVR equivalence conditions

## Open Questions

- Formal equivalence proof for PW-OLS ≈ CoefProj (scalar X case)
- Why doesn't CoefProj advantage monotonically increase with p?

## Safety

- Don't exfiltrate private data
- `trash` > `rm`
- Ask before external actions (emails, API calls, etc.)

---

*Last updated: 2026-02-10*
