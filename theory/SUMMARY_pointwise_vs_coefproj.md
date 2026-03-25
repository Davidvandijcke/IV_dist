# Pointwise FIVR vs Coefficient Projection: Summary

**Date:** 2025-01-27

## The Question

In the linear quantile model Q(u|X) = α(u) + β(u)(X - μ_X):

- **Pointwise FIVR:** Projects Q̂(u|X_j) to be monotone at each observed X_j, then regresses
- **Coefficient Projection:** Projects endpoint QFs γ₋, γ₊ to be monotone, recovers β directly

Both enforce monotonicity. When are they equivalent?

---

## Main Finding: Practical Equivalence

**Despite theoretical differences, the methods produce essentially identical results.**

| DGP | Divergence (L2) | MSE Difference | Both Viols |
|-----|-----------------|----------------|------------|
| Oscillating β (freq=2) | 0.020 | +0.07% | 93% |
| Oscillating β (freq=4) | 0.014 | −0.03% | 74% |
| Positive X | 0.069 | −0.47% | 66% |
| Symmetric X | 0.012 | +0.05% | 45% |
| Hetero β, weak IV | 0.017 | +0.03% | 40% |
| Hetero β, strong IV | 0.002 | +0.49% | 1% |
| Extreme weak IV | 0.040 | +0.17% | 78% |

All MSE differences are **<0.5%** of 2SLS baseline. Methods are indistinguishable.

---

## Theoretical Conditions

### When Exactly Equivalent
1. **No violations:** Both γ₋ and γ₊ already monotone → both = unprojected 2SLS
2. **Same pooling structure:** PAVA pools at the same intervals for all X
3. **Only endpoint observations:** No intermediate X values to project

### When They (Theoretically) Diverge
1. **Different violation patterns:** γ₋ violates at [0.3,0.5], γ₊ at [0.6,0.8]
2. **Mixed-sign X:** Symmetric support around μ_X
3. **Oscillating β(u):** Creates asymmetric violations

The core issue: **PAVA is nonlinear.** π_Q(λf + (1-λ)g) ≠ λπ_Q(f) + (1-λ)π_Q(g)

---

## Why Divergence Doesn't Matter

Three mechanisms:

1. **Small violations:** Actual decrements are small → PAVA smooths similarly
2. **Approximate linearity:** For near-monotone functions, PAVA respects convex combinations
3. **Regression averaging:** OLS on projected values averages out small differences

---

## Recommendation

**Use Coefficient Projection as the default method.**

| Criterion | Coefficient Projection | Pointwise FIVR |
|-----------|------------------------|----------------|
| **Computation** | 2 PAVA operations | M PAVA operations |
| **Constraint** | Valid QF for ALL X | Valid QF for observed X only |
| **Multivariate** | 2^p vertices | M projections |
| **MSE** | Essentially equal | Essentially equal |

Coefficient Projection is simpler, cleaner, and just as good.

---

## Files Created

- `theory/pointwise_vs_coefproj_equivalence.tex` — Full LaTeX writeup with proofs
- `code/simulations_divergence_test.R` — Simulation code testing divergence
- `code/simulation_results_divergence.rds` — Simulation results
- `code/plot_divergence_test.png` — Visualization

---

## Next Steps

1. ✅ Formalized conditions for equivalence
2. ✅ Identified (theoretical) divergence conditions
3. ✅ Ran simulations confirming practical equivalence
4. ✅ Wrote up findings

**Open questions:**
- Inference/bootstrap for projection estimators?
- Extension to other functional forms beyond linear?
