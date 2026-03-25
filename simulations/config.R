# =============================================================================
# Simulation Configuration
# =============================================================================
#
# Four experiments, all built on the Melly-Pons DGP framework.
#
#   1. iv_strength    - Varying instrument strength (the main result)
#   2. sample_size    - Varying (M, N) at pi_Z=1 (MP-comparable)
#   3. heavy_tails    - Heavy-tailed base distribution + small N (independent channel)
#   4. heterogeneity  - Varying beta_slope at moderate IV (amplification channel)
#
# =============================================================================

Q_GRID <- seq(0.05, 0.95, by = 0.05)   # 19 points

ESTIMATOR_REGISTRY <- list(
  "2sls" = "estimate_2sls",
  "div"  = "estimate_div"
)


EXPERIMENTS <- list(

  # -----------------------------------------------------------------------
  # 1. IV strength: the key result.
  # -----------------------------------------------------------------------
  iv_strength = list(
    dgps = "dgp_mp",
    param_grid = data.frame(
      M             = 50,
      N             = 50,
      endogenous    = TRUE,
      heterogeneity = TRUE,
      pi_Z          = c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0),
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  ),

  # -----------------------------------------------------------------------
  # 2. Sample size at pi_Z=1 (direct comparison with MP/CLP).
  # -----------------------------------------------------------------------
  sample_size = list(
    dgps = "dgp_mp",
    param_grid = expand.grid(
      M             = c(25, 50, 100, 200),
      N             = c(25, 50, 100, 200),
      endogenous    = TRUE,
      heterogeneity = TRUE,
      pi_Z          = 1.0,
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  ),

  # -----------------------------------------------------------------------
  # 3. Heavy-tailed base distribution: independent channel via noisy group
  #    QF estimation. Strong instrument (pi_Z=1).
  #    The heavy tails are built into beta0(u) (the base quantile function),
  #    so the population target gamma(u) = sqrt(u) is unchanged.
  # -----------------------------------------------------------------------
  heavy_tails = list(
    dgps = "dgp_mp",
    param_grid = expand.grid(
      M             = 50,
      N             = c(10, 25, 50),
      endogenous    = TRUE,
      heterogeneity = TRUE,
      pi_Z          = 1.0,
      base_dist     = c("normal", "t3", "t2", "lognormal"),
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  ),

  # -----------------------------------------------------------------------
  # 4. Treatment effect heterogeneity at moderate IV (pi_Z=0.5).
  #    Shows that beta_slope amplifies weak-IV gains.
  #    Also run at pi_Z=1 to confirm no effect with strong instruments.
  # -----------------------------------------------------------------------
  heterogeneity = list(
    dgps = "dgp_mp",
    param_grid = expand.grid(
      M             = 50,
      N             = 50,
      endogenous    = TRUE,
      heterogeneity = TRUE,
      pi_Z          = c(0.5, 1.0),
      beta_slope    = c(0.0, 0.1, 0.2, 0.3),
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  )
)
