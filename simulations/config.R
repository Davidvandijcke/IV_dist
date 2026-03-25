# =============================================================================
# Simulation Configuration
# =============================================================================
#
# Four experiments, all built on the Melly-Pons DGP framework.
#
#   1. iv_strength    - Varying instrument strength (the main result)
#   2. sample_size    - Varying (M, N) at pi_Z=1 (MP-comparable)
#   3. heavy_tails    - Heavy-tailed errors + small N (independent channel)
#   4. heterogeneity  - Varying treatment effect heterogeneity
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
  #    Endogenous MP DGP, varying pi_Z.
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
  # 2. Sample size: vary M and N at pi_Z=1 (direct comparison with MP/CLP).
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
  # 3. Heavy tails: independent channel via noisy group QF estimation.
  #    Strong instrument (pi_Z=1) to isolate the effect.
  # -----------------------------------------------------------------------
  heavy_tails = list(
    dgps = "dgp_mp",
    param_grid = expand.grid(
      M             = 50,
      N             = c(10, 25, 50),
      endogenous    = TRUE,
      heterogeneity = TRUE,
      pi_Z          = 1.0,
      error_df      = c(Inf, 3, 2),
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  ),

  # -----------------------------------------------------------------------
  # 4. Treatment effect heterogeneity: vary beta_slope at pi_Z=1.
  #    gamma(u) = sqrt(u) + beta_slope * Phi^{-1}(u)
  # -----------------------------------------------------------------------
  heterogeneity = list(
    dgps = "dgp_mp",
    param_grid = data.frame(
      M             = 50,
      N             = 50,
      endogenous    = TRUE,
      heterogeneity = TRUE,
      pi_Z          = 1.0,
      beta_slope    = c(0.0, 0.05, 0.1, 0.15, 0.2, 0.3),
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  )
)
