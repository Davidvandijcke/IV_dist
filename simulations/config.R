# =============================================================================
# Simulation Configuration
# =============================================================================
#
# Four experiments, all built on the Melly-Pons DGP framework.
#
#   1. baseline     - MP DGP 1-3 at reference sample sizes (validate implementation)
#   2. iv_strength  - Endogenous MP DGP, varying instrument strength
#   3. sample_size  - Endogenous MP DGP, varying (M, N)
#   4. heterogeneity - Endogenous MP DGP, varying treatment effect shape
#
# =============================================================================

Q_GRID <- seq(0.05, 0.95, by = 0.05)   # 19 points

ESTIMATOR_REGISTRY <- list(
  "2sls"         = "estimate_2sls",
  "div"          = "estimate_div",
  "div_endpoint" = "estimate_div_endpoint"
)


EXPERIMENTS <- list(

  # -----------------------------------------------------------------------
  # 1. Baseline: replicate the three MP DGP variants at reference (m,n)
  #    Shows D-IV works correctly on the established benchmark.
  # -----------------------------------------------------------------------
  baseline = list(
    dgps = "dgp_mp",
    param_grid = expand.grid(
      M             = c(25, 200),
      N             = c(25, 200),
      endogenous    = c(FALSE, FALSE, TRUE),
      heterogeneity = c(FALSE, TRUE,  TRUE),
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 1000
  ),

  # -----------------------------------------------------------------------
  # 2. IV strength: the key result.
  #    Endogenous MP DGP with group heterogeneity, varying pi_Z.
  #    Demonstrates D-IV gains grow as instrument weakens.
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
  # 3. Sample size: vary M (groups) and N (within-group).
  #    Endogenous MP DGP with moderate IV strength.
  #    Shows D-IV helps most when M is small.
  # -----------------------------------------------------------------------
  sample_size = list(
    dgps = "dgp_mp",
    param_grid = expand.grid(
      M             = c(25, 50, 100, 200),
      N             = c(25, 50, 100, 200),
      endogenous    = TRUE,
      heterogeneity = TRUE,
      pi_Z          = 0.5,
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  ),

  # -----------------------------------------------------------------------
  # 4. Treatment effect heterogeneity: vary beta_slope.
  #    gamma(u) = sqrt(u) + beta_slope * Phi^{-1}(u)
  #    beta_slope = 0 is homogeneous (original MP); higher = more spread.
  #    Moderate IV (pi_Z = 0.3), moderate samples.
  # -----------------------------------------------------------------------
  heterogeneity = list(
    dgps = "dgp_mp",
    param_grid = data.frame(
      M             = 50,
      N             = 50,
      endogenous    = TRUE,
      heterogeneity = TRUE,
      pi_Z          = 0.3,
      beta_slope    = c(0.0, 0.05, 0.1, 0.15, 0.2, 0.3),
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  )
)
