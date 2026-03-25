# =============================================================================
# Simulation Configuration
# =============================================================================
#
# Experiments using both the MP DGP (dgp_mp) and centered-x₂ DGP (dgp_centered).
#
# =============================================================================

Q_GRID <- seq(0.05, 0.95, by = 0.05)   # 19 points

ESTIMATOR_REGISTRY <- list(
  "2sls" = "estimate_2sls",
  "div"  = "estimate_div"
)

EXPERIMENTS <- list(

  # -----------------------------------------------------------------------
  # 1. IV strength: centered x₂, vary pi_Z
  # -----------------------------------------------------------------------
  iv_strength = list(
    dgps = "dgp_centered",
    param_grid = data.frame(
      M    = 50, N = 50,
      pi_Z = c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0),
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  ),

  # -----------------------------------------------------------------------
  # 2. Sample size: centered x₂, F≈12
  # -----------------------------------------------------------------------
  sample_size = list(
    dgps = "dgp_centered",
    param_grid = expand.grid(
      M    = c(25, 50, 100, 200),
      N    = c(25, 50, 100, 200),
      pi_Z = 0.5,
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  ),

  # -----------------------------------------------------------------------
  # 3. Controls + heterogeneous first stage
  # -----------------------------------------------------------------------
  controls = list(
    dgps = "dgp_centered",
    param_grid = expand.grid(
      M         = 50, N = 50,
      pi_Z      = 0.5,
      p         = c(1, 2, 3, 5),
      hetero_fs = c(0, 0.5, 1.0),
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  ),

  # -----------------------------------------------------------------------
  # 4. Heavy tails (centered x₂)
  # -----------------------------------------------------------------------
  heavy_tails = list(
    dgps = "dgp_centered",
    param_grid = expand.grid(
      M         = 50,
      N         = c(25, 50),
      pi_Z      = 0.5,
      base_dist = c("normal", "t5", "t3"),
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  ),

  # -----------------------------------------------------------------------
  # 5. Realistic combination
  #    p=3, hetero_fs=0.5, t5 base, beta_slope=0.2, F≈12
  # -----------------------------------------------------------------------
  realistic = list(
    dgps = "dgp_centered",
    param_grid = data.frame(
      M          = 50, N = 50,
      pi_Z       = 0.5,
      p          = 3,
      hetero_fs  = 0.5,
      base_dist  = "t5",
      beta_slope = 0.2,
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  ),

  # -----------------------------------------------------------------------
  # 6. MP baseline comparison (pi_Z=1, original DGP)
  # -----------------------------------------------------------------------
  mp_baseline = list(
    dgps = "dgp_mp",
    param_grid = expand.grid(
      M             = c(25, 50, 100),
      N             = c(25, 50),
      endogenous    = TRUE,
      heterogeneity = TRUE,
      pi_Z          = 1.0,
      stringsAsFactors = FALSE
    ),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID,
    n_reps     = 500
  )
)
