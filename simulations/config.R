# =============================================================================
# Simulation Study Configuration
# =============================================================================
#
# Defines all experiments as a list. Each experiment specifies:
#   dgps        - character vector of DGP function names
#   param_grid  - data.frame of parameter combinations to sweep
#   estimators  - character vector of estimator names
#   q_grid      - quantile grid
#   n_reps      - number of Monte Carlo replications
#
# =============================================================================

# --- Quantile grids ---
Q_GRID_DEFAULT <- seq(0.1, 0.9, by = 0.05)  # 17 points
Q_GRID_SPARSE  <- c(0.1, 0.5, 0.9)          # 3 points (for MP comparison)

# --- Estimator dispatch table ---
# Maps short names to function names defined in R/estimators.R
ESTIMATOR_REGISTRY <- list(
  "2sls"         = "estimate_2sls",
  "div"          = "estimate_div",
  "div_endpoint" = "estimate_div_endpoint"
)

# --- Experiments ---
EXPERIMENTS <- list(

  # A. Validate against Melly-Pons paper (DGPs 1-3, sparse quantiles)
  mp_validation = list(
    dgps       = c("dgp_mp1", "dgp_mp2", "dgp_mp3"),
    param_grid = expand.grid(M = c(25, 200), N = c(25, 200)),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID_SPARSE,
    n_reps     = 1000
  ),

  # B. IV strength study (main paper result)
  iv_strength = list(
    dgps       = "dgp_mp3_weak_iv",
    param_grid = data.frame(M = 30, N = 25,
                            pi_Z = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.5)),
    estimators = c("2sls", "div", "div_endpoint"),
    q_grid     = Q_GRID_DEFAULT,
    n_reps     = 500
  ),

  # C. Sample size study
  sample_size = list(
    dgps       = "dgp_mp3",
    param_grid = expand.grid(M = c(15, 25, 50, 100),
                             N = c(10, 25, 50, 100)),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID_DEFAULT,
    n_reps     = 500
  ),

  # D. Treatment effect heterogeneity
  heterogeneity = list(
    dgps       = "dgp_mp3_hetero",
    param_grid = data.frame(M = 30, N = 25, pi_Z = 0.2,
                            beta_slope = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID_DEFAULT,
    n_reps     = 500
  ),

  # E. Reproduce current paper tables -- IV strength sweep
  paper_v6_iv = list(
    dgps       = "dgp_adhoc",
    param_grid = data.frame(M = 30, N = 25,
                            pi_Z = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.5),
                            beta_slope = 0.3),
    estimators = c("2sls", "div", "div_endpoint"),
    q_grid     = Q_GRID_DEFAULT,
    n_reps     = 500
  ),

  # F. Reproduce current paper tables -- sample size sweep
  paper_v6_sample = list(
    dgps       = "dgp_adhoc",
    param_grid = expand.grid(M = c(15, 25, 50, 100),
                             N = c(10, 25, 50, 100),
                             pi_Z = 0.25,
                             beta_slope = 0.3),
    estimators = c("2sls", "div"),
    q_grid     = Q_GRID_DEFAULT,
    n_reps     = 500
  )
)
