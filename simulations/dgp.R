# =============================================================================
# DGP for D-IV Simulation Study
# =============================================================================
#
# Single unified DGP based on Melly & Pons (2025) / Chetverikov et al. (2016).
#
# Model (Skorohod representation):
#   y_ij = beta0(u_ij) + x_{1ij} * (u_ij/2) + x_{2j} * gamma(u_ij) + alpha_j(u_ij)
#
# where u_ij ~ U(0,1) is the within-group rank.
#
# Parameters:
#   endogenous    - whether x2 is endogenous (DGP 3) or exogenous (DGP 1/2)
#   heterogeneity - whether group effects alpha_j are present (DGP 2/3) or not
#   pi_Z          - instrument strength (1 = original MP)
#   beta_slope    - heterogeneity in gamma: gamma(u) = sqrt(u) + beta_slope*Phi^{-1}(u)
#   base_dist     - base distribution for beta0(u):
#                   "normal" = Phi^{-1}(u) (original MP)
#                   "t3"     = qt(u, df=3)
#                   "t2"     = qt(u, df=2)
#                   "lognormal" = qlnorm(u, 0, 1.5) - mean (skewed, heavy right tail)
#
# The Skorohod map h_j(u) = beta0(u) + ... is always monotone in u (since
# beta0 is monotone and gamma = sqrt(u) + beta_slope*Phi^{-1}(u) is monotone
# for moderate beta_slope). So u_ij IS the within-group rank, the group
# quantile function IS h_j(u), and gamma(u) IS the true population coefficient.
#
# Base R only. No package dependencies.
# =============================================================================


#' Unified Melly-Pons DGP
#'
#' @param M             Number of groups
#' @param N             Number of individuals per group
#' @param q_grid        Quantile grid in (0,1)
#' @param endogenous    If TRUE, x2 is endogenous via shared eta.
#' @param heterogeneity If TRUE, group effects alpha_j(u) are present.
#' @param pi_Z          Instrument strength in (0, 1].
#' @param beta_slope    gamma(u) = sqrt(u) + beta_slope * Phi^{-1}(u).
#' @param base_dist     Base distribution: "normal", "t3", "t2", or "lognormal".
#' @return Standardized DGP list
dgp_mp <- function(M, N, q_grid,
                   endogenous    = TRUE,
                   heterogeneity = TRUE,
                   pi_Z          = 1.0,
                   beta_slope    = 0.0,
                   base_dist     = "normal") {

  # --- Base distribution quantile function ---
  # "normal" = no added beta0 (matches original MP DGP exactly)
  # Others add a heavy-tailed component to the group QF
  beta0_fn <- switch(base_dist,
    normal    = function(u) 0,
    t3        = function(u) qt(u, df = 3),
    t2        = function(u) qt(u, df = 2),
    lognormal = function(u) qlnorm(u, meanlog = 0, sdlog = 1.5) - exp(1.5^2 / 2),
    stop("Unknown base_dist: ", base_dist)
  )

  # --- Group heterogeneity ---
  if (heterogeneity) {
    eta <- runif(M)
  } else {
    eta <- rep(0.5, M)
  }

  # --- Instrument and treatment ---
  z <- exp(0.25 * rnorm(M))

  if (endogenous) {
    nu <- exp(0.25 * rnorm(M))
    x2 <- pi_Z * z + eta + sqrt(max(0, 1 - pi_Z^2)) * nu
  } else {
    x2 <- exp(0.25 * rnorm(M))
    z  <- x2
  }

  # --- True treatment effect ---
  true_gamma <- sqrt(q_grid) + beta_slope * qnorm(q_grid)

  # --- Generate individual outcomes via Skorohod representation ---
  y_list  <- vector("list", M)
  x1_list <- vector("list", M)
  mu_x <- mean(x2)  # needed for centering in the structural model

  for (j in seq_len(M)) {
    u_ij  <- runif(N)
    x1_ij <- exp(0.25 * rnorm(N))

    gamma_u <- sqrt(u_ij) + beta_slope * qnorm(u_ij)
    alpha_j <- u_ij * eta[j] - u_ij / 2

    y_ij <- beta0_fn(u_ij) +          # base distribution
            x1_ij * (u_ij / 2) +      # individual covariate effect
            x2[j] * gamma_u +         # group treatment effect
            alpha_j                    # group heterogeneity

    y_list[[j]]  <- y_ij
    x1_list[[j]] <- x1_ij
  }

  # --- Descriptive name ---
  parts <- "MP"
  if (!heterogeneity)       parts <- paste0(parts, "_nohetero")
  if (!endogenous)          parts <- paste0(parts, "_exog")
  if (pi_Z < 1)             parts <- paste0(parts, sprintf("_piZ%.2f", pi_Z))
  if (beta_slope != 0)      parts <- paste0(parts, sprintf("_slope%.1f", beta_slope))
  if (base_dist != "normal") parts <- paste0(parts, "_", base_dist)

  list(
    y_list     = y_list,
    x1_list    = x1_list,
    x2         = x2,
    z          = z,
    q_grid     = q_grid,
    true_gamma = true_gamma,
    dgp_name   = parts,
    params     = list(M = M, N = N, endogenous = endogenous,
                      heterogeneity = heterogeneity, pi_Z = pi_Z,
                      beta_slope = beta_slope, base_dist = base_dist)
  )
}


# =============================================================================
# Centered-x₂ DGP: group-level model with normal instruments
# =============================================================================
#
# Model (Skorohod representation):
#   y_ij = sigma * QF_base(u_ij) + x_{2j} * gamma(u_ij) + alpha_j(u_ij)
#
# where:
#   QF_base = Phi^{-1} (normal), qt(df=5) (t_5), qt(df=3) (t_3)
#   gamma(u) = u + beta_slope * Phi^{-1}(u)  [bounded derivative gamma'=1+...]
#   alpha_j(u) = u * eta_j - u/2,  eta_j ~ U(0,1)
#   x_{2j} = (pi_Z + delta*(eta_j-0.5)) * z_j + eta_j + sqrt(1-pi_Z^2) * nu_j
#   z_j, nu_j ~ N(0,1)  [CENTERED — x₂ can be negative]
#
# Optional controls: W_k ~ N(0,1) with effects u/(k+2)
#
# The key differences from dgp_mp:
#   1. Normal (not lognormal) instruments → centered x₂ with negative values
#   2. gamma(u) = u (bounded derivative), not sqrt(u) (unbounded near 0)
#   3. sigma * QF_base provides the intercept (tunable steepness)
#   4. Heterogeneous first stage via delta parameter
#   5. No individual covariates x1 (cleaner group-level model)
#
# Returns eta vector for W₂² computation.
# =============================================================================

#' Centered-x₂ DGP with normal instruments
#'
#' @param M             Number of groups
#' @param N             Number of individuals per group
#' @param q_grid        Quantile grid in (0,1)
#' @param pi_Z          Instrument strength in (0, 1].
#' @param sigma_beta0   Scale of base QF (steepness of intercept, default 2).
#' @param base_dist     "normal", "t5", or "t3".
#' @param hetero_fs     First-stage heterogeneity delta.
#'                      x2 = (pi_Z + delta*(eta-0.5))*z + eta + noise.
#' @param p             Number of regressors (1 = endogenous only, >1 adds controls).
#' @param beta_slope    gamma(u) = u + beta_slope * Phi^{-1}(u).
#' @return Standardized DGP list (includes eta for W₂² computation)
dgp_centered <- function(M, N, q_grid,
                          pi_Z         = 0.5,
                          sigma_beta0  = 3,
                          base_dist    = "normal",
                          hetero_fs    = 0,
                          p            = 1L,
                          beta_slope   = 0) {

  beta0_fn <- switch(base_dist,
    normal      = function(u) sigma_beta0 * qnorm(u),
    t5          = function(u) sigma_beta0 * qt(u, df = 5),
    t3          = function(u) sigma_beta0 * qt(u, df = 3),
    lognormal   = function(u) sigma_beta0 * (qlnorm(u, 0, 1) - exp(0.5)),
    exponential = function(u) sigma_beta0 * (qexp(u, 1) - 1),
    chisq3      = function(u) sigma_beta0 * (qchisq(u, 3) - 3),
    stop("Unknown base_dist: ", base_dist)
  )

  eta <- runif(M)
  z_raw <- rnorm(M)
  nu <- rnorm(M)
  x_endog <- (pi_Z + hetero_fs * (eta - 0.5)) * z_raw +
             eta + sqrt(max(0, 1 - pi_Z^2)) * nu

  # Controls (if p > 1)
  if (p > 1L) {
    W <- matrix(rnorm(M * (p - 1L)), M, p - 1L)
    X <- cbind(x_endog, W)
    Z <- cbind(z_raw, W)  # controls are own instruments
  } else {
    X <- x_endog
    Z <- z_raw
  }

  # Treatment effect: gamma(u) = u + beta_slope * sin(2*pi*u)
  # The sine perturbation makes gamma non-monotone for beta_slope > 0.
  # gamma' = 1 + 2*pi*beta_slope*cos(2*pi*u) > 0 for beta_slope < 1/(2*pi) ~ 0.16,
  # but the QF stays monotone even for larger beta_slope because sigma*Phi^{-1}'
  # dominates. For beta_slope = 0 this reduces to gamma(u) = u.
  true_gamma <- q_grid + beta_slope * sin(2 * pi * q_grid)

  # Control effect functions: odd controls have monotone effects (u),
  # even controls have non-monotone effects (0.1 * sin(k*pi*u)).
  # This mix is more realistic than all-monotone controls.
  ctrl_fn <- function(u, k) {
    if (k %% 2 == 1) u else 0.1 * sin(k * pi * u)
  }

  # True intercept: E[Q_{Y_j}(u)]
  mu_x <- mean(x_endog)
  true_beta0 <- beta0_fn(q_grid) + mu_x * true_gamma
  if (p > 1L) {
    for (k in seq_len(p - 1L)) {
      true_beta0 <- true_beta0 + mean(W[, k]) * ctrl_fn(q_grid, k)
    }
  }

  # alpha_j(u) at grid points (for W₂² computation)
  alpha_mat <- outer(eta, q_grid) - outer(rep(0.5, M), q_grid)  # M x Q

  y_list <- vector("list", M)
  for (j in seq_len(M)) {
    u_ij <- runif(N)
    gamma_u <- u_ij + beta_slope * sin(2 * pi * u_ij)
    alpha_u <- u_ij * eta[j] - u_ij / 2
    y_ij <- beta0_fn(u_ij) + x_endog[j] * gamma_u + alpha_u
    if (p > 1L) {
      for (k in seq_len(p - 1L)) {
        y_ij <- y_ij + W[j, k] * ctrl_fn(u_ij, k)
      }
    }
    y_list[[j]] <- y_ij
  }

  # Descriptive name
  parts <- "centered"
  if (pi_Z != 0.5)     parts <- paste0(parts, sprintf("_piZ%.2f", pi_Z))
  if (hetero_fs != 0)   parts <- paste0(parts, sprintf("_hfs%.1f", hetero_fs))
  if (p > 1)            parts <- paste0(parts, sprintf("_p%d", p))
  if (beta_slope != 0)  parts <- paste0(parts, sprintf("_slope%.1f", beta_slope))
  if (base_dist != "normal") parts <- paste0(parts, "_", base_dist)

  list(
    y_list     = y_list,
    x1_list    = NULL,
    x2         = X,
    z          = Z,
    q_grid     = q_grid,
    true_gamma = true_gamma,
    true_beta0 = true_beta0,
    eta        = eta,
    alpha_mat  = alpha_mat,
    dgp_name   = parts,
    params     = list(M = M, N = N, pi_Z = pi_Z, sigma_beta0 = sigma_beta0,
                      base_dist = base_dist, hetero_fs = hetero_fs,
                      p = p, beta_slope = beta_slope)
  )
}


# =============================================================================
# DGP variant with monotone control effects of varying magnitude
# =============================================================================
#
# Same as dgp_centered but control effects are gamma_k(u) = c_k * u
# where c_k ~ N(0, ctrl_magnitude) rather than u/(k+2).
# This produces MONOTONE population QFs but NOISY estimated coefficient
# functions that oscillate, driving violations in the fitted psi curves.
# Matches the pattern seen in the CLP empirical application.
#
#' @param ctrl_magnitude SD of control effect coefficients c_k. Larger values
#'   produce noisier estimated control coefficients and more violations.
dgp_centered_ctrl <- function(M, N, q_grid,
                               pi_Z          = 0.5,
                               sigma_beta0   = 3,
                               base_dist     = "normal",
                               hetero_fs     = 0,
                               p             = 1L,
                               beta_slope    = 0,
                               ctrl_magnitude = 1.0) {

  beta0_fn <- switch(base_dist,
    normal      = function(u) sigma_beta0 * qnorm(u),
    t5          = function(u) sigma_beta0 * qt(u, df = 5),
    t3          = function(u) sigma_beta0 * qt(u, df = 3),
    lognormal   = function(u) sigma_beta0 * (qlnorm(u, 0, 1) - exp(0.5)),
    exponential = function(u) sigma_beta0 * (qexp(u, 1) - 1),
    chisq3      = function(u) sigma_beta0 * (qchisq(u, 3) - 3),
    stop("Unknown base_dist: ", base_dist)
  )

  eta <- runif(M)
  z_raw <- rnorm(M)
  nu <- rnorm(M)
  x_endog <- (pi_Z + hetero_fs * (eta - 0.5)) * z_raw +
             eta + sqrt(max(0, 1 - pi_Z^2)) * nu

  # Controls with random monotone effects
  if (p > 1L) {
    n_ctrl <- p - 1L
    W <- matrix(rnorm(M * n_ctrl), M, n_ctrl)
    X <- cbind(x_endog, W)
    Z <- cbind(z_raw, W)
    # Control effect coefficients with magnitude scaling
    c_k <- rnorm(n_ctrl, 0, ctrl_magnitude)
    # Control effect functions: odd = monotone (u), even = non-monotone (0.1*sin)
    ctrl_fn_c <- function(u, k) {
      if (k %% 2 == 1) u else 0.1 * sin(k * pi * u)
    }
  } else {
    X <- x_endog
    Z <- z_raw
    c_k <- numeric(0)
  }

  true_gamma <- q_grid + beta_slope * sin(2 * pi * q_grid)

  # True intercept
  mu_x <- mean(x_endog)
  true_beta0 <- beta0_fn(q_grid) + mu_x * true_gamma
  if (p > 1L) {
    for (k in seq_len(p - 1L)) {
      true_beta0 <- true_beta0 + mean(W[, k]) * c_k[k] * ctrl_fn_c(q_grid, k)
    }
  }

  # alpha at grid points
  alpha_mat <- outer(eta, q_grid) - outer(rep(0.5, M), q_grid)

  y_list <- vector("list", M)
  for (j in seq_len(M)) {
    u_ij <- runif(N)
    gamma_u <- u_ij + beta_slope * sin(2 * pi * u_ij)
    alpha_u <- u_ij * eta[j] - u_ij / 2
    y_ij <- beta0_fn(u_ij) + x_endog[j] * gamma_u + alpha_u
    if (p > 1L) {
      for (k in seq_len(p - 1L)) {
        y_ij <- y_ij + W[j, k] * c_k[k] * ctrl_fn_c(u_ij, k)
      }
    }
    y_list[[j]] <- y_ij
  }

  parts <- "centered_ctrl"
  if (pi_Z != 0.5)        parts <- paste0(parts, sprintf("_piZ%.2f", pi_Z))
  if (hetero_fs != 0)      parts <- paste0(parts, sprintf("_hfs%.1f", hetero_fs))
  if (p > 1)               parts <- paste0(parts, sprintf("_p%d", p))
  if (ctrl_magnitude != 1) parts <- paste0(parts, sprintf("_cm%.1f", ctrl_magnitude))

  list(
    y_list     = y_list,
    x1_list    = NULL,
    x2         = X,
    z          = Z,
    q_grid     = q_grid,
    true_gamma = true_gamma,
    true_beta0 = true_beta0,
    eta        = eta,
    alpha_mat  = alpha_mat,
    dgp_name   = parts,
    params     = list(M = M, N = N, pi_Z = pi_Z, sigma_beta0 = sigma_beta0,
                      base_dist = base_dist, hetero_fs = hetero_fs,
                      p = p, beta_slope = beta_slope, ctrl_magnitude = ctrl_magnitude)
  )
}
