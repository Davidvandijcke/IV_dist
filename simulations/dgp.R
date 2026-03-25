# =============================================================================
# DGP (Data Generating Processes) for D-IV Simulation Study
# =============================================================================
#
# Each DGP returns a standardized list:
#   y_list     - M-length list of N-vectors (individual outcomes)
#   x1_list    - M-length list (individual covariates); NULL if none
#   x2         - M-vector of group-level treatment (scalar)
#   z          - M-vector of group-level instrument (scalar)
#   q_grid     - quantile grid
#   true_gamma - true gamma(u) at each quantile
#   dgp_name   - character string identifying the DGP
#   params     - list of DGP parameters for reference
#
# DGPs:
#   1. dgp_mp1           - Melly-Pons DGP 1: no group heterogeneity
#   2. dgp_mp2           - Melly-Pons DGP 2: exogenous group heterogeneity
#   3. dgp_mp3           - Melly-Pons DGP 3: endogenous treatment
#   4. dgp_adhoc         - Ad-hoc DGP from current paper
#   5. dgp_mp3_weak_iv   - DGP 3 with controllable IV strength
#   6. dgp_mp3_hetero    - DGP 3 with heterogeneous gamma
#
# Base R only. No package dependencies.
# =============================================================================


# =============================================================================
# DGP 1: Melly-Pons DGP 1 -- No group heterogeneity
# =============================================================================
#
# From Melly & Pons (2025) Section 4.2, matching Chetverikov et al. (2016).
#
# Model (Skorohod representation):
#   y_ij = beta_0(u_ij) + x_{1ij} * beta_1(u_ij) + x_{2j} * gamma(u_ij)
#
# where:
#   u_ij ~ U(0,1)
#   x_{1ij} ~ exp(0.25 * N(0,1))     (log-normal individual covariate)
#   x_{2j}  ~ exp(0.25 * N(0,1))     (log-normal group treatment, EXOGENOUS)
#   gamma(u) = sqrt(u)
#   beta_1(u) = u / 2
#   beta_0(u) = 0
#   alpha_j(u) = 0                    (no group heterogeneity)
#   z_j = x_{2j}                      (instrument = treatment)
#
# So: y_ij = x_{1ij} * (u_ij / 2) + x_{2j} * sqrt(u_ij)
#
#' @param M Number of groups
#' @param N Number of individuals per group
#' @param q_grid Quantile grid in (0,1)
#' @return Standardized DGP list
dgp_mp1 <- function(M, N, q_grid) {

  # Group-level treatment (exogenous, log-normal)
  x2 <- exp(0.25 * rnorm(M))

  # Instrument = treatment (exogenous case)
  z <- x2

  # Generate individual outcomes and covariates
  y_list  <- vector("list", M)
  x1_list <- vector("list", M)

  for (j in seq_len(M)) {
    u_ij   <- runif(N)
    x1_ij  <- exp(0.25 * rnorm(N))

    y_ij <- x1_ij * (u_ij / 2) + x2[j] * sqrt(u_ij)

    y_list[[j]]  <- y_ij
    x1_list[[j]] <- x1_ij
  }

  # True gamma: sqrt(u)
  true_gamma <- sqrt(q_grid)

  list(
    y_list     = y_list,
    x1_list    = x1_list,
    x2         = x2,
    z          = z,
    q_grid     = q_grid,
    true_gamma = true_gamma,
    dgp_name   = "MP1_no_heterogeneity",
    params     = list(M = M, N = N)
  )
}


# =============================================================================
# DGP 2: Melly-Pons DGP 2 -- Exogenous group heterogeneity
# =============================================================================
#
# Same as DGP 1, but with group-level heterogeneity:
#   alpha_j(u) = u * eta_j - u / 2,   eta_j ~ U(0,1)
#   E[alpha_j(u)] = u * 0.5 - u/2 = 0  (mean zero)
#   x_{2j} still exogenous (uncorrelated with eta_j)
#
# So: y_ij = x_{1ij} * (u_ij / 2) + x_{2j} * sqrt(u_ij) + u_ij * eta_j - u_ij / 2
#
#' @param M Number of groups
#' @param N Number of individuals per group
#' @param q_grid Quantile grid in (0,1)
#' @return Standardized DGP list
dgp_mp2 <- function(M, N, q_grid) {

  # Group-level treatment (exogenous, log-normal)
  x2 <- exp(0.25 * rnorm(M))

  # Instrument = treatment (exogenous case)
  z <- x2

  # Group heterogeneity
  eta <- runif(M)

  # Generate individual outcomes and covariates
  y_list  <- vector("list", M)
  x1_list <- vector("list", M)

  for (j in seq_len(M)) {
    u_ij  <- runif(N)
    x1_ij <- exp(0.25 * rnorm(N))

    y_ij <- x1_ij * (u_ij / 2) +
            x2[j] * sqrt(u_ij) +
            u_ij * eta[j] - u_ij / 2

    y_list[[j]]  <- y_ij
    x1_list[[j]] <- x1_ij
  }

  # True gamma: sqrt(u)
  true_gamma <- sqrt(q_grid)

  list(
    y_list     = y_list,
    x1_list    = x1_list,
    x2         = x2,
    z          = z,
    q_grid     = q_grid,
    true_gamma = true_gamma,
    dgp_name   = "MP2_exogenous_heterogeneity",
    params     = list(M = M, N = N)
  )
}


# =============================================================================
# DGP 3: Melly-Pons DGP 3 -- Endogenous treatment
# =============================================================================
#
# Same heterogeneity as DGP 2, but x_{2j} is endogenous:
#   eta_j ~ U(0,1)
#   z_j   ~ exp(0.25 * N(0,1))      (log-normal instrument)
#   nu_j  ~ exp(0.25 * N(0,1))      (additional noise)
#   x_{2j} = z_j + eta_j + nu_j     (ENDOGENOUS: eta_j in both x2 and alpha_j)
#   alpha_j(u) = u * eta_j - u / 2
#
# So: y_ij = x_{1ij} * (u_ij / 2) + x_{2j} * sqrt(u_ij) + u_ij * eta_j - u_ij / 2
#
#' @param M Number of groups
#' @param N Number of individuals per group
#' @param q_grid Quantile grid in (0,1)
#' @return Standardized DGP list
dgp_mp3 <- function(M, N, q_grid) {

  # Group heterogeneity (source of endogeneity)
  eta <- runif(M)

  # Instrument (log-normal)
  z <- exp(0.25 * rnorm(M))

  # Additional noise
  nu <- exp(0.25 * rnorm(M))

  # Endogenous treatment: depends on eta
  x2 <- z + eta + nu

  # Generate individual outcomes and covariates
  y_list  <- vector("list", M)
  x1_list <- vector("list", M)

  for (j in seq_len(M)) {
    u_ij  <- runif(N)
    x1_ij <- exp(0.25 * rnorm(N))

    y_ij <- x1_ij * (u_ij / 2) +
            x2[j] * sqrt(u_ij) +
            u_ij * eta[j] - u_ij / 2

    y_list[[j]]  <- y_ij
    x1_list[[j]] <- x1_ij
  }

  # True gamma: sqrt(u)
  true_gamma <- sqrt(q_grid)

  list(
    y_list     = y_list,
    x1_list    = x1_list,
    x2         = x2,
    z          = z,
    q_grid     = q_grid,
    true_gamma = true_gamma,
    dgp_name   = "MP3_endogenous",
    params     = list(M = M, N = N)
  )
}


# =============================================================================
# DGP 4: Ad-hoc DGP from current paper
# =============================================================================
#
# Preserves the existing simulation from the FIVR paper.
# No individual covariates (x1_list = NULL).
#
# Model:
#   Z_j ~ N(0,1)
#   v_j ~ N(0, sigma_v)
#   X_j = pi_Z * Z_j + 0.7 * v_j + N(0, 0.3)
#   beta_0(u) = Phi^{-1}(u)
#   beta_1(u) = beta_base + safe_slope * Phi^{-1}(u)
#     where safe_slope = min(beta_slope, 0.9 / max|X_j - mu_X|)
#   Y_ij = Phi^{-1}(U_ij) + beta_1(U_ij) * (X_j - mu_X) + 0.5 * v_j
#     where U_ij ~ U(0,1)
#
#' @param M Number of groups
#' @param N Number of individuals per group
#' @param q_grid Quantile grid in (0,1)
#' @param pi_Z Instrument strength (default 0.5)
#' @param beta_base Base treatment effect (default 1.0)
#' @param beta_slope Heterogeneity slope (default 0.3)
#' @param sigma_v Confounder std dev (default 0.8)
#' @return Standardized DGP list
dgp_adhoc <- function(M, N, q_grid, pi_Z = 0.5, beta_base = 1.0,
                      beta_slope = 0.3, sigma_v = 0.8) {

  z  <- rnorm(M)
  v  <- rnorm(M, 0, sigma_v)
  x2 <- pi_Z * z + 0.7 * v + rnorm(M, 0, 0.3)

  mu_x    <- mean(x2)
  max_dev <- max(abs(x2 - mu_x))

  # Bound heterogeneity to ensure valid quantile functions
  safe_slope <- min(beta_slope, 0.9 / max_dev)

  # True gamma = beta_1(u)
  true_gamma <- beta_base + safe_slope * qnorm(q_grid)

  # Generate individual outcomes
  y_list <- vector("list", M)

  for (j in seq_len(M)) {
    U_ij <- runif(N)
    x_dev <- x2[j] - mu_x

    y_ij <- qnorm(U_ij) +
            (beta_base + safe_slope * qnorm(U_ij)) * x_dev +
            0.5 * v[j]

    y_list[[j]] <- y_ij
  }

  list(
    y_list     = y_list,
    x1_list    = NULL,
    x2         = x2,
    z          = z,
    q_grid     = q_grid,
    true_gamma = true_gamma,
    dgp_name   = "adhoc",
    params     = list(M = M, N = N, pi_Z = pi_Z, beta_base = beta_base,
                      beta_slope = safe_slope, sigma_v = sigma_v)
  )
}


# =============================================================================
# DGP 5: Weak IV extension of DGP 3
# =============================================================================
#
# Like DGP 3 but with controllable instrument strength:
#   z_j  ~ exp(0.25 * N(0,1))
#   nu_j ~ exp(0.25 * N(0,1))
#   x_{2j} = pi_Z * z_j + eta_j + sqrt(1 - pi_Z^2) * nu_j
#
# pi_Z controls first-stage strength; pi_Z = 1 recovers DGP 3 approximately.
#
#' @param M Number of groups
#' @param N Number of individuals per group
#' @param q_grid Quantile grid in (0,1)
#' @param pi_Z Instrument strength in (0,1] (default 0.5)
#' @return Standardized DGP list
dgp_mp3_weak_iv <- function(M, N, q_grid, pi_Z = 0.5) {

  # Group heterogeneity (source of endogeneity)
  eta <- runif(M)

  # Instrument (log-normal)
  z <- exp(0.25 * rnorm(M))

  # Additional noise (log-normal)
  nu <- exp(0.25 * rnorm(M))

  # Endogenous treatment with controllable IV strength
  x2 <- pi_Z * z + eta + sqrt(1 - pi_Z^2) * nu

  # Generate individual outcomes and covariates
  y_list  <- vector("list", M)
  x1_list <- vector("list", M)

  for (j in seq_len(M)) {
    u_ij  <- runif(N)
    x1_ij <- exp(0.25 * rnorm(N))

    y_ij <- x1_ij * (u_ij / 2) +
            x2[j] * sqrt(u_ij) +
            u_ij * eta[j] - u_ij / 2

    y_list[[j]]  <- y_ij
    x1_list[[j]] <- x1_ij
  }

  # True gamma: sqrt(u)
  true_gamma <- sqrt(q_grid)

  list(
    y_list     = y_list,
    x1_list    = x1_list,
    x2         = x2,
    z          = z,
    q_grid     = q_grid,
    true_gamma = true_gamma,
    dgp_name   = "MP3_weak_iv",
    params     = list(M = M, N = N, pi_Z = pi_Z)
  )
}


# =============================================================================
# DGP 6: Heterogeneous treatment effect extension
# =============================================================================
#
# Like DGP 5 but with heterogeneous gamma:
#   gamma(u) = sqrt(u) + beta_slope * qnorm(u)
#
# beta_slope = 0 recovers standard DGP 3/5.
# Everything else same as DGP 5.
#
#' @param M Number of groups
#' @param N Number of individuals per group
#' @param q_grid Quantile grid in (0,1)
#' @param pi_Z Instrument strength in (0,1] (default 0.5)
#' @param beta_slope Additional heterogeneity in gamma (default 0.3)
#' @return Standardized DGP list
dgp_mp3_hetero <- function(M, N, q_grid, pi_Z = 0.5, beta_slope = 0.3) {

  # Group heterogeneity (source of endogeneity)
  eta <- runif(M)

  # Instrument (log-normal)
  z <- exp(0.25 * rnorm(M))

  # Additional noise (log-normal)
  nu <- exp(0.25 * rnorm(M))

  # Endogenous treatment with controllable IV strength
  x2 <- pi_Z * z + eta + sqrt(1 - pi_Z^2) * nu

  # Generate individual outcomes and covariates
  y_list  <- vector("list", M)
  x1_list <- vector("list", M)

  for (j in seq_len(M)) {
    u_ij  <- runif(N)
    x1_ij <- exp(0.25 * rnorm(N))

    # gamma(u) = sqrt(u) + beta_slope * qnorm(u)
    gamma_u <- sqrt(u_ij) + beta_slope * qnorm(u_ij)

    y_ij <- x1_ij * (u_ij / 2) +
            x2[j] * gamma_u +
            u_ij * eta[j] - u_ij / 2

    y_list[[j]]  <- y_ij
    x1_list[[j]] <- x1_ij
  }

  # True gamma: sqrt(u) + beta_slope * Phi^{-1}(u)
  true_gamma <- sqrt(q_grid) + beta_slope * qnorm(q_grid)

  list(
    y_list     = y_list,
    x1_list    = x1_list,
    x2         = x2,
    z          = z,
    q_grid     = q_grid,
    true_gamma = true_gamma,
    dgp_name   = "MP3_hetero_gamma",
    params     = list(M = M, N = N, pi_Z = pi_Z, beta_slope = beta_slope)
  )
}
