# =============================================================================
# DGP for D-IV Simulation Study
# =============================================================================
#
# Single unified DGP based on Melly & Pons (2025) / Chetverikov et al. (2016).
# All simulation experiments are parameterized variants of this DGP.
#
# Base model (Skorohod representation):
#   y_ij = x_{1ij} * beta_1(u_ij) + x_{2j} * gamma(u_ij) + alpha_j(u_ij)
#
# where u_ij ~ U(0,1), x_{1ij} ~ Lognormal, x_{2j} group-level treatment.
#
# Parameters control:
#   endogenous    - whether x2 is endogenous (DGP 3) or exogenous (DGP 1/2)
#   heterogeneity - whether group effects alpha_j are present (DGP 2/3) or not (DGP 1)
#   pi_Z          - instrument strength (1 = original, <1 = weaker)
#   beta_slope    - additional heterogeneity in gamma
#   error_df      - degrees of freedom for t-distributed individual-level noise
#                   (Inf = no added noise, matching the original MP DGP)
#   error_scale   - scale of the added noise
#
# Base R only. No package dependencies.
# =============================================================================


#' Unified Melly-Pons DGP
#'
#' @param M             Number of groups
#' @param N             Number of individuals per group
#' @param q_grid        Quantile grid in (0,1)
#' @param endogenous    If TRUE, x2 is endogenous via shared eta (DGP 3).
#' @param heterogeneity If TRUE, group effects alpha_j(u) are present (DGP 2/3).
#' @param pi_Z          Instrument strength in (0, 1].
#' @param beta_slope    gamma(u) = sqrt(u) + beta_slope * Phi^{-1}(u).
#' @param error_df      Degrees of freedom for added t-distributed noise.
#'                      Inf = no added noise (original MP). 2-3 = heavy tails.
#' @param error_scale   Scale multiplier for the added noise (default 1).
#' @return Standardized DGP list
dgp_mp <- function(M, N, q_grid,
                   endogenous    = TRUE,
                   heterogeneity = TRUE,
                   pi_Z          = 1.0,
                   beta_slope    = 0.0,
                   error_df      = Inf,
                   error_scale   = 1.0) {

  # --- Group heterogeneity ---
  if (heterogeneity) {
    eta <- runif(M)  # eta_j ~ U(0,1)
  } else {
    eta <- rep(0.5, M)  # alpha_j(u) = u*0.5 - u/2 = 0
  }

  # --- Instrument and treatment ---
  z <- exp(0.25 * rnorm(M))  # log-normal instrument

  if (endogenous) {
    nu <- exp(0.25 * rnorm(M))
    x2 <- pi_Z * z + eta + sqrt(max(0, 1 - pi_Z^2)) * nu
  } else {
    x2 <- exp(0.25 * rnorm(M))
    z  <- x2
  }

  # --- True treatment effect function ---
  true_gamma <- sqrt(q_grid) + beta_slope * qnorm(q_grid)

  # --- Generate individual outcomes ---
  y_list  <- vector("list", M)
  x1_list <- vector("list", M)

  for (j in seq_len(M)) {
    u_ij  <- runif(N)
    x1_ij <- exp(0.25 * rnorm(N))

    gamma_u <- sqrt(u_ij) + beta_slope * qnorm(u_ij)
    alpha_j <- u_ij * eta[j] - u_ij / 2

    y_ij <- x1_ij * (u_ij / 2) +
            x2[j] * gamma_u +
            alpha_j

    # Add heavy-tailed noise if requested
    if (is.finite(error_df)) {
      y_ij <- y_ij + error_scale * rt(N, df = error_df)
    }

    y_list[[j]]  <- y_ij
    x1_list[[j]] <- x1_ij
  }

  # --- Build descriptive name ---
  parts <- "MP"
  if (!heterogeneity)       parts <- paste0(parts, "_nohetero")
  if (!endogenous)          parts <- paste0(parts, "_exog")
  if (pi_Z < 1)             parts <- paste0(parts, sprintf("_piZ%.2f", pi_Z))
  if (beta_slope != 0)      parts <- paste0(parts, sprintf("_slope%.1f", beta_slope))
  if (is.finite(error_df))  parts <- paste0(parts, sprintf("_t%g", error_df))

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
                      beta_slope = beta_slope, error_df = error_df,
                      error_scale = error_scale)
  )
}
