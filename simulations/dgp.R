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
#   beta_slope    - additional heterogeneity in gamma: gamma(u) = sqrt(u) + beta_slope * Phi^{-1}(u)
#
# Returns a standardized list:
#   y_list     - M-length list of N-vectors
#   x1_list    - M-length list of N-vectors (individual covariates)
#   x2         - M-vector (group treatment)
#   z          - M-vector (instrument)
#   q_grid     - quantile grid
#   true_gamma - true gamma(u) at each quantile
#   dgp_name   - descriptive name
#   params     - all parameters
#
# Base R only. No package dependencies.
# =============================================================================


#' Unified Melly-Pons DGP
#'
#' Generates grouped data with individual-level outcomes following the
#' quantile model of Melly & Pons (2025, Section 4.2).
#'
#' @param M             Number of groups
#' @param N             Number of individuals per group
#' @param q_grid        Quantile grid in (0,1)
#' @param endogenous    If TRUE, x2 is endogenous via shared eta (DGP 3).
#'                      If FALSE, x2 is exogenous (DGP 1 or 2).
#' @param heterogeneity If TRUE, group effects alpha_j(u) are present (DGP 2/3).
#'                      If FALSE, alpha_j = 0 (DGP 1).
#' @param pi_Z          Instrument strength in (0, 1]. Controls the coefficient
#'                      on z in the treatment equation. pi_Z = 1 gives the
#'                      original MP specification.
#' @param beta_slope    Additional heterogeneity in treatment effect across
#'                      quantiles. gamma(u) = sqrt(u) + beta_slope * Phi^{-1}(u).
#'                      beta_slope = 0 gives the original sqrt(u).
#' @return Standardized DGP list
dgp_mp <- function(M, N, q_grid,
                   endogenous    = TRUE,
                   heterogeneity = TRUE,
                   pi_Z          = 1.0,
                   beta_slope    = 0.0) {

  # --- Group heterogeneity ---
  if (heterogeneity) {
    eta <- runif(M)  # eta_j ~ U(0,1)
  } else {
    eta <- rep(0.5, M)  # alpha_j(u) = u*0.5 - u/2 = 0
  }

  # --- Instrument and treatment ---
  z <- exp(0.25 * rnorm(M))  # log-normal instrument

  if (endogenous) {
    # x2 depends on eta (endogeneity source)
    nu <- exp(0.25 * rnorm(M))
    x2 <- pi_Z * z + eta + sqrt(max(0, 1 - pi_Z^2)) * nu
  } else {
    # x2 independent of eta (exogenous)
    x2 <- exp(0.25 * rnorm(M))
    z  <- x2  # instrument = treatment when exogenous
  }

  # --- True treatment effect function ---
  # gamma(u) = sqrt(u) + beta_slope * Phi^{-1}(u)
  true_gamma <- sqrt(q_grid) + beta_slope * qnorm(q_grid)

  # --- Generate individual outcomes ---
  y_list  <- vector("list", M)
  x1_list <- vector("list", M)

  for (j in seq_len(M)) {
    u_ij  <- runif(N)
    x1_ij <- exp(0.25 * rnorm(N))  # log-normal individual covariate

    # Treatment effect at individual ranks
    gamma_u <- sqrt(u_ij) + beta_slope * qnorm(u_ij)

    # Group heterogeneity (mean zero by construction)
    alpha_j <- u_ij * eta[j] - u_ij / 2

    # Outcome
    y_ij <- x1_ij * (u_ij / 2) +   # individual covariate effect: beta_1(u) = u/2
            x2[j] * gamma_u +        # group treatment effect
            alpha_j                   # group heterogeneity

    y_list[[j]]  <- y_ij
    x1_list[[j]] <- x1_ij
  }

  # --- Build descriptive name ---
  parts <- "MP"
  if (!heterogeneity)       parts <- paste0(parts, "_nohetero")
  if (!endogenous)          parts <- paste0(parts, "_exog")
  if (pi_Z < 1)             parts <- paste0(parts, sprintf("_piZ%.2f", pi_Z))
  if (beta_slope != 0)      parts <- paste0(parts, sprintf("_slope%.1f", beta_slope))

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
                      beta_slope = beta_slope)
  )
}
