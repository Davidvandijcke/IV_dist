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
