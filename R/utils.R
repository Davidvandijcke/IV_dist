# =============================================================================
# D-IV Estimator: Core Utility Functions
# =============================================================================
#
# Required packages: stats (base)
#
# Functions:
#   pava()                 - Pool-Adjacent-Violators (isotonic regression)
#   compute_group_quantiles() - Empirical quantile functions per group
#   compute_2sls_coefs()   - Unconstrained 2SLS coefficient functions
#
# =============================================================================


#' Pool-Adjacent-Violators Algorithm for isotonic regression
#'
#' Enforces monotone non-decreasing order on a numeric vector using
#' weighted least-squares isotonic regression via isoreg().
#'
#' @param y Numeric vector to project onto monotone non-decreasing cone.
#' @return Numeric vector of same length, monotone non-decreasing.
pava <- function(y) {
  n <- length(y)
  if (n <= 1) return(y)
  if (any(is.na(y))) return(rep(NA_real_, n))
  isoreg(seq_len(n), y)$yf
}


#' Compute empirical quantile functions for each group
#'
#' Given a list of individual-level outcome vectors (one per group), computes
#' the empirical quantile function evaluated at a common grid.
#'
#' @param y_list List of length M, each element a numeric vector of individual
#'   outcomes for that group.
#' @param q_grid Numeric vector of quantile probabilities in (0, 1).
#' @return M x Q matrix where row j is the empirical quantile function for
#'   group j, evaluated at q_grid.
compute_group_quantiles <- function(y_list, q_grid) {
  M <- length(y_list)
  Q <- length(q_grid)
  Q_mat <- matrix(NA_real_, nrow = M, ncol = Q)

  for (j in seq_len(M)) {
    y_j <- y_list[[j]]
    if (length(y_j) == 0 || all(is.na(y_j))) {
      Q_mat[j, ] <- NA_real_
    } else {
      Q_mat[j, ] <- quantile(y_j, probs = q_grid, type = 7, names = FALSE)
    }
  }

  Q_mat
}


#' Compute unconstrained 2SLS coefficient functions
#'
#' Estimates beta_0(u) and beta_1(u) via two-stage least squares applied to
#' group-level quantile functions. The model is:
#'
#'   Q_{Y_j}(u) = beta_0(u) + beta_1(u)' (X_j - mu_X) + epsilon_j(u)
#'   E[epsilon_j(u) | Z_j] = 0
#'
#' Handles both scalar X (M-vector) and multivariate X (M x p matrix).
#' For multivariate X with a single instrument, the first column is treated
#' as endogenous and remaining columns as exogenous controls (partialled out).
#'
#' @param Q_Yk M x Q matrix of group quantile functions.
#' @param X M-vector (scalar case) or M x p matrix (multivariate case) of
#'   endogenous/exogenous covariates.
#' @param Z M-vector or M x l matrix of instruments.
#' @param endog_idx Integer index of the endogenous variable in X when X is a
#'   matrix. Default: 1. Ignored when X is a vector.
#' @return Named list:
#'   \item{beta0}{Q-vector of intercept coefficients (mean quantile function).}
#'   \item{beta1}{Q-vector (scalar X) or p x Q matrix (multivariate X) of
#'     slope coefficients.}
#'   \item{beta_endog}{Q-vector of the endogenous variable's coefficients
#'     (same as beta1 for scalar X).}
compute_2sls_coefs <- function(Q_Yk, X, Z, endog_idx = 1L) {
  M <- nrow(Q_Yk)
  Q <- ncol(Q_Yk)

  # Intercept: always the column mean of quantile functions
  beta0 <- colMeans(Q_Yk)

  # Detect scalar vs matrix X
  if (is.null(dim(X))) {
    # ------ Scalar X, scalar Z ------
    X_centered <- X - mean(X)
    Z_centered <- Z - mean(Z)

    # First stage: X_hat = proj_Z(X) (demeaned)
    Sigma_ZZ <- mean(Z_centered^2)
    Sigma_ZX <- mean(Z_centered * X_centered)

    if (abs(Sigma_ZX) < 1e-10) {
      return(list(beta0 = beta0,
                  beta1 = rep(NA_real_, Q),
                  beta_endog = rep(NA_real_, Q)))
    }

    X_hat <- Z_centered * (Sigma_ZX / Sigma_ZZ)
    denom <- mean(X_hat * X_centered)  # = Sigma_ZX^2 / Sigma_ZZ

    # Second stage: beta_1(u) = cov(X_hat, Q_Y(u)) / cov(X_hat, X)
    beta1 <- as.vector((t(Q_Yk) %*% X_hat / M) / denom)

    return(list(beta0 = beta0, beta1 = beta1, beta_endog = beta1))
  }

  # ------ Matrix X (M x p) ------
  X <- as.matrix(X)
  Z <- if (is.null(dim(Z))) as.matrix(Z) else as.matrix(Z)
  p <- ncol(X)
  l <- ncol(Z)

  if (p == 1) {
    # Single column matrix: treat as scalar
    return(compute_2sls_coefs(Q_Yk, X[, 1], Z[, 1], endog_idx = 1L))
  }

  # Partition into endogenous and exogenous
  X_endog <- X[, endog_idx]
  X_exog <- X[, -endog_idx, drop = FALSE]

  # Build exogenous projection matrix (with intercept)
  W <- cbind(1, X_exog)  # M x (1 + p_exog)
  WtW_inv <- solve(crossprod(W))
  P_W <- W %*% WtW_inv %*% t(W)
  M_W <- diag(M) - P_W  # Annihilator

  # Partial out exogenous variables
  X_endog_resid <- as.vector(M_W %*% X_endog)
  Z_resid <- M_W %*% Z  # M x l

  # First stage on residualized variables
  Sigma_ZZ_resid <- crossprod(Z_resid) / M
  Sigma_ZX_resid <- as.vector(crossprod(Z_resid, X_endog_resid) / M)

  # Check instrument relevance
  if (all(abs(Sigma_ZX_resid) < 1e-10)) {
    return(list(beta0 = beta0,
                beta1 = matrix(NA_real_, nrow = p, ncol = Q),
                beta_endog = rep(NA_real_, Q)))
  }

  # 2SLS for endogenous coefficient
  # X_hat_resid = Z_resid * (Z'Z)^{-1} Z'X_resid (projected)
  pi_hat <- solve(Sigma_ZZ_resid, Sigma_ZX_resid)  # l-vector
  X_hat_resid <- as.vector(Z_resid %*% pi_hat)

  denom_endog <- mean(X_hat_resid * X_endog_resid)
  if (abs(denom_endog) < 1e-10) {
    return(list(beta0 = beta0,
                beta1 = matrix(NA_real_, nrow = p, ncol = Q),
                beta_endog = rep(NA_real_, Q)))
  }

  # Second stage for endogenous: beta_endog(u)
  beta_endog <- numeric(Q)
  for (iq in seq_len(Q)) {
    Q_resid <- as.vector(M_W %*% Q_Yk[, iq])
    beta_endog[iq] <- mean(X_hat_resid * Q_resid) / denom_endog
  }

  # OLS for exogenous coefficients (after removing endogenous effect)
  p_exog <- ncol(X_exog)
  X_exog_centered <- scale(X_exog, center = TRUE, scale = FALSE)
  beta_exog <- matrix(NA_real_, nrow = p_exog, ncol = Q)

  for (iq in seq_len(Q)) {
    # Remove intercept and endogenous effect
    resid_iq <- Q_Yk[, iq] - beta0[iq] - beta_endog[iq] * (X_endog - mean(X_endog))
    # OLS of residual on centered exogenous
    beta_exog[, iq] <- solve(crossprod(X_exog_centered), crossprod(X_exog_centered, resid_iq))
  }

  # Assemble full beta1 matrix (p x Q)
  beta1 <- matrix(NA_real_, nrow = p, ncol = Q)
  beta1[endog_idx, ] <- beta_endog
  exog_indices <- seq_len(p)[-endog_idx]
  for (j in seq_along(exog_indices)) {
    beta1[exog_indices[j], ] <- beta_exog[j, ]
  }

  list(beta0 = beta0, beta1 = beta1, beta_endog = beta_endog)
}
