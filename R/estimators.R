# =============================================================================
# D-IV Estimators
# =============================================================================
#
# All estimator functions for the D-IV (Distribution-valued IV) project.
#
# Each estimator takes the same core inputs:
#   Q_Yk  - M x Q matrix of group quantile functions
#   X     - M-vector (scalar) or M x p matrix (multivariate) of covariates
#   Z     - M-vector or M x l matrix of instruments
#   q_grid - Q-vector of quantile levels
#
# Each returns:
#   list(beta0 = Q-vector, beta1 = Q-vector or p x Q matrix, method = string)
#
# Required packages: none (base R only)
# =============================================================================

# Source utility functions (pava).
# When sourced from another script, resolve path relative to this file's location.
local({
  this_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) NULL
  )
  if (is.null(this_dir) || this_dir == "") {
    this_dir <- "R"
  }
  source(file.path(this_dir, "utils.R"), local = FALSE)
})


# =============================================================================
# Estimator 1: Unconstrained 2SLS
# =============================================================================

#' Standard unconstrained 2SLS on group quantile functions.
#'
#' For each quantile u, runs 2SLS with design matrix cbind(1, X - mu_X)
#' and instrument matrix cbind(1, Z - mu_Z). In the just-identified scalar
#' case this reduces to the Wald estimator: beta_1(u) = cov(Z, Q_Y(u)) / cov(Z, X).
#'
#' @param Q_Yk M x Q matrix of group quantile functions
#' @param X     M-vector or M x p matrix of group-level covariates
#' @param Z     M-vector or M x l matrix of instruments
#' @param q_grid Q-vector of quantile levels
#' @param weights Optional M-vector of analytic weights (default: uniform)
#' @return list(beta0, beta1, method) where beta1 is Q-vector (scalar X) or p x Q matrix
estimate_2sls <- function(Q_Yk, X, Z, q_grid, weights = NULL) {
  M  <- nrow(Q_Yk)
  Q  <- length(q_grid)

  # Coerce to matrix form
  X_mat <- if (is.matrix(X)) X else matrix(X, ncol = 1)
  Z_mat <- if (is.matrix(Z)) Z else matrix(Z, ncol = 1)
  p <- ncol(X_mat)
  l <- ncol(Z_mat)

  # Weights: normalize so mean(w) = 1, then scale by sqrt(w)
  if (is.null(weights)) {
    sw <- rep(1, M)
    mu_X <- colMeans(X_mat)
    mu_Z <- colMeans(Z_mat)
  } else {
    w <- weights / mean(weights)
    sw <- sqrt(w)
    mu_X <- colSums(w * X_mat) / sum(w)
    mu_Z <- colSums(w * Z_mat) / sum(w)
  }

  # Centered design and instrument matrices (with intercept)
  X_design <- cbind(1, sweep(X_mat, 2, mu_X))   # M x (p+1)
  Z_design <- cbind(1, sweep(Z_mat, 2, mu_Z))   # M x (l+1)

  # Scale by sqrt(weights) for weighted LS
  X_w <- X_design * sw   # M x (p+1)
  Z_w <- Z_design * sw   # M x (l+1)
  Q_w <- Q_Yk * sw       # M x Q

  # Instrument projection (memory-efficient: avoid forming M x M matrix P_Z)
  # X_hat = Z_w (Z_w'Z_w)^{-1} Z_w' X_w = Z_w %*% (ZtZ_inv %*% Z_w'X_w)
  ZtZ <- crossprod(Z_w)                  # (l+1) x (l+1)
  ZtZ_inv <- tryCatch(solve(ZtZ), error = function(e) NULL)
  if (is.null(ZtZ_inv)) {
    return(list(
      beta0  = rep(NA_real_, Q),
      beta1  = if (p == 1) rep(NA_real_, Q) else matrix(NA_real_, p, Q),
      method = "2sls"
    ))
  }

  ZtX <- crossprod(Z_w, X_w)             # (l+1) x (p+1)
  ZtQ <- crossprod(Z_w, Q_w)             # (l+1) x Q

  # X_hat'X_hat = (ZtX)' ZtZ_inv ZtX  — small (p+1) x (p+1) matrix
  XhXh     <- crossprod(ZtX, ZtZ_inv %*% ZtX)
  XhXh_inv <- tryCatch(solve(XhXh), error = function(e) NULL)
  if (is.null(XhXh_inv)) {
    return(list(
      beta0  = rep(NA_real_, Q),
      beta1  = if (p == 1) rep(NA_real_, Q) else matrix(NA_real_, p, Q),
      method = "2sls"
    ))
  }

  # beta = (X_hat'X_hat)^{-1} X_hat'Q_w = XhXh_inv %*% ZtX' ZtZ_inv ZtQ
  # All operations on small matrices — no M x M allocation
  beta_all <- XhXh_inv %*% crossprod(ZtX, ZtZ_inv %*% ZtQ)  # (p+1) x Q

  beta0 <- as.vector(beta_all[1, ])
  beta1 <- if (p == 1) {
    as.vector(beta_all[2, ])
  } else {
    beta_all[2:(p + 1), , drop = FALSE]  # p x Q
  }

  list(beta0 = beta0, beta1 = beta1, method = "2sls")
}


# =============================================================================
# Estimator 2: D-IV (pointwise projection then OLS)
# =============================================================================

#' D-IV estimator: the paper's primary estimator (Section 3.3).
#'
#' Four steps:
#'   1. Compute unconstrained 2SLS coefficients.
#'   2. Evaluate fitted curves psi_hat(X_j, u) at each observed X_j.
#'   3. Project each curve to monotonicity via PAVA.
#'   4. Recover coefficients by OLS regression of projected curves on X.
#'
#' @inheritParams estimate_2sls
#' @param return_internals If TRUE, also return unconstrained coefficients,
#'   fitted curves, and projected curves (useful for inference).
#' @return list(beta0, beta1, method) plus optional internals
estimate_div <- function(Q_Yk, X, Z, q_grid, weights = NULL,
                         return_internals = FALSE) {
  M <- nrow(Q_Yk)
  Q <- length(q_grid)

  X_mat <- if (is.matrix(X)) X else matrix(X, ncol = 1)
  p <- ncol(X_mat)

  # Weighted or unweighted means for centering
  if (is.null(weights)) {
    mu_X <- colMeans(X_mat)
    sw <- rep(1, M)
  } else {
    w <- weights / mean(weights)
    sw <- sqrt(w)
    mu_X <- colSums(w * X_mat) / sum(w)
  }

  # Step 1: unconstrained 2SLS
  fit_unc <- estimate_2sls(Q_Yk, X, Z, q_grid, weights = weights)
  if (any(is.na(fit_unc$beta0))) {
    return(list(
      beta0  = rep(NA_real_, Q),
      beta1  = if (p == 1) rep(NA_real_, Q) else matrix(NA_real_, p, Q),
      method = "div"
    ))
  }

  beta0_tilde <- fit_unc$beta0  # Q-vector
  # Ensure beta1_tilde is p x Q
  beta1_tilde <- if (p == 1) {
    matrix(fit_unc$beta1, nrow = 1)
  } else {
    fit_unc$beta1
  }

  # Step 2: evaluate fitted curves at each observed X_j
  # psi_hat(X_j, u) = beta0(u) + (X_j - mu_X)' beta1(u)
  X_centered <- sweep(X_mat, 2, mu_X)  # M x p
  # psi_mat = M x Q: each row is beta0 + X_centered[j,] %*% beta1_tilde
  psi_mat <- matrix(rep(beta0_tilde, each = M), nrow = M) +
             X_centered %*% beta1_tilde  # M x Q

  # Step 3: project each row to monotonicity via PAVA
  Q_hat <- t(apply(psi_mat, 1, pava))  # M x Q

  # Step 4: recover coefficients by weighted OLS (vectorized across all quantiles)
  X_tilde <- cbind(1, X_centered)  # M x (p+1)
  X_w <- X_tilde * sw
  Q_w <- Q_hat * sw
  XtX_inv <- solve(crossprod(X_w))
  beta_all <- XtX_inv %*% crossprod(X_w, Q_w)  # (p+1) x Q

  beta0 <- as.vector(beta_all[1, ])
  beta1 <- if (p == 1) {
    as.vector(beta_all[2, ])
  } else {
    beta_all[2:(p + 1), , drop = FALSE]
  }

  result <- list(beta0 = beta0, beta1 = beta1, method = "div")

  if (return_internals) {
    # Unconstrained coefficients as (p+1) x Q matrix
    result$beta_unc_all <- rbind(beta0_tilde, beta1_tilde)
    result$psi_mat <- psi_mat  # M x Q unconstrained fitted curves
    result$Q_hat   <- Q_hat    # M x Q projected curves
  }

  result
}
