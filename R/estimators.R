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
# Required packages: quadprog (for multivariate endpoint projection)
# =============================================================================

# Source utility functions (pava, compute_2sls_coefs).
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
#' @return list(beta0, beta1, method) where beta1 is Q-vector (scalar X) or p x Q matrix
estimate_2sls <- function(Q_Yk, X, Z, q_grid) {
  M  <- nrow(Q_Yk)
  Q  <- length(q_grid)

  # Coerce to matrix form
  X_mat <- if (is.matrix(X)) X else matrix(X, ncol = 1)
  Z_mat <- if (is.matrix(Z)) Z else matrix(Z, ncol = 1)
  p <- ncol(X_mat)
  l <- ncol(Z_mat)

  # Centered design and instrument matrices (with intercept)
  mu_X <- colMeans(X_mat)
  mu_Z <- colMeans(Z_mat)
  X_design <- cbind(1, sweep(X_mat, 2, mu_X))   # M x (p+1)
  Z_design <- cbind(1, sweep(Z_mat, 2, mu_Z))   # M x (l+1)

  # Instrument projection: P_Z = Z (Z'Z)^{-1} Z'
  ZtZ <- crossprod(Z_design)
  ZtZ_inv <- tryCatch(solve(ZtZ), error = function(e) NULL)
  if (is.null(ZtZ_inv)) {
    return(list(
      beta0  = rep(NA_real_, Q),
      beta1  = if (p == 1) rep(NA_real_, Q) else matrix(NA_real_, p, Q),
      method = "2sls"
    ))
  }

  P_Z   <- Z_design %*% ZtZ_inv %*% t(Z_design)
  X_hat <- P_Z %*% X_design  # M x (p+1)

  # Check rank of X_hat'X_hat
  XhXh     <- crossprod(X_hat)
  XhXh_inv <- tryCatch(solve(XhXh), error = function(e) NULL)
  if (is.null(XhXh_inv)) {
    return(list(
      beta0  = rep(NA_real_, Q),
      beta1  = if (p == 1) rep(NA_real_, Q) else matrix(NA_real_, p, Q),
      method = "2sls"
    ))
  }

  # Compute all coefficients at once: (p+1) x Q matrix
  # beta_all[, iq] = (X_hat'X_hat)^{-1} X_hat' Q_Yk[, iq]
  beta_all <- XhXh_inv %*% crossprod(X_hat, Q_Yk)  # (p+1) x Q

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
#' @return list(beta0, beta1, method)
estimate_div <- function(Q_Yk, X, Z, q_grid) {
  M <- nrow(Q_Yk)
  Q <- length(q_grid)

  X_mat <- if (is.matrix(X)) X else matrix(X, ncol = 1)
  p <- ncol(X_mat)
  mu_X <- colMeans(X_mat)

  # Step 1: unconstrained 2SLS
  fit_unc <- estimate_2sls(Q_Yk, X, Z, q_grid)
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

  # Step 4: recover coefficients by OLS (vectorized across all quantiles)
  X_tilde <- cbind(1, X_centered)  # M x (p+1)
  XtX_inv <- solve(crossprod(X_tilde))
  beta_all <- XtX_inv %*% crossprod(X_tilde, Q_hat)  # (p+1) x Q

  beta0 <- as.vector(beta_all[1, ])
  beta1 <- if (p == 1) {
    as.vector(beta_all[2, ])
  } else {
    beta_all[2:(p + 1), , drop = FALSE]
  }

  list(beta0 = beta0, beta1 = beta1, method = "div")
}


# =============================================================================
# Estimator 3: D-IV Endpoint Projection
# =============================================================================

#' D-IV endpoint projection estimator.
#'
#' For scalar X: projects the two endpoint quantile functions gamma_{-} and
#' gamma_{+} via PAVA, then recovers coefficients analytically.
#'
#' For multivariate X (p > 1): solves a QP that enforces monotonicity at all
#' 2^p vertices of the support hypercube.
#'
#' @inheritParams estimate_2sls
#' @return list(beta0, beta1, method)
estimate_div_endpoint <- function(Q_Yk, X, Z, q_grid) {
  X_mat <- if (is.matrix(X)) X else matrix(X, ncol = 1)
  p <- ncol(X_mat)

  if (p == 1) {
    return(estimate_div_endpoint_scalar(Q_Yk, X_mat[, 1], Z, q_grid))
  } else {
    return(estimate_div_endpoint_multi(Q_Yk, X_mat, Z, q_grid))
  }
}


# -----------------------------------------------------------------------------
# Scalar endpoint projection
# -----------------------------------------------------------------------------

#' Endpoint projection for scalar X.
#'
#' 1. Compute unconstrained 2SLS: beta_tilde_0(u), beta_tilde_1(u)
#' 2. Form endpoint curves at min(X) and max(X)
#' 3. Project each via PAVA
#' 4. Recover coefficients from the two projected endpoints
#'
#' @param Q_Yk M x Q matrix
#' @param x    M-vector (scalar covariate)
#' @param Z    M-vector or M x l matrix of instruments
#' @param q_grid Q-vector
#' @return list(beta0, beta1, method)
estimate_div_endpoint_scalar <- function(Q_Yk, x, Z, q_grid) {
  Q <- length(q_grid)

  # Unconstrained 2SLS
  fit_unc <- estimate_2sls(Q_Yk, x, Z, q_grid)
  if (any(is.na(fit_unc$beta0))) {
    return(list(beta0 = rep(NA_real_, Q), beta1 = rep(NA_real_, Q),
                method = "div_endpoint"))
  }

  beta0_tilde <- fit_unc$beta0
  beta1_tilde <- fit_unc$beta1  # Q-vector

  mu_x   <- mean(x)
  a_minus <- min(x) - mu_x
  a_plus  <- max(x) - mu_x

  # Endpoint curves
  gamma_minus <- beta0_tilde + a_minus * beta1_tilde
  gamma_plus  <- beta0_tilde + a_plus  * beta1_tilde

  # Project to monotonicity
  gamma_minus_proj <- pava(gamma_minus)
  gamma_plus_proj  <- pava(gamma_plus)

  # Recover coefficients
  denom <- a_plus - a_minus
  beta1 <- (gamma_plus_proj - gamma_minus_proj) / denom
  beta0 <- (a_plus * gamma_minus_proj - a_minus * gamma_plus_proj) / denom

  list(beta0 = beta0, beta1 = beta1, method = "div_endpoint")
}


# -----------------------------------------------------------------------------
# Multivariate endpoint projection via QP
# -----------------------------------------------------------------------------

#' Endpoint projection for multivariate X via quadratic programming.
#'
#' Stacks all coefficients: theta = c(beta_0, beta_1_1, ..., beta_1_p), length
#' (p+1)*Q. Generates 2^p vertices of the support hypercube. For each vertex v
#' and consecutive quantile pair, enforces gamma_v(u_{q+1}) >= gamma_v(u_q).
#' Solves the QP: min ||theta - theta_tilde||^2 subject to A theta >= 0.
#'
#' @param Q_Yk M x Q matrix
#' @param X    M x p matrix
#' @param Z    M-vector or M x l matrix
#' @param q_grid Q-vector
#' @return list(beta0, beta1 as p x Q matrix, method)
estimate_div_endpoint_multi <- function(Q_Yk, X, Z, q_grid) {
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("Package 'quadprog' is required for multivariate endpoint projection.")
  }

  M <- nrow(Q_Yk)
  Q <- length(q_grid)
  p <- ncol(X)
  mu_X <- colMeans(X)

  # Step 1: unconstrained 2SLS coefficients
  fit_unc <- estimate_2sls(Q_Yk, X, Z, q_grid)
  if (any(is.na(fit_unc$beta0))) {
    return(list(beta0 = rep(NA_real_, Q),
                beta1 = matrix(NA_real_, p, Q),
                method = "div_endpoint"))
  }

  # Stack theta_tilde = c(beta0, beta1_1, ..., beta1_p), length (p+1)*Q
  beta1_mat <- if (p == 1) matrix(fit_unc$beta1, nrow = 1) else fit_unc$beta1
  theta_tilde <- c(fit_unc$beta0)
  for (j in seq_len(p)) {
    theta_tilde <- c(theta_tilde, beta1_mat[j, ])
  }
  n_vars <- (p + 1) * Q

  # Step 2: generate vertices of support hypercube (centered)
  x_min <- apply(X, 2, min)
  x_max <- apply(X, 2, max)
  vertices_raw <- as.matrix(
    expand.grid(lapply(seq_len(p), function(j) c(x_min[j], x_max[j])))
  )
  # Center vertices: v - mu_X
  vertices <- sweep(vertices_raw, 2, mu_X)  # 2^p x p
  n_vertices <- nrow(vertices)

  # Step 3: build differencing matrix D (first-differences, Q-1 x Q)
  D <- matrix(0, Q - 1, Q)
  for (i in seq_len(Q - 1)) {
    D[i, i]     <- -1
    D[i, i + 1] <-  1
  }

  # Step 4: build constraint matrix A_ineq (n_constraints x n_vars)
  # For vertex v, constraint is: D %*% gamma_v >= 0
  # gamma_v(u) = beta0(u) + sum_j (v_j - mu_j) * beta1_j(u)
  # In stacked form: gamma_v = [I, v1*I, v2*I, ..., vp*I] %*% theta
  # So A_ineq block for vertex v = D %*% [I, v1*I, ..., vp*I]
  n_constraints <- n_vertices * (Q - 1)
  A_ineq <- matrix(0, n_constraints, n_vars)

  for (v in seq_len(n_vertices)) {
    row_idx <- ((v - 1) * (Q - 1) + 1):(v * (Q - 1))

    # beta0 block: columns 1:Q
    A_ineq[row_idx, 1:Q] <- D

    # beta1_j blocks: columns (j*Q+1):((j+1)*Q)
    for (j in seq_len(p)) {
      col_idx <- (j * Q + 1):((j + 1) * Q)
      A_ineq[row_idx, col_idx] <- D * vertices[v, j]
    }
  }

  # Step 5: solve QP
  # minimize 0.5 * ||theta - theta_tilde||^2
  # quadprog solves: min 0.5 * x'Dx - d'x  s.t. A'x >= b
  Dmat <- diag(n_vars) + 1e-8 * diag(n_vars)
  dvec <- theta_tilde
  Amat <- t(A_ineq)
  bvec <- rep(0, n_constraints)

  result <- tryCatch(
    quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 0),
    error = function(e) {
      # Fallback: return unprojected
      list(solution = theta_tilde)
    }
  )

  theta_proj <- result$solution

  # Extract coefficients
  beta0 <- theta_proj[1:Q]
  beta1 <- matrix(NA_real_, p, Q)
  for (j in seq_len(p)) {
    beta1[j, ] <- theta_proj[(j * Q + 1):((j + 1) * Q)]
  }

  list(beta0 = beta0, beta1 = beta1, method = "div_endpoint")
}
