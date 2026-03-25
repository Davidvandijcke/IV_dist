# =============================================================================
# D-IV Estimator: First-Stage Functions
# =============================================================================
#
# Required packages: stats (base), quantreg (for compute_qr_fitted only)
#
# Computes group-level quantile functions from individual-level data.
# These serve as the "outcome" in the second-stage IV regression.
#
# Functions:
#   compute_sample_quantiles() - Empirical quantile functions (D-IV method)
#   compute_qr_fitted()        - Within-group quantile regression (CLP/MP)
#
# =============================================================================


#' Compute sample quantile functions per group
#'
#' For each group, computes the empirical quantile function evaluated at a
#' common grid of quantile probabilities. This is the first-stage estimator
#' used in D-IV.
#'
#' @param y_list List of length M, each element a numeric vector of individual
#'   outcomes for that group.
#' @param q_grid Numeric vector of quantile probabilities in (0, 1).
#' @return M x Q matrix where row j is the empirical quantile function for
#'   group j.
compute_sample_quantiles <- function(y_list, q_grid) {
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


#' Within-group quantile regression fitted values
#'
#' For each group j and quantile level u, runs a linear quantile regression
#' of the outcome on within-group covariates:
#'
#'   Q_{Y|X1}(u | x1) = alpha_j(u) + gamma_j(u) * x1
#'
#' Returns both the group-level intercepts (used by CLP-type methods) and
#' the full matrix of fitted values (used by MP-type methods).
#'
#' Falls back to sample quantiles for any group where QR fails (common with
#' small within-group sample sizes or collinear data).
#'
#' @param y_list List of length M, each element a numeric vector of individual
#'   outcomes for group j.
#' @param x1_list List of length M, each element a numeric vector of
#'   within-group covariates for group j. Must have same length as the
#'   corresponding element of y_list.
#' @param q_grid Numeric vector of quantile probabilities in (0, 1).
#' @return Named list:
#'   \item{intercepts}{M x Q matrix of group-level intercepts alpha_j(u).
#'     For groups where QR failed, this is the sample quantile.}
#'   \item{fitted_values}{List of length M. Each element is an N_j x Q matrix
#'     of fitted quantile values at observed within-group covariate values.
#'     For groups where QR failed, each row is the sample quantile (constant
#'     across individuals).}
#'   \item{fallback}{Logical vector of length M indicating which groups fell
#'     back to sample quantiles.}
compute_qr_fitted <- function(y_list, x1_list, q_grid) {
  M <- length(y_list)
  Q <- length(q_grid)

  intercepts <- matrix(NA_real_, nrow = M, ncol = Q)
  fitted_values <- vector("list", M)
  fallback <- logical(M)

  for (j in seq_len(M)) {
    y_j <- y_list[[j]]
    x1_j <- x1_list[[j]]
    N_j <- length(y_j)

    # Fallback: sample quantiles (constant across individuals)
    sq_j <- quantile(y_j, probs = q_grid, type = 7, names = FALSE)

    if (N_j < 3 || length(unique(x1_j)) < 2) {
      # Too few observations or no covariate variation: use sample quantiles
      intercepts[j, ] <- sq_j
      fitted_values[[j]] <- matrix(rep(sq_j, each = N_j), nrow = N_j, ncol = Q)
      fallback[j] <- TRUE
      next
    }

    # Build data frame for rq
    df_j <- data.frame(y = y_j, x1 = x1_j)

    # Run QR across all quantiles at once
    fit <- tryCatch({
      quantreg::rq(y ~ x1, tau = q_grid, data = df_j)
    }, error = function(e) NULL,
       warning = function(w) {
         # Suppress common QR warnings but still try to use the result
         tryCatch(
           suppressWarnings(quantreg::rq(y ~ x1, tau = q_grid, data = df_j)),
           error = function(e) NULL
         )
       })

    if (is.null(fit)) {
      # QR failed entirely: use sample quantiles
      intercepts[j, ] <- sq_j
      fitted_values[[j]] <- matrix(rep(sq_j, each = N_j), nrow = N_j, ncol = Q)
      fallback[j] <- TRUE
      next
    }

    # Extract coefficients: 2 x Q matrix (intercept, slope)
    coefs <- coef(fit)
    if (!is.matrix(coefs)) {
      # Single quantile case: coef returns a named vector
      coefs <- matrix(coefs, ncol = 1)
    }

    intercepts[j, ] <- coefs[1, ]

    # Compute fitted values: N_j x Q
    X_design <- cbind(1, x1_j)  # N_j x 2
    fv_j <- X_design %*% coefs  # N_j x Q
    fitted_values[[j]] <- fv_j
    fallback[j] <- FALSE
  }

  list(intercepts = intercepts,
       fitted_values = fitted_values,
       fallback = fallback)
}
