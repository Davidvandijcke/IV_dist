# =============================================================================
# D-IV Inference: Pointwise and Uniform Confidence Bands
# =============================================================================
#
# Implements inference for the D-IV estimator based on Theorems 3-4 of the
# paper. Under Average Strictness, the unconstrained 2SLS and D-IV estimators
# share the same asymptotic Gaussian process, so inference uses the same
# influence function for both.
#
# Public functions:
#   div_pointwise_ci()      - Pointwise confidence intervals
#   div_uniform_cb()        - Uniform confidence bands (multiplier bootstrap)
#   div_bootstrap_process() - Raw bootstrap draws for custom inference
#
# Internal:
#   .inference_setup()      - Shared computation (residuals, SEs, 2SLS operator)
#
# Required packages: none (base R only)
# Depends on: R/utils.R (pava), R/estimators.R (estimate_2sls, estimate_div)
# =============================================================================

# Source dependencies when loaded directly
local({
  this_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) NULL
  )
  if (is.null(this_dir) || this_dir == "") {
    this_dir <- "R"
  }
  source(file.path(this_dir, "utils.R"), local = FALSE)
  source(file.path(this_dir, "estimators.R"), local = FALSE)
})


# =============================================================================
# Internal: shared inference setup
# =============================================================================

#' Compute all shared quantities needed for inference
#'
#' Builds the 2SLS operator, unconstrained residuals, influence matrix,
#' and pointwise standard errors. Used by all three public functions.
#'
#' @param Q_Yk M x Q matrix of group quantile functions
#' @param X    M-vector or M x p matrix of covariates
#' @param Z    M-vector or M x l matrix of instruments
#' @param q_grid Q-vector of quantile levels
#' @param weights Optional M-vector of analytic weights
#' @param cluster Optional M-vector of cluster identifiers for cluster-robust
#'   inference. When provided, the meat of the sandwich sums influence
#'   functions within clusters, and the bootstrap draws one weight per cluster.
#' @return Named list of inference components
.inference_setup <- function(Q_Yk, X, Z, q_grid, weights = NULL, cluster = NULL) {
  M <- nrow(Q_Yk)
  Q <- length(q_grid)

  X_mat <- if (is.matrix(X)) X else matrix(X, ncol = 1)
  Z_mat <- if (is.matrix(Z)) Z else matrix(Z, ncol = 1)
  p <- ncol(X_mat)
  l <- ncol(Z_mat)

  # --- Weights ---
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

  # --- Augmented centered design matrices (unweighted) ---
  X_aug <- cbind(1, sweep(X_mat, 2, mu_X))  # M x (p+1)
  Z_aug <- cbind(1, sweep(Z_mat, 2, mu_Z))  # M x (l+1)

  # --- Weighted versions ---
  X_w <- X_aug * sw  # M x (p+1)
  Z_w <- Z_aug * sw  # M x (l+1)
  Q_w <- Q_Yk * sw   # M x Q

  # --- Sample moment matrices ---
  Sig_ZZ <- crossprod(Z_w) / M      # (l+1) x (l+1)
  Sig_ZX <- crossprod(Z_w, X_w) / M # (l+1) x (p+1)

  Sig_ZZ_inv <- tryCatch(solve(Sig_ZZ), error = function(e) NULL)
  if (is.null(Sig_ZZ_inv)) {
    stop("Sigma_ZZ is singular; instruments are collinear")
  }

  # --- 2SLS operator: S_hat = (Sig_ZX' Sig_ZZ^{-1} Sig_ZX)^{-1} Sig_ZX' Sig_ZZ^{-1} ---
  # Dimension: (p+1) x (l+1)
  A <- crossprod(Sig_ZX, Sig_ZZ_inv %*% Sig_ZX)  # (p+1) x (p+1)
  A_inv <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(A_inv)) {
    stop("2SLS projection matrix is singular; check instrument relevance")
  }
  S_hat <- A_inv %*% crossprod(Sig_ZX, Sig_ZZ_inv)  # (p+1) x (l+1)

  # --- Unconstrained 2SLS coefficients ---
  # beta_unc_all = S_hat %*% (1/M) Z_w' Q_w  =  (p+1) x Q
  beta_unc_all <- S_hat %*% (crossprod(Z_w, Q_w) / M)

  # --- Residuals ---
  # Xi[j,q] = Q_Yj(u_q) - X_aug_j' beta_unc(u_q)  (unweighted)
  Xi <- Q_Yk - X_aug %*% beta_unc_all  # M x Q
  # Xi_w[j,q] = sqrt(w_j) * xi_j(q)  (weighted residuals for the WLS influence function)
  # In the transformed model y* = X* beta + e*, the residual is e* = sqrt(w) * xi.
  # The HC0 sandwich uses (e*)^2, so the score is Phi_j * e*_j = Phi_j * sw_j * xi_j.
  Xi_w <- Xi * sw  # M x Q

  # --- Influence matrix: Phi[j,k] = Z_w_j' S_hat[k,]' ---
  # Phi = Z_w %*% t(S_hat), dimension M x (p+1)
  # Note: Phi already contains sqrt(w) through Z_w.
  Phi <- Z_w %*% t(S_hat)

  # --- Clustering ---
  # With clustering, the effective "number of observations" for the CLT is the
  # number of clusters C, not M. The meat sums influence functions within
  # clusters: Meat_cl(u) = (1/M) sum_c (sum_{j in c} Phi_j xi_j(u))^2
  # and SE = sqrt(Omega_cl / M), scaled by sqrt(C/(C-1)) for small-sample
  # correction (analogous to HC1 for clusters).
  if (!is.null(cluster)) {
    cluster <- as.factor(cluster)
    cluster_ids <- levels(cluster)
    C <- length(cluster_ids)
    cluster_idx <- as.integer(cluster)  # integer mapping j -> cluster index

    # Cluster-robust SEs: sum scores within clusters, then compute variance.
    # Score = Phi[j,k] * Xi_w[j,q]. Small-sample correction: C/(C-1).
    se_mat <- matrix(NA_real_, p + 1, Q)

    for (k in seq_len(p + 1)) {
      # Score: Phi[j,k] * Xi_w[j,q] = Phi[j,k] * sqrt(w_j) * xi_j(q)
      phi_k_Xi_w <- Phi[, k] * Xi_w  # M x Q
      # Sum within clusters: cluster_scores[c, q]
      cluster_scores <- rowsum(phi_k_Xi_w, cluster_idx)  # C x Q
      # Variance: (1/M^2) * sum_c cluster_scores[c,q]^2 * C/(C-1)
      se_mat[k, ] <- sqrt(colSums(cluster_scores^2) * C / ((C - 1) * M^2))
    }
  } else {
    cluster_ids <- NULL
    cluster_idx <- NULL
    C <- M

    # --- Pointwise standard errors (HC0, vectorized over Q) ---
    # Score = Phi[j,k] * Xi_w[j,q] where Xi_w = sqrt(w) * xi (weighted residual)
    # se[k,q] = sqrt( mean(Phi[,k]^2 * Xi_w[,q]^2) / M )
    se_mat <- matrix(NA_real_, p + 1, Q)
    for (k in seq_len(p + 1)) {
      phi_k_sq <- Phi[, k]^2
      se_mat[k, ] <- sqrt(colMeans(phi_k_sq * Xi_w^2) / M)
    }
  }

  list(
    M = M, Q = Q, p = p, l = l, C = C,
    X_aug = X_aug, Z_aug = Z_aug,
    X_w = X_w, Z_w = Z_w, Q_w = Q_w,
    sw = sw, mu_X = mu_X, mu_Z = mu_Z,
    S_hat = S_hat,
    beta_unc_all = beta_unc_all,
    Xi = Xi, Xi_w = Xi_w,
    Phi = Phi,
    se_mat = se_mat,
    q_grid = q_grid,
    cluster_idx = cluster_idx,
    cluster_ids = cluster_ids
  )
}


# =============================================================================
# Function 1: Pointwise confidence intervals
# =============================================================================

#' Pointwise confidence intervals for D-IV coefficient functions
#'
#' Computes sandwich-based pointwise CIs for both the unconstrained 2SLS
#' and D-IV coefficient functions. Under Average Strictness (Assumption 3),
#' both estimators share the same asymptotic variance, so the SEs are
#' identical for both.
#'
#' @param Q_Yk  M x Q matrix of group quantile functions
#' @param X     M-vector or M x p matrix of covariates
#' @param Z     M-vector or M x l matrix of instruments
#' @param q_grid Q-vector of quantile levels
#' @param weights Optional M-vector of analytic weights
#' @param cluster Optional M-vector of cluster identifiers for cluster-robust SEs
#' @param alpha  Significance level (default 0.05 for 95% CIs)
#' @return Named list with components:
#'   \item{q_grid}{Quantile grid}
#'   \item{alpha}{Significance level}
#'   \item{se}{(p+1) x Q matrix of pointwise standard errors}
#'   \item{ci_unc}{List with est, lo, hi: (p+1) x Q matrices for unconstrained 2SLS}
#'   \item{ci_div}{List with est, lo, hi: (p+1) x Q matrices for D-IV}
div_pointwise_ci <- function(Q_Yk, X, Z, q_grid, weights = NULL,
                             cluster = NULL, alpha = 0.05) {
  setup <- .inference_setup(Q_Yk, X, Z, q_grid, weights, cluster = cluster)
  z_crit <- qnorm(1 - alpha / 2)

  # Point estimates: unconstrained (already in setup)
  beta_unc <- setup$beta_unc_all  # (p+1) x Q

  # Point estimates: D-IV
  fit_div <- estimate_div(Q_Yk, X, Z, q_grid, weights = weights)
  p <- setup$p
  beta_div <- rbind(
    fit_div$beta0,
    if (p == 1) matrix(fit_div$beta1, nrow = 1) else fit_div$beta1
  )  # (p+1) x Q

  # CIs use the same SEs for both (same asymptotic distribution)
  se <- setup$se_mat
  half_width <- z_crit * se

  list(
    q_grid = q_grid,
    alpha  = alpha,
    se     = se,
    ci_unc = list(
      est = beta_unc,
      lo  = beta_unc - half_width,
      hi  = beta_unc + half_width
    ),
    ci_div = list(
      est = beta_div,
      lo  = beta_div - half_width,
      hi  = beta_div + half_width
    )
  )
}


# =============================================================================
# Function 2: Uniform confidence bands (multiplier bootstrap)
# =============================================================================

#' Uniform confidence bands for D-IV coefficient functions
#'
#' Uses a multiplier bootstrap to construct uniform confidence bands over
#' the quantile grid. Two bootstrap variants are available:
#'
#' 1. **Unprojected** (projected=FALSE): Bootstraps the unconstrained
#'    2SLS process. Valid for both unconstrained and D-IV coefficients
#'    (same limiting process under Average Strictness).
#'
#' 2. **Projected** (projected=TRUE): Additionally runs PAVA inside each
#'    bootstrap draw. Can give tighter bands in finite samples, especially
#'    with weak instruments. Only valid for D-IV coefficients.
#'
#' @inheritParams div_pointwise_ci
#' @param B         Number of bootstrap replications (default 1000)
#' @param projected If TRUE, also compute projected bootstrap bands for D-IV
#' @param multiplier Bootstrap weight distribution: "gaussian" (default) or "rademacher"
#' @param seed      Optional random seed for reproducibility
#' @return Named list with components:
#'   \item{q_grid, alpha, B, se}{As in div_pointwise_ci}
#'   \item{c_alpha}{(p+1)-vector of unprojected critical values}
#'   \item{ucb_unc}{List with est, lo, hi for unconstrained uniform bands}
#'   \item{ucb_div}{List with est, lo, hi for D-IV uniform bands (unprojected)}
#'   \item{c_alpha_proj}{(p+1)-vector of projected critical values (if projected=TRUE)}
#'   \item{ucb_div_proj}{List with est, lo, hi for D-IV projected uniform bands (if projected=TRUE)}
#'   \item{se_proj}{(p+1) x Q projected bootstrap SEs (if projected=TRUE)}
#'   \item{ci_div_proj_pw}{Normal-based pointwise CIs using projected bootstrap SEs (if projected=TRUE)}
#'   \item{ci_div_proj_pctl}{Percentile-based pointwise CIs from projected bootstrap (if projected=TRUE)}
div_uniform_cb <- function(Q_Yk, X, Z, q_grid, weights = NULL,
                           cluster = NULL,
                           alpha = 0.05, B = 1000L, projected = FALSE,
                           multiplier = c("gaussian", "rademacher"),
                           seed = NULL) {
  multiplier <- match.arg(multiplier)
  setup <- .inference_setup(Q_Yk, X, Z, q_grid, weights, cluster = cluster)

  M   <- setup$M
  Q   <- setup$Q
  p   <- setup$p
  Phi <- setup$Phi          # M x (p+1)
  Xi  <- setup$Xi           # M x Q
  se  <- setup$se_mat       # (p+1) x Q

  # Point estimates
  beta_unc <- setup$beta_unc_all  # (p+1) x Q

  fit_div <- estimate_div(Q_Yk, X, Z, q_grid, weights = weights)
  beta_div <- rbind(
    fit_div$beta0,
    if (p == 1) matrix(fit_div$beta1, nrow = 1) else fit_div$beta1
  )

  # Guard against zero SEs (would cause division by zero in studentization)
  se_safe <- pmax(se, .Machine$double.eps)

  # --- Unprojected bootstrap ---
  if (!is.null(seed)) set.seed(seed)

  T_mat <- matrix(NA_real_, B, p + 1)  # sup stats per coefficient

  # Projected bootstrap setup (if requested)
  if (projected) {
    psi_hat    <- setup$X_aug %*% beta_unc  # M x Q: unconstrained fitted curves
    XtX_w_inv  <- solve(crossprod(setup$X_w))
    T_proj_mat <- matrix(NA_real_, B, p + 1)
    # Accumulate per-quantile draws for pointwise inference
    delta_div_sq <- matrix(0, p + 1, Q)        # running sum of squares
    delta_div_draws <- array(NA_real_, dim = c(B, p + 1, Q))  # for percentile CIs
  }

  # Cluster bootstrap: draw one weight per cluster, expand to observation level
  clustered <- !is.null(setup$cluster_idx)
  C <- setup$C

  for (b in seq_len(B)) {
    # Draw multiplier weights (one per cluster, or one per obs if unclustered)
    if (clustered) {
      omega_cl <- if (multiplier == "gaussian") {
        rnorm(C)
      } else {
        sample(c(-1, 1), C, replace = TRUE)
      }
      omega <- omega_cl[setup$cluster_idx]  # expand to M-vector
    } else {
      omega <- if (multiplier == "gaussian") {
        rnorm(M)
      } else {
        sample(c(-1, 1), M, replace = TRUE)
      }
    }

    # Bootstrap draw: G_b = t(Phi) %*% (omega * Xi_w) / sqrt(M)
    # Uses Xi_w (weighted residuals) for correct HC0 sandwich in WLS
    # Dimension: (p+1) x Q
    G_b <- crossprod(Phi, omega * setup$Xi_w) / sqrt(M)

    # Studentized sup statistic for each coefficient
    for (k in seq_len(p + 1)) {
      T_mat[b, k] <- max(abs(G_b[k, ]) / (sqrt(M) * se_safe[k, ]))
    }

    # --- Projected bootstrap ---
    if (projected) {
      # Perturbed unconstrained curves: psi_star = psi_hat + X_aug %*% G_b / sqrt(M)
      # G_b is already the (p+1) x Q bootstrap perturbation
      psi_star <- psi_hat + setup$X_aug %*% G_b / sqrt(M)  # M x Q

      # Project each curve to monotonicity
      Q_star <- t(apply(psi_star, 1, pava))  # M x Q

      # OLS to recover D-IV bootstrap coefficients
      beta_star <- XtX_w_inv %*% crossprod(setup$X_w, Q_star * setup$sw)  # (p+1) x Q

      # Centered bootstrap statistic
      delta_div <- sqrt(M) * (beta_star - beta_div)  # (p+1) x Q

      for (k in seq_len(p + 1)) {
        T_proj_mat[b, k] <- max(abs(delta_div[k, ]) / (sqrt(M) * se_safe[k, ]))
      }

      # Accumulate for pointwise inference
      delta_div_sq <- delta_div_sq + delta_div^2
      delta_div_draws[b, , ] <- delta_div
    }
  }

  # --- Critical values ---
  c_alpha <- apply(T_mat, 2, quantile, probs = 1 - alpha)  # (p+1)-vector

  # --- Uniform bands ---
  result <- list(
    q_grid  = q_grid,
    alpha   = alpha,
    B       = B,
    se      = se,
    c_alpha = c_alpha,
    ucb_unc = list(
      est = beta_unc,
      lo  = beta_unc - outer(c_alpha, rep(1, Q)) * se,
      hi  = beta_unc + outer(c_alpha, rep(1, Q)) * se
    ),
    ucb_div = list(
      est = beta_div,
      lo  = beta_div - outer(c_alpha, rep(1, Q)) * se,
      hi  = beta_div + outer(c_alpha, rep(1, Q)) * se
    )
  )

  if (projected) {
    c_alpha_proj <- apply(T_proj_mat, 2, quantile, probs = 1 - alpha)
    result$c_alpha_proj <- c_alpha_proj
    result$ucb_div_proj <- list(
      est = beta_div,
      lo  = beta_div - outer(c_alpha_proj, rep(1, Q)) * se,
      hi  = beta_div + outer(c_alpha_proj, rep(1, Q)) * se
    )

    # --- Projected bootstrap pointwise CIs ---
    # Bootstrap SE: sqrt( (1/B) sum_b delta_div_b^2 / M )
    se_proj <- sqrt(delta_div_sq / (B * M))  # (p+1) x Q

    # Normal-based pointwise CI using bootstrap SE
    z_crit <- qnorm(1 - alpha / 2)
    result$se_proj <- se_proj
    result$ci_div_proj_pw <- list(
      est = beta_div,
      lo  = beta_div - z_crit * se_proj,
      hi  = beta_div + z_crit * se_proj
    )

    # Percentile-based pointwise CI (captures asymmetry from PAVA)
    lo_pctl <- matrix(NA_real_, p + 1, Q)
    hi_pctl <- matrix(NA_real_, p + 1, Q)
    for (k in seq_len(p + 1)) {
      for (q in seq_len(Q)) {
        qs <- quantile(delta_div_draws[, k, q],
                       probs = c(alpha / 2, 1 - alpha / 2))
        hi_pctl[k, q] <- beta_div[k, q] - qs[1] / sqrt(M)
        lo_pctl[k, q] <- beta_div[k, q] - qs[2] / sqrt(M)
      }
    }
    result$ci_div_proj_pctl <- list(
      est = beta_div,
      lo  = lo_pctl,
      hi  = hi_pctl
    )
  }

  result
}


# =============================================================================
# Function 3: Raw bootstrap process
# =============================================================================

#' Bootstrap process for D-IV inference
#'
#' Returns the full bootstrap draws, useful for custom tests or
#' visualization of the bootstrap distribution.
#'
#' @inheritParams div_uniform_cb
#' @return Named list with components:
#'   \item{q_grid}{Quantile grid}
#'   \item{B}{Number of bootstrap draws}
#'   \item{se}{(p+1) x Q matrix of pointwise SEs}
#'   \item{beta_unc}{(p+1) x Q unconstrained point estimates}
#'   \item{beta_div}{(p+1) x Q D-IV point estimates}
#'   \item{G_star}{List of B matrices, each (p+1) x Q: unprojected bootstrap draws}
#'   \item{G_div_star}{List of B matrices (if projected=TRUE): projected bootstrap draws}
div_bootstrap_process <- function(Q_Yk, X, Z, q_grid, weights = NULL,
                                  cluster = NULL,
                                  B = 1000L, projected = FALSE,
                                  multiplier = c("gaussian", "rademacher"),
                                  seed = NULL) {
  multiplier <- match.arg(multiplier)
  setup <- .inference_setup(Q_Yk, X, Z, q_grid, weights, cluster = cluster)

  M   <- setup$M
  Q   <- setup$Q
  p   <- setup$p
  Phi <- setup$Phi
  Xi  <- setup$Xi

  beta_unc <- setup$beta_unc_all

  fit_div <- estimate_div(Q_Yk, X, Z, q_grid, weights = weights)
  beta_div <- rbind(
    fit_div$beta0,
    if (p == 1) matrix(fit_div$beta1, nrow = 1) else fit_div$beta1
  )

  if (!is.null(seed)) set.seed(seed)

  G_star <- vector("list", B)

  if (projected) {
    psi_hat    <- setup$X_aug %*% beta_unc
    XtX_w_inv  <- solve(crossprod(setup$X_w))
    G_div_star <- vector("list", B)
  }

  clustered <- !is.null(setup$cluster_idx)
  C <- setup$C

  for (b in seq_len(B)) {
    if (clustered) {
      omega_cl <- if (multiplier == "gaussian") rnorm(C) else sample(c(-1, 1), C, replace = TRUE)
      omega <- omega_cl[setup$cluster_idx]
    } else {
      omega <- if (multiplier == "gaussian") rnorm(M) else sample(c(-1, 1), M, replace = TRUE)
    }

    G_b <- crossprod(Phi, omega * setup$Xi_w) / sqrt(M)  # (p+1) x Q
    G_star[[b]] <- G_b

    if (projected) {
      psi_star   <- psi_hat + setup$X_aug %*% G_b / sqrt(M)
      Q_star     <- t(apply(psi_star, 1, pava))
      beta_star  <- XtX_w_inv %*% crossprod(setup$X_w, Q_star * setup$sw)
      G_div_star[[b]] <- sqrt(M) * (beta_star - beta_div)
    }
  }

  result <- list(
    q_grid   = q_grid,
    B        = B,
    se       = setup$se_mat,
    beta_unc = beta_unc,
    beta_div = beta_div,
    G_star   = G_star
  )

  if (projected) {
    result$G_div_star <- G_div_star
  }

  result
}
