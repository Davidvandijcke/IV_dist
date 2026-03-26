# =============================================================================
# D-IV Estimator: Core Utility Functions
# =============================================================================
#
# Required packages: stats (base)
#
# Functions:
#   pava()                    - Pool-Adjacent-Violators (isotonic regression)
#   compute_group_quantiles() - Empirical quantile functions per group
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
