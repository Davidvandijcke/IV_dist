# =============================================================================
# Melly & Pons (2025) Replication + D-IV Comparison
# =============================================================================
#
# Replicates the empirical application from Section 6 of Melly & Pons (2025):
# "The Effect of the Food Stamps Program on Birth Weight"
#
# Groups = county × trimester × year × race
# Outcome = birth weight quantile functions (19 points, 0.05 to 0.95)
# Treatment = fsp (binary: food stamp program in place 3+ months before birth)
# Controls = REIS transfers/income + 1960 county characteristics × time
# Fixed effects = county (absorbed by demeaning) + state×year + time
# Clustering = county level
#
# Usage:
#   cd IV_dist && Rscript empirical/mp_replication.R
#
# =============================================================================

# --- Resolve project root ---
PROJECT_ROOT <- tryCatch({
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}, error = function(e) getwd())

# --- Dependencies ---
library(haven)
library(ggplot2)

source(file.path(PROJECT_ROOT, "R", "utils.R"))
source(file.path(PROJECT_ROOT, "R", "estimators.R"))
source(file.path(PROJECT_ROOT, "R", "inference.R"))

# --- Configuration ---
DATA_PATH <- file.path(PROJECT_ROOT, "data", "out", "mp_analysis_data.dta")
FIG_DIR <- file.path(PROJECT_ROOT, "empirical", "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

Q_GRID <- seq(0.05, 0.95, by = 0.05)
Q_LABELS <- paste0("q_", as.integer(Q_GRID * 100))

RACE_CODES <- c(black = 2, white = 1)
RACE_LABELS <- c(black = "Black Mothers", white = "White Mothers")


# =============================================================================
# Load data
# =============================================================================

cat("Loading MP analysis data...\n")
d <- read_dta(DATA_PATH)
cat(sprintf("  %d groups, %d columns\n", nrow(d), ncol(d)))
cat(sprintf("  Black groups: %d\n", sum(d$mrace == 2)))
cat(sprintf("  White groups: %d\n", sum(d$mrace == 1)))


# =============================================================================
# Helper: within-county demeaning (absorb county FEs)
# =============================================================================

#' Demean columns within groups (fixed effects transformation)
#' @param mat Matrix or vector to demean
#' @param group_id Integer vector of group identifiers
#' @return Demeaned matrix/vector
demean_within <- function(mat, group_id) {
  mat <- as.matrix(mat)
  groups <- as.integer(as.factor(group_id))
  n_groups <- max(groups)
  group_sums <- rowsum(mat, groups)
  group_counts <- tabulate(groups, nbins = n_groups)
  group_means <- group_sums / group_counts
  mat - group_means[groups, , drop = FALSE]
}

#' Absorb multiple fixed effects via iterative demeaning (Gauss-Seidel)
#' @param mat Matrix to demean
#' @param ... Vectors of group identifiers (one per FE dimension)
#' @param tol Convergence tolerance
#' @param maxiter Maximum iterations
#' @return Demeaned matrix with all FEs absorbed
demean_multi_fe <- function(mat, ..., tol = 1e-8, maxiter = 100) {
  fe_groups <- list(...)
  mat <- as.matrix(mat)
  for (iter in seq_len(maxiter)) {
    mat_old <- mat
    for (g in fe_groups) {
      mat <- demean_within(mat, g)
    }
    # Check convergence
    delta <- max(abs(mat - mat_old))
    if (delta < tol) {
      if (iter > 1) cat(sprintf("    FE absorption converged in %d iterations\n", iter))
      break
    }
  }
  mat
}


# =============================================================================
# Run D-IV for each race
# =============================================================================

cat("\n=== Running D-IV estimation ===\n")

run_mp_div <- function(data, race_name) {
  race_code <- RACE_CODES[race_name]
  sub <- data[data$mrace == race_code, ]
  M <- nrow(sub)
  cat(sprintf("\n  %s: %d groups\n", RACE_LABELS[race_name], M))

  # --- Build Q_Yk matrix (M x 19) ---
  Q_Yk <- as.matrix(sub[, Q_LABELS])

  # Check for NAs
  na_rows <- rowSums(is.na(Q_Yk)) > 0
  if (any(na_rows)) {
    cat(sprintf("    Dropping %d groups with NA quantiles\n", sum(na_rows)))
    sub <- sub[!na_rows, ]
    Q_Yk <- Q_Yk[!na_rows, ]
    M <- nrow(sub)
  }

  # --- Build treatment and controls ---
  fsp <- sub$fsp

  # Group-level controls (REIS transfers/income + 1960 chars × time)
  control_cols <- c()
  reis_cols <- intersect(names(sub),
    c("tranpcret", "tranpcmed", "tranpcpub", "tranpcvet", "tranpcoth",
      "ripc", "afdcpcret", "fspcret", "ssipcrret", "medpcret"))
  control_cols <- c(control_cols, reis_cols)
  time_cols <- intersect(names(sub),
    c("black60_time", "urban60_time", "farmlandpct60_time", "lnpop60_time"))
  control_cols <- c(control_cols, time_cols)
  cat(sprintf("    Controls: %s\n", paste(control_cols, collapse = ", ")))

  controls_mat <- as.matrix(sub[, control_cols, drop = FALSE])
  controls_mat[is.na(controls_mat)] <- 0

  # X = fsp + controls only (FEs absorbed by iterative demeaning)
  X_raw <- cbind(fsp, controls_mat)
  cat(sprintf("    X columns (before FE absorption): %d\n", ncol(X_raw)))

  # --- Absorb county + state×year + time FEs via iterative demeaning ---
  county_id    <- sub$county_id
  state_year   <- sub$state_year
  time_id      <- sub$time

  cat("    Absorbing FEs (county + state×year + time)...\n")
  Q_Yk_dm <- demean_multi_fe(Q_Yk, county_id, state_year, time_id)
  X_dm    <- demean_multi_fe(X_raw, county_id, state_year, time_id)

  # Drop near-zero variance columns (if any control is fully absorbed)
  col_var <- apply(X_dm, 2, var)
  keep_cols <- col_var > 1e-10
  if (any(!keep_cols)) {
    cat(sprintf("    Dropping %d absorbed columns\n", sum(!keep_cols)))
    X_dm <- X_dm[, keep_cols, drop = FALSE]
  }

  # Z = X (exogenous case)
  Z_dm <- X_dm

  # --- Estimate ---
  cat("    Running 2SLS (exogenous)...\n")
  fit_2sls <- estimate_2sls(Q_Yk_dm, X_dm, Z_dm, Q_GRID)

  cat("    Running D-IV...\n")
  fit_div <- estimate_div(Q_Yk_dm, X_dm, Z_dm, Q_GRID)

  # Extract FSP coefficient
  # After filtering, fsp should still be the first column of X_dm (column 1 of beta1)
  beta1_2sls <- if (is.matrix(fit_2sls$beta1)) fit_2sls$beta1[1, ] else fit_2sls$beta1
  beta1_div  <- if (is.matrix(fit_div$beta1))  fit_div$beta1[1, ]  else fit_div$beta1

  if (any(is.na(beta1_2sls))) {
    cat("    WARNING: 2SLS estimates contain NAs\n")
  }

  cat(sprintf("    FSP effect at u=0.05 (2SLS): %.1f grams\n", beta1_2sls[1]))
  cat(sprintf("    FSP effect at u=0.05 (D-IV): %.1f grams\n", beta1_div[1]))
  cat(sprintf("    FSP effect at u=0.50 (2SLS): %.1f grams\n", beta1_2sls[10]))
  cat(sprintf("    FSP effect at u=0.50 (D-IV): %.1f grams\n", beta1_div[10]))

  # --- Inference (clustered by county) ---
  cat("    Computing confidence bands (B=500, clustered by county)...\n")
  ci_pw <- div_pointwise_ci(Q_Yk_dm, X_dm, Z_dm, Q_GRID,
                             cluster = county_id, alpha = 0.05)
  # Use unprojected bootstrap only (projected is too memory-intensive for 80K groups)
  ci_ub <- div_uniform_cb(Q_Yk_dm, X_dm, Z_dm, Q_GRID,
                           cluster = county_id,
                           alpha = 0.05, B = 500L, projected = FALSE,
                           seed = 42L)

  # FSP coefficient is row 2 of (p+1) x Q matrix (row 1 = intercept)
  k <- 2L

  data.frame(
    quantile   = Q_GRID,
    beta1_2sls = beta1_2sls,
    beta1_div  = beta1_div,
    div_pw_lo  = ci_pw$ci_div$lo[k, ],
    div_pw_hi  = ci_pw$ci_div$hi[k, ],
    div_ub_lo  = ci_ub$ucb_div$lo[k, ],
    div_ub_hi  = ci_ub$ucb_div$hi[k, ],
    div_se     = ci_pw$se[k, ]
  )
}

# Process one race at a time to limit memory usage
mp_results <- list()
for (r in names(RACE_CODES)) {
  mp_results[[r]] <- run_mp_div(d, r)
  gc()  # free memory between races
}


# =============================================================================
# Figures: Replicate MP Figure 2
# =============================================================================

cat("\n=== Generating figures ===\n")

theme_mp <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.margin = margin(t = -5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

for (race in names(RACE_CODES)) {
  res <- mp_results[[race]]

  p <- ggplot(res, aes(x = quantile)) +
    # Uniform confidence band (lighter)
    geom_ribbon(aes(ymin = div_ub_lo, ymax = div_ub_hi),
                fill = "grey30", alpha = 0.12) +
    # Sandwich pointwise CI (darker)
    geom_ribbon(aes(ymin = div_pw_lo, ymax = div_pw_hi),
                fill = "grey30", alpha = 0.25) +
    # Zero reference line
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    # 2SLS point estimates
    geom_line(aes(y = beta1_2sls, color = "2SLS"), linewidth = 0.5) +
    geom_point(aes(y = beta1_2sls, color = "2SLS", shape = "2SLS"), size = 2) +
    # D-IV point estimates
    geom_line(aes(y = beta1_div, color = "D-IV"), linewidth = 0.5) +
    geom_point(aes(y = beta1_div, color = "D-IV", shape = "D-IV"), size = 2) +
    # Scales
    scale_color_manual(
      name = NULL,
      values = c("2SLS" = "grey50", "D-IV" = "black")
    ) +
    scale_shape_manual(
      name = NULL,
      values = c("2SLS" = 1, "D-IV" = 16)
    ) +
    scale_x_continuous(
      breaks = seq(0.1, 0.9, by = 0.1),
      limits = c(0.04, 0.96), expand = c(0, 0)
    ) +
    labs(x = "Quantiles",
         y = "Point Estimates",
         title = RACE_LABELS[race]) +
    theme_mp

  fig_path <- file.path(FIG_DIR, sprintf("mp_div_%s.pdf", race))
  ggsave(fig_path, p, width = 7, height = 5)
  cat(sprintf("  Saved: %s\n", fig_path))
}


# =============================================================================
# Summary
# =============================================================================

cat("\n=== Summary ===\n")
for (race in names(RACE_CODES)) {
  res <- mp_results[[race]]
  cat(sprintf("\n  %s:\n", RACE_LABELS[race]))
  cat(sprintf("    FSP effect at u=0.05: %.1f grams (2SLS) / %.1f grams (D-IV)\n",
              res$beta1_2sls[1], res$beta1_div[1]))
  cat(sprintf("    FSP effect at u=0.50: %.1f grams (2SLS) / %.1f grams (D-IV)\n",
              res$beta1_2sls[10], res$beta1_div[10]))
  cat(sprintf("    FSP effect at u=0.95: %.1f grams (2SLS) / %.1f grams (D-IV)\n",
              res$beta1_2sls[19], res$beta1_div[19]))
  cat(sprintf("    Max |D-IV - 2SLS|: %.1f grams\n",
              max(abs(res$beta1_div - res$beta1_2sls))))
  cat(sprintf("    Mean SE (sandwich): %.1f grams\n",
              mean(res$div_se)))
}

cat("\nDone.\n")
