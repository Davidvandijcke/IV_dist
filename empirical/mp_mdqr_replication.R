# =============================================================================
# Replicate MP (2025) Figure 2 using their mdqr package
# =============================================================================
#
# Uses the mdqr R package (github.com/martinapons/mdqr) to run the exact
# MP estimator: first-stage within-group QR, second-stage OLS/GMM with FEs.
#
# Formula: bweight ~ sex + mom_age + mom_age_sq + legit_bin |
#          fsp |
#          0 |
#          county_id + state_year + time |
#          group_id
#
# =============================================================================

PROJECT_ROOT <- tryCatch({
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}, error = function(e) getwd())

library(haven)
library(mdqr)
library(ggplot2)

DATA_DIR <- file.path(PROJECT_ROOT, "data", "out")
FIG_DIR <- file.path(PROJECT_ROOT, "empirical", "figures")

Q_GRID <- seq(0.05, 0.95, by = 0.05)

# =============================================================================
# Load individual-level data (black mothers only)
# =============================================================================

cat("Loading individual-level data (black mothers)...\n")
d <- read_dta(file.path(DATA_DIR, "mp_individual_black.dta"))
cat(sprintf("  %d observations, %d groups\n", nrow(d), length(unique(d$group_id))))

# =============================================================================
# Run mdqr — both stages at once
# =============================================================================
#
# Formula: y ~ exo | endo | instruments | FEs | group_ID
# mdqr uses fixest::feols internally for the second stage → FEs absorbed efficiently
#
# First stage: within-group QR of bweight on (sex, mom_age, mom_age_sq, legit_bin)
# Second stage: feols of fitted values on FSP + controls | county + state_year + time
#
# =============================================================================

cat("\nRunning mdqr (first + second stage)...\n")
n_cores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("  Cores: %d\n", n_cores))

fit <- mdqr(
  bweight ~ sex + mom_age + mom_age_sq + legit_bin + fsp +
    tranpcret + tranpcmed + tranpcvet + ripc +
    black60_time + urban60_time + farmlandpct60_time + lnpop60_time |
    0 | 0 |
    county_id + state_year + time |
    group_id,
  data = d,
  method = "ols",
  quantiles = Q_GRID,
  clustervar = "county_id",
  cores = n_cores,
  n_small = 25,
  run_time = TRUE
)

# =============================================================================
# Extract FSP coefficient
# =============================================================================

cat("\nExtracting results...\n")

# fit$results is a fixest_multi object — one fixest per quantile
fsp_coef <- numeric(length(Q_GRID))
fsp_se   <- numeric(length(Q_GRID))

for (i in seq_along(Q_GRID)) {
  res_i <- fit$results[[i]]
  ct <- fixest::coeftable(res_i)
  fsp_idx <- which(rownames(ct) == "fsp")
  if (length(fsp_idx) > 0) {
    fsp_coef[i] <- ct[fsp_idx, "Estimate"]
    fsp_se[i]   <- ct[fsp_idx, "Std. Error"]
  }
}

cat(sprintf("\nFSP effect (MP estimator, black mothers):\n"))
for (i in seq_along(Q_GRID)) {
  cat(sprintf("  u=%.2f: %6.1f grams (SE %5.1f)\n", Q_GRID[i], fsp_coef[i], fsp_se[i]))
}

# =============================================================================
# Figure: MP replication
# =============================================================================

cat("\n=== Generating figure ===\n")

plot_df <- data.frame(
  quantile = Q_GRID,
  estimate = fsp_coef,
  ci_lo = fsp_coef - 1.96 * fsp_se,
  ci_hi = fsp_coef + 1.96 * fsp_se
)

p <- ggplot(plot_df, aes(x = quantile)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "grey70", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(y = estimate), linewidth = 0.8) +
  geom_point(aes(y = estimate), size = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1),
                     limits = c(0.04, 0.96), expand = c(0, 0)) +
  labs(x = "Quantiles", y = "Point Estimates",
       title = "Black Mothers (MP Estimator)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggsave(file.path(FIG_DIR, "mp_mdqr_black.pdf"), p, width = 6, height = 4.5)
cat("  Saved: mp_mdqr_black.pdf\n")

cat("\nDone.\n")
