# =============================================================================
# Full Decomposition: β(u) = δ(u) + composition + re-ranking
# =============================================================================
#
# Uses mdqr fitted values to compute Q^⊕_j(u) exactly, then decomposes
# the D-IV estimate β(u) into direct + composition + re-ranking.
#
# Strategy:
#   1. Run mdqr first stage → get individual fitted values ŷ_ij(u)
#   2. Q^⊕_j(u) = mean of ŷ_ij(u) within each group j (composition-fixed QF)
#   3. Q_j(u) = empirical group quantile (what D-IV uses)
#   4. Δ_j(u) = Q_j(u) - Q^⊕_j(u) (re-ranking gap)
#   5. Regress Q^⊕ and Δ on FSP with FEs via fixest → decomposition coefficients
#   6. δ(u) from mdqr second stage (= MP's estimator)
#   7. Composition = coef_on_Q^⊕ - δ(u)
#
# =============================================================================

PROJECT_ROOT <- tryCatch({
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}, error = function(e) getwd())

library(haven)
library(mdqr)
library(fixest)
library(ggplot2)

DATA_DIR <- file.path(PROJECT_ROOT, "data", "out")
FIG_DIR <- file.path(PROJECT_ROOT, "empirical", "figures")

Q_GRID <- seq(0.05, 0.95, by = 0.05)
n_q <- length(Q_GRID)
Q_LABELS <- paste0("q_", as.integer(Q_GRID * 100))


# =============================================================================
# Step 1: Load data and run mdqr (first + second stage, save fitted values)
# =============================================================================

cat("Loading individual data (black mothers)...\n")
d <- read_dta(file.path(DATA_DIR, "mp_individual_black.dta"))
cat(sprintf("  %d obs, %d groups\n", nrow(d), length(unique(d$group_id))))

cat("\nRunning mdqr (first + second stage, saving fitted values)...\n")
n_cores <- max(1L, parallel::detectCores() - 1L)

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
  run_time = TRUE,
  save_data = TRUE
)

# Extract δ(u) = MP's FSP coefficient
delta_mp <- numeric(n_q)
delta_mp_se <- numeric(n_q)
for (i in seq_along(Q_GRID)) {
  ct <- fixest::coeftable(fit$results[[i]])
  delta_mp[i] <- ct["fsp", "Estimate"]
  delta_mp_se[i] <- ct["fsp", "Std. Error"]
}
cat(sprintf("  delta_MP(0.05) = %.1f, delta_MP(0.50) = %.1f\n", delta_mp[1], delta_mp[10]))


# =============================================================================
# Step 2: Compute Q^⊕_j(u) = group mean of fitted values
# =============================================================================

cat("\nComputing Q^oplus_j and Q_j from saved data...\n")

# The saved data has: bweight, fitted_*, group, county_id, state_year, time, fsp, n
fv_data <- fit$data
fv_cols <- paste0("fitted_", Q_GRID)

gid <- fv_data$group
groups <- sort(unique(gid))
M <- length(groups)
gid_int <- match(gid, groups)
cat(sprintf("  Groups: %d, Individuals: %d\n", M, nrow(fv_data)))

fv_mat <- as.matrix(fv_data[, fv_cols])
bw_vec <- fv_data$bweight

# Q^⊕_j(u) = group mean of fitted values (composition-fixed quantile)
Q_oplus <- matrix(NA_real_, M, n_q)
for (q in seq_len(n_q)) {
  Q_oplus[, q] <- tapply(fv_mat[, q], gid_int, mean, na.rm = TRUE)
}

# Q_j(u) = group quantile of actual birth weights (realized quantile)
Q_j <- matrix(NA_real_, M, n_q)
for (q in seq_len(n_q)) {
  Q_j[, q] <- tapply(bw_vec, gid_int, quantile, probs = Q_GRID[q], type = 7, na.rm = TRUE)
}

# Group-level variables (constant within group)
group_fsp       <- as.numeric(tapply(fv_data$fsp, gid_int, function(x) x[1]))
group_county    <- as.numeric(tapply(fv_data$county_id, gid_int, function(x) x[1]))
group_stateyear <- as.numeric(tapply(fv_data$state_year, gid_int, function(x) x[1]))
group_time      <- as.numeric(tapply(fv_data$time, gid_int, function(x) x[1]))

matched <- M
cat(sprintf("  Q^oplus mean(0.05)=%.0f, Q_j mean(0.05)=%.0f\n",
            mean(Q_oplus[, 1]), mean(Q_j[, 1])))

rm(fv_data, fv_mat, bw_vec, fit); gc()

# Δ_j(u) = Q_j(u) - Q^⊕_j(u)
Delta_j <- Q_j - Q_oplus

cat(sprintf("  Mean Delta(0.05) = %.1f, Delta(0.50) = %.1f\n",
            mean(Delta_j[, 1], na.rm = TRUE), mean(Delta_j[, 10], na.rm = TRUE)))


# =============================================================================
# Step 4: Regress Q^⊕ and Δ on FSP with FEs via fixest
# =============================================================================

cat("\nRegressing Q^oplus and Delta on FSP...\n")

# Build group-level data frame
gdf <- data.frame(
  fsp = group_fsp,
  county_id = as.factor(group_county),
  state_year = as.factor(group_stateyear),
  time = as.factor(group_time)
)

good <- !is.na(Q_j[, 1]) & !is.na(Q_oplus[, 1]) & group_county > 0
gdf <- gdf[good, ]
Q_j_g <- Q_j[good, ]
Q_oplus_g <- Q_oplus[good, ]
Delta_g <- Delta_j[good, ]

cat(sprintf("  Groups for regression: %d\n", nrow(gdf)))

# β(u): D-IV-style regression of Q_j on FSP
beta_u <- numeric(n_q)
# coef on Q^⊕: direct + composition
coef_oplus <- numeric(n_q)
# coef on Δ: re-ranking
coef_delta <- numeric(n_q)

for (q in seq_len(n_q)) {
  gdf$y_q <- Q_j_g[, q]
  gdf$y_oplus <- Q_oplus_g[, q]
  gdf$y_delta <- Delta_g[, q]

  fit_q <- feols(y_q ~ fsp | county_id + state_year + time,
                 data = gdf, cluster = ~county_id)
  fit_o <- feols(y_oplus ~ fsp | county_id + state_year + time,
                 data = gdf, cluster = ~county_id)
  fit_d <- feols(y_delta ~ fsp | county_id + state_year + time,
                 data = gdf, cluster = ~county_id)

  beta_u[q] <- coef(fit_q)["fsp"]
  coef_oplus[q] <- coef(fit_o)["fsp"]
  coef_delta[q] <- coef(fit_d)["fsp"]
}

# Composition = coef_oplus - δ_MP
composition <- coef_oplus - delta_mp

# Verify: β ≈ coef_oplus + coef_delta
check <- max(abs(beta_u - coef_oplus - coef_delta))
cat(sprintf("  max|beta - (coef_Qoplus + coef_Delta)| = %.4f\n", check))


# =============================================================================
# Step 5: Summary table
# =============================================================================

cat(sprintf("\n=== Full Decomposition (Black Mothers) ===\n"))
cat(sprintf("  %-30s  u=0.05   u=0.50   u=0.95\n", ""))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "beta(u) [D-IV total]", beta_u[1], beta_u[10], beta_u[19]))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "delta(u) [MP direct]", delta_mp[1], delta_mp[10], delta_mp[19]))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "coef on Q^oplus [direct+comp]", coef_oplus[1], coef_oplus[10], coef_oplus[19]))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "Composition [Q^oplus - delta]", composition[1], composition[10], composition[19]))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "Re-ranking [coef on Delta]", coef_delta[1], coef_delta[10], coef_delta[19]))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "Gap [beta - delta]",
            beta_u[1] - delta_mp[1], beta_u[10] - delta_mp[10], beta_u[19] - delta_mp[19]))


# =============================================================================
# Figure: Decomposition
# =============================================================================

cat("\n=== Generating figures ===\n")

decomp_df <- data.frame(
  quantile = rep(Q_GRID, 4),
  value = c(beta_u, delta_mp, composition, coef_delta),
  component = factor(rep(c("Total (D-IV)", "Direct (MP)",
                           "Composition", "Re-ranking"),
                         each = n_q),
                     levels = c("Total (D-IV)", "Direct (MP)",
                                "Composition", "Re-ranking"))
)

p <- ggplot(decomp_df, aes(x = quantile, y = value,
                            color = component, linetype = component)) +
  geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Total (D-IV)" = "black",
                                "Direct (MP)" = "#2980b9",
                                "Composition" = "#27ae60",
                                "Re-ranking" = "#c0392b")) +
  scale_linetype_manual(values = c("Total (D-IV)" = "solid",
                                   "Direct (MP)" = "solid",
                                   "Composition" = "dashed",
                                   "Re-ranking" = "dotdash")) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1),
                     limits = c(0.04, 0.96), expand = c(0, 0)) +
  labs(x = "Quantile", y = "Effect on Birth Weight (grams)",
       title = "Decomposition of FSP Effect (Black Mothers)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "mp_full_decomposition_black.pdf"), p, width = 7, height = 5)
cat("  Saved: mp_full_decomposition_black.pdf\n")

cat("\nDone.\n")
