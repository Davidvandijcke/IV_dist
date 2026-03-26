# =============================================================================
# Full Decomposition: β(u) = δ(u) + composition + re-ranking
# =============================================================================
#
# Binary type split: Young (<24) vs Old (>=24)
# Better coverage than 4-way split — most county-trimester groups have both.
#
# Exercise 1: Decompose β(u) via Q^⊕ and Δ regressions
# Exercise 2: Type-specific contributions to direct effect
#
# =============================================================================

PROJECT_ROOT <- tryCatch({
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}, error = function(e) getwd())

library(haven)
library(ggplot2)

source(file.path(PROJECT_ROOT, "R", "utils.R"))
source(file.path(PROJECT_ROOT, "R", "estimators.R"))

demean_within <- function(mat, group_id) {
  mat <- as.matrix(mat)
  groups <- as.integer(as.factor(group_id))
  n_groups <- max(groups)
  group_sums <- rowsum(mat, groups)
  group_counts <- tabulate(groups, nbins = n_groups)
  group_means <- group_sums / group_counts
  mat - group_means[groups, , drop = FALSE]
}
demean_multi_fe <- function(mat, ..., tol = 1e-8, maxiter = 100) {
  fe_groups <- list(...)
  mat <- as.matrix(mat)
  for (iter in seq_len(maxiter)) {
    mat_old <- mat
    for (g in fe_groups) mat <- demean_within(mat, g)
    if (max(abs(mat - mat_old)) < tol) break
  }
  mat
}

DATA_DIR <- file.path(PROJECT_ROOT, "data", "out")
FIG_DIR <- file.path(PROJECT_ROOT, "empirical", "figures")

Q_GRID <- seq(0.05, 0.95, by = 0.05)
Q_LABELS <- paste0("q_", as.integer(Q_GRID * 100))
n_q <- length(Q_GRID)

TYPE_LABELS <- c("1" = "Old (>=24)", "2" = "Young (<24)")

control_cols <- c("tranpcret", "tranpcmed", "tranpcvet", "ripc",
                  "black60_time", "urban60_time", "farmlandpct60_time", "lnpop60_time")


# =============================================================================
# Load data
# =============================================================================

cat("Loading data...\n")
d_orig  <- read_dta(file.path(DATA_DIR, "mp_analysis_data.dta"))
d_types <- read_dta(file.path(DATA_DIR, "mp_binary_type_groups.dta"))

orig  <- d_orig[d_orig$mrace == 2, ]
types <- d_types[d_types$mrace == 2, ]

cat(sprintf("  Original groups: %d\n", nrow(orig)))
cat(sprintf("  Binary-type groups: %d\n", nrow(types)))
cat(sprintf("    Old (>=24): %d\n", sum(types$btype == 1)))
cat(sprintf("    Young (<24): %d\n", sum(types$btype == 2)))

rm(d_orig, d_types); gc()


# =============================================================================
# Exercise 1: β(u) on full sample
# =============================================================================

cat("\n=== β(u) from full sample ===\n")

Q_Yk <- as.matrix(orig[, Q_LABELS])
avail_ctrl <- intersect(control_cols, names(orig))
ctrls <- as.matrix(orig[, avail_ctrl, drop = FALSE]); ctrls[is.na(ctrls)] <- 0
X_raw <- cbind(orig$fsp, ctrls)

Q_dm <- demean_multi_fe(Q_Yk, orig$county_id, orig$state_year, orig$time)
X_dm <- demean_multi_fe(X_raw, orig$county_id, orig$state_year, orig$time)
cv <- apply(X_dm, 2, var); X_dm <- X_dm[, cv > 1e-10, drop = FALSE]

fit_beta <- estimate_2sls(Q_dm, X_dm, X_dm, Q_GRID)
beta_u <- if (is.matrix(fit_beta$beta1)) fit_beta$beta1[1, ] else fit_beta$beta1
cat(sprintf("  beta(0.05)=%.1f, beta(0.50)=%.1f, beta(0.95)=%.1f\n",
            beta_u[1], beta_u[10], beta_u[19]))


# =============================================================================
# Exercise 1: Compute Q^⊕_j and Δ_j using binary types
# =============================================================================

cat("\n=== Computing Q^oplus and Delta ===\n")

# Match type-specific QFs to original groups
orig$gkey  <- paste(orig$stfips, orig$countyfips, orig$year, orig$trimester, sep = "_")
types$gkey <- paste(types$stfips, types$countyfips, types$year, types$trimester, sep = "_")

Q_oplus  <- matrix(NA_real_, nrow(orig), n_q)
coverage <- rep(0, nrow(orig))

for (i in seq_len(nrow(orig))) {
  gk <- orig$gkey[i]
  type_rows <- types[types$gkey == gk, ]
  if (nrow(type_rows) == 0) next

  matched_births <- sum(type_rows$nbirths_type)
  coverage[i] <- matched_births / orig$nbirths[i]

  Q_oplus[i, ] <- 0
  for (r in seq_len(nrow(type_rows))) {
    wt <- type_rows$nbirths_type[r] / matched_births
    Q_oplus[i, ] <- Q_oplus[i, ] + wt * as.numeric(type_rows[r, Q_LABELS])
  }
}

# Filter: keep groups with >=70% coverage
good <- !is.na(Q_oplus[, 1]) & coverage >= 0.70
cat(sprintf("  Groups with >=70%% coverage: %d / %d (%.0f%%)\n",
            sum(good), nrow(orig), 100 * mean(good)))
cat(sprintf("  Mean coverage: %.1f%%\n", 100 * mean(coverage[good])))

# Restrict
orig_r    <- orig[good, ]
Q_Yk_r    <- Q_Yk[good, ]
Q_oplus_r <- Q_oplus[good, ]
ctrls_r   <- ctrls[good, ]
X_raw_r   <- cbind(orig_r$fsp, ctrls_r)

# Δ_j = Q_j - Q^⊕_j
Delta_r <- Q_Yk_r - Q_oplus_r

cat(sprintf("  Mean Delta(0.05)=%.1f, Delta(0.50)=%.1f\n",
            mean(Delta_r[, 1]), mean(Delta_r[, 10])))

# FE absorption on restricted sample
Q_dm_r     <- demean_multi_fe(Q_Yk_r, orig_r$county_id, orig_r$state_year, orig_r$time)
Q_oplus_dm <- demean_multi_fe(Q_oplus_r, orig_r$county_id, orig_r$state_year, orig_r$time)
Delta_dm   <- demean_multi_fe(Delta_r, orig_r$county_id, orig_r$state_year, orig_r$time)
X_dm_r     <- demean_multi_fe(X_raw_r, orig_r$county_id, orig_r$state_year, orig_r$time)
cv <- apply(X_dm_r, 2, var); X_dm_r <- X_dm_r[, cv > 1e-10, drop = FALSE]

# Regressions
fit_beta_r  <- estimate_2sls(Q_dm_r, X_dm_r, X_dm_r, Q_GRID)
fit_oplus_r <- estimate_2sls(Q_oplus_dm, X_dm_r, X_dm_r, Q_GRID)
fit_delta_r <- estimate_2sls(Delta_dm, X_dm_r, X_dm_r, Q_GRID)

beta_r     <- if (is.matrix(fit_beta_r$beta1)) fit_beta_r$beta1[1, ] else fit_beta_r$beta1
coef_oplus <- if (is.matrix(fit_oplus_r$beta1)) fit_oplus_r$beta1[1, ] else fit_oplus_r$beta1
coef_delta <- if (is.matrix(fit_delta_r$beta1)) fit_delta_r$beta1[1, ] else fit_delta_r$beta1

check <- max(abs(beta_r - coef_oplus - coef_delta))
cat(sprintf("  beta_r(0.05)=%.1f (restricted sample)\n", beta_r[1]))
cat(sprintf("  max|beta - (coef_Qoplus + coef_Delta)| = %.4f\n", check))


# =============================================================================
# Exercise 2: Type-specific D-IV
# =============================================================================

cat("\n=== Type-specific D-IV ===\n")

delta_k <- list()
n_k <- c()

for (tid in sort(unique(types$btype))) {
  sub_t <- types[types$btype == tid, ]
  M_t <- nrow(sub_t)
  n_k[as.character(tid)] <- M_t

  Q_t <- as.matrix(sub_t[, Q_LABELS])
  ctrl_t <- as.matrix(sub_t[, intersect(avail_ctrl, names(sub_t)), drop = FALSE])
  ctrl_t[is.na(ctrl_t)] <- 0
  X_t <- cbind(sub_t$fsp, ctrl_t)

  Q_t_dm <- demean_multi_fe(Q_t, sub_t$county_id, sub_t$state_year, sub_t$time)
  X_t_dm <- demean_multi_fe(X_t, sub_t$county_id, sub_t$state_year, sub_t$time)
  cv <- apply(X_t_dm, 2, var); X_t_dm <- X_t_dm[, cv > 1e-10, drop = FALSE]

  fit_t <- tryCatch(estimate_2sls(Q_t_dm, X_t_dm, X_t_dm, Q_GRID), error = function(e) NULL)
  if (!is.null(fit_t)) {
    d_t <- if (is.matrix(fit_t$beta1)) fit_t$beta1[1, ] else fit_t$beta1
    delta_k[[as.character(tid)]] <- d_t
    cat(sprintf("  %s: %d grps, delta(0.05)=%.1f, delta(0.50)=%.1f, delta(0.95)=%.1f\n",
                TYPE_LABELS[as.character(tid)], M_t, d_t[1], d_t[10], d_t[19]))
  }
}

valid_types <- names(delta_k)
pi_k <- n_k[valid_types] / sum(n_k[valid_types])
delta_pop <- rep(0, n_q)
for (tid in valid_types) delta_pop <- delta_pop + pi_k[tid] * delta_k[[tid]]

cat(sprintf("  Pop-avg delta(0.05)=%.1f, delta(0.50)=%.1f\n", delta_pop[1], delta_pop[10]))

# Composition = coef_oplus - delta_pop
composition_u <- coef_oplus - delta_pop


# =============================================================================
# Summary table
# =============================================================================

cat(sprintf("\n=== Full Decomposition (restricted sample, binary types) ===\n"))
cat(sprintf("  %-30s  u=0.05   u=0.50   u=0.95\n", ""))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "beta(u) [total, restricted]", beta_r[1], beta_r[10], beta_r[19]))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "coef Q^oplus [direct+comp]", coef_oplus[1], coef_oplus[10], coef_oplus[19]))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "  delta(u) [pop-avg direct]", delta_pop[1], delta_pop[10], delta_pop[19]))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "  Composition", composition_u[1], composition_u[10], composition_u[19]))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "coef Delta [re-ranking]", coef_delta[1], coef_delta[10], coef_delta[19]))
cat(sprintf("  %-30s  %6.1f   %6.1f   %6.1f\n", "Sum (should = beta_r)",
            delta_pop[1] + composition_u[1] + coef_delta[1],
            delta_pop[10] + composition_u[10] + coef_delta[10],
            delta_pop[19] + composition_u[19] + coef_delta[19]))
cat(sprintf("\n  %-30s  %6.1f   %6.1f   %6.1f\n", "beta(u) [full sample]", beta_u[1], beta_u[10], beta_u[19]))


# =============================================================================
# Figures
# =============================================================================

cat("\n=== Generating figures ===\n")

# Figure 1: Decomposition
decomp_df <- data.frame(
  quantile = rep(Q_GRID, 4),
  value = c(beta_r, delta_pop, composition_u, coef_delta),
  component = factor(rep(c("Total (D-IV)", "Direct (within-type)",
                           "Composition", "Re-ranking"),
                         each = n_q),
                     levels = c("Total (D-IV)", "Direct (within-type)",
                                "Composition", "Re-ranking"))
)

p1 <- ggplot(decomp_df, aes(x = quantile, y = value,
                              color = component, linetype = component)) +
  geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Total (D-IV)" = "black",
                                "Direct (within-type)" = "#2980b9",
                                "Composition" = "#27ae60",
                                "Re-ranking" = "#c0392b")) +
  scale_linetype_manual(values = c("Total (D-IV)" = "solid",
                                   "Direct (within-type)" = "solid",
                                   "Composition" = "dashed",
                                   "Re-ranking" = "dotdash")) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1),
                     limits = c(0.04, 0.96), expand = c(0, 0)) +
  labs(x = "Quantile", y = "Effect on Birth Weight (grams)",
       title = "Decomposition of FSP Effect (Black Mothers, binary type split)",
       color = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "mp_full_decomposition_black.pdf"), p1, width = 7, height = 5)
cat("  Saved: mp_full_decomposition_black.pdf\n")

# Figure 2: Type-specific effects
type_df <- data.frame(
  quantile = rep(Q_GRID, length(valid_types) + 1),
  value = c(do.call(c, delta_k[valid_types]), delta_pop),
  type = factor(c(rep(TYPE_LABELS[valid_types], each = n_q), rep("Pop average", n_q)),
                levels = c(TYPE_LABELS[valid_types], "Pop average"))
)

p2 <- ggplot(type_df, aes(x = quantile, y = value, color = type, linewidth = type)) +
  geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
  geom_line() +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Old (>=24)" = "#e74c3c", "Young (<24)" = "#3498db",
                                "Pop average" = "black")) +
  scale_linewidth_manual(values = c("Old (>=24)" = 0.7, "Young (<24)" = 0.7,
                                    "Pop average" = 1.1)) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1),
                     limits = c(0.04, 0.96), expand = c(0, 0)) +
  labs(x = "Quantile", y = "FSP Effect (grams)",
       title = "Type-Specific Direct Effects (Black Mothers)",
       color = NULL, linewidth = NULL) +
  guides(linewidth = "none") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "mp_type_effects_black.pdf"), p2, width = 7, height = 5)
cat("  Saved: mp_type_effects_black.pdf\n")

cat("\nDone.\n")
