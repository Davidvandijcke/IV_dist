# =============================================================================
# Type Decomposition of Direct Effect + Rank Invariance Diagnostic
# =============================================================================
#
# 1. Run mdqr separately by type (young/old × legit/illegit) to get
#    type-specific δ_k(u). Show which types drive the large left-tail effect.
#
# 2. Compute density weights: at each quantile u, which types are marginal
#    at the group quantile cutoff?
#
# 3. Rank invariance diagnostic: check for non-monotonicity in type-specific
#    QTEs, which would indicate rank violations.
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
FIG_DIR  <- file.path(PROJECT_ROOT, "empirical", "figures")

Q_GRID <- seq(0.05, 0.95, by = 0.05)
n_q    <- length(Q_GRID)

cat("Loading individual data (black mothers)...\n")
d <- read_dta(file.path(DATA_DIR, "mp_individual_black.dta"))

# Create type: young (<24) × illegit
d$mom_young <- as.integer(d$mom_age < 24)
d$illegit   <- as.integer(d$legit_bin == 0)
d$type_id   <- 1L + d$mom_young + 2L * d$illegit
# 1=Old+Legit, 2=Young+Legit, 3=Old+Illegit, 4=Young+Illegit

TYPE_LABELS <- c("Old+Legit", "Young+Legit", "Old+Illegit", "Young+Illegit")
names(TYPE_LABELS) <- as.character(1:4)

type_counts <- table(d$type_id)
cat("Type distribution:\n")
for (tid in 1:4) {
  cat(sprintf("  %d (%s): %d (%.1f%%)\n",
              tid, TYPE_LABELS[as.character(tid)],
              type_counts[as.character(tid)],
              100 * type_counts[as.character(tid)] / nrow(d)))
}

n_cores <- max(1L, parallel::detectCores() - 1L)


# =============================================================================
# Part 1: Type-specific mdqr estimates δ_k(u)
# =============================================================================

cat("\n=== Part 1: Type-specific δ_k(u) via mdqr ===\n")

delta_k    <- list()
delta_k_se <- list()
n_k        <- c()

for (tid in 1:4) {
  d_t <- d[d$type_id == tid, ]
  n_t <- nrow(d_t)
  n_groups_t <- length(unique(d_t$group_id))
  n_k[as.character(tid)] <- n_t

  if (n_groups_t < 500) {
    cat(sprintf("  Type %d (%s): %d obs, %d groups — too few, skipping\n",
                tid, TYPE_LABELS[as.character(tid)], n_t, n_groups_t))
    next
  }

  cat(sprintf("  Type %d (%s): %d obs, %d groups — running mdqr...\n",
              tid, TYPE_LABELS[as.character(tid)], n_t, n_groups_t))

  # For type-specific QR, the first stage has NO individual covariates
  # (since we're already conditioning on type). Just intercept per group.
  # Second stage: regress intercepts on fsp + controls with FEs.
  fit_t <- tryCatch(
    mdqr(
      bweight ~ fsp +
        tranpcret + tranpcmed + tranpcvet + ripc +
        black60_time + urban60_time + farmlandpct60_time + lnpop60_time |
        0 | 0 |
        county_id + state_year + time |
        group_id,
      data = d_t,
      method = "ols",
      quantiles = Q_GRID,
      clustervar = "county_id",
      cores = n_cores,
      n_small = 25,
      run_time = TRUE
    ),
    error = function(e) { cat(sprintf("    Error: %s\n", e$message)); NULL }
  )

  if (!is.null(fit_t)) {
    coef_t <- se_t <- numeric(n_q)
    for (i in seq_along(Q_GRID)) {
      ct <- fixest::coeftable(fit_t$results[[i]])
      if ("fsp" %in% rownames(ct)) {
        coef_t[i] <- ct["fsp", "Estimate"]
        se_t[i]   <- ct["fsp", "Std. Error"]
      }
    }
    delta_k[[as.character(tid)]]    <- coef_t
    delta_k_se[[as.character(tid)]] <- se_t
    cat(sprintf("    δ(0.05)=%.1f, δ(0.50)=%.1f, δ(0.95)=%.1f\n",
                coef_t[1], coef_t[10], coef_t[19]))
  }

  rm(d_t, fit_t); gc()
}


# =============================================================================
# Part 2: Density weights — which types drive the pooled 5th percentile?
# =============================================================================

cat("\n=== Part 2: Density weights ===\n")

# Pooled quantiles (from all black mothers)
pooled_q <- quantile(d$bweight, probs = Q_GRID, type = 7)
bw_h <- bw.nrd0(d$bweight)

# KDE per type, evaluate at pooled quantiles
valid_types <- names(delta_k)
pi_k <- n_k[valid_types] / sum(n_k[valid_types])

w_k <- matrix(0, length(valid_types), n_q)
rownames(w_k) <- valid_types

for (i in seq_along(valid_types)) {
  tid <- as.numeric(valid_types[i])
  bw_t <- d$bweight[d$type_id == tid]
  kde_t <- density(bw_t, bw = bw_h, n = 2048,
                   from = min(d$bweight) - 500, to = max(d$bweight) + 500)
  f_t <- approx(kde_t$x, kde_t$y, xout = pooled_q, rule = 2)$y
  w_k[i, ] <- pi_k[valid_types[i]] * f_t
}

# Normalize
for (q in seq_len(n_q)) {
  s <- sum(w_k[, q])
  if (s > 0) w_k[, q] <- w_k[, q] / s
}

cat("Density weights at u=0.05, 0.50, 0.95:\n")
cat(sprintf("  %-15s  %6s  %7s  %7s  %7s\n", "Type", "pi", "w(0.05)", "w(0.50)", "w(0.95)"))
for (i in seq_along(valid_types)) {
  tid <- valid_types[i]
  cat(sprintf("  %-15s  %6.3f  %7.3f  %7.3f  %7.3f\n",
              TYPE_LABELS[tid], pi_k[tid], w_k[i, 1], w_k[i, 10], w_k[i, 19]))
}

# Density-weighted direct effect
delta_density_wt <- rep(0, n_q)
delta_pop_avg    <- rep(0, n_q)
for (i in seq_along(valid_types)) {
  delta_density_wt <- delta_density_wt + w_k[i, ] * delta_k[[valid_types[i]]]
  delta_pop_avg    <- delta_pop_avg + pi_k[valid_types[i]] * delta_k[[valid_types[i]]]
}

cat(sprintf("\nDensity-wtd δ(0.05) = %.1f, Pop-avg δ(0.05) = %.1f\n",
            delta_density_wt[1], delta_pop_avg[1]))

rm(d); gc()


# =============================================================================
# Part 3: Rank invariance diagnostic
# =============================================================================

cat("\n=== Part 3: Rank invariance diagnostic ===\n")

# Under rank invariance, δ_k(u) should be monotone non-increasing in u
# (or at least smooth). Non-monotonicity suggests rank violations.
cat("Checking monotonicity of type-specific QTEs:\n")
for (tid in valid_types) {
  d_t <- delta_k[[tid]]
  # Count sign changes in first differences
  diffs <- diff(d_t)
  n_increase <- sum(diffs > 0)
  n_decrease <- sum(diffs < 0)
  # Correlation with quantile index (should be negative if decreasing)
  cor_tau <- cor(Q_GRID, d_t)
  cat(sprintf("  %s: cor(tau, delta)=%.2f, increases=%d, decreases=%d\n",
              TYPE_LABELS[tid], cor_tau, n_increase, n_decrease))
}


# =============================================================================
# Figures
# =============================================================================

cat("\n=== Generating figures ===\n")

# --- Figure 1: Type-specific QTEs ---
type_df <- do.call(rbind, lapply(valid_types, function(tid) {
  data.frame(
    quantile = Q_GRID,
    estimate = delta_k[[tid]],
    ci_lo = delta_k[[tid]] - 1.96 * delta_k_se[[tid]],
    ci_hi = delta_k[[tid]] + 1.96 * delta_k_se[[tid]],
    type = TYPE_LABELS[tid],
    stringsAsFactors = FALSE
  )
}))
type_df$type <- factor(type_df$type, levels = TYPE_LABELS[valid_types])

p1 <- ggplot(type_df, aes(x = quantile, color = type, fill = type)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(y = estimate), linewidth = 0.8) +
  geom_point(aes(y = estimate), size = 1.5) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1),
                     limits = c(0.04, 0.96), expand = c(0, 0)) +
  labs(x = "Quantile", y = "FSP Effect (grams)",
       title = "Type-Specific Conditional QTE (Black Mothers)",
       color = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "mp_type_qte_black.pdf"), p1, width = 7, height = 5)
cat("  Saved: mp_type_qte_black.pdf\n")


# --- Figure 2: Density-weighted contributions ---
contrib_df <- do.call(rbind, lapply(valid_types, function(tid) {
  i <- which(valid_types == tid)
  data.frame(
    quantile = Q_GRID,
    contribution = w_k[i, ] * delta_k[[tid]],
    weight = w_k[i, ],
    type = TYPE_LABELS[tid],
    stringsAsFactors = FALSE
  )
}))
contrib_df$type <- factor(contrib_df$type, levels = TYPE_LABELS[valid_types])

p2 <- ggplot(contrib_df, aes(x = quantile, y = contribution, color = type)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  # Add total as black line
  geom_line(data = data.frame(quantile = Q_GRID, contribution = delta_density_wt),
            aes(x = quantile, y = contribution),
            inherit.aes = FALSE, linewidth = 1.1, color = "black") +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1),
                     limits = c(0.04, 0.96), expand = c(0, 0)) +
  labs(x = "Quantile", y = "Density-Weighted Contribution (grams)",
       title = "Type Contributions to Direct Effect (Black Mothers)",
       subtitle = "Black line = total density-weighted direct effect",
       color = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "mp_type_contributions_black.pdf"), p2, width = 7, height = 5)
cat("  Saved: mp_type_contributions_black.pdf\n")


# --- Figure 3: Density weights across quantiles ---
w_df <- do.call(rbind, lapply(seq_along(valid_types), function(i) {
  tid <- valid_types[i]
  data.frame(
    quantile = Q_GRID,
    weight = w_k[i, ],
    pop_share = pi_k[tid],
    type = TYPE_LABELS[tid],
    stringsAsFactors = FALSE
  )
}))
w_df$type <- factor(w_df$type, levels = TYPE_LABELS[valid_types])

p3 <- ggplot(w_df, aes(x = quantile)) +
  geom_line(aes(y = weight, color = type), linewidth = 0.8) +
  geom_line(aes(y = pop_share, color = type), linetype = "dashed", linewidth = 0.4) +
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.2),
                     limits = c(0.04, 0.96), expand = c(0, 0)) +
  labs(x = "Quantile", y = "Density Weight",
       title = "Density Weights by Type (Black Mothers)",
       subtitle = "Solid = density weight; Dashed = population share",
       color = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "mp_density_weights_types_black.pdf"), p3, width = 7, height = 5)
cat("  Saved: mp_density_weights_types_black.pdf\n")


# --- Summary ---
cat("\n=== Summary ===\n")
cat(sprintf("  Overall MP δ(0.05) ≈ 27g\n"))
cat(sprintf("  Pop-avg of type-specific δ_k(0.05) = %.1f\n", delta_pop_avg[1]))
cat(sprintf("  Density-wtd δ_k(0.05) = %.1f\n", delta_density_wt[1]))
cat(sprintf("  Dominant type at u=0.05: %s (w=%.2f, δ=%.1f)\n",
            TYPE_LABELS[valid_types[which.max(w_k[, 1])]],
            max(w_k[, 1]),
            delta_k[[valid_types[which.max(w_k[, 1])]]][1]))

cat("\nDone.\n")
