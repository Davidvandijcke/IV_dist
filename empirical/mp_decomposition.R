# =============================================================================
# Decomposition: D-IV vs MP Estimates
# =============================================================================
#
# Decomposes the gap between our D-IV estimate β(u) and MP's within-type
# estimate δ(u) using the paper's Appendix C framework:
#
#   β(u) = δ(u) + composition + re-ranking
#
# Components:
#   1. δ(u)  — type-specific D-IV, averaged across types (direct effect)
#   2. Composition — FSP changes who gives birth (type shares shift)
#   3. Re-ranking  — density weighting: types contribute differently at
#                    different quantiles of the pooled distribution
#
# Usage:
#   cd IV_dist && Rscript empirical/mp_decomposition.R
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
source(file.path(PROJECT_ROOT, "R", "inference.R"))

# Reuse the FE absorption function from mp_replication.R
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

# --- Paths ---
DATA_DIR <- file.path(PROJECT_ROOT, "data", "out")
FIG_DIR <- file.path(PROJECT_ROOT, "empirical", "figures")

Q_GRID <- seq(0.05, 0.95, by = 0.05)
Q_LABELS <- paste0("q_", as.integer(Q_GRID * 100))

TYPE_LABELS <- c("1" = "Old+Legit", "2" = "Young+Legit",
                 "3" = "Old+Illegit", "4" = "Young+Illegit")


# =============================================================================
# Load data
# =============================================================================

cat("Loading datasets...\n")

# Original D-IV estimates (from mp_analysis_data.dta via mp_replication.R)
d_orig <- read_dta(file.path(DATA_DIR, "mp_analysis_data.dta"))

# Type-specific groups
d_types <- read_dta(file.path(DATA_DIR, "mp_type_groups.dta"))

# Type composition
d_comp <- read_dta(file.path(DATA_DIR, "mp_type_composition.dta"))

# Density data (10% sample of individuals)
d_dens <- read_dta(file.path(DATA_DIR, "mp_type_density.dta"))

cat(sprintf("  Original groups: %d\n", nrow(d_orig)))
cat(sprintf("  Type-specific groups: %d\n", nrow(d_types)))
cat(sprintf("  Composition records: %d\n", nrow(d_comp)))
cat(sprintf("  Density sample: %d individuals\n", nrow(d_dens)))


# =============================================================================
# Run decomposition for each race
# =============================================================================

run_decomposition <- function(race_code, race_label) {
  cat(sprintf("\n=== %s ===\n", race_label))

  # --- Original D-IV estimate β(u) ---
  sub_orig <- d_orig[d_orig$mrace == race_code, ]
  Q_Yk_orig <- as.matrix(sub_orig[, Q_LABELS])

  control_cols <- intersect(names(sub_orig),
    c("tranpcret", "tranpcmed", "tranpcvet", "ripc",
      "black60_time", "urban60_time", "farmlandpct60_time", "lnpop60_time"))
  controls_orig <- as.matrix(sub_orig[, control_cols, drop = FALSE])
  controls_orig[is.na(controls_orig)] <- 0
  X_orig <- cbind(sub_orig$fsp, controls_orig)

  Q_dm <- demean_multi_fe(Q_Yk_orig, sub_orig$county_id, sub_orig$state_year, sub_orig$time)
  X_dm <- demean_multi_fe(X_orig, sub_orig$county_id, sub_orig$state_year, sub_orig$time)
  col_var <- apply(X_dm, 2, var); X_dm <- X_dm[, col_var > 1e-10, drop = FALSE]

  fit_div_orig <- estimate_div(Q_dm, X_dm, X_dm, Q_GRID)
  beta_u <- if (is.matrix(fit_div_orig$beta1)) fit_div_orig$beta1[1, ] else fit_div_orig$beta1
  cat(sprintf("  β(u=0.05) = %.1f, β(u=0.50) = %.1f\n", beta_u[1], beta_u[10]))

  rm(Q_Yk_orig, Q_dm); gc()

  # --- Type-specific D-IV estimates δ_k(u) ---
  sub_types <- d_types[d_types$mrace == race_code, ]
  types_present <- sort(unique(sub_types$type_id))
  cat(sprintf("  Types present: %s\n", paste(types_present, collapse = ", ")))

  delta_k <- list()  # type_id -> Q-vector of D-IV estimates
  n_k <- c()         # type counts

  for (tid in types_present) {
    sub_t <- sub_types[sub_types$type_id == tid, ]
    M_t <- nrow(sub_t)
    n_k[as.character(tid)] <- M_t

    if (M_t < 50) {
      cat(sprintf("    Type %d (%s): only %d groups, skipping\n",
                  tid, TYPE_LABELS[as.character(tid)], M_t))
      next
    }

    Q_t <- as.matrix(sub_t[, Q_LABELS])
    ctrl_t <- as.matrix(sub_t[, intersect(control_cols, names(sub_t)), drop = FALSE])
    ctrl_t[is.na(ctrl_t)] <- 0
    X_t <- cbind(sub_t$fsp, ctrl_t)

    Q_t_dm <- demean_multi_fe(Q_t, sub_t$county_id, sub_t$state_year, sub_t$time)
    X_t_dm <- demean_multi_fe(X_t, sub_t$county_id, sub_t$state_year, sub_t$time)
    cv <- apply(X_t_dm, 2, var); X_t_dm <- X_t_dm[, cv > 1e-10, drop = FALSE]

    fit_t <- tryCatch(
      estimate_div(Q_t_dm, X_t_dm, X_t_dm, Q_GRID),
      error = function(e) { cat(sprintf("    Type %d: error - %s\n", tid, e$message)); NULL }
    )

    if (!is.null(fit_t)) {
      d_t <- if (is.matrix(fit_t$beta1)) fit_t$beta1[1, ] else fit_t$beta1
      delta_k[[as.character(tid)]] <- d_t
      cat(sprintf("    Type %d (%s): %d groups, δ(0.05) = %.1f, δ(0.50) = %.1f\n",
                  tid, TYPE_LABELS[as.character(tid)], M_t, d_t[1], d_t[10]))
    }

    rm(Q_t, Q_t_dm, X_t, X_t_dm); gc()
  }

  # --- Population-weighted average δ(u) = Σ_k π_k δ_k(u) ---
  valid_types <- names(delta_k)
  pi_k <- n_k[valid_types] / sum(n_k[valid_types])
  delta_avg <- rep(0, length(Q_GRID))
  for (tid in valid_types) {
    delta_avg <- delta_avg + pi_k[tid] * delta_k[[tid]]
  }
  cat(sprintf("  Population-avg δ(0.05) = %.1f, δ(0.50) = %.1f\n",
              delta_avg[1], delta_avg[10]))

  # --- Density weights w_k(u) ---
  cat("  Computing density weights...\n")
  dens_sub <- d_dens[d_dens$mrace == race_code, ]

  # Pooled quantiles (from original group QFs, average across groups)
  pooled_quantiles <- colMeans(as.matrix(d_orig[d_orig$mrace == race_code, Q_LABELS]))

  # KDE for each type
  density_at_quantile <- matrix(0, length(valid_types), length(Q_GRID))
  rownames(density_at_quantile) <- valid_types

  for (i in seq_along(valid_types)) {
    tid <- valid_types[i]
    bw_vals <- dens_sub$bweight[dens_sub$type_id == as.numeric(tid)]
    if (length(bw_vals) < 100) next
    kde <- density(bw_vals, n = 1024, from = min(bw_vals) - 100, to = max(bw_vals) + 100)
    # Evaluate at pooled quantile values
    density_at_quantile[i, ] <- approx(kde$x, kde$y, xout = pooled_quantiles, rule = 2)$y
  }

  # Density weights: w_k(u) = π_k f_k(Q(u)) / Σ_l π_l f_l(Q(u))
  w_k <- matrix(0, length(valid_types), length(Q_GRID))
  for (q in seq_along(Q_GRID)) {
    numerator <- pi_k[valid_types] * density_at_quantile[, q]
    denom <- sum(numerator)
    if (denom > 0) w_k[, q] <- numerator / denom
  }
  rownames(w_k) <- valid_types

  cat("  Density weights at u=0.05:\n")
  for (i in seq_along(valid_types)) {
    cat(sprintf("    %s: π=%.3f, w(0.05)=%.3f, w(0.50)=%.3f, w(0.95)=%.3f\n",
                TYPE_LABELS[valid_types[i]], pi_k[valid_types[i]],
                w_k[i, 1], w_k[i, 10], w_k[i, 19]))
  }

  # --- Density-weighted effect: Σ_k w_k(u) δ_k(u) ---
  delta_density_weighted <- rep(0, length(Q_GRID))
  for (i in seq_along(valid_types)) {
    delta_density_weighted <- delta_density_weighted + w_k[i, ] * delta_k[[valid_types[i]]]
  }
  cat(sprintf("  Density-weighted δ(0.05) = %.1f, δ(0.50) = %.1f\n",
              delta_density_weighted[1], delta_density_weighted[10]))

  # --- Composition effect ---
  # Regress type shares on FSP (with FEs) to get E[W̄(1) - W̄(0)]
  comp_sub <- d_comp[d_comp$mrace == race_code, ]
  share_cols <- paste0("type_share", types_present)
  avail_share_cols <- intersect(share_cols, names(comp_sub))

  composition_effect <- rep(0, length(Q_GRID))
  if (length(avail_share_cols) > 0) {
    for (sc in avail_share_cols) {
      tid_num <- gsub("type_share", "", sc)
      if (!(tid_num %in% valid_types)) next

      # Regress type share on fsp (with FEs)
      share_vec <- comp_sub[[sc]]
      share_vec[is.na(share_vec)] <- 0
      fsp_vec <- comp_sub$fsp

      share_dm <- demean_multi_fe(
        cbind(share_vec, fsp_vec),
        comp_sub$county_id, comp_sub$state_year, comp_sub$time
      )

      # Simple regression of demeaned share on demeaned fsp
      fsp_dm <- share_dm[, 2]
      share_dm_dep <- share_dm[, 1]
      fsp_effect_on_share <- sum(fsp_dm * share_dm_dep) / sum(fsp_dm^2)

      # γ_k(u) ≈ δ_k(u) here (return to being type k at quantile u)
      # The composition effect for this type: Δπ_k × δ_k(u)
      if (tid_num %in% names(delta_k)) {
        composition_effect <- composition_effect +
          fsp_effect_on_share * delta_k[[tid_num]]
        cat(sprintf("    FSP effect on share of %s: %.4f\n",
                    TYPE_LABELS[tid_num], fsp_effect_on_share))
      }
    }
  }

  # --- Re-ranking residual ---
  reranking <- beta_u - delta_avg - composition_effect

  # --- Summary ---
  cat(sprintf("\n  Decomposition at u=0.05:\n"))
  cat(sprintf("    β(u)         = %.1f (total D-IV effect)\n", beta_u[1]))
  cat(sprintf("    δ(u) pop-avg = %.1f (direct, pop-weighted)\n", delta_avg[1]))
  cat(sprintf("    δ(u) density = %.1f (direct, density-weighted)\n", delta_density_weighted[1]))
  cat(sprintf("    Composition  = %.1f\n", composition_effect[1]))
  cat(sprintf("    Re-ranking   = %.1f\n", reranking[1]))

  cat(sprintf("\n  Decomposition at u=0.50:\n"))
  cat(sprintf("    β(u)         = %.1f\n", beta_u[10]))
  cat(sprintf("    δ(u) pop-avg = %.1f\n", delta_avg[10]))
  cat(sprintf("    δ(u) density = %.1f\n", delta_density_weighted[10]))
  cat(sprintf("    Composition  = %.1f\n", composition_effect[10]))
  cat(sprintf("    Re-ranking   = %.1f\n", reranking[10]))

  # Return results for plotting
  data.frame(
    quantile            = Q_GRID,
    beta_total          = beta_u,
    delta_pop_avg       = delta_avg,
    delta_density_wt    = delta_density_weighted,
    composition         = composition_effect,
    reranking           = reranking
  )
}


# =============================================================================
# Run for blacks (main result)
# =============================================================================

results_black <- run_decomposition(2, "Black Mothers")
gc()


# =============================================================================
# Decomposition figure
# =============================================================================

cat("\n=== Generating decomposition figure ===\n")

library(tidyr)

plot_df <- results_black
plot_long <- data.frame(
  quantile = rep(Q_GRID, 4),
  value = c(plot_df$beta_total, plot_df$delta_pop_avg,
            plot_df$delta_density_wt, plot_df$composition),
  component = rep(c("Total (D-IV)", "Direct (pop-avg)",
                     "Direct (density-wtd)", "Composition"),
                  each = length(Q_GRID))
)
plot_long$component <- factor(plot_long$component,
  levels = c("Total (D-IV)", "Direct (density-wtd)", "Direct (pop-avg)", "Composition"))

p <- ggplot(plot_long, aes(x = quantile, y = value, color = component, linetype = component)) +
  geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c(
    "Total (D-IV)" = "black",
    "Direct (density-wtd)" = "#c0392b",
    "Direct (pop-avg)" = "#2980b9",
    "Composition" = "#27ae60"
  )) +
  scale_linetype_manual(values = c(
    "Total (D-IV)" = "solid",
    "Direct (density-wtd)" = "solid",
    "Direct (pop-avg)" = "dashed",
    "Composition" = "dotted"
  )) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1),
                     limits = c(0.04, 0.96), expand = c(0, 0)) +
  labs(x = "Quantiles", y = "Effect on Birth Weight (grams)",
       title = "Decomposition: D-IV vs Within-Type Effect (Black Mothers)",
       color = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 11),
    legend.position = "bottom"
  )

fig_path <- file.path(FIG_DIR, "mp_decomposition_black.pdf")
ggsave(fig_path, p, width = 7, height = 5)
cat(sprintf("  Saved: %s\n", fig_path))


# --- Density weights figure ---
# Show how density weights vary across quantiles

if (exists("results_black")) {
  cat("\nDone.\n")
}
