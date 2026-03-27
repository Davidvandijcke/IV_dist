# =============================================================================
# CLP (2016) Replication + D-IV Comparison
# =============================================================================
#
# Replicates Figures 1-3 of Chetverikov, Larsen, Palmer (2016, Econometrica)
# and applies the D-IV estimator to the same data.
#
# CLP estimated distributional effects of Chinese import competition on
# U.S. wages using quantile-by-quantile 2SLS. We replicate their results
# and compare with D-IV, which projects fitted quantile curves to
# monotonicity before recovering coefficients by OLS.
#
# Usage:
#   cd IV_dist && Rscript empirical/clp_replication.R
#
# =============================================================================

# --- Resolve project root ---
PROJECT_ROOT <- tryCatch({
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}, error = function(e) {
  getwd()
})

# --- Dependencies ---
library(haven)
library(ivreg)
library(sandwich)
library(lmtest)
library(ggplot2)

source(file.path(PROJECT_ROOT, "R", "utils.R"))
source(file.path(PROJECT_ROOT, "R", "estimators.R"))
source(file.path(PROJECT_ROOT, "R", "inference.R"))

# --- Configuration ---
DATA_PATH <- file.path(PROJECT_ROOT, "data", "in",
                       "CLP_supplmentary material",
                       "CLP_empirical_application",
                       "CLP_regression_data.dta")
FIG_DIR <- file.path(PROJECT_ROOT, "empirical", "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

Q_GRID <- seq(0.05, 0.95, by = 0.05)  # 19 quantiles
Q_LABELS <- sprintf("d_p%d_lnwkwage", as.integer(Q_GRID * 100))

CONTROLS <- c("l_shind_manuf_cbp", "l_sh_popedu_c", "l_sh_popfborn",
              "l_sh_empl_f", "l_sh_routine33", "l_task_outsource")
REG_DUMMIES <- c("reg_midatl", "reg_encen", "reg_wncen", "reg_satl",
                 "reg_escen", "reg_wscen", "reg_mount", "reg_pacif")
ALL_CONTROLS <- c(CONTROLS, REG_DUMMIES, "t2")

CLASSES <- c("all", "f", "m")
CLASS_LABELS <- c(all = "Full Sample", f = "Females", m = "Males")

# =============================================================================
# Load data
# =============================================================================

cat("Loading CLP data...\n")
d <- read_dta(DATA_PATH)
cat(sprintf("  %d observations, %d columns\n", nrow(d), ncol(d)))

# Fix trailing underscore in ADH column name
if ("adh_d_avg_lnwkwage_" %in% names(d)) {
  names(d)[names(d) == "adh_d_avg_lnwkwage_"] <- "adh_d_avg_lnwkwage"
}


# =============================================================================
# Part 1: Replicate CLP quantile-by-quantile 2SLS
# =============================================================================

cat("\n=== Part 1: Replicating CLP 2SLS ===\n")

# Build regression formula (same for all quantiles, swap LHS)
rhs_endo <- "d_tradeusch_pw"
rhs_inst <- "d_tradeotch_pw_lag"
rhs_controls <- paste(ALL_CONTROLS, collapse = " + ")
iv_formula_rhs <- sprintf("%s + %s | %s + %s",
                          rhs_endo, rhs_controls, rhs_inst, rhs_controls)

run_clp_2sls <- function(data, class_name) {
  sub <- data[data$class == class_name, ]
  n <- nrow(sub)
  cat(sprintf("  Class '%s': %d obs\n", class_name, n))

  # Quantile-by-quantile 2SLS
  qr_results <- data.frame(
    quantile = Q_GRID,
    estimate = NA_real_,
    se       = NA_real_,
    ci_lo    = NA_real_,
    ci_hi    = NA_real_
  )

  for (i in seq_along(Q_GRID)) {
    lhs <- Q_LABELS[i]
    fml <- as.formula(paste(lhs, "~", iv_formula_rhs))

    fit <- ivreg(fml, data = sub, weights = timepwt48)
    vcov_cl <- vcovCL(fit, cluster = sub$statefip)
    ct <- coeftest(fit, vcov. = vcov_cl)

    idx <- which(rownames(ct) == rhs_endo)
    qr_results$estimate[i] <- ct[idx, "Estimate"]
    qr_results$se[i]       <- ct[idx, "Std. Error"]
    qr_results$ci_lo[i]    <- ct[idx, "Estimate"] - 1.96 * ct[idx, "Std. Error"]
    qr_results$ci_hi[i]    <- ct[idx, "Estimate"] + 1.96 * ct[idx, "Std. Error"]
  }

  # ADH mean regression
  fml_mean <- as.formula(paste("adh_d_avg_lnwkwage ~", iv_formula_rhs))
  fit_mean <- ivreg(fml_mean, data = sub, weights = timepwt48)
  vcov_mean <- vcovCL(fit_mean, cluster = sub$statefip)
  ct_mean <- coeftest(fit_mean, vcov. = vcov_mean)
  idx_mean <- which(rownames(ct_mean) == rhs_endo)

  adh <- list(
    estimate = ct_mean[idx_mean, "Estimate"],
    se       = ct_mean[idx_mean, "Std. Error"],
    ci_lo    = ct_mean[idx_mean, "Estimate"] - 1.96 * ct_mean[idx_mean, "Std. Error"],
    ci_hi    = ct_mean[idx_mean, "Estimate"] + 1.96 * ct_mean[idx_mean, "Std. Error"]
  )

  list(qr = qr_results, adh = adh)
}

clp_results <- lapply(CLASSES, function(cl) run_clp_2sls(d, cl))
names(clp_results) <- CLASSES


# =============================================================================
# Part 2: Apply D-IV estimator
# =============================================================================

cat("\n=== Part 2: Applying D-IV ===\n")

run_div <- function(data, class_name) {
  sub <- data[data$class == class_name, ]
  n <- nrow(sub)

  # Build Q_Yk matrix: n x 19
  Q_Yk <- as.matrix(sub[, Q_LABELS])

  # Build X: n x (1 + n_controls) — endogenous first, then controls
  controls_mat <- as.matrix(sub[, ALL_CONTROLS])
  X <- cbind(sub$d_tradeusch_pw, controls_mat)

  # Build Z: n x (1 + n_controls) — instrument first, then same controls
  Z <- cbind(sub$d_tradeotch_pw_lag, controls_mat)

  # Weights
  w <- sub$timepwt48

  # Unconstrained 2SLS (sanity check against ivreg)
  fit_2sls <- estimate_2sls(Q_Yk, X, Z, Q_GRID, weights = w)

  # D-IV (project then OLS)
  fit_div <- estimate_div(Q_Yk, X, Z, Q_GRID, weights = w)

  # Extract the endogenous variable's coefficient (first row of beta1)
  beta1_2sls <- if (is.matrix(fit_2sls$beta1)) fit_2sls$beta1[1, ] else fit_2sls$beta1
  beta1_div  <- if (is.matrix(fit_div$beta1))  fit_div$beta1[1, ]  else fit_div$beta1

  # Sanity check: compare our 2SLS with ivreg
  clp_est <- clp_results[[class_name]]$qr$estimate
  max_diff <- max(abs(beta1_2sls - clp_est))
  cat(sprintf("  Class '%s': max |2SLS - ivreg| = %.6f\n", class_name, max_diff))

  # D-IV confidence bands
  # Endogenous variable is in column 2 of (p+1) x Q coefficient matrix
  # (row 1 = intercept, row 2 = endogenous, rows 3+ = controls)
  # Cluster by state (same as CLP)
  cl <- sub$statefip

  cat(sprintf("  Computing D-IV confidence bands (B=2000, clustered by state)...\n"))
  ci_pw <- div_pointwise_ci(Q_Yk, X, Z, Q_GRID, weights = w,
                             cluster = cl, alpha = 0.05)
  ci_ub <- div_uniform_cb(Q_Yk, X, Z, Q_GRID, weights = w,
                           cluster = cl,
                           alpha = 0.05, B = 2000L, projected = TRUE,
                           seed = 42L)

  # Row 2 = endogenous variable coefficient
  k <- 2L

  cat(sprintf("  SE ratio (proj boot / sandwich): %.3f\n",
              mean(ci_ub$se_proj[k, ] / ci_ub$se[k, ])))

  data.frame(
    quantile           = Q_GRID,
    beta1_2sls         = beta1_2sls,
    beta1_div          = beta1_div,
    div_pw_lo          = ci_pw$ci_div$lo[k, ],
    div_pw_hi          = ci_pw$ci_div$hi[k, ],
    div_proj_pw_lo     = ci_ub$ci_div_proj_pw$lo[k, ],
    div_proj_pw_hi     = ci_ub$ci_div_proj_pw$hi[k, ],
    div_ub_lo          = ci_ub$ucb_div$lo[k, ],
    div_ub_hi          = ci_ub$ucb_div$hi[k, ],
    div_ub_proj_lo     = ci_ub$ucb_div_proj$lo[k, ],
    div_ub_proj_hi     = ci_ub$ucb_div_proj$hi[k, ],
    div_se             = ci_pw$se[k, ],
    div_se_proj        = ci_ub$se_proj[k, ]
  )
}

div_results <- lapply(CLASSES, function(cl) run_div(d, cl))
names(div_results) <- CLASSES


# =============================================================================
# Part 3: Figures
# =============================================================================

cat("\n=== Part 3: Generating figures ===\n")

# Matching CLP (2016) Econometrica figure style
theme_clp <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 2),
    legend.spacing.y = unit(4, "pt"),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.4, "cm"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

for (cl in CLASSES) {
  clp <- clp_results[[cl]]
  div <- div_results[[cl]]

  plot_df <- data.frame(
    quantile          = Q_GRID,
    clp_est           = clp$qr$estimate,
    clp_lo            = clp$qr$ci_lo,
    clp_hi            = clp$qr$ci_hi,
    div_est           = div$beta1_div,
    div_pw_lo         = div$div_pw_lo,
    div_pw_hi         = div$div_pw_hi,
    div_proj_pw_lo    = div$div_proj_pw_lo,
    div_proj_pw_hi    = div$div_proj_pw_hi,
    div_ub_lo         = div$div_ub_lo,
    div_ub_hi         = div$div_ub_hi,
    div_ub_proj_lo    = div$div_ub_proj_lo,
    div_ub_proj_hi    = div$div_ub_proj_hi
  )

  # Data-driven y-axis limits (encompass all plotted elements with 5% buffer)
  all_y <- c(plot_df$clp_lo, plot_df$clp_hi,
             plot_df$div_ub_lo, plot_df$div_ub_hi,
             clp$adh$ci_lo, clp$adh$ci_hi)
  y_range <- range(all_y)
  y_pad <- 0.05 * diff(y_range)
  y_lim <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  # Round breaks to nearest 0.5
  y_breaks <- seq(floor(y_lim[1] * 2) / 2, ceiling(y_lim[2] * 2) / 2, by = 0.5)

  p <- ggplot(plot_df, aes(x = quantile)) +
    # Zero line (thin, distinct from ADH lines)
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    # Uniform bands (blue shading, mapped to fill for legend)
    geom_ribbon(aes(ymin = div_ub_lo, ymax = div_ub_hi,
                    fill = "Uniform CB"),
                alpha = 0.18) +
    geom_ribbon(aes(ymin = div_ub_proj_lo, ymax = div_ub_proj_hi,
                    fill = "Uniform CB (proj.)"),
                alpha = 0.28) +
    # Pointwise bands (red shading, mapped to fill for legend)
    geom_ribbon(aes(ymin = div_pw_lo, ymax = div_pw_hi,
                    fill = "Pointwise CI"),
                alpha = 0.20) +
    geom_ribbon(aes(ymin = div_proj_pw_lo, ymax = div_proj_pw_hi,
                    fill = "Pointwise CI (proj.)"),
                alpha = 0.35) +
    # ADH mean 95% CI (dotted lines)
    geom_hline(aes(yintercept = clp$adh$ci_lo, linetype = "ADH 95% CI"),
               color = "grey40", linewidth = 0.5) +
    geom_hline(aes(yintercept = clp$adh$ci_hi, linetype = "ADH 95% CI"),
               color = "grey40", linewidth = 0.5) +
    # ADH mean estimate (solid horizontal)
    geom_hline(aes(yintercept = clp$adh$estimate, linetype = "ADH Estimate"),
               color = "grey30", linewidth = 0.6) +
    # CLP 2SLS 95% CI (dashed)
    geom_line(aes(y = clp_lo, linetype = "CLP 95% CI"),
              color = "black", linewidth = 0.4) +
    geom_line(aes(y = clp_hi, linetype = "CLP 95% CI"),
              color = "black", linewidth = 0.4) +
    # CLP 2SLS point estimates (black connected dots)
    geom_line(aes(y = clp_est, color = "CLP"), linewidth = 0.5) +
    geom_point(aes(y = clp_est, color = "CLP", shape = "CLP"), size = 2.2) +
    # D-IV estimates (red connected triangles)
    geom_line(aes(y = div_est, color = "D-IV"), linewidth = 0.5) +
    geom_point(aes(y = div_est, color = "D-IV", shape = "D-IV"), size = 2.2) +
    # Scales
    scale_fill_manual(
      name = NULL,
      values = c("Pointwise CI"         = "#e8b4b0",
                 "Pointwise CI (proj.)" = "#c0392b",
                 "Uniform CB"           = "#b3d4e8",
                 "Uniform CB (proj.)"   = "#2980b9"),
      guide = guide_legend(order = 2, nrow = 1,
                           override.aes = list(alpha = 1),
                           keywidth = unit(0.7, "cm"))
    ) +
    scale_color_manual(
      name = NULL,
      values = c("CLP" = "black", "D-IV" = "#c0392b"),
      guide = guide_legend(order = 1, nrow = 1,
                           keywidth = unit(1, "cm"))
    ) +
    scale_shape_manual(
      name = NULL,
      values = c("CLP" = 16, "D-IV" = 17),
      guide = guide_legend(order = 1, nrow = 1,
                           keywidth = unit(1, "cm"))
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("ADH Estimate" = "solid",
                  "ADH 95% CI" = "dotted",
                  "CLP 95% CI" = "dashed"),
      guide = guide_legend(order = 1, nrow = 1,
                           keywidth = unit(1, "cm"))
    ) +
    scale_x_continuous(
      breaks = seq(0.1, 0.9, by = 0.1),
      labels = sprintf("%.1f", seq(0.1, 0.9, by = 0.1)),
      limits = c(0.04, 0.96), expand = c(0, 0)
    ) +
    scale_y_continuous(breaks = y_breaks, limits = y_lim) +
    labs(x = "Quantile", y = "Coefficient (log points)",
         title = CLASS_LABELS[cl]) +
    theme_clp

  fig_path <- file.path(FIG_DIR, sprintf("clp_div_comparison_%s.pdf", cl))
  ggsave(fig_path, p, width = 7, height = 5.5)
  cat(sprintf("  Saved: %s\n", fig_path))
}


# =============================================================================
# Summary statistics
# =============================================================================

cat("\n=== Summary ===\n")
for (cl in CLASSES) {
  clp <- clp_results[[cl]]
  div <- div_results[[cl]]

  diff <- div$beta1_div - div$beta1_2sls
  pct_changed <- mean(abs(diff) > 1e-6) * 100
  max_diff <- max(abs(diff))
  mean_diff <- mean(diff)

  cat(sprintf("\n  %s:\n", CLASS_LABELS[cl]))
  cat(sprintf("    ADH mean effect: %.3f (SE %.3f)\n",
              clp$adh$estimate, clp$adh$se))
  cat(sprintf("    Observations with PAVA adjustment: %.0f%%\n", pct_changed))
  cat(sprintf("    Max |D-IV - 2SLS|: %.4f\n", max_diff))
  cat(sprintf("    Mean (D-IV - 2SLS): %.4f\n", mean_diff))
}

cat("\nDone.\n")
