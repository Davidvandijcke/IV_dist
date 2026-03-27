# =============================================================================
# D-IV vs CLP: Subsampling Exercise + Fitted Distribution Plots
# =============================================================================
#
# Demonstrates D-IV advantages over unconstrained 2SLS using CLP data:
#   1. Subsampling exercise: IMSE comparison at various sample sizes
#   2. Spaghetti plot: visual instability of 2SLS vs D-IV
#   3. Fitted quantile functions for specific CZs (non-monotonicity examples)
#
# Usage:
#   cd IV_dist && Rscript empirical/clp_div_advantages.R
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

FIG_DIR <- file.path(PROJECT_ROOT, "empirical", "figures")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# --- Configuration ---
Q_GRID <- seq(0.05, 0.95, by = 0.05)
Q_LABELS <- sprintf("d_p%d_lnwkwage", as.integer(Q_GRID * 100))
CONTROLS <- c("l_shind_manuf_cbp", "l_sh_popedu_c", "l_sh_popfborn",
              "l_sh_empl_f", "l_sh_routine33", "l_task_outsource")
REG_DUMMIES <- c("reg_midatl", "reg_encen", "reg_wncen", "reg_satl",
                 "reg_escen", "reg_wscen", "reg_mount", "reg_pacif")
ALL_CONTROLS <- c(CONTROLS, REG_DUMMIES, "t2")

N_REPS <- 500L
M_SUBS <- c(50L, 75L, 100L, 150L, 200L, 350L, 500L)
N_SPAGHETTI <- 30L  # number of subsample curves to plot

# --- Load data ---
cat("Loading CLP data...\n")
DATA_PATH <- file.path(PROJECT_ROOT, "data", "in",
                       "CLP_supplmentary material",
                       "CLP_empirical_application",
                       "CLP_regression_data.dta")
d <- read_dta(DATA_PATH)
sub <- d[d$class == "all", ]

czones <- unique(sub$czone)
M_full <- length(czones)
cat(sprintf("  %d observations, %d unique CZs\n", nrow(sub), M_full))

# --- Full-sample estimates (pseudo-truth) ---
cat("Computing full-sample estimates...\n")
Q_Yk_full <- as.matrix(sub[, Q_LABELS])
controls_full <- as.matrix(sub[, ALL_CONTROLS])
X_full <- cbind(sub$d_tradeusch_pw, controls_full)
Z_full <- cbind(sub$d_tradeotch_pw_lag, controls_full)
w_full <- sub$timepwt48

fit_2sls_full <- estimate_2sls(Q_Yk_full, X_full, Z_full, Q_GRID, weights = w_full)
fit_div_full  <- estimate_div(Q_Yk_full, X_full, Z_full, Q_GRID, weights = w_full,
                               return_internals = TRUE)

b1_2sls_truth <- fit_2sls_full$beta1[1, ]
b1_div_truth  <- fit_div_full$beta1[1, ]


# =============================================================================
# Exercise 1: Subsampling IMSE comparison
# =============================================================================

cat("\n=== Exercise 1: Subsampling ===\n")

run_subsample_rep <- function(seed, m_sub, sub_data, czones, q_grid, q_labels,
                              all_controls, b1_target) {
  set.seed(seed)
  cz_sample <- sample(czones, m_sub, replace = FALSE)
  idx <- sub_data$czone %in% cz_sample
  s <- sub_data[idx, ]

  Q_Yk <- as.matrix(s[, q_labels])
  ctrls <- as.matrix(s[, all_controls])
  X <- cbind(s$d_tradeusch_pw, ctrls)
  Z <- cbind(s$d_tradeotch_pw_lag, ctrls)
  w <- s$timepwt48

  fit_2sls <- tryCatch(estimate_2sls(Q_Yk, X, Z, q_grid, weights = w),
                       error = function(e) NULL)
  fit_div  <- tryCatch(estimate_div(Q_Yk, X, Z, q_grid, weights = w,
                                     return_internals = TRUE),
                       error = function(e) NULL)

  if (is.null(fit_2sls) || any(is.na(fit_2sls$beta1))) {
    return(list(valid = FALSE))
  }

  b1_2sls <- if (is.matrix(fit_2sls$beta1)) fit_2sls$beta1[1, ] else fit_2sls$beta1
  b1_div  <- if (is.matrix(fit_div$beta1))  fit_div$beta1[1, ]  else fit_div$beta1

  # Fraction of fitted curves that are non-monotone
  psi_mat <- fit_div$psi_mat  # M x Q unconstrained fitted curves
  n_nonmono <- sum(apply(psi_mat, 1, function(r) any(diff(r) < -1e-10)))
  frac_nonmono <- n_nonmono / nrow(psi_mat)

  # Max monotonicity violation (largest backward step across all CZs)
  max_viol <- max(apply(psi_mat, 1, function(r) max(0, -min(diff(r)))))

  # Return raw coefficients for bias/variance decomposition
  list(b1_2sls = b1_2sls, b1_div = b1_div,
       frac_nonmono = frac_nonmono, max_violation = max_viol,
       valid = TRUE)
}

n_cores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("  Using %d cores, %d reps per size\n", n_cores, N_REPS))

subsample_results <- list()

for (m_sub in M_SUBS) {
  cat(sprintf("  M=%d ... ", m_sub))
  t0 <- proc.time()[3]

  reps <- parallel::mclapply(seq_len(N_REPS), function(s) {
    run_subsample_rep(s, m_sub, sub, czones, Q_GRID, Q_LABELS,
                      ALL_CONTROLS, b1_div_truth)
  }, mc.cores = n_cores)

  # Filter valid reps and stack coefficient matrices
  valid_reps <- Filter(function(r) isTRUE(r$valid), reps)
  n_valid <- length(valid_reps)
  Q <- length(Q_GRID)

  b1_2sls_mat <- matrix(NA_real_, n_valid, Q)  # R x Q
  b1_div_mat  <- matrix(NA_real_, n_valid, Q)
  frac_nm <- numeric(n_valid)
  max_viol <- numeric(n_valid)

  for (i in seq_len(n_valid)) {
    b1_2sls_mat[i, ] <- valid_reps[[i]]$b1_2sls
    b1_div_mat[i, ]  <- valid_reps[[i]]$b1_div
    frac_nm[i]       <- valid_reps[[i]]$frac_nonmono
    max_viol[i]      <- valid_reps[[i]]$max_violation
  }

  # Bias-variance decomposition: IMSE = IBias² + IVar
  # IBias²(u) = (E[β̂(u)] - β*(u))², IVar(u) = Var(β̂(u))
  # Integrated over u
  mean_2sls <- colMeans(b1_2sls_mat)
  mean_div  <- colMeans(b1_div_mat)

  ibias2_2sls <- mean((mean_2sls - b1_div_truth)^2)
  ibias2_div  <- mean((mean_div  - b1_div_truth)^2)

  ivar_2sls <- mean(apply(b1_2sls_mat, 2, var))
  ivar_div  <- mean(apply(b1_div_mat,  2, var))

  imse_2sls <- ibias2_2sls + ivar_2sls
  imse_div  <- ibias2_div  + ivar_div

  subsample_results[[as.character(m_sub)]] <- data.frame(
    m_sub          = m_sub,
    n_obs          = 2L * m_sub,
    imse_2sls      = imse_2sls,
    imse_div       = imse_div,
    ibias2_2sls    = ibias2_2sls,
    ibias2_div     = ibias2_div,
    ivar_2sls      = ivar_2sls,
    ivar_div       = ivar_div,
    improvement    = 1 - imse_div / imse_2sls,
    frac_nonmono   = mean(frac_nm),
    max_violation   = mean(max_viol),
    n_valid        = n_valid,
    stringsAsFactors = FALSE
  )

  elapsed <- round(proc.time()[3] - t0, 1)
  r <- subsample_results[[as.character(m_sub)]]
  cat(sprintf("%.1fs | IMSE: 2SLS=%.4f (b²=%.4f v=%.4f), D-IV=%.4f (b²=%.4f v=%.4f), gain=%.1f%%\n",
              elapsed, r$imse_2sls, r$ibias2_2sls, r$ivar_2sls,
              r$imse_div, r$ibias2_div, r$ivar_div,
              100 * r$improvement))
}

ss_df <- do.call(rbind, subsample_results)


# =============================================================================
# Figure 1: IMSE comparison
# =============================================================================

cat("\n=== Generating figures ===\n")

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 10))

# Panel (a): Bias-variance decomposition (stacked bars, drop M=50)
ss_plot <- ss_df[ss_df$m_sub >= 75, ]

bv_long <- rbind(
  data.frame(m_sub = ss_plot$m_sub, value = ss_plot$ibias2_2sls,
             component = "Bias²", method = "2SLS"),
  data.frame(m_sub = ss_plot$m_sub, value = ss_plot$ivar_2sls,
             component = "Variance", method = "2SLS"),
  data.frame(m_sub = ss_plot$m_sub, value = ss_plot$ibias2_div,
             component = "Bias²", method = "D-IV"),
  data.frame(m_sub = ss_plot$m_sub, value = ss_plot$ivar_div,
             component = "Variance", method = "D-IV")
)
bv_long$component <- factor(bv_long$component, levels = c("Variance", "Bias²"))
# Offset bars side-by-side
bv_long$x_pos <- bv_long$m_sub + ifelse(bv_long$method == "2SLS", -8, 8)
# Scale offsets proportionally for large M values
offsets <- c("75" = 3, "100" = 4, "150" = 6, "200" = 8, "350" = 14, "500" = 20)
bv_long$x_pos <- bv_long$m_sub +
  ifelse(bv_long$method == "2SLS", -1, 1) * offsets[as.character(bv_long$m_sub)]
bar_widths <- 2 * offsets[as.character(bv_long$m_sub)] * 0.85

p1a <- ggplot(bv_long, aes(x = x_pos, y = value, fill = interaction(component, method))) +
  geom_col(width = bar_widths, position = "stack") +
  scale_fill_manual(
    name = NULL,
    values = c("Variance.2SLS" = "grey50", "Bias².2SLS" = "grey80",
               "Variance.D-IV" = "#c0392b", "Bias².D-IV" = "#e8a0a0"),
    labels = c("Variance.2SLS" = "2SLS Variance", "Bias².2SLS" = "2SLS Bias²",
               "Variance.D-IV" = "D-IV Variance", "Bias².D-IV" = "D-IV Bias²")
  ) +
  scale_x_continuous(breaks = ss_plot$m_sub) +
  labs(x = "Number of Commuting Zones",
       y = "IMSE",
       title = "Bias-Variance Decomposition") +
  theme_paper +
  theme(legend.key.size = unit(0.4, "cm"))

# Panel (b): IMSE ratio (D-IV / 2SLS)
p1b <- ggplot(ss_plot) +
  geom_hline(yintercept = 1, color = "grey60", linetype = "dashed") +
  geom_line(aes(x = m_sub, y = 1 - improvement), color = "#c0392b",
            linewidth = 0.8) +
  geom_point(aes(x = m_sub, y = 1 - improvement), color = "#c0392b",
             size = 2.5, shape = 17) +
  annotate("text", x = max(ss_plot$m_sub), y = 1.02, label = "No improvement",
           hjust = 1, size = 3, color = "grey50") +
  scale_x_continuous(breaks = ss_plot$m_sub) +
  scale_y_continuous(limits = c(0.5, 1.05),
                     labels = function(x) sprintf("%.0f%%", x * 100)) +
  labs(x = "Number of Commuting Zones",
       y = "IMSE Ratio (D-IV / 2SLS)",
       title = "D-IV Relative Efficiency") +
  theme_paper

fig1_path <- file.path(FIG_DIR, "clp_subsampling_imse.pdf")
pdf(fig1_path, width = 12, height = 5)
gridExtra::grid.arrange(p1a, p1b, ncol = 2)
dev.off()
cat(sprintf("  Saved: %s\n", fig1_path))


# =============================================================================
# Exercise 2: Spaghetti plot (visual instability)
# =============================================================================

cat("\n=== Exercise 2: Spaghetti plot ===\n")

# Pick a small M for dramatic effect, and a moderate M for comparison
spaghetti_ms <- c(75L, 200L)

spaghetti_data <- list()
for (m_sub in spaghetti_ms) {
  cat(sprintf("  Generating %d subsample curves at M=%d...\n", N_SPAGHETTI, m_sub))
  curves_2sls <- matrix(NA, N_SPAGHETTI, length(Q_GRID))
  curves_div  <- matrix(NA, N_SPAGHETTI, length(Q_GRID))

  for (i in seq_len(N_SPAGHETTI)) {
    set.seed(10000L + i)
    cz_sample <- sample(czones, m_sub, replace = FALSE)
    idx <- sub$czone %in% cz_sample
    s <- sub[idx, ]

    Q_Yk <- as.matrix(s[, Q_LABELS])
    ctrls <- as.matrix(s[, ALL_CONTROLS])
    X <- cbind(s$d_tradeusch_pw, ctrls)
    Z <- cbind(s$d_tradeotch_pw_lag, ctrls)
    w <- s$timepwt48

    f2 <- tryCatch(estimate_2sls(Q_Yk, X, Z, Q_GRID, weights = w),
                   error = function(e) NULL)
    fd <- tryCatch(estimate_div(Q_Yk, X, Z, Q_GRID, weights = w),
                   error = function(e) NULL)

    if (!is.null(f2) && !any(is.na(f2$beta1))) {
      curves_2sls[i, ] <- if (is.matrix(f2$beta1)) f2$beta1[1, ] else f2$beta1
      curves_div[i, ]  <- if (is.matrix(fd$beta1)) fd$beta1[1, ] else fd$beta1
    }
  }

  spaghetti_data[[as.character(m_sub)]] <- list(
    curves_2sls = curves_2sls,
    curves_div  = curves_div
  )
}

# Build data frame for plotting
spag_list <- list()
for (m_sub in spaghetti_ms) {
  sd <- spaghetti_data[[as.character(m_sub)]]
  for (i in seq_len(N_SPAGHETTI)) {
    if (!any(is.na(sd$curves_2sls[i, ]))) {
      spag_list[[length(spag_list) + 1]] <- data.frame(
        quantile = Q_GRID, coef = sd$curves_2sls[i, ],
        method = "2SLS", draw = i, m_sub = m_sub, stringsAsFactors = FALSE)
      spag_list[[length(spag_list) + 1]] <- data.frame(
        quantile = Q_GRID, coef = sd$curves_div[i, ],
        method = "D-IV", draw = i, m_sub = m_sub, stringsAsFactors = FALSE)
    }
  }
}
spag_df <- do.call(rbind, spag_list)
spag_df$panel <- sprintf("M = %d CZs", spag_df$m_sub)
spag_df$panel <- factor(spag_df$panel, levels = sprintf("M = %d CZs", spaghetti_ms))

# Full-sample reference
ref_df <- rbind(
  data.frame(quantile = Q_GRID, coef = b1_2sls_truth, method = "2SLS"),
  data.frame(quantile = Q_GRID, coef = b1_div_truth,  method = "D-IV")
)

p2 <- ggplot(spag_df, aes(x = quantile, y = coef)) +
  geom_line(aes(group = interaction(draw, method), color = method),
            alpha = 0.12, linewidth = 0.4) +
  geom_line(data = ref_df, aes(color = method), linewidth = 1.1) +
  facet_grid(method ~ panel) +
  scale_color_manual(name = NULL, values = c("2SLS" = "black", "D-IV" = "#c0392b")) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.2)) +
  labs(x = "Quantile", y = "Coefficient (log points)",
       title = "Subsample Coefficient Estimates (Full-Sample Estimate in Bold)") +
  theme_paper +
  theme(legend.position = "none",
        strip.text = element_text(size = 10))

fig2_path <- file.path(FIG_DIR, "clp_spaghetti.pdf")
ggsave(fig2_path, p2, width = 9, height = 6)
cat(sprintf("  Saved: %s\n", fig2_path))


# =============================================================================
# Exercise 3: CZ-level fitted distributions (non-monotonicity examples)
# =============================================================================

cat("\n=== Exercise 3: Fitted quantile function examples ===\n")

# Compute fitted curves for all CZs using full-sample coefficients
beta0_full <- fit_div_full$beta_unc_all[1, ]  # intercept
beta1_full <- fit_div_full$beta_unc_all[-1, ]  # p x Q slope matrix

wn <- w_full / mean(w_full)
mu_X <- colSums(wn * X_full) / sum(wn)
X_centered <- sweep(X_full, 2, mu_X)

# Fitted unconstrained curves: psi_j(u) = beta0(u) + X_j' beta1(u)
psi_all <- matrix(rep(beta0_full, each = nrow(sub)), nrow = nrow(sub)) +
           X_centered %*% beta1_full  # N x Q

# Project each to monotonicity
psi_proj <- t(apply(psi_all, 1, pava))

# Monotonicity violation measure: sum of (negative differences)^2
violation <- apply(psi_all, 1, function(r) {
  diffs <- diff(r)
  sum(pmin(diffs, 0)^2)
})

# Pick top 4 most-violated CZs (use first decade for cleaner display)
decade1 <- sub$t2 == 0
viol_d1 <- violation
viol_d1[!decade1] <- 0  # zero out decade 2
top_idx <- head(order(viol_d1, decreasing = TRUE), 4)

cat(sprintf("  Top 4 CZs by violation (CZ, violation):\n"))
for (j in top_idx) {
  cat(sprintf("    CZ %d: violation=%.4f, trade_shock=%.2f\n",
              sub$czone[j], violation[j], sub$d_tradeusch_pw[j]))
}

# Build plot data
cz_list <- list()
for (j in top_idx) {
  cz_list[[length(cz_list) + 1]] <- data.frame(
    quantile = Q_GRID,
    value    = psi_all[j, ],
    type     = "2SLS (unconstrained)",
    cz       = sprintf("CZ %d", sub$czone[j]),
    stringsAsFactors = FALSE
  )
  cz_list[[length(cz_list) + 1]] <- data.frame(
    quantile = Q_GRID,
    value    = psi_proj[j, ],
    type     = "D-IV (projected)",
    cz       = sprintf("CZ %d", sub$czone[j]),
    stringsAsFactors = FALSE
  )
}
cz_df <- do.call(rbind, cz_list)

p3 <- ggplot(cz_df, aes(x = quantile, y = value, color = type, linetype = type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.3) +
  facet_wrap(~ cz, scales = "free_y", ncol = 2) +
  scale_color_manual(name = NULL,
                     values = c("2SLS (unconstrained)" = "black",
                                "D-IV (projected)" = "#c0392b")) +
  scale_linetype_manual(name = NULL,
                        values = c("2SLS (unconstrained)" = "dashed",
                                   "D-IV (projected)" = "solid")) +
  labs(x = "Quantile", y = "Fitted Quantile Function",
       title = "Fitted Wage Change Distributions: Non-Monotonicity Examples") +
  theme_paper +
  theme(strip.text = element_text(size = 10))

fig3_path <- file.path(FIG_DIR, "clp_fitted_distributions.pdf")
ggsave(fig3_path, p3, width = 9, height = 6)
cat(sprintf("  Saved: %s\n", fig3_path))


# =============================================================================
# Print summary table
# =============================================================================

cat("\n=== Subsampling Results (Bias-Variance Decomposition) ===\n")
cat(sprintf("  %-5s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-8s\n",
            "M_CZ", "IMSE_2SLS", "Bias²", "Var", "IMSE_DIV", "Bias²", "Var", "Gain(%)"))
for (i in seq_len(nrow(ss_df))) {
  cat(sprintf("  %-5d  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-8.1f\n",
              ss_df$m_sub[i],
              ss_df$imse_2sls[i], ss_df$ibias2_2sls[i], ss_df$ivar_2sls[i],
              ss_df$imse_div[i], ss_df$ibias2_div[i], ss_df$ivar_div[i],
              100 * ss_df$improvement[i]))
}

cat("\nDone.\n")
