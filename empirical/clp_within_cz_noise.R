# =============================================================================
# D-IV vs 2SLS: Effect of Within-CZ Sample Size
# =============================================================================
#
# Simulates smaller within-CZ samples by adding calibrated quantile-process
# noise to the CLP pre-computed quantile outcomes. The noise follows the
# asymptotic distribution of sample quantiles:
#
#   Cov(Q̂(u), Q̂(v)) = min(u,v)(1-max(u,v)) / (N * f(Q(u)) * f(Q(v)))
#
# where f is estimated from the quantile grid by finite differencing.
#
# Usage:
#   cd IV_dist && Rscript empirical/clp_within_cz_noise.R
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
N_WITHIN <- c(50L, 100L, 200L, 500L, 1000L, 5000L)
N_SPAGHETTI <- 30L

# --- Load data ---
cat("Loading CLP data...\n")
DATA_PATH <- file.path(PROJECT_ROOT, "data", "in",
                       "CLP_supplmentary material",
                       "CLP_empirical_application",
                       "CLP_regression_data.dta")
d <- read_dta(DATA_PATH)
sub <- d[d$class == "all", ]
n_obs <- nrow(sub)
n_q <- length(Q_GRID)
cat(sprintf("  %d observations\n", n_obs))

# --- Prepare full-sample matrices ---
Q_Yk_clean <- as.matrix(sub[, Q_LABELS])
controls_mat <- as.matrix(sub[, ALL_CONTROLS])
X <- cbind(sub$d_tradeusch_pw, controls_mat)
Z <- cbind(sub$d_tradeotch_pw_lag, controls_mat)
w <- sub$timepwt48

# --- Full-sample truth (from the clean data) ---
cat("Computing full-sample estimates (truth)...\n")
fit_2sls_truth <- estimate_2sls(Q_Yk_clean, X, Z, Q_GRID, weights = w)
fit_div_truth  <- estimate_div(Q_Yk_clean, X, Z, Q_GRID, weights = w)
b1_2sls_truth <- fit_2sls_truth$beta1[1, ]
b1_div_truth  <- fit_div_truth$beta1[1, ]


# =============================================================================
# Estimate density at each CZ x quantile (for noise calibration)
# =============================================================================

cat("Estimating within-CZ densities for noise calibration...\n")

# f_j(Q_j(u_k)) ≈ delta_u / delta_Q  (central differences, forward/backward at edges)
delta_u <- diff(Q_GRID)  # length n_q - 1, all = 0.05

f_hat <- matrix(NA_real_, n_obs, n_q)
for (j in seq_len(n_obs)) {
  q_j <- Q_Yk_clean[j, ]
  dq <- diff(q_j)  # Q_j(u_{k+1}) - Q_j(u_k), length n_q - 1

  # Central differences for interior points
  for (k in 2:(n_q - 1)) {
    denom <- q_j[k + 1] - q_j[k - 1]
    if (abs(denom) > 1e-10) {
      f_hat[j, k] <- (Q_GRID[k + 1] - Q_GRID[k - 1]) / denom
    }
  }

  # Forward/backward at edges
  if (abs(dq[1]) > 1e-10) f_hat[j, 1] <- delta_u[1] / dq[1]
  if (abs(dq[n_q - 1]) > 1e-10) f_hat[j, n_q] <- delta_u[n_q - 1] / dq[n_q - 1]
}

# Clamp density to positive minimum (avoid division by zero / extreme noise)
f_hat[is.na(f_hat)] <- 0
f_hat <- pmax(f_hat, 0.01)

# --- Precompute Cholesky factors for each observation's noise covariance ---
# Cov(u_k, u_l) = min(u_k, u_l)(1 - max(u_k, u_l)) / (f(Q(u_k)) * f(Q(u_l)))
# (without the 1/N factor, which we scale at draw time)

cat("Precomputing Cholesky factors...\n")
chol_list <- vector("list", n_obs)
for (j in seq_len(n_obs)) {
  f_j <- f_hat[j, ]
  Sigma_j <- matrix(NA_real_, n_q, n_q)
  for (k in seq_len(n_q)) {
    for (l in seq_len(n_q)) {
      Sigma_j[k, l] <- min(Q_GRID[k], Q_GRID[l]) *
                        (1 - max(Q_GRID[k], Q_GRID[l])) /
                        (f_j[k] * f_j[l])
    }
  }
  # Ensure positive definite (numerical regularization)
  Sigma_j <- Sigma_j + 1e-8 * diag(n_q)
  chol_list[[j]] <- tryCatch(chol(Sigma_j), error = function(e) {
    # Fallback: diagonal noise
    diag(sqrt(diag(Sigma_j)))
  })
}


# =============================================================================
# Exercise: IMSE vs within-CZ sample size
# =============================================================================

cat("\n=== IMSE vs Within-CZ Sample Size ===\n")

run_one_noisy_rep <- function(seed, n_within, Q_Yk_clean, X, Z, w, q_grid,
                              chol_list, b1_2sls_target, b1_div_target) {
  set.seed(seed)
  n_obs <- nrow(Q_Yk_clean)
  n_q <- ncol(Q_Yk_clean)

  # Add calibrated noise: Q_noisy = Q_clean + (1/sqrt(N)) * L * z
  Q_noisy <- Q_Yk_clean
  for (j in seq_len(n_obs)) {
    z_j <- rnorm(n_q)
    noise_j <- as.vector(crossprod(chol_list[[j]], z_j)) / sqrt(n_within)
    Q_noisy[j, ] <- Q_Yk_clean[j, ] + noise_j
  }

  fit_2sls <- tryCatch(estimate_2sls(Q_noisy, X, Z, q_grid, weights = w),
                       error = function(e) NULL)
  fit_div  <- tryCatch(estimate_div(Q_noisy, X, Z, q_grid, weights = w,
                                     return_internals = TRUE),
                       error = function(e) NULL)

  if (is.null(fit_2sls) || any(is.na(fit_2sls$beta1))) {
    return(data.frame(imse_2sls = NA, imse_div = NA,
                      frac_nonmono = NA, stringsAsFactors = FALSE))
  }

  b1_2sls <- if (is.matrix(fit_2sls$beta1)) fit_2sls$beta1[1, ] else fit_2sls$beta1
  b1_div  <- if (is.matrix(fit_div$beta1))  fit_div$beta1[1, ]  else fit_div$beta1

  # Each method measured against its own clean-data estimate (pure noise effect)
  imse_2sls <- mean((b1_2sls - b1_2sls_target)^2)
  imse_div  <- mean((b1_div  - b1_div_target)^2)

  # Fraction of fitted curves that are non-monotone
  psi_mat <- fit_div$psi_mat
  n_nonmono <- sum(apply(psi_mat, 1, function(r) any(diff(r) < -1e-10)))

  data.frame(imse_2sls = imse_2sls, imse_div = imse_div,
             frac_nonmono = n_nonmono / nrow(psi_mat),
             stringsAsFactors = FALSE)
}

n_cores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("  Using %d cores, %d reps per N_within\n", n_cores, N_REPS))

noise_results <- list()

for (n_w in N_WITHIN) {
  cat(sprintf("  N_within=%d ... ", n_w))
  t0 <- proc.time()[3]

  reps <- parallel::mclapply(seq_len(N_REPS), function(s) {
    run_one_noisy_rep(s, n_w, Q_Yk_clean, X, Z, w, Q_GRID,
                      chol_list, b1_2sls_truth, b1_div_truth)
  }, mc.cores = n_cores)

  rep_df <- do.call(rbind, reps)
  valid <- !is.na(rep_df$imse_2sls)

  noise_results[[as.character(n_w)]] <- data.frame(
    n_within       = n_w,
    imse_2sls      = mean(rep_df$imse_2sls[valid]),
    imse_div       = mean(rep_df$imse_div[valid]),
    improvement    = 1 - mean(rep_df$imse_div[valid]) / mean(rep_df$imse_2sls[valid]),
    frac_nonmono   = mean(rep_df$frac_nonmono[valid]),
    n_valid        = sum(valid),
    stringsAsFactors = FALSE
  )

  elapsed <- round(proc.time()[3] - t0, 1)
  r <- noise_results[[as.character(n_w)]]
  cat(sprintf("%.1fs | IMSE: 2SLS=%.4f, D-IV=%.4f, gain=%.1f%%, non-mono=%.0f%%\n",
              elapsed, r$imse_2sls, r$imse_div, 100 * r$improvement,
              100 * r$frac_nonmono))
}

nr_df <- do.call(rbind, noise_results)


# =============================================================================
# Figure 1: IMSE and improvement vs within-CZ N
# =============================================================================

cat("\n=== Generating figures ===\n")

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 10))

# Panel (a): IMSE vs N_within
imse_long <- rbind(
  data.frame(n_within = nr_df$n_within, imse = nr_df$imse_2sls, method = "2SLS"),
  data.frame(n_within = nr_df$n_within, imse = nr_df$imse_div,  method = "D-IV")
)

p1a <- ggplot(imse_long, aes(x = n_within, y = imse, color = method, shape = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_color_manual(name = NULL, values = c("2SLS" = "black", "D-IV" = "#c0392b")) +
  scale_shape_manual(name = NULL, values = c("2SLS" = 16, "D-IV" = 17)) +
  scale_x_log10(breaks = N_WITHIN) +
  scale_y_log10() +
  labs(x = "Within-CZ Sample Size (N)",
       y = "IMSE (log scale)",
       title = "Estimation Accuracy vs Within-CZ Sample Size") +
  theme_paper

# Panel (b): Improvement ratio and non-monotonicity
# Use log-scale x positions for even bar spacing
nr_df$log_nw <- log10(nr_df$n_within)

p1b <- ggplot(nr_df, aes(x = log_nw)) +
  geom_col(aes(y = improvement * 100), fill = "#c0392b", alpha = 0.7,
           width = 0.18) +
  geom_line(aes(y = frac_nonmono * 100), color = "grey30",
            linewidth = 0.7, linetype = "dashed") +
  geom_point(aes(y = frac_nonmono * 100), color = "grey30", size = 2) +
  scale_x_continuous(breaks = log10(N_WITHIN),
                     labels = N_WITHIN) +
  scale_y_continuous(
    name = "D-IV IMSE Improvement (%)",
    sec.axis = sec_axis(~ ., name = "Fitted Curves Non-Monotone (%)")
  ) +
  labs(x = "Within-CZ Sample Size (N)",
       title = "D-IV Improvement vs Within-CZ Precision") +
  theme_paper

fig1_path <- file.path(FIG_DIR, "clp_within_cz_imse.pdf")
pdf(fig1_path, width = 12, height = 5)
gridExtra::grid.arrange(p1a, p1b, ncol = 2)
dev.off()
cat(sprintf("  Saved: %s\n", fig1_path))


# =============================================================================
# Figure 2: Spaghetti plot at two noise levels
# =============================================================================

cat("  Generating spaghetti plots...\n")
spag_nw <- c(100L, 1000L)
spag_list <- list()

for (n_w in spag_nw) {
  for (i in seq_len(N_SPAGHETTI)) {
    set.seed(20000L + i)
    Q_noisy <- Q_Yk_clean
    for (j in seq_len(n_obs)) {
      z_j <- rnorm(n_q)
      Q_noisy[j, ] <- Q_Yk_clean[j, ] + as.vector(crossprod(chol_list[[j]], z_j)) / sqrt(n_w)
    }

    f2 <- tryCatch(estimate_2sls(Q_noisy, X, Z, Q_GRID, weights = w),
                   error = function(e) NULL)
    fd <- tryCatch(estimate_div(Q_noisy, X, Z, Q_GRID, weights = w),
                   error = function(e) NULL)

    if (!is.null(f2) && !any(is.na(f2$beta1))) {
      b2 <- if (is.matrix(f2$beta1)) f2$beta1[1, ] else f2$beta1
      bd <- if (is.matrix(fd$beta1)) fd$beta1[1, ] else fd$beta1
      spag_list[[length(spag_list) + 1]] <- data.frame(
        quantile = Q_GRID, coef = b2, method = "2SLS",
        draw = i, n_within = n_w, stringsAsFactors = FALSE)
      spag_list[[length(spag_list) + 1]] <- data.frame(
        quantile = Q_GRID, coef = bd, method = "D-IV",
        draw = i, n_within = n_w, stringsAsFactors = FALSE)
    }
  }
}

spag_df <- do.call(rbind, spag_list)
spag_df$panel <- sprintf("N = %d workers/CZ", spag_df$n_within)
spag_df$panel <- factor(spag_df$panel, levels = sprintf("N = %d workers/CZ", spag_nw))

ref_df <- rbind(
  data.frame(quantile = Q_GRID, coef = b1_2sls_truth, method = "2SLS"),
  data.frame(quantile = Q_GRID, coef = b1_div_truth,  method = "D-IV")
)

p2 <- ggplot(spag_df, aes(x = quantile, y = coef)) +
  geom_line(aes(group = interaction(draw, method), color = method),
            alpha = 0.15, linewidth = 0.4) +
  geom_line(data = ref_df, aes(color = method), linewidth = 1.1) +
  facet_grid(method ~ panel) +
  scale_color_manual(name = NULL, values = c("2SLS" = "black", "D-IV" = "#c0392b")) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.2)) +
  labs(x = "Quantile", y = "Coefficient (log points)",
       title = "Effect of Within-CZ Sample Size on Estimation Stability") +
  theme_paper +
  theme(legend.position = "none",
        strip.text = element_text(size = 10))

fig2_path <- file.path(FIG_DIR, "clp_within_cz_spaghetti.pdf")
ggsave(fig2_path, p2, width = 9, height = 6)
cat(sprintf("  Saved: %s\n", fig2_path))


# =============================================================================
# Summary table
# =============================================================================

cat("\n=== Results ===\n")
cat(sprintf("  %-10s  %-10s  %-10s  %-10s  %-10s\n",
            "N_within", "IMSE_2SLS", "IMSE_DIV", "Gain(%)", "NonMono(%)"))
for (i in seq_len(nrow(nr_df))) {
  cat(sprintf("  %-10d  %-10.4f  %-10.4f  %-10.1f  %-10.0f\n",
              nr_df$n_within[i], nr_df$imse_2sls[i], nr_df$imse_div[i],
              100 * nr_df$improvement[i], 100 * nr_df$frac_nonmono[i]))
}

cat("\nDone.\n")
