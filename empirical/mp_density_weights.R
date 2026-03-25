# =============================================================================
# Estimate Density Weights Directly (Continuous Covariates)
# =============================================================================
#
# Estimates w(x₁; u) = f_{Y|X₁}(Q(u) | x₁) / f_Y(Q(u)) using kernel
# density estimation with continuous mom_age (+ discrete sex, legit).
#
# The density-weighted average of type-specific effects approximates MP's
# conditional quantile regression estimand, while the unweighted average
# gives D-IV's unconditional estimand.
#
# Usage:
#   cd IV_dist && Rscript empirical/mp_density_weights.R
#
# =============================================================================

PROJECT_ROOT <- tryCatch({
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}, error = function(e) getwd())

library(haven)
library(ggplot2)

# --- Paths ---
DATA_DIR <- file.path(PROJECT_ROOT, "data", "out")
FIG_DIR <- file.path(PROJECT_ROOT, "empirical", "figures")

Q_GRID <- seq(0.05, 0.95, by = 0.05)


# =============================================================================
# Load data
# =============================================================================

cat("Loading enhanced density sample...\n")
d <- read_dta(file.path(DATA_DIR, "mp_density_enhanced.dta"))
cat(sprintf("  %d individuals\n", nrow(d)))

# Also load original group-level QFs for the pooled quantiles
d_orig <- read_dta(file.path(DATA_DIR, "mp_analysis_data.dta"))
Q_LABELS <- paste0("q_", as.integer(Q_GRID * 100))


# =============================================================================
# Estimate density weights for black mothers
# =============================================================================

cat("\n=== Black Mothers ===\n")

sub <- d[d$mrace == 2, ]
cat(sprintf("  Individuals: %d\n", nrow(sub)))
cat(sprintf("  Mom age range: %d-%d, median: %d\n",
            min(sub$mom_age), max(sub$mom_age), median(sub$mom_age)))

# Pooled quantiles (average group QFs across all black groups)
orig_black <- d_orig[d_orig$mrace == 2, ]
pooled_q <- colMeans(as.matrix(orig_black[, Q_LABELS]))
cat(sprintf("  Pooled Q(0.05)=%.0f, Q(0.50)=%.0f, Q(0.95)=%.0f\n",
            pooled_q[1], pooled_q[10], pooled_q[19]))

# --- Estimate conditional density f(bw | mom_age) at each pooled quantile ---
# Strategy: for each quantile u, evaluate f(Q(u) | age) for a grid of ages
# using Nadaraya-Watson kernel density estimation.
#
# f(y | age) = [Σᵢ K_h1(ageᵢ - age) K_h2(bwᵢ - y)] / [Σᵢ K_h1(ageᵢ - age)]
#
# For efficiency, bin into 1-year age groups and compute KDE of bweight
# within each bin. This is exact for a discrete kernel on age.

bw_vec <- sub$bweight
age_vec <- sub$mom_age

# Bandwidth for birth weight kernel
bw_h <- bw.nrd0(bw_vec)  # Silverman's rule
cat(sprintf("  KDE bandwidth for bweight: %.0f grams\n", bw_h))

# Age groups: 1-year bins from 14 to 45
age_bins <- 14:45
n_ages <- length(age_bins)
n_q <- length(Q_GRID)

# For each age bin, estimate f_Y|age(y) via KDE and evaluate at pooled quantiles
# density_matrix[age_idx, q_idx] = f(Q(u) | age)
density_matrix <- matrix(0, n_ages, n_q)
age_counts <- rep(0, n_ages)

cat("  Computing conditional densities by age...\n")
for (a in seq_along(age_bins)) {
  age <- age_bins[a]
  bw_age <- bw_vec[age_vec == age]
  age_counts[a] <- length(bw_age)

  if (length(bw_age) < 50) next

  # KDE of birth weight for this age group
  kde <- density(bw_age, bw = bw_h, n = 1024,
                 from = min(bw_vec) - 500, to = max(bw_vec) + 500)
  # Evaluate at pooled quantiles
  density_matrix[a, ] <- approx(kde$x, kde$y, xout = pooled_q, rule = 2)$y
}

# Marginal density f_Y(Q(u)) — pooled across all ages
kde_pooled <- density(bw_vec, bw = bw_h, n = 1024,
                      from = min(bw_vec) - 500, to = max(bw_vec) + 500)
f_marginal <- approx(kde_pooled$x, kde_pooled$y, xout = pooled_q, rule = 2)$y

# Population share of each age
pi_age <- age_counts / sum(age_counts)

# --- Density weights: w(age; u) = π(age) · f(Q(u)|age) / f(Q(u)) ---
# weight_matrix[age_idx, q_idx]
weight_matrix <- matrix(0, n_ages, n_q)
for (q in seq_len(n_q)) {
  num <- pi_age * density_matrix[, q]
  denom <- sum(num)
  if (denom > 0) weight_matrix[, q] <- num / denom
}

# --- Visualize density weights across quantiles ---
cat("\n  Density weights by age group at selected quantiles:\n")
age_groups <- list(
  "15-19" = 15:19,
  "20-24" = 20:24,
  "25-29" = 25:29,
  "30-34" = 30:34,
  "35+"   = 35:45
)

weight_by_group <- matrix(0, length(age_groups), n_q)
share_by_group <- rep(0, length(age_groups))
for (g in seq_along(age_groups)) {
  idx <- age_bins %in% age_groups[[g]]
  weight_by_group[g, ] <- colSums(weight_matrix[idx, , drop = FALSE])
  share_by_group[g] <- sum(pi_age[idx])
}

cat(sprintf("  %-8s  %-8s  %-8s  %-8s  %-8s\n",
            "Group", "π", "w(0.05)", "w(0.50)", "w(0.95)"))
for (g in seq_along(age_groups)) {
  cat(sprintf("  %-8s  %.3f     %.3f     %.3f     %.3f\n",
              names(age_groups)[g], share_by_group[g],
              weight_by_group[g, 1], weight_by_group[g, 10], weight_by_group[g, 19]))
}

# --- Key diagnostic: how much do density weights deviate from population shares? ---
# If weights ≈ shares, then density weighting doesn't matter (homogeneous densities)
# If weights differ, it tells us which types are over/under-represented at each quantile
max_deviation <- max(abs(weight_by_group - share_by_group))
cat(sprintf("\n  Max |w(u) - π|: %.3f\n", max_deviation))
cat(sprintf("  This measures how much density weighting differs from population averaging.\n"))


# =============================================================================
# Density weight figure
# =============================================================================

cat("\n=== Generating density weight figures ===\n")

# Figure 1: Density weights across quantiles by age group
plot_df <- data.frame(
  quantile = rep(Q_GRID, length(age_groups)),
  weight = as.vector(t(weight_by_group)),
  group = rep(names(age_groups), each = n_q)
)
# Add population share as horizontal reference
plot_df$pop_share <- rep(share_by_group, each = n_q)

plot_df$group <- factor(plot_df$group, levels = names(age_groups))

p1 <- ggplot(plot_df, aes(x = quantile)) +
  geom_line(aes(y = weight, color = group), linewidth = 0.8) +
  geom_line(aes(y = pop_share, color = group), linetype = "dashed", linewidth = 0.4) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.2),
                     limits = c(0.04, 0.96), expand = c(0, 0)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Quantile", y = "Density Weight w(age; u)",
       title = "Density Weights by Mother's Age (Black Mothers)",
       subtitle = "Solid = density weight at each quantile; Dashed = population share",
       color = "Age group") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "mp_density_weights_black.pdf"), p1, width = 7, height = 5)
cat("  Saved: mp_density_weights_black.pdf\n")


# Figure 2: Conditional density f(bw | age) at different ages
# Show how the birth weight distribution shifts with age
cat("  Generating conditional density plot...\n")
sample_ages <- c(16, 20, 25, 30, 35)
dens_df <- data.frame()
for (age in sample_ages) {
  bw_age <- bw_vec[age_vec == age]
  if (length(bw_age) < 50) next
  kde <- density(bw_age, bw = bw_h, n = 512, from = 500, to = 5500)
  dens_df <- rbind(dens_df, data.frame(
    bweight = kde$x, density = kde$y,
    age = factor(paste0("Age ", age))
  ))
}

p2 <- ggplot(dens_df, aes(x = bweight, y = density, color = age)) +
  geom_line(linewidth = 0.7) +
  # Add vertical lines at pooled quantiles
  geom_vline(xintercept = pooled_q[c(1, 10, 19)],
             linetype = "dotted", color = "grey50", linewidth = 0.4) +
  annotate("text", x = pooled_q[1] - 50, y = max(dens_df$density) * 0.95,
           label = "Q(0.05)", hjust = 1, size = 3) +
  annotate("text", x = pooled_q[10] + 50, y = max(dens_df$density) * 0.95,
           label = "Q(0.50)", hjust = 0, size = 3) +
  annotate("text", x = pooled_q[19] + 50, y = max(dens_df$density) * 0.95,
           label = "Q(0.95)", hjust = 0, size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Birth Weight (grams)", y = "Density",
       title = "Birth Weight Distribution by Mother's Age (Black Mothers)",
       color = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "mp_conditional_density_black.pdf"), p2, width = 7, height = 5)
cat("  Saved: mp_conditional_density_black.pdf\n")

cat("\nDone.\n")
