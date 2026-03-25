# =============================================================================
# D-IV Simulation Analysis: Paper Output Pipeline
# =============================================================================
#
# Reads simulation results from simulations/results/*.rds and produces:
#   1. tabs/tab_simulations.tex  — combined simulation table
#   2. figs/fig_iv_strength.pdf  — % IMSE improvement vs IV strength
#
# Output goes to the synced Overleaf folder.
#
# Usage:
#   Rscript simulations/analysis.R
#
# =============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# --- Paths ---
PROJECT_ROOT <- tryCatch({
  normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = TRUE)
}, error = function(e) getwd())

RESULTS_DIR <- file.path(PROJECT_ROOT, "simulations", "results")

OVERLEAF <- file.path(
  Sys.getenv("HOME"),
  "University of Michigan Dropbox/David Van Dijcke/Apps/Overleaf/IV_distributions"
)
TABS_DIR <- file.path(OVERLEAF, "tabs")
FIGS_DIR <- file.path(OVERLEAF, "figs")


# =============================================================================
# Helpers
# =============================================================================

fmt <- function(x, digits = 2) {
  ifelse(abs(x) >= 100, formatC(x, format = "f", digits = 0, big.mark = ","),
  ifelse(abs(x) >= 1,   formatC(x, format = "f", digits = digits),
                         formatC(x, format = "f", digits = digits + 1)))
}

pw <- function(df, ...) {
  df %>%
    filter(estimator %in% c("2sls", "div")) %>%
    select(..., estimator, imse) %>%
    pivot_wider(names_from = estimator, values_from = imse, names_prefix = "imse_") %>%
    mutate(gain = 100 * (1 - imse_div / imse_2sls))
}


# =============================================================================
# Generate paper table
# =============================================================================

generate_table <- function() {
  res_iv  <- readRDS(file.path(RESULTS_DIR, "iv_strength.rds"))
  res_ss  <- readRDS(file.path(RESULTS_DIR, "sample_size.rds"))
  res_ht  <- readRDS(file.path(RESULTS_DIR, "heavy_tails.rds"))

  # --- Panel A: IV strength ---
  pA <- pw(res_iv, pi_Z) %>% arrange(pi_Z)

  rows_A <- sapply(seq_len(nrow(pA)), function(i) {
    r <- pA[i, ]
    sprintf("%.2f & %s & %s & %.1f",
            r$pi_Z, fmt(r$imse_2sls), fmt(r$imse_div), r$gain)
  })

  # --- Panel B: Sample size (selected rows) ---
  pB <- pw(res_ss, M, N) %>% arrange(M, N)
  pB_sel <- pB %>% filter(
    (M == 25  & N == 25)  |
    (M == 25  & N == 200) |
    (M == 50  & N == 50)  |
    (M == 100 & N == 50)  |
    (M == 200 & N == 200)
  )

  rows_B <- sapply(seq_len(nrow(pB_sel)), function(i) {
    r <- pB_sel[i, ]
    sprintf("(%d, %d) & %s & %s & %.1f",
            r$M, r$N, fmt(r$imse_2sls), fmt(r$imse_div), r$gain)
  })

  # --- Panel C: Heavy tails (selected rows) ---
  pC <- pw(res_ht, N, error_df) %>%
    mutate(error_label = ifelse(is.infinite(error_df), "None", sprintf("$t_{%g}$", error_df))) %>%
    arrange(desc(error_df), N)  # None first, then t_3, then t_2

  # Select representative rows
  pC_sel <- pC %>% filter(
    (is.infinite(error_df) & N == 10) |
    (is.infinite(error_df) & N == 50) |
    (error_df == 3 & N == 10)  |
    (error_df == 3 & N == 25)  |
    (error_df == 2 & N == 10)  |
    (error_df == 2 & N == 25)
  )

  rows_C <- sapply(seq_len(nrow(pC_sel)), function(i) {
    r <- pC_sel[i, ]
    sprintf("%s, $N = %d$ & %s & %s & %.1f",
            r$error_label, r$N, fmt(r$imse_2sls), fmt(r$imse_div), r$gain)
  })

  # --- Assemble LaTeX ---
  tex <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{Monte Carlo results: IMSE of unconstrained 2SLS and \\est. ",
    "IMSE is $\\int \\|\\hat{\\beta}_1(u) - \\beta_1(u)\\|^2\\, du$ averaged over 500 replications. ",
    "``\\% Gain'' is the percentage reduction in IMSE from the projection step. ",
    "The DGP follows \\citet{melly2025minimum} with $\\gamma(u) = \\sqrt{u}$.}\n",
    "\\label{tab:simulations}\n",
    "\\begin{tabular}{lccc}\n",
    "\\toprule\n",
    " & IMSE (2SLS) & IMSE (\\est) & \\% Gain \\\\\n",
    "\\midrule\n",
    "\\multicolumn{4}{l}{\\textit{Panel A: Instrument strength} ($n = 50$, $N = 50$)} \\\\\n",
    "\\addlinespace[2pt]\n",
    "$\\pi_Z$ & & & \\\\\n",
    paste(rows_A, collapse = " \\\\\n"), " \\\\\n",
    "\\addlinespace[6pt]\n",
    "\\multicolumn{4}{l}{\\textit{Panel B: Sample size} ($\\pi_Z = 1$)} \\\\\n",
    "\\addlinespace[2pt]\n",
    "$(n, N)$ & & & \\\\\n",
    paste(rows_B, collapse = " \\\\\n"), " \\\\\n",
    "\\addlinespace[6pt]\n",
    "\\multicolumn{4}{l}{\\textit{Panel C: Heavy-tailed within-group errors} ($n = 50$, $\\pi_Z = 1$)} \\\\\n",
    "\\addlinespace[2pt]\n",
    "Errors & & & \\\\\n",
    paste(rows_C, collapse = " \\\\\n"), " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )

  if (!dir.exists(TABS_DIR)) dir.create(TABS_DIR, recursive = TRUE)
  out_path <- file.path(TABS_DIR, "tab_simulations.tex")
  writeLines(tex, out_path)
  cat("Table written to:", out_path, "\n")

  invisible(tex)
}


# =============================================================================
# Generate paper figure
# =============================================================================

generate_figure <- function() {
  res_iv <- readRDS(file.path(RESULTS_DIR, "iv_strength.rds"))

  df <- pw(res_iv, pi_Z) %>% arrange(pi_Z)

  p <- ggplot(df, aes(x = pi_Z, y = gain)) +
    geom_line(linewidth = 0.9, color = "black") +
    geom_point(size = 2.5, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    scale_x_continuous(breaks = df$pi_Z) +
    labs(
      x = expression(paste("Instrument strength (", pi[Z], ")")),
      y = "% IMSE reduction"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 10, 5, 5)
    )

  if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR, recursive = TRUE)
  out_path <- file.path(FIGS_DIR, "fig_iv_strength.pdf")
  ggsave(out_path, p, width = 5.5, height = 3.5, device = cairo_pdf)
  cat("Figure written to:", out_path, "\n")

  invisible(p)
}


# =============================================================================
# Main
# =============================================================================

generate_paper_outputs <- function() {
  cat("Generating paper outputs...\n")
  generate_table()
  generate_figure()
  cat("Done.\n")
}

if (sys.nframe() == 0) {
  generate_paper_outputs()
}
