# =============================================================================
# D-IV Simulation Analysis: Paper Output Pipeline
# =============================================================================
#
# Reads simulation results and produces LaTeX tables + figures for Overleaf.
#
# Usage:  Rscript simulations/analysis.R
# =============================================================================

suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})

PROJECT_ROOT <- tryCatch({
  normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = TRUE)
}, error = function(e) getwd())

RESULTS_DIR <- file.path(PROJECT_ROOT, "simulations", "results")
OVERLEAF <- file.path(Sys.getenv("HOME"),
  "University of Michigan Dropbox/David Van Dijcke/Apps/Overleaf/IV_distributions")
TABS_DIR <- file.path(OVERLEAF, "tabs")
FIGS_DIR <- file.path(OVERLEAF, "figs")

fmt <- function(x, d=2) {
  ifelse(abs(x) >= 100, formatC(x, format="f", digits=0, big.mark=","),
  ifelse(abs(x) >= 1,   formatC(x, format="f", digits=d),
                         formatC(x, format="f", digits=d+1)))
}

# =============================================================================
# Table: Main simulation results
# =============================================================================

generate_table <- function() {
  r_iv   <- readRDS(file.path(RESULTS_DIR, "iv_strength.rds"))
  r_ctrl <- readRDS(file.path(RESULTS_DIR, "controls.rds"))
  r_real <- readRDS(file.path(RESULTS_DIR, "realistic.rds"))
  r_mp   <- readRDS(file.path(RESULTS_DIR, "mp_baseline.rds"))

  pw <- function(df, ...) {
    cols <- enquos(...)
    df %>% filter(estimator %in% c("2sls","div")) %>%
      select(!!!cols, estimator, imse, w2_sq, frac_invalid) %>%
      pivot_wider(names_from=estimator,
                  values_from=c(imse, w2_sq, frac_invalid),
                  names_glue="{.value}_{estimator}") %>%
      mutate(b1_gain = 100*(1 - imse_div/imse_2sls),
             w2_gain = ifelse(!is.na(w2_sq_2sls) & w2_sq_2sls > 0,
                              100*(1 - w2_sq_div/w2_sq_2sls), NA))
  }

  # --- Panel A: IV strength ---
  # Map pi_Z to approximate median F-stat (precomputed)
  f_map <- c("0.1"=1, "0.2"=3, "0.3"=6, "0.5"=17, "0.7"=44, "1"=563)
  pA <- pw(r_iv, pi_Z) %>% arrange(pi_Z) %>%
    mutate(F_approx = f_map[as.character(pi_Z)])
  rows_A <- sapply(seq_len(nrow(pA)), function(i) {
    r <- pA[i,]
    sprintf("%.0f & %s & %s & %.1f & %s & %s & %.1f & %.1f",
      r$F_approx, fmt(r$imse_2sls), fmt(r$imse_div), r$b1_gain,
      fmt(r$w2_sq_2sls), fmt(r$w2_sq_div), r$w2_gain,
      100*r$frac_invalid_2sls)
  })

  # --- Panel B: Controls + hetero FS (selected rows) ---
  pB <- pw(r_ctrl, p, hetero_fs) %>% arrange(p, hetero_fs)
  pB_sel <- pB %>% filter(
    (p==1 & hetero_fs==0) | (p==2 & hetero_fs==0) | (p==5 & hetero_fs==0) |
    (p==1 & hetero_fs==0.5) | (p==2 & hetero_fs==0.5) | (p==5 & hetero_fs==0.5) |
    (p==2 & hetero_fs==1) | (p==5 & hetero_fs==1)
  )
  rows_B <- sapply(seq_len(nrow(pB_sel)), function(i) {
    r <- pB_sel[i,]
    sprintf("$p=%d$, $\\delta=%.1f$ & %s & %s & %.1f & %s & %s & %.1f & %.1f",
      r$p, r$hetero_fs, fmt(r$imse_2sls), fmt(r$imse_div), r$b1_gain,
      fmt(r$w2_sq_2sls), fmt(r$w2_sq_div), r$w2_gain,
      100*r$frac_invalid_2sls)
  })

  # --- Panel C: Realistic ---
  pC <- pw(r_real, p)
  rows_C <- sprintf("$p=3$, $\\delta=0.5$, $t_5$, $\\beta_{\\text{slope}}=0.2$ & %s & %s & %.1f & %s & %s & %.1f & %.1f",
    fmt(pC$imse_2sls), fmt(pC$imse_div), pC$b1_gain,
    fmt(pC$w2_sq_2sls), fmt(pC$w2_sq_div), pC$w2_gain,
    100*pC$frac_invalid_2sls)

  # --- Panel D: MP baseline ---
  pD <- r_mp %>% filter(estimator %in% c("2sls","div")) %>%
    select(M, N, estimator, imse, frac_invalid) %>%
    pivot_wider(names_from=estimator, values_from=c(imse, frac_invalid),
                names_glue="{.value}_{estimator}") %>%
    mutate(b1_gain = 100*(1-imse_div/imse_2sls)) %>%
    filter((M==25 & N==25) | (M==25 & N==50) | (M==50 & N==50))
  rows_D <- sapply(seq_len(nrow(pD)), function(i) {
    r <- pD[i,]
    sprintf("$(n,N)=(%d,%d)$ & %s & %s & %.1f & --- & --- & --- & %.1f",
      r$M, r$N, fmt(r$imse_2sls), fmt(r$imse_div), r$b1_gain,
      100*r$frac_invalid_2sls)
  })

  tex <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\\small\n",
    "\\caption{Monte Carlo results (500 replications). ",
    "IMSE is $\\int \\|\\hat{\\beta}_1(u) - \\beta_1(u)\\|^2\\, du$ averaged over replications. ",
    "$W_2^2$ is the average squared Wasserstein distance between estimated and true conditional quantile functions. ",
    "``Inv.'' is the fraction of groups with non-monotone fitted $\\hat{\\psi}_{X_j}$.}\n",
    "\\label{tab:simulations}\n",
    "\\begin{tabular}{lccccccr}\n",
    "\\toprule\n",
    " & \\multicolumn{3}{c}{Coefficient IMSE} & \\multicolumn{3}{c}{$W_2^2$} & \\\\\n",
    "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}\n",
    " & 2SLS & \\est & \\% Gain & 2SLS & \\est & \\% Gain & Inv.~(\\%) \\\\\n",
    "\\midrule\n",
    "\\multicolumn{8}{l}{\\textit{Panel A: Instrument strength} ($n=50$, $N=50$, centered $x_2$)} \\\\\n",
    "\\addlinespace[2pt]\n",
    "$F$ & & & & & & & \\\\\n",
    paste(rows_A, collapse=" \\\\\n"), " \\\\\n",
    "\\addlinespace[4pt]\n",
    "\\multicolumn{8}{l}{\\textit{Panel B: Controls \\& heterogeneous first stage} ($n=50$, $N=50$, $F \\approx 13$--$17$)} \\\\\n",
    "\\addlinespace[2pt]\n",
    paste(rows_B, collapse=" \\\\\n"), " \\\\\n",
    "\\addlinespace[4pt]\n",
    "\\multicolumn{8}{l}{\\textit{Panel C: Realistic combination} ($n=50$, $N=50$, $F \\approx 15$)} \\\\\n",
    "\\addlinespace[2pt]\n",
    paste(rows_C, collapse=" \\\\\n"), " \\\\\n",
    "\\addlinespace[4pt]\n",
    "\\multicolumn{8}{l}{\\textit{Panel D: \\citet{melly2025minimum} baseline} ($\\pi_Z=1$, log-normal $x_2$)} \\\\\n",
    "\\addlinespace[2pt]\n",
    paste(rows_D, collapse=" \\\\\n"), " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )

  if (!dir.exists(TABS_DIR)) dir.create(TABS_DIR, recursive=TRUE)
  out <- file.path(TABS_DIR, "tab_simulations.tex")
  writeLines(tex, out)
  cat("Table:", out, "\n")
  invisible(tex)
}

# =============================================================================
# Figure: IMSE improvement vs IV strength
# =============================================================================

generate_figure <- function() {
  r_iv <- readRDS(file.path(RESULTS_DIR, "iv_strength.rds"))

  f_map <- c("0.1"=1, "0.2"=3, "0.3"=6, "0.5"=17, "0.7"=44, "1"=563)

  df <- r_iv %>%
    filter(estimator %in% c("2sls","div")) %>%
    select(pi_Z, estimator, imse, w2_sq) %>%
    pivot_wider(names_from=estimator, values_from=c(imse, w2_sq), names_glue="{.value}_{estimator}") %>%
    mutate(F_stat = f_map[as.character(pi_Z)],
           b1_gain = 100*(1-imse_div/imse_2sls),
           w2_gain = 100*(1-w2_sq_div/w2_sq_2sls)) %>%
    filter(F_stat <= 50) %>%  # drop F=563 for readability
    select(F_stat, b1_gain, w2_gain) %>%
    pivot_longer(c(b1_gain, w2_gain), names_to="metric", values_to="gain") %>%
    mutate(metric = factor(metric, levels=c("b1_gain","w2_gain"),
                           labels=c("Coefficient IMSE", expression(W[2]^2))))

  p <- ggplot(df, aes(x=F_stat, y=gain, linetype=metric, shape=metric)) +
    geom_line(linewidth=0.8, color="black") +
    geom_point(size=2.5, color="black") +
    geom_hline(yintercept=0, linetype="dashed", color="gray60") +
    scale_x_continuous(breaks=unique(df$F_stat), trans="log10",
                       labels=function(x) round(x)) +
    scale_linetype_manual(values=c("solid","dashed")) +
    scale_shape_manual(values=c(16, 17)) +
    labs(x="First-stage F-statistic",
         y="% improvement over 2SLS",
         linetype="Metric", shape="Metric") +
    theme_minimal(base_size=11) +
    theme(panel.grid.minor=element_blank(),
          legend.position="bottom",
          plot.margin=margin(5,10,5,5))

  if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR, recursive=TRUE)
  out <- file.path(FIGS_DIR, "fig_iv_strength.pdf")
  ggsave(out, p, width=5.5, height=4, device=cairo_pdf)
  cat("Figure:", out, "\n")
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

if (sys.nframe() == 0) generate_paper_outputs()
