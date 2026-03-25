# =============================================================================
# D-IV Simulation Analysis: Post-Processing and Tables
# =============================================================================
#
# Functions for loading results, computing improvements, formatting tables,
# and creating figures. Uses tidyverse for data manipulation.
#
# Usage:
#   source("simulations/analysis.R")
#   df <- load_results("iv_strength")
#   df <- compute_improvements(df)
#   create_latex_table(df, caption = "IV Strength Study", label = "tab:iv")
#   create_plots(df, "iv_strength", "simulations/results/")
#
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# --- Resolve paths ---
ANALYSIS_ROOT <- tryCatch({
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}, error = function(e) getwd())

RESULTS_DIR <- file.path(ANALYSIS_ROOT, "simulations", "results")


# =============================================================================
# Loading
# =============================================================================

#' Load saved experiment results
#'
#' @param experiment_name Character name of the experiment
#' @param results_dir     Directory containing .rds files
#' @return data.frame of results
load_results <- function(experiment_name, results_dir = RESULTS_DIR) {
  path <- file.path(results_dir, paste0(experiment_name, ".rds"))
  if (!file.exists(path)) stop("Results not found: ", path)
  readRDS(path)
}


#' Load all available results
#'
#' @param results_dir Directory containing .rds files
#' @return Named list of data.frames
load_all_results <- function(results_dir = RESULTS_DIR) {
  combined_path <- file.path(results_dir, "all_results.rds")
  if (file.exists(combined_path)) return(readRDS(combined_path))

  # Fallback: load individual files
  rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
  rds_files <- rds_files[!grepl("all_results", rds_files)]
  results <- lapply(rds_files, readRDS)
  names(results) <- tools::file_path_sans_ext(basename(rds_files))
  results
}


# =============================================================================
# Compute improvements
# =============================================================================

#' Add columns for % MSE improvement of each estimator over 2SLS
#'
#' For each (dgp, parameter combo), computes:
#'   pct_imse_improvement = 100 * (1 - imse / imse_2sls)
#'
#' @param results_df data.frame from load_results()
#' @return data.frame with added columns
compute_improvements <- function(results_df) {
  # Drop any previous improvement columns to make this idempotent
  improvement_cols <- c("imse_2sls", "ibias_2sls",
                        "pct_imse_improvement", "pct_ibias_improvement")
  results_df <- results_df %>% select(-any_of(improvement_cols))

  # Identify grouping columns: everything except estimator and metric columns
  metric_cols <- c("estimator", "imse", "ibias", "iabs_bias", "n_valid")
  group_cols  <- setdiff(names(results_df), metric_cols)

  # Get 2SLS baseline
  baseline <- results_df %>%
    filter(estimator == "2sls") %>%
    select(all_of(group_cols), imse_2sls = imse, ibias_2sls = ibias)

  # Join and compute improvement
  results_df %>%
    left_join(baseline, by = group_cols) %>%
    mutate(
      pct_imse_improvement = 100 * (1 - imse / imse_2sls),
      pct_ibias_improvement = 100 * (1 - abs(ibias) / pmax(abs(ibias_2sls), 1e-12))
    )
}


# =============================================================================
# Formatting
# =============================================================================

#' Create a clean summary table
#'
#' @param results_df   data.frame (after compute_improvements)
#' @param vary_col     Character name of the column that varies across rows
#' @param digits       Number of decimal places for metrics
#' @return data.frame formatted for display
format_table <- function(results_df, vary_col, digits = 3) {
  results_df %>%
    select(
      all_of(vary_col), estimator,
      imse, ibias, iabs_bias,
      any_of("pct_imse_improvement")
    ) %>%
    mutate(
      across(c(imse, ibias, iabs_bias), ~ round(.x, digits)),
      across(any_of("pct_imse_improvement"), ~ round(.x, 1))
    ) %>%
    arrange(across(all_of(vary_col)), estimator)
}


# =============================================================================
# LaTeX tables
# =============================================================================

#' Generate LaTeX tabular code
#'
#' Produces a self-contained tabular environment (not a full table float).
#' Each row is a parameter setting; columns are estimators with their IMSE
#' and % improvement.
#'
#' @param results_df  data.frame (after compute_improvements)
#' @param vary_col    Column that varies across rows
#' @param caption     Table caption
#' @param label       LaTeX label (e.g., "tab:iv_strength")
#' @param estimators  Which estimators to include (default: all present)
#' @return Character string of LaTeX code (also cat'd to console)
create_latex_table <- function(results_df, vary_col, caption, label,
                               estimators = NULL) {
  if (is.null(estimators)) {
    estimators <- unique(results_df$estimator)
  }

  df <- results_df %>%
    filter(estimator %in% estimators) %>%
    compute_improvements()

  # Build header
  n_est <- length(estimators)
  col_spec <- paste0("l", paste(rep("rr", n_est), collapse = ""))
  header_names <- paste(
    sapply(estimators, function(e) sprintf("\\multicolumn{2}{c}{%s}", e)),
    collapse = " & "
  )
  sub_header <- paste(
    rep("IMSE & \\%Impr.", n_est),
    collapse = " & "
  )

  # Build rows
  param_vals <- sort(unique(df[[vary_col]]))
  rows <- sapply(param_vals, function(pv) {
    cells <- sapply(estimators, function(est) {
      row_data <- df %>% filter(.data[[vary_col]] == pv, estimator == est)
      if (nrow(row_data) == 0) return("--- & ---")

      imse_str <- formatC(row_data$imse[1], format = "f", digits = 3)
      impr_str <- if (est == "2sls") {
        "---"
      } else {
        formatC(row_data$pct_imse_improvement[1], format = "f", digits = 1)
      }
      paste(imse_str, "&", impr_str)
    })
    paste(c(pv, cells), collapse = " & ")
  })

  # Assemble
  latex <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    sprintf("\\caption{%s}\n", caption),
    sprintf("\\label{%s}\n", label),
    sprintf("\\begin{tabular}{%s}\n", col_spec),
    "\\toprule\n",
    sprintf("%s & %s \\\\\n", vary_col, header_names),
    sprintf(" & %s \\\\\n", sub_header),
    "\\midrule\n",
    paste(rows, collapse = " \\\\\n"),
    " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )

  cat(latex)
  invisible(latex)
}


# =============================================================================
# Plots
# =============================================================================

#' Create standard plots for an experiment
#'
#' Produces three plots:
#'   1. IMSE comparison across estimators
#'   2. % improvement vs. the varying parameter
#'   3. Bias comparison
#'
#' @param results_df      data.frame (raw or after compute_improvements)
#' @param experiment_name Character name (used in titles and filenames)
#' @param vary_col        Column name for the x-axis variable
#' @param output_dir      Directory to save plots (NULL = don't save)
#' @return List of ggplot objects (invisibly)
create_plots <- function(results_df, experiment_name, vary_col,
                         output_dir = RESULTS_DIR) {

  df <- compute_improvements(results_df)

  # Shared theme
  theme_div <- theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 13))

  # 1. IMSE comparison
  p_imse <- ggplot(df, aes(x = .data[[vary_col]], y = imse,
                            color = estimator, shape = estimator)) +
    geom_point(size = 2.5) +
    geom_line(linewidth = 0.7) +
    scale_y_log10() +
    labs(title = paste0(experiment_name, ": IMSE comparison"),
         x = vary_col, y = "IMSE (log scale)",
         color = "Estimator", shape = "Estimator") +
    theme_div

  # 2. % improvement
  df_impr <- df %>% filter(estimator != "2sls")

  p_improvement <- ggplot(df_impr, aes(x = .data[[vary_col]],
                                        y = pct_imse_improvement,
                                        color = estimator, shape = estimator)) +
    geom_point(size = 2.5) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = paste0(experiment_name, ": % IMSE improvement over 2SLS"),
         x = vary_col, y = "% improvement",
         color = "Estimator", shape = "Estimator") +
    theme_div

  # 3. Bias comparison
  p_bias <- ggplot(df, aes(x = .data[[vary_col]], y = iabs_bias,
                            color = estimator, shape = estimator)) +
    geom_point(size = 2.5) +
    geom_line(linewidth = 0.7) +
    labs(title = paste0(experiment_name, ": integrated absolute bias"),
         x = vary_col, y = "Integrated |bias|",
         color = "Estimator", shape = "Estimator") +
    theme_div

  plots <- list(imse = p_imse, improvement = p_improvement, bias = p_bias)

  # Save
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    for (pname in names(plots)) {
      fname <- sprintf("plot_%s_%s.pdf", experiment_name, pname)
      fpath <- file.path(output_dir, fname)
      ggsave(fpath, plots[[pname]], width = 7, height = 5)
      cat(sprintf("  Saved: %s\n", fpath))
    }
  }

  invisible(plots)
}
