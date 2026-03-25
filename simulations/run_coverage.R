# =============================================================================
# D-IV Coverage Simulation
# =============================================================================
#
# Monte Carlo study of confidence band coverage for D-IV.
# Tests both pointwise and uniform coverage across DGP configurations.
#
# Usage:
#   Rscript simulations/run_coverage.R
#
# Or source interactively:
#   source("simulations/run_coverage.R")
#   results <- run_coverage_experiment(n_reps = 200, n_cores = 4)
#
# =============================================================================

# --- Resolve project root ---
PROJECT_ROOT <- tryCatch({
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}, error = function(e) getwd())

# --- Source dependencies ---
source(file.path(PROJECT_ROOT, "R", "utils.R"))
source(file.path(PROJECT_ROOT, "R", "first_stage.R"))

# Source estimators.R with the self-sourcing block suppressed
local({
  lines <- readLines(file.path(PROJECT_ROOT, "R", "estimators.R"))
  in_block <- FALSE
  brace_depth <- 0L
  for (i in seq_along(lines)) {
    if (grepl("^local\\(\\{", lines[i])) {
      in_block <- TRUE; brace_depth <- 1L; lines[i] <- ""; next
    }
    if (in_block) {
      brace_depth <- brace_depth +
        nchar(gsub("[^{]", "", lines[i])) - nchar(gsub("[^}]", "", lines[i]))
      lines[i] <- ""
      if (brace_depth <= 0) in_block <- FALSE
    }
  }
  eval(parse(text = lines), envir = globalenv())
})

# Source inference.R with the same treatment
local({
  lines <- readLines(file.path(PROJECT_ROOT, "R", "inference.R"))
  in_block <- FALSE
  brace_depth <- 0L
  for (i in seq_along(lines)) {
    if (grepl("^local\\(\\{", lines[i])) {
      in_block <- TRUE; brace_depth <- 1L; lines[i] <- ""; next
    }
    if (in_block) {
      brace_depth <- brace_depth +
        nchar(gsub("[^{]", "", lines[i])) - nchar(gsub("[^}]", "", lines[i]))
      lines[i] <- ""
      if (brace_depth <= 0) in_block <- FALSE
    }
  }
  eval(parse(text = lines), envir = globalenv())
})

source(file.path(PROJECT_ROOT, "simulations", "dgp.R"))
source(file.path(PROJECT_ROOT, "simulations", "config.R"))

RESULTS_DIR <- file.path(PROJECT_ROOT, "simulations", "results")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)


# =============================================================================
# Single MC replication with coverage checking
# =============================================================================

#' Run one MC replication and check confidence band coverage
#'
#' @param seed       Random seed
#' @param dgp_args   Named list of DGP arguments
#' @param q_grid     Quantile grid
#' @param alpha      Significance level
#' @param B          Bootstrap replications for uniform bands
#' @param projected  Whether to compute projected bootstrap
#' @return data.frame with coverage results
run_one_rep_coverage <- function(seed, dgp_args, q_grid,
                                 alpha = 0.05, B = 500L,
                                 projected = TRUE) {
  set.seed(seed)

  # Generate data
  data <- do.call(dgp_mp, dgp_args)
  Q_Yk <- compute_sample_quantiles(data$y_list, q_grid)

  true_gamma <- data$true_gamma  # Q-vector: true slope coefficient gamma(u)
  Q <- length(q_grid)

  # Run inference (tryCatch in case of singular matrices)
  ci_pw <- tryCatch(
    div_pointwise_ci(Q_Yk, data$x2, data$z, q_grid, alpha = alpha),
    error = function(e) NULL
  )
  ci_ub <- tryCatch(
    div_uniform_cb(Q_Yk, data$x2, data$z, q_grid,
                   alpha = alpha, B = B, projected = projected),
    error = function(e) NULL
  )

  if (is.null(ci_pw) || is.null(ci_ub)) {
    return(data.frame(
      estimator         = c("2sls", "div", "div_proj"),
      pw_coverage_rate  = NA_real_,
      ub_covers         = NA_integer_,
      avg_pw_width      = NA_real_,
      avg_ub_width      = NA_real_,
      stringsAsFactors  = FALSE
    ))
  }

  # Slope coefficient is row 2 of the (p+1) x Q matrices
  k <- 2L

  # --- Pointwise coverage (fraction of quantiles where true gamma is in CI) ---
  pw_in_unc <- (true_gamma >= ci_pw$ci_unc$lo[k, ]) &
               (true_gamma <= ci_pw$ci_unc$hi[k, ])
  pw_in_div <- (true_gamma >= ci_pw$ci_div$lo[k, ]) &
               (true_gamma <= ci_pw$ci_div$hi[k, ])

  # --- Uniform coverage (does entire true curve lie within band?) ---
  ub_covers_unc <- all(true_gamma >= ci_ub$ucb_unc$lo[k, ] &
                       true_gamma <= ci_ub$ucb_unc$hi[k, ])
  ub_covers_div <- all(true_gamma >= ci_ub$ucb_div$lo[k, ] &
                       true_gamma <= ci_ub$ucb_div$hi[k, ])

  # --- Band widths ---
  pw_width_unc <- mean(ci_pw$ci_unc$hi[k, ] - ci_pw$ci_unc$lo[k, ])
  pw_width_div <- mean(ci_pw$ci_div$hi[k, ] - ci_pw$ci_div$lo[k, ])
  ub_width_unc <- mean(ci_ub$ucb_unc$hi[k, ] - ci_ub$ucb_unc$lo[k, ])
  ub_width_div <- mean(ci_ub$ucb_div$hi[k, ] - ci_ub$ucb_div$lo[k, ])

  results <- data.frame(
    estimator        = c("2sls", "div"),
    pw_coverage_rate = c(mean(pw_in_unc), mean(pw_in_div)),
    ub_covers        = c(as.integer(ub_covers_unc), as.integer(ub_covers_div)),
    avg_pw_width     = c(pw_width_unc, pw_width_div),
    avg_ub_width     = c(ub_width_unc, ub_width_div),
    stringsAsFactors = FALSE
  )

  # Add projected bootstrap results if available
  if (projected && !is.null(ci_ub$ucb_div_proj)) {
    ub_covers_proj <- all(true_gamma >= ci_ub$ucb_div_proj$lo[k, ] &
                          true_gamma <= ci_ub$ucb_div_proj$hi[k, ])
    ub_width_proj <- mean(ci_ub$ucb_div_proj$hi[k, ] - ci_ub$ucb_div_proj$lo[k, ])

    results <- rbind(results, data.frame(
      estimator        = "div_proj",
      pw_coverage_rate = mean(pw_in_div),  # same pointwise CIs
      ub_covers        = as.integer(ub_covers_proj),
      avg_pw_width     = pw_width_div,
      avg_ub_width     = ub_width_proj,
      stringsAsFactors = FALSE
    ))
  }

  results
}


# =============================================================================
# Coverage experiment runner
# =============================================================================

#' Run the full coverage simulation study
#'
#' @param n_reps   Number of MC replications
#' @param alpha    Significance level
#' @param B        Bootstrap replications per MC draw
#' @param n_cores  Cores for mclapply
#' @return data.frame with aggregated coverage results
run_coverage_experiment <- function(n_reps = 500L, alpha = 0.05, B = 500L,
                                    n_cores = max(1L, parallel::detectCores() - 1L)) {

  q_grid <- Q_GRID

  # Parameter configurations to test
  param_grid <- data.frame(
    M    = c(50,  50,  100, 200, 50),
    N    = c(50,  50,  50,  50,  50),
    pi_Z = c(0.3, 0.5, 0.5, 0.5, 0.7),
    stringsAsFactors = FALSE
  )

  cat(sprintf("D-IV Coverage Simulation\n"))
  cat(sprintf("  Reps: %d | Bootstrap: %d | Alpha: %.2f | Cores: %d\n",
              n_reps, B, alpha, n_cores))

  all_results <- list()

  for (row_i in seq_len(nrow(param_grid))) {
    params <- param_grid[row_i, ]
    dgp_args <- list(
      M = params$M, N = params$N, q_grid = q_grid,
      endogenous = TRUE, heterogeneity = TRUE,
      pi_Z = params$pi_Z, beta_slope = 0
    )

    param_str <- sprintf("M=%d, N=%d, pi_Z=%.1f", params$M, params$N, params$pi_Z)
    cat(sprintf("\n  [%d] %s ... ", row_i, param_str))
    t0 <- proc.time()[3]

    # Run replications in parallel
    seeds <- seq_len(n_reps)
    rep_results <- parallel::mclapply(seeds, function(s) {
      run_one_rep_coverage(s, dgp_args, q_grid, alpha = alpha, B = B)
    }, mc.cores = n_cores)

    rep_df <- do.call(rbind, rep_results)

    # Aggregate by estimator
    agg <- do.call(rbind, lapply(unique(rep_df$estimator), function(est) {
      sub <- rep_df[rep_df$estimator == est, ]
      valid <- !is.na(sub$pw_coverage_rate)
      data.frame(
        M                = params$M,
        N                = params$N,
        pi_Z             = params$pi_Z,
        estimator        = est,
        pw_coverage_mean = mean(sub$pw_coverage_rate[valid]),
        ub_coverage      = mean(sub$ub_covers[valid]),
        avg_pw_width     = mean(sub$avg_pw_width[valid]),
        avg_ub_width     = mean(sub$avg_ub_width[valid]),
        n_valid          = sum(valid),
        stringsAsFactors = FALSE
      )
    }))

    elapsed <- round(proc.time()[3] - t0, 1)
    cat(sprintf("done (%.1fs)\n", elapsed))

    # Print results immediately
    for (i in seq_len(nrow(agg))) {
      cat(sprintf("    %-10s  PW cov: %.3f  UB cov: %.3f  PW width: %.3f  UB width: %.3f\n",
                  agg$estimator[i], agg$pw_coverage_mean[i], agg$ub_coverage[i],
                  agg$avg_pw_width[i], agg$avg_ub_width[i]))
    }

    all_results[[row_i]] <- agg
  }

  result_df <- do.call(rbind, all_results)
  rownames(result_df) <- NULL

  # Save
  out_path <- file.path(RESULTS_DIR, "coverage.rds")
  saveRDS(result_df, out_path)
  cat(sprintf("\nResults saved: %s\n", out_path))

  result_df
}


# =============================================================================
# Main execution
# =============================================================================

if (sys.nframe() == 0) {
  n_cores <- max(1L, parallel::detectCores() - 1L)
  cat(sprintf("Detected %d cores, using %d.\n", parallel::detectCores(), n_cores))
  run_coverage_experiment(n_reps = 500L, B = 500L, n_cores = n_cores)
}
