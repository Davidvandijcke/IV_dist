# =============================================================================
# D-IV Simulation Runner
# =============================================================================
#
# Main entry point for Monte Carlo simulations.
#
# Usage:
#   Rscript simulations/run_simulations.R                 # run all experiments
#   Rscript simulations/run_simulations.R iv_strength     # run one experiment
#   Rscript simulations/run_simulations.R iv_strength,sample_size  # run several
#
# Or source interactively:
#   source("simulations/run_simulations.R")
#   results <- run_experiment("iv_strength")
#
# =============================================================================

# --- Resolve project root (works whether sourced or Rscript'd) ---
PROJECT_ROOT <- tryCatch({
  script_dir <- dirname(sys.frame(1)$ofile)
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}, error = function(e) {
  # Fallback: assume working directory is project root
  getwd()
})

# --- Source dependencies ---
# Source utils.R and first_stage.R first.
source(file.path(PROJECT_ROOT, "R", "utils.R"))
source(file.path(PROJECT_ROOT, "R", "first_stage.R"))

# estimators.R has an internal local() block that re-sources utils.R via a
# path derived from sys.frame(1)$ofile. In nested source() calls, that frame
# can point to THIS file instead, producing a wrong path. Since we already
# loaded utils.R above, we suppress any error from the redundant re-source
# inside estimators.R by reading the file, replacing its self-sourcing block
# with a no-op, and evaluating the result.
local({
  lines <- readLines(file.path(PROJECT_ROOT, "R", "estimators.R"))
  # Blank out the local({...source(...)...}) block (lines that re-source utils.R)
  in_block <- FALSE
  brace_depth <- 0L
  for (i in seq_along(lines)) {
    if (grepl("^local\\(\\{", lines[i])) {
      in_block <- TRUE
      brace_depth <- 1L
      lines[i] <- ""
      next
    }
    if (in_block) {
      brace_depth <- brace_depth +
        nchar(gsub("[^{]", "", lines[i])) -
        nchar(gsub("[^}]", "", lines[i]))
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
# Core: single MC replication
# =============================================================================

#' Run one Monte Carlo replication
#'
#' @param seed         Integer seed for reproducibility
#' @param dgp_fn       DGP function (e.g., dgp_mp3)
#' @param dgp_args     Named list of arguments to dgp_fn (must include M, N, q_grid)
#' @param estimator_names Character vector of estimator short names (keys in ESTIMATOR_REGISTRY)
#' @param q_grid       Quantile grid
#' @return data.frame with one row per estimator: estimator, imse, ibias, iabs_bias
run_one_rep <- function(seed, dgp_fn, dgp_args, estimator_names, q_grid) {
  set.seed(seed)

  # Generate data
  data <- do.call(dgp_fn, dgp_args)

  # Compute group-level quantile functions
  Q_Yk <- compute_sample_quantiles(data$y_list, q_grid)

  true_gamma <- data$true_gamma
  Q <- length(q_grid)
  M <- length(data$y_list)

  # Run both estimators
  fit_2sls <- tryCatch(
    estimate_2sls(Q_Yk, data$x2, data$z, q_grid),
    error = function(e) NULL
  )
  fit_div <- tryCatch(
    estimate_div(Q_Yk, data$x2, data$z, q_grid),
    error = function(e) NULL
  )

  if (is.null(fit_2sls) || any(is.na(fit_2sls$beta1))) {
    return(data.frame(
      estimator = c("2sls", "div"),
      imse = NA_real_, ibias = NA_real_, iabs_bias = NA_real_,
      w2_sq = NA_real_, frac_invalid = NA_real_, avg_correction = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  # Extract slope coefficients
  b1_2sls <- if (is.matrix(fit_2sls$beta1)) fit_2sls$beta1[1, ] else fit_2sls$beta1
  b1_div  <- if (is.matrix(fit_div$beta1))  fit_div$beta1[1, ]  else fit_div$beta1

  # Coefficient IMSE
  imse_2sls <- mean((b1_2sls - true_gamma)^2)
  imse_div  <- mean((b1_div  - true_gamma)^2)

  # W₂², fraction invalid, avg correction
  # Compute fitted psi curves and their projections
  x_vec <- if (is.matrix(data$x2)) data$x2[, 1] else data$x2
  mu_x  <- mean(x_vec)

  w2_unc <- 0; w2_prj <- 0; n_invalid <- 0; total_corr <- 0

  # True group QF at grid points (if available from DGP)
  has_alpha <- !is.null(data$alpha_mat)
  has_beta0 <- !is.null(data$true_beta0)

  for (j in seq_len(M)) {
    psi_j  <- fit_2sls$beta0 + b1_2sls * (x_vec[j] - mu_x)
    proj_j <- pava(psi_j)

    # Correction magnitude
    total_corr <- total_corr + mean((proj_j - psi_j)^2)

    # Invalid group?
    if (any(diff(psi_j) < 0)) n_invalid <- n_invalid + 1

    # W₂²: need true Q_{Y_j}(u) at grid points
    if (has_alpha && has_beta0) {
      Q_true_j <- data$true_beta0 + true_gamma * (x_vec[j] - mu_x) + data$alpha_mat[j, ]
      w2_unc <- w2_unc + mean((psi_j  - Q_true_j)^2)
      w2_prj <- w2_prj + mean((proj_j - Q_true_j)^2)
    }
  }

  frac_invalid <- n_invalid / M
  avg_correction <- total_corr / M

  # W₂² (set to NA if not computable)
  if (has_alpha && has_beta0) {
    w2_2sls <- w2_unc / M
    w2_div  <- w2_prj / M
  } else {
    w2_2sls <- NA_real_
    w2_div  <- NA_real_
  }

  data.frame(
    estimator      = c("2sls", "div"),
    imse           = c(imse_2sls, imse_div),
    ibias          = c(mean(b1_2sls - true_gamma), mean(b1_div - true_gamma)),
    iabs_bias      = c(mean(abs(b1_2sls - true_gamma)), mean(abs(b1_div - true_gamma))),
    w2_sq          = c(w2_2sls, w2_div),
    frac_invalid   = c(frac_invalid, 0),  # D-IV always produces valid QFs
    avg_correction = c(avg_correction, 0),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Run one experiment
# =============================================================================

#' Run a named experiment from EXPERIMENTS
#'
#' @param experiment_name Character name matching a key in EXPERIMENTS
#' @param n_reps_override Optional integer to override the experiment's n_reps
#' @param n_cores         Number of cores for mclapply (default: detectCores() - 1)
#' @return data.frame with columns: experiment, dgp, <param cols>, estimator,
#'         imse, ibias, iabs_bias, n_valid
run_experiment <- function(experiment_name,
                           n_reps_override = NULL,
                           n_cores = max(1L, parallel::detectCores() - 1L)) {

  exp <- EXPERIMENTS[[experiment_name]]
  if (is.null(exp)) stop("Unknown experiment: ", experiment_name)

  n_reps     <- if (!is.null(n_reps_override)) n_reps_override else exp$n_reps
  q_grid     <- exp$q_grid
  estimators <- exp$estimators
  dgp_names  <- exp$dgps
  param_grid <- exp$param_grid

  cat(sprintf("\n=== Experiment: %s ===\n", experiment_name))
  cat(sprintf("  DGPs: %s\n", paste(dgp_names, collapse = ", ")))
  cat(sprintf("  Estimators: %s\n", paste(estimators, collapse = ", ")))
  cat(sprintf("  Parameter combos: %d | Reps: %d | Cores: %d\n",
              nrow(param_grid), n_reps, n_cores))

  all_results <- list()
  combo_idx <- 0L

  for (dgp_name in dgp_names) {
    dgp_fn <- match.fun(dgp_name)

    for (row_i in seq_len(nrow(param_grid))) {
      combo_idx <- combo_idx + 1L
      params <- as.list(param_grid[row_i, , drop = FALSE])

      # Build DGP argument list: always include M, N, q_grid
      dgp_args <- list(M = params$M, N = params$N, q_grid = q_grid)
      # Pass through any extra parameters the DGP accepts
      extra_params <- setdiff(names(params), c("M", "N"))
      for (ep in extra_params) {
        dgp_args[[ep]] <- params[[ep]]
      }

      param_str <- paste(
        sapply(names(params), function(nm) sprintf("%s=%s", nm, params[[nm]])),
        collapse = ", "
      )
      cat(sprintf("  [%d] %s | %s ... ", combo_idx, dgp_name, param_str))

      t0 <- proc.time()[3]

      # Run replications in parallel
      seeds <- seq_len(n_reps)
      rep_results <- parallel::mclapply(seeds, function(s) {
        run_one_rep(s, dgp_fn, dgp_args, estimators, q_grid)
      }, mc.cores = n_cores)

      # Aggregate
      rep_df <- do.call(rbind, rep_results)

      agg <- do.call(rbind, lapply(estimators, function(est) {
        sub <- rep_df[rep_df$estimator == est, ]
        valid <- !is.na(sub$imse)
        row <- data.frame(
          experiment     = experiment_name,
          dgp            = dgp_name,
          estimator      = est,
          imse           = mean(sub$imse[valid]),
          ibias          = mean(sub$ibias[valid]),
          iabs_bias      = mean(sub$iabs_bias[valid]),
          n_valid        = sum(valid),
          stringsAsFactors = FALSE
        )
        # Add W₂² and diagnostic columns if present
        if ("w2_sq" %in% names(sub)) {
          w2_valid <- valid & !is.na(sub$w2_sq)
          row$w2_sq          <- if (any(w2_valid)) mean(sub$w2_sq[w2_valid]) else NA_real_
          row$frac_invalid   <- mean(sub$frac_invalid[valid])
          row$avg_correction <- mean(sub$avg_correction[valid])
        }
        row
      }))

      # Attach parameter columns
      for (nm in names(params)) {
        agg[[nm]] <- params[[nm]]
      }

      elapsed <- round(proc.time()[3] - t0, 1)
      cat(sprintf("done (%.1fs)\n", elapsed))

      all_results[[combo_idx]] <- agg
    }
  }

  result_df <- do.call(rbind, all_results)
  rownames(result_df) <- NULL

  # Save to disk
  out_path <- file.path(RESULTS_DIR, paste0(experiment_name, ".rds"))
  saveRDS(result_df, out_path)
  cat(sprintf("  Saved: %s\n", out_path))

  result_df
}


# =============================================================================
# Run all experiments
# =============================================================================

#' Run all (or selected) experiments and save results
#'
#' @param experiments Character vector of experiment names (default: all)
#' @param n_reps_override Optional integer to override n_reps for all experiments
#' @param n_cores Number of cores for mclapply
#' @return Named list of result data.frames
run_all <- function(experiments = names(EXPERIMENTS),
                    n_reps_override = NULL,
                    n_cores = max(1L, parallel::detectCores() - 1L)) {

  cat(sprintf("D-IV Simulation Study\n"))
  cat(sprintf("Experiments to run: %s\n", paste(experiments, collapse = ", ")))
  cat(sprintf("Cores: %d\n", n_cores))

  t_start <- proc.time()[3]
  results <- list()

  for (exp_name in experiments) {
    results[[exp_name]] <- run_experiment(
      exp_name,
      n_reps_override = n_reps_override,
      n_cores = n_cores
    )
  }

  elapsed_total <- round(proc.time()[3] - t_start, 1)
  cat(sprintf("\nAll done. Total time: %.1fs\n", elapsed_total))

  # Save combined results
  combined_path <- file.path(RESULTS_DIR, "all_results.rds")
  saveRDS(results, combined_path)
  cat(sprintf("Combined results saved: %s\n", combined_path))

  results
}


# =============================================================================
# Main execution (when run as a script)
# =============================================================================

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) > 0) {
    experiments <- unlist(strsplit(args[1], ","))
    # Validate names
    bad <- setdiff(experiments, names(EXPERIMENTS))
    if (length(bad) > 0) {
      stop("Unknown experiment(s): ", paste(bad, collapse = ", "),
           "\nAvailable: ", paste(names(EXPERIMENTS), collapse = ", "))
    }
  } else {
    experiments <- names(EXPERIMENTS)
  }

  n_cores <- max(1L, parallel::detectCores() - 1L)
  cat(sprintf("Detected %d cores, using %d.\n", parallel::detectCores(), n_cores))

  run_all(experiments = experiments, n_cores = n_cores)
}
