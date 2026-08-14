# ============================================================================
# MONTE CARLO SIMULATION
#   Differential Censoring Is Not Informative Censoring: A Simulation Study
#   of Censoring Bias in Survival Analyses -- manuscript programming code share.
#
#   Sources simulation_ipcw_github.R for all data-generating and estimation
#   logic, runs R Monte Carlo replicates across the 4 censoring scenarios
#   (S1-S4), and summarizes bias / coverage / MSE for the naive vs. IPCW
#   estimators against the counterfactual truth.
#
#   Run this script with the repository root (where both .R files live) as
#   the working directory, so that "simulation_ipcw_github.R" and
#   "results/" resolve correctly.
# ============================================================================

library(tidyverse)
library(survival)
library(riskRegression)
library(survRM2)
library(data.table)
library(parallel)

output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

source("simulation_ipcw_github.R")

# ============================================================================
# USER-CONTROLLED SETTINGS
# ============================================================================

R         <- 200       # number of Monte Carlo replicates (manuscript used 1000)
tmax      <- 24
max_fu    <- 24
base_seed <- 20260606
n_boot    <- 200       # bootstrap iterations per replicate (0 = skip; manuscript used 500)

n_s1 <- 2000
n_s2 <- 2000
n_s3 <- 5000
n_s4 <- 5000

n_cores <- max(1L, parallel::detectCores() - 1L)   # lower this if RAM-constrained

# ============================================================================
# ONE REPLICATE
# ============================================================================

run_one_replicate <- function(r) {
  set.seed(base_seed + 1000L + r)

  dat_s1 <- simulate_scenario1(n = n_s1, dt = 1, max_fu = max_fu)
  dat_s2 <- simulate_scenario2(n = n_s2, dt = 1, max_fu = max_fu)
  dat_s3 <- simulate_scenario3(n = n_s3, dt = 1, max_fu = max_fu)
  dat_s4 <- simulate_scenario4(n = n_s4, dt = 1, max_fu = max_fu)

  res_s1 <- analyze_coxph_pt(dat_s1$dat_pt, dat_s1$dat, "S1_noninf_nondiff",     tmax, n_boot = n_boot)
  res_s2 <- analyze_coxph_pt(dat_s2$dat_pt, dat_s2$dat, "S2_noninf_diff_nausea", tmax, n_boot = n_boot)
  res_s3 <- analyze_coxph_pt(dat_s3$dat_pt, dat_s3$dat, "S3_inf_nondiff",        tmax, n_boot = n_boot)
  res_s4 <- analyze_coxph_pt(dat_s4$dat_pt, dat_s4$dat, "S4_inf_diff",           tmax, n_boot = n_boot)

  bind_rows(res_s1$numeric, res_s2$numeric, res_s3$numeric, res_s4$numeric) %>%
    mutate(replicate = r)
}

# A single bad replicate should not kill the whole run.
run_replicate_safe <- function(r) {
  tryCatch(run_one_replicate(r), error = function(e) {
    message("Replicate ", r, " failed: ", conditionMessage(e)); NULL
  })
}

# ============================================================================
# RUN ALL REPLICATES (parallel across replicates)
# ============================================================================
# NOTE: on Unix, bootstrap_ipcw_rd_rmstd() also parallelises internally via
# mclapply when n_boot is large. Combined with replicate-level parallelism
# here, that oversubscribes cores (n_cores x internal bootstrap cores) rather
# than causing a correctness problem -- reduce n_cores if runs feel slower
# than expected on a Unix machine.

cat(sprintf("Running %d Monte Carlo replicates on %d core(s)...\n", R, n_cores))

if (n_cores > 1L) {
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterSetRNGStream(cl, base_seed)
  parallel::clusterExport(cl, c("n_s1", "n_s2", "n_s3", "n_s4", "tmax", "max_fu",
                                "n_boot", "base_seed", "run_one_replicate"))
  parallel::clusterEvalQ(cl, source("simulation_ipcw_github.R"))
  mc_results_list <- tryCatch(
    parallel::parLapply(cl, seq_len(R), run_replicate_safe),
    finally = parallel::stopCluster(cl)
  )
} else {
  mc_results_list <- lapply(seq_len(R), run_replicate_safe)
}

mc_results_raw <- bind_rows(mc_results_list)

n_ok <- length(unique(mc_results_raw$replicate))
if (n_ok < R) {
  warning(sprintf("Only %d of %d replicates completed successfully.", n_ok, R))
}
if (n_ok == 0L) stop("No replicates completed -- aborting.")

saveRDS(mc_results_raw, file.path(output_dir, "mc_results_raw.rds"))
write_csv(mc_results_raw, file.path(output_dir, "mc_results_raw.csv"))
cat(sprintf("Raw results saved (%d rows, %d replicates).\n", nrow(mc_results_raw), n_ok))

# ============================================================================
# PERFORMANCE METRIC FUNCTION
# ============================================================================

# Computes bias, MSE, empirical SE, empirical CI of estimates, and (optionally)
# coverage + mean width of the model/bootstrap CIs supplied via lower_ci/upper_ci.
# Coverage/MSE are evaluated against a single fixed scalar benchmark
# (truth_mean = mean of the per-replicate counterfactual truths) rather than
# each replicate's own correlated truth, which would otherwise understate
# variance and mask bias via spurious over-coverage.
calculate_performance <- function(estimates, lower_ci, upper_ci, truth) {
  valid   <- !is.na(estimates) & !is.na(truth)
  est_v   <- estimates[valid]
  lo_v    <- lower_ci[valid]
  hi_v    <- upper_ci[valid]
  truth_v <- truth[valid]
  n_valid <- length(est_v)

  if (n_valid == 0L) {
    return(tibble(
      n_converged   = 0L,
      mean          = NA_real_, empirical_se  = NA_real_, mse        = NA_real_,
      bias          = NA_real_, rel_bias_pct  = NA_real_,
      empirical_lo  = NA_real_, empirical_hi  = NA_real_,
      coverage      = NA_real_, mean_ci_width = NA_real_
    ))
  }

  truth_mean <- mean(truth_v)
  mean_est   <- mean(est_v)
  bias       <- mean_est - truth_mean
  rel_bias   <- if (truth_mean != 0) 100 * bias / truth_mean else NA_real_

  ci_valid  <- !is.na(lo_v) & !is.na(hi_v)
  coverage  <- if (any(ci_valid)) {
    mean(lo_v[ci_valid] <= truth_mean & hi_v[ci_valid] >= truth_mean)
  } else NA_real_
  ci_width  <- if (any(ci_valid)) mean(hi_v[ci_valid] - lo_v[ci_valid]) else NA_real_

  tibble(
    n_converged   = n_valid,
    mean          = mean_est,
    empirical_se  = sd(est_v),
    mse           = mean((est_v - truth_mean)^2),
    bias          = bias,
    rel_bias_pct  = rel_bias,
    empirical_lo  = as.numeric(quantile(est_v, 0.025)),
    empirical_hi  = as.numeric(quantile(est_v, 0.975)),
    coverage      = coverage,
    mean_ci_width = ci_width
  )
}

# ============================================================================
# METRIC SPECIFICATIONS
# ============================================================================
# Each row defines one performance calculation:
#   est_col   — column with the point estimate
#   lo_col    — column with CI lower bound (NA = no CI)
#   hi_col    — column with CI upper bound (NA = no CI)
#   truth_col — column with the truth value
#   ci_type   — "analytical", "bootstrap", or "none"

metric_specs <- tribble(
  ~metric,              ~est_col,          ~lo_col,              ~hi_col,              ~truth_col,    ~ci_type,
  "HR_naive",           "HR_naive_est",    "HR_naive_lo",        "HR_naive_hi",        "HR_cf",       "analytical",
  "HR_ipcw",            "HR_ipcw_est",     "HR_ipcw_lo",         "HR_ipcw_hi",         "HR_cf",       "analytical",
  "RD_naive",           "RD_naive_est",    "RD_naive_lo",        "RD_naive_hi",        "RD_cf",       "analytical",
  "RD_ipcw_boot",       "RD_ipcw_est",     "RD_ipcw_boot_lo",    "RD_ipcw_boot_hi",    "RD_cf",       "bootstrap",
  "RMSTD_naive",        "RMSTD_naive_est", "RMSTD_naive_lo",     "RMSTD_naive_hi",     "RMSTD_cf",    "analytical",
  "RMSTD_ipcw_boot",    "RMSTD_ipcw_est",  "RMSTD_ipcw_boot_lo", "RMSTD_ipcw_boot_hi", "RMSTD_cf",    "bootstrap",
  "Risk_A0_naive",      "Risk_A0_naive",   "Risk_A0_naive_lo",   "Risk_A0_naive_hi",   "Risk_A0_cf",  "analytical",
  "Risk_A1_naive",      "Risk_A1_naive",   "Risk_A1_naive_lo",   "Risk_A1_naive_hi",   "Risk_A1_cf",  "analytical",
  "Risk_A0_ipcw",       "Risk_A0_ipcw",    "Risk_A0_ipcw_lo",    "Risk_A0_ipcw_hi",    "Risk_A0_cf",  "analytical",
  "Risk_A1_ipcw",       "Risk_A1_ipcw",    "Risk_A1_ipcw_lo",    "Risk_A1_ipcw_hi",    "Risk_A1_cf",  "analytical",
  # IRR — analytical sandwich SE CIs
  "IRR_naive",          "IRR_naive_est",   "IRR_naive_lo",       "IRR_naive_hi",       "IRR_cf",    "analytical",
  "IRR_ipcw",           "IRR_ipcw_est",    "IRR_ipcw_lo",        "IRR_ipcw_hi",        "IRR_cf",    "analytical",
  # Arm-specific IRs — point estimates only
  "IR_A0_naive",        "IR_obs_A0",       "IR_obs_A0_lo",       "IR_obs_A0_hi",       "IR_cf_A0",  "analytical",
  "IR_A1_naive",        "IR_obs_A1",       "IR_obs_A1_lo",       "IR_obs_A1_hi",       "IR_cf_A1",  "analytical",
  "IR_A0_ipcw",         "IR_ipcw_A0",      "IR_ipcw_A0_lo",      "IR_ipcw_A0_hi",      "IR_cf_A0",  "analytical",
  "IR_A1_ipcw",         "IR_ipcw_A1",      "IR_ipcw_A1_lo",      "IR_ipcw_A1_hi",      "IR_cf_A1",  "analytical"
)

# ============================================================================
# COMPUTE PERFORMANCE SUMMARY
# ============================================================================

mc_summary_long <- mc_results_raw %>%
  group_by(scenario) %>%
  group_modify(function(.x, .y) {
    n_rep <- nrow(.x)
    map_dfr(seq_len(nrow(metric_specs)), function(i) {
      spec <- metric_specs[i, ]
      lo_v <- if (!is.na(spec$lo_col)) .x[[spec$lo_col]] else rep(NA_real_, n_rep)
      hi_v <- if (!is.na(spec$hi_col)) .x[[spec$hi_col]] else rep(NA_real_, n_rep)

      perf <- calculate_performance(
        estimates = .x[[spec$est_col]],
        lower_ci  = lo_v,
        upper_ci  = hi_v,
        truth     = .x[[spec$truth_col]]
      )

      bind_cols(
        tibble(
          metric      = spec$metric,
          ci_type     = spec$ci_type,
          truth_mean  = mean(.x[[spec$truth_col]], na.rm = TRUE)
        ),
        perf
      )
    })
  }) %>%
  ungroup()

# Weight diagnostics averaged across replicates
weight_diag <- mc_results_raw %>%
  group_by(scenario) %>%
  summarise(
    w_mean_avg          = mean(w_ipcw_mean,      na.rm = TRUE),
    w_max_avg           = mean(w_ipcw_max,       na.rm = TRUE),
    w_p95_avg           = mean(w_ipcw_p95,       na.rm = TRUE),
    n_events_A0_avg     = mean(n_events_A0,      na.rm = TRUE),
    n_events_A1_avg     = mean(n_events_A1,      na.rm = TRUE),
    pct_censored_A0_avg = mean(pct_censored_A0,  na.rm = TRUE),
    pct_censored_A1_avg = mean(pct_censored_A1,  na.rm = TRUE),
    .groups = "drop"
  )

# Formatted summary table
mc_summary_table <- mc_summary_long %>%
  mutate(
    empirical_result = if_else(
      !is.na(empirical_lo),
      sprintf("%.4f (%.4f, %.4f)", mean, empirical_lo, empirical_hi),
      sprintf("%.4f", mean)
    )
  ) %>%
  select(scenario, metric, ci_type, truth_mean, n_converged,
         empirical_result, bias, rel_bias_pct, empirical_se, mse,
         coverage, mean_ci_width)

# ============================================================================
# EXPORT
# ============================================================================

write_csv(mc_summary_long,  file.path(output_dir, "mc_summary_performance_long.csv"))
write_csv(mc_summary_table, file.path(output_dir, "mc_summary_performance.csv"))
write_csv(weight_diag,      file.path(output_dir, "mc_weight_diagnostics.csv"))

cat("Summary results exported.\n")
print(mc_summary_table, n = Inf, width = Inf)
