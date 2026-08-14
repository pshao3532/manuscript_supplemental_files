# ============================================================================
# DIAGNOSTIC FUNCTION — which covariates drive BOTH the event AND censoring?
#
# A covariate that predicts the Outcome (event) *and* non-administrative censoring is a
# shared cause => it has the potential to make censoring INFORMATIVE (and, if its censoring effect
# differs by arm, DIFFERENTIAL). For a given dataset this function fits two Cox
# models on the same covariates:
#     (1) Event model :  Surv(time, outcome_flag) ~ exposure + confounders
#     (2) Censor model:  Surv(time, censor_flag)  ~ exposure + confounders
# and reports HR (95% CI) for each covariate in each model (rows = {Event,
# Non-administrative censor}; cols = covariates). A covariate whose CI excludes
# 1 in BOTH models is flagged as an informative-censoring driver.
#
# When interaction = TRUE, the model also includes exposure:confounder terms
# (one per confounder) and those interaction HRs are exported alongside the
# main effects -- a formal CI test of whether an effect is DIFFERENTIAL by arm.
# ============================================================================

suppressMessages({
  library(survival)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# Needed only for the EXAMPLE USE block below (simulate_scenario1..4()).
# Not required if you're sourcing this file just for analyze_censor_outcome_drivers()
# and supplying your own data.
source("simulation_ipcw_github.R")

# ----------------------------------------------------------------------------
# analyze_censor_outcome_drivers()
#
# input_table             data.frame / tibble (one row per subject)
# exposure                <chr> column name of the exposure / treatment (e.g. "A")
# outcome_flag            <chr> column name of the event indicator (1 = event)
# censor_flag             <chr> column name of the non-administrative censoring
#                               indicator (1 = LTFU censored, not admin)
# time_to_event           <chr> column name of the follow-up time
# censoring_confounders   <chr> vector of candidate censoring-confounder columns
# interaction             logical; if TRUE add exposure:confounder terms and
#                               export their HRs too
# weight                  <chr> column name of weights to use (e.g. "IPTW") for
#                               a weighted Cox model; FALSE (default) fits an
#                               unweighted Cox model
# output_folder           <chr> directory to write the CSV(s) into
#
# Returns (invisibly) a list with $drivers and $driver_flags, and writes
#   <output_folder>/censor_outcome_drivers.csv
#   <output_folder>/censor_outcome_driver_flags.csv
# ----------------------------------------------------------------------------
analyze_censor_outcome_drivers <- function(input_table,
                                           exposure,
                                           outcome_flag,
                                           censor_flag,
                                           time_to_event,
                                           censoring_confounders,
                                           interaction   = FALSE,
                                           weight        = FALSE,
                                           output_folder) {

# ---- validate -----------------------------------------------------------
  needed <- c(exposure, outcome_flag, censor_flag, time_to_event,
              censoring_confounders)
  if (!isFALSE(weight))
    needed <- c(needed, weight)
  missing_cols <- setdiff(needed, names(input_table))
  if (length(missing_cols))
    stop("input_table is missing column(s): ", paste(missing_cols, collapse = ", "))
  if (!dir.exists(output_folder))
    dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# ---- helper: pull HR + 95% CI for one term out of a coxph fit -----------
  hr_ci <- function(fit, term) {
    if (is.null(fit))
      return(tibble(hr = NA_real_, lo = NA_real_, hi = NA_real_,
                    sig = NA, ci_txt = NA_character_))
    ci <- summary(fit)$conf.int
    if (!term %in% rownames(ci))
      return(tibble(hr = NA_real_, lo = NA_real_, hi = NA_real_,
                    sig = NA, ci_txt = NA_character_))
    hr  <- ci[term, "exp(coef)"]
    lo  <- ci[term, "lower .95"]
    hi  <- ci[term, "upper .95"]
    sig <- (lo > 1) | (hi < 1)                       # CI excludes 1?
    tibble(hr = hr, lo = lo, hi = hi, sig = sig,
           ci_txt = sprintf("%.2f (%.2f, %.2f)%s", hr, lo, hi,
                            ifelse(isTRUE(sig), "*", "")))
  }
  fit_cox <- function(form)
    if (isFALSE(weight))
      tryCatch(coxph(form, data = input_table, ties = "efron"),
               error = function(e) NULL)
    else
      tryCatch(coxph(form, data = input_table, weights = input_table[[weight]],
                     ties = "efron"),
               error = function(e) NULL)

# ---- build the two model formulas ---------------------------------------
  main_terms <- c(exposure, censoring_confounders)
  int_terms  <- if (isTRUE(interaction))
    paste0(exposure, ":", censoring_confounders) else character(0)
  rhs <- paste(c(main_terms, int_terms), collapse = " + ")

  ev_form <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_to_event, outcome_flag, rhs))
  cn_form <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_to_event, censor_flag,  rhs))
  fit_ev  <- fit_cox(ev_form)
  fit_cn  <- fit_cox(cn_form)

# ---- assemble HR table over all reported terms --------------------------
  report_terms <- c(main_terms, int_terms)
  term_type    <- c(rep("main", length(main_terms)),
                    rep("interaction", length(int_terms)))

  drivers <- do.call(rbind, lapply(seq_along(report_terms), function(j) {
    term <- report_terms[j]
    rbind(
      cbind(term = term, term_type = term_type[j],
            model = "Event",                       hr_ci(fit_ev, term)),
      cbind(term = term, term_type = term_type[j],
            model = "Non-administrative censor",   hr_ci(fit_cn, term))
    )
  }))
  drivers <- as_tibble(drivers)

# ---- driver flag: CI excludes 1 in BOTH models (main effects only) ------
  driver_flags <- drivers %>%
    filter(term_type == "main") %>%
    select(term, model, sig) %>%
    pivot_wider(names_from = model, values_from = sig) %>%
    rename(sig_event = Event, sig_censor = `Non-administrative censor`) %>%
    mutate(informative_censoring_driver = dplyr::coalesce(sig_event, FALSE) &
                                          dplyr::coalesce(sig_censor, FALSE))

# ---- report -------------------------------------------------------------
  cat("\n================= HR (95% CI)  ('*' = CI excludes 1) =================\n")
  shell <- drivers %>%
    filter(term_type == "main") %>%
    select(model, term, ci_txt) %>%
    pivot_wider(names_from = term, values_from = ci_txt)
  print(as.data.frame(shell[, c("model", main_terms)]), row.names = FALSE)

  cat("\n=============== INFORMATIVE-CENSORING DRIVERS ========================\n")
  cat("(TRUE = covariate predicts the outcome AND non-admin censoring)\n\n")
  print(as.data.frame(driver_flags), row.names = FALSE)

  if (isTRUE(interaction) && length(int_terms)) {
    cat("\n=============== DIFFERENTIAL EFFECT (exposure x confounder) ==========\n")
    cat("CI excluding 1 => that effect differs by exposure arm.\n\n")
    print(as.data.frame(
      drivers %>% filter(term_type == "interaction") %>%
        select(model, term, ci_txt)), row.names = FALSE)
  }

# ---- save ---------------------------------------------------------------
  table_label <- deparse(substitute(input_table)) 
  drivers_path <- file.path(output_folder, paste0("censor_outcome_drivers_", table_label,".csv") )
  flags_path   <- file.path(output_folder, paste0("censor_outcome_driver_flags_",table_label,".csv") )
  write.csv(drivers,      drivers_path, row.names = FALSE)
  write.csv(driver_flags, flags_path,   row.names = FALSE)
  cat(sprintf("\nSaved: %s\n       %s\n", drivers_path, flags_path))

  invisible(list(drivers = drivers, driver_flags = driver_flags))
}

# ----------------------------------------------------------------------------
# EXAMPLE USE:
# Generates one draw of each scenario (same n and max_fu as the manuscript)
# via simulate_scenario1..4() from simulation_ipcw_github.R, then checks
# which covariates predict both the event and non-administrative censoring.
#

set.seed(20260606)
dat_s1 <- simulate_scenario1(n = 20000, dt = 1, max_fu = 24)$dat
dat_s2 <- simulate_scenario2(n = 20000, dt = 1, max_fu = 24)$dat
dat_s3 <- simulate_scenario3(n = 20000, dt = 1, max_fu = 24)$dat
dat_s4 <- simulate_scenario4(n = 20000, dt = 1, max_fu = 24)$dat

  analyze_censor_outcome_drivers(
    input_table           = dat_s1,
    exposure              = "A",
    outcome_flag          = "event1",
    censor_flag           = "censor_event",
    time_to_event         = "time",
    censoring_confounders = c("L"),
    interaction           = TRUE,
    output_folder         = "results"
  )
# ----------------------------------------------------------------------------
analyze_censor_outcome_drivers(
  input_table           = dat_s2,
  exposure              = "A",
  outcome_flag          = "event1",
  censor_flag           = "censor_event",
  time_to_event         = "time",
  censoring_confounders = c("L_c1", "L_c2", "L"),
  interaction           = TRUE,
  output_folder         = "results"
)

analyze_censor_outcome_drivers(
  input_table           = dat_s3,
  exposure              = "A",
  outcome_flag          = "event1",
  censor_flag           = "censor_event",
  time_to_event         = "time",
  censoring_confounders = c("L_c1", "L_c2", "L"),
  interaction           = TRUE,
  output_folder         = "results"
)

analyze_censor_outcome_drivers(
  input_table           = dat_s4,
  exposure              = "A",
  outcome_flag          = "event1",
  censor_flag           = "censor_event",
  time_to_event         = "time",
  censoring_confounders = c("L_c1", "L_c2", "L"),
  interaction           = TRUE,
  output_folder         = "results"
)
