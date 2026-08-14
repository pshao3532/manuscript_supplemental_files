# ============================================================================
# Censoring SIMULATION
#   Differential Censoring Is Not Informative Censoring: A Simulation Study
#   of Censoring Bias in Survival Analyses -- manuscript programming code share.
# 
# This script was create to simulate 4 censoring scenarios and related analysis
# ============================================================================
# Author: Ping
# Date: 2026-05-03
# ============================================================================

library(tidyverse)
library(survival)
library(riskRegression)
library(survRM2)
library(data.table)
library(sandwich)
library(parallel)   


# Set the uniform coefficient for exposure in event generation models
betaA       <- log(3.0)
# ============================================================================
# SECTION 1: DATA GENERATION
# ============================================================================

generate_baseline <- function(n = 2000) {
  L      <- rbinom(n, 1, 0.5)
  alpha0 <- -0.5
  alphaL <- 1.2
  prob_A <- plogis(alpha0 + alphaL * L)
  A      <- rbinom(n, 1, prob_A)
  tibble(id = 1:n, L = L, A = A)
}

make_person_time <- function(dat, dt_ = 1, tmax = NULL) {
  stopifnot(all(c("id", "time", "event1") %in% names(dat)))
  if (is.null(tmax)) tmax <- max(dat$time, na.rm = TRUE)

  has_counterfactuals <- all(c("T_counterfactual", "C_counterfactual") %in% names(dat))
  dt <- data.table::as.data.table(dat)

  if (!"censor_event" %in% names(dt)) dt[, censor_event := 1L - event1]

  dt[, t_obs       := pmin(time, tmax)]
  dt[, n_intervals := pmax(1L, ceiling(t_obs / dt_))]

  expand_idx   <- rep(seq_len(nrow(dt)), dt$n_intervals)
  cols_to_keep <- setdiff(names(dt), c("time", "T_counterfactual", "C_counterfactual", "n_intervals"))
  pt           <- dt[expand_idx, ..cols_to_keep]

  pt[, k     := sequence(dt$n_intervals)]
  pt[, start := (k - 1L) * dt_]
  pt[, stop  := pmin(k * dt_, t_obs)]
  pt <- pt[stop > start]
  pt[, cens_interval  := as.integer(censor_event == 1L & stop == t_obs)]
  pt[, event_interval := as.integer(event1 == 1L & stop == t_obs)]

  if (has_counterfactuals) {
    cf_data <- dt[, .(
      id,
      T_counterfactual = pmin(T_counterfactual, tmax + 0.0001),
      time,
      C_counterfactual = pmin(C_counterfactual, tmax)
    )]
    pt <- cf_data[pt, on = "id"]
  }

  tibble::as_tibble(pt)
}

# Scenario 1: Non-informative, non-differential censoring (no censoring risk factors)
simulate_scenario1 <- function(n = 2000, dt = 1, max_fu = 24) {
  dat <- generate_baseline(n)

  baseline_h1 <- 0.03
  betaL       <- 0.0
  dat$T_counterfactual <- rexp(n, rate = baseline_h1 * exp(betaL * dat$L + betaA * dat$A))

  baseline_hC <- 0.07
  gammaL      <- 0.0
  dat$C_ltfu  <- rexp(n, rate = baseline_hC * exp(gammaL * dat$L))
  dat$C_counterfactual <- pmin(dat$C_ltfu, max_fu)

  dat <- dat %>%
    mutate(
      time         = pmin(T_counterfactual, C_counterfactual),
      event1       = as.integer(T_counterfactual <= pmin(C_ltfu, max_fu)),
      censor_event = as.integer(C_ltfu < pmin(T_counterfactual, max_fu)),
      admin_cens   = as.integer(time >= max_fu & event1 == 0 & censor_event == 0),
      scenario     = "S1_noninf_nondiff"
    )
  dat_pt<- make_person_time(dat, dt_ = dt)
  list(dat = dat,dat_pt= dat_pt)
}

# Scenario 2: Non-informative (not associated with event risk), DIFFERENTIAL censoring.

simulate_scenario2 <- function(n = 10000, dt = 1, max_fu = 24) {
  dat <- generate_baseline(n)

  dat$L_c1 <- rbinom(n, 1, ifelse(dat$A == 1, 0.30, 0.02))
  dat$L_c2 <- rbinom(n, 1, ifelse(dat$A == 1, 0.20, 0.02))

  baseline_h1 <- 0.02
  betaL       <- 0.0
  dat$T_counterfactual <- rexp(n, rate = baseline_h1 * exp(betaL * dat$L + betaA * dat$A))

  baseline_hC <- 0.04
  gammaN1     <- 0.4
  gammaN2     <- 0.8
  dat$C_ltfu  <- rexp(n, rate = baseline_hC * exp(gammaN1 * dat$L_c1 + gammaN2 * dat$L_c2))
  dat$C_counterfactual <- pmin(dat$C_ltfu, max_fu)

  dat <- dat %>%
    mutate(
      time         = pmin(T_counterfactual, C_counterfactual),
      event1       = as.integer(T_counterfactual <= pmin(C_ltfu, max_fu)),
      censor_event = as.integer(C_ltfu < pmin(T_counterfactual, max_fu)),
      admin_cens   = as.integer(time >= max_fu & event1 == 0 & censor_event == 0),
      scenario     = "S2_noninf_diff_nausea"
    )
  dat_pt<- make_person_time(dat, dt_ = dt)
  list(dat = dat,dat_pt= dat_pt)
}

# Scenario 3: INFORMATIVE, non-differential censoring.

simulate_scenario3 <- function(n = 10000, dt = 1, max_fu = 24, betaA = log(3.0)) {
  dat <- generate_baseline(n)

  dat$L_c1 <- rbinom(n, 1, 0.35)
  dat$L_c2 <- rbinom(n, 1, 0.25)

  baseline_h1 <- 0.015
  betaL       <- 0
  betaL_c1    <- log(4) * 0.4   # 40% share of log-HR; HR for L_c1=1 alone: 
  betaL_c2    <- log(4) * 0.6   # 60% share of log-HR; HR for L_c2=1 alone:

  dat$T_counterfactual <- rexp(
    n,
    rate = baseline_h1 * exp(
      betaA * dat$A + betaL * dat$L +
        betaL_c1 * dat$L_c1 + betaL_c2 * dat$L_c2
    )
  )

  baseline_hC <- 0.025
  gammaL_c1   <- log(6) * 0.4   # 40% share; HR for L_c1=1 
  gammaL_c2   <- log(6) * 0.6   # 60% share; HR for L_c2=1 
  dat$C_ltfu  <- rexp(n, rate = baseline_hC * exp(gammaL_c1 * dat$L_c1 + gammaL_c2 * dat$L_c2))
  dat$C_counterfactual <- pmin(dat$C_ltfu, max_fu)

  dat <- dat %>%
    mutate(
      time         = pmin(T_counterfactual, C_counterfactual),
      event1       = as.integer(T_counterfactual <= pmin(C_ltfu, max_fu)),
      censor_event = as.integer(C_ltfu < pmin(T_counterfactual, max_fu)),
      admin_cens   = as.integer(time >= max_fu & event1 == 0 & censor_event == 0),
      scenario     = "S3_inf_nondiff"
    )

  dat_pt<- make_person_time(dat, dt_ = dt)
  list(dat = dat,dat_pt= dat_pt)
}

# Scenario 4: INFORMATIVE AND DIFFERENTIAL censoring.
# A:L_c interaction on censoring kept for L_c1 only.
simulate_scenario4 <- function(n = 10000, dt = 1, max_fu = 24, betaA = log(3.0)) {
  L    <- rbinom(n, 1, 0.5)
  L_c1 <- rbinom(n, 1, 0.3)
  L_c2 <- rbinom(n, 1, 0.4)
  alpha0 <- -0.5
  alphaL <- 1.2
  prob_A <- plogis(alpha0 + alphaL * L)
  A      <- rbinom(n, 1, prob_A)
  dat    <- tibble(id = 1:n, L = L, A = A, L_c1 = L_c1, L_c2 = L_c2)

  baseline_h1 <- 0.015
  betaL       <- 0
  betaL_c1    <- log(4) * 0.4
  betaL_c2    <- log(4) * 0.6
  dat$T_counterfactual <- rexp(
    n,
    rate = baseline_h1 * exp(
      betaA * dat$A + betaL * dat$L +
        betaL_c1 * dat$L_c1 + betaL_c2 * dat$L_c2
    )
  )

  baseline_hC  <- 0.01
  gammaL_c1    <- log(4) * 0.4
  gammaL_c2    <- log(4) * 0.6
  gamma_A      <- log(2.5)
  gammaA_L_c1  <- log(3)          # A:L_c1 interaction; L_c2 has no interaction term
  dat$C_ltfu   <- rexp(
    n,
    rate = baseline_hC * exp(
      gamma_A     * dat$A +
        gammaL_c1   * dat$L_c1 +
        gammaL_c2   * dat$L_c2 +
        gammaA_L_c1 * dat$A * dat$L_c1
    )
  )
  dat$C_counterfactual <- pmin(dat$C_ltfu, max_fu)

  dat <- dat %>%
    mutate(
      time         = pmin(T_counterfactual, C_counterfactual),
      event1       = as.integer(T_counterfactual <= pmin(C_ltfu, max_fu)),
      censor_event = as.integer(C_ltfu < pmin(T_counterfactual, max_fu)),
      admin_cens   = as.integer(time >= max_fu & event1 == 0 & censor_event == 0),
      scenario     = "S4_inf_diff"
    )

  dat_pt<- make_person_time(dat, dt_ = dt)
  list(dat = dat,dat_pt= dat_pt)
}

# ============================================================================
# SECTION 2: PERSON-LEVEL DATA & IPCW VIA riskRegression::ipcw()
# ============================================================================
# Build the correctly-specified censoring model formula for each scenario.
# The denominator includes all variables that predict LTFU censoring in the DGP.

censor_surv_formula <- function(scenario_id, time_var, status_var) {
  lhs <- paste0("survival::Surv(", time_var, ", ", status_var, ")")

  rhs <- switch(
    scenario_id,
    S1_noninf_nondiff     = "1",
    S2_noninf_diff_nausea = "A+ L_c1 + L_c2",          # A acts through L_c1/L_c2, not directly
    S3_inf_nondiff        = "A + L_c1 + L_c2",          # A has no effect on censoring in S3 DGP
    S4_inf_diff           = "L_c1 + L_c2 + A + L_c1:A",
    stop("Unknown scenario_id: ", scenario_id)
  )

  as.formula(paste(lhs, rhs, sep = " ~ "))
}

# Compute stabilised IPCW weights via riskRegression::ipcw().
# dat_pt: counting-process person-time data (start/stop format).
# dat:    person-level data (one row per person, from simulate_scenario*()$dat).
# Denominator: full censoring model (all predictors from DGP).
# Numerator:   marginal for S1; arm-conditional for S2–S4 (standard stabilisation).
estimate_ipcw_weights_riskregression <- function(dat_pt =dat_s3$dat_pt,
                                                 dat= dat_s3$dat,
                                                 scenario_id = "S3_inf_nondiff",
                                                 tmax = 36,
                                                 eps        = 1e-6,
                                                 truncate   = NULL,
                                                 eval_times = NULL) {  
  # R1: data.table filtering is 3-5x faster than dplyr
  pt <- data.table::as.data.table(dat_pt)
  pt[, ptime := stop - start]
  pt <- pt[ptime > 0 & start < tmax]
  data.table::setorder(pt, id, start)

  # ipcw() uses status=0 meaning censored, then internally reverses so that
  # status=0 becomes the "censoring event". We want LTFU (censor_event=1) to
  # be that event, so ltfu_status = 1 - censor_event gives 0 for LTFU (the
  # event we want to model) and 1 for clinical event / admin censoring.
  # dat is kept unmodified for use in counterfactual_truth / calc_naive_estimates.
  dat_cens <- dat %>% mutate(ltfu_status = 1L - censor_event)

  formula_den <- censor_surv_formula(scenario_id, "time", "ltfu_status")

  # Numerator: arm-conditional (~ A) for all scenarios with differential or
  # arm-dependent censoring (S2–S4); marginal (~ 1) only for S1.
  formula_num <- if (scenario_id == "S1_noninf_nondiff") {
    survival::Surv(time, ltfu_status) ~ 1
  } else {
    survival::Surv(time, ltfu_status) ~ A
  }

  method_den <- if (scenario_id == "S1_noninf_nondiff") "marginal" else "cox"
  method_num <- if (scenario_id == "S1_noninf_nondiff") "marginal" else "cox"

  # R1: reuse pre-computed eval_times from caller when available (saves sort/unique in bootstrap)
  if (is.null(eval_times)) eval_times <- sort(unique(pt$start))

  # lag = 1 gives G(t-) = P(C_ltfu >= t), the left-continuous version required by IPCW.
  fit_den <- riskRegression::ipcw(
    formula       = formula_den,
    data          = dat_cens,
    method        = method_den,
    times         = eval_times,
    subject.times = dat_cens$time,
    lag           = 1,
    keep          = c("fit")
  )

  fit_num <- riskRegression::ipcw(
    formula       = formula_num,
    data          = dat_cens,
    method        = method_num,
    times         = eval_times,
    subject.times = dat_cens$time,
    lag           = 1,
    keep          = c("fit")
  )

  id_to_row   <- match(pt$id, dat_cens$id)
  time_to_col <- match(pt$start, eval_times)

  get_ipcw_at_start <- function(fit) {
    if (is.null(dim(fit$IPCW.times))) {
      fit$IPCW.times[time_to_col]
    } else {
      fit$IPCW.times[cbind(id_to_row, time_to_col)]
    }
  }

  G_den <- get_ipcw_at_start(fit_den)
  G_num <- get_ipcw_at_start(fit_num)
  w     <- G_num / pmax(G_den, eps)
  w[!is.finite(w) | is.na(w)] <- 1

  if (!is.null(truncate)) {
    qs <- quantile(w, probs = truncate, na.rm = TRUE)
    w  <- pmin(pmax(w, qs[1]), qs[2])
  }

  pt$w_ipcw <- w
  pt
}

estimate_ipcw_weights <- function(dat = dat_s3,
                                  scenario_id,
                                  tmax = NULL,
                                  eps = 1e-6,
                                  truncate = NULL) {
  
  pt <- dat%>%
    mutate(ptime = stop - start) %>% 
    filter(ptime >0) 
  
  # Define censoring models by scenario
  if (scenario_id == "S1_noninf_nondiff") {
    formula_den <- cens_interval ~ 1 + offset(log(ptime))
    formula_num <- cens_interval ~ 1 + offset(log(ptime))
  } else if (scenario_id == "S2_noninf_diff_nausea") {
    formula_den <- cens_interval ~ L_c1 + L_c2+  + A + offset(log(ptime))
    formula_num <- cens_interval ~ A + offset(log(ptime))
  } else if (scenario_id == "S3_inf_nondiff") {
    # L_c is baseline covariate, include in denominator
    formula_den <- cens_interval ~ L_c1 + L_c2 + A+ offset(log(ptime))
    formula_num <- cens_interval ~  A + offset(log(ptime))
  } else if (scenario_id == "S4_inf_diff") {
    # Differential censoring: include A and A:L_c
    formula_den <- cens_interval ~  L_c1 + L_c2 + A  + L_c1:A + offset(log(ptime))
    formula_num <- cens_interval ~  A+ offset(log(ptime))
  } else {
    stop("Unknown scenario_id")
  }
  
  pt_fit <- pt %>% filter(event_interval == 0)  # fit censoring model here
  # Fit pooled logistic models
  fit_den <- glm(formula_den, data = pt_fit, family = binomial(link = "cloglog"))
  fit_num <- glm(formula_num, data = pt_fit, family = binomial(link = "cloglog"))
  
  # Predicted probabilities
  p_den <- predict(fit_den, newdata = pt, type = "response")
  p_num <- predict(fit_num, newdata = pt, type = "response")
  
  p_den <- pmin(pmax(p_den, eps), 1 - eps)
  p_num <- pmin(pmax(p_num, eps), 1 - eps)
  
  # Calculate weights
  pt_weighted <- pt %>%
    arrange(id, start) %>%
    mutate(
      q_den = 1 - p_den,
      q_num = 1 - p_num
    ) %>%
    group_by(id) %>%
    mutate(
      G_den_stop = cumprod(q_den),
      G_num_stop = cumprod(q_num),
      G_den_start = lag(G_den_stop, default = 1),
      G_num_start = lag(G_num_stop, default = 1),
      w_ipcw = G_num_start / G_den_start
    ) %>%
    ungroup()
  
  # Clean up weights
  w <- pt_weighted$w_ipcw
  w[!is.finite(w) | is.na(w)] <- 1
  
  if (!is.null(truncate)) {
    qs <- quantile(w, probs = truncate, na.rm = TRUE)
    w <- pmin(pmax(w, qs[1]), qs[2])
  }
  
  pt_weighted$w_ipcw <- w
  
  # Counterfactuals are already in pt_weighted from make_person_time()
  return(pt_weighted)
}

# ============================================================================
# SECTION 3: COUNTERFACTUAL TRUTH
# ============================================================================

calc_rmst_from_survfit <- function(fit, tau, arm) {
  s           <- summary(fit)
  strata_name <- paste0("A=", arm)
  idx         <- which(s$strata == strata_name)
  times       <- s$time[idx]
  surv        <- s$surv[idx]
  keep        <- times <= tau

  times_full  <- c(0, times[keep], tau)
  surv_full   <- c(1, surv[keep], if (any(keep)) tail(surv[keep], 1) else 1)
  keep_unique <- !duplicated(times_full)
  sum(diff(times_full[keep_unique]) * head(surv_full[keep_unique], -1))
}

# dat: person-level data from simulate_scenario*()$dat; contains T_counterfactual.
counterfactual_truth <- function(dat = dat_s3$dat, tmax = 36) {
  dat_cf <- dat %>%
    transmute(
      id       = id,
      A        = A,
      time_cf  = pmin(T_counterfactual, tmax),
      event_cf = as.integer(T_counterfactual <= tmax)
    )

  fit_cf_cox <- survival::coxph(survival::Surv(time_cf, event_cf) ~ A, data = dat_cf,
                                control = survival::coxph.control(timefix = FALSE))
  fit_cf_km  <- survival::survfit(survival::Surv(time_cf, event_cf) ~ A, data = dat_cf,
                                  timefix = FALSE)

  surv_cf <- summary(fit_cf_km, times = tmax, extend = TRUE)
  s0_cf   <- surv_cf$surv[surv_cf$strata == "A=0"]
  s1_cf   <- surv_cf$surv[surv_cf$strata == "A=1"]

  list(
    hr         = exp(coef(fit_cf_cox)[["A"]]),
    rd         = (1 - s1_cf) - (1 - s0_cf),
    rmstd      = calc_rmst_from_survfit(fit_cf_km, tmax, 1) -
                   calc_rmst_from_survfit(fit_cf_km, tmax, 0),
    risk_a0_cf = 1 - s0_cf,
    risk_a1_cf = 1 - s1_cf
  )
}

# ============================================================================
# SECTION 4: NAIVE (UNWEIGHTED) ESTIMATES
# ============================================================================

# dat: person-level data from simulate_scenario*()$dat.
calc_naive_estimates <- function(dat, tau) {

  # ---- HR: Cox with robust sandwich SE clustered by id ----
  fit_hr_cox <- survival::coxph(
    survival::Surv(time, event1) ~ A + cluster(id),
    data    = dat,
    control = survival::coxph.control(timefix = FALSE)
  )
  beta_hr <- coef(fit_hr_cox)[["A"]]
  se_hr   <- sqrt(vcov(fit_hr_cox)["A", "A"])
  hr <- list(
    hr = exp(beta_hr),
    lo = exp(beta_hr - 1.96 * se_hr),
    hi = exp(beta_hr + 1.96 * se_hr)
  )

  # ---- RD: Kaplan-Meier with Greenwood SE at tau ----
  fit_km <- survival::survfit(survival::Surv(time, event1) ~ A, data = dat,
                              timefix = FALSE)
  s_tau  <- summary(fit_km, times = tau, extend = TRUE)
  s0     <- s_tau$surv[s_tau$strata == "A=0"]
  s1     <- s_tau$surv[s_tau$strata == "A=1"]
  se0    <- s_tau$std.err[s_tau$strata == "A=0"]
  se1    <- s_tau$std.err[s_tau$strata == "A=1"]

  rd_est <- (1 - s1) - (1 - s0)
  rd_se  <- sqrt(se0^2 + se1^2)
  rd <- list(
    rd = rd_est,
    lo = rd_est - 1.96 * rd_se,
    hi = rd_est + 1.96 * rd_se
  )

  # ---- RMSTD: survRM2::rmst2() ----
  dat_rmst <- dat %>%
    select(time, status = event1, A) %>%
    as.data.frame()

  rmst2_fit <- tryCatch(
    survRM2::rmst2(
      time   = dat_rmst$time,
      status = dat_rmst$status,
      arm    = dat_rmst$A,
      tau    = tau
    ),
    error = function(e) NULL
  )

  if (!is.null(rmst2_fit)) {
    res_row <- rmst2_fit$unadjusted.result[1L, ]   # row 1 = arm1 - arm0 difference
    rmstd <- list(
      rmstd = unname(res_row[1L]),   # Est.
      lo    = unname(res_row[2L]),   # lower .95
      hi    = unname(res_row[3L])    # upper .95
    )
  } else {
    rmstd <- list(rmstd = NA_real_, lo = NA_real_, hi = NA_real_)
  }

  list(
    hr               = hr,
    rd               = rd,
    rmstd            = rmstd,
    risk_a0_naive    = 1 - s0,
    risk_a0_naive_lo = pmax(0, (1 - s0) - 1.96 * se0),
    risk_a0_naive_hi = pmin(1, (1 - s0) + 1.96 * se0),
    risk_a1_naive    = 1 - s1,
    risk_a1_naive_lo = pmax(0, (1 - s1) - 1.96 * se1),
    risk_a1_naive_hi = pmin(1, (1 - s1) + 1.96 * se1)
  )
}

# ============================================================================
# SECTION 5: RD and RMSTD ESTIMATES VIA WEIGHTED KM
# ============================================================================

# IPCW-weighted KM risk difference at tau -- POINT ESTIMATE ONLY.
# v11: the Greenwood plug-in SE/CI (which used raw, unweighted counts and
# understated the weighted estimator's variance) has been removed. lo/hi/se0/se1
# are returned as NA so downstream code keeps a stable structure; RD inference
# now comes from the bootstrap (RD_ipcw_boot_*). fit_raw is unused (kept for a
# stable signature).
extract_rd_ci_corrected <- function(fit_weighted, fit_raw, tau) {
  s_w <- summary(fit_weighted, times = tau, extend = TRUE)
  s0  <- s_w$surv[s_w$strata == "A=0"]
  s1  <- s_w$surv[s_w$strata == "A=1"]
  rd  <- (1 - s1) - (1 - s0)
  list(rd = rd, lo = NA_real_, hi = NA_real_,
       se0 = NA_real_, se1 = NA_real_)
}

# IPCW-weighted RMST up to tau -- POINT ESTIMATE ONLY.
# v11: the Greenwood-style plug-in variance (raw, unweighted counts) has been
# removed because it understated the weighted estimator's variance. se is
# returned as NA; RMSTD inference now comes from the bootstrap. fit_raw is
# unused (kept for a stable signature).
calc_rmst_with_ci_corrected <- function(fit_weighted, fit_raw, tau, arm) {
  strata_name <- paste0("A=", arm)

  sw      <- summary(fit_weighted)
  idx_w   <- which(sw$strata == strata_name)
  times_w <- sw$time[idx_w]
  surv_w  <- sw$surv[idx_w]

  keep       <- times_w <= tau
  times_k    <- times_w[keep]
  surv_k     <- surv_w[keep]

  times_full <- c(0, times_k, tau)
  surv_full  <- c(1, surv_k, if (length(surv_k)) tail(surv_k, 1) else 1)
  dup        <- duplicated(times_full)
  rmst       <- sum(diff(times_full[!dup]) * head(surv_full[!dup], -1))

  list(rmst = rmst, se = NA_real_)
}

# HR with IPCW: weighted Cox with robust sandwich SE clustered by id.
# pt_weighted is pre-computed in analyze_coxph_pt() and shared with calc_rd_rmstd_ipcw().
calc_hr_ipcw <- function(pt_weighted = pt, tau) {
  dat_cox <- pt_weighted %>%
    filter(stop > start, start < tau) %>%
    mutate(event = as.integer(event_interval))

  fit <- survival::coxph(
    survival::Surv(start, stop, event) ~ A + cluster(id),
    data    = dat_cox,
    weights = w_ipcw,
    control = survival::coxph.control(timefix = FALSE)
  )

  beta <- coef(fit)[["A"]]
  se   <- sqrt(vcov(fit)["A", "A"])
  list(hr = exp(beta), lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se))
}

# RD, RMSTD, and arm-specific absolute risks from IPCW-weighted KM.
# Reuses pt_weighted (same LTFU-correct weights used for HR) — no ate() call.
# CI uses corrected Greenwood: weighted KM survival, raw counts in variance.
calc_rd_rmstd_ipcw <- function(pt_weighted, dat, tau) {
  # Drop zero-/near-zero-length intervals so survival's aeqSurv (timefix) does
  # not throw "an interval has effective length 0"; disable timefix to match
  # calc_hr_ipcw (interval endpoints are exact by construction).
  pt_weighted <- pt_weighted %>% filter(stop > start)

  # IPCW-weighted KM (counting-process format for time-varying weights)
  fit_ipcw_km <- survival::survfit(
    survival::Surv(start, stop, event_interval) ~ A,
    data    = pt_weighted,
    weights = w_ipcw,
    timefix = FALSE
  )

  # Naive KM on person-level data supplies raw n.risk / n.event for Greenwood
  fit_naive_km <- survival::survfit(
    survival::Surv(time, event1) ~ A,
    data    = dat,
    timefix = FALSE
  )

  rd   <- tryCatch(extract_rd_ci_corrected(fit_ipcw_km, fit_naive_km, tau),
                   error = function(e) list(rd = NA, lo = NA, hi = NA))

  rmst_a0 <- tryCatch(calc_rmst_with_ci_corrected(fit_ipcw_km, fit_naive_km, tau, 0),
                      error = function(e) list(rmst = NA, se = NA))
  rmst_a1 <- tryCatch(calc_rmst_with_ci_corrected(fit_ipcw_km, fit_naive_km, tau, 1),
                      error = function(e) list(rmst = NA, se = NA))

  rmstd_est <- rmst_a1$rmst - rmst_a0$rmst

  # Arm-specific IPCW absolute risks at tau from weighted KM
  s_tau        <- summary(fit_ipcw_km, times = tau, extend = TRUE)
  risk_a0_ipcw <- 1 - s_tau$surv[s_tau$strata == "A=0"]
  risk_a1_ipcw <- 1 - s_tau$surv[s_tau$strata == "A=1"]

  # v11: Greenwood plug-in CIs removed for all IPCW absolute estimands.
  # Point estimates only; analytical CI bounds are NA. RD/RMSTD inference is via
  # the bootstrap (RD_ipcw_boot_*, RMSTD_ipcw_boot_*); the arm-specific IPCW
  # risks have no analytical CI (coverage reported as NA in the summary).
  list(
    rd    = rd,
    rmstd = list(
      rmstd = rmstd_est,
      lo    = NA_real_,
      hi    = NA_real_
    ),
    risk_a0_ipcw    = risk_a0_ipcw,
    risk_a0_ipcw_lo = NA_real_,
    risk_a0_ipcw_hi = NA_real_,
    risk_a1_ipcw    = risk_a1_ipcw,
    risk_a1_ipcw_lo = NA_real_,
    risk_a1_ipcw_hi = NA_real_
  )
}

# ============================================================================
# SECTION 5b: INCIDENCE RATE RATIO VIA POISSON REGRESSION
# ============================================================================

# Computes IRR_naive and IRR_ipcw from interval Poisson models with clustered
# sandwich SE. Arm-specific IRs are also returned (observed, IPCW-weighted,
# counterfactual). dat is passed in to avoid re-deriving from pt_weighted.

calc_irr <- function(dat_pois, dat, tmax) {

  f_pois <- y ~ A + offset(log(ptime))

  # R4: NAIVE (unweighted) Poisson fit on PERSON-AGGREGATED data.
  # NOTE: the IPCW fit CANNOT be aggregated the same way because w_ipcw varies
  # within person across intervals — it stays on the interval-level data.
  dat_pois_dt   <- data.table::as.data.table(dat_pois)
  dat_naive_agg <- dat_pois_dt[, .(y = sum(y), ptime = sum(ptime)),
                               by = .(id, A)]

  fit_naive <- glm(f_pois, family = poisson(), data = dat_naive_agg)
  fit_ipcw  <- glm(f_pois, family = poisson(), data = dat_pois, weights = w_ipcw)

  V_naive    <- sandwich::vcovCL(fit_naive, cluster = dat_naive_agg$id)
  V_ipcw     <- sandwich::vcovCL(fit_ipcw,  cluster = dat_pois$id)
  coef_naive <- coef(fit_naive)["A"]; se_naive <- sqrt(V_naive["A", "A"])
  coef_ipcw  <- coef(fit_ipcw)["A"];  se_ipcw  <- sqrt(V_ipcw["A",  "A"])

  # Crude (unweighted) IR and observed person-time by arm
  ir_obs <- dat_pois %>%
    group_by(A) %>%
    summarise(events = sum(y), ptime = sum(ptime), .groups = "drop") %>%
    mutate(IR = events / ptime)

  # IPCW-weighted IR by arm
  ir_ipcw_arm <- dat_pois %>%
    group_by(A) %>%
    summarise(
      events = sum(y * w_ipcw),
      ptime  = sum(ptime * w_ipcw),
      .groups = "drop"
    ) %>%
    mutate(IR = events / ptime)

  # Counterfactual IR from dat (person-level; avoids group_by/slice on pt_weighted)
  ir_cf <- dat %>%
    transmute(A, ptime_cf = pmin(T_counterfactual, tmax),
              event_cf = as.numeric(T_counterfactual <= tmax)) %>%
    group_by(A) %>%
    summarise(events = sum(event_cf), ptime = sum(ptime_cf), .groups = "drop") %>%
    mutate(IR = events / ptime)

  # Exact Poisson (chi-squared) CIs for observed IR by arm
  e0 <- ir_obs$events[ir_obs$A == 0]; pt0 <- ir_obs$ptime[ir_obs$A == 0]
  e1 <- ir_obs$events[ir_obs$A == 1]; pt1 <- ir_obs$ptime[ir_obs$A == 1]

  # Delta-method CIs for IPCW-weighted IR by arm from fit_ipcw / V_ipcw
  int_ipcw <- coef(fit_ipcw)["(Intercept)"]
  se_log_A0 <- sqrt(V_ipcw["(Intercept)", "(Intercept)"])
  se_log_A1 <- sqrt(V_ipcw["(Intercept)", "(Intercept)"] +
                    V_ipcw["A", "A"] +
                    2 * V_ipcw["(Intercept)", "A"])

  list(
    IRR_naive_est  = exp(coef_naive),
    IRR_naive_lo   = exp(coef_naive - 1.96 * se_naive),
    IRR_naive_hi   = exp(coef_naive + 1.96 * se_naive),
    IRR_ipcw_est   = exp(coef_ipcw),
    IRR_ipcw_lo    = exp(coef_ipcw  - 1.96 * se_ipcw),
    IRR_ipcw_hi    = exp(coef_ipcw  + 1.96 * se_ipcw),
    # Observed IR by arm + exact Poisson CI
    IR_obs_A0      = e0 / pt0,
    IR_obs_A0_lo   = qchisq(0.025, 2 * e0)       / (2 * pt0),
    IR_obs_A0_hi   = qchisq(0.975, 2 * (e0 + 1)) / (2 * pt0),
    IR_obs_A1      = e1 / pt1,
    IR_obs_A1_lo   = qchisq(0.025, 2 * e1)       / (2 * pt1),
    IR_obs_A1_hi   = qchisq(0.975, 2 * (e1 + 1)) / (2 * pt1),
    # IPCW-weighted IR by arm + delta-method CI (log scale, from fit_ipcw)
    IR_ipcw_A0     = exp(int_ipcw),
    IR_ipcw_A0_lo  = exp(int_ipcw - 1.96 * se_log_A0),
    IR_ipcw_A0_hi  = exp(int_ipcw + 1.96 * se_log_A0),
    IR_ipcw_A1     = exp(int_ipcw + coef_ipcw),
    IR_ipcw_A1_lo  = exp(int_ipcw + coef_ipcw - 1.96 * se_log_A1),
    IR_ipcw_A1_hi  = exp(int_ipcw + coef_ipcw + 1.96 * se_log_A1),
    IR_cf_A0       = ir_cf$IR[ir_cf$A == 0],
    IR_cf_A1       = ir_cf$IR[ir_cf$A == 1],
    IRR_cf         = ir_cf$IR[ir_cf$A == 1] / ir_cf$IR[ir_cf$A == 0],
    ptime_obs_A0   = pt0,
    ptime_obs_A1   = pt1
  )
}

# ============================================================================
# SECTION 6: BOOTSTRAP (patient-level) FOR IPCW RD, RMSTD, AND ABSOLUTE RISKS
# ============================================================================

# R1: eval_times passed in to avoid recomputing sort/unique on every bootstrap draw.
# R1: dat_pt_b kept as data.table (no tibble conversion) throughout hot path.
# R2: bootstrap iterations run in parallel via mclapply (Unix) or parLapply (Windows).
bootstrap_ipcw_rd_rmstd <- function(dat, dat_pt, scenario_id, tmax, n_boot = 20,
                                     eval_times = NULL) {
  ids   <- dat$id
  n     <- length(ids)
  pt_dt <- data.table::as.data.table(dat_pt)

  # R1: pre-compute eval_times once if not supplied
  if (is.null(eval_times)) {
    eval_times <- sort(unique(pt_dt[stop - start > 0 & start < tmax, start]))
  }

  # Single bootstrap draw — factored out so it can be called from lapply/parLapply
  one_boot <- function(b) {
    samp_ids <- sample(ids, n, replace = TRUE)

    dat_b    <- dat[match(samp_ids, ids), ]
    dat_b$id <- seq_len(n)

    id_dt    <- data.table::data.table(orig_id = samp_ids, new_id = seq_len(n))
    dat_pt_b <- pt_dt[id_dt, on = .(id = orig_id), allow.cartesian = TRUE, nomatch = 0L]
    dat_pt_b[, id := new_id][, new_id := NULL]
    # R1: no tibble conversion — pass data.table directly

    pt_w <- tryCatch(
      estimate_ipcw_weights_riskregression(dat_pt_b, dat_b, scenario_id, tmax,
                                           eval_times = eval_times),
      error = function(e) NULL
    )
    res_b <- if (!is.null(pt_w)) {
      tryCatch(calc_rd_rmstd_ipcw(pt_w, dat_b, tmax), error = function(e) NULL)
    } else NULL

    if (!is.null(res_b)) {
      list(rd = res_b$rd$rd, rmstd = res_b$rmstd$rmstd)
    } else {
      list(rd = NA_real_, rmstd = NA_real_)
    }
  }

  # R2: parallel bootstrap — mclapply on Unix (fork), parLapply on Windows (PSOCK).
  # Falls back to sequential if parallel is unavailable or n_boot is small.
  n_cores_boot <- if (.Platform$OS.type == "unix") {
    max(1L, parallel::detectCores() - 1L)
  } else {
    1L   # Windows: parallelise at the MC-replicate level instead (see monte_carlo_run_v10)
  }

  results <- if (n_cores_boot > 1L && n_boot >= 4L) {
    parallel::mclapply(seq_len(n_boot), one_boot,
                       mc.cores = n_cores_boot, mc.set.seed = TRUE)
  } else {
    lapply(seq_len(n_boot), one_boot)
  }

  rd_v    <- vapply(results, function(x) x$rd,    NA_real_)
  rmstd_v <- vapply(results, function(x) x$rmstd, NA_real_)

  pct <- c(0.025, 0.975)
  list(
    RD_ipcw_boot_lo    = unname(quantile(rd_v,    pct[1], na.rm = TRUE)),
    RD_ipcw_boot_hi    = unname(quantile(rd_v,    pct[2], na.rm = TRUE)),
    RMSTD_ipcw_boot_lo = unname(quantile(rmstd_v, pct[1], na.rm = TRUE)),
    RMSTD_ipcw_boot_hi = unname(quantile(rmstd_v, pct[2], na.rm = TRUE))
  )
}

# ============================================================================
# SECTION 7: MAIN ANALYSIS
# ============================================================================
# Allow others the option to use the package to calculate ipcw or my manual way. Results are the same.
analyze_coxph_pt <- function(dat_pt, dat, scenario_id, tmax, ipcw= "package", truncate = NULL, n_boot = 20) {
  # Compute IPCW weights once; shared by HR, RD, and RMSTD.
  if (ipcw == "package") {
    pt_weighted <- estimate_ipcw_weights_riskregression(
      dat_pt      = dat_pt,
      dat         = dat,
      scenario_id = scenario_id,
      tmax        = tmax,
      truncate    = truncate)
  } else if (ipcw == "manual") {
    pt_weighted <- estimate_ipcw_weights(
      dat         = dat_pt,
      scenario_id = scenario_id,
      tmax        = tmax,
      truncate    = truncate)
  } else {
    stop("ipcw must be either 'package' or 'manual'")
  }
  
  truth    <- counterfactual_truth(dat, tmax)
  naive    <- calc_naive_estimates(dat, tmax)
  hr_ipcw  <- calc_hr_ipcw(pt_weighted, tmax)
  ipcw_est <- calc_rd_rmstd_ipcw(pt_weighted, dat, tmax)

  # R1: pre-compute dat_pois ONCE here 
  #     and reuses the same object for eval_times_cache below.
  dat_pois <- pt_weighted %>%
    mutate(ptime = stop - start, y = as.integer(event_interval)) %>%
    filter(ptime > 0, start < tmax)

  # R1: pre-compute the time grid once; reused by bootstrap to skip sort/unique
  eval_times_cache <- sort(unique(dat_pois$start))

  irr_est  <- calc_irr(dat_pois, dat, tmax)

  # ---- Bootstrap CIs (RD and RMSTD only) ----
  boot_ci <- if (n_boot > 0L) {
    bootstrap_ipcw_rd_rmstd(dat, dat_pt, scenario_id, tmax, n_boot,
                             eval_times = eval_times_cache)
  } else {
    list(
      RD_ipcw_boot_lo    = NA_real_, RD_ipcw_boot_hi    = NA_real_,
      RMSTD_ipcw_boot_lo = NA_real_, RMSTD_ipcw_boot_hi = NA_real_
    )
  }

  # ---- Diagnostics from dat and pt_weighted ----
  diag_arm <- dat %>%
    group_by(A) %>%
    summarise(
      n_events = sum(event1),
      n_ltfu   = sum(censor_event),
      n_admin  = sum(admin_cens),
      n_total  = dplyr::n(),
      .groups  = "drop"
    )
  gd <- function(a, col) diag_arm[[col]][diag_arm$A == a]

  n_events_A0     <- gd(0, "n_events");  n_events_A1     <- gd(1, "n_events")
  n_ltfu_A0       <- gd(0, "n_ltfu");    n_ltfu_A1       <- gd(1, "n_ltfu")
  n_admin_A0      <- gd(0, "n_admin");   n_admin_A1      <- gd(1, "n_admin")
  pct_censored_A0 <- n_ltfu_A0 / gd(0, "n_total")
  pct_censored_A1 <- n_ltfu_A1 / gd(1, "n_total")

  w_ipcw_mean <- mean(pt_weighted$w_ipcw, na.rm = TRUE)
  w_ipcw_max  <- max(pt_weighted$w_ipcw,  na.rm = TRUE)
  w_ipcw_p95  <- unname(quantile(pt_weighted$w_ipcw, 0.95, na.rm = TRUE))

  scen <- paste0(scenario_id, "_pt")

  # ---- Formatted tibble (unchanged output for human review) ----
  fmt_hr    <- function(x) sprintf("%.2f (%.2f, %.2f)", x$hr, x$lo, x$hi)
  fmt_rd    <- function(x) sprintf("%.1f%% (%.1f%%, %.1f%%)", x$rd * 100, x$lo * 100, x$hi * 100)
  fmt_rmstd <- function(x) sprintf("%.2f (%.2f, %.2f)", x$rmstd, x$lo, x$hi)
  fmt_risk  <- function(x) sprintf("%.1f%%", x * 100)

  formatted <- tibble(
    scenario      = scen,
    HR_cf         = sprintf("%.2f", truth$hr),
    RD_cf         = sprintf("%.1f%%", truth$rd * 100),
    RMSTD_cf      = sprintf("%.2f", truth$rmstd),
    Risk_A0_cf    = fmt_risk(truth$risk_a0_cf),
    Risk_A1_cf    = fmt_risk(truth$risk_a1_cf),
    HR_naive      = fmt_hr(naive$hr),
    RD_naive      = fmt_rd(naive$rd),
    RMSTD_naive   = fmt_rmstd(naive$rmstd),
    Risk_A0_naive = fmt_risk(naive$risk_a0_naive),
    Risk_A1_naive = fmt_risk(naive$risk_a1_naive),
    HR_ipcw       = fmt_hr(hr_ipcw),
    RD_ipcw       = sprintf("%.1f%%", ipcw_est$rd$rd * 100),   # v11: point only (no Greenwood CI)
    RMSTD_ipcw    = sprintf("%.2f",  ipcw_est$rmstd$rmstd),    # v11: point only (no Greenwood CI)
    Risk_A0_ipcw  = fmt_risk(ipcw_est$risk_a0_ipcw),
    Risk_A1_ipcw  = fmt_risk(ipcw_est$risk_a1_ipcw),
    IRR_naive     = sprintf("%.2f (%.2f, %.2f)", irr_est$IRR_naive_est, irr_est$IRR_naive_lo, irr_est$IRR_naive_hi),
    IRR_ipcw      = sprintf("%.2f (%.2f, %.2f)", irr_est$IRR_ipcw_est,  irr_est$IRR_ipcw_lo,  irr_est$IRR_ipcw_hi)
  )

  # ---- Numeric tibble (machine-readable for Monte Carlo) ----
  numeric_out <- tibble(
    scenario              = scen,
    # Truth (scalars, no CI)
    HR_cf                 = truth$hr,
    RD_cf                 = truth$rd,
    RMSTD_cf              = truth$rmstd,
    Risk_A0_cf            = truth$risk_a0_cf,
    Risk_A1_cf            = truth$risk_a1_cf,
    # Naive estimates + analytical CIs
    HR_naive_est          = naive$hr$hr,
    HR_naive_lo           = naive$hr$lo,
    HR_naive_hi           = naive$hr$hi,
    RD_naive_est          = naive$rd$rd,
    RD_naive_lo           = naive$rd$lo,
    RD_naive_hi           = naive$rd$hi,
    RMSTD_naive_est       = naive$rmstd$rmstd,
    RMSTD_naive_lo        = naive$rmstd$lo,
    RMSTD_naive_hi        = naive$rmstd$hi,
    Risk_A0_naive         = naive$risk_a0_naive,
    Risk_A0_naive_lo      = naive$risk_a0_naive_lo,
    Risk_A0_naive_hi      = naive$risk_a0_naive_hi,
    Risk_A1_naive         = naive$risk_a1_naive,
    Risk_A1_naive_lo      = naive$risk_a1_naive_lo,
    Risk_A1_naive_hi      = naive$risk_a1_naive_hi,
    # IPCW estimates + analytical CIs
    HR_ipcw_est           = hr_ipcw$hr,
    HR_ipcw_lo            = hr_ipcw$lo,
    HR_ipcw_hi            = hr_ipcw$hi,
    RD_ipcw_est           = ipcw_est$rd$rd,
    RD_ipcw_lo            = ipcw_est$rd$lo,
    RD_ipcw_hi            = ipcw_est$rd$hi,
    RMSTD_ipcw_est        = ipcw_est$rmstd$rmstd,
    RMSTD_ipcw_lo         = ipcw_est$rmstd$lo,
    RMSTD_ipcw_hi         = ipcw_est$rmstd$hi,
    Risk_A0_ipcw          = ipcw_est$risk_a0_ipcw,
    Risk_A0_ipcw_lo       = ipcw_est$risk_a0_ipcw_lo,
    Risk_A0_ipcw_hi       = ipcw_est$risk_a0_ipcw_hi,
    Risk_A1_ipcw          = ipcw_est$risk_a1_ipcw,
    Risk_A1_ipcw_lo       = ipcw_est$risk_a1_ipcw_lo,
    Risk_A1_ipcw_hi       = ipcw_est$risk_a1_ipcw_hi,
    # IPCW bootstrap CIs (RD and RMSTD only)
    RD_ipcw_boot_lo       = boot_ci$RD_ipcw_boot_lo,
    RD_ipcw_boot_hi       = boot_ci$RD_ipcw_boot_hi,
    RMSTD_ipcw_boot_lo    = boot_ci$RMSTD_ipcw_boot_lo,
    RMSTD_ipcw_boot_hi    = boot_ci$RMSTD_ipcw_boot_hi,
    # Diagnostics
    n_events_A0           = n_events_A0,
    n_events_A1           = n_events_A1,
    n_ltfu_A0             = n_ltfu_A0,
    n_ltfu_A1             = n_ltfu_A1,
    n_admin_A0            = n_admin_A0,
    n_admin_A1            = n_admin_A1,
    pct_censored_A0       = pct_censored_A0,
    pct_censored_A1       = pct_censored_A1,
    w_ipcw_mean           = w_ipcw_mean,
    w_ipcw_max            = w_ipcw_max,
    w_ipcw_p95            = w_ipcw_p95,
    # IRR truth (counterfactual, no CI)
    IRR_cf                = irr_est$IRR_cf,
    IR_cf_A0              = irr_est$IR_cf_A0,
    IR_cf_A1              = irr_est$IR_cf_A1,
    # IRR naive + analytical sandwich CI
    IRR_naive_est         = irr_est$IRR_naive_est,
    IRR_naive_lo          = irr_est$IRR_naive_lo,
    IRR_naive_hi          = irr_est$IRR_naive_hi,
    # IRR IPCW + analytical sandwich CI
    IRR_ipcw_est          = irr_est$IRR_ipcw_est,
    IRR_ipcw_lo           = irr_est$IRR_ipcw_lo,
    IRR_ipcw_hi           = irr_est$IRR_ipcw_hi,
    # Arm-specific IRs + analytical CIs
    IR_obs_A0             = irr_est$IR_obs_A0,
    IR_obs_A0_lo          = irr_est$IR_obs_A0_lo,
    IR_obs_A0_hi          = irr_est$IR_obs_A0_hi,
    IR_obs_A1             = irr_est$IR_obs_A1,
    IR_obs_A1_lo          = irr_est$IR_obs_A1_lo,
    IR_obs_A1_hi          = irr_est$IR_obs_A1_hi,
    IR_ipcw_A0            = irr_est$IR_ipcw_A0,
    IR_ipcw_A0_lo         = irr_est$IR_ipcw_A0_lo,
    IR_ipcw_A0_hi         = irr_est$IR_ipcw_A0_hi,
    IR_ipcw_A1            = irr_est$IR_ipcw_A1,
    IR_ipcw_A1_lo         = irr_est$IR_ipcw_A1_lo,
    IR_ipcw_A1_hi         = irr_est$IR_ipcw_A1_hi,
    # Observed person-time by arm
    ptime_obs_A0          = irr_est$ptime_obs_A0,
    ptime_obs_A1          = irr_est$ptime_obs_A1
  )

  list(formatted = formatted, numeric = numeric_out)
}

# ============================================================================
# SECTION 8: DESCRIPTIVE SUMMARY
# ============================================================================

# dat: person-level data from simulate_scenario*()$dat (one row per person).
summarize_followup <- function(dat) {
  summarise_by <- function(group_var) {
    dat %>%
      group_by(A, !!sym(group_var)) %>%
      summarise(
        pts       = n_distinct(id),
        events    = sum(event1),
        py        = sum(time),
        median_fu = median(time),
        IR        = events / py,
        IR_lower  = ifelse(events == 0, 0, qchisq(0.025, 2 * events) / (2 * py)),
        IR_upper  = qchisq(0.975, 2 * (events + 1)) / (2 * py),
        IR_CI     = sprintf("%.3f (%.3f, %.3f)", IR, IR_lower, IR_upper),
        .groups   = "drop"
      ) %>%
      select(A, !!sym(group_var), pts, median_fu, IR_CI)
  }

  result <- list()
  if ("L_c1" %in% names(dat)) result$by_l_c1 <- summarise_by("L_c1")
  if ("L_c2" %in% names(dat)) result$by_l_c2 <- summarise_by("L_c2")
  result
}


