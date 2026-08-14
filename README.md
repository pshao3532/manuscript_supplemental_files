# Censoring Bias Simulation — Code Companion

Code accompanying **"Differential Censoring Is Not Informative Censoring: A
Simulation Study of Censoring Bias in Survival Analyses."**

Two scripts are included:

| File | Role |
|---|---|
| [`simulation_ipcw_github.R`](simulation_ipcw_github.R) | Defines everything: the 4 data-generating scenarios, the naive and IPCW estimators, the counterfactual truth, and the bootstrap. Sourced by the other script — not meant to be run on its own. |
| [`monte_carlo_run_github.R`](monte_carlo_run_github.R) | The entry point. Sources the file above, runs `R` Monte Carlo replicates across all 4 scenarios in parallel, and writes out bias / coverage / MSE summary tables. |

You only need to run `monte_carlo_run_github.R`.

---

## 1. Requirements

- R ≥ 4.3 (developed under R 4.5.1)
- Packages:

```r
install.packages(c(
  "tidyverse", "survival", "riskRegression",
  "survRM2", "data.table", "sandwich"
))
# "parallel" ships with base R -- no install needed
```

## 2. Folder layout

Both scripts live at the repository root and assume they're run with the
**repository root as the working directory**:

```
manuscript_supplemental_files/
├── simulation_ipcw_github.R
├── monte_carlo_run_github.R
└── results/            # created automatically on first run
```

If you move the scripts into a subfolder, update the `source(...)` call
at the top of `monte_carlo_run_github.R` to match.

## 3. Quick start

```r
# from an R session with the working directory set to the repo root
# (e.g. setwd() to where you cloned this repo, or open it as an RStudio project):
source("monte_carlo_run_github.R")
```

Or from a terminal:

```bash
Rscript monte_carlo_run_github.R
```

The default settings run quickly (a few minutes on a laptop) so you can
confirm everything works before scaling up to the manuscript's full
settings (see §6).

## 4. The four scenarios

All four scenarios share the same exposure-generating model
(`L ~ Bernoulli(0.5)`, `A | L ~ Bernoulli(plogis(-0.5 + 1.2*L))`) and a fixed
administrative end of follow-up (`max_fu`). They differ only in how the
loss-to-follow-up (LTFU) censoring hazard is generated:

| Scenario | Name | Censoring hazard depends on... | Interpretation |
|---|---|---|---|
| S1 | `S1_noninf_nondiff` | nothing (constant) | Non-informative, non-differential |
| S2 | `S2_noninf_diff_nausea` | `L_c1`, `L_c2` (whose *prevalence* depends on arm `A`), but not on event risk | Non-informative, differential |
| S3 | `S3_inf_nondiff` | `L_c1`, `L_c2` (which also drive event risk), not on `A` directly | Informative, non-differential |
| S4 | `S4_inf_diff` | `L_c1`, `L_c2`, `A`, and an `A:L_c1` interaction | Informative, differential |

"Informative" = the censoring hazard shares risk factors with the event
hazard (censoring is related to prognosis). "Differential" = the censoring
hazard depends directly on treatment arm `A`. All event and censoring times
are generated as independent exponentials (given covariates); administrative
censoring is applied at `max_fu`.

Each replicate produces, per scenario:

- **Counterfactual truth** (`_cf` columns) — estimated from `T_counterfactual`
  with censoring turned off (only administrative censoring at `max_fu`
  remains).
- **Naive estimates** (`_naive` columns) — standard unweighted Cox/KM/RMST
  on the observed (censored) data.
- **IPCW estimates** (`_ipcw` columns) — the same estimands, but using
  inverse-probability-of-censoring weights fit with the *correctly specified*
  censoring model for that scenario (see `censor_surv_formula()` in
  `simulation_ipcw_github.R`).

Estimands reported: hazard ratio (HR), risk difference (RD), restricted mean
survival time difference (RMSTD), arm-specific risk at `tmax`, incidence
rate ratio (IRR), and arm-specific incidence rate (IR).

## 5. User-controlled settings

At the top of `monte_carlo_run_github.R`:

```r
R         <- 200        # number of Monte Carlo replicates
tmax      <- 24         # analysis horizon
max_fu    <- 24         # administrative end of follow-up
base_seed <- 20260606   # RNG seed (reproducibility)
n_boot    <- 200        # bootstrap iterations per replicate (0 = skip)

n_s1 <- 2000             # sample size per scenario, per replicate
n_s2 <- 2000
n_s3 <- 5000
n_s4 <- 5000

n_cores <- max(1L, parallel::detectCores() - 1L)
```

## 6. Reproducing the manuscript's numbers

The manuscript results used `R = 1000` replicates and `n_boot = 500`. This
is substantially slower (expect a few hours depending on cores/RAM) —
start with the defaults above to confirm your setup works, then increase
`R` and `n_boot` for a full run. If you're RAM-constrained, lower `n_cores`.

## 7. Output

Written to `results/`:

| File | Contents |
|---|---|
| `mc_results_raw.csv` / `.rds` | One row per replicate × scenario: every point estimate, CI, and diagnostic (event/censoring counts, IPCW weight summary) produced by that replicate. |
| `mc_summary_performance_long.csv` | One row per scenario × metric: mean estimate, bias, relative bias, empirical SE, MSE, CI coverage, and mean CI width, computed across all replicates. |
| `mc_summary_performance.csv` | Same as above, reformatted as a compact table with `"estimate (2.5%, 97.5%)"` strings — closest to what's reported in the manuscript. |
| `mc_weight_diagnostics.csv` | Per-scenario average IPCW weight (mean/max/95th percentile) and average censoring/event counts by arm — useful for spotting extreme-weight problems. |

**Coverage/MSE note**: these are computed against a single fixed benchmark —
the *mean* of each replicate's own counterfactual truth — rather than each
replicate's individually correlated truth, which would otherwise inflate
apparent coverage. See the comment above `calculate_performance()` in
`monte_carlo_run_github.R`.

## 8. A note on parallelism

`monte_carlo_run_github.R` parallelizes **across replicates** using a PSOCK
cluster (portable across Windows/Mac/Linux). On Unix, the bootstrap inside
each replicate (`bootstrap_ipcw_rd_rmstd()`, in `simulation_ipcw_github.R`)
*also* parallelizes internally via `mclapply` when it detects it's on Unix.
Running both at once can oversubscribe your cores (not incorrect, just
slower than expected) — if that happens, lower `n_cores`.
