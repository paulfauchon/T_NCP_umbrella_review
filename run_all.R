# ============================================================================
#  run_all.R  —  One-click, self-contained run of the whole pipeline.
#
#  HOW TO RUN
#   1. Open  clean_pipeline.Rproj  in RStudio (this sets the working directory
#      and anchors all paths; a `.here` file also anchors it for plain Rscript).
#   2. source("run_all.R")
#
#  The two switches below control scope. Leave QUICK = TRUE for a first check
#  (a few minutes); set FALSE for the full analysis. RUN_HEAVY adds the
#  recovery simulation + full sensitivity battery (can take hours).
# ============================================================================

# Switches — defaults here, overridable from the terminal via env vars, e.g.
#   QUICK=FALSE RUN_HEAVY=TRUE Rscript run_all.R
envflag <- function(name, default)
  toupper(Sys.getenv(name, unset = default)) %in% c("TRUE", "1", "YES", "T")

QUICK       <- envflag("QUICK",       "TRUE")   # TRUE: fast smoke test (primary model, 300 iters)
RUN_HEAVY   <- envflag("RUN_HEAVY",   "FALSE")  # TRUE: also run 06 recovery + 08 sensitivity
BUILD_SIGMA <- envflag("BUILD_SIGMA", "FALSE")  # TRUE: rebuild data/sigma.csv from the .xlsx

Sys.setenv(TRAITNCP_QUICK = if (QUICK) "1" else "0")
message("Run mode: QUICK=", QUICK, " RUN_HEAVY=", RUN_HEAVY, " BUILD_SIGMA=", BUILD_SIGMA)

# ---- 0. One-time dependency setup ------------------------------------------
setup <- function() {
  repos <- c("https://cloud.r-project.org")
  pkgs  <- c("here", "dplyr", "tidyr", "stringr", "tibble", "posterior", "loo",
             "ordinal", "ggplot2", "lme4", "svglite")   # last three: figures
  if (BUILD_SIGMA) pkgs <- c(pkgs, "readxl")
  miss <- setdiff(pkgs, rownames(installed.packages()))
  if (length(miss)) install.packages(miss, repos = repos)

  # cmdstanr (not on CRAN) + the CmdStan toolchain
  if (!requireNamespace("cmdstanr", quietly = TRUE))
    install.packages("cmdstanr",
                     repos = c("https://mc-stan.org/r-packages/", repos))
  ok <- tryCatch(!is.null(cmdstanr::cmdstan_path()), error = function(e) FALSE)
  if (!ok) {
    message("Installing CmdStan (one-time, a few minutes)...")
    cmdstanr::install_cmdstan()
  }
}
setup()

# ---- 1. Pipeline in order --------------------------------------------------
run <- function(f) { message("\n=== ", f, " ==="); source(here::here("R", f)) }

if (BUILD_SIGMA) run("00b_build_sigma.R")   # else use the shipped data/sigma.csv
run("01_prepare_data.R")
run("02_fit_models.R")
run("03_diagnostics.R")
run("04_variance_decomposition.R")
run("05_grouping_sensitivity.R")            # skips models not fitted in QUICK mode
run("07_results_pairs.R")
run("09_frequentist_check.R")               # framework-independent check (no MCMC)

if (RUN_HEAVY) {
  run("08_sensitivity_battery.R")
  run("06_recovery_simulation.R")
}

# ---- 2a. Section-2 figure (review disagreement) ----------------------------
# Needs only data/data_model.csv, so it runs in every mode.
fig2 <- here::here("figures", "fig2_review_disagreement.R")
if (file.exists(fig2)) {
  message("\n=== figures/fig2_review_disagreement.R ===")
  tryCatch(source(fig2),
           error = function(e)
             message("  (figure step failed: ", conditionMessage(e), ")"))
}

# ---- 2a-bis. Section-3 figure (review methodology signature + GLMM) --------
# Needs data/data_model.csv + data/meta-data2.csv; fits a study-level GLMM
# (requires lme4) and writes outputs/fig3_glmm.csv. Runs in every mode.
fig3 <- here::here("figures", "fig3_review_methodology.R")
if (file.exists(fig3)) {
  message("\n=== figures/fig3_review_methodology.R ===")
  tryCatch(source(fig3),
           error = function(e)
             message("  (figure step failed: ", conditionMessage(e),
                     " — needs lme4; install.packages('lme4').)"))
}

# ---- 2b. Main figure (pair-level trade-off matrix) -------------------------
# Needs the pair-only 'simple' fit: check_reversals.R writes
# outputs/reversal_check.csv, then figures/fig_main_tradeoff.R renders the
# main (13-trait) + supplementary (19-trait) figures into outputs/.
# (A stdlib Python twin, figures/fig_main_tradeoff.py, produces the same figure
#  if you prefer to run it outside R.)
# Skipped automatically in QUICK mode (only the primary model is fitted there).
if (!QUICK && file.exists(here::here("outputs", "fits", "fit_simple.rds"))) {
  run("check_reversals.R")
  fig <- here::here("figures", "fig_main_tradeoff.R")
  if (file.exists(fig)) {
    message("\n=== figures/fig_main_tradeoff.R ===")
    tryCatch(source(fig),
             error = function(e)
               message("  (figure step failed: ", conditionMessage(e),
                       " — run figures/fig_main_tradeoff.R manually.)"))
  }
}

message("\nDONE. Outputs in outputs/ (fits/, *.csv, fig_*_tradeoff.svg/png). ",
        if (QUICK) "QUICK run — set QUICK<-FALSE for the real analysis."
        else "Full run complete.")

