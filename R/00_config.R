# ============================================================================
#  00_config.R  —  Global configuration (paths, seed, model specifications)
#
#  Sourced first by every other script. Nothing here reads data or fits models;
#  it only defines WHERE things live and WHICH models we run.
# ============================================================================

# ---- Packages --------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)        # project-root-relative paths (no absolute paths anywhere)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(posterior)
  library(cmdstanr)
  library(loo)
})

# ---- Paths (all derived from the project root) -----------------------------
BASE_DIR   <- here::here()
DATA_DIR   <- file.path(BASE_DIR, "data")
STAN_DIR   <- file.path(BASE_DIR, "stan")
OUT_DIR    <- file.path(BASE_DIR, "outputs")
FIT_DIR    <- file.path(OUT_DIR, "fits")
dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)

STAN_FILE  <- file.path(STAN_DIR, "model_traitNCP.stan")

# ---- Reproducibility -------------------------------------------------------
SEED <- 20260722L

# ---- Sampler settings (one place, reused everywhere) -----------------------
MCMC <- list(
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1500,
  iter_sampling   = 1500,
  adapt_delta     = 0.97,
  max_treedepth   = 12,
  refresh         = 500
)

# ---- Model specifications --------------------------------------------------
#  Every model in the paper is one row here: a name + the two Stan flags.
#  trait_structure: 0 pair | 1 service | 2 trait+service | 3 organ+trait+service
#  quality_mode   : 0 none | 1 likelihood-weight | 2 precision (quality->tau_j)
#
#  PRIMARY model for the paper = "qp"  (full structure, quality as precision).
MODEL_SPECS <- tribble(
  ~name,     ~label,                              ~trait_structure, ~quality_mode,
  "simple",  "Pair only",                          0L,               0L,
  "ncp",     "NCP only + quality-weight",          1L,               1L,
  "full",    "Trait+NCP + quality-weight",         2L,               1L,
  "nq",      "Trait+NCP, no quality",              2L,               0L,
  "qp",      "Trait+NCP, quality-precision (PRIMARY)", 2L,           2L,
  "qp_organ","Organ+Trait+NCP, quality-precision", 3L,               2L
)

PRIMARY_MODEL <- "qp"

# ---- Quick smoke-test mode -------------------------------------------------
#  Set env TRAITNCP_QUICK=1 (run_all.R does this when QUICK=TRUE) for a fast
#  end-to-end check: tiny sampler + primary model only. Re-applied on every
#  source() of this file, so it survives the per-script re-sourcing.
if (Sys.getenv("TRAITNCP_QUICK") == "1") {
  MCMC$iter_warmup   <- 300
  MCMC$iter_sampling <- 300
  MODEL_SPECS <- MODEL_SPECS %>% filter(name == PRIMARY_MODEL)
  message("QUICK mode ON: 300+300 iterations, primary model only.")
}

message("Config loaded. Primary model = '", PRIMARY_MODEL,
        "'. ", nrow(MODEL_SPECS), " model specification(s).")
