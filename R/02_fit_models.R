# ============================================================================
#  02_fit_models.R  —  Fit every model specification from MODEL_SPECS
#
#  Compiles the single Stan model once, then fits each specification by only
#  changing the two flags (trait_structure, quality_mode). Saves one fit per
#  model to outputs/fits/fit_<name>.rds.
#
#  Run 01_prepare_data.R first (or it is sourced automatically below).
# ============================================================================

source(here::here("R", "00_config.R"))

prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
stan_data <- prep$stan_data

# ---- Compile once ----------------------------------------------------------
mod <- cmdstan_model(STAN_FILE)

# Resume-safe: skip a model whose fit is already on disk. Force a full refit
# from the terminal with  FORCE_REFIT=TRUE Rscript run_all.R
FORCE_REFIT <- toupper(Sys.getenv("FORCE_REFIT", "FALSE")) %in%
  c("TRUE", "1", "YES", "T")

# ---- Fit each specification ------------------------------------------------
fit_one <- function(spec) {
  fit_path <- file.path(FIT_DIR, paste0("fit_", spec$name, ".rds"))
  if (file.exists(fit_path) && !FORCE_REFIT) {
    message("Skipping '", spec$name, "' — fit already on disk ",
            "(FORCE_REFIT=TRUE to redo).")
    return(invisible(readRDS(fit_path)))
  }
  message("\n=== Fitting '", spec$name, "'  (", spec$label, ") ===")
  sd <- stan_data
  sd$trait_structure <- spec$trait_structure
  sd$quality_mode    <- spec$quality_mode

  fit <- mod$sample(
    data            = sd,
    seed            = SEED,
    chains          = MCMC$chains,
    parallel_chains = MCMC$parallel_chains,
    iter_warmup     = MCMC$iter_warmup,
    iter_sampling   = MCMC$iter_sampling,
    adapt_delta     = MCMC$adapt_delta,
    max_treedepth   = MCMC$max_treedepth,
    refresh         = MCMC$refresh
  )

  # Persist as a portable .rds (draws + diagnostics travel with it)
  fit$save_object(file.path(FIT_DIR, paste0("fit_", spec$name, ".rds")))

  # Quick convergence readout
  key <- c("sigma_alpha", "sigma_trait", "sigma_service", "tau", "gamma",
           "phi", "c[1]", "c[2]")
  key <- key[key %in% fit$metadata()$stan_variables |
               key %in% dimnames(fit$draws())$variable]
  print(fit$summary(key))
  div <- sum(fit$diagnostic_summary()$num_divergent)
  message("Divergences: ", div,
          " | max Rhat: ",
          round(max(fit$summary()$rhat, na.rm = TRUE), 3))
  invisible(fit)
}

for (i in seq_len(nrow(MODEL_SPECS))) fit_one(MODEL_SPECS[i, ])

message("\nAll models fitted. Fits in outputs/fits/.")
