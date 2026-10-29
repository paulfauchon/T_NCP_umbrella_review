# ============================================================================
#  03_diagnostics.R  —  Convergence, posterior predictive checks, model ranking
#
#  For each fitted model: R-hat / ESS / divergences, a 3-category posterior
#  predictive check (observed vs. replicated proportions), and PSIS-LOO.
#  Produces a single model-comparison table (outputs/model_comparison.csv).
#
#  This is where the "5 models converge" robustness argument is quantified.
# ============================================================================

source(here::here("R", "00_config.R"))
prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
stan_data <- prep$stan_data

load_fit <- function(name)
  readRDS(file.path(FIT_DIR, paste0("fit_", name, ".rds")))

# ---- 1. Convergence summary per model --------------------------------------
convergence_row <- function(name) {
  fit <- load_fit(name)
  s   <- fit$summary()
  tibble(
    model        = name,
    max_rhat     = round(max(s$rhat,     na.rm = TRUE), 3),
    min_ess_bulk = round(min(s$ess_bulk, na.rm = TRUE)),
    divergences  = sum(fit$diagnostic_summary()$num_divergent)
  )
}
convergence <- bind_rows(lapply(MODEL_SPECS$name, convergence_row))
print(convergence)

# ---- 2. Posterior predictive check (aggregate direction proportions) -------
ppc_row <- function(name) {
  fit    <- load_fit(name)
  N      <- stan_data$N
  y_obs  <- stan_data$counts
  p_obs  <- colSums(y_obs) / sum(y_obs)

  yr <- fit$draws("y_rep", format = "matrix")
  prop_rep <- sapply(1:3, function(d) {
    cols <- paste0("y_rep[", seq_len(N), ",", d, "]")
    rowSums(yr[, cols])
  })
  prop_rep <- prop_rep / rowSums(prop_rep)
  tibble(
    model = name,
    ppc_neg  = sprintf("%.3f obs / %.3f rep", p_obs[1], median(prop_rep[, 1])),
    ppc_null = sprintf("%.3f obs / %.3f rep", p_obs[2], median(prop_rep[, 2])),
    ppc_pos  = sprintf("%.3f obs / %.3f rep", p_obs[3], median(prop_rep[, 3]))
  )
}
ppc <- bind_rows(lapply(MODEL_SPECS$name, ppc_row))
print(ppc)

# ---- 3. PSIS-LOO and model ranking -----------------------------------------
#  log_lik is UNWEIGHTED by construction (see the Stan generated quantities),
#  so LOO compares fit, not data quality.
loo_list <- setNames(
  lapply(MODEL_SPECS$name, function(nm)
    loo(load_fit(nm)$draws("log_lik", format = "matrix"), cores = 4)),
  MODEL_SPECS$name
)

if (length(loo_list) >= 2) {
  loo_tab <- as.data.frame(loo_compare(loo_list))
  loo_tab$model <- rownames(loo_tab)
  loo_tab <- loo_tab[, c("model", "elpd_diff", "se_diff",
                         "elpd_loo", "se_elpd_loo")]
} else {                                   # single model (e.g. QUICK mode)
  e <- loo_list[[1]]$estimates
  loo_tab <- data.frame(model = names(loo_list)[1],
                        elpd_diff = 0, se_diff = 0,
                        elpd_loo = e["elpd_loo", "Estimate"],
                        se_elpd_loo = e["elpd_loo", "SE"])
}
print(loo_tab, row.names = FALSE)

# ---- 4. Save unified comparison table --------------------------------------
comparison <- convergence %>%
  left_join(loo_tab, by = "model") %>%
  left_join(ppc, by = "model") %>%
  left_join(MODEL_SPECS %>% select(name, label), by = c("model" = "name"))

write.csv(comparison, file.path(OUT_DIR, "model_comparison.csv"),
          row.names = FALSE)
message("Saved outputs/model_comparison.csv")

