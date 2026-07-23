# ============================================================================
#  run_lowrank.R  —  Fit and evaluate the reduced-rank (K=2) trait x NCP model
#
#  BOUNDED PROTOTYPE. Standalone: reads outputs/prepared_data.rds and the qp /
#  simple fits (for LOO comparison). Refits nothing else.
#     Rscript -e 'source("R/run_lowrank.R")'
#
#  Decision rule:
#   * clean convergence (0 divergences, Rhat < 1.01) AND better LOO than qp/simple
#       -> promote as the single main model (variance + trade-off directions).
#   * poor convergence or no LOO gain -> keep the two-model design; this stays a
#       supplementary robustness check.
#
#  OUTPUTS  outputs/fits/fit_lowrank.rds, outputs/lowrank_loadings.csv,
#           outputs/pair_results_lowrank.csv
# ============================================================================

source(here::here("R", "00_config.R"))
prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))

sd <- prep$stan_data
sd$K <- 2L                                  # latent rank
sd$quality_mode <- 2L; sd$link <- 0L; sd$use_dm <- 1L; sd$estimate_kappa <- 1L

mod <- cmdstan_model(file.path(STAN_DIR, "model_lowrank.stan"))
fit <- mod$sample(
  data = sd, seed = SEED, chains = 4, parallel_chains = 4,
  iter_warmup = 2000, iter_sampling = 2000,
  adapt_delta = 0.99, max_treedepth = 12, refresh = 500, init = 0.3)
fit$save_object(file.path(FIT_DIR, "fit_lowrank.rds"))

# ---- 1. Convergence ---------------------------------------------------------
cat("\n=== CONVERGENCE (low-rank K=2) ===\n")
cat("Divergences:", sum(fit$diagnostic_summary()$num_divergent),
    "| max Rhat:", round(max(fit$summary()$rhat, na.rm = TRUE), 3),
    "| min ESS:", round(min(fit$summary()$ess_bulk, na.rm = TRUE)), "\n")
print(fit$summary(c("sigma_alpha", "tau", "gamma", "kappa", "phi")))

# ---- 2. LOO vs qp and simple ------------------------------------------------
loo_of <- function(f) loo(readRDS(f)$draws("log_lik", format = "matrix"), cores = 4)
loo_lr     <- loo(fit$draws("log_lik", format = "matrix"), cores = 4)
loo_qp     <- loo_of(file.path(FIT_DIR, "fit_qp.rds"))
loo_simple <- loo_of(file.path(FIT_DIR, "fit_simple.rds"))
cat("\n=== MODEL COMPARISON (higher ELPD = better) ===\n")
print(loo_compare(list(lowrank = loo_lr, qp = loo_qp, simple = loo_simple)))

# ---- 3. Latent NCP loadings F (the interpretable syndrome axes) -------------
service_levels <- sort(unique(prep$data_model$ncp))
trait_levels   <- sort(unique(prep$data_model$trait_grouped))
Fsum <- fit$summary("F", "mean")
Fmat <- matrix(Fsum$mean, nrow = length(service_levels), byrow = FALSE)
rownames(Fmat) <- service_levels; colnames(Fmat) <- paste0("factor", seq_len(sd$K))
cat("\n=== NCP loadings on the latent axes (F) ===\n")
print(round(Fmat, 2))

Lsum <- fit$summary("L", "mean")
Lmat <- matrix(Lsum$mean, nrow = length(trait_levels), byrow = FALSE)
rownames(Lmat) <- trait_levels; colnames(Lmat) <- paste0("factor", seq_len(sd$K))
for (f in seq_len(sd$K)) {
  cat("\nTop traits on factor", f, "(+ and -):\n")
  o <- order(Lmat[, f])
  print(round(head(Lmat[o, f], 5), 2)); print(round(tail(Lmat[o, f], 5), 2))
}
write.csv(cbind(component = "NCP", label = service_levels, Fmat),
          file.path(OUT_DIR, "lowrank_loadings.csv"), row.names = FALSE)

# ---- 4. Pair-level directions from the low-rank model -----------------------
P_lr <- apply(fit$draws("alpha", format = "matrix"), 2, function(x) mean(x > 0))
pair <- prep$data_model %>% distinct(ts_id, trait_grouped, ncp) %>% arrange(ts_id) %>%
  mutate(prob_positive = P_lr[ts_id])
write.csv(pair, file.path(OUT_DIR, "pair_results_lowrank.csv"), row.names = FALSE)
cat("\nSaved pair_results_lowrank.csv, lowrank_loadings.csv, fit_lowrank.rds\n")
cat("\nRead the NCP loadings: if factor 1 separates provisioning (+) from",
    "regulation (-),\nthe model has recovered the trade-off axis you saw by hand.\n")
