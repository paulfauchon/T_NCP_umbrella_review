# ============================================================================
#  08_sensitivity_battery.R  —  Is the headline a property of the data or of
#                               the model? (the reviewer-proofing script)
#
#  Re-fits the primary model under a battery of perturbations and reports the
#  genuine between-review heterogeneity share each time. A conclusion that
#  survives is robust; one that swings is an assumption, not a result.
#
#  Battery (each answers a specific reviewer objection):
#   - kappa fixed {0.5, 1, 2}          : Sigma-scale sensitivity
#   - kappa estimated                  : + posterior corr(tau, kappa, gamma)
#   - Dirichlet-multinomial vs multinomial : does overdispersion phi matter?
#   - logit vs probit link             : link-family sensitivity
#   - wide priors on sigma/tau/gamma   : prior sensitivity
#   - Sigma diagonal +/-20%            : Sigma-is-estimated (not known) sensitivity
#  Plus prior-vs-posterior contraction for tau/gamma/kappa/phi and the count of
#  records that can inform phi (n>1).
#
#  OUTPUT  outputs/sensitivity_battery.csv, outputs/param_correlations.csv
# ============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "_helpers.R"))

prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
base <- prep$stan_data
spec <- MODEL_SPECS %>% filter(name == PRIMARY_MODEL)
base$trait_structure <- spec$trait_structure
base$quality_mode    <- spec$quality_mode
use_cat <- spec$trait_structure >= 3
mod  <- cmdstan_model(STAN_FILE)

# Refit with a set of overrides; return heterogeneity share + diagnostics
run_variant <- function(name, overrides = list(), Sigma_used = prep$Sigma) {
  sd <- modifyList(base, overrides)
  sd$Sigma <- Sigma_used
  fit <- mod$sample(data = sd, seed = SEED, chains = 4, parallel_chains = 4,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.97, refresh = 0)
  vd  <- variance_decomposition(fit, Sigma_used, use_cat = use_cat)
  het_full   <- vd %>% filter(component == "review_heterogeneity", basis == "D_full")
  het_struct <- vd %>% filter(component == "review_heterogeneity", basis == "D_struct")
  tibble(
    variant           = name,
    het_share_Dfull   = het_full$median,
    het_Dfull_lo      = het_full$q05,
    het_Dfull_hi      = het_full$q95,
    het_share_Dstruct = het_struct$median,
    max_rhat          = round(max(fit$summary()$rhat, na.rm = TRUE), 3),
    divergences       = sum(fit$diagnostic_summary()$num_divergent)
  )
}

# ---- Sigma perturbations (diagonal +/-20%: Sigma is estimated, not known) ---
Sigma_hi <- prep$Sigma; diag(Sigma_hi) <- diag(Sigma_hi) * 1.2
Sigma_lo <- prep$Sigma; diag(Sigma_lo) <- diag(Sigma_lo) * 0.8

battery <- bind_rows(
  run_variant("primary (kappa estimated, DM, logit)"),
  run_variant("kappa fixed = 0.5", list(estimate_kappa = 0L, kappa_fixed = 0.5)),
  run_variant("kappa fixed = 1",   list(estimate_kappa = 0L, kappa_fixed = 1.0)),
  run_variant("kappa fixed = 2",   list(estimate_kappa = 0L, kappa_fixed = 2.0)),
  run_variant("multinomial (no phi)", list(use_dm = 0L)),
  run_variant("probit link",      list(link = 1L)),
  run_variant("wide priors",      list(ps_sigma = 3, ps_tau = 3, ps_gamma = 2)),
  run_variant("Sigma diag +20%",  Sigma_used = Sigma_hi),
  run_variant("Sigma diag -20%",  Sigma_used = Sigma_lo),
  # Sigma ~ 0: drop the overlap structure entirely (tiny PD matrix so the
  # cov_matrix constraint still holds). Does review dominance survive with no
  # dependence structure at all? (reviewer request)
  run_variant("Sigma = 0 (no overlap)",
              Sigma_used = diag(1e-8, nrow(prep$Sigma)))
)
write.csv(battery, file.path(OUT_DIR, "sensitivity_battery.csv"), row.names = FALSE)

cat("\n=============== SENSITIVITY BATTERY ===============\n")
print(as.data.frame(battery %>%
        mutate(across(starts_with("het"), ~sprintf("%.1f%%", 100 * .)))),
      row.names = FALSE)
cat("\nIf het_share_Dfull is stable across ALL rows, 'variance is mostly at the\n",
    "review level' is a data property. If it swings, report it as conditional.\n", sep = "")

# ---- Posterior correlations corr(tau, kappa, gamma) from the kappa-free fit --
fit_kf <- mod$sample(
  data = modifyList(base, list(estimate_kappa = 1L)),
  seed = SEED, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.99, refresh = 0)

corr_draws <- fit_kf$draws(c("tau", "kappa", "gamma"), format = "draws_matrix")
cor_mat <- cor(as.matrix(corr_draws))
write.csv(round(cor_mat, 3), file.path(OUT_DIR, "param_correlations.csv"))

cat("\n--- Posterior correlations (kappa-free fit) ---\n")
print(round(cor_mat, 3))
cat("Correlations near +/-1 => the split is a ridge (not identified by the data),\n",
    "which is exactly why kappa is FIXED to 1 in the primary model.\n", sep = "")

# ---- Prior vs posterior contraction + phi information -----------------------
n_draw <- 4000
prior_draws <- list(
  tau   = abs(rnorm(n_draw, 0, base$ps_tau)),
  gamma = abs(rnorm(n_draw, 0, base$ps_gamma)),
  kappa = rlnorm(n_draw, 0, base$ps_kappa_ln),
  phi   = exp(rnorm(n_draw, log(10), 1.5))
)
contraction <- bind_rows(lapply(names(prior_draws), function(p) {
  post <- as.numeric(fit_kf$draws(p, format = "draws_matrix"))
  tibble(parameter = p,
         contraction  = round(1 - sd(post) / sd(prior_draws[[p]]), 3),
         prior_median = round(median(prior_draws[[p]]), 3),
         post_median  = round(median(post), 3))
}))
cat("\n--- Prior vs posterior contraction (0 = not learned) ---\n")
print(as.data.frame(contraction), row.names = FALSE)
cat(sprintf("\nRecords with n>1 (the only ones informing phi): %d / %d\n",
            sum(base$n_trials > 1), base$N))
