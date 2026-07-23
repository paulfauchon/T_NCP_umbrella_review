# ============================================================================
#  06_recovery_simulation.R  —  Does the model recover the variance components?
#
#  Refutes the objection "the small trait-NCP variance is a shrinkage artefact".
#  We simulate from the generative model on the REAL design (same pairs, reviews,
#  sample sizes, Sigma) and check the variance decomposition recovers the truth.
#
#  Corrections after statistical review:
#   * truth and estimator now use the SAME per-effect basis (mean_j tau_j^2);
#   * MANY replicate datasets per scenario -> report bias AND interval coverage
#     (two single draws proved nothing);
#   * a dedicated scenario exercises the CENSORED marginalisation path.
#
#  Heavy: N_REP * n_scenarios model fits. Lower N_REP for a quick check.
#  OUTPUT  outputs/recovery_simulation.csv
# ============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "_helpers.R"))

prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
sd0  <- prep$stan_data
mod  <- cmdstan_model(STAN_FILE)

N_REP <- 10L            # replicate datasets per scenario (raise to 20+ for final)

# ---- Generative helpers ----------------------------------------------------
rdirmult <- function(n_trials, prob, phi) {
  a <- phi * prob + 1e-6
  g <- rgamma(length(a), shape = a, rate = 1)
  as.integer(rmultinom(1, n_trials, g / sum(g)))
}

# True per-effect structural shares implied by a "truth" (common basis with
# the estimator: review heterogeneity = mean_j tau_j^2, no link term here since
# we compare STRUCTURAL shares).
true_structural_shares <- function(truth) {
  Vt <- truth$sigma_trait^2; Vs <- truth$sigma_service^2
  Vp <- truth$sigma_alpha^2
  Vh <- mean((truth$tau * exp(truth$gamma * (1 - sd0$quality_j)))^2)   # per-effect mean
  tot <- Vt + Vs + Vp + Vh
  c(trait = Vt/tot, service_NCP = Vs/tot,
    pair_residual = Vp/tot, review_heterogeneity = Vh/tot)
}

# Simulate one dataset. censor = TRUE re-censors the originally-censored records
# (using each record's own simulated counts) to exercise the marginalisation.
simulate_dataset <- function(truth, censor = FALSE) {
  J <- sd0$J
  tau_j <- truth$tau * exp(truth$gamma * (1 - sd0$quality_j))
  Cov   <- sd0$Sigma; diag(Cov) <- diag(Cov) + tau_j^2
  u     <- as.numeric(t(chol(Cov)) %*% rnorm(J))

  mu_trait   <- rnorm(sd0$K_trait,   0, truth$sigma_trait)
  mu_service <- rnorm(sd0$K_service, 0, truth$sigma_service)
  eps        <- rnorm(sd0$K_ts,      0, truth$sigma_alpha)
  alpha <- mu_trait[sd0$trait_id] + mu_service[sd0$service_id] + eps

  ts_all <- c(sd0$ts_id, sd0$ts_id_cens)
  j_all  <- c(sd0$j_id,  sd0$j_id_cens)
  n_cens <- vapply(seq_len(sd0$N_cens),
                   function(i) sum(sd0$all_configs[sd0$config_start[i], ]),
                   integer(1))
  n_all  <- c(sd0$n_trials, n_cens)

  counts <- t(vapply(seq_along(ts_all), function(i) {
    eta <- alpha[ts_all[i]] + u[j_all[i]]
    p   <- c(plogis(truth$c1 - eta),
             plogis(truth$c2 - eta) - plogis(truth$c1 - eta),
             1 - plogis(truth$c2 - eta))
    rdirmult(n_all[i], p, truth$phi)
  }, integer(3)))

  if (!censor) {
    return(list(counts = counts, ts_id = ts_all, j_id = j_all, n_trials = n_all,
                N_cens = 0L))
  }

  # Keep the first N complete; re-censor the last N_cens using their own counts
  N <- sd0$N; Nc <- sd0$N_cens
  comp_idx <- seq_len(N); cens_idx <- N + seq_len(Nc)
  cfg_list <- vector("list", Nc); cs <- integer(Nc); cc <- integer(Nc); tot <- 0L
  for (i in seq_len(Nc)) {
    ci <- counts[cens_idx[i], ]; n <- sum(ci)
    d  <- (ci[3] - ci[1]) / max(n, 1)
    # disc_min = 0: no discrepancy filter, so the true (incl. all-null) triple is
    # always inside its own enumerated set (avoids a silent likelihood mis-spec).
    cfg <- enumerate_configs("pd", n, pd_max = d, pd_min = d, disc_min = 0)
    cfg_list[[i]] <- cfg; cs[i] <- tot + 1L; cc[i] <- nrow(cfg); tot <- tot + nrow(cfg)
  }
  list(counts = counts[comp_idx, , drop = FALSE],
       ts_id = ts_all[comp_idx], j_id = j_all[comp_idx], n_trials = n_all[comp_idx],
       N_cens = Nc, ts_id_cens = ts_all[cens_idx], j_id_cens = j_all[cens_idx],
       all_configs = do.call(rbind, cfg_list), config_start = cs,
       config_count = cc, total_configs = tot)
}

fit_sim <- function(sim) {
  sd <- sd0
  sd$trait_structure <- 2L; sd$quality_mode <- 2L
  sd$N <- length(sim$ts_id); sd$counts <- sim$counts
  sd$n_trials <- sim$n_trials; sd$ts_id <- sim$ts_id; sd$j_id <- sim$j_id
  if (sim$N_cens == 0L) {
    sd$N_cens <- 0L; sd$total_configs <- 1L
    sd$all_configs <- matrix(c(0L, 1L, 0L), nrow = 1)
    sd$config_start <- integer(0); sd$config_count <- integer(0)
    sd$ts_id_cens <- integer(0);  sd$j_id_cens <- integer(0)
  } else {
    sd$N_cens <- sim$N_cens; sd$total_configs <- sim$total_configs
    sd$all_configs <- sim$all_configs; sd$config_start <- sim$config_start
    sd$config_count <- sim$config_count
    sd$ts_id_cens <- sim$ts_id_cens; sd$j_id_cens <- sim$j_id_cens
  }
  # study weights must be re-sized to the simulated N / N_cens (default 1)
  sd$study_w      <- rep(1, sd$N)
  sd$study_w_cens <- rep(1, sd$N_cens)
  mod$sample(data = sd, seed = SEED, chains = 4, parallel_chains = 4,
             iter_warmup = 800, iter_sampling = 800,
             adapt_delta = 0.97, refresh = 0)
}

# ---- Scenarios (variance components on the logit scale) --------------------
#  The KEY test (reviewer objection #1): under a trait/pair-dominated truth with
#  review ~ 0, does the model recover a SMALL review share, or does sparse lower
#  levels spuriously inflate review toward 80%? Scenario C is that stress test.
scenarios <- list(
  A_strong_trait  = list(sigma_trait = 0.9, sigma_service = 0.9, sigma_alpha = 0.3,
                         tau = 0.2, gamma = 0.1, c1 = -1.0, c2 = 0.5, phi = 10),
  B_strong_review = list(sigma_trait = 0.2, sigma_service = 0.2, sigma_alpha = 0.2,
                         tau = 1.2, gamma = 0.3, c1 = -1.0, c2 = 0.5, phi = 10),
  C_zero_review   = list(sigma_trait = 0.7, sigma_service = 0.7, sigma_alpha = 0.6,
                         tau = 0.05, gamma = 0.05, c1 = -1.0, c2 = 0.5, phi = 10)
)

# ---- Run replicated recovery (complete-data path) --------------------------
run_replicates <- function(nm, truth, censor = FALSE) {
  truth_shares <- true_structural_shares(truth)
  out <- vector("list", N_REP)
  for (r in seq_len(N_REP)) {
    set.seed(SEED + r)
    fit <- fit_sim(simulate_dataset(truth, censor = censor))
    vd  <- variance_decomposition(fit, prep$Sigma, use_cat = FALSE) %>%
      filter(basis == "D_struct",
             component %in% names(truth_shares))
    vd$truth   <- truth_shares[vd$component]
    vd$rep     <- r
    vd$covered <- vd$truth >= vd$q05 & vd$truth <= vd$q95
    out[[r]]   <- vd
  }
  bind_rows(out) %>% mutate(scenario = nm, censored = censor)
}

recovery <- bind_rows(
  run_replicates("A_strong_trait",  scenarios$A_strong_trait),
  run_replicates("B_strong_review", scenarios$B_strong_review),
  run_replicates("C_zero_review",   scenarios$C_zero_review),
  run_replicates("A_censored_path", scenarios$A_strong_trait, censor = TRUE)
)
write.csv(recovery, file.path(OUT_DIR, "recovery_simulation.csv"), row.names = FALSE)

# ---- Summarise: bias and 90% coverage per component ------------------------
summary_tab <- recovery %>%
  group_by(scenario, component) %>%
  summarise(
    true_share    = mean(truth),
    est_share     = mean(median),
    bias          = mean(median - truth),
    coverage_90   = mean(covered),
    .groups = "drop"
  ) %>%
  mutate(across(c(true_share, est_share, bias), ~round(., 3)))

cat("\n--- Recovery: bias & 90% coverage over", N_REP, "replicates ---\n")
print(as.data.frame(summary_tab), row.names = FALSE)

# HEADLINE: does the model inflate review when the truth is trait/pair-dominated?
cat("\n=== ANTI-SHRINKAGE TEST (reviewer objection #1) ===\n")
rev_row <- summary_tab %>% filter(component == "review_heterogeneity")
for (i in seq_len(nrow(rev_row))) {
  r <- rev_row[i, ]
  cat(sprintf("  %-16s TRUE review share = %.0f%%  ->  RECOVERED = %.0f%%%s\n",
              r$scenario, 100 * r$true_share, 100 * r$est_share,
              if (abs(r$bias) < 0.1) "   (no inflation)" else "   (CHECK: biased)"))
}
cat("If the model recovers a LOW review share when the truth is low (esp. scenario\n",
    "C_zero_review), the 85% headline is NOT a shrinkage artefact of sparse levels.\n", sep = "")
cat("\nGood recovery = bias near 0 and coverage near 0.90 in BOTH scenarios,\n",
    "including the censored-path scenario. This shows the decomposition is not\n",
    "biased toward 'review dominates' and recovers trait-NCP variance when present.\n",
    sep = "")
