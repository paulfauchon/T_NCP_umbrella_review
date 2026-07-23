# ============================================================================
#  ppc_by_pair.R  —  Ecology-oriented posterior predictive check (corrected)
#
#  Proper Bayesian PPC: for each summary statistic T, compare the OBSERVED value
#  to the posterior-predictive DISTRIBUTION of T across y_rep draws, and report a
#  Bayesian p-value = P(T_rep >= T_obs). Well-calibrated => p in ~(0.05, 0.95);
#  values near 0 or 1 flag a systematic misfit of that feature.
#
#  (A record-level 90% interval coverage is uninformative here because the median
#  record has n=1, so a discrete [0,1] interval trivially contains the observed
#  0/1 — we therefore use summary-statistic PPCs, and report coverage only for
#  the n>1 records where it is meaningful.)
#
#  STANDALONE (no MCMC): reads fit_qp.rds (y_rep) + prepared_data.rds. Seconds.
#     Rscript -e 'source("R/ppc_by_pair.R")'
#  OUTPUT  outputs/ppc_by_pair.csv
# ============================================================================

source(here::here("R", "00_config.R"))
prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
sd   <- prep$stan_data
fit  <- readRDS(file.path(FIT_DIR, "fit_qp.rds"))

N     <- sd$N
obs   <- sd$counts                                  # N x 3 (neg, null, pos)
n_tr  <- rowSums(obs)
yr    <- fit$draws("y_rep", format = "matrix")      # draws x (N*3)
D     <- nrow(yr)
cols  <- function(d) paste0("y_rep[", seq_len(N), ",", d, "]")
Yneg <- yr[, cols(1)]; Ynull <- yr[, cols(2)]; Ypos <- yr[, cols(3)]  # draws x N

# ---- Summary statistics: observed value + posterior-predictive draws --------
tot <- sum(n_tr)
stat <- function(obs_val, rep_vec, name) {
  p <- mean(rep_vec >= obs_val)
  data.frame(statistic = name, observed = round(obs_val, 3),
             rep_median = round(median(rep_vec), 3),
             rep_q05 = round(quantile(rep_vec, .05, names = FALSE), 3),
             rep_q95 = round(quantile(rep_vec, .95, names = FALSE), 3),
             bayes_p = round(p, 3))
}

# 1-3: overall category proportions
res <- bind_rows(
  stat(sum(obs[, 1]) / tot, rowSums(Yneg)  / tot, "prop_negative"),
  stat(sum(obs[, 2]) / tot, rowSums(Ynull) / tot, "prop_neutral"),
  stat(sum(obs[, 3]) / tot, rowSums(Ypos)  / tot, "prop_positive"),
  # 4: dispersion of the per-record positive fraction (does the model reproduce
  #    the bimodality/unanimity of single-study records?)
  {
    sd_obs <- sd(obs[, 3] / n_tr)
    sd_rep <- apply(Ypos, 1, function(y) sd(y / n_tr))
    stat(sd_obs, sd_rep, "sd_positive_fraction")
  },
  # 5: fraction of records that are unanimously positive
  {
    u_obs <- mean(obs[, 3] == n_tr)
    u_rep <- rowMeans(Ypos == matrix(n_tr, D, N, byrow = TRUE))
    stat(u_obs, u_rep, "frac_unanimous_positive")
  }
)

cat("\n=== POSTERIOR PREDICTIVE CHECK (Bayesian p-values) ===\n")
cat("Well-calibrated: bayes_p in ~(0.05, 0.95). Near 0/1 = that feature misfit.\n\n")
print(as.data.frame(res), row.names = FALSE)

# Coverage on the informative (n>1) records only
idx <- which(n_tr > 1)
if (length(idx) > 0) {
  cov <- sapply(1:3, function(d) {
    M <- yr[, cols(d)][, idx, drop = FALSE]
    lo <- apply(M, 2, quantile, .05); hi <- apply(M, 2, quantile, .95)
    mean(obs[idx, d] >= lo & obs[idx, d] <= hi)
  })
  cat(sprintf("\n90%% interval coverage on n>1 records (%d of %d): neg %.2f, null %.2f, pos %.2f\n",
              length(idx), N, cov[1], cov[2], cov[3]))
}

write.csv(res, file.path(OUT_DIR, "ppc_by_pair.csv"), row.names = FALSE)
cat("\nSaved outputs/ppc_by_pair.csv\n")
