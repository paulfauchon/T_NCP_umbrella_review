# ============================================================================
#  10_biome_moderator.R  —  Separate ecological CONTEXT from residual review
#                          heterogeneity (the artefact-vs-context question)
#
#  The primary model leaves the review effect u_j unexplained, so "review
#  variance" conflates ecological context (biome/climate) with synthesis
#  choices. Here we add review-level moderators (ecosystem, climate) to the
#  review effect and split its variance into:
#     * CONTEXT      = var over reviews of X_rev %*% beta_rev  (explained by biome)
#     * RESIDUAL     = mean_j tau_j^2                          (unexplained: synthesis
#                                                               + unmeasured context)
#  If context explains little, "review-level dominance is not just biome" becomes
#  a result. If it explains a lot, the honest reading shifts toward context.
#
#  STANDALONE: reads outputs/prepared_data.rds + data/meta-data2.csv; refits the
#  primary (qp) model once with covariates. Compare to the base fit_qp.rds.
#     Rscript -e 'source("R/10_biome_moderator.R")'
#
#  OUTPUT  outputs/biome_moderator.csv
# ============================================================================

source(here::here("R", "00_config.R"))

prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
sd   <- prep$stan_data
review_levels <- prep$review_levels
J <- sd$J

# ---- 1. Review-level moderators, in the canonical review order -------------
meta <- read.csv(file.path(DATA_DIR, "meta-data2.csv"), sep = ";",
                 fileEncoding = "latin1", stringsAsFactors = FALSE)
meta$source_id <- trimws(as.character(meta$source_id))
idx <- match(review_levels, meta$source_id)
stopifnot(!any(is.na(idx)))
eco <- factor(meta$ecosystem[idx])
cli <- factor(meta$climate[idx])
cat("Ecosystem per review:\n"); print(table(eco))
cat("Climate per review:\n");   print(table(cli))

# Design matrix (treatment contrasts; reference level folded into the thresholds).
# Toggle: ecosystem only, or ecosystem + climate.
X_rev <- model.matrix(~ eco)[, -1, drop = FALSE]         # ecosystem only
# X_rev <- model.matrix(~ eco + cli)[, -1, drop = FALSE]  # + climate
colnames(X_rev) <- make.names(colnames(X_rev))
cat("\nModerator design matrix: ", J, " reviews x ", ncol(X_rev),
    " covariates (", paste(colnames(X_rev), collapse = ", "), ")\n", sep = "")

# ---- 2. Refit the primary model WITH the moderators ------------------------
spec <- MODEL_SPECS %>% filter(name == PRIMARY_MODEL)
sd$trait_structure <- spec$trait_structure
sd$quality_mode    <- spec$quality_mode
sd$P     <- ncol(X_rev)
sd$X_rev <- X_rev

mod <- cmdstan_model(STAN_FILE)
fit <- mod$sample(data = sd, seed = SEED, chains = 4, parallel_chains = 4,
                  iter_warmup = 1500, iter_sampling = 1500,
                  adapt_delta = 0.97, refresh = 500)
fit$save_object(file.path(FIT_DIR, "fit_qp_biome.rds"))
cat("\nDivergences:", sum(fit$diagnostic_summary()$num_divergent),
    "| max Rhat:", round(max(fit$summary()$rhat, na.rm = TRUE), 3), "\n")

# ---- 3. Split review variance into CONTEXT vs RESIDUAL ---------------------
d <- fit$draws(c("beta_rev", "tau_j", "sigma_trait", "sigma_service",
                 "sigma_alpha", "kappa"), format = "draws_matrix")
beta <- as.matrix(d[, paste0("beta_rev[", seq_len(ncol(X_rev)), "]")])
tauj2 <- as.matrix(d[, paste0("tau_j[", seq_len(J), "]")])^2

ctx_var <- apply(beta, 1, function(b) {              # var over reviews of X%*%b
  v <- as.numeric(X_rev %*% b); mean((v - mean(v))^2)
})
V_context  <- ctx_var
V_resid    <- rowMeans(tauj2)
V_trait    <- as.numeric(d[, "sigma_trait"])^2
V_service  <- as.numeric(d[, "sigma_service"])^2
V_pair     <- as.numeric(d[, "sigma_alpha"])^2

q <- function(x) sprintf("%.1f%% [%.1f-%.1f]",
                         100 * median(x), 100 * quantile(x, .05),
                         100 * quantile(x, .95))
struct <- V_trait + V_service + V_pair + V_context + V_resid

cat("\n============ REVIEW VARIANCE: CONTEXT vs RESIDUAL ============\n")
cat("Share of the review-level variance explained by biome (context):\n")
cat("  context / (context + residual) =",
    q(V_context / (V_context + V_resid)), "\n")
cat("\nStructural decomposition WITH moderators (share of structural variance):\n")
cat("  ecological context (biome) :", q(V_context / struct), "\n")
cat("  residual review heterogen. :", q(V_resid   / struct), "\n")
cat("  trait                      :", q(V_trait   / struct), "\n")
cat("  service (NCP)              :", q(V_service / struct), "\n")
cat("  pair residual              :", q(V_pair    / struct), "\n")

# Review share of STRUCTURAL variance, without vs with biome (shares, comparable)
vd_path <- file.path(OUT_DIR, "variance_decomposition.csv")
if (file.exists(vd_path)) {
  vd <- read.csv(vd_path)
  base_rev <- vd$median[vd$component == "review_heterogeneity" & vd$basis == "D_struct"]
  cat("\nReview share of structural variance (comparable, share-based):\n")
  cat("  without biome (all review) :", sprintf("%.1f%%", 100 * base_rev), "\n")
  cat("  with biome                 : context",
      sprintf("%.1f%%", 100 * median(V_context / struct)),
      "+ residual", sprintf("%.1f%%", 100 * median(V_resid / struct)), "\n")
}

cat("\nEcosystem coefficients (logit scale; + = more positive trait-NCP effects):\n")
print(round(fit$summary("beta_rev", "mean", ~quantile(.x, c(.05, .95)))[, -1], 2))

out <- data.frame(
  component = c("context_biome", "residual_review", "trait", "service", "pair"),
  share_of_structural = c(median(V_context/struct), median(V_resid/struct),
                          median(V_trait/struct), median(V_service/struct),
                          median(V_pair/struct)))
write.csv(out, file.path(OUT_DIR, "biome_moderator.csv"), row.names = FALSE)
cat("\nSaved outputs/biome_moderator.csv\n")
cat("\nINTERPRETATION: if 'context (biome)' is small and 'residual review' stays\n",
    "large, the review-level dominance is NOT reducible to biome -> supports the\n",
    "synthesis-level reading. If context is large, report context-dependence.\n", sep = "")
