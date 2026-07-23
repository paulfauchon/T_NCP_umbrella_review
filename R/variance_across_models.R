# ============================================================================
#  variance_across_models.R  —  Is review dominance model-dependent?
#
#  Answers the "two-model = cherry-picking" objection by computing the review
#  share of structural variance in a PARAMETERISATION-FREE way for every model:
#     review_share = mean_j(tau_j^2) / ( mean_j(tau_j^2) + Var_k(alpha_k) )
#  where Var_k(alpha_k) is the realised variance of the pair effects (whatever
#  their internal structure: additive means, none, or latent factors). If review
#  dominates in qp AND simple AND low-rank AND the frequentist clmm, the split is
#  robust, not a modelling choice.
#
#  STANDALONE: reads the saved fits + outputs/frequentist_check.csv.
#     Rscript -e 'source("R/variance_across_models.R")'
#  OUTPUT  outputs/variance_across_models.csv
# ============================================================================

source(here::here("R", "00_config.R"))

rev_share_bayes <- function(name) {
  fp <- file.path(FIT_DIR, paste0("fit_", name, ".rds"))
  if (!file.exists(fp)) return(NULL)
  fit <- readRDS(fp)
  J   <- length(grep("^tau_j\\[", fit$metadata()$model_params, value = TRUE))
  if (J == 0) J <- sum(grepl("^tau_j\\[", dimnames(fit$draws())$variable))
  a  <- fit$draws("alpha", format = "matrix")                 # draws x K_ts
  tj <- fit$draws("tau_j", format = "matrix")                 # draws x J
  V_pair <- apply(a, 1, var)                                  # per-draw pair variance
  V_rev  <- rowMeans(tj^2)                                    # per-draw review het.
  share  <- V_rev / (V_rev + V_pair)
  data.frame(model = name,
             review_share = median(share),
             q05 = quantile(share, .05, names = FALSE),
             q95 = quantile(share, .95, names = FALSE))
}

tab <- bind_rows(lapply(c("qp", "simple", "lowrank"), rev_share_bayes))

# Frequentist clmm (from 09): review variance / total structural variance
fq <- file.path(OUT_DIR, "frequentist_check.csv")
if (file.exists(fq)) {
  cl <- read.csv(fq)
  struct <- cl$variance[cl$component != "link_residual"]
  names(struct) <- cl$component[cl$component != "link_residual"]
  tab <- bind_rows(tab, data.frame(
    model = "clmm (frequentist)",
    review_share = struct["review"] / sum(struct),
    q05 = NA, q95 = NA))
}

tab$review_share <- round(tab$review_share, 3)
write.csv(tab, file.path(OUT_DIR, "variance_across_models.csv"), row.names = FALSE)

cat("\n=== Review share of structural variance, across model families ===\n")
print(transform(tab,
        review_share = sprintf("%.0f%%", 100 * review_share),
        q05 = ifelse(is.na(q05), "", sprintf("%.0f%%", 100 * q05)),
        q95 = ifelse(is.na(q95), "", sprintf("%.0f%%", 100 * q95))),
      row.names = FALSE)
cat("\nIf review dominates in ALL rows, the two-model design is 'most interpretable\n",
    "parameterisation per question', not cherry-picking. Report this table.\n", sep = "")
