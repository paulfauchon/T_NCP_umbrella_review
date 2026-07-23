# ============================================================================
#  07_results_pairs.R  —  Per-pair posterior support vs. vote counting
#
#  For the primary model, computes for every trait-NCP pair:
#    - posterior P(latent effect > 0)  and  P(< 0)
#    - the raw vote-counting proportions (positive / neutral / negative)
#  so the manuscript can contrast Bayesian support with naive aggregation.
#
#  OUTPUT  outputs/pair_results.csv
# ============================================================================

source(here::here("R", "00_config.R"))
prep       <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
data_model <- prep$data_model

fit   <- readRDS(file.path(FIT_DIR, paste0("fit_", PRIMARY_MODEL, ".rds")))
alpha <- fit$draws("alpha", format = "matrix")     # [draws x K_ts]

# ---- Posterior directional support -----------------------------------------
posterior_support <- tibble(
  ts_id         = seq_len(ncol(alpha)),
  mean_alpha    = apply(alpha, 2, mean),
  q05_alpha     = apply(alpha, 2, quantile, 0.05),
  q95_alpha     = apply(alpha, 2, quantile, 0.95),
  prob_positive = apply(alpha, 2, function(x) mean(x > 0)),
  prob_negative = apply(alpha, 2, function(x) mean(x < 0))
)

# ---- Raw vote counting (per pair, across all primary studies) --------------
vote_counts <- data_model %>%
  group_by(ts_id) %>%
  summarise(
    n_reviews = n_distinct(source_id),
    n_studies = sum(ntot, na.rm = TRUE),
    vote_pos  = sum(Npos,  na.rm = TRUE) / sum(ntot, na.rm = TRUE),
    vote_neg  = sum(Nneg,  na.rm = TRUE) / sum(ntot, na.rm = TRUE),
    vote_neu  = sum(Nnull, na.rm = TRUE) / sum(ntot, na.rm = TRUE),
    .groups   = "drop"
  )

pair_labels <- data_model %>%
  distinct(ts_id, ts, trait_category, trait_grouped, ncp_category, ncp) %>%
  arrange(ts_id)

pair_results <- pair_labels %>%
  left_join(posterior_support, by = "ts_id") %>%
  left_join(vote_counts,       by = "ts_id") %>%
  arrange(desc(prob_positive))

write.csv(pair_results, file.path(OUT_DIR, "pair_results.csv"),
          row.names = FALSE)

cat("\nMost strongly supported POSITIVE pairs:\n")
print(head(pair_results %>%
             select(ts, prob_positive, vote_pos, n_studies) %>%
             mutate(across(where(is.numeric), ~round(., 3))), 10),
      row.names = FALSE)

cat("\nMost strongly supported NEGATIVE pairs:\n")
print(head(pair_results %>% arrange(desc(prob_negative)) %>%
             select(ts, prob_negative, vote_neg, n_studies) %>%
             mutate(across(where(is.numeric), ~round(., 3))), 10),
      row.names = FALSE)

cat(sprintf("\nPairs with P>0.90 in either direction: %d / %d\n",
            sum(pair_results$prob_positive > 0.90 |
                  pair_results$prob_negative > 0.90),
            nrow(pair_results)))
