# ============================================================================
#  05_grouping_sensitivity.R  —  Is the trait grouping driving the result?
#
#  Tests whether partial pooling by trait_grouped ("borrowing of strength")
#  is biased by a too-coarse grouping. If traits pooled together react
#  systematically differently, the trait-NCP variance is artificially deflated
#  and the review share artificially inflated.
#
#  Two complementary checks:
#   (1) RESOLUTION sensitivity: re-run the variance decomposition at three
#       trait resolutions (organ / trait_grouped / trait_grouped+organ nested)
#       and verify the conclusion is stable. Uses fits already produced:
#         - "ncp"      : service only (no trait pooling)      -> coarsest
#         - "qp"       : trait + service                       -> primary
#         - "qp_organ" : organ + trait + service (nested)      -> finest
#   (2) EXCHANGEABILITY diagnostic: within each trait_grouped, is the spread of
#       pair-level effects (alpha) explained by the underlying raw trait
#       (trait_uni)? A high within-group R^2 by trait_uni flags a bad grouping.
#
#  OUTPUT  outputs/grouping_sensitivity.csv
# ============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "_helpers.R"))

prep       <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
data_model <- prep$data_model

# ---- (1) Variance decomposition across trait resolutions -------------------
resolutions <- tribble(
  ~fit_name,   ~resolution,               ~use_cat, ~use_trait,
  "ncp",       "service only (no trait)",  FALSE,    FALSE,
  "qp",        "trait + service",          FALSE,    TRUE,
  "qp_organ",  "organ + trait + service",  TRUE,     TRUE
)

decomp_by_res <- bind_rows(lapply(seq_len(nrow(resolutions)), function(i) {
  r   <- resolutions[i, ]
  fp  <- file.path(FIT_DIR, paste0("fit_", r$fit_name, ".rds"))
  if (!file.exists(fp)) return(NULL)
  variance_decomposition(readRDS(fp), prep$Sigma, use_cat = r$use_cat,
                         use_trait = r$use_trait) %>%
    filter(basis == "D_struct") %>%
    mutate(resolution = r$resolution)
}))
write.csv(decomp_by_res, file.path(OUT_DIR, "grouping_sensitivity.csv"),
          row.names = FALSE)

cat("\n--- Review-heterogeneity share across trait resolutions ---\n")
print(decomp_by_res %>%
        filter(component == "review_heterogeneity") %>%
        mutate(across(median:q95, ~sprintf("%.1f%%", 100 * .))) %>%
        select(resolution, median, q05, q95), row.names = FALSE)
cat("If this share is stable across resolutions, the trait grouping is NOT\n",
    "driving the review dominance. If it drops sharply at finer resolution,\n",
    "the coarse grouping was masking real trait-NCP signal.\n", sep = "")

# ---- (2) Within-group exchangeability diagnostic ---------------------------
#  Needs the primary fit's pair-level posterior means + the raw trait (trait_uni).
fit_primary <- readRDS(file.path(FIT_DIR, paste0("fit_", PRIMARY_MODEL, ".rds")))
alpha_mean  <- fit_primary$summary("alpha", "mean")$mean   # length K_ts, pair order

pair_map <- data_model %>%
  distinct(ts_id, trait_grouped, trait_uni) %>%
  arrange(ts_id) %>%
  mutate(alpha = alpha_mean[ts_id])

# For groups containing >1 raw trait and >1 pair: variance of alpha explained
# by trait_uni (one-way ANOVA R^2). High R^2 = raw trait matters within group.
group_diag <- pair_map %>%
  group_by(trait_grouped) %>%
  filter(n() >= 3, n_distinct(trait_uni) >= 2) %>%
  summarise(
    n_pairs   = n(),
    n_rawtrait = n_distinct(trait_uni),
    r2_by_rawtrait = summary(lm(alpha ~ trait_uni))$r.squared,
    .groups = "drop"
  ) %>%
  arrange(desc(r2_by_rawtrait))

cat("\n--- Groups where the RAW trait explains within-group spread (top 10) ---\n")
print(head(group_diag, 10), row.names = FALSE)
cat("\nGroups with high r2_by_rawtrait are candidates for a too-coarse grouping\n",
    "and should be reviewed manually before finalising.\n", sep = "")

write.csv(group_diag, file.path(OUT_DIR, "grouping_exchangeability.csv"),
          row.names = FALSE)
