# ============================================================================
#  check_reversals.R  —  Are "reversed" relationships real, or an artefact of
#                        cross-NCP borrowing through the trait/service means?
#
#  For every pair, compares the posterior direction under:
#   * qp        (trait + NCP means + quality precision -> borrowing across NCPs)
#   * simple (pair-only + quality precision -> NO cross-NCP borrowing)
#  A reversal/deflation that PERSISTS in `simple` is data-driven (real). One that
#  DISAPPEARS in `simple` was produced by the trait main effect dragging the pair
#  toward the trait's behaviour on other NCPs (borrowing artefact) — those must
#  not be presented as "the model overturns the literature".
#
#  STANDALONE: needs outputs/fits/fit_qp.rds and fit_simple.rds.
#     Rscript -e 'source("R/check_reversals.R")'
#  OUTPUT  outputs/reversal_check.csv
#
#  ADDED: test for positivity bias
# ============================================================================

source(here::here("R", "00_config.R"))
prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
dm   <- prep$data_model

# P(alpha > 0) per pair, in ts_id order (same mapping in both fits)
p_dir <- function(name) {
  a <- readRDS(file.path(FIT_DIR, paste0("fit_", name, ".rds")))$draws(
    "alpha", format = "matrix")
  apply(a, 2, function(x) mean(x > 0))
}
P_qp     <- p_dir("qp")
P_simple <- p_dir("simple")

pair_map <- dm %>% distinct(ts_id, trait_grouped, ncp) %>% arrange(ts_id)
votes <- dm %>% group_by(ts_id) %>%
  summarise(n_studies = sum(ntot, na.rm = TRUE),
            n_reviews = n_distinct(source_id),
            vpos = sum(Npos, na.rm = TRUE) / sum(ntot, na.rm = TRUE),
            vneg = sum(Nneg, na.rm = TRUE) / sum(ntot, na.rm = TRUE),
            .groups = "drop")

tab <- pair_map %>%
  mutate(P_qp = P_qp[ts_id], P_simple = P_simple[ts_id]) %>%
  left_join(votes, by = "ts_id") %>%
  mutate(
    apparent      = 0.5 + 0.5 * (vpos - vneg),          # vote-counting positivity
    shift_qp      = P_qp - apparent,                    # how far the model moved it
    borrowing_gap = P_qp - P_simple,                    # qp minus no-borrowing
    # qp changes the apparent direction (deflation past 0.5 or reversal)
    qp_flips_vote = (apparent > 0.6 & P_qp < 0.5) | (apparent < 0.4 & P_qp > 0.5),
    # is that change also present WITHOUT borrowing?
    verdict = case_when(
      !qp_flips_vote                       ~ "consistent with vote",
      (P_qp < 0.5) == (P_simple < 0.5)     ~ "REAL (holds without borrowing)",
      TRUE                                 ~ "BORROWING ARTEFACT (gone without means)"
    ))

# Focus on well-studied pairs where the model departs from the vote
flagged <- tab %>%
  filter(n_reviews >= 2, n_studies >= 8, qp_flips_vote) %>%
  arrange(verdict, desc(abs(borrowing_gap))) %>%
  mutate(across(c(apparent, P_qp, P_simple, borrowing_gap), ~round(., 2)))

cat("\n=== Pairs where qp departs from vote-counting (well-studied) ===\n")
cat("apparent = vote positivity | P_qp = full model | P_simple = no-borrowing\n\n")
print(as.data.frame(flagged %>%
        select(trait_grouped, ncp, n_studies, apparent, P_qp, P_simple,
               borrowing_gap, verdict)), row.names = FALSE)

write.csv(tab, file.path(OUT_DIR, "reversal_check.csv"), row.names = FALSE)
cat("\nSaved outputs/reversal_check.csv (all pairs).\n")
cat("\nUse in the main figure: keep only departures marked REAL; drop or caveat\n",
    "the BORROWING ARTEFACT ones (they reflect the trait mean, not the pair's data).\n",
    sep = "")

#% of same direction relationships between reviews that shared a pair

# Direction of each review's relationship within each pair
review_directions <- dm %>%
  group_by(ts_id, source_id) %>%
  summarise(
    npos = sum(Npos, na.rm = TRUE),
    nneg = sum(Nneg, na.rm = TRUE),
    direction = case_when(
      npos > nneg ~ "positive",
      nneg > npos ~ "negative",
      TRUE        ~ "tie"
    ),
    .groups = "drop"
  )

review_agreement <- review_directions %>%
  filter(direction != "tie") %>%
  group_by(ts_id) %>%
  summarise(
    n_directional_reviews = n(),
    n_review_comparisons = choose(n(), 2),
    n_same_direction = if (n() < 2) {
      0L
    } else {
      sum(combn(
        direction, 2,
        FUN = function(x) x[1] == x[2]
      ))
    },
    pct_same_direction = if (n_review_comparisons > 0) {
      100 * n_same_direction / n_review_comparisons
    } else {
      NA_real_
    },
    .groups = "drop"
  )

# Overall percentage across all comparable review pairs
overall_agreement <- review_agreement %>%
  summarise(
    total_review_pairs = sum(n_review_comparisons),
    total_same_direction = sum(n_same_direction),
    pct_same_direction = if (total_review_pairs > 0) {
      100 * total_same_direction / total_review_pairs
    } else {
      NA_real_
    }
  )

print(overall_agreement)


perfect_agreement <- review_agreement %>%
  filter(n_review_comparisons > 0) %>%
  summarise(
    n_shared_relationships = n(),
    n_perfect_agreement = sum(n_same_direction == n_review_comparisons),
    pct_perfect_agreement = 100 * n_perfect_agreement / n_shared_relationships
  )

print(perfect_agreement)

#test for positivity bias

names(dm)
dplyr::glimpse(dm)

dm %>%
  summarise(
    n_rows = n(),
    n_relationships = n_distinct(ts_id),
    n_reviews = n_distinct(source_id),
    total_ntot = sum(ntot, na.rm = TRUE),
    total_positive = sum(Npos, na.rm = TRUE),
    total_negative = sum(Nneg, na.rm = TRUE)
  )

dm %>%
  count(ts_id, source_id) %>%
  filter(n > 1)

library(dplyr)
library(ggplot2)
library(scales)

relationship_trend <- dm %>%
  group_by(ts_id, trait_grouped, ncp) %>%
  summarise(
    n_studies = sum(ntot, na.rm = TRUE),
    n_positive = sum(Npos, na.rm = TRUE),
    n_negative = sum(Nneg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_directional = n_positive + n_negative,
    prop_positive = if_else(
      n_directional > 0,
      n_positive / n_directional,
      NA_real_
    ),
    relationship = paste(trait_grouped, ncp, sep = " – ")
  ) %>%
  filter(
    n_studies > 0,
    n_directional > 0
  )

p <- ggplot(
  relationship_trend,
  aes(x = n_studies, y = prop_positive)
) +
  geom_point(alpha = 0.65, size = 2) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = label_percent()
  ) +
  scale_x_log10() +
  labs(
    x = "Number of studies or outcomes per trait–NCP relationship",
    y = "Proportion of positive directional outcomes"
  ) +
  theme_classic()

p

m_trend <- glm(
  cbind(n_positive, n_negative) ~ log(n_studies),
  family = binomial(),
  data = relationship_trend
)

summary(m_trend)
confint(m_trend)

