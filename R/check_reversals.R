# ============================================================================
#  check_reversals.R  —  Are "reversed" relationships real, or an artefact of
#                        cross-NCP borrowing through the trait/service means?
#
#  For every pair, compares the posterior direction under:
#    * qp     (full model: trait + service means -> borrowing across NCPs)
#    * simple (pair-only : NO trait/service means -> NO borrowing)
#  A reversal/deflation that PERSISTS in `simple` is data-driven (real). One that
#  DISAPPEARS in `simple` was produced by the trait main effect dragging the pair
#  toward the trait's behaviour on other NCPs (borrowing artefact) — those must
#  not be presented as "the model overturns the literature".
#
#  STANDALONE: needs outputs/fits/fit_qp.rds and fit_simple.rds.
#     Rscript -e 'source("R/check_reversals.R")'
#  OUTPUT  outputs/reversal_check.csv
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
