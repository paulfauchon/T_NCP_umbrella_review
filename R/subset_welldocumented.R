# ============================================================================
#  subset_welldocumented.R  —  Does review dominance survive on rich data?
#
#  Reviewer request: the 85% review share might be a sparsity artefact (41% of
#  pairs rest on 1 study). Here we DROP the sparse pairs (keep only pairs with
#  >= N_MIN primary studies), refit the primary (qp) model on that richer subset,
#  and check whether the review share still dominates. If it holds when the lower
#  levels are better replicated, the sparsity-artefact objection is closed
#  (complements the fake-data recovery in 06).
#
#  STANDALONE: reads prepared_data.rds; refits qp on the filtered data only.
#     Rscript -e 'source("R/subset_welldocumented.R")'
#  OUTPUT  outputs/subset_welldocumented.csv
# ============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "_helpers.R"))

N_MIN <- 3L                       # keep pairs with >= N_MIN primary studies (try 5 too)

prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
dm   <- prep$data_model

# ---- 1. Keep only well-documented pairs ------------------------------------
pair_n <- dm %>% group_by(ts) %>% summarise(n = sum(ntot, na.rm = TRUE), .groups = "drop")
keep_ts <- pair_n$ts[pair_n$n >= N_MIN]
d <- dm[dm$ts %in% keep_ts, ]
message("Pairs kept (>= ", N_MIN, " studies): ", length(keep_ts), " / ",
        nrow(pair_n), " | records: ", nrow(d), " / ", nrow(dm))

# ---- 2. Re-index everything on the subset (mirrors 01_prepare_data) ---------
review_levels <- intersect(prep$review_levels, unique(d$source_id))
ts_levels     <- sort(unique(d$ts))
trait_levels  <- sort(unique(d$trait_grouped))
serv_levels   <- sort(unique(d$ncp))
cat_levels    <- sort(unique(d$trait_category))
J <- length(review_levels)

d <- d %>% mutate(
  j_id = match(source_id, review_levels),
  ts_id = match(ts, ts_levels),
  .t = match(trait_grouped, trait_levels),
  .s = match(ncp, serv_levels),
  .c = match(trait_category, cat_levels))
ts_h <- d %>% distinct(ts_id, .t, .s, .c) %>% arrange(ts_id)

Sigma <- prep$Sigma[review_levels, review_levels]
qual  <- prep$quality_j[match(review_levels, prep$review_levels)]

comp <- d[is.na(d$cens_type) | trimws(d$cens_type) == "", ]
cens <- d[!(is.na(d$cens_type) | trimws(d$cens_type) == ""), ]

# censored configs
cfg <- vector("list", nrow(cens)); cs <- integer(nrow(cens)); cc <- integer(nrow(cens)); tot <- 0L
for (i in seq_len(nrow(cens))) {
  m <- enumerate_configs(cens$cens_type[i], cens$ntot[i], cens$pd_max[i],
                         cens$pd_min[i], cens$discrepancy_min[i])
  cfg[[i]] <- m; cs[i] <- tot + 1L; cc[i] <- nrow(m); tot <- tot + nrow(m)
}
all_cfg <- if (nrow(cens) > 0) do.call(rbind, cfg) else matrix(c(0L,1L,0L), 1)

sd <- list(
  N = nrow(comp), K_ts = length(ts_levels), K_trait = length(trait_levels),
  K_service = length(serv_levels), K_cat = length(cat_levels), J = J, P = 0L,
  trait_structure = 2L, quality_mode = 2L, estimate_kappa = 1L, kappa_fixed = 1,
  link = 0L, use_dm = 1L, ps_sigma = 1, ps_tau = 1, ps_gamma = 0.5, ps_kappa_ln = 0.5,
  trait_id = ts_h$.t, service_id = ts_h$.s, cat_id = ts_h$.c,
  counts = as.matrix(comp[, c("Nneg","Nnull","Npos")]), n_trials = comp$ntot,
  ts_id = comp$ts_id, j_id = comp$j_id, study_w = rep(1, nrow(comp)),
  N_cens = nrow(cens), ts_id_cens = cens$ts_id, j_id_cens = cens$j_id,
  study_w_cens = rep(1, nrow(cens)), total_configs = max(tot, 1L),
  all_configs = all_cfg, config_start = if (nrow(cens)>0) cs else integer(0),
  config_count = if (nrow(cens)>0) cc else integer(0),
  quality_j = as.numeric(qual), Sigma = Sigma, X_rev = matrix(0, J, 0))

# ---- 3. Refit and decompose ------------------------------------------------
mod <- cmdstan_model(STAN_FILE)
fit <- mod$sample(data = sd, seed = SEED, chains = 4, parallel_chains = 4,
                  iter_warmup = 1500, iter_sampling = 1500,
                  adapt_delta = 0.97, refresh = 500)
cat("\nDivergences:", sum(fit$diagnostic_summary()$num_divergent),
    "| max Rhat:", round(max(fit$summary()$rhat, na.rm = TRUE), 3), "\n")

vd <- variance_decomposition(fit, Sigma, use_cat = FALSE)
het <- vd %>% filter(component == "review_heterogeneity", basis == "D_struct")

cat(sprintf("\n=== Review share on well-documented subset (>= %d studies) ===\n", N_MIN))
cat(sprintf("  pairs %d, reviews %d, records %d\n",
            length(ts_levels), J, nrow(comp) + nrow(cens)))
cat(sprintf("  review share of structural variance = %.0f%% [%.0f-%.0f]\n",
            100*het$median, 100*het$q05, 100*het$q95))
cat("  (full data was ~85%. If it stays high here, review dominance is NOT a\n",
    "   sparsity artefact of single-study pairs.)\n", sep = "")
write.csv(vd, file.path(OUT_DIR, "subset_welldocumented.csv"), row.names = FALSE)

