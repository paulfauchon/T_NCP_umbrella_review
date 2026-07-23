# ============================================================================
#  01_prepare_data.R  —  Build the Stan data list from the raw extraction
#
#  INPUT   data/data_model.csv   (one row = one trait-NCP record within a review)
#          data/sigma.csv        (review x review sampling+overlap covariance)
#          data/study_scores.csv (review methodological-quality scores)
#  OUTPUT  outputs/stan_data.rds  (+ objects kept in memory for later scripts)
#
#  KEY FIX vs. the original scripts: reviews are aligned to ONE canonical order
#  used consistently for j_id, quality_j and Sigma. In the original pipeline
#  j_id came from unique(source_id) while quality/Sigma came from rownames(Sigma);
#  if those orders differ the review effects, quality and covariance are silently
#  mismatched. Here everything is keyed on `review_levels` explicitly.
# ============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "_helpers.R"))

# ---- 1. Load ---------------------------------------------------------------
data_model   <- read.csv(file.path(DATA_DIR, "data_model.csv"),
                         stringsAsFactors = FALSE)
Sigma_raw    <- as.matrix(read.csv(file.path(DATA_DIR, "sigma.csv"),
                                   check.names = FALSE))
study_scores <- read.csv(file.path(DATA_DIR, "study_scores.csv"),
                         stringsAsFactors = FALSE)

clean_id <- function(x) sub("^X", "", as.character(x))
colnames(Sigma_raw) <- clean_id(colnames(Sigma_raw))
rownames(Sigma_raw) <- colnames(Sigma_raw)

# ---- 2. Canonical review order (single source of truth) --------------------
review_levels <- colnames(Sigma_raw)                 # define the order once
stopifnot(all(unique(data_model$source_id) %in% review_levels))

Sigma <- Sigma_raw[review_levels, review_levels]     # explicit alignment
Sigma <- (Sigma + t(Sigma)) / 2                      # enforce exact symmetry
min_eig <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
if (min_eig <= 0) {                                  # CSV rounding can break PD
  Sigma <- Sigma + diag(abs(min_eig) + 1e-8, nrow(Sigma))
  warning("Sigma was not positive-definite (min eig=", signif(min_eig, 3),
          "); nudged to PD.")
}
J <- length(review_levels)

# NOTE (from statistical review): diag(Sigma) is treated as a FIXED sampling
# floor and MUST be supplied on the latent logit scale (delta-method from n_j
# and typical category probabilities), NOT as a raw 1/n_j count heuristic.
# The construction of sigma.csv should be shipped and documented; see README.

# ---- 3. Neutral column: make Nnull explicit --------------------------------
#  The model needs counts = (Nneg, Nnull, Npos). "Nno" (no direction) is the
#  neutral category. Alias once so Stan and figures share the same source.
if (!"Nnull" %in% names(data_model) && "Nno" %in% names(data_model))
  data_model$Nnull <- data_model$Nno
stopifnot(all(c("Nneg", "Nnull", "Npos", "ntot") %in% names(data_model)))

# ---- 4. Pair / trait / service / organ indices -----------------------------
ts_levels      <- sort(unique(data_model$ts))
trait_levels   <- sort(unique(data_model$trait_grouped))
service_levels <- sort(unique(data_model$ncp))
cat_levels     <- sort(unique(data_model$trait_category))

data_model <- data_model %>%
  mutate(
    j_id       = match(source_id,      review_levels),
    ts_id      = match(ts,             ts_levels),
    .trait_id  = match(trait_grouped,  trait_levels),
    .serv_id   = match(ncp,            service_levels),
    .cat_id    = match(trait_category, cat_levels)
  )

K_ts <- length(ts_levels); K_trait <- length(trait_levels)
K_service <- length(service_levels); K_cat <- length(cat_levels)

# One (trait, service, organ) triplet per pair, in pair-id order
ts_hierarchy <- data_model %>%
  distinct(ts_id, .trait_id, .serv_id, .cat_id) %>%
  arrange(ts_id)
stopifnot(nrow(ts_hierarchy) == K_ts)   # each pair maps to exactly one triplet

# ---- 5. Review quality: scale to (0,1] on the canonical order --------------
Q_raw <- study_scores$prop_score[match(review_levels, study_scores$source_id)]
if (any(is.na(Q_raw))) {
  warning("Reviews without a quality score set to the median.")
  Q_raw[is.na(Q_raw)] <- median(Q_raw, na.rm = TRUE)
}
quality_j <- (Q_raw - min(Q_raw)) / (max(Q_raw) - min(Q_raw))   # min-max -> [0,1]
quality_j <- pmax(quality_j, 0.05)                              # avoid exact 0

# ---- 6. Split complete vs. censored records --------------------------------
is_cens        <- !is.na(data_model$cens_type) & trimws(data_model$cens_type) != ""
data_complete  <- data_model[!is_cens, ]
data_cens      <- data_model[ is_cens, ]
message("Complete records: ", nrow(data_complete),
        " | Censored: ", nrow(data_cens))

# ---- 7. Enumerate compatible configurations for censored records -----------
cfg_list  <- vector("list", nrow(data_cens))
cfg_start <- integer(nrow(data_cens))
cfg_count <- integer(nrow(data_cens))
total     <- 0L
for (i in seq_len(nrow(data_cens))) {
  cfg <- enumerate_configs(
    cens_type = data_cens$cens_type[i],
    ntot      = data_cens$ntot[i],
    pd_max    = data_cens$pd_max[i],
    pd_min    = data_cens$pd_min[i],
    disc_min  = data_cens$discrepancy_min[i]
  )
  cfg_list[[i]] <- cfg
  cfg_start[i]  <- total + 1L
  cfg_count[i]  <- nrow(cfg)
  total         <- total + nrow(cfg)
}
all_configs <- do.call(rbind, cfg_list)
message("Total censored configurations enumerated: ", total)

# ---- 7b. Integrity checks --------------------------------------------------
stopifnot(all(rowSums(as.matrix(
  data_complete[, c("Nneg", "Nnull", "Npos")])) == data_complete$ntot))

# ---- 7c. Review-level covariates (optional; default: none) -----------------
#  To separate synthesis heterogeneity from ecological context (statistical
#  review point I1), plug an ecosystem/climate/year design matrix here, one row
#  per review in `review_levels` order. Default P = 0 (residual heterogeneity
#  only). Example to enable later:
#    X_rev <- model.matrix(~ ecosystem + climate, review_meta[review_levels, ])[, -1]
X_rev <- matrix(0, nrow = J, ncol = 0)

# ---- 8. Assemble the Stan data list ----------------------------------------
#  trait_structure / quality_mode are placeholders here; 02_fit_models.R
#  overwrites them per model specification.
stan_data <- list(
  N         = nrow(data_complete),
  K_ts      = K_ts, K_trait = K_trait, K_service = K_service, K_cat = K_cat,
  J         = J,
  P         = ncol(X_rev),
  trait_structure = 2L, quality_mode = 2L,          # overwritten per spec
  # Sigma is now a UNITLESS dependence STRUCTURE (built by 00b_build_sigma.R),
  # so kappa carries the per-study logit-scale sampling variance and is ESTIMATED.
  # It is identified because n_j spans ~12-213 (the heteroscedastic kappa/n_j term
  # is distinguishable from the flat tau). 08 keeps fixed-kappa robustness checks.
  estimate_kappa  = 1L,
  kappa_fixed     = 1.0,                            # used only in 08 sensitivity
  link            = 0L,                             # 0 logit (primary), 1 probit
  use_dm          = 1L,                             # 1 Dirichlet-multinomial (primary)
  # Prior scales as data (prior-sensitivity = data change)
  ps_sigma        = 1.0,
  ps_tau          = 1.0,
  ps_gamma        = 0.5,
  ps_kappa_ln     = 0.5,                            # per-study logit variance, weakly informative

  trait_id   = ts_hierarchy$.trait_id,
  service_id = ts_hierarchy$.serv_id,
  cat_id     = ts_hierarchy$.cat_id,

  counts   = as.matrix(data_complete[, c("Nneg", "Nnull", "Npos")]),
  n_trials = data_complete$ntot,
  ts_id    = data_complete$ts_id,
  j_id     = data_complete$j_id,
  # Fractional study weights (1/multiplicity) to blunt double-counting of shared
  # primary studies. Default 1 (no downweighting) until primary-study IDs are
  # available; see README limitation #2 and the study-effect TODO.
  study_w      = rep(1, nrow(data_complete)),

  N_cens        = nrow(data_cens),
  ts_id_cens    = data_cens$ts_id,
  j_id_cens     = data_cens$j_id,
  study_w_cens  = rep(1, nrow(data_cens)),
  total_configs = total,
  all_configs   = all_configs,
  config_start  = cfg_start,
  config_count  = cfg_count,

  quality_j = as.numeric(quality_j),
  Sigma     = Sigma,
  X_rev     = X_rev
)

saveRDS(list(stan_data = stan_data, data_model = data_model,
             review_levels = review_levels, ts_levels = ts_levels,
             quality_j = quality_j, Sigma = Sigma),
        file.path(OUT_DIR, "prepared_data.rds"))

message("Saved outputs/prepared_data.rds  (N=", stan_data$N,
        ", K_ts=", K_ts, ", J=", J, ")")
