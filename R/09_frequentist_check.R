# ============================================================================
#  09_frequentist_check.R  —  Framework-independent robustness (no priors)
#
#  Shows the variance-partition conclusion (review >> NCP > trait) is NOT an
#  artefact of the Bayesian framework or the priors, by refitting a STANDARD
#  frequentist multilevel model and reading off its variance components.
#
#  STANDALONE: reads outputs/prepared_data.rds; needs no Bayesian fit and does
#  not re-run anything. Run on its own with
#     Rscript -e 'source("R/09_frequentist_check.R")'
#
#  Model: cumulative link mixed model (ordinal::clmm) on the COMPLETE records,
#  each primary-study direction expanded to one ordinal outcome (neg<null<pos),
#  with crossed random intercepts for review / trait / service / pair. This is
#  the exact frequentist analogue of the Bayesian cumulative-logit model.
#
#  HONEST LIMITS (state them): the frequentist check drops the 410 censored
#  records, ignores the Dirichlet overdispersion and the fixed Sigma overlap
#  covariance. So its "review" variance conflates genuine heterogeneity with
#  sampling (compare it to the Bayesian D_full review share, not D_struct). It
#  is a coarser, approximate parallel — the point is only that the variance
#  ORDER is preserved without any Bayesian machinery.
#
#  OUTPUT  outputs/frequentist_check.csv
# ============================================================================

source(here::here("R", "00_config.R"))

if (!requireNamespace("ordinal", quietly = TRUE))
  install.packages("ordinal", repos = "https://cloud.r-project.org")
library(ordinal)

prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
dm   <- prep$data_model

# ---- 1. Complete records, expanded to one row per primary-study direction --
comp <- dm[is.na(dm$cens_type) | trimws(dm$cens_type) == "", ]
long <- do.call(rbind, lapply(seq_len(nrow(comp)), function(i) {
  r <- comp[i, ]
  dirs <- c(rep("neg", r$Nneg), rep("null", r$Nnull), rep("pos", r$Npos))
  if (length(dirs) == 0) return(NULL)
  data.frame(resp = dirs, review = r$source_id, trait = r$trait_grouped,
             service = r$ncp, pair = r$ts, stringsAsFactors = FALSE)
}))
long$resp <- factor(long$resp, levels = c("neg", "null", "pos"), ordered = TRUE)
for (col in c("review", "trait", "service", "pair")) long[[col]] <- factor(long[[col]])
message("Frequentist check on ", nrow(long),
        " expanded ordinal outcomes (", nrow(comp), " complete records).")

# ---- 2. Fit clmm, with graceful fallback if the full crossed model struggles
fit_clmm <- function(form) tryCatch(
  clmm(form, data = long, Hess = FALSE),
  error = function(e) { message("  clmm failed: ", conditionMessage(e)); NULL })

message("Fitting clmm (review + trait + service + pair)...")
m <- fit_clmm(resp ~ 1 + (1 | review) + (1 | trait) + (1 | service) + (1 | pair))
if (is.null(m)) {
  message("Falling back to review + service + pair (dropping trait).")
  m <- fit_clmm(resp ~ 1 + (1 | review) + (1 | service) + (1 | pair))
}
stopifnot(!is.null(m))

# ---- 3. Variance components (each RE is scalar -> ST is a 1x1 SD) -----------
re_sd  <- sapply(m$ST, function(s) as.numeric(s))     # named by random term
re_var <- re_sd^2
V_link <- pi^2 / 3                                     # logistic link residual

comp_tab <- data.frame(
  component = names(re_var),
  variance  = as.numeric(re_var),
  stringsAsFactors = FALSE
)
comp_tab <- rbind(comp_tab, data.frame(component = "link_residual", variance = V_link))
comp_tab$share_of_total <- comp_tab$variance / sum(comp_tab$variance)

cat("\n=========== FREQUENTIST multilevel variance (clmm) ===========\n")
print(transform(comp_tab, variance = round(variance, 3),
                share_of_total = sprintf("%.1f%%", 100 * share_of_total)),
      row.names = FALSE)

# ---- 4. Side-by-side with the Bayesian decomposition (if available) ---------
vd_path <- file.path(OUT_DIR, "variance_decomposition.csv")
if (file.exists(vd_path)) {
  vd <- read.csv(vd_path)
  bayes_full <- vd[vd$basis == "D_full", c("component", "median")]
  cat("\n--- Bayesian shares of full latent variance (D_full), for comparison ---\n")
  print(transform(bayes_full, median = sprintf("%.1f%%", 100 * median)),
        row.names = FALSE)
  cat("\nCompare the frequentist 'review' share to the Bayesian review block\n",
      "(heterogeneity + sampling floor): both should place most variance at the\n",
      "review level. Matching ORDER = the result is not a Bayesian/prior artefact.\n",
      sep = "")
}

write.csv(comp_tab, file.path(OUT_DIR, "frequentist_check.csv"), row.names = FALSE)
message("\nSaved outputs/frequentist_check.csv")

# ---- Optional second opinion (metafor 3-level meta-analysis) ---------------
#  If you also want the field-standard meta-analytic tool, on a per-record
#  directional score d = (Npos - Nneg)/ntot with binomial-type variance:
#    library(metafor)
#    comp$d  <- (comp$Npos - comp$Nneg) / comp$ntot
#    comp$vd <- (1 - comp$d^2) / comp$ntot            # crude sampling variance
#    m2 <- rma.mv(d, vd, random = ~ 1 | review/pair, data = comp)
#    # sigma2[1] = review-level, sigma2[2] = pair-level variance components.
