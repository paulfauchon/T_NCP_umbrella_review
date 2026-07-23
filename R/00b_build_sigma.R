# ============================================================================
#  00b_build_sigma.R  —  Construct the review dependence matrix Sigma from the
#                        primary-paper lists (fixes the "wrong-scale Sigma"
#                        limitation raised in review).
#
#  INPUT   data/2nd order SR Traits - Services.xlsx  (tabs E1_primary_papers,
#          E2, E3): per review (source_id) the DOIs of its primary studies.
#  OUTPUT  data/sigma.csv  (22 x 22, review order = `reviews` below)
#
#  ----------------------------------------------------------------------------
#  WHAT Sigma IS NOW (and why this is the correct fix)
#  ----------------------------------------------------------------------------
#  Sigma is a UNITLESS DEPENDENCE STRUCTURE, not a variance on an arbitrary
#  scale. With B the review x unique-study incidence matrix (B[j,d]=1 iff review
#  j includes primary study d) and n_j = number of unique studies in review j:
#
#     Sigma[j,j] = 1 / n_j                         (per-study sampling term)
#     Sigma[j,k] = shared(j,k) / (n_j * n_k)       (shared-study covariance)
#
#  i.e. Sigma = D^{-1} (B B^T) D^{-1}, D = diag(n_j). Because B B^T is a Gram
#  matrix, Sigma is positive semi-definite by construction.
#
#  The LOGIT-SCALE per-study sampling variance is carried by the model parameter
#  `kappa` (Var(u_j) = kappa/n_j + tau_j^2). This is the key correction: the raw
#  1/n_j is only the *structure*; the *scale* is estimated, and it is identified
#  because n_j spans ~12-213 (an 18x spread) so the heteroscedastic kappa/n_j
#  term is distinguishable from the flat tau. Hence, unlike the previous fixed
#  1/n_j-as-variance matrix, nothing on the wrong scale is baked in.
# ============================================================================

source(here::here("R", "00_config.R"))
suppressPackageStartupMessages(library(readxl))

XLSX <- file.path(DATA_DIR, "2nd order SR Traits - Services.xlsx")

# Canonical review order (must match study_scores.csv / data_model source_id)
reviews <- c("3","171","235","276","B","244","54","165","161","177","124",
             "332","314","A","125","151","55","382","50","238","232","35")

# Normalise a DOI from any messy cell
norm_doi <- function(x) {
  x <- tolower(trimws(as.character(x)))
  m <- regmatches(x, regexpr("10\\.\\d{4,9}/[^\\s\"\\)\\]]+", x, perl = TRUE))
  if (length(m) == 0) return(NA_character_)
  sub("\\.$", "", m)
}

# Pull (source_id, DOI) pairs from one sheet given candidate DOI columns
pull <- function(sheet, doi_cols) {
  df <- suppressMessages(read_excel(XLSX, sheet = sheet, col_types = "text"))
  if (!"source_id" %in% names(df)) return(NULL)
  present <- doi_cols[doi_cols %in% names(df)]
  out <- lapply(seq_len(nrow(df)), function(i) {
    doi <- NA_character_
    for (cc in present) { d <- norm_doi(df[[cc]][i]); if (!is.na(d)) { doi <- d; break } }
    if (is.na(doi)) NULL else data.frame(source_id = trimws(df$source_id[i]),
                                         doi = doi, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

pairs <- unique(rbind(
  pull("E1_primary_papers", c("DOI", "doi_article", "crossref_doi", "crossref_result_DOI")),
  pull("E2",                c("doi")),
  pull("E3",                c("DOI"))
))
pairs <- pairs[pairs$source_id %in% reviews & !is.na(pairs$doi), ]

# DOI set per review
doi_sets <- lapply(reviews, function(j) unique(pairs$doi[pairs$source_id == j]))
names(doi_sets) <- reviews
n_j <- vapply(doi_sets, length, integer(1))
if (any(n_j == 0)) stop("Reviews with no DOI: ", paste(reviews[n_j == 0], collapse = ", "))

# Shared-study counts C, then Sigma = D^{-1} C D^{-1}
J <- length(reviews)
C <- matrix(0, J, J, dimnames = list(reviews, reviews))
for (a in seq_len(J)) for (b in a:J) {
  s <- length(intersect(doi_sets[[a]], doi_sets[[b]]))
  C[a, b] <- s; C[b, a] <- s
}
Dinv  <- diag(1 / n_j)
Sigma <- Dinv %*% C %*% Dinv
dimnames(Sigma) <- list(reviews, reviews)
Sigma <- (Sigma + t(Sigma)) / 2

# Validate PSD
ev <- min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
if (ev <= 0) { Sigma <- Sigma + diag(abs(ev) + 1e-8, J)
  warning("Sigma nudged to PD (min eig=", signif(ev, 3), ").") }

write.csv(Sigma, file.path(DATA_DIR, "sigma.csv"), row.names = FALSE)

message("Built data/sigma.csv from ", nrow(pairs), " (review, DOI) pairs.")
message("n_j spread: ", min(n_j), "-", max(n_j),
        " (", round(max(n_j) / min(n_j), 1), "x -> kappa is identified).")
message("Overlap ratio sum(Sigma)/sum(diag): ",
        round(sum(Sigma) / sum(diag(Sigma)), 3),
        " | min eigenvalue: ", signif(min(eigen(Sigma, symmetric = TRUE,
                                                 only.values = TRUE)$values), 3))
