# ============================================================================
#  _helpers.R  —  Reusable functions (censoring enumeration, variance decomp)
# ============================================================================

# ----------------------------------------------------------------------------
#  enumerate_configs(): all (Nneg, Nnull, Npos) triples compatible with a
#  censored record. "pd" = predominant-direction only; "sig" = significant but
#  undirected (Nnull = 0). Returns a 3-column integer matrix (>=1 row) or errors
#  loudly if an unsatisfiable window yields zero configurations.
# ----------------------------------------------------------------------------
enumerate_configs <- function(cens_type, ntot, pd_max, pd_min, disc_min) {
  cens_type <- trimws(as.character(cens_type))
  out <- list()
  if (cens_type == "pd") {
    diff_min <- floor(pd_min * ntot)
    diff_max <- ceiling(pd_max * ntot)
    for (npos in 0:ntot) for (nneg in 0:ntot) {
      nnull <- ntot - npos - nneg
      if (nnull < 0) next
      d <- npos - nneg
      if (d < diff_min || d > diff_max) next
      if (disc_min == 1 && npos < disc_min && nneg < disc_min) next
      if (disc_min == 2 && (npos < disc_min || nneg < disc_min)) next
      out[[length(out) + 1]] <- c(nneg, nnull, npos)
    }
  } else if (cens_type == "sig") {
    for (npos in 0:ntot) out[[length(out) + 1]] <- c(ntot - npos, 0L, npos)
  }
  m <- do.call(rbind, out)
  if (is.null(m) || nrow(m) < 1)
    stop("enumerate_configs: no compatible configuration for cens_type='",
         cens_type, "', ntot=", ntot, " (check pd_min/pd_max/disc_min).")
  m
}

# ----------------------------------------------------------------------------
#  variance_decomposition(): latent-scale variance-components partition.
#
#  ** Corrected after statistical review (per-effect footing). **
#  Every component is a PER-EFFECT latent (logit-scale) variance, so trait,
#  service, pair and review sit on a common footing (design-independent). The
#  earlier n_j/N weighting made the sampling term collapse to J/N and mixed
#  per-observation with per-effect quantities — removed.
#
#   V_trait   = sigma_trait^2
#   V_service = sigma_service^2
#   V_cat     = sigma_cat^2               (only if use_cat)
#   V_pair    = sigma_alpha^2
#   V_rev_het = mean_j tau_j^2            genuine residual between-review heterogeneity
#   V_rev_smp = kappa * mean_j Sigma_jj   fixed sampling floor (kappa scales it; has
#                                          a posterior when kappa is estimated)
#   V_link    = pi^2/3                    logistic-link residual (latent-scale ICC)
#
#  Two denominators, both reported:
#   D_struct = V_cat+V_trait+V_service+V_pair+V_rev_het        (synthesis-driven)
#   D_full   = D_struct + V_rev_smp + V_link                   (full latent variance)
#
#  phi (DM overdispersion) is a COUNT-scale quantity and is NOT folded into the
#  latent denominator; report it separately. Sigma off-diagonal overlap is a
#  reported diagnostic, not a variance component (a single obs sees one u_j).
#
#  Returns a tidy tibble: component, basis (D_struct/D_full), median, q05, q95.
# ----------------------------------------------------------------------------
variance_decomposition <- function(fit, Sigma, use_cat = FALSE,
                                   use_trait = TRUE, use_service = TRUE) {
  vars_wanted <- c("sigma_trait", "sigma_service", "sigma_cat",
                   "sigma_alpha", "kappa")
  J <- nrow(Sigma)
  d <- fit$draws(
    variables = c(vars_wanted, paste0("tau_j[", seq_len(J), "]")),
    format = "draws_matrix")
  n_draw <- nrow(d)

  g   <- function(v) as.numeric(d[, v])
  tau_j2 <- as.matrix(d[, paste0("tau_j[", seq_len(J), "]")])^2

  # Gate dormant levels: a level not used by the model samples its prior, so its
  # sigma must NOT be counted as variance (this otherwise inflates, e.g., the
  # "trait" share of a service-only model with the trait prior).
  V_trait   <- if (use_trait)   g("sigma_trait")^2   else rep(0, n_draw)
  V_service <- if (use_service) g("sigma_service")^2 else rep(0, n_draw)
  V_cat     <- if (use_cat)     g("sigma_cat")^2     else rep(0, n_draw)
  V_pair    <- g("sigma_alpha")^2
  V_rev_het <- rowMeans(tau_j2)                      # per-effect (1/J) sum tau_j^2
  V_rev_smp <- g("kappa") * mean(diag(Sigma))        # fixed floor, scaled by kappa
  V_link    <- pi^2 / 3

  D_struct <- V_cat + V_trait + V_service + V_pair + V_rev_het
  D_full   <- D_struct + V_rev_smp + V_link

  q  <- function(x) c(median = median(x),
                      q05 = quantile(x, .05, names = FALSE),
                      q95 = quantile(x, .95, names = FALSE))
  add <- function(comp, basis, num, den)
    c(comp, basis, q(num / den))

  rows <- list(
    add("organ",                "D_struct", V_cat,     D_struct),
    add("trait",                "D_struct", V_trait,   D_struct),
    add("service_NCP",          "D_struct", V_service, D_struct),
    add("pair_residual",        "D_struct", V_pair,    D_struct),
    add("review_heterogeneity", "D_struct", V_rev_het, D_struct),

    add("organ",                "D_full",   V_cat,     D_full),
    add("trait",                "D_full",   V_trait,   D_full),
    add("service_NCP",          "D_full",   V_service, D_full),
    add("pair_residual",        "D_full",   V_pair,    D_full),
    add("review_heterogeneity", "D_full",   V_rev_het, D_full),
    add("review_sampling_floor","D_full",   V_rev_smp, D_full),
    add("link_residual",        "D_full",   V_link,    D_full)
  )

  tibble::tibble(
    component = vapply(rows, `[`, "", 1),
    basis     = vapply(rows, `[`, "", 2),
    # rows are character vectors (comp/basis are strings), so extract as "" then coerce
    median    = as.numeric(vapply(rows, `[`, "", 3)),
    q05       = as.numeric(vapply(rows, `[`, "", 4)),
    q95       = as.numeric(vapply(rows, `[`, "", 5))
  )
}

# ----------------------------------------------------------------------------
#  overlap_diagnostic(): summarise Sigma off-diagonal mass (reviewer request).
#  ratio > 1 signals shared-primary-study overlap not captured by the diagonal.
# ----------------------------------------------------------------------------
overlap_diagnostic <- function(Sigma) {
  list(
    overlap_ratio     = sum(Sigma) / sum(diag(Sigma)),
    mean_offdiag_corr = {
      D <- sqrt(diag(Sigma)); R <- Sigma / outer(D, D)
      mean(R[upper.tri(R)])
    }
  )
}
