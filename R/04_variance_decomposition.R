# ============================================================================
#  04_variance_decomposition.R  —  THE central result (corrected after review)
#
#  Latent-scale variance-components partition on a common PER-EFFECT footing.
#  Separates genuine residual between-review heterogeneity (tau_j^2) from the
#  fixed sampling floor (kappa * diag Sigma) and the logistic link residual
#  (pi^2/3). Reports shares over the full latent variance (D_full, headline)
#  and over the synthesis-driven structural variance (D_struct).
#
#  Honesty notes surfaced by the review and echoed to the console:
#   - the sampling floor is a FIXED input scaled by kappa; report kappa's
#     posterior (08 runs the fixed-kappa sensitivity);
#   - phi (count-scale overdispersion) is reported separately, not folded in;
#   - Sigma off-diagonal overlap is a diagnostic, not a variance component.
#
#  OUTPUT  outputs/variance_decomposition.csv
# ============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "_helpers.R"))

prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
fit  <- readRDS(file.path(FIT_DIR, paste0("fit_", PRIMARY_MODEL, ".rds")))
spec <- MODEL_SPECS %>% filter(name == PRIMARY_MODEL)
use_cat <- spec$trait_structure >= 3

vd <- variance_decomposition(fit, prep$Sigma, use_cat = use_cat)
write.csv(vd, file.path(OUT_DIR, "variance_decomposition.csv"), row.names = FALSE)

pct <- function(x) sprintf("%.1f%%", 100 * x)
get <- function(comp, basis) vd %>% filter(component == comp, basis == !!basis)

het_full   <- get("review_heterogeneity", "D_full")
het_struct <- get("review_heterogeneity", "D_struct")
smp_full   <- get("review_sampling_floor", "D_full")
lnk_full   <- get("link_residual", "D_full")

# Supporting posterior summaries
kappa <- fit$summary("kappa")
phi   <- fit$summary("phi")
gamma <- fit$summary("gamma")
ov    <- overlap_diagnostic(prep$Sigma)

cat("\n================ VARIANCE DECOMPOSITION (corrected) ================\n")
cat("\nShares of FULL latent variance (D_full = structural + sampling + link):\n")
print(vd %>% filter(basis == "D_full") %>%
        mutate(across(median:q95, pct)), row.names = FALSE)

cat("\nShares of SYNTHESIS-DRIVEN structural variance (D_struct):\n")
print(vd %>% filter(basis == "D_struct") %>%
        mutate(across(median:q95, pct)), row.names = FALSE)

cat("\n--- Headline numbers ---\n")
cat(sprintf("Genuine between-review heterogeneity: %s of D_full (90%% CI %s-%s)\n",
            pct(het_full$median), pct(het_full$q05), pct(het_full$q95)))
cat(sprintf("                                  ... %s of D_struct (90%% CI %s-%s)\n",
            pct(het_struct$median), pct(het_struct$q05), pct(het_struct$q95)))
cat(sprintf("Fixed sampling floor: %s of D_full | link residual: %s of D_full\n",
            pct(smp_full$median), pct(lnk_full$median)))

cat("\n--- Identification / honesty checks (report these) ---\n")
cat(sprintf("kappa (Sigma scale): median %.2f (90%% CI %.2f-%.2f) [prior median 1]\n",
            kappa$median, kappa$q5, kappa$q95))
cat(sprintf("gamma (quality->precision): median %.2f (90%% CI %.2f-%.2f) [prior ~half-N(0,0.5)]\n",
            gamma$median, gamma$q5, gamma$q95))
cat(sprintf("phi (count overdispersion, reported separately): median %.1f\n", phi$median))
cat(sprintf("Sigma overlap ratio (off-diagonal mass): %.3f | mean off-diag corr: %.3f\n",
            ov$overlap_ratio, ov$mean_offdiag_corr))
cat("===================================================================\n")

