# ============================================================================
#  anatomy_reviews.R  —  What IS the review-level effect made of?
#
#  Answers the reviewers' demand: "if biome explains only ~10%, what are the
#  other ~75%?" Dissects the 22 review deviations u_j (from the primary qp fit)
#  against synthesis-level metadata: synthesis method (vote-counting vs none),
#  #databases searched, explicit trait/NCP definitions, scope restriction,
#  ecological framework, PRISMA screening, year, #primary studies, quality.
#
#  EXPLORATORY by necessity (n = 22 reviews, many candidate moderators, low
#  power, multiple comparisons) — report effect sizes, not just p-values.
#
#  STANDALONE (no MCMC): reads outputs/prepared_data.rds, fit_qp.rds,
#  data/meta-data2.csv. Runs in seconds; safe to run while other fits are going.
#     Rscript -e 'source("R/anatomy_reviews.R")'
#  OUTPUT  outputs/review_anatomy.csv
# ============================================================================

source(here::here("R", "00_config.R"))
prep <- readRDS(file.path(OUT_DIR, "prepared_data.rds"))
rl   <- prep$review_levels

# ---- 1. Review deviations u_j (posterior mean), in review order ------------
u <- readRDS(file.path(FIT_DIR, "fit_qp.rds"))$summary("u", "mean")$mean
stopifnot(length(u) == length(rl))

# ---- 2. Synthesis-level metadata, recoded ----------------------------------
meta <- read.csv(file.path(DATA_DIR, "meta-data2.csv"), sep = ";",
                 fileEncoding = "latin1", stringsAsFactors = FALSE)
meta$source_id <- trimws(as.character(meta$source_id))
m <- meta[match(rl, meta$source_id), ]
g <- function(x) trimws(as.character(x))
n_db <- g(m$Number_of_Databases_Queried)

studies <- prep$data_model %>% group_by(source_id) %>%
  summarise(n_studies = sum(ntot, na.rm = TRUE),
            n_records = n(), .groups = "drop")

df <- data.frame(
  u                = u,
  year             = suppressWarnings(as.numeric(g(m$year))),
  n_databases      = ifelse(n_db == "1", 1, ifelse(n_db == "02-mars", 2.5,
                       ifelse(n_db == "5+", 5, NA))),          # Excel mangled "2-3"
  method           = ifelse(grepl("vote", g(m$Statistical_Method_to_Assess_the_Effect),
                                  ignore.case = TRUE), "vote-counting", "none"),
  def_traits       = g(m$Specific_Traits_Definitions),
  def_services     = g(m$Specific_Services_Definitions),
  context          = g(m$Research_resticted_to_a_certain_context),
  eco_standpoint   = g(m$Ecological_standpoint),
  other_sources    = g(m$Other_sources_used),
  screening        = ifelse(g(m$Screening_Process_Detailled) == "PRISMA", "PRISMA", "other"),
  scope_ncp        = ifelse(g(m$number_of_ncp) == "multiple", "multiple", "single"),
  ecosystem        = g(m$ecosystem),
  climate          = g(m$climate),
  stringsAsFactors = FALSE
)
df$source_id <- rl
df <- left_join(df, studies, by = "source_id")
df$quality <- prep$quality_j

# ---- 3. Univariate screening of each moderator vs u_j -----------------------
assoc <- function(mod, u) {
  if (is.numeric(mod)) {
    ok <- !is.na(mod)
    ct <- suppressWarnings(cor.test(mod[ok], u[ok], method = "spearman"))
    data.frame(type = "numeric", effect = unname(ct$estimate), p = ct$p.value,
               detail = sprintf("Spearman rho = %.2f", ct$estimate))
  } else {
    f <- factor(mod)
    if (nlevels(f) < 2) return(data.frame(type = "const", effect = NA, p = NA, detail = "constant"))
    kt <- suppressWarnings(kruskal.test(u ~ f))
    gm <- tapply(u, f, mean)
    data.frame(type = "categorical", effect = diff(range(gm)), p = kt$p.value,
               detail = paste(sprintf("%s=%.2f", names(gm), gm), collapse = " | "))
  }
}

mods <- c("method", "def_services", "def_traits", "context", "eco_standpoint",
          "other_sources", "screening", "scope_ncp", "ecosystem", "climate",
          "year", "n_databases", "n_studies", "quality")
res <- bind_rows(lapply(mods, function(nm) cbind(moderator = nm, assoc(df[[nm]], df$u))))
res <- res %>% arrange(p) %>% mutate(effect = round(effect, 3), p = round(p, 3))

cat("\n=== ANATOMY OF THE REVIEW EFFECT (u_j) — exploratory, n = 22 ===\n")
cat("effect = Spearman rho (numeric) or spread of group means (categorical).\n")
cat("u_j > 0 = review reports MORE positive trait-NCP effects.\n\n")
print(as.data.frame(res), row.names = FALSE)

# ---- 4. How much do the strongest moderators JOINTLY explain? ---------------
top <- head(res$moderator[!is.na(res$p)], 3)
fml <- as.formula(paste("u ~", paste(top, collapse = " + ")))
r2  <- summary(lm(fml, data = df))$r.squared
cat(sprintf("\nJoint R^2 of u_j ~ %s : %.2f\n", paste(top, collapse = " + "), r2))
cat("(exploratory; n=22 so this over-fits — read it as an upper bound.)\n")

write.csv(res, file.path(OUT_DIR, "review_anatomy.csv"), row.names = FALSE)
cat("\nSaved outputs/review_anatomy.csv\n")
cat("\nEither reading is publishable:\n",
    " - a synthesis characteristic stands out -> you NAME what fabricates generality;\n",
    " - nothing stands out -> even obvious methodological axes don't explain the\n",
    "   between-review heterogeneity, sharpening 'the evidence is not synthesizable'.\n",
    sep = "")

