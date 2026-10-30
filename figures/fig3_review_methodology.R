# ============================================================================
#  fig3_review_methodology.R  —  Section-3 figure (R / ggplot2) + GLMM stats
#     "How a review is done predicts what it concludes"
#
#  Each dot is one review, placed by its positivity (share of directional
#  evidence reported as positive), split by three synthesis-rigour attributes.
#  Group means + the shift between them, annotated with the p-value from a
#  study-level binomial GLMM (random review & trait-NCP-pair effects), which
#  properly weights by evidence and controls for which relationships each review
#  covers. Canonical twin of figures/fig3_review_methodology.py (the .py reads
#  the p-values this script writes to outputs/fig3_glmm.csv).
#
#  INPUT   data/data_model.csv, data/meta-data2.csv
#  OUTPUT  outputs/fig3_glmm.csv         (attribute, OR, CI, p)
#          outputs/fig3_review_methodology.(png|svg)
#  RUN     Rscript figures/fig3_review_methodology.R
#  DEPS    dplyr, ggplot2, lme4. SVG needs svglite (optional).
#
#  ADDED:  Colinearity test / phi correlations bewteen 3 quality criteria
#          Script for alternative figure: Fig. S4 (number of weaknesses)
#         (see below)
# ============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(lme4) })

root <- tryCatch(here::here(), error = function(e) getwd())
DM   <- file.path(root, "data", "data_model.csv")
META <- file.path(root, "data", "meta-data2.csv")
for (p in c(DM, META)) if (!file.exists(p)) stop("Missing ", p)

dm <- read.csv(DM, stringsAsFactors = FALSE)
for (k in c("Npos", "Nneg")) dm[[k]] <- suppressWarnings(as.numeric(dm[[k]]))
dm[is.na(dm)] <- 0
# robust latin1 read (native reader truncates on invalid bytes): convert whole
# file to UTF-8 in memory, then parse with the proper CSV reader.
.raw <- readChar(META, file.info(META)$size, useBytes = TRUE)
.tmp <- tempfile(); writeChar(iconv(.raw, "latin1", "UTF-8", sub = ""), .tmp, eos = NULL)
meta <- read.csv(.tmp, sep = ";", stringsAsFactors = FALSE)
meta$source_id <- trimws(as.character(meta$source_id))
meta <- meta[!duplicated(meta$source_id), ]   # one row per review (attributes are review-level)
m   <- meta[match(dm$source_id, meta$source_id), ]
low <- function(x) trimws(tolower(as.character(x)))
dm$trait_def <- low(m$Specific_Traits_Definitions)
dm$method    <- ifelse(grepl("vote", low(m$Statistical_Method_to_Assess_the_Effect)), "yes_quant", "narrative")
dm$eco       <- low(m$Ecological_standpoint)
dmd <- dm[(dm$Npos + dm$Nneg) > 0, ]
#From  here J=21 because review #382 provide only significant observations but without direction given, these are then marginalized

# ---- (1) study-level binomial GLMM per attribute ---------------------------
glmm_one <- function(var, rig, loose) {
  d <- dmd[dmd[[var]] %in% c(rig, loose), ]
  d$g <- factor(ifelse(d[[var]] == loose, "loose", "rigorous"), levels = c("rigorous", "loose"))
  f <- glmer(cbind(Npos, Nneg) ~ g + (1 | source_id) + (1 | ts),
             data = d, family = binomial, control = glmerControl(optimizer = "bobyqa"))
  co <- summary(f)$coefficients["gloose", ]
  ci <- tryCatch(suppressMessages(confint(f, parm = "gloose", method = "Wald")),
                 error = function(e) c(NA, NA))
  data.frame(OR = exp(fixef(f)["gloose"]), OR_lo = exp(ci[1]), OR_hi = exp(ci[2]),
             p = co["Pr(>|z|)"], n_reviews = length(unique(d$source_id)),
             n_outcomes = sum(d$Npos + d$Nneg))
}
glmm <- bind_rows(
  cbind(attribute = "method",    glmm_one("method",    "yes_quant", "narrative")),
  cbind(attribute = "trait_def", glmm_one("trait_def", "yes",       "no")),
  cbind(attribute = "eco",       glmm_one("eco",       "yes",       "no")))
write.csv(glmm, file.path(root, "outputs", "fig3_glmm.csv"), row.names = FALSE)
plab <- function(a) { p <- glmm$p[glmm$attribute == a]
  if (is.na(p)) "" else if (p < 0.10) sprintf("p = %.2f", p) else "n.s." }

# ---- (2) per-review positivity ---------------------------------------------
rev <- dm %>% group_by(source_id) %>%
  summarise(p = sum(Npos), n = sum(Nneg),
            trait_def = first(trait_def), method = first(method), eco = first(eco),
            .groups = "drop") %>%
  filter(p + n > 0) %>% mutate(pos = p / (p + n))

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; INK <- "#0b0b0b"; SUB <- "#52514e"
spec <- tibble::tribble(
  ~block,                 ~rowlab,        ~y, ~sel,     ~col,   ~pl,
  "a Synthesis method",     "Quantitative",  6, "m_q",   BLUE,   "",
  "a Synthesis method",     "Narrative only",5, "m_n",   ORANGE, plab("method"),
  "b Trait definitions",    "Explicit",      4, "t_y",   BLUE,   "",
  "b Trait definitions",    "None",          3, "t_n",   ORANGE, plab("trait_def"),
  "c Ecological framework", "Explicit",      2, "e_y",   BLUE,   "",
  "c Ecological framework", "None",          1, "e_n",   ORANGE, plab("eco"))
SF <- list(m_q=function(d) d$method %in% "yes_quant", m_n=function(d) d$method %in% "narrative",
            t_y=function(d) d$trait_def %in% "yes",    t_n=function(d) d$trait_def %in% "no",
            e_y=function(d) d$eco %in% "yes",          e_n=function(d) d$eco %in% "no")

grp_vals <- lapply(spec$sel, function(k) rev$pos[SF[[k]](rev)])

# dots with deterministic beeswarm (stack within x-bins)
dots <- bind_rows(lapply(seq_len(nrow(spec)), function(i) {
  v <- grp_vals[[i]]
  if (!length(v)) return(NULL)
  data.frame(pos = v, col = spec$col[i], y0 = spec$y[i]) %>% arrange(pos) %>%
    mutate(bin = round(pos / 0.02)) %>% group_by(bin) %>%
    mutate(k = row_number()-1, y = y0 + (k %% 5 - 2) * 0.11) %>% ungroup()
}))
means <- spec %>% mutate(mp = sapply(grp_vals, function(v) if (length(v)) mean(v) else NA_real_),
                         nn = sapply(grp_vals, length))
# shift arrows: one per block (rigorous -> loose)
blk <- means %>% group_by(block) %>%
  summarise(x_rig = mp[col == BLUE], x_loose = mp[col == ORANGE],
            y_rig = y[col == BLUE], y_loose = y[col == ORANGE],
            pl = pl[col == ORANGE], .groups = "drop") %>%
  mutate(ymid = (y_rig + y_loose)/2, dpp = round((x_loose - x_rig)*100))

p <- ggplot() +
  geom_vline(xintercept = seq(.5,1,.1), color = "#f1f0eb", linewidth = .5) +
  geom_rect(data = blk, aes(xmin = x_rig, xmax = x_loose, ymin = y_loose, ymax = y_rig),
            fill = ORANGE, alpha = .06) +
  geom_point(data = dots, aes(pos, y, color = col), size = 2.3, alpha = .5, stroke = .3) +
  scale_color_identity() +
  geom_segment(data = blk, aes(x = x_rig, xend = x_loose, y = ymid, yend = ymid),
               arrow = arrow(length = unit(6,"pt"), type="closed"), linewidth = .7, color = INK) +
  geom_text(data = blk, aes(x = (x_rig+x_loose)/2, y = ymid + .30,
            label = paste0("+", dpp, " pts")), fontface = 2, size = 3.3, color = INK) +
  geom_text(data = blk, aes(x = (x_rig+x_loose)/2, y = ymid - .30, label = pl),
            fontface = 2, size = 3.0, color = "#0b0b0b") +
  # mean diamonds + labels
  geom_point(data = means, aes(mp, y, fill = col), shape = 23, size = 3.4,
             color = "white", stroke = .8) + scale_fill_identity() +
  geom_text(data = means, aes(mp, y + ifelse(col==BLUE, .34, -.40),
            label = paste0(round(mp*100), "%"), color = col), fontface = 2, size = 3.1) +
  # row labels
  geom_text(data = spec, aes(x = .485, y = y, label = rowlab), hjust = 1, size = 3.4, color = INK) +
  geom_text(data = means, aes(x = .485, y = y - .26, label = paste0("n=", nn)),
            hjust = 1, size = 2.5, color = SUB) +
  geom_text(data = distinct(spec, block, .keep_all = TRUE),
            aes(x = .27, y = y - .5, label = block), hjust = 0, fontface = 2, size = 4, color = INK) +
  #annotate("text", x = .58, y = 6.9, label = "← more rigorous, less positive", color = BLUE, fontface = 2, size = 3.2) +
  #annotate("text", x = .93, y = 6.9, label = "less rigorous, more positive →", color = ORANGE, fontface = 2, size = 3.2) +
  scale_x_continuous(limits = c(.24, 1.02), breaks = seq(.5,1,.1),
                     labels = paste0(seq(50,100,10), "%"), position = "top") +
  scale_y_continuous(limits = c(.3, 6.3)) +
  # labs(title = "How a review is done predicts what it concludes",
  #      subtitle = paste0("Each dot is one systematic review, placed by the share of its directional ",
  #        "evidence reported as positive.\nReviews using less rigorous methods report systematically ",
  #        "more positive trait-NCP relationships."),
  #      x = NULL, y = NULL,
  #      caption = paste0("Group means shown as diamonds; p from a study-level binomial GLMM ",
  #        "(2,203 outcomes; random review & trait-NCP-pair effects).\nBetween-review differences track ",
  #        "how evidence is synthesised, not the ecosystem studied (biome ~10%; recovery + subset checks in Supplementary).\n(n = ", nrow(rev), ", one review did not provide any directional evidence).")) +
  theme_void(base_size = 11) +
  theme(axis.text.x.top = element_text(size = 9, color = SUB),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 10, color = SUB, margin = margin(b = 10)),
        plot.caption = element_text(hjust = 0, size = 8.5, color = SUB),
        plot.margin = margin(16, 24, 14, 24),
        plot.background = element_rect(fill = "#ffffff", color = NA))

out <- file.path(root, "outputs", "fig3_review_methodology")
ggsave(paste0(out, ".png"), p, width = 11, height = 8, dpi = 200, bg = "#fcfcfb")
ok <- requireNamespace("svglite", quietly = TRUE)
if (ok) ggsave(paste0(out, ".svg"), p, width = 11, height = 8, bg = "#fcfcfb")
message("wrote ", out, ".png", if (ok) " (+ .svg)" else "", " ; outputs/fig3_glmm.csv")

# ---- Export and print figure values ----------------------------------------
figure_values <- means %>%
  select(block, rowlab, col, mp, nn) %>%
  mutate(
    positivity_pct = round(mp * 100, 1),
    n_reviews = nn
  ) %>%
  select(block, rowlab, positivity_pct, n_reviews) %>%
  left_join(
    blk %>%
      select(block, x_rig, x_loose, dpp, pl),
    by = "block"
  )

print(figure_values)

# One row per methodological attribute: raw group means and difference
differences <- blk %>%
  transmute(
    attribute = block,
    rigorous_pct = round(x_rig * 100, 1),
    less_rigorous_pct = round(x_loose * 100, 1),
    difference_pp = round((x_loose - x_rig) * 100, 1),
    p_label = pl
  )

print(differences)

write.csv(
  differences,
  file.path(root, "outputs", "fig3_group_differences.csv"),
  row.names = FALSE
)

# ============================================================================
# Colinearity
# ============================================================================

# ---- Collinearity among the three methodological weaknesses ---------------
review_flags <- rev %>%
  distinct(source_id, trait_def, method, eco) %>%
  mutate(
    weak_trait = case_when(
      trait_def == "no"  ~ 1L,
      trait_def == "yes" ~ 0L,
      TRUE ~ NA_integer_
    ),
    narrative_weakness = case_when(
      method == "narrative" ~ 1L,
      method == "yes_quant" ~ 0L,
      TRUE ~ NA_integer_
    ),
    weak_eco = case_when(
      eco == "no"  ~ 1L,
      eco == "yes" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

weak_vars <- c("weak_trait", "narrative_weakness", "weak_eco")
flags_complete <- review_flags %>%
  filter(if_all(all_of(weak_vars), ~ !is.na(.x)))

cat("Reviews with complete methodology data:",
    nrow(flags_complete), "\n")

# Pearson correlations between binary variables = phi correlations
cor_mat <- cor(flags_complete[, weak_vars])
print(round(cor_mat, 2))

write.csv(
  cor_mat,
  file.path(root, "outputs", "fig3_weakness_correlations.csv")
)

# Cross-tabulation: number of reviews with each combination of weaknesses
combos <- flags_complete %>%
  count(
    weak_trait, narrative_weakness, weak_eco,
    name = "n_reviews"
  )

all_combos <- expand.grid(
  weak_trait = 0:1,
  narrative_weakness = 0:1,
  weak_eco = 0:1
)

combos <- all_combos %>%
  left_join(
    combos,
    by = c("weak_trait", "narrative_weakness", "weak_eco")
  ) %>%
  mutate(
    n_reviews = ifelse(is.na(n_reviews), 0L, n_reviews),
    combination = paste0(
      "Trait undefined=", weak_trait,
      "; Narrative=", narrative_weakness,
      "; No ecological framework=", weak_eco
    )
  )

print(combos)

write.csv(
  combos,
  file.path(root, "outputs", "fig3_weakness_combinations.csv"),
  row.names = FALSE
)

# ---- Variance reduction associated with the rigour gradient ---------------
# Aggregate directional counts by review x trait-NCP pair.

model_var <- dm %>%
  group_by(source_id, ts) %>%
  summarise(
    Npos = sum(Npos, na.rm = TRUE),
    Nneg = sum(Nneg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Npos + Nneg > 0) %>%
  left_join(
    review_flags %>%
      select(source_id, all_of(weak_vars)),
    by = "source_id"
  ) %>%
  filter(if_all(all_of(weak_vars), ~ !is.na(.x))) %>%
  mutate(
    weaknesses = weak_trait + narrative_weakness + weak_eco,
    source_id = factor(source_id),
    ts = factor(ts)
  )

# Null model: no methodological predictor
fit_null <- glmer(
  cbind(Npos, Nneg) ~ 1 +
    (1 | source_id) + (1 | ts),
  data = model_var,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

# Rigour-gradient model: categorical count of weaknesses (0 to 3)
fit_rigour <- glmer(
  cbind(Npos, Nneg) ~ factor(weaknesses) +
    (1 | source_id) + (1 | ts),
  data = model_var,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

var_null <- as.data.frame(VarCorr(fit_null)) %>%
  filter(grp == "source_id") %>%
  pull(vcov)

var_rigour <- as.data.frame(VarCorr(fit_rigour)) %>%
  filter(grp == "source_id") %>%
  pull(vcov)

variance_reduction_pct <- 100 * (var_null - var_rigour) / var_null

cat("Review-level variance, null model:", var_null, "\n")
cat("Review-level variance, rigour model:", var_rigour, "\n")
cat("Reduction in review-level variance:",
    round(variance_reduction_pct, 1), "%\n")

# ============================================================================
# fig3_review_methodology.R
# Fig. 3: methodological weaknesses vs positive outcomes
# Supplementary Fig. S: the three weaknesses considered separately -original
#
# Inputs:
#   data/data_model.csv
#   data/meta-data2.csv
#
# Outputs:
#   outputs/fig3_methodology.svg
#   outputs/fig3_methodology.png
#   outputs/figS_review_methodology.svg -original figure, vertical
#   outputs/fig3_alt_glmm.csv
#
# GLMM: binomial, random intercepts for review and trait-NCP pair.
# ============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

# ---- Packages ---------------------------------------------------------------
required <- c("lme4", "ggplot2", "dplyr", "tidyr")
missing <- required[!vapply(required, requireNamespace,
                            logical(1), quietly = TRUE)]

if (length(missing)) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing, collapse = ", "),
    "\nRun: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "),
    "))"
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(emmeans)

# ---- Colours ---------------------------------------------------------------
col_blue   <- "#377EB8"
col_orange <- "#E68632"
col_grey   <- "#D9DDE2"

weakness_cols <- c(
  "0" = col_blue,
  "1" = "#FAD9B8",
  "2" = "#E9A15B",
  "3" = col_orange
)

# ---- Paths ------------------------------------------------------------------
wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
root <- if (basename(wd) == "figures") dirname(wd) else wd

DM   <- file.path(root, "data", "data_model.csv")
META <- file.path(root, "data", "meta-data2.csv")
OUT  <- file.path(root, "outputs")

if (!file.exists(DM)) stop("Missing input: ", DM)
if (!file.exists(META)) stop("Missing input: ", META)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----------------------------------------------------------------
clean_text <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[is.na(x)] <- ""
  x
}

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# Map metadata values to explicit yes/no categories.
yes_no <- function(x) {
  z <- clean_text(x)
  out <- rep(NA_character_, length(z))
  out[z %in% c("yes", "y", "oui", "1", "true")] <- "yes"
  out[z %in% c("no", "n", "non", "0", "false")] <- "no"
  out
}

# ---- Read data --------------------------------------------------------------
dm <- read.csv(
  DM, check.names = FALSE, na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

meta <- read.csv(
  META, sep = ";", fileEncoding = "latin1",
  check.names = FALSE, na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

required_dm <- c("source_id", "Npos", "Nneg", "ts_id")
required_meta <- c(
  "source_id",
  "Specific_Traits_Definitions",
  "Statistical_Method_to_Assess_the_Effect",
  "Ecological_standpoint"
)

miss_dm <- setdiff(required_dm, names(dm))
miss_meta <- setdiff(required_meta, names(meta))

if (length(miss_dm)) {
  stop("Missing columns in data_model.csv: ",
       paste(miss_dm, collapse = ", "))
}
if (length(miss_meta)) {
  stop("Missing columns in meta-data2.csv: ",
       paste(miss_meta, collapse = ", "))
}

dm <- dm %>%
  mutate(
    source_id = trimws(as.character(source_id)),
    ts_id = as.character(ts_id),
    Npos = to_num(Npos),
    Nneg = to_num(Nneg)
  ) %>%
  filter(
    !is.na(source_id), source_id != "",
    !is.na(ts_id), ts_id != ""
  )

dm$Npos[is.na(dm$Npos)] <- 0
dm$Nneg[is.na(dm$Nneg)] <- 0

if (any(dm$Npos < 0 | dm$Nneg < 0)) {
  stop("Npos/Nneg contain negative counts.")
}

meta$source_id <- trimws(as.character(meta$source_id))
meta <- meta[!is.na(meta$source_id) & meta$source_id != "", ]

# Match the original analysis: first metadata row per review.
meta <- meta[!duplicated(meta$source_id), ]

# ---- Define the three weaknesses --------------------------------------------
meta2 <- meta %>%
  transmute(
    source_id,
    trait_defined = yes_no(Specific_Traits_Definitions),
    method_text = clean_text(Statistical_Method_to_Assess_the_Effect),
    eco_defined = yes_no(Ecological_standpoint)
  ) %>%
  mutate(
    # A method explicitly mentioning "vote" is quantitative vote counting.
    narrative = case_when(
      grepl("vote", method_text, fixed = TRUE) ~ 0L,
      method_text != "" ~ 1L,
      TRUE ~ NA_integer_
    ),
    weak_trait = case_when(
      trait_defined == "no" ~ 1L,
      trait_defined == "yes" ~ 0L,
      TRUE ~ NA_integer_
    ),
    weak_eco = case_when(
      eco_defined == "no" ~ 1L,
      eco_defined == "yes" ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(source_id, weak_trait, narrative, weak_eco)

# ---- Review-level dots ------------------------------------------------------
# Each dot is one review.
# Positivity = sum(Npos) / sum(Npos + Nneg).
# Dot size = number of trait-NCP outcomes with directional evidence.

review_data <- dm %>%
  group_by(source_id) %>%
  summarise(
    Npos = sum(Npos),
    Nneg = sum(Nneg),
    n_outcomes = n_distinct(ts_id[(Npos + Nneg) > 0]),
    .groups = "drop"
  ) %>%
  mutate(
    total = Npos + Nneg,
    positivity = ifelse(total > 0, Npos / total, NA_real_)
  ) %>%
  filter(total > 0, !is.na(positivity)) %>%
  left_join(meta2, by = "source_id") %>%
  mutate(
    weaknesses = weak_trait + narrative + weak_eco
  )

# Complete cases for the three weakness attributes.
analysis_reviews <- review_data %>%
  filter(
    !is.na(weak_trait),
    !is.na(narrative),
    !is.na(weak_eco),
    !is.na(weaknesses)
  ) %>%
  mutate(weaknesses = as.integer(weaknesses))

if (nrow(analysis_reviews) < 5) {
  stop("Too few reviews with complete metadata to fit the model.")
}

if (n_distinct(analysis_reviews$source_id) < 2 ||
    n_distinct(dm$ts_id) < 2) {
  stop("Insufficient levels for the review/pair random effects.")
}

# ---- Build study/pair-level binomial observations ---------------------------
# Aggregate repeated rows within each review x trait-NCP pair.
model_data <- dm %>%
  group_by(source_id, ts_id) %>%
  summarise(
    Npos = sum(Npos),
    Nneg = sum(Nneg),
    .groups = "drop"
  ) %>%
  filter(Npos + Nneg > 0) %>%
  left_join(meta2, by = "source_id") %>%
  mutate(
    weaknesses = weak_trait + narrative + weak_eco
  ) %>%
  filter(
    !is.na(weak_trait),
    !is.na(narrative),
    !is.na(weak_eco),
    !is.na(weaknesses)
  )

model_data$source_id <- factor(model_data$source_id)
model_data$ts_id <- factor(model_data$ts_id)

# Treat the number of weaknesses as a categorical predictor so each of
# levels 0, 1, 2 and 3 gets its own adjusted estimate.
model_data$weakness_factor <- factor(
  model_data$weaknesses, levels = 0:3
)

observed_levels <- sort(unique(model_data$weaknesses))
missing_levels <- setdiff(0:3, observed_levels)

if (length(missing_levels)) {
  stop(
    "Cannot estimate all four levels (0-3 weaknesses). Missing level(s): ",
    paste(missing_levels, collapse = ", "),
    ". Check the metadata categories and review counts."
  )
}

# Ensure the model has enough reviews and pairs.
if (n_distinct(model_data$source_id) < 2 ||
    n_distinct(model_data$ts_id) < 2) {
  stop("Model data do not contain enough review/pair levels.")
}

fit <- lme4::glmer(
  cbind(Npos, Nneg) ~ weakness_factor +
    (1 | source_id) + (1 | ts_id),
  data = model_data,
  family = binomial(link = "logit"),
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)
  )
)

# ---- Fixed-effect estimates and 95% CIs on the probability scale ------------
# Predictions set both random intercepts to zero. These are adjusted
# conditional estimates, not raw group means.
beta <- lme4::fixef(fit)
V <- as.matrix(vcov(fit))

terms <- paste0("weakness_factor", 1:3)
if (!all(terms %in% names(beta))) {
  stop("Unexpected GLMM coefficient names: ",
       paste(names(beta), collapse = ", "))
}

eta <- c(beta["(Intercept)"],
         beta["(Intercept)"] + beta[terms])

se_eta <- numeric(4)
se_eta[1] <- sqrt(V["(Intercept)", "(Intercept)"])

for (i in 1:3) {
  term <- terms[i]
  se_eta[i + 1] <- sqrt(
    V["(Intercept)", "(Intercept)"] +
      V[term, term] +
      2 * V["(Intercept)", term]
  )
}

estimates <- data.frame(
  weaknesses = 0:3,
  estimate = plogis(eta),
  lower = plogis(eta - 1.96 * se_eta),
  upper = plogis(eta + 1.96 * se_eta)
)

counts <- analysis_reviews %>%
  count(weaknesses, name = "n_reviews")

estimates <- estimates %>%
  left_join(counts, by = "weaknesses") %>%
  mutate(
    n_reviews = ifelse(is.na(n_reviews), 0L, n_reviews),
    n_outcomes = vapply(
      weaknesses,
      function(k) sum(
        model_data$Npos[model_data$weaknesses == k] +
          model_data$Nneg[model_data$weaknesses == k] > 0
      ),
      numeric(1)
    )
  )

write.csv(
  estimates,
  file.path(OUT, "fig3_alt_glmm.csv"),
  row.names = FALSE
)

# ---- Pairwise comparisons and compact-letter display -----------------------

emm <- emmeans(fit, ~ weakness_factor)

# Pairwise comparisons with Tukey adjustment
pw <- pairs(emm, adjust = "tukey")
print(pw)

# Compact letter display: levels sharing a letter are not
# significantly different at the chosen alpha level.
cld_results <- multcomp::cld(
  emm,
  adjust = "tukey",
  Letters = letters,
  alpha = 0.05
)

letters_df <- as.data.frame(cld_results) %>%
  transmute(
    weaknesses = as.integer(as.character(weakness_factor)),
    sig_letter = gsub(" ", "", .group)
  )

print(letters_df)


# ---- Main Figure 3 ----------------------------------------------------------
# One dot per review, with area scaled by number of directional outcomes.
# Slight horizontal jitter is only for display.
set.seed(2026)

# Position letters just above the upper confidence limit

letters_df <- letters_df %>%
  mutate(label_y = 1.04)

plot_reviews <- analysis_reviews %>%
  mutate(
    x_jitter = weaknesses +
      runif(n(), min = -0.09, max = 0.09)
  )

x_counts <- estimates %>%
  mutate(
    label = paste0("n = ", n_reviews)
  )

p3 <- ggplot() +
  geom_point(
    data = plot_reviews,
    aes(
      x = x_jitter,
      y = positivity,
      size = n_outcomes,
      fill = factor(weaknesses)
    ),
    shape = 21,
    colour = "grey25",
    alpha = 0.85,
    stroke = 0.4
  ) +
  scale_fill_manual(
    values = weakness_cols,
    name = "Methodological\nweaknesses",
    breaks = as.character(0:3)
  ) +
  geom_errorbar(
    data = estimates,
    aes(x = weaknesses, ymin = lower, ymax = upper),
    width = 0.10, linewidth = 0.8, colour = "black"
  ) +
  geom_point(
    data = estimates,
    aes(x = weaknesses, y = estimate),
    shape = 18, size = 3.4, colour = "black"
  ) +
  geom_text(
    data = letters_df,
    aes(
      x = weaknesses,
      y = label_y,
      label = sig_letter
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5,
    colour = "black"
  ) +
  geom_text(
    data = x_counts,
    aes(x = weaknesses, y = 0.455, label = label),
    size = 3.6, colour = "black"
  ) +
  scale_x_continuous(
    breaks = 0:3,
    limits = c(-0.35, 3.35),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = seq(0.5, 1, by = 0.1),
    labels = function(x) paste0(round(100 * x), "%"),
    limits = c(0.44, 1.08),
    expand = expansion(mult = 0)
  ) +
  scale_size_area(
    max_size = 8,
    name = "Directional outcomes"
  ) +
  labs(
    #title = "Less rigorous reviews report more positive relationships",
    x = "Number of methodological weaknesses",
    y = "Positive outcomes"
    #,
    #caption = paste0("Black diamonds and bars: binomial GLMM estimates and 95% CIs. ", nrow(model_data)," review–trait–NCP observations; random intercepts for review ","and trait–NCP pair. Grey dots: individual reviews.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    plot.caption = element_text(size = 9, hjust = 0),
    plot.margin = margin(12, 18, 12, 12)
  )

ggsave(
  file.path(OUT, "fig3_methodology.svg"),
  p3, width = 9.5, height = 6.3, device = "svg"
)
ggsave(
  file.path(OUT, "fig3_methodology.png"),
  p3, width = 9.5, height = 6.3, dpi = 320
)

# ---- Supplementary Figure S: each weakness separately ----------------------
# Group labels and n are shown explicitly. Points are reviews; black diamonds
# show the observed group means. These panels are descriptive, not GLMM CIs.

supp_long <- analysis_reviews %>%
  select(
    source_id, positivity,
    weak_trait, narrative, weak_eco
  ) %>%
  pivot_longer(
    cols = c(weak_trait, narrative, weak_eco),
    names_to = "attribute",
    values_to = "weak"
  ) %>%
  mutate(
    attribute = recode(
      attribute,
      weak_trait = "Traits not defined",
      narrative = "Narrative synthesis",
      weak_eco = "No ecological framework"
    ),
    group = factor(
      ifelse(weak == 1, "Weakness present", "Weakness absent"),
      levels = c("Weakness absent", "Weakness present")
    )
  )

supp_counts <- supp_long %>%
  count(attribute, group, name = "n") %>%
  mutate(label = paste0("n = ", n))

supp_means <- supp_long %>%
  group_by(attribute, group) %>%
  summarise(
    mean_pos = mean(positivity, na.rm = TRUE),
    .groups = "drop"
  )

pS <- ggplot(
  supp_long,
  aes(x = group, y = positivity)
) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 2026),
    alpha = 0.65, size = 2.2, colour = "grey45"
  ) +
  geom_point(
    data = supp_means,
    aes(x = group, y = mean_pos),
    inherit.aes = FALSE,
    shape = 18, size = 3.5, colour = "black"
  ) +
  geom_text(
    data = supp_counts,
    aes(x = group, y = 0.455, label = label),
    inherit.aes = FALSE,
    size = 3.3
  ) +
  facet_wrap(~attribute, nrow = 1) +
  scale_y_continuous(
    breaks = seq(0.5, 1, by = 0.1),
    labels = function(x) paste0(round(100 * x), "%"),
    limits = c(0.44, 1.02),
    expand = expansion(mult = 0)
  ) +
  labs(
    title = "Review positivity by individual methodological weakness",
    x = NULL,
    y = "Positive outcomes",
    caption = "Grey dots: individual reviews. Black diamonds: unadjusted group means."
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 15, hjust = 1),
    plot.caption = element_text(size = 9, hjust = 0)
  )

ggsave(
  file.path(OUT, "figS_review_methodology.svg"),
  pS, width = 12, height = 5.8, device = "svg"
)

# ---- Report -----------------------------------------------------------------
cat("\n=== FIGURE 3: BINOMIAL GLMM ===\n")
cat("Reviews included:", nrow(analysis_reviews), "\n")
cat("Review-pair observations in GLMM:", nrow(model_data), "\n\n")
print(
  estimates %>%
    mutate(
      across(c(estimate, lower, upper), ~round(100 * ., 1))
    ) %>%
    select(
      weaknesses, n_reviews, n_outcomes,
      estimate, lower, upper
    )
)

if (all(c(0, 3) %in% estimates$weaknesses)) {
  e0 <- estimates$estimate[estimates$weaknesses == 0]
  e3 <- estimates$estimate[estimates$weaknesses == 3]
  
  cat(
    "\nAdjusted probability, 0 weaknesses:",
    round(100 * e0, 1), "%\n"
  )
  cat(
    "Adjusted probability, 3 weaknesses:",
    round(100 * e3, 1), "%\n"
  )
  cat(
    "Adjusted difference:",
    round(100 * (e3 - e0), 1), "percentage points\n"
  )
}

cat("\nSaved files in:", OUT, "\n")
cat(" - fig3_methodology.svg / .png\n")
cat(" - figS_review_methodology.svg\n")
cat(" - fig3_glmm.csv\n")
cat(
  "\nNote: check the reported outcome count against the intended ",
  "554 observations before finalizing the manuscript.\n",
  sep = ""
)


