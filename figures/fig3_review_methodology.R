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
  "Synthesis method",     "Quantitative",  6, "m_q",   BLUE,   "",
  "Synthesis method",     "Narrative only",5, "m_n",   ORANGE, plab("method"),
  "Trait definitions",    "Explicit",      4, "t_y",   BLUE,   "",
  "Trait definitions",    "None",          3, "t_n",   ORANGE, plab("trait_def"),
  "Ecological framework", "Explicit",      2, "e_y",   BLUE,   "",
  "Ecological framework", "None",          1, "e_n",   ORANGE, plab("eco"))
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
  annotate("text", x = .58, y = 6.9, label = "← more rigorous, less positive",
           color = BLUE, fontface = 2, size = 3.2) +
  annotate("text", x = .93, y = 6.9, label = "less rigorous, more positive →",
           color = ORANGE, fontface = 2, size = 3.2) +
  scale_x_continuous(limits = c(.24, 1.02), breaks = seq(.5,1,.1),
                     labels = paste0(seq(50,100,10), "%"), position = "top") +
  scale_y_continuous(limits = c(.3, 7.2)) +
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

