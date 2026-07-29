# ============================================================================
#  fig2_review_disagreement.R  —  Section-2 figure (R / ggplot2)
#     "The same relationship, opposite verdicts across reviews"
#
#  The human-readable face of the between-review result: for well-covered
#  trait-NCP relationships, each dot is ONE systematic review's overall verdict
#  (negative / no clear direction / positive). Reviews contradict one another
#  pair by pair. (Full variance decomposition + robustness battery: text /
#  Supplementary.)  Twin of figures/fig2_review_disagreement.py.
#
#  INPUT   data/data_model.csv
#  OUTPUT  outputs/fig2_review_disagreement.(png|svg)
#  RUN     Rscript figures/fig2_review_disagreement.R
#  DEPS    dplyr, tidyr, ggplot2. SVG needs svglite (optional).
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

MIN_REVIEWS <- 5L
N_PAIRS     <- 10L

root <- tryCatch(here::here(), error = function(e) getwd())
dm_path <- file.path(root, "data", "data_model.csv")
if (!file.exists(dm_path)) stop("Missing ", dm_path)
dm <- read.csv(dm_path, stringsAsFactors = FALSE)
for (k in c("Npos", "Nneg", "Nnull"))
  dm[[k]] <- suppressWarnings(as.numeric(dm[[k]])); dm[is.na(dm)] <- 0

SHORT_TRAIT <- c("specific leaf area"="Specific leaf area","leaf nitrogen content"="Leaf nitrogen",
 "leaf dry matter content"="Leaf dry-matter content","height"="Plant height","wood density"="Wood density",
 "root diameter"="Root diameter","aboveground biomass"="Aboveground biomass","leaf area"="Leaf area",
 "leaf phosphorus content"="Leaf phosphorus","crown size"="Crown size","leaf carbon content"="Leaf carbon")
SHORT_NCP <- c("regulation of climate"="climate regulation","formation and protection of soils"="soil formation",
 "regulation of freshwater quantity"="freshwater regulation","food and feed"="food & feed",
 "regulation of hazards and extreme events"="hazard regulation","regulation of detrimental organisms"="pest regulation")
nice <- function(ts) {
  parts <- strsplit(ts, "__", fixed = TRUE)[[1]]
  t <- ifelse(parts[1] %in% names(SHORT_TRAIT), SHORT_TRAIT[parts[1]],
              paste0(toupper(substring(parts[1],1,1)), substring(parts[1],2)))
  n <- ifelse(parts[2] %in% names(SHORT_NCP), SHORT_NCP[parts[2]], parts[2])
  paste0(t, "  ->  ", n)
}

# ---- one direction per UNIQUE review per pair -------------------------------
per_rev <- dm %>%
  group_by(ts, source_id) %>%
  summarise(p = sum(Npos), n = sum(Nneg), z = sum(Nnull), .groups = "drop") %>%
  mutate(dir = case_when(p > n ~ "pos", n > p ~ "neg", TRUE ~ "null"))

pair_sum <- per_rev %>%
  group_by(ts) %>%
  summarise(nrev = n(),
            npos = sum(dir == "pos"), nneg = sum(dir == "neg"),
            nnull = sum(dir == "null"),
            disagree = pmin(npos, nneg), .groups = "drop") %>%
  filter(nrev >= MIN_REVIEWS) %>%
  arrange(desc(disagree), desc(nrev)) %>%
  slice_head(n = N_PAIRS)

# ---- explicit dot coordinates (beeswarm-style, 3 columns) -------------------
col_x <- c(neg = 1, null = 2, pos = 3)
per   <- 6; gap <- 0.13
dots <- per_rev %>%
  filter(ts %in% pair_sum$ts) %>%
  mutate(pair = factor(ts, levels = rev(pair_sum$ts))) %>%
  group_by(pair, dir) %>%
  mutate(k = row_number() - 1,
         line = k %/% per, idx = k %% per,
         inrow = pmin(n() - line*per, per),
         nlines = ceiling(n()/per),
         x = col_x[dir] + (idx - (inrow - 1)/2) * gap,
         y = as.integer(pair) + (line - (nlines - 1)/2) * 0.16) %>%
  ungroup()

lab_df <- pair_sum %>% mutate(pair = factor(ts, levels = rev(pair_sum$ts)),
                              y = as.integer(pair),
                              label = vapply(ts, nice, character(1)))

pal <- c(neg = "#d1495b", null = "#b9bcc2", pos = "#2a78d6")
ny  <- nrow(pair_sum)

p <- ggplot() +
  geom_vline(xintercept = c(1.5, 2.5), color = "#ecebe6", linewidth = 0.5) +
  geom_point(data = dots, aes(x, y, color = dir), size = 3.2,
             stroke = 0.5, show.legend = FALSE) +
  scale_color_manual(values = pal) +
  # column headers
  annotate("text", x = 1, y = ny + 0.75, label = "reported NEGATIVE",
           color = pal["neg"], fontface = 2, size = 3.4) +
  annotate("text", x = 2, y = ny + 0.75, label = "no clear direction",
           color = "#8f9298", fontface = 2, size = 3.4) +
  annotate("text", x = 3, y = ny + 0.75, label = "reported POSITIVE",
           color = pal["pos"], fontface = 2, size = 3.4) +
  # pair labels
  geom_text(data = lab_df, aes(x = 0.45, y = y, label = label),
            hjust = 1, fontface = 2, size = 3.6, color = "#0b0b0b") +
  geom_text(data = lab_df, aes(x = 0.45, y = y - 0.28,
            label = paste0(nrev, " reviews")),
            hjust = 1, size = 2.7, color = "#52514e") +
  scale_x_continuous(limits = c(-1.7, 3.6)) +
  scale_y_continuous(limits = c(0.3, ny + 1.1)) +
  # labs(title = "The same relationship, opposite verdicts across reviews",
  #      subtitle = paste0("Each dot is one systematic review's overall conclusion for that ",
  #                        "trait-NCP relationship (well-covered pairs, >=", MIN_REVIEWS, " reviews)."),
  #      caption = paste0("Reviews contradict one another pair by pair: the reported direction depends\n",
  #                       "more on which review reports it than on the trait or the service involved.")) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0),
        plot.subtitle = element_text(size = 10, color = "#52514e", hjust = 0,
                                     margin = margin(b = 8)),
        plot.caption = element_text(face = "bold", size = 9.5, color = "#0b0b0b",
                                    hjust = 0, margin = margin(t = 10)),
        plot.margin = margin(16, 20, 14, 20),
        plot.background = element_rect(fill = "#ffffff", color = NA))

out <- file.path(root, "outputs", "fig2_review_disagreement")
ggsave(paste0(out, ".png"), p, width = 10.5, height = 2.6 + ny*0.52, dpi = 200, bg = "#fcfcfb")
ok_svg <- requireNamespace("svglite", quietly = TRUE)
if (ok_svg) ggsave(paste0(out, ".svg"), p, width = 10.5, height = 2.6 + ny*0.52, bg = "#fcfcfb")
message("wrote ", out, ".png", if (ok_svg) " (+ .svg)" else " (install 'svglite' for .svg)",
        "  [", ny, " pairs]")

