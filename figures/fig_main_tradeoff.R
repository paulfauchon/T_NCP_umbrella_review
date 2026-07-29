# ============================================================================
#  fig_main_tradeoff.R  —  Main figure: "Same trait, opposite contributions"
#                          (R / ggplot2 version — same figure as the .py)
#
#  Trait x NCP posterior effect-direction matrix from the PAIR-LEVEL (simple)
#  model (no cross-NCP borrowing) — the correct parameterisation for reading
#  per-pair directions. Fill = continuous P(effect > 0) on a diverging RdBu
#  ramp; cell number = primary studies. Two conventional reading aids overlaid
#  on the continuous scale:
#     bold black outline  = strong      (P > 0.90 or < 0.10)
#     thin grey outline   = suggestive  (P > 0.75 or < 0.25)
#  Cells with < 3 primary studies are faded (under-evidenced).
#
#  INPUT   outputs/reversal_check.csv   (from R/check_reversals.R)
#          data/data_model.csv          (trait/NCP hierarchy)
#  OUTPUT  outputs/fig_main_tradeoff.(png|svg)     (13 curated traits)
#          outputs/fig_supp_tradeoff.(png|svg)     (19 traits, wider set)
#
#  RUN     Rscript figures/fig_main_tradeoff.R      (or source() it)
#  DEPS    ggplot2, dplyr, tidyr, tibble, forcats. SVG needs svglite (optional).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(forcats); library(ggplot2)
})

# ---- Resolve pipeline paths (script lives in <root>/figures/) --------------
root <- tryCatch(here::here(), error = function(e) getwd())
OUT_DIR  <- file.path(root, "outputs")
DATA_DIR <- file.path(root, "data")
rev_path <- file.path(OUT_DIR, "reversal_check.csv")
if (!file.exists(rev_path))
  stop("Missing ", rev_path, " — run R/check_reversals.R first.")

rev <- read.csv(rev_path, stringsAsFactors = FALSE)

# ---- Trait rows: short label -> trait_grouped, in super-groups --------------
main_traits <- tribble(
  ~group,                      ~label,           ~trait,
  "Leaf economics & size",     "SLA",            "specific leaf area",
  "Leaf economics & size",     "LDMC",           "leaf dry matter content",
  "Leaf economics & size",     "Leaf N",         "leaf nitrogen content",
  "Leaf economics & size",     "Leaf P",         "leaf phosphorus content",
  "Leaf economics & size",     "Leaf area",      "leaf area",
  "Leaf economics & size",     "LAI",            "leaf area index",
  "Whole-plant",  "Height",         "height",
  "Whole-plant",  "Crown size",     "crown size",
  "Whole-plant",  "Wood density",   "wood density",
  "Whole-plant",  "Foliage density","foliage density",
  "Root",                      "Root diameter",  "root diameter",
  "Root",                      "Root length",    "root length",
  "Root",                      "Rooting depth",  "rooting depth")

supp_traits <- tribble(
  ~group,                      ~label,            ~trait,
  "Leaf economics & size",     "SLA",             "specific leaf area",
  "Leaf economics & size",     "LDMC",            "leaf dry matter content",
  "Leaf economics & size",     "Leaf N",          "leaf nitrogen content",
  "Leaf economics & size",     "Leaf P",          "leaf phosphorus content",
  "Leaf economics & size",     "Leaf C",          "leaf carbon content",
  "Leaf economics & size",     "Leaf area",       "leaf area",
  "Leaf economics & size",     "LAI",             "leaf area index",
  "Leaf economics & size",     "Leaf thickness",  "leaf thickness",
  "Leaf economics & size",     "Leaf toughness",  "leaf toughness",
  "Whole-plant",  "Height",          "height",
  "Whole-plant",  "Crown size",      "crown size",
  "Whole-plant",  "Wood density",    "wood density",
  "Whole-plant",  "Foliage density", "foliage density",
  "Whole-plant",  "Aboveg. biomass", "aboveground biomass",
  "Whole-plant",  "Aboveg. litter",  "aboveground litter",
  "Root",                      "Root diameter",   "root diameter",
  "Root",                      "Root length",     "root length",
  "Root",                      "Rooting depth",   "rooting depth",
  "Reprod.",                   "Seed mass",       "seed mass")

# ---- Columns: short -> full NCP name, grouped by NCP category ----------------
ncp_tab <- tribble(
  ~cat,          ~label,            ~ncp,
  "REGULATION",  "Soil formation",  "formation and protection of soils",
  "REGULATION",  "Climate",         "regulation of climate",
  "REGULATION",  "Freshwater qty",  "regulation of freshwater quantity",
  "REGULATION",  "Freshwater qual", "regulation of freshwater quality",
  "REGULATION",  "Hazards",         "regulation of hazards and extreme events",
  "REGULATION",  "Pest reg.",       "regulation of detrimental organisms",
  "REGULATION",  "Air quality",     "regulation of air quality",
  "REGULATION",  "Pollination",     "pollination and seed dispersal",
  "REGULATION",  "Habitat",         "habitat creation and maintenance",
  "MATERIAL",    "Food & feed",     "food and feed",
  "MATERIAL",    "Materials",       "materials and assistance",
  "MATERIAL",    "Medicinal",       "medicinal and genetic ressources",
  "NON-MAT.",    "Experiences",     "experiences")

# Diverging RdBu ramp (low P = negative = red, high P = positive = blue)
RDBU <- c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7",
          "#d1e5f0","#92c5de","#4393c3","#2166ac")

make_fig <- function(traits, out_base, title) {
  ncp2 <- ncp_tab %>% rename(ncp_label = label)   # avoid clash with traits$label
  d <- tidyr::crossing(traits, ncp2) %>%
    left_join(rev %>% transmute(trait = trait_grouped, ncp,
                                P = P_simple, n = n_studies),
              by = c("trait", "ncp")) %>%
    mutate(
      documented = !is.na(P) & n >= 3,
      Pf   = ifelse(documented, P, NA_real_),
      tier = case_when(!documented                 ~ "none",
                       P > 0.90 | P < 0.10          ~ "strong",
                       P > 0.75 | P < 0.25          ~ "sugg",
                       TRUE                         ~ "weak"),
      nlab = ifelse(documented, as.character(round(n)), ""),
      tcol = ifelse(documented & (P < 0.16 | P > 0.84), "white", "#52514e"),
      label     = factor(label,     levels = rev(traits$label)),        # top -> bottom
      ncp_label = factor(ncp_label, levels = ncp_tab$label),            # left -> right
      ncat      = factor(cat,       levels = unique(ncp_tab$cat)),
      grp       = factor(group,     levels = unique(traits$group)))

  p <- ggplot(d, aes(ncp_label, label)) +
    geom_tile(aes(fill = Pf), color = "white", linewidth = 0.4,
              width = 0.94, height = 0.94) +
    geom_tile(data = dplyr::filter(d, tier == "sugg"),
              fill = NA, color = "#6f6f6f", linewidth = 0.4,
              width = 0.94, height = 0.94) +
    geom_tile(data = dplyr::filter(d, tier == "strong"),
              fill = NA, color = "#1a1a1a", linewidth = 0.9,
              width = 0.94, height = 0.94) +
    geom_text(aes(label = nlab, color = tcol), size = 2.5) +
    scale_color_identity() +
    scale_fill_gradientn(colours = RDBU, limits = c(0, 1),
                         na.value = "#f7f6f3", name = "P(effect > 0)",
                         breaks = c(0, 0.5, 1),
                         labels = c("negative", "uncertain", "positive")) +
    scale_x_discrete(position = "top") +
    facet_grid(grp ~ ncat, scales = "free", space = "free", switch = "y") +
    labs(title = title,
         caption = paste0(
           "Pair-level model. Fill = continuous P(effect>0). ",
           "Bold outline = strong (P>0.90); thin outline = suggestive (P>0.75); ",
           "faded = <3 studies. Cell number = primary studies.")) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title = element_blank(),
      axis.text.x.top = element_text(angle = 45, hjust = 0, size = 8, color = "#52514e"),
      axis.text.y = element_text(size = 9, face = "bold", color = "#0b0b0b"),
      panel.grid = element_blank(),
      panel.spacing = unit(3, "pt"),
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 8.5, color = "#52514e"),
      strip.text.x = element_text(face = "bold", size = 9, color = "#52514e"),
      strip.placement = "outside",
      plot.title = element_text(face = "bold", size = 15),
      plot.caption = element_text(hjust = 0, size = 7.5, color = "#898781"),
      legend.position = "bottom", legend.key.width = unit(28, "pt"),
      plot.background = element_rect(fill = "#ffffff", color = NA))

  nrow_t <- nrow(traits)
  ggsave(paste0(out_base, ".png"), p, width = 11, height = 2.4 + nrow_t * 0.34,
         dpi = 200, bg = "#ffffff")
  ok_svg <- requireNamespace("svglite", quietly = TRUE)
  if (ok_svg)
    ggsave(paste0(out_base, ".svg"), p, width = 11, height = 2.4 + nrow_t * 0.34,
           bg = "#ffffff")
  message("wrote ", out_base, ".png",
          if (ok_svg) " (+ .svg)" else " (install 'svglite' for .svg)",
          "  [", nrow_t, " traits]")
}

make_fig(main_traits, file.path(OUT_DIR, "fig_main_tradeoff"),
         "Same trait, opposite contributions")
make_fig(supp_traits, file.path(OUT_DIR, "fig_supp_tradeoff"),
         "Same trait, opposite contributions (extended set)")

