# =============================================================================
# Fst / Dxy genome-wide distribution with QTL window averages
# Extracted from Population_genetics_vill_alle.qmd (windowed avgs section)
# =============================================================================
#
# Computes average f_dM, Dxy, and Fst in 1 cM windows (excluding NAs from means),
# plots histograms of genome-wide distributions, and overlays QTL credible-
# interval window averages as vertical lines.
#
# Requires: pixy_wide (from pixy import; columns: chrN, chromosome, window_pos_1,
#   window_pos_2, avg_dxy, avg_wc_fst, chr_position).
# f_dM: Dsuite Fstats .txt files are read in-script (path below).
# Paths: Dsuite dir, lepmap dir, all_QTL_CI CSV; ggsave outputs in working directory.
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# ---- Dsuite f_dM: load window stats and build pos_f_dM_clean ----
dsuite_dir <- "/Users/juliaharencar/Documents/Github/Population_genetics_vill_alle/Dsuite/Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120"
file_list <- list.files(path = dsuite_dir, pattern = "*.txt", full.names = TRUE)

LASIout_w50s25_all9 <- data.frame()
for (file in file_list) {
  temp_df <- read.delim(file, sep = "\t")
  LASIout_w50s25_all9 <- rbind(LASIout_w50s25_all9, temp_df)
}

pos_f_dM_chr <- LASIout_w50s25_all9 %>%
  select(chr, windowStart, windowEnd, f_dM) %>%
  mutate(f_dM_pos = ifelse(f_dM > 0, f_dM, NA_real_))

pos_f_dM_clean <- pos_f_dM_chr[is.finite(pos_f_dM_chr$f_dM_pos), ] %>%
  mutate(chrN = as.numeric(str_extract(chr, "\\d+")))

# ---- cM map: load lepmap and build bp -> cM interpolation ----
lepmap_dir <- "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/best_iter_ordermarkers/"
lepmap_files <- list.files(lepmap_dir,
  pattern = "*_clean_LodLim11_theta0.001_allF2s_0.5covCO_withRef.txt",
  full.names = TRUE
)

lepmap_df <- do.call(rbind, lapply(lepmap_files, function(f) {
  markers <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  names(markers) <- c("marker_num", "cM_pos", "chrom", "bp_pos")
  return(markers)
}))

bp_to_cm <- lepmap_df %>%
  filter(chrom != 4) %>%
  filter(!(chrom == 2 & bp_pos > 1.37e+08 & bp_pos < 1.4e+08)) %>%
  select(chrom, bp_pos, cM_pos)

chrom_invert <- c(2, 3, 7, 9)
bp_to_cm_fixed <- bp_to_cm %>%
  group_by(chrom) %>%
  mutate(cM_pos = ifelse(chrom %in% chrom_invert, max(cM_pos) - cM_pos, cM_pos)) %>%
  ungroup()

apply_iso <- function(df, bp_col = "bp_pos", cm_col = "cM_pos") {
  df <- df %>% arrange(.data[[bp_col]])
  y <- df[[cm_col]]
  iso <- isoreg(y)
  df$cM_iso <- iso$yf
  df
}

bp_to_cm_iso <- bp_to_cm_fixed %>%
  group_by(chrom) %>%
  group_modify(~ apply_iso(.x)) %>%
  ungroup()

# Optional: plot to check cM map
# ggplot(bp_to_cm_iso, aes(x = bp_pos)) +
#   geom_point(aes(y = cM_pos), size = 0.8, alpha = 0.4) +
#   geom_line(aes(y = cM_iso), color = "blue", size = 0.7) +
#   facet_wrap(~ chrom, scales = "free_x") +
#   labs(title = "Original cM (points) vs Isotonic monotone cM (blue line)", x = "bp_pos", y = "cM") +
#   theme_bw()

bp_to_cm_funs <- bp_to_cm_iso %>%
  group_by(chrom) %>%
  group_map(~ approxfun(.x$bp_pos, .x$cM_iso, rule = 2))
names(bp_to_cm_funs) <- unique(bp_to_cm_iso$chrom)

# ---- QTL and trait colors (shared) ----
qtl_ci_path <- "/Users/juliaharencar/Documents/Github/alle_vill_QTL/all_QTL_CIs_bp_cm.251130.csv"
all_QTL_CI <- read.csv(qtl_ci_path, header = TRUE)

all_QTL_CI <- all_QTL_CI %>%
  filter(!chr == 2 | !trait == "Nreabs") %>%
  filter(!chr == 5 | !trait == "GR")

trait_colors <- c(
  SD = "#0072B2",
  GR = "#CC79A7",
  tough = "#8250C4",
  Nreabs = "#E69F00",
  Chlreabs = "#009E73"
)

# ---- f_dM: 1 cM window averages and histogram ----
# Add cM position and drop chr 4 for genetic map consistency
pos_f_dM_clean <- pos_f_dM_clean %>%
  filter(chrN != 4) %>%
  rowwise() %>%
  mutate(cM_start = bp_to_cm_funs[[as.character(chrN)]](windowStart)) %>%
  ungroup()

# Per–1 cM bin: mean of f_dM_pos (positive f_dM only); NAs excluded in mean()
f_dM_summary <- pos_f_dM_clean %>%
  mutate(cM_bin = floor(cM_start)) %>%
  group_by(chrN, cM_bin) %>%
  summarise(
    avg_f_dM_per_cM = mean(f_dM_pos, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  )

# QTL interval average: mean f_dM in windows overlapping each QTL CI (NAs excluded)
qtl_f_dM_bp <- all_QTL_CI %>%
  rowwise() %>%
  mutate(
    avg_f_dM_QTL = mean(
      pos_f_dM_clean$f_dM[
        pos_f_dM_clean$chrN == chr &
          pos_f_dM_clean$windowStart <= upper95CI_bp &
          pos_f_dM_clean$windowEnd >= lower95CI_bp
      ],
      na.rm = TRUE
    )
  ) %>%
  ungroup()

f_dM_windows_hist_w.QTL.CIs <- ggplot() +
  geom_histogram(
    data = f_dM_summary,
    aes(x = avg_f_dM_per_cM, y = after_stat(density)),
    bins = 50,
    fill = "grey80",
    color = "black"
  ) +
  geom_vline(
    data = qtl_f_dM_bp,
    aes(xintercept = avg_f_dM_QTL, color = trait),
    linetype = "solid",
    linewidth = 1
  ) +
  scale_color_manual(values = trait_colors) +
  labs(x = "average f_dM", y = "frequency", color = "Trait") +
  theme_bw()

# ---- Dxy: 1 cM window averages and histogram ----
# Exclude NA and Chrom4; then mean per cM bin excludes NAs (na.rm = TRUE)
dxy <- pixy_wide %>%
  filter(is.finite(avg_dxy), chromosome != "Chrom4") %>%
  mutate(window_center = chr_position * 10e5) %>%
  select(chrN, chromosome, window_pos_1, window_pos_2, avg_dxy, window_center) %>%
  rowwise() %>%
  mutate(cM_window_start = bp_to_cm_funs[[as.character(chrN)]](window_pos_1)) %>%
  ungroup()

dxy_summary <- dxy %>%
  mutate(cM_bin = floor(cM_window_start)) %>%
  group_by(chrN, cM_bin) %>%
  summarise(
    avg_dxy_per_cM = mean(avg_dxy, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  )

# Re-load/filter QTL for Dxy (same as above; chr6 SD filter applied here as in qmd)
all_QTL_CI_dxy <- read.csv(qtl_ci_path, header = TRUE) %>%
  filter(!chr == 2 | !trait == "Nreabs") %>%
  filter(!chr == 6 | !trait == "SD") %>%
  filter(!chr == 5 | !trait == "GR")

qtl_dxy_bp <- all_QTL_CI_dxy %>%
  rowwise() %>%
  mutate(
    avg_dxy_QTL = mean(
      dxy$avg_dxy[
        dxy$chrN == chr &
          dxy$window_pos_1 <= upper95CI_bp &
          dxy$window_pos_2 >= lower95CI_bp
      ],
      na.rm = TRUE
    )
  ) %>%
  ungroup()

dxy_windows_hist_w.QTL.CIs <- ggplot() +
  geom_histogram(
    data = dxy_summary,
    aes(x = avg_dxy_per_cM, y = after_stat(density)),
    bins = 50,
    fill = "grey80",
    color = "black"
  ) +
  geom_vline(
    data = qtl_dxy_bp,
    aes(xintercept = avg_dxy_QTL, color = trait),
    linetype = "solid",
    linewidth = 1,
    show.legend = TRUE,
    key_glyph = "point"
  ) +
  scale_color_manual(values = trait_colors) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 15, size = 5, linetype = 0)
    )
  ) +
  labs(x = "average Dxy", y = "frequency", color = "Trait") +
  theme_bw() +
  theme(legend.position = "right", legend.justification = "center")

# ---- Fst: 1 cM window averages and histogram ----
# Exclude NA and Chrom4; mean per cM bin uses na.rm = TRUE
fst <- pixy_wide %>%
  filter(is.finite(avg_wc_fst), chromosome != "Chrom4") %>%
  mutate(window_center = chr_position * 10e5) %>%
  select(chrN, chromosome, window_pos_1, window_pos_2, avg_wc_fst, window_center) %>%
  rowwise() %>%
  mutate(cM_window_start = bp_to_cm_funs[[as.character(chrN)]](window_pos_1)) %>%
  ungroup()

fst_summary <- fst %>%
  mutate(cM_bin = floor(cM_window_start)) %>%
  group_by(chrN, cM_bin) %>%
  summarise(
    avg_fst_per_cM = mean(avg_wc_fst, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  )

all_QTL_CI_fst <- read.csv(qtl_ci_path, header = TRUE) %>%
  filter(!chr == 2 | !trait == "Nreabs") %>%
  filter(!chr == 6 | !trait == "SD") %>%
  filter(!chr == 5 | !trait == "GR")

trait_labels <- c(
  SD        = "Seed dormancy",
  GR        = "Growth rate",
  tough     = "Leaf toughness",
  Nreabs    = "Nitrogen reabsorption during drought",
  Chlreabs  = "Chlorophyll reabsorption during drought"
)

fst_cm <- fst %>%
  arrange(chrN, cM_window_start) %>%
  group_by(chrN) %>%
  mutate(cM_window_end = lead(cM_window_start)) %>%
  ungroup()

# QTL average in cM: mean Fst in windows overlapping CI (NAs excluded)
qtl_fst_cm <- all_QTL_CI_fst %>%
  rowwise() %>%
  mutate(
    avg_fst_QTL = mean(
      fst_cm$avg_wc_fst[
        fst_cm$chrN == chr &
          fst_cm$cM_window_start <= upper95CI_cm &
          fst_cm$cM_window_end >= lower95CI_cm
      ],
      na.rm = TRUE
    )
  ) %>%
  ungroup()

fst_windows_hist_w.QTL.CIs <- ggplot() +
  geom_histogram(
    data = fst_summary,
    aes(x = avg_fst_per_cM, y = after_stat(density)),
    binwidth = 0.01,
    fill = "grey80",
    color = "black"
  ) +
  geom_vline(
    data = qtl_fst_cm,
    aes(xintercept = avg_fst_QTL, color = trait),
    linetype = "solid",
    linewidth = 1,
    show.legend = TRUE,
    key_glyph = "point"
  ) +
  scale_color_manual(values = trait_colors, labels = trait_labels) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 15, size = 5, linetype = 0)
    )
  ) +
  labs(
    x = "average differentiation (Fst in 1 cM windows)",
    y = "probability density",
    color = "QTL window averages\nTrait:"
  ) +
  theme_bw() +
  theme(legend.position = "right", legend.justification = "center")

# ---- Figure outputs ----
# 1) Solo figure: Fst genome-wide distribution and Fst values from QTL windows
ggsave("Fst_genomewide_distribution_w_QTL_windows.png", fst_windows_hist_w.QTL.CIs,
  device = "png", width = 17, height = 6.5, units = "cm"
)

# 2) Three-panel figure: Fst, Dxy, f_dM stacked (top to bottom), shared legend
three_panel_hist <- fst_windows_hist_w.QTL.CIs / dxy_windows_hist_w.QTL.CIs / f_dM_windows_hist_w.QTL.CIs +
  plot_layout(guides = "collect")
ggsave("Fst_Dxy_fdM_genomewide_distribution_w_QTL_3panel.png", three_panel_hist,
  device = "png", width = 17, height = 18, units = "cm"
)
