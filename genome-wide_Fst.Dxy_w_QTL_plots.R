# =============================================================================
# genome-wide_Fst.Dxy_w_QTL_plots.R
# Genome-wide Fst, Dxy, and f_dM with QTL locations and distribution histograms
# =============================================================================
#
# Single script: (1) imports pixy + Dsuite f_dM once, (2) along-genome Fst/Dxy
# plots (500 kb window averages + QTL peaks/CIs), (3) 1 cM windowed distribution
# histograms (Fst, Dxy, f_dM) with QTL window averages, density curves, tail shading.
#
# Dependencies: tidyverse, ggplot2, patchwork. Paths below for pixy, Dsuite,
# QTL CSV, lepmap; ggsave outputs in "pixy/", "Dsuite/", and working directory.
# =============================================================================

# ---- load dependencies ----
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

theme_set(theme_bw(base_size = 11))

# ---- pixy: import data ----
# function from the pixy website (https://pixy.readthedocs.io/en/latest/plotting.html)
pixy_to_long <- function(pixy_files){

  pixy_df <- list()

  for(i in seq_along(pixy_files)){

    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])

    if(stat_file_type == "pi"){

      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)

      pixy_df[[i]] <- df


    } else{

      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df

    }

  }

  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)

}

# Paths (edit as needed)
# Project root: where pixy/ and Dsuite/ output folders live (script's repo)
project_root <- "/Users/juliaharencar/Documents/Github/Costus_popgen_clean"
pixy_folder <- "/Users/juliaharencar/Documents/Github/Population_genetics_vill_alle/pixy/output_251124_LASIref"
dsuite_dir <- "/Users/juliaharencar/Documents/Github/Population_genetics_vill_alle/Dsuite/Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120"
qtl_ci_path <- "/Users/juliaharencar/Documents/Github/alle_vill_QTL/all_QTL_CIs_bp_cm.251130.csv"
lepmap_dir <- "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/best_iter_ordermarkers/"

# ---- Pixy: import ----
pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)
pixy_wide <- pivot_wider(pixy_df, names_from = statistic, values_from = value) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000)
# make chromosome number only col:
pixy_wide <- pixy_wide %>%
  mutate(chrN = as.numeric(str_extract(chromosome, "\\d+")))

# ---- Dsuite f_dM: load window stats ----
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

# ---- Genome-wide f_dM figure (positive values only, 99.9th % cutoff) ----
data_cum <- LASIout_w50s25_all9 %>%
  group_by(chr) %>%
  summarise(max_windowStart = max(windowStart)) %>%
  mutate(windowStart_add = lag(cumsum(max_windowStart), default = 0)) %>%
  select(chr, windowStart_add)

LASIout_w50s25_all9 <- LASIout_w50s25_all9 %>%
  inner_join(data_cum, by = "chr") %>%
  mutate(windowStart_cum = windowStart + windowStart_add,
         chr_label = gsub("Chrom", "", chr))

axis_set_dsuite <- LASIout_w50s25_all9 %>%
  group_by(chr) %>%
  summarize(center = mean(windowStart_cum), chr_label = first(chr_label))

cutoff99.9thPercentile <- quantile(abs(LASIout_w50s25_all9$f_dM[LASIout_w50s25_all9$f_dM < 0]), 0.999)

LASIout_positive <- LASIout_w50s25_all9 %>% filter(f_dM > 0)

LASI_99.9thper_w50s25_f_dM <- ggplot(data = LASIout_positive, aes(x = windowStart_cum, y = f_dM, color = as_factor(chr))) +
  geom_point(size = 3) +
  coord_cartesian(xlim = c(24398, 1120430367)) +
  geom_hline(yintercept = cutoff99.9thPercentile, linetype = "dashed", color = "red") +
  scale_color_manual(values = rep(c("black", "darkgray"), ceiling(length(unique(axis_set_dsuite$chr)) / 2))) +
  guides(color = "none") +
  geom_point(data = subset(LASIout_positive, f_dM > cutoff99.9thPercentile), color = "red") +
  scale_x_continuous(label = axis_set_dsuite$chr_label, breaks = axis_set_dsuite$center) +
  labs(x = NULL, y = "f_dM") +
  theme(plot.margin = unit(c(2, 1, 1, 1), "cm"))

if (!dir.exists(file.path(project_root, "Dsuite"))) dir.create(file.path(project_root, "Dsuite"))
ggsave(file.path(project_root, "Dsuite", "LASI_99.9thper_w50s25_f_dM_pos_only_all_chroms.png"),
       plot = LASI_99.9thper_w50s25_f_dM,
       width = 24, height = 5.5, units = "in", dpi = 300)

# ---- QTL and trait colors (used by along-genome plots and histograms) ----
all_QTL_CI <- read.csv(qtl_ci_path, header = TRUE) %>%
  filter(!chr == 2 | !trait == "Nreabs") %>%
  filter(!chr == 5 | !trait == "GR")
trait_colors <- c(
  SD = "#0072B2",
  GR = "#CC79A7",
  tough = "#8250C4",
  Nreabs = "#E69F00",
  Chlreabs = "#009E73"
)
all_QTL_CI_dxy <- all_QTL_CI %>% filter(!chr == 6 | !trait == "SD")
all_QTL_CI_fst <- all_QTL_CI %>% filter(!chr == 6 | !trait == "SD")

# ---- cM map (for 1 cM windowed histograms) ----
lepmap_files <- list.files(lepmap_dir, pattern = "*_clean_LodLim11_theta0.001_allF2s_0.5covCO_withRef.txt", full.names = TRUE)
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
bp_to_cm_funs <- bp_to_cm_iso %>%
  group_by(chrom) %>%
  group_map(~ approxfun(.x$bp_pos, .x$cM_iso, rule = 2))
names(bp_to_cm_funs) <- unique(bp_to_cm_iso$chrom)

# ---- Along-genome: 500 kb window averages and QTL ----
# One row per 500 kb bin per chromosome: WindowMid (bp), MeanY (mean of raw values).
# Bins with >50% of 10 kb values missing are excluded (no line segment there).
window_size_bp <- 500000
min_frac_valid <- 0.25  # require at least 25% of 10 kb windows in bin to have valid value (raise to 0.5 for stricter)

windowed_Dxy_avg <- pixy_wide %>%
  mutate(bin_mid = (floor(window_pos_1 / window_size_bp) + 0.5) * window_size_bp) %>%
  group_by(chrN, WindowMid = bin_mid) %>%
  summarise(
    n_total = n(),
    n_valid = sum(is.finite(avg_dxy)),
    MeanY = mean(avg_dxy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_total > 0, n_valid / n_total >= min_frac_valid) %>%
  select(chrN, WindowMid, MeanY)

windowed_Fst_avg <- pixy_wide %>%
  mutate(bin_mid = (floor(window_pos_1 / window_size_bp) + 0.5) * window_size_bp) %>%
  group_by(chrN, WindowMid = bin_mid) %>%
  summarise(
    n_total = n(),
    n_valid = sum(is.finite(avg_wc_fst)),
    MeanY = mean(avg_wc_fst, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_total > 0, n_valid / n_total >= min_frac_valid) %>%
  select(chrN, WindowMid, MeanY)

# Split by chromosome for plotting (list with chr1, ..., chr9; empty tibble if no passing bins)
windowed_Dxy_avg_list <- lapply(1:9, function(i) windowed_Dxy_avg %>% filter(chrN == i))
names(windowed_Dxy_avg_list) <- paste0("chr", 1:9)
windowed_Fst_avg_list <- lapply(1:9, function(i) windowed_Fst_avg %>% filter(chrN == i))
names(windowed_Fst_avg_list) <- paste0("chr", 1:9)

# ---- Dxy: plot with QTL ----
# Scale bar: 2 Mb at top left, 1 Mb black + 1 Mb white
scale_bar_2Mb <- function(raw_df, win_df, y_col = "avg_dxy") {
  x_range <- range(c(raw_df$window_pos_1, win_df$WindowMid), na.rm = TRUE)
  y_vals <- c(0, pmax(raw_df[[y_col]], 0, na.rm = TRUE), win_df$MeanY)
  y_range <- range(y_vals, na.rm = TRUE)
  bar_x_left <- x_range[1] + 0.01 * diff(x_range)
  bar_x_mid <- bar_x_left + 1e6
  bar_x_right <- bar_x_left + 2e6
  bar_ht <- 0.05 * diff(y_range)
  bar_y_top <- y_range[2] - 0.02 * diff(y_range)
  bar_y_bottom <- bar_y_top - bar_ht
  scale_df <- data.frame(
    xmin = c(bar_x_left, bar_x_mid),
    xmax = c(bar_x_mid, bar_x_right),
    ymin = bar_y_bottom,
    ymax = bar_y_top,
    fill = c("black", "white")
  )
  list(
    geom_rect(data = scale_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              colour = "black", linewidth = 0.4, inherit.aes = FALSE),
    scale_fill_identity(),
    annotate("text", x = bar_x_left + 1e6, y = bar_y_bottom - 0.015 * diff(y_range),
             label = "2 Mb", size = 3.2, hjust = 0.5)
  )
}

plot_chr_combined <- function(raw_df, win_df) {
  sb <- scale_bar_2Mb(raw_df, win_df, "avg_dxy")
  ggplot(raw_df, aes(x = window_pos_1)) +
    geom_point(
      aes(y = pmax(avg_dxy, 0)),
      size = 1,
      alpha = 0.4,
      color = "grey40"
    ) +
    geom_line(
      data = win_df,
      aes(x = WindowMid, y = MeanY),
      color = "darkred",
      linewidth = 0.7
    ) +
    sb +
    labs(
      x = "Position (bp)",
      y = "Dxy"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)))
}

chroms <- 1:9

plots <- lapply(chroms, function(i) {

  raw_df_i <- pixy_wide %>% filter(chrN == i, is.finite(avg_dxy))
  win_df_i <- windowed_Dxy_avg_list[[paste0("chr", i)]]

  plot_chr_combined(
    raw_df = raw_df_i,
    win_df = win_df_i
  )
})
names(plots) <- paste0("chr", chroms)

chrom_list <- sort(unique(pixy_wide$chrN))

add_qtl_layers <- function(base_plot, qtl_rows, trait_colors) {

  if (nrow(qtl_rows) == 0) return(base_plot)

  base_plot +
    geom_rect(
      data = qtl_rows,
      aes(
        xmin = lower95CI_bp,
        xmax = upper95CI_bp,
        ymin = -Inf,
        ymax = Inf,
        fill = trait
      ),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    geom_vline(
      data = qtl_rows,
      aes(xintercept = peak_bp, color = trait),
      linewidth = 1
    ) +
    scale_fill_manual(values = trait_colors) +
    scale_color_manual(values = trait_colors)
}

final_plots <- list()

for (i in seq_along(chrom_list)) {

  chr_num <- chrom_list[i]
  base_plot <- plots[[i]]

  qtl_here <- all_QTL_CI %>%
    filter(chr == chr_num)

  plot_i <- add_qtl_layers(base_plot, qtl_here, trait_colors)
  plot_i <- plot_i + ggtitle(paste("Chromosome", chr_num))
  final_plots[[i]] <- plot_i
}
names(final_plots) <- paste0("Chrom", chrom_list)

# 2-panel figure: chromosome 1 above chromosome 3 (Dxy)
Dxy_chr1_chr3_2panel <- (final_plots$Chrom1 / final_plots$Chrom3) &
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
ggsave(file.path(project_root, "pixy", "Dxy_chr1_chr3_2panel.png"), Dxy_chr1_chr3_2panel,
       device = "png", width = 22, height = 10, units = "cm")

Dxy_LASI_sep_chrs_w.QTL <- final_plots$Chrom1 + final_plots$Chrom2 + final_plots$Chrom3 +
  final_plots$Chrom5 + final_plots$Chrom6 + final_plots$Chrom7 + final_plots$Chrom8 + final_plots$Chrom9

ggsave(file.path(project_root, "pixy", "Dxy_10kb_LASIout_separate_chrs_w.500kb.avg_7QTL_narrowed.CIs.png"),
       Dxy_LASI_sep_chrs_w.QTL, device = "png", width = 50, height = 20, units = "cm")

# ---- Fst: plot with QTL ----
plot_chr_combined <- function(raw_df, win_df) {
  sb <- scale_bar_2Mb(raw_df, win_df, "avg_wc_fst")
  ggplot(raw_df, aes(x = window_pos_1)) +
    geom_point(
      aes(y = pmax(avg_wc_fst, 0)),
      size = 1,
      alpha = 0.4,
      color = "grey40"
    ) +
    geom_line(
      data = win_df,
      aes(x = WindowMid, y = MeanY),
      color = "darkred",
      linewidth = 0.7
    ) +
    sb +
    labs(
      x = "Position (bp)",
      y = "Fst"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)))
}

plots <- lapply(chroms, function(i) {

  raw_df_i <- pixy_wide %>% filter(chrN == i, is.finite(avg_wc_fst))
  win_df_i <- windowed_Fst_avg_list[[paste0("chr", i)]]

  plot_chr_combined(
    raw_df = raw_df_i,
    win_df = win_df_i
  )
})
names(plots) <- paste0("chr", chroms)

# Reuse QTL data and helper (already in env)
final_plots <- list()

for (i in seq_along(chrom_list)) {

  chr_num <- chrom_list[i]
  base_plot <- plots[[i]]

  qtl_here <- all_QTL_CI %>%
    filter(chr == chr_num)

  plot_i <- add_qtl_layers(base_plot, qtl_here, trait_colors)
  plot_i <- plot_i + ggtitle(paste("Chromosome", chr_num))
  final_plots[[i]] <- plot_i
}
names(final_plots) <- paste0("Chrom", chrom_list)

# 2-panel figure: chromosome 1 above chromosome 3 (Fst)
Fst_chr1_chr3_2panel <- (final_plots$Chrom1 / final_plots$Chrom3) &
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
ggsave(file.path(project_root, "pixy", "Fst_chr1_chr3_2panel.png"), Fst_chr1_chr3_2panel,
       device = "png", width = 22, height = 10, units = "cm")

Fst_LASI_sep_chrs_w.QTL <- final_plots$Chrom1 + final_plots$Chrom2 + final_plots$Chrom3 +
  final_plots$Chrom5 + final_plots$Chrom6 + final_plots$Chrom7 + final_plots$Chrom8 + final_plots$Chrom9

ggsave(file.path(project_root, "pixy", "Fst_10kb_LASIout_separate_chrs_w.500kb.avg_7QTL_narrowed.CIs.png"),
       Fst_LASI_sep_chrs_w.QTL, device = "png", width = 50, height = 20, units = "cm")

# ---- Genome-wide distribution: 1 cM window histograms with QTL averages ----
# f_dM: add cM position (drop chr 4), 1 cM bins, QTL window averages
pos_f_dM_clean <- pos_f_dM_clean %>%
  filter(chrN != 4) %>%
  rowwise() %>%
  mutate(cM_start = bp_to_cm_funs[[as.character(chrN)]](windowStart)) %>%
  ungroup()
f_dM_summary <- pos_f_dM_clean %>%
  mutate(cM_bin = floor(cM_start)) %>%
  group_by(chrN, cM_bin) %>%
  summarise(avg_f_dM_per_cM = mean(f_dM_pos, na.rm = TRUE), n_windows = n(), .groups = "drop")
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
# Density and 5th/95th percentile tails for f_dM (shade outer 5% each tail)
d_fdM <- density(f_dM_summary$avg_f_dM_per_cM)
p05_fdM <- quantile(f_dM_summary$avg_f_dM_per_cM, 0.05)
p95_fdM <- quantile(f_dM_summary$avg_f_dM_per_cM, 0.95)
tail_left_fdM <- data.frame(x = d_fdM$x, y = d_fdM$y) %>% filter(x <= p05_fdM)
tail_right_fdM <- data.frame(x = d_fdM$x, y = d_fdM$y) %>% filter(x >= p95_fdM)

f_dM_windows_hist_w.QTL.CIs <- ggplot() +
  geom_histogram(data = f_dM_summary, aes(x = avg_f_dM_per_cM, y = after_stat(density)), bins = 50, fill = "grey80", color = "black") +
  geom_ribbon(data = tail_left_fdM, aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.35) +
  geom_ribbon(data = tail_right_fdM, aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.35) +
  geom_density(data = f_dM_summary, aes(x = avg_f_dM_per_cM, y = after_stat(density)), color = "darkblue", linewidth = 0.8) +
  geom_vline(data = qtl_f_dM_bp, aes(xintercept = avg_f_dM_QTL, color = trait), linetype = "solid", linewidth = 1) +
  scale_color_manual(values = trait_colors) +
  labs(x = "average f_dM", y = "frequency", color = "Trait")

# Dxy: 1 cM bins, QTL window averages (all_QTL_CI_dxy)
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
  summarise(avg_dxy_per_cM = mean(avg_dxy, na.rm = TRUE), n_windows = n(), .groups = "drop")
qtl_dxy_bp <- all_QTL_CI_dxy %>%
  rowwise() %>%
  mutate(
    avg_dxy_QTL = mean(
      dxy$avg_dxy[dxy$chrN == chr & dxy$window_pos_1 <= upper95CI_bp & dxy$window_pos_2 >= lower95CI_bp],
      na.rm = TRUE
    )
  ) %>%
  ungroup()
# Density and 5th/95th percentile tails for Dxy (shade outer 5% each tail)
d_dxy <- density(dxy_summary$avg_dxy_per_cM)
p05_dxy <- quantile(dxy_summary$avg_dxy_per_cM, 0.05)
p95_dxy <- quantile(dxy_summary$avg_dxy_per_cM, 0.95)
tail_left_dxy <- data.frame(x = d_dxy$x, y = d_dxy$y) %>% filter(x <= p05_dxy)
tail_right_dxy <- data.frame(x = d_dxy$x, y = d_dxy$y) %>% filter(x >= p95_dxy)

dxy_windows_hist_w.QTL.CIs <- ggplot() +
  geom_histogram(data = dxy_summary, aes(x = avg_dxy_per_cM, y = after_stat(density)), bins = 50, fill = "grey80", color = "black") +
  geom_ribbon(data = tail_left_dxy, aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.35) +
  geom_ribbon(data = tail_right_dxy, aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.35) +
  geom_density(data = dxy_summary, aes(x = avg_dxy_per_cM, y = after_stat(density)), color = "darkblue", linewidth = 0.8) +
  geom_vline(data = qtl_dxy_bp, aes(xintercept = avg_dxy_QTL, color = trait), linetype = "solid", linewidth = 1, show.legend = TRUE, key_glyph = "point") +
  scale_color_manual(values = trait_colors) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5, linetype = 0))) +
  labs(x = "average Dxy", y = "frequency", color = "Trait") +
  theme(legend.position = "right", legend.justification = "center")

# Fst: 1 cM bins, QTL window averages in cM (all_QTL_CI_fst)
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
  summarise(avg_fst_per_cM = mean(avg_wc_fst, na.rm = TRUE), n_windows = n(), .groups = "drop")
trait_labels <- c(
  SD = "Seed dormancy", GR = "Growth rate", tough = "Leaf toughness",
  Nreabs = "Nitrogen reabsorption\nduring drought", Chlreabs = "Chlorophyll reabsorption\nduring drought"
)
fst_cm <- fst %>%
  arrange(chrN, cM_window_start) %>%
  group_by(chrN) %>%
  mutate(cM_window_end = lead(cM_window_start)) %>%
  ungroup()
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
# Density and 5th/95th percentile tails for Fst (shade outer 5% each tail)
d_fst <- density(fst_summary$avg_fst_per_cM)
p05_fst <- quantile(fst_summary$avg_fst_per_cM, 0.05)
p95_fst <- quantile(fst_summary$avg_fst_per_cM, 0.95)
tail_left_fst <- data.frame(x = d_fst$x, y = d_fst$y) %>% filter(x <= p05_fst)
tail_right_fst <- data.frame(x = d_fst$x, y = d_fst$y) %>% filter(x >= p95_fst)

fst_windows_hist_w.QTL.CIs <- ggplot() +
  geom_histogram(data = fst_summary, aes(x = avg_fst_per_cM, y = after_stat(density)), binwidth = 0.01, fill = "grey80", color = "black") +
  geom_ribbon(data = tail_left_fst, aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.35) +
  geom_ribbon(data = tail_right_fst, aes(x = x, ymin = 0, ymax = y), fill = "red", alpha = 0.35) +
  geom_density(data = fst_summary, aes(x = avg_fst_per_cM, y = after_stat(density)), color = "darkblue", linewidth = 0.8) +
  geom_vline(data = qtl_fst_cm, aes(xintercept = avg_fst_QTL, color = trait), linetype = "solid", linewidth = 1, show.legend = TRUE, key_glyph = "point") +
  scale_color_manual(values = trait_colors, labels = trait_labels) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5, linetype = 0))) +
  labs(x = "Fst (1 cM windows)", y = "probability density", color = "QTL Trait:") +
  theme(legend.position = "left", legend.justification = "center") +
  scale_x_continuous(expand = c(-.01, -.01))

# Solo Fst distribution figure
ggsave(file.path(project_root, "Fst_genomewide_distribution_w_QTL_windows.png"), fst_windows_hist_w.QTL.CIs,
  device = "png", width = 15, height = 5.5, units = "cm")
# Three-panel: Fst, Dxy, f_dM stacked
three_panel_hist <- fst_windows_hist_w.QTL.CIs /
  (dxy_windows_hist_w.QTL.CIs + guides(color = "none")) /
  (f_dM_windows_hist_w.QTL.CIs + guides(color = "none")) +
  plot_layout(guides = "collect")
ggsave(file.path(project_root, "Fst_Dxy_fdM_genomewide_distribution_w_QTL_3panel.png"), three_panel_hist,
  device = "png", width = 17, height = 18, units = "cm")

# ---- Recombination rate vs Dsuite stats (f_dM, f_d, d_f) from cM–bp relationship ----
# Recombination rate = |ΔcM| / window_bp * 1e6 (cM/Mb) per Dsuite window; one dataset for all three stats
dsuite_with_recomb <- LASIout_w50s25_all9 %>%
  select(chr, windowStart, windowEnd, f_dM, f_d, d_f) %>%
  mutate(chrN = as.numeric(str_extract(chr, "\\d+"))) %>%
  filter(chrN != 4) %>%
  rowwise() %>%
  mutate(
    cM_start = bp_to_cm_funs[[as.character(chrN)]](windowStart),
    cM_end = bp_to_cm_funs[[as.character(chrN)]](windowEnd),
    window_bp = windowEnd - windowStart,
    recomb_cM_Mb = if_else(window_bp > 0, abs(cM_end - cM_start) / window_bp * 1e6, NA_real_)
  ) %>%
  ungroup() %>%
  filter(is.finite(recomb_cM_Mb))

# f_dM: positive values only (same as histogram logic)
recomb_vs_fdM_plot <- ggplot(dsuite_with_recomb %>% filter(f_dM > 0), aes(x = recomb_cM_Mb, y = f_dM)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, color = "darkred", linewidth = 0.8) +
  labs(x = "Recombination rate (cM/Mb)", y = "f_dM") +
  theme(panel.grid.minor = element_blank())

# f_d
recomb_vs_f_d_plot <- ggplot(dsuite_with_recomb %>% filter(is.finite(f_d)), aes(x = recomb_cM_Mb, y = f_d)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, color = "darkred", linewidth = 0.8) +
  labs(x = "Recombination rate (cM/Mb)", y = "f_d") +
  theme(panel.grid.minor = element_blank())

# d_f
recomb_vs_d_f_plot <- ggplot(dsuite_with_recomb %>% filter(is.finite(d_f)), aes(x = recomb_cM_Mb, y = d_f)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, color = "darkred", linewidth = 0.8) +
  labs(x = "Recombination rate (cM/Mb)", y = "d_f") +
  theme(panel.grid.minor = element_blank())

# Three-panel: recombination rate vs f_dM, f_d, d_f
recomb_by_stat_3panel <- recomb_vs_fdM_plot / recomb_vs_f_d_plot / recomb_vs_d_f_plot

if (!dir.exists(file.path(project_root, "Dsuite"))) dir.create(file.path(project_root, "Dsuite"))
ggsave(file.path(project_root, "Dsuite", "recombination_rate_vs_f_dM.png"),
       plot = recomb_vs_fdM_plot,
       device = "png", width = 14, height = 10, units = "cm", dpi = 300)
ggsave(file.path(project_root, "Dsuite", "recombination_rate_vs_f_d.png"),
       plot = recomb_vs_f_d_plot,
       device = "png", width = 14, height = 10, units = "cm", dpi = 300)
ggsave(file.path(project_root, "Dsuite", "recombination_rate_vs_d_f.png"),
       plot = recomb_vs_d_f_plot,
       device = "png", width = 14, height = 10, units = "cm", dpi = 300)
ggsave(file.path(project_root, "Dsuite", "recombination_rate_vs_Dsuite_stats_3panel.png"),
       plot = recomb_by_stat_3panel,
       device = "png", width = 14, height = 24, units = "cm", dpi = 300)
