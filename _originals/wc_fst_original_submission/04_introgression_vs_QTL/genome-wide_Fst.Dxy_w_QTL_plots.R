# =============================================================================
# genome-wide_Fst.Dxy_w_QTL_plots.R
# Genome-wide Fst, Dxy, and f_dM with QTL locations and distribution histograms
# =============================================================================
#
# Single script: (1) imports pixy + Dsuite f_dM once, (2) along-genome Fst/Dxy
# plots (500 kb window averages + QTL peaks/CIs), (3) 1 cM windowed distribution
# histograms (Fst, Dxy, f_dM) with QTL window averages, density curves, tail shading;
# (4) first-pass simplified 1 cM Fst density (black curve + dark-grey QTL CI means);
# (5) Fst at 10 kb (native pixy) with QTL lines at max Fst within each 95% CI.
#
# Dependencies: tidyverse, ggplot2, patchwork, ggrastr (rasterized points on Fst panels).
# QTL CSV, lepmap; ggsave outputs in "pixy/", "Dsuite/", and working directory.
# =============================================================================

# ---- load dependencies ----
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggrastr)

theme_set(theme_bw(base_size = 11))
fst_raster_dpi <- 300

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
for (p in c(
  file.path("scripts", "project_paths.R"),
  file.path("..", "scripts", "project_paths.R"),
  file.path("..", "..", "scripts", "project_paths.R")
)) {
  if (file.exists(p)) { source(p); break }
}
paths <- project_paths()
project_root <- paths$project_root
pixy_folder <- paths$pixy_folder
dsuite_dir <- paths$dsuite_dir
qtl_ci_path <- paths$qtl_ci_path
lepmap_dir <- paths$lepmap_dir
if (!file.exists(qtl_ci_path)) {
  qtl_ci_path <- file.path(project_root, "QTL", "all_QTL_CIs_bp_cm.251130.csv")
}
final_figures_dir <- file.path(
  Sys.getenv("HOME"),
  "Library/CloudStorage/GoogleDrive-jgharencar@ucdavis.edu",
  "My Drive/Costus/PopGen/Figures&Tables/Final_figures"
)

# ---- Pixy: import ----
pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)
pixy_wide <- pivot_wider(pixy_df, names_from = statistic, values_from = value) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000)
# make chromosome number only col:
pixy_wide <- pixy_wide %>%
  mutate(chrN = as.numeric(str_extract(chromosome, "\\d+")))

# ---- Genome-wide mean Fst (W&C from pixy) ----
mean_fst_windows <- mean(pixy_wide$avg_wc_fst, na.rm = TRUE)
mean_fst_snp_weighted <- weighted.mean(
  pixy_wide$avg_wc_fst,
  w = pixy_wide$no_snps,
  na.rm = TRUE
)
cat(
  "Genome-wide mean Fst (alle vs vill):\n",
  "  mean across 10 kb windows (unweighted) = ", sprintf("%.6f", mean_fst_windows), "\n",
  "  SNP-weighted mean (weights = no_snps per window) = ", sprintf("%.6f", mean_fst_snp_weighted), "\n",
  sep = ""
)

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

# x-axis: data limits + slight expand so first/last windows are not clipped at panel edges
cutoff99.9thPercentile <- quantile(abs(LASIout_w50s25_all9$f_dM[LASIout_w50s25_all9$f_dM < 0]), 0.999)

LASIout_positive <- LASIout_w50s25_all9 %>% filter(f_dM > 0)

fdM_xlim <- range(LASIout_positive$windowStart_cum, na.rm = TRUE)

LASI_99.9thper_w50s25_f_dM <- ggplot(data = LASIout_positive, aes(x = windowStart_cum, y = f_dM, color = as_factor(chr))) +
  geom_point(size = 3) +
  geom_hline(yintercept = cutoff99.9thPercentile, linetype = "dashed", color = "red") +
  scale_color_manual(values = rep(c("black", "darkgray"), ceiling(length(unique(axis_set_dsuite$chr)) / 2))) +
  guides(color = "none") +
  geom_point(
    data = subset(LASIout_positive, f_dM > cutoff99.9thPercentile),
    color = "red",
    size = 4.5
  ) +
  scale_x_continuous(
    label = axis_set_dsuite$chr_label,
    breaks = axis_set_dsuite$center,
    limits = fdM_xlim,
    expand = expansion(mult = 0.018, add = 0)
  ) +
  scale_y_continuous(expand = expansion(mult = 0.06, add = 0)) +
  labs(x = NULL, y = "f_dM") +
  theme(
    plot.margin = unit(c(2, 1, 1, 1), "cm"),
    axis.text = element_text(size = 20 ),
    axis.title = element_text(size = 22)
  )

if (!dir.exists(file.path(project_root, "Dsuite"))) dir.create(file.path(project_root, "Dsuite"))
ggsave(file.path(project_root, "Dsuite", "LASI_99.9thper_w50s25_f_dM_pos_only_all_chroms.png"),
       plot = LASI_99.9thper_w50s25_f_dM,
       width = 24, height = 3.8, units = "in", dpi = 300)

# ---- Genome-wide f_dM figure (all finite values; negatives included) ----
LASIout_fd <- LASIout_w50s25_all9 %>% filter(is.finite(f_dM))
fdM_xlim_all <- range(LASIout_fd$windowStart_cum, na.rm = TRUE)
cutoff99.9_abs_fdM <- quantile(abs(LASIout_fd$f_dM), 0.999, na.rm = TRUE)

LASI_f_dM_incl_negative <- ggplot(data = LASIout_fd, aes(x = windowStart_cum, y = f_dM, color = as_factor(chr))) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_point(size = 3) +
  scale_color_manual(values = rep(c("black", "darkgray"), ceiling(length(unique(axis_set_dsuite$chr)) / 2))) +
  guides(color = "none") +
  geom_point(
    data = LASIout_fd %>% filter(abs(f_dM) > cutoff99.9_abs_fdM),
    aes(x = windowStart_cum, y = f_dM),
    inherit.aes = FALSE,
    color = "red",
    size = 4.5
  ) +
  scale_x_continuous(
    label = axis_set_dsuite$chr_label,
    breaks = axis_set_dsuite$center,
    limits = fdM_xlim_all,
    expand = expansion(mult = 0.018, add = 0)
  ) +
  scale_y_continuous(expand = expansion(mult = 0.06, add = 0)) +
  labs(x = NULL, y = "f_dM") +
  theme(plot.margin = unit(c(2, 1, 1, 1), "cm"))

ggsave(file.path(project_root, "Dsuite", "LASI_99.9thper_w50s25_f_dM_incl_negative_all_chroms.png"),
       plot = LASI_f_dM_incl_negative,
       width = 24, height = 5.5, units = "in", dpi = 300)

# ---- QTL and trait colors (used by along-genome plots and histograms) ----
all_QTL_CI <- read.csv(qtl_ci_path, header = TRUE) %>%
  filter(!chr == 2 | !trait == "Nreabs") %>%
  filter(!chr == 5 | !trait == "GR")
trait_colors <- c(
  SD = "#0072B2",
  GR = "#E69F00",
  tough = "#8250C4",
  Nreabs = "#CC79A7",
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
cutoff99.9_dxy <- {
  x <- pixy_wide$avg_dxy
  x <- pmax(x, 0)
  quantile(x[is.finite(x)], 0.999)
}

plot_chr_combined <- function(raw_df, win_df, cutoff99.9) {
  ggplot(raw_df, aes(x = window_pos_1)) +
    geom_point(
      data = raw_df %>% filter(pmax(avg_dxy, 0) <= cutoff99.9),
      aes(y = pmax(avg_dxy, 0)),
      size = 1,
      alpha = 0.4,
      color = "grey40"
    ) +
    geom_point(
      data = raw_df %>% filter(pmax(avg_dxy, 0) > cutoff99.9),
      aes(y = pmax(avg_dxy, 0)),
      size = 2,
      alpha = 0.9,
      color = "red"
    ) +
    geom_line(
      data = win_df,
      aes(x = WindowMid, y = MeanY),
      color = "darkred",
      linewidth = 0.7
    ) +
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
    win_df = win_df_i,
    cutoff99.9 = cutoff99.9_dxy
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

# Fst supplement: merged fill+color legend (CI band + peak line), collected once below grid
add_qtl_layers_fst <- function(base_plot, qtl_rows, trait_colors) {
  if (nrow(qtl_rows) == 0) {
    return(base_plot + theme(legend.position = "none"))
  }

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
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = trait_colors, name = NULL, drop = FALSE) +
    scale_color_manual(values = trait_colors, name = NULL, drop = FALSE) +
    theme(legend.position = "none")
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

# 2-panel figure: chromosome 1 above chromosome 6 (Dxy)
Dxy_chr1_chr6_2panel <- (final_plots$Chrom1 / final_plots$Chrom6) &
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
ggsave(file.path(project_root, "pixy", "Dxy_chr1_chr6_2panel.png"), Dxy_chr1_chr6_2panel,
       device = "png", width = 22, height = 10, units = "cm")

Dxy_LASI_sep_chrs_w.QTL <- final_plots$Chrom1 + final_plots$Chrom2 + final_plots$Chrom3 +
  final_plots$Chrom5 + final_plots$Chrom6 + final_plots$Chrom7 + final_plots$Chrom8 + final_plots$Chrom9

ggsave(file.path(project_root, "pixy", "Dxy_10kb_LASIout_separate_chrs_w.500kb.avg_7QTL_narrowed.CIs.png"),
       Dxy_LASI_sep_chrs_w.QTL, device = "png", width = 50, height = 20, units = "cm")

# ---- Fst: plot with QTL ----
cutoff99.9_fst <- {
  x <- pixy_wide$avg_wc_fst
  x <- pmax(x, 0)
  quantile(x[is.finite(x)], 0.999)
}

plot_chr_combined <- function(raw_df, win_df, cutoff99.9) {
  ggplot(raw_df, aes(x = window_pos_1)) +
    geom_point_rast(
      data = raw_df %>% filter(pmax(avg_wc_fst, 0) <= cutoff99.9),
      aes(y = pmax(avg_wc_fst, 0)),
      size = 1,
      alpha = 0.4,
      color = "grey40",
      raster.dpi = fst_raster_dpi
    ) +
    geom_point_rast(
      data = raw_df %>% filter(pmax(avg_wc_fst, 0) > cutoff99.9),
      aes(y = pmax(avg_wc_fst, 0)),
      size = 2,
      alpha = 0.9,
      color = "red",
      raster.dpi = fst_raster_dpi
    ) +
    geom_line(
      data = win_df,
      aes(x = WindowMid, y = MeanY),
      color = "darkred",
      linewidth = 0.7
    ) +
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
    win_df = win_df_i,
    cutoff99.9 = cutoff99.9_fst
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

  plot_i <- add_qtl_layers_fst(base_plot, qtl_here, trait_colors)
  plot_i <- plot_i + ggtitle(paste("Chr", chr_num))
  final_plots[[i]] <- plot_i
}
names(final_plots) <- paste0("Chrom", chrom_list)

# 2-panel figure: chromosome 1 above chromosome 3 (Fst)
Fst_chr1_chr3_2panel <- (final_plots$Chrom1 / final_plots$Chrom3) &
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
ggsave(file.path(project_root, "pixy", "Fst_chr1_chr3_2panel.png"), Fst_chr1_chr3_2panel,
       device = "png", width = 22, height = 10, units = "cm")

# 2-panel figure: chromosome 1 above chromosome 6 (Fst)
Fst_chr1_chr6_2panel <- (final_plots$Chrom1 / final_plots$Chrom6) &
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
ggsave(file.path(project_root, "pixy", "Fst_chr1_chr6_2panel.png"), Fst_chr1_chr6_2panel,
       device = "png", width = 22, height = 10, units = "cm")

# 8 chromosomes in 2x4 grid; one shared legend (not collected per panel — avoids duplicate traits)
fst_sup_theme <- theme(
  plot.title = element_text(size = 10, hjust = 0.5),
  axis.title = element_text(size = 9),
  axis.text = element_text(size = 8),
  legend.position = "none"
)
for (nm in names(final_plots)) {
  final_plots[[nm]] <- final_plots[[nm]] + fst_sup_theme
}

qtl_legend_df <- tibble(
  trait = factor(names(trait_colors), levels = names(trait_colors)),
  lower95CI_bp = seq_along(trait_colors),
  upper95CI_bp = seq_along(trait_colors) + 1,
  peak_bp = seq_along(trait_colors) + 0.5
)
fst_qtl_legend_plot <- ggplot() +
  geom_rect(
    data = qtl_legend_df,
    aes(
      xmin = lower95CI_bp,
      xmax = upper95CI_bp,
      ymin = -Inf,
      ymax = Inf,
      fill = trait
    ),
    alpha = 0.2
  ) +
  geom_vline(
    data = qtl_legend_df,
    aes(xintercept = peak_bp, color = trait),
    linewidth = 1
  ) +
  scale_fill_manual(values = trait_colors, name = NULL, drop = FALSE) +
  scale_color_manual(values = trait_colors, name = NULL, drop = FALSE) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.title = element_blank()
  )
fst_legend_grob <- get_legend(fst_qtl_legend_plot)

fst_panels <- (
  final_plots$Chrom1 + final_plots$Chrom2 + final_plots$Chrom3 + final_plots$Chrom5 +
  final_plots$Chrom6 + final_plots$Chrom7 + final_plots$Chrom8 + final_plots$Chrom9
) +
  plot_layout(ncol = 2, nrow = 4) &
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.position = "none"
  )

Fst_LASI_sep_chrs_w.QTL <- fst_panels / wrap_elements(full = fst_legend_grob) +
  plot_layout(heights = c(1, 0.06), guides = "keep")

fst_fig_w_cm <- 22
fst_fig_h_cm <- 28
fst_pixy_pdf <- file.path(
  project_root, "pixy",
  "Fst_10kb_LASIout_separate_chrs_w.500kb.avg_7QTL_narrowed.CIs.pdf"
)
fst_pixy_png <- sub("\\.pdf$", ".png", fst_pixy_pdf)

ggsave(fst_pixy_pdf, Fst_LASI_sep_chrs_w.QTL, device = "pdf", width = fst_fig_w_cm, height = fst_fig_h_cm, units = "cm")
ggsave(fst_pixy_png, Fst_LASI_sep_chrs_w.QTL, device = "png", width = fst_fig_w_cm, height = fst_fig_h_cm, units = "cm", dpi = 200)

if (dir.exists(final_figures_dir)) {
  ggsave(
    file.path(final_figures_dir, "sup_genomewide_Fst_QTL_CIs.pdf"),
    Fst_LASI_sep_chrs_w.QTL,
    device = "pdf",
    width = fst_fig_w_cm,
    height = fst_fig_h_cm,
    units = "cm"
  )
  message("Wrote sup_genomewide_Fst_QTL_CIs to ", final_figures_dir)
}

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
# 5th and 95th percentile cutoffs (red dashed vertical lines)
p05_fdM <- quantile(f_dM_summary$avg_f_dM_per_cM, 0.05)
p95_fdM <- quantile(f_dM_summary$avg_f_dM_per_cM, 0.95)

f_dM_windows_hist_w.QTL.CIs <- ggplot() +
  geom_histogram(data = f_dM_summary, aes(x = avg_f_dM_per_cM, y = after_stat(density)), bins = 50, fill = "grey80", color = "black") +
  geom_vline(xintercept = p05_fdM, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = p95_fdM, linetype = "dashed", color = "red", linewidth = 0.8) +
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
# 5th and 95th percentile cutoffs (red dashed vertical lines)
p05_dxy <- quantile(dxy_summary$avg_dxy_per_cM, 0.05)
p95_dxy <- quantile(dxy_summary$avg_dxy_per_cM, 0.95)

dxy_windows_hist_w.QTL.CIs <- ggplot() +
  geom_histogram(data = dxy_summary, aes(x = avg_dxy_per_cM, y = after_stat(density)), bins = 50, fill = "grey80", color = "black") +
  geom_vline(xintercept = p05_dxy, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = p95_dxy, linetype = "dashed", color = "red", linewidth = 0.8) +
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
# Cache 1 cM Fst objects so the simplified plot chunk can be re-run alone
saveRDS(
  list(fst_summary = fst_summary, qtl_fst_cm = qtl_fst_cm, project_root = project_root),
  file.path(project_root, "Fst_1cM_window_cache.rds")
)
# 5th and 95th percentile cutoffs (red dashed vertical lines)
p05_fst <- quantile(fst_summary$avg_fst_per_cM, 0.05)
p95_fst <- quantile(fst_summary$avg_fst_per_cM, 0.95)

fst_windows_hist_w.QTL.CIs <- ggplot() +
  geom_histogram(data = fst_summary, aes(x = avg_fst_per_cM, y = after_stat(density)), binwidth = 0.01, fill = "grey80", color = "black") +
  geom_vline(xintercept = p05_fst, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = p95_fst, linetype = "dashed", color = "red", linewidth = 0.8) +
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

# First-pass simplified 1 cM Fst: density only, dark-grey QTL CI means, no histogram/percentiles
# This block is re-runnable alone if Fst_1cM_window_cache.rds exists (written above on a full run).
if (!exists("fst_summary") || !exists("qtl_fst_cm")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  library(ggplot2)
  cache_candidates <- c(
    if (exists("project_root")) file.path(project_root, "Fst_1cM_window_cache.rds") else NULL,
    "Fst_1cM_window_cache.rds",
    file.path("..", "Fst_1cM_window_cache.rds"),
    file.path("..", "..", "Fst_1cM_window_cache.rds")
  )
  cache_path <- cache_candidates[file.exists(cache_candidates)][1]
  if (is.na(cache_path) || !nzchar(cache_path)) {
    stop(
      "fst_summary/qtl_fst_cm not in memory and Fst_1cM_window_cache.rds not found. ",
      "Run the script once through the 1 cM Fst summary (needs pixy + lepmap + QTL CIs)."
    )
  }
  cache <- readRDS(cache_path)
  fst_summary <- cache$fst_summary
  qtl_fst_cm <- cache$qtl_fst_cm
  if (!exists("project_root")) project_root <- cache$project_root
}
fst_windows_density_simplified <- ggplot() +
  geom_density(
    data = fst_summary,
    aes(x = avg_fst_per_cM, y = after_stat(density)),
    color = "black",
    linewidth = 0.8
  ) +
  geom_vline(
    data = qtl_fst_cm,
    aes(xintercept = avg_fst_QTL),
    color = "grey60",
    linetype = "solid",
    linewidth = 1
  ) +
  labs(x = "Distribution of Fst (1cM windows)", y = "probability density") +
  scale_x_continuous(expand = c(-.01, -.01)) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA)
  )
ggsave(
  file.path(project_root, "Fst_genomewide_distribution_w_QTL_windows_simplified.png"),
  fst_windows_density_simplified,
  device = "png",
  width = 7.5,
  height = 4.5,
  units = "cm"
)

# ---- Fst: 10 kb window histogram; QTL lines = max Fst window in CI ----
# Blue curve = geom_density() (kernel-smoothed genome-wide 10 kb Fst PDF).
# Spike at Fst ~ 0: many windows with Fst = 0 or << 0.01 (low/no allelic differentiation,
# often few SNPs); binwidth 0.01 stacks those in the leftmost histogram bin.
max_fst_in_qtl_bp <- function(qtl_rows, fst_df) {
  qtl_rows %>%
    rowwise() %>%
    mutate(
      max_fst_QTL = max(
        fst_df$avg_wc_fst[
          fst_df$chrN == chr &
            fst_df$window_pos_1 <= upper95CI_bp &
            fst_df$window_pos_2 >= lower95CI_bp
        ],
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

mean_fst_in_qtl_bp <- function(qtl_rows, fst_df) {
  qtl_rows %>%
    rowwise() %>%
    mutate(
      avg_fst_QTL = mean(
        fst_df$avg_wc_fst[
          fst_df$chrN == chr &
            fst_df$window_pos_1 <= upper95CI_bp &
            fst_df$window_pos_2 >= lower95CI_bp
        ],
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

make_fst_bp_hist <- function(summary_df, x_col, qtl_df, xlab, outfile) {
  p05 <- quantile(summary_df[[x_col]], 0.05, na.rm = TRUE)
  p95 <- quantile(summary_df[[x_col]], 0.95, na.rm = TRUE)
  p <- ggplot() +
    geom_histogram(
      data = summary_df,
      aes(x = .data[[x_col]], y = after_stat(density)),
      binwidth = 0.01,
      fill = "grey80",
      color = "black"
    ) +
    geom_vline(xintercept = p05, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_vline(xintercept = p95, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_density(
      data = summary_df,
      aes(x = .data[[x_col]], y = after_stat(density)),
      color = "darkblue",
      linewidth = 0.8
    ) +
    geom_vline(
      data = qtl_df,
      aes(xintercept = max_fst_QTL, color = trait),
      linetype = "solid",
      linewidth = 1,
      show.legend = TRUE,
      key_glyph = "point"
    ) +
    scale_color_manual(values = trait_colors, labels = trait_labels) +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 5, linetype = 0))) +
    labs(x = xlab, y = "probability density", color = "QTL Trait:") +
    theme(legend.position = "left", legend.justification = "center") +
    scale_x_continuous(expand = c(-0.01, -0.01))
  ggsave(outfile, p, device = "png", width = 15, height = 5.5, units = "cm")
  invisible(p)
}

# 10 kb: native pixy windows (histogram = per-window Fst; QTL = max window in CI)
fst_10kb_summary <- fst %>% select(avg_fst = avg_wc_fst)
qtl_fst_10kb_max <- max_fst_in_qtl_bp(all_QTL_CI_fst, fst)
fst_10kb_hist <- make_fst_bp_hist(
  fst_10kb_summary,
  "avg_fst",
  qtl_fst_10kb_max,
  "Fst (10 kb windows)",
  file.path(project_root, "Fst_genomewide_distribution_10kb_w_QTL_max_windows.png")
)

# 10 kb: mean Fst across windows in each QTL CI (no genome-wide histogram)
qtl_fst_10kb_avg <- mean_fst_in_qtl_bp(all_QTL_CI_fst, fst) %>%
  mutate(qtl_label = paste0(trait, " (chr", chr, ")"))
fst_10kb_qtl_mean_plot <- ggplot(qtl_fst_10kb_avg, aes(x = avg_fst_QTL, y = reorder(qtl_label, avg_fst_QTL), color = trait)) +
  geom_segment(aes(x = 0, xend = avg_fst_QTL, y = qtl_label, yend = qtl_label), linewidth = 0.9, color = "grey70") +
  geom_point(size = 4) +
  scale_color_manual(values = trait_colors, labels = trait_labels) +
  labs(
    x = "Mean Fst (10 kb windows in QTL 95% CI)",
    y = NULL,
    color = "QTL Trait:"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "left", legend.justification = "center")
ggsave(
  file.path(project_root, "Fst_10kb_QTL_CI_mean_only.png"),
  fst_10kb_qtl_mean_plot,
  device = "png",
  width = 12,
  height = 6,
  units = "cm"
)

# Three-panel: Fst, Dxy, f_dM stacked
three_panel_hist <- fst_windows_hist_w.QTL.CIs /
  (dxy_windows_hist_w.QTL.CIs + guides(color = "none")) /
  (f_dM_windows_hist_w.QTL.CIs + guides(color = "none")) +
  plot_layout(guides = "collect")
ggsave(file.path(project_root, "Fst_Dxy_fdM_genomewide_distribution_w_QTL_3panel.png"), three_panel_hist,
  device = "png", width = 17, height = 18, units = "cm")
