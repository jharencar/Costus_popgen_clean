# =============================================================================
# combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R
#
# Same stacked layout as Population_genetics_vill_alle.qmd combo_Dsuite_pixy_5panel
# (f_dM / d_f / Dxy / pi / Fst), but with four panels only: f_dM, d_f, Dxy, Fst (no pi).
# f_dM is plotted on the full signed scale (positive and negative); cutoffs are symmetric
# from the 99th percentile of |f_dM| among negative windows (qmd convention).
# QTL: 95% CI ribbons + peak vlines from QTL/all_QTL_CIs_bp_cm.251130.csv (same filters as
# genome-wide_*_w_QTL_plots.R); x positions use cumulative bp (Dsuite offsets for f_dM/d_f,
# pixy offsets for Dxy/Fst). CI layers are drawn last so they appear above the points.
# Point layers use ggrastr (geom_point_rast) to keep PDF size manageable.
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggrastr)
library(grid)

# ---- Paths ----
for (p in c(
  file.path("scripts", "project_paths.R"),
  file.path("..", "scripts", "project_paths.R"),
  file.path("..", "..", "scripts", "project_paths.R")
)) {
  if (file.exists(p)) { source(p); break }
}
paths <- project_paths()
project_root <- paths$project_root
dsuite_dir <- paths$dsuite_dir
pixy_folder <- paths$pixy_folder
qtl_ci_path <- paths$qtl_ci_path
if (!file.exists(qtl_ci_path)) {
  qtl_ci_path <- file.path(project_root, "QTL", "all_QTL_CIs_bp_cm.251130.csv")
}
final_figures_dir <- file.path(
  Sys.getenv("HOME"),
  "Library/CloudStorage/GoogleDrive-jgharencar@ucdavis.edu",
  "My Drive/Costus/PopGen/Figures&Tables/Final_figures"
)
if (nzchar(Sys.getenv("PLOT_OUT_ROOT", ""))) {
  project_root <- Sys.getenv("PLOT_OUT_ROOT")
  final_figures_dir <- project_root
  dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
}
point_size <- 3
combo_raster_dpi <- 300
legend_text_size <- 40

# Same palette as genome-wide_Fst.Dxy_w_QTL_plots.R (SD, GR, tough, Nreabs, Chlreabs)
trait_colors <- c(
  SD = "#0072B2",
  GR = "#E69F00",
  tough = "#8250C4",
  Nreabs = "#CC79A7",
  Chlreabs = "#009E73"
)

# Same row filters as genome-wide_Fst.Dxy_w_QTL_plots.R / genome-wide_*_w_QTL_plots.R
# Opposite of expected species direction: dashed peaks. GR chromosome 5 is listed
# and is already removed by the CI filter.
opposite_direction_qtl <- tibble::tribble(
  ~trait, ~chr,
  "SD", 6L,
  "GR", 5L,
  "tough", 3L,
  "tough", 8L
)
all_QTL_CI <- read.csv(qtl_ci_path, header = TRUE) %>%
  filter(!chr == 2 | !trait == "Nreabs") %>%
  filter(!chr == 5 | !trait == "GR") %>%
  mutate(
    opposite_direction = paste(trait, as.integer(chr)) %in%
      paste(opposite_direction_qtl$trait, opposite_direction_qtl$chr),
    qtl_linetype = ifelse(opposite_direction, "dashed", "solid")
  )

# =============================================================================
# Dsuite: f_dM and d_f (from Population_genetics_vill_alle.qmd "Plotting f_dM and d_f")
# =============================================================================
LASIout_w50s25_all9 <- data.frame()
file_list <- list.files(path = dsuite_dir, pattern = "*.txt", full.names = TRUE)
for (file in file_list) {
  temp_df <- read.delim(file, sep = "\t")
  LASIout_w50s25_all9 <- rbind(LASIout_w50s25_all9, temp_df)
}

cutoff99thPercentile <- quantile(
  abs(LASIout_w50s25_all9$f_dM[LASIout_w50s25_all9$f_dM < 0]),
  0.99
)

data_cum <- LASIout_w50s25_all9 %>%
  group_by(chr) %>%
  summarise(max_windowStart = max(windowStart)) %>%
  mutate(windowStart_add = lag(cumsum(max_windowStart), default = 0)) %>%
  select(chr, windowStart_add)

LASIout_w50s25_all9 <- LASIout_w50s25_all9 %>%
  inner_join(data_cum, by = "chr") %>%
  mutate(
    windowStart_cum = windowStart + windowStart_add,
    chr_label = gsub("Chrom", "", chr)
  )

axis_set_dsuite <- LASIout_w50s25_all9 %>%
  group_by(chr) %>%
  summarize(center = mean(windowStart_cum), chr_label = first(chr_label))

dsuite_chr_add <- LASIout_w50s25_all9 %>%
  distinct(chr, windowStart_add) %>%
  mutate(chrN = as.numeric(gsub("Chrom", "", chr)))

qtl_dsuite <- all_QTL_CI %>%
  mutate(chrN = as.numeric(as.character(chr))) %>%
  left_join(dsuite_chr_add, by = "chrN") %>%
  mutate(
    lower95CI_cum = lower95CI_bp + windowStart_add,
    upper95CI_cum = upper95CI_bp + windowStart_add,
    peak_cum = peak_bp + windowStart_add
  )

dsuite_xlim <- range(LASIout_w50s25_all9$windowStart_cum, na.rm = TRUE)

# =============================================================================
# Pixy: Dxy and Fst (from qmd "My plotting to match Dsuite plots"; pi omitted)
# =============================================================================
pixy_to_long <- function(pixy_files) {
  pixy_df <- list()
  for (i in seq_along(pixy_files)) {
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    if (stat_file_type == "pi") {
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
          key = "statistic", value = "value"
        ) %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      pixy_df[[i]] <- df
    } else {
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
          key = "statistic", value = "value"
        )
      pixy_df[[i]] <- df
    }
  }
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

if (is.na(pixy_folder) || !nzchar(pixy_folder) || !dir.exists(pixy_folder)) {
  stop("10 kb pixy folder not found: ", pixy_folder, "\nSet PIXY_FOLDER to the Hudson run.")
}
pixy_files <- list.files(pixy_folder, pattern = "_(pi|dxy|fst)\\.txt$", full.names = TRUE)
if (length(pixy_files) == 0) {
  stop("No pixy *_pi.txt / *_dxy.txt / *_fst.txt files in ", pixy_folder)
}
pixy_df <- pixy_to_long(pixy_files)
pixy_wide <- pivot_wider(pixy_df, names_from = statistic, values_from = value) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000)
if (!"avg_hudson_fst" %in% names(pixy_wide)) {
  if ("avg_wc_fst" %in% names(pixy_wide)) {
    stop("Pixy folder is Weir–Cockerham (avg_wc_fst). Use the Hudson run, not output_251124_LASIref.")
  }
  stop("Pixy folder is missing avg_hudson_fst: ", pixy_folder)
}
pixy_wide <- pixy_wide %>%
  mutate(across(
    any_of(c("avg_hudson_fst", "avg_dxy", "avg_pi", "no_snps", "no_sites")),
    ~ suppressWarnings(as.numeric(.x))
  ))

pixy_data_cum <- pixy_wide %>%
  group_by(chromosome) %>%
  summarise(max_windowStart = max(window_pos_1)) %>%
  mutate(windowStart_add = lag(cumsum(max_windowStart), default = 0)) %>%
  select(chromosome, windowStart_add)

pixy_wide <- pixy_wide %>%
  inner_join(pixy_data_cum, by = "chromosome") %>%
  mutate(windowStart_cum = window_pos_1 + windowStart_add)

axis_set_pixy <- pixy_wide %>%
  group_by(chromosome) %>%
  summarize(center = mean(windowStart_cum))

pixy_chr_add <- pixy_wide %>%
  distinct(chromosome, windowStart_add) %>%
  mutate(chrN = as.numeric(gsub("Chrom", "", chromosome)))

qtl_pixy <- all_QTL_CI %>%
  mutate(chrN = as.numeric(as.character(chr))) %>%
  left_join(pixy_chr_add, by = "chrN") %>%
  mutate(
    lower95CI_cum = lower95CI_bp + windowStart_add,
    upper95CI_cum = upper95CI_bp + windowStart_add,
    peak_cum = peak_bp + windowStart_add
  )

pixy_xlim <- range(pixy_wide$windowStart_cum, na.rm = TRUE)
combo_xlim <- c(
  min(dsuite_xlim[1], pixy_xlim[1]),
  max(dsuite_xlim[2], pixy_xlim[2])
)
# Data limits + small symmetric buffer (tighter than genome-wide 0.018)
x_expand <- expansion(mult = 0.012, add = 0)

LASI_99thper_w50s25_f_dM <- ggplot(
  data = LASIout_w50s25_all9,
  aes(x = windowStart_cum, y = f_dM, color = as_factor(chr))
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_point_rast(size = point_size, raster.dpi = combo_raster_dpi) +
  geom_hline(
    yintercept = c(cutoff99thPercentile, -cutoff99thPercentile),
    linetype = "dashed",
    color = "red"
  ) +
  scale_color_manual(values = rep(c("black", "darkgray"), unique(length(axis_set_dsuite$chr)))) +
  guides(color = "none") +
  geom_point_rast(
    data = subset(
      LASIout_w50s25_all9,
      f_dM > cutoff99thPercentile | f_dM < -cutoff99thPercentile
    ),
    color = "red",
    size = point_size,
    raster.dpi = combo_raster_dpi
  ) +
  scale_x_continuous(
    label = axis_set_dsuite$chr_label,
    breaks = axis_set_dsuite$center,
    limits = combo_xlim,
    expand = x_expand
  ) +
  labs(x = NULL, y = "f_dM") +
  geom_rect(
    data = qtl_dsuite,
    aes(xmin = lower95CI_cum, xmax = upper95CI_cum, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.2,
    fill = trait_colors[as.character(qtl_dsuite$trait)]
  ) +
  geom_vline(
    data = qtl_dsuite,
    aes(xintercept = peak_cum, linetype = qtl_linetype),
    color = trait_colors[as.character(qtl_dsuite$trait)],
    linewidth = 0.7
  ) +
  scale_linetype_identity(guide = "none") +
  theme(
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 30),
    plot.margin = unit(c(2, 1, 1, 1), "cm")
  )

cutoff99thPercentile_d_f <- quantile(
  abs(LASIout_w50s25_all9$d_f[which(LASIout_w50s25_all9$d_f < 0)]),
  0.99
)

LASI_99thper_w50s25_d_f <- ggplot(
  data = LASIout_w50s25_all9,
  aes(x = windowStart_cum, y = d_f, color = as_factor(chr))
) +
  geom_point_rast(size = point_size, raster.dpi = combo_raster_dpi) +
  geom_hline(
    yintercept = c(cutoff99thPercentile_d_f, -cutoff99thPercentile_d_f),
    linetype = "dashed",
    color = "red"
  ) +
  scale_color_manual(values = rep(c("black", "darkgray"), unique(length(axis_set_dsuite$chr)))) +
  guides(color = "none") +
  geom_point_rast(
    data = subset(
      LASIout_w50s25_all9,
      d_f > cutoff99thPercentile_d_f | d_f < -cutoff99thPercentile_d_f
    ),
    color = "red",
    size = point_size,
    raster.dpi = combo_raster_dpi
  ) +
  scale_x_continuous(
    label = axis_set_dsuite$chr_label,
    breaks = axis_set_dsuite$center,
    limits = combo_xlim,
    expand = x_expand
  ) +
  labs(x = NULL, y = "d_f") +
  geom_rect(
    data = qtl_dsuite,
    aes(xmin = lower95CI_cum, xmax = upper95CI_cum, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.2,
    fill = trait_colors[as.character(qtl_dsuite$trait)]
  ) +
  geom_vline(
    data = qtl_dsuite,
    aes(xintercept = peak_cum, linetype = qtl_linetype),
    color = trait_colors[as.character(qtl_dsuite$trait)],
    linewidth = 0.7
  ) +
  scale_linetype_identity(guide = "none") +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40))

# Same 10 kb floors as genome-wide_Fst.Dxy_w_QTL_plots.R. Axis limits stay on the
# unfiltered chromosome span so this panel still aligns with the Dsuite panels.
pixy_dxy <- pixy_wide %>%
  filter(is.finite(avg_dxy), is.finite(no_sites), no_sites >= 1000)
pixy_fst <- pixy_wide %>%
  filter(is.finite(avg_hudson_fst), is.finite(no_snps), no_snps >= 10)
cutoff99thPercentile.dxy <- quantile(pixy_dxy$avg_dxy, 0.99, na.rm = TRUE)
Dxy_plot <- ggplot(data = pixy_dxy, aes(x = windowStart_cum, y = avg_dxy, color = as_factor(chromosome))) +
  geom_point_rast(size = point_size, raster.dpi = combo_raster_dpi) +
  scale_color_manual(values = rep(c("black", "darkgray"), unique(length(axis_set_pixy$chromosome)))) +
  guides(color = "none") +
  scale_x_continuous(
    label = axis_set_pixy$chromosome,
    breaks = axis_set_pixy$center,
    limits = combo_xlim,
    expand = x_expand
  ) +
  labs(x = NULL, y = "Dxy") +
  geom_hline(yintercept = cutoff99thPercentile.dxy, linetype = "dashed", color = "red") +
  geom_point_rast(
    data = subset(pixy_dxy, avg_dxy > cutoff99thPercentile.dxy),
    color = "red",
    size = point_size,
    raster.dpi = combo_raster_dpi
  ) +
  geom_rect(
    data = qtl_pixy,
    aes(xmin = lower95CI_cum, xmax = upper95CI_cum, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.2,
    fill = trait_colors[as.character(qtl_pixy$trait)]
  ) +
  geom_vline(
    data = qtl_pixy,
    aes(xintercept = peak_cum, linetype = qtl_linetype),
    color = trait_colors[as.character(qtl_pixy$trait)],
    linewidth = 0.7
  ) +
  scale_linetype_identity(guide = "none") +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40))

upper_99thPercentile.fst <- quantile(pixy_fst$avg_hudson_fst, 0.99, na.rm = TRUE)
lower_99thPercentile.fst <- quantile(pixy_fst$avg_hudson_fst, 0.01, na.rm = TRUE)
hudson_fst_plot <- ggplot(data = pixy_fst, aes(x = windowStart_cum, y = avg_hudson_fst, color = as_factor(chromosome))) +
  geom_point_rast(size = point_size, raster.dpi = combo_raster_dpi) +
  scale_color_manual(values = rep(c("black", "darkgray"), unique(length(axis_set_pixy$chromosome)))) +
  guides(color = "none") +
  scale_x_continuous(
    label = axis_set_pixy$chromosome,
    breaks = axis_set_pixy$center,
    limits = combo_xlim,
    expand = x_expand
  ) +
  labs(x = NULL, y = "Fst") +
  geom_hline(
    yintercept = c(upper_99thPercentile.fst, lower_99thPercentile.fst),
    linetype = "dashed",
    color = "red"
  ) +
  geom_point_rast(
    data = subset(
      pixy_fst,
      avg_hudson_fst > upper_99thPercentile.fst | avg_hudson_fst < lower_99thPercentile.fst
    ),
    color = "red",
    size = point_size,
    raster.dpi = combo_raster_dpi
  ) +
  geom_rect(
    data = qtl_pixy,
    aes(xmin = lower95CI_cum, xmax = upper95CI_cum, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.2,
    fill = trait_colors[as.character(qtl_pixy$trait)]
  ) +
  geom_vline(
    data = qtl_pixy,
    aes(xintercept = peak_cum, linetype = qtl_linetype),
    color = trait_colors[as.character(qtl_pixy$trait)],
    linewidth = 0.7
  ) +
  scale_linetype_identity(guide = "none") +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40))

# =============================================================================
# Four-panel combo (5-panel from qmd minus pi)
# =============================================================================
combo_panel_theme <- theme(
  plot.background = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
  axis.line = element_blank(),
  panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
  panel.grid.minor = element_blank(),
  plot.tag.position = c(0, 1),
  plot.tag = element_text(size = 40, hjust = 0, vjust = 0)
)

# Shared QTL legend (CI band + peak line); one row per trait — matches sup_genomewide_Fst_QTL_CIs
qtl_legend_df <- tibble(
  trait = factor(names(trait_colors), levels = names(trait_colors)),
  lower95CI_bp = seq_along(trait_colors),
  upper95CI_bp = seq_along(trait_colors) + 1,
  peak_bp = seq_along(trait_colors) + 0.5
)
qtl_legend_plot <- ggplot() +
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
    legend.text = element_text(size = legend_text_size),
    legend.key.height = unit(1.1, "cm"),
    legend.key.width = unit(1.6, "cm"),
    legend.title = element_blank()
  )
qtl_legend_grob <- get_legend(qtl_legend_plot)

# Tag only the four data panels (A–D). Explicit tags so the QTL legend row is not "E)".
combo_Dsuite_pixy_4panel_fdm_df_dxy_fst <-
  LASI_99thper_w50s25_f_dM +
  LASI_99thper_w50s25_d_f +
  Dxy_plot +
  hudson_fst_plot +
  wrap_elements(full = qtl_legend_grob) +
  plot_layout(
    design = "
      AAAA
      BBBB
      CCCC
      DDDD
      EEEE
    ",
    heights = c(0.55, 0.55, 1.45, 1.45, 0.08),
    guides = "keep"
  ) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D")), tag_suffix = ")") &
  combo_panel_theme

out_png <- file.path(project_root, "combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.png")
out_pdf <- file.path(project_root, "combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.pdf")
ggsave(
  out_png,
  combo_Dsuite_pixy_4panel_fdm_df_dxy_fst,
  device = "png",
  width = 120,
  height = 80,
  units = "cm",
  bg = "white"
)
ggsave(
  out_pdf,
  combo_Dsuite_pixy_4panel_fdm_df_dxy_fst,
  device = "pdf",
  width = 100,
  height = 70,
  units = "cm",
  bg = "white"
)

message("Wrote ", out_png)
message("Wrote ", out_pdf)

if (dir.exists(final_figures_dir)) {
  sup_pdf <- file.path(final_figures_dir, "sup_4panel_fdm_df_dxy_fst.pdf")
  ggsave(
    sup_pdf,
    combo_Dsuite_pixy_4panel_fdm_df_dxy_fst,
    device = "pdf",
    width = 100,
    height = 70,
    units = "cm",
    bg = "white"
  )
  message("Wrote ", sup_pdf)
}
