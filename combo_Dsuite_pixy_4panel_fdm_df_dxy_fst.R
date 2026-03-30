# =============================================================================
# combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R
#
# Same stacked layout as Population_genetics_vill_alle.qmd combo_Dsuite_pixy_5panel
# (f_dM / d_f / Dxy / pi / Fst), but with four panels only: f_dM, d_f, Dxy, Fst (no pi).
# f_dM is plotted on the full signed scale (positive and negative); cutoffs are symmetric
# from the 99.9th percentile of |f_dM| among negative windows (qmd convention).
# QTL: 95% CI ribbons + peak vlines from QTL/all_QTL_CIs_bp_cm.251130.csv (same filters as
# genome-wide_*_w_QTL_plots.R); x positions use cumulative bp (Dsuite offsets for f_dM/d_f,
# pixy offsets for Dxy/Fst). CI layers are drawn last so they appear above the points.
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(grid)

# ---- Paths ----
project_root <- "/Users/admin-jgharenc/Documents/Github/Costus_popgen_clean"
dsuite_dir <- file.path(project_root, "Dsuite", "Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120")
pixy_folder <- file.path(project_root, "pixy", "output_251124_LASIref")
qtl_ci_path <- file.path(project_root, "QTL", "all_QTL_CIs_bp_cm.251130.csv")
point_size <- 3

trait_colors <- c(
  SD = "#0072B2",
  GR = "#E69F00",
  tough = "#8250C4",
  Nreabs = "#CC79A7",
  Chlreabs = "#009E73"
)

# Same row filters as genome-wide_Fst.Dxy_w_QTL_plots.R / genome-wide_*_w_QTL_plots.R
all_QTL_CI <- read.csv(qtl_ci_path, header = TRUE) %>%
  filter(!chr == 2 | !trait == "Nreabs") %>%
  filter(!chr == 5 | !trait == "GR")

# =============================================================================
# Dsuite: f_dM and d_f (from Population_genetics_vill_alle.qmd "Plotting f_dM and d_f")
# =============================================================================
LASIout_w50s25_all9 <- data.frame()
file_list <- list.files(path = dsuite_dir, pattern = "*.txt", full.names = TRUE)
for (file in file_list) {
  temp_df <- read.delim(file, sep = "\t")
  LASIout_w50s25_all9 <- rbind(LASIout_w50s25_all9, temp_df)
}

cutoff99.9thPercentile <- quantile(
  abs(LASIout_w50s25_all9$f_dM[LASIout_w50s25_all9$f_dM < 0]),
  0.999
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

LASI_99.9thper_w50s25_f_dM <- ggplot(
  data = LASIout_w50s25_all9,
  aes(x = windowStart_cum, y = f_dM, color = as_factor(chr))
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_point(size = point_size) +
  coord_cartesian(xlim = c(24398, 1120430367)) +
  geom_hline(
    yintercept = c(cutoff99.9thPercentile, -cutoff99.9thPercentile),
    linetype = "dashed",
    color = "red"
  ) +
  scale_color_manual(values = rep(c("black", "darkgray"), unique(length(axis_set_dsuite$chr)))) +
  guides(color = "none") +
  geom_point(
    data = subset(
      LASIout_w50s25_all9,
      f_dM > cutoff99.9thPercentile | f_dM < -cutoff99.9thPercentile
    ),
    color = "red",
    size = point_size
  ) +
  scale_x_continuous(label = axis_set_dsuite$chr_label, breaks = axis_set_dsuite$center) +
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
    aes(xintercept = peak_cum),
    color = trait_colors[as.character(qtl_dsuite$trait)],
    linewidth = 0.7
  ) +
  theme(
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 30),
    plot.margin = unit(c(2, 1, 1, 1), "cm")
  )

cutoff99.9thPercentile_d_f <- quantile(
  abs(LASIout_w50s25_all9$d_f[which(LASIout_w50s25_all9$d_f < 0)]),
  0.999
)

LASI_99.9thper_w50s25_d_f <- ggplot(
  data = LASIout_w50s25_all9,
  aes(x = windowStart_cum, y = d_f, color = as_factor(chr))
) +
  geom_point(size = point_size) +
  coord_cartesian(xlim = c(24398, 1120430367)) +
  geom_hline(
    yintercept = c(cutoff99.9thPercentile_d_f, -cutoff99.9thPercentile_d_f),
    linetype = "dashed",
    color = "red"
  ) +
  scale_color_manual(values = rep(c("black", "darkgray"), unique(length(axis_set_dsuite$chr)))) +
  guides(color = "none") +
  geom_point(
    data = subset(
      LASIout_w50s25_all9,
      d_f > cutoff99.9thPercentile_d_f | d_f < -cutoff99.9thPercentile_d_f
    ),
    color = "red",
    size = point_size
  ) +
  scale_x_continuous(label = axis_set_dsuite$chr, breaks = axis_set_dsuite$center) +
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
    aes(xintercept = peak_cum),
    color = trait_colors[as.character(qtl_dsuite$trait)],
    linewidth = 0.7
  ) +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40))

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

pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)
pixy_wide <- pivot_wider(pixy_df, names_from = statistic, values_from = value) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000)

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

cutoff99.9thPercentile.dxy <- quantile(pixy_wide$avg_dxy, 0.999, na.rm = TRUE)
Dxy_plot <- ggplot(data = pixy_wide, aes(x = windowStart_cum, y = avg_dxy, color = as_factor(chromosome))) +
  geom_point(size = point_size) +
  scale_color_manual(values = rep(c("black", "darkgray"), unique(length(axis_set_pixy$chromosome)))) +
  guides(color = "none") +
  scale_x_continuous(label = axis_set_pixy$chromosome, breaks = axis_set_pixy$center) +
  labs(x = NULL, y = "Dxy") +
  geom_hline(yintercept = cutoff99.9thPercentile.dxy, linetype = "dashed", color = "red") +
  geom_point(
    data = subset(pixy_wide, avg_dxy > cutoff99.9thPercentile.dxy),
    color = "red",
    size = point_size
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
    aes(xintercept = peak_cum),
    color = trait_colors[as.character(qtl_pixy$trait)],
    linewidth = 0.7
  ) +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40))

upper_99.9thPercentile.fst <- quantile(pixy_wide$avg_wc_fst, 0.999, na.rm = TRUE)
lower_99.9thPercentile.fst <- quantile(pixy_wide$avg_wc_fst, 0.001, na.rm = TRUE)
wc_fst_plot <- ggplot(data = pixy_wide, aes(x = windowStart_cum, y = avg_wc_fst, color = as_factor(chromosome))) +
  geom_point(size = point_size) +
  scale_color_manual(values = rep(c("black", "darkgray"), unique(length(axis_set_pixy$chromosome)))) +
  guides(color = "none") +
  scale_x_continuous(label = axis_set_pixy$chromosome, breaks = axis_set_pixy$center) +
  labs(x = NULL, y = "Fst") +
  geom_hline(
    yintercept = c(upper_99.9thPercentile.fst, lower_99.9thPercentile.fst),
    linetype = "dashed",
    color = "red"
  ) +
  geom_point(
    data = subset(
      pixy_wide,
      avg_wc_fst > upper_99.9thPercentile.fst | avg_wc_fst < lower_99.9thPercentile.fst
    ),
    color = "red",
    size = point_size
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
    aes(xintercept = peak_cum),
    color = trait_colors[as.character(qtl_pixy$trait)],
    linewidth = 0.7
  ) +
  theme(axis.text = element_text(size = 40), axis.title = element_text(size = 40))

# =============================================================================
# Four-panel combo (5-panel from qmd minus pi)
# =============================================================================
combo_Dsuite_pixy_4panel_fdm_df_dxy_fst <-
  LASI_99.9thper_w50s25_f_dM /
  LASI_99.9thper_w50s25_d_f /
  Dxy_plot /
  wc_fst_plot +
  plot_layout(heights = c(0.55, 0.55, 1.45, 1.45)) +
  plot_annotation(tag_levels = "A", tag_suffix = ")") &
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.line = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size = 40, hjust = 0, vjust = 0)
  )

out_png <- file.path(project_root, "combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.png")
ggsave(
  out_png,
  combo_Dsuite_pixy_4panel_fdm_df_dxy_fst,
  device = "png",
  width = 120,
  height = 80,
  units = "cm",
  bg = "white"
)

message("Wrote ", out_png)
