# =============================================================================
# combo_4panel_chem_alk_flav_QTL_vlines.R
#
# Same f_dM / d_f / Dxy / Fst layout + red outlier windows as
# combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R, but replaces ecological-trait QTL
# CI ribbons with vertical lines for alkaloid vs flavonoid chemistry QTL peaks
# from the WatersMass first-pass screen (10% FDR), clustered within 15 cM:
#   - solid  = shared peak region (≥2 traits)
#   - dashed = independent peak (1 trait)
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

source(file.path(project_root, "04_introgression_vs_QTL", "chem_alk_flav_QTL_regions.R"))
chem_colors <- chem_colors_alk_flav
point_size <- 3
combo_raster_dpi <- 300
legend_text_size <- 32

# =============================================================================
# Chemistry peaks
# =============================================================================
chem_regions <- load_chem_alk_flav_regions(project_root)
message("Chemistry QTL regions for overlay:")
print(as.data.frame(chem_regions))
out_tbl <- file.path(project_root, "04_introgression_vs_QTL", "chem_alk_flav_QTL_regions_for_4panel.csv")
write.csv(chem_regions, out_tbl, row.names = FALSE)

# =============================================================================
# Dsuite + pixy (same as combo 4-panel)
# =============================================================================
LASIout_w50s25_all9 <- data.frame()
file_list <- list.files(path = dsuite_dir, pattern = "*.txt", full.names = TRUE)
stopifnot(length(file_list) > 0)
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
  summarise(max_windowStart = max(windowStart), .groups = "drop") %>%
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
  summarize(center = mean(windowStart_cum), chr_label = first(chr_label), .groups = "drop")

dsuite_chr_add <- LASIout_w50s25_all9 %>%
  distinct(chr, windowStart_add) %>%
  mutate(chrN = as.numeric(gsub("Chrom", "", chr)))

chem_dsuite <- chem_regions %>%
  left_join(dsuite_chr_add, by = "chrN") %>%
  mutate(peak_cum = peak_bp + windowStart_add)

dsuite_xlim <- range(LASIout_w50s25_all9$windowStart_cum, na.rm = TRUE)

pixy_to_long <- function(pixy_files) {
  pixy_df <- list()
  for (i in seq_along(pixy_files)) {
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    if (stat_file_type == "pi") {
      df <- read_delim(pixy_files[i], delim = "\t", show_col_types = FALSE)
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
          key = "statistic", value = "value"
        ) %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      pixy_df[[i]] <- df
    } else {
      df <- read_delim(pixy_files[i], delim = "\t", show_col_types = FALSE)
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
stopifnot(length(pixy_files) > 0)
pixy_df <- pixy_to_long(pixy_files)
pixy_wide <- pivot_wider(pixy_df, names_from = statistic, values_from = value) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000)

pixy_data_cum <- pixy_wide %>%
  group_by(chromosome) %>%
  summarise(max_windowStart = max(window_pos_1), .groups = "drop") %>%
  mutate(windowStart_add = lag(cumsum(max_windowStart), default = 0)) %>%
  select(chromosome, windowStart_add)

pixy_wide <- pixy_wide %>%
  inner_join(pixy_data_cum, by = "chromosome") %>%
  mutate(windowStart_cum = window_pos_1 + windowStart_add)

axis_set_pixy <- pixy_wide %>%
  group_by(chromosome) %>%
  summarize(center = mean(windowStart_cum), .groups = "drop")

pixy_chr_add <- pixy_wide %>%
  distinct(chromosome, windowStart_add) %>%
  mutate(chrN = as.numeric(gsub("Chrom", "", chromosome)))

chem_pixy <- chem_regions %>%
  left_join(pixy_chr_add, by = "chrN") %>%
  mutate(peak_cum = peak_bp + windowStart_add)

pixy_xlim <- range(pixy_wide$windowStart_cum, na.rm = TRUE)
combo_xlim <- c(
  min(dsuite_xlim[1], pixy_xlim[1]),
  max(dsuite_xlim[2], pixy_xlim[2])
)
x_expand <- expansion(mult = 0.012, add = 0)

# Hard-coded colors per group so we do not overwrite chromosome point colors
add_chem_vlines <- function(p, chem_df) {
  for (g in names(chem_colors)) {
    sub <- chem_df %>% filter(as.character(group) == g)
    if (nrow(sub) == 0) next
    p <- p +
      geom_vline(
        data = sub,
        aes(xintercept = peak_cum, linetype = sharing),
        color = chem_colors[[g]],
        linewidth = 0.9,
        inherit.aes = FALSE
      )
  }
  p +
    scale_linetype_manual(
      values = c(shared = "solid", independent = "dashed"),
      drop = FALSE,
      guide = "none"
    )
}

panel_theme_big <- theme(
  axis.text = element_text(size = 30),
  axis.title = element_text(size = 30),
  plot.margin = unit(c(2, 1, 1, 1), "cm")
)
panel_theme <- theme(
  axis.text = element_text(size = 40),
  axis.title = element_text(size = 40)
)

LASI_99.9thper_w50s25_f_dM <- ggplot(
  data = LASIout_w50s25_all9,
  aes(x = windowStart_cum, y = f_dM, color = as_factor(chr))
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_point_rast(size = point_size, raster.dpi = combo_raster_dpi) +
  geom_hline(
    yintercept = c(cutoff99.9thPercentile, -cutoff99.9thPercentile),
    linetype = "dashed",
    color = "red"
  ) +
  scale_color_manual(values = rep(c("black", "darkgray"), unique(length(axis_set_dsuite$chr)))) +
  guides(color = "none") +
  geom_point_rast(
    data = subset(
      LASIout_w50s25_all9,
      f_dM > cutoff99.9thPercentile | f_dM < -cutoff99.9thPercentile
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
  panel_theme_big
LASI_99.9thper_w50s25_f_dM <- add_chem_vlines(LASI_99.9thper_w50s25_f_dM, chem_dsuite)

cutoff99.9thPercentile_d_f <- quantile(
  abs(LASIout_w50s25_all9$d_f[which(LASIout_w50s25_all9$d_f < 0)]),
  0.999
)

LASI_99.9thper_w50s25_d_f <- ggplot(
  data = LASIout_w50s25_all9,
  aes(x = windowStart_cum, y = d_f, color = as_factor(chr))
) +
  geom_point_rast(size = point_size, raster.dpi = combo_raster_dpi) +
  geom_hline(
    yintercept = c(cutoff99.9thPercentile_d_f, -cutoff99.9thPercentile_d_f),
    linetype = "dashed",
    color = "red"
  ) +
  scale_color_manual(values = rep(c("black", "darkgray"), unique(length(axis_set_dsuite$chr)))) +
  guides(color = "none") +
  geom_point_rast(
    data = subset(
      LASIout_w50s25_all9,
      d_f > cutoff99.9thPercentile_d_f | d_f < -cutoff99.9thPercentile_d_f
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
  panel_theme
LASI_99.9thper_w50s25_d_f <- add_chem_vlines(LASI_99.9thper_w50s25_d_f, chem_dsuite)

cutoff99.9thPercentile.dxy <- quantile(pixy_wide$avg_dxy, 0.999, na.rm = TRUE)
Dxy_plot <- ggplot(data = pixy_wide, aes(x = windowStart_cum, y = avg_dxy, color = as_factor(chromosome))) +
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
  geom_hline(yintercept = cutoff99.9thPercentile.dxy, linetype = "dashed", color = "red") +
  geom_point_rast(
    data = subset(pixy_wide, avg_dxy > cutoff99.9thPercentile.dxy),
    color = "red",
    size = point_size,
    raster.dpi = combo_raster_dpi
  ) +
  panel_theme
Dxy_plot <- add_chem_vlines(Dxy_plot, chem_pixy)

upper_99.9thPercentile.fst <- quantile(pixy_wide$avg_wc_fst, 0.999, na.rm = TRUE)
lower_99.9thPercentile.fst <- quantile(pixy_wide$avg_wc_fst, 0.001, na.rm = TRUE)
wc_fst_plot <- ggplot(data = pixy_wide, aes(x = windowStart_cum, y = avg_wc_fst, color = as_factor(chromosome))) +
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
    yintercept = c(upper_99.9thPercentile.fst, lower_99.9thPercentile.fst),
    linetype = "dashed",
    color = "red"
  ) +
  geom_point_rast(
    data = subset(
      pixy_wide,
      avg_wc_fst > upper_99.9thPercentile.fst | avg_wc_fst < lower_99.9thPercentile.fst
    ),
    color = "red",
    size = point_size,
    raster.dpi = combo_raster_dpi
  ) +
  panel_theme
wc_fst_plot <- add_chem_vlines(wc_fst_plot, chem_pixy)

# =============================================================================
# Legend: color = chemical group; linetype = shared vs independent
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

legend_df <- tidyr::crossing(
  group = factor(c("Alkaloid", "Flavonoid"), levels = c("Alkaloid", "Flavonoid")),
  sharing = factor(c("shared", "independent"), levels = c("shared", "independent"))
) %>%
  mutate(
    label = paste(as.character(group), as.character(sharing), sep = " / "),
    x = as.numeric(group) + (as.numeric(sharing) - 1) * 0.15
  )

chem_legend_plot <- ggplot(legend_df, aes(x = x, y = 0, color = group, linetype = sharing)) +
  geom_vline(aes(xintercept = x, color = group, linetype = sharing), linewidth = 1.4) +
  scale_color_manual(values = chem_colors, name = NULL) +
  scale_linetype_manual(
    values = c(shared = "solid", independent = "dashed"),
    name = NULL,
    labels = c(shared = "shared peak", independent = "independent peak")
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 1.4, linetype = "solid")),
    linetype = guide_legend(order = 2, override.aes = list(color = "black", linewidth = 1.4))
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = legend_text_size),
    legend.key.height = unit(1.0, "cm"),
    legend.key.width = unit(2.2, "cm"),
    legend.title = element_blank(),
    legend.box = "horizontal"
  )
chem_legend_grob <- get_legend(chem_legend_plot)

combo_chem <-
  LASI_99.9thper_w50s25_f_dM +
  LASI_99.9thper_w50s25_d_f +
  Dxy_plot +
  wc_fst_plot +
  wrap_elements(full = chem_legend_grob) +
  plot_layout(
    design = "
      AAAA
      BBBB
      CCCC
      DDDD
      EEEE
    ",
    heights = c(0.55, 0.55, 1.45, 1.45, 0.1),
    guides = "keep"
  ) +
  plot_annotation(tag_levels = "A", tag_suffix = ")") &
  combo_panel_theme

out_dir <- file.path(project_root, "04_introgression_vs_QTL")
out_png <- file.path(out_dir, "sup_4panel_fdm_df_dxy_fst_chem_alk_flav_QTL.png")
out_pdf <- file.path(out_dir, "sup_4panel_fdm_df_dxy_fst_chem_alk_flav_QTL.pdf")
# Also write informal copies at repo root for easy find
out_png_root <- file.path(project_root, "sup_4panel_fdm_df_dxy_fst_chem_alk_flav_QTL.png")
out_pdf_root <- file.path(project_root, "sup_4panel_fdm_df_dxy_fst_chem_alk_flav_QTL.pdf")

for (f in c(out_png, out_png_root)) {
  ggsave(f, combo_chem, device = "png", width = 120, height = 80, units = "cm", bg = "white")
  message("Wrote ", f)
}
for (f in c(out_pdf, out_pdf_root)) {
  ggsave(f, combo_chem, device = "pdf", width = 100, height = 70, units = "cm", bg = "white")
  message("Wrote ", f)
}
message("Region table: ", out_tbl)
