suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

theme_set(theme_bw(base_size = 11))

# ---- Paths ----
project_root <- "/Users/juliaharencar/Documents/Github/Costus_popgen_clean"
dsuite_input_dir <- "/Users/juliaharencar/Documents/Github/Population_genetics_vill_alle/Dsuite/Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120"
lepmap_dir <- "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/best_iter_ordermarkers/"

# ---- Parameters ----
recomb_window_bp <- 1e6

# ---- LepMap: read map and build bp -> cM interpolation ----
lepmap_files <- list.files(
  lepmap_dir,
  pattern = "*_clean_LodLim11_theta0.001_allF2s_0.5covCO_withRef.txt",
  full.names = TRUE
)

lepmap_df <- do.call(rbind, lapply(lepmap_files, function(f) {
  markers <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  names(markers) <- c("marker_num", "cM_pos", "chrom", "bp_pos")
  markers
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
  df <- df[order(df[[bp_col]]), , drop = FALSE]
  iso <- isoreg(df[[cm_col]])
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

chrom_bounds <- bp_to_cm_iso %>%
  group_by(chrom) %>%
  summarise(
    chr_min_bp = min(bp_pos),
    chr_max_bp = max(bp_pos),
    .groups = "drop"
  ) %>%
  rename(chrN = chrom)

# ---- Dsuite: read window stats ----
dsuite_files <- list.files(dsuite_input_dir, pattern = "*.txt", full.names = TRUE)

dsuite_stats <- do.call(rbind, lapply(dsuite_files, function(file) {
  read.delim(file, sep = "\t")
})) %>%
  mutate(
    chrN = as.numeric(str_extract(chr, "\\d+")),
    window_mid_bp = (windowStart + windowEnd) / 2,
    recomb_bin_start = floor(window_mid_bp / recomb_window_bp) * recomb_window_bp
  ) %>%
  filter(chrN != 4)

# ---- Fixed 1 Mb recombination windows ----
# Assign each Dsuite window to a fixed 1 Mb physical bin by midpoint, then
# estimate recombination rate from the cM span across that full bin.
recomb_bins <- dsuite_stats %>%
  distinct(chrN, recomb_bin_start) %>%
  left_join(chrom_bounds, by = "chrN") %>%
  mutate(
    bin_start_bp = pmax(recomb_bin_start, chr_min_bp),
    bin_end_bp = pmin(recomb_bin_start + recomb_window_bp, chr_max_bp)
  ) %>%
  rowwise() %>%
  mutate(
    cM_start = bp_to_cm_funs[[as.character(chrN)]](bin_start_bp),
    cM_end = bp_to_cm_funs[[as.character(chrN)]](bin_end_bp),
    bin_width_bp = bin_end_bp - bin_start_bp,
    recomb_cM_Mb = if_else(
      bin_width_bp > 0,
      abs(cM_end - cM_start) / bin_width_bp * 1e6,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  filter(is.finite(recomb_cM_Mb))

dsuite_with_recomb <- dsuite_stats %>%
  left_join(
    recomb_bins %>% select(chrN, recomb_bin_start, recomb_cM_Mb),
    by = c("chrN", "recomb_bin_start")
  ) %>%
  filter(is.finite(recomb_cM_Mb))

dsuite_with_recomb <- dsuite_with_recomb %>%
  mutate(
    recomb_decile_num = ntile(recomb_cM_Mb, 10),
    recomb_decile = factor(
      recomb_decile_num,
      levels = 1:10,
      labels = paste0("D", 1:10)
    )
  )

selected_chroms <- c(1, 3, 5, 7, 8, 9)

dsuite_with_recomb_selected <- dsuite_with_recomb %>%
  filter(chrN %in% selected_chroms) %>%
  mutate(
    recomb_decile_num = ntile(recomb_cM_Mb, 10),
    recomb_decile = factor(
      recomb_decile_num,
      levels = 1:10,
      labels = paste0("D", 1:10)
    )
  )

dsuite_with_recomb_selected_by_chr_decile <- dsuite_with_recomb_selected %>%
  group_by(chrN) %>%
  mutate(
    recomb_decile_chr_num = ntile(recomb_cM_Mb, 10),
    recomb_decile_chr = factor(
      recomb_decile_chr_num,
      levels = 1:10,
      labels = paste0("D", 1:10)
    )
  ) %>%
  ungroup()

# ---- Diagnostic summaries ----
recomb_diagnostic_overall <- recomb_bins %>%
  summarise(
    scope = "all_chromosomes",
    n_bins = n(),
    min_cM_Mb = min(recomb_cM_Mb),
    p05_cM_Mb = quantile(recomb_cM_Mb, 0.05)[[1]],
    median_cM_Mb = median(recomb_cM_Mb),
    mean_cM_Mb = mean(recomb_cM_Mb),
    p95_cM_Mb = quantile(recomb_cM_Mb, 0.95)[[1]],
    max_cM_Mb = max(recomb_cM_Mb)
  )

recomb_diagnostic_by_chr <- recomb_bins %>%
  group_by(chrN) %>%
  summarise(
    scope = paste0("chr", first(chrN)),
    n_bins = n(),
    min_cM_Mb = min(recomb_cM_Mb),
    p05_cM_Mb = quantile(recomb_cM_Mb, 0.05)[[1]],
    median_cM_Mb = median(recomb_cM_Mb),
    mean_cM_Mb = mean(recomb_cM_Mb),
    p95_cM_Mb = quantile(recomb_cM_Mb, 0.95)[[1]],
    max_cM_Mb = max(recomb_cM_Mb),
    .groups = "drop"
  ) %>%
  select(-chrN)

recomb_diagnostic_table <- bind_rows(
  recomb_diagnostic_overall,
  recomb_diagnostic_by_chr
)

# ---- Plots ----
recomb_vs_fdM_plot <- ggplot(
  dsuite_with_recomb %>% filter(f_dM > 0),
  aes(x = recomb_cM_Mb, y = f_dM)
) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  labs(x = "Recombination rate (cM/Mb; fixed 1 Mb window)", y = "f_dM") +
  theme(panel.grid.minor = element_blank())

recomb_vs_fdM_all_plot <- ggplot(
  dsuite_with_recomb %>% filter(is.finite(f_dM)),
  aes(x = recomb_cM_Mb, y = f_dM)
) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  labs(x = "Recombination rate (cM/Mb; fixed 1 Mb window)", y = "f_dM") +
  theme(panel.grid.minor = element_blank())

recomb_vs_f_d_plot <- ggplot(
  dsuite_with_recomb %>% filter(is.finite(f_d)),
  aes(x = recomb_cM_Mb, y = f_d)
) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  labs(x = "Recombination rate (cM/Mb; fixed 1 Mb window)", y = "f_d") +
  theme(panel.grid.minor = element_blank())

recomb_vs_d_f_plot <- ggplot(
  dsuite_with_recomb %>% filter(is.finite(d_f)),
  aes(x = recomb_cM_Mb, y = d_f)
) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  labs(x = "Recombination rate (cM/Mb; fixed 1 Mb window)", y = "d_f") +
  theme(panel.grid.minor = element_blank())

recomb_by_stat_3panel <- recomb_vs_fdM_plot / recomb_vs_f_d_plot / recomb_vs_d_f_plot

plot_stat_vs_recomb <- function(df, stat_col, y_lab, positive_only = FALSE) {
  plot_df <- df[is.finite(df[[stat_col]]), , drop = FALSE]

  if (positive_only) {
    plot_df <- plot_df[plot_df[[stat_col]] > 0, , drop = FALSE]
  }

  ggplot(
    plot_df,
    aes(x = .data[["recomb_cM_Mb"]], y = .data[[stat_col]])
  ) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
    labs(x = "Recombination rate (cM/Mb; fixed 1 Mb window)", y = y_lab) +
    theme(panel.grid.minor = element_blank())
}

plot_stat_by_recomb_decile <- function(df, stat_col, y_lab) {
  plot_df <- df[is.finite(df[[stat_col]]), , drop = FALSE]

  ggplot(
    plot_df,
    aes(x = .data[["recomb_decile"]], y = .data[[stat_col]])
  ) +
    geom_boxplot(outlier.alpha = 0.2, width = 0.7) +
    labs(
      x = "Recombination rate decile",
      y = y_lab
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

plot_stat_by_recomb_decile_col <- function(df, stat_col, y_lab, decile_col) {
  plot_df <- df[is.finite(df[[stat_col]]), , drop = FALSE]

  ggplot(
    plot_df,
    aes(x = .data[[decile_col]], y = .data[[stat_col]])
  ) +
    geom_boxplot(outlier.alpha = 0.2, width = 0.7) +
    labs(
      x = "Recombination rate decile",
      y = y_lab
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

fdM_by_recomb_decile_boxplot <- plot_stat_by_recomb_decile(
  dsuite_with_recomb,
  "f_dM",
  "f_dM"
)

f_d_by_recomb_decile_boxplot <- plot_stat_by_recomb_decile(
  dsuite_with_recomb,
  "f_d",
  "f_d"
)

d_f_by_recomb_decile_boxplot <- plot_stat_by_recomb_decile(
  dsuite_with_recomb,
  "d_f",
  "d_f"
)

recomb_decile_boxplots_3panel <- fdM_by_recomb_decile_boxplot /
  f_d_by_recomb_decile_boxplot /
  d_f_by_recomb_decile_boxplot

# ---- Selected chromosomes only: 1, 3, 5, 7, 8, 9 ----
recomb_vs_fdM_plot_selected <- plot_stat_vs_recomb(
  dsuite_with_recomb_selected,
  "f_dM",
  "f_dM",
  positive_only = TRUE
)

recomb_vs_fdM_all_plot_selected <- plot_stat_vs_recomb(
  dsuite_with_recomb_selected,
  "f_dM",
  "f_dM"
)

recomb_vs_f_d_plot_selected <- plot_stat_vs_recomb(
  dsuite_with_recomb_selected,
  "f_d",
  "f_d"
)

recomb_vs_d_f_plot_selected <- plot_stat_vs_recomb(
  dsuite_with_recomb_selected,
  "d_f",
  "d_f"
)

recomb_by_stat_3panel_selected <- recomb_vs_fdM_plot_selected /
  recomb_vs_f_d_plot_selected /
  recomb_vs_d_f_plot_selected

fdM_by_recomb_decile_boxplot_selected <- plot_stat_by_recomb_decile(
  dsuite_with_recomb_selected,
  "f_dM",
  "f_dM"
)

f_d_by_recomb_decile_boxplot_selected <- plot_stat_by_recomb_decile(
  dsuite_with_recomb_selected,
  "f_d",
  "f_d"
)

d_f_by_recomb_decile_boxplot_selected <- plot_stat_by_recomb_decile(
  dsuite_with_recomb_selected,
  "d_f",
  "d_f"
)

recomb_decile_boxplots_3panel_selected <- fdM_by_recomb_decile_boxplot_selected /
  f_d_by_recomb_decile_boxplot_selected /
  d_f_by_recomb_decile_boxplot_selected

# ---- Selected chromosomes only: f_dM per chromosome multipanel PDFs ----
fdM_raw_by_chr_plots <- lapply(selected_chroms, function(chr_num) {
  plot_stat_vs_recomb(
    dsuite_with_recomb_selected[dsuite_with_recomb_selected$chrN == chr_num, , drop = FALSE],
    "f_dM",
    "f_dM"
  ) +
    ggtitle(paste("Chromosome", chr_num))
})

fdM_raw_by_chr_pdf <- wrap_plots(fdM_raw_by_chr_plots, ncol = 2)

fdM_decile_by_chr_plots <- lapply(selected_chroms, function(chr_num) {
  plot_stat_by_recomb_decile_col(
    dsuite_with_recomb_selected_by_chr_decile[
      dsuite_with_recomb_selected_by_chr_decile$chrN == chr_num,
      ,
      drop = FALSE
    ],
    "f_dM",
    "f_dM",
    "recomb_decile_chr"
  ) +
    ggtitle(paste("Chromosome", chr_num))
})

fdM_decile_by_chr_pdf <- wrap_plots(fdM_decile_by_chr_plots, ncol = 2)

# ---- Save outputs ----
if (!dir.exists(file.path(project_root, "Dsuite"))) {
  dir.create(file.path(project_root, "Dsuite"))
}

ggsave(
  file.path(project_root, "Dsuite", "recombination_rate_vs_f_dM.png"),
  plot = recomb_vs_fdM_plot,
  device = "png", width = 7, height = 5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "recombination_rate_vs_f_dM_all_values.png"),
  plot = recomb_vs_fdM_all_plot,
  device = "png", width = 7, height = 5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "recombination_rate_vs_f_d.png"),
  plot = recomb_vs_f_d_plot,
  device = "png", width = 7, height = 5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "recombination_rate_vs_d_f.png"),
  plot = recomb_vs_d_f_plot,
  device = "png", width = 7, height = 5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "recombination_rate_vs_Dsuite_stats_3panel.png"),
  plot = recomb_by_stat_3panel,
  device = "png", width = 7, height = 15, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "f_dM_by_recombination_rate_decile_boxplot.png"),
  plot = fdM_by_recomb_decile_boxplot,
  device = "png", width = 8, height = 5.5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "f_d_by_recombination_rate_decile_boxplot.png"),
  plot = f_d_by_recomb_decile_boxplot,
  device = "png", width = 8, height = 5.5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "d_f_by_recombination_rate_decile_boxplot.png"),
  plot = d_f_by_recomb_decile_boxplot,
  device = "png", width = 8, height = 5.5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "introgression_stats_by_recombination_rate_decile_3panel.png"),
  plot = recomb_decile_boxplots_3panel,
  device = "png", width = 8, height = 16.5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "recombination_rate_vs_f_dM_chr135789.png"),
  plot = recomb_vs_fdM_plot_selected,
  device = "png", width = 7, height = 5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "recombination_rate_vs_f_dM_all_values_chr135789.png"),
  plot = recomb_vs_fdM_all_plot_selected,
  device = "png", width = 7, height = 5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "recombination_rate_vs_f_d_chr135789.png"),
  plot = recomb_vs_f_d_plot_selected,
  device = "png", width = 7, height = 5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "recombination_rate_vs_d_f_chr135789.png"),
  plot = recomb_vs_d_f_plot_selected,
  device = "png", width = 7, height = 5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "recombination_rate_vs_Dsuite_stats_3panel_chr135789.png"),
  plot = recomb_by_stat_3panel_selected,
  device = "png", width = 7, height = 15, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "f_dM_by_recombination_rate_decile_boxplot_chr135789.png"),
  plot = fdM_by_recomb_decile_boxplot_selected,
  device = "png", width = 8, height = 5.5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "f_d_by_recombination_rate_decile_boxplot_chr135789.png"),
  plot = f_d_by_recomb_decile_boxplot_selected,
  device = "png", width = 8, height = 5.5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "d_f_by_recombination_rate_decile_boxplot_chr135789.png"),
  plot = d_f_by_recomb_decile_boxplot_selected,
  device = "png", width = 8, height = 5.5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "introgression_stats_by_recombination_rate_decile_3panel_chr135789.png"),
  plot = recomb_decile_boxplots_3panel_selected,
  device = "png", width = 8, height = 16.5, units = "cm", dpi = 300
)

ggsave(
  file.path(project_root, "Dsuite", "f_dM_raw_by_chromosome_recombination_chr135789.pdf"),
  plot = fdM_raw_by_chr_pdf,
  device = "pdf", width = 18, height = 24, units = "cm"
)

ggsave(
  file.path(project_root, "Dsuite", "f_dM_deciles_by_chromosome_recombination_chr135789.pdf"),
  plot = fdM_decile_by_chr_pdf,
  device = "pdf", width = 18, height = 24, units = "cm"
)

write.csv(
  recomb_diagnostic_table,
  file.path(project_root, "Dsuite", "recombination_rate_1Mb_diagnostic_summary.csv"),
  row.names = FALSE
)
