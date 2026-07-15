# =============================================================================
# genomewide_Fst_chem_alk_flav_QTL_vlines.R
#
# Same layout as sup_genomewide_Fst_QTL_CIs.pdf (per-chr Fst 10 kb points +
# 500 kb mean line + red 99.9th outliers; 2x4 grid, no chr4), but replaces
# ecological QTL CI ribbons with alkaloid / flavonoid chemistry QTL peak lines:
#   solid = shared (≥2 traits within 15 cM); dashed = independent.
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggrastr)

theme_set(theme_bw(base_size = 11))
fst_raster_dpi <- 300

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

source(file.path(project_root, "04_introgression_vs_QTL", "chem_alk_flav_QTL_regions.R"))

chem_regions <- load_chem_alk_flav_regions(project_root)
message("Chemistry QTL regions:")
print(as.data.frame(chem_regions))
write.csv(
  chem_regions,
  file.path(project_root, "04_introgression_vs_QTL", "chem_alk_flav_QTL_regions_for_Fst_panels.csv"),
  row.names = FALSE
)
if (any(chem_regions$chrN == 4)) {
  message(
    "NOTE: chr4 alkaloid peak(s) exist but this supplement layout omits chr4 ",
    "(same as published Fst QTL CI figure). See region table."
  )
}

# ---- pixy ----
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
pixy_wide <- pivot_wider(pixy_to_long(pixy_files), names_from = statistic, values_from = value) %>%
  mutate(
    chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000,
    chrN = as.numeric(str_extract(chromosome, "\\d+"))
  )

window_size_bp <- 500000
min_frac_valid <- 0.25
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

windowed_Fst_avg_list <- lapply(1:9, function(i) windowed_Fst_avg %>% filter(chrN == i))
names(windowed_Fst_avg_list) <- paste0("chr", 1:9)

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
    labs(x = "Position (bp)", y = "Fst") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)))
}

# Same chromosomes as published Fst QTL CI supplement (exclude chr4)
chrom_list <- c(1, 2, 3, 5, 6, 7, 8, 9)

final_plots <- list()
for (chr_num in chrom_list) {
  raw_df_i <- pixy_wide %>% filter(chrN == chr_num, is.finite(avg_wc_fst))
  win_df_i <- windowed_Fst_avg_list[[paste0("chr", chr_num)]]
  base_plot <- plot_chr_combined(raw_df_i, win_df_i, cutoff99.9_fst)
  chem_here <- chem_regions %>% filter(chrN == chr_num)
  plot_i <- add_chem_vlines_bp(base_plot, chem_here) + ggtitle(paste("Chr", chr_num))
  final_plots[[paste0("Chrom", chr_num)]] <- plot_i
}

fst_sup_theme <- theme(
  plot.title = element_text(size = 10, hjust = 0.5),
  axis.title = element_text(size = 9),
  axis.text = element_text(size = 8),
  legend.position = "none"
)
for (nm in names(final_plots)) {
  final_plots[[nm]] <- final_plots[[nm]] + fst_sup_theme
}

# Legend: color = group; linetype = shared vs independent
legend_df <- tidyr::crossing(
  group = factor(c("Alkaloid", "Flavonoid"), levels = c("Alkaloid", "Flavonoid")),
  sharing = factor(c("shared", "independent"), levels = c("shared", "independent"))
) %>%
  mutate(x = as.numeric(group) + (as.numeric(sharing) - 1) * 0.12)

chem_legend_plot <- ggplot(legend_df, aes(x = x, y = 0, color = group, linetype = sharing)) +
  geom_vline(aes(xintercept = x, color = group, linetype = sharing), linewidth = 1.2) +
  scale_color_manual(values = chem_colors_alk_flav, name = NULL) +
  scale_linetype_manual(
    values = c(shared = "solid", independent = "dashed"),
    name = NULL,
    labels = c(shared = "shared peak", independent = "independent peak")
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(linewidth = 1.2, linetype = "solid")),
    linetype = guide_legend(order = 2, override.aes = list(color = "black", linewidth = 1.2))
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.title = element_blank(),
    legend.box = "horizontal"
  )
fst_legend_grob <- get_legend(chem_legend_plot)

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

Fst_chem_QTL <- fst_panels / wrap_elements(full = fst_legend_grob) +
  plot_layout(heights = c(1, 0.07), guides = "keep")

fst_fig_w_cm <- 22
fst_fig_h_cm <- 28
out_dir <- file.path(project_root, "04_introgression_vs_QTL")
out_pdf <- file.path(out_dir, "sup_genomewide_Fst_chem_alk_flav_QTL.pdf")
out_png <- sub("\\.pdf$", ".png", out_pdf)
out_pdf_root <- file.path(project_root, "sup_genomewide_Fst_chem_alk_flav_QTL.pdf")
out_png_root <- sub("\\.pdf$", ".png", out_pdf_root)

for (f in c(out_pdf, out_pdf_root)) {
  ggsave(f, Fst_chem_QTL, device = "pdf", width = fst_fig_w_cm, height = fst_fig_h_cm, units = "cm")
  message("Wrote ", f)
}
for (f in c(out_png, out_png_root)) {
  ggsave(f, Fst_chem_QTL, device = "png", width = fst_fig_w_cm, height = fst_fig_h_cm, units = "cm", dpi = 200)
  message("Wrote ", f)
}
