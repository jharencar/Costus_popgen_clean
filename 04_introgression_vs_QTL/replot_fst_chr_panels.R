# Replot Fst chromosome panels only (styling iteration helper).
# Usage (from repo root):
#   Rscript 04_introgression_vs_QTL/replot_fst_chr_panels.R
#
# Writes the same Fst chromosome outputs as genome-wide_Fst.Dxy_w_QTL_plots.R
# without loading Dsuite, lepmap, or building the remaining figures.

library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggrastr)

theme_set(theme_bw(base_size = 11))
fst_raster_dpi <- 300
min_snps_per_10kb <- 10

pixy_to_long <- function(pixy_files) {
  pixy_df <- list()
  for (i in seq_along(pixy_files)) {
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    if (stat_file_type == "pi") {
      df <- read_delim(pixy_files[i], delim = "\t", show_col_types = FALSE)
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      pixy_df[[i]] <- df
    } else {
      df <- read_delim(pixy_files[i], delim = "\t", show_col_types = FALSE)
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df
    }
  }
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

read_pixy_wide <- function(folder, label) {
  if (is.na(folder) || !nzchar(folder) || !dir.exists(folder)) {
    stop(label, " pixy folder not found: ", folder)
  }
  pixy_files <- list.files(folder, pattern = "_(pi|dxy|fst)\\.txt$", full.names = TRUE)
  if (length(pixy_files) == 0) stop("No pixy tables in ", folder)
  wide <- pixy_to_long(pixy_files) %>%
    pivot_wider(names_from = statistic, values_from = value) %>%
    mutate(
      chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000,
      chrN = as.numeric(str_extract(chromosome, "\\d+"))
    )
  if (!"avg_hudson_fst" %in% names(wide)) {
    stop(label, " is missing avg_hudson_fst.")
  }
  num_cols <- intersect(
    c("avg_hudson_fst", "hudson_fst_num", "hudson_fst_den", "avg_dxy", "avg_pi", "no_snps", "no_sites"),
    names(wide)
  )
  wide %>% mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x))))
}

line_500kb <- function(wide_500, wide_10, value_col, count_col, min_count) {
  valid_10kb <- wide_10 %>%
    filter(
      is.finite(.data[[value_col]]),
      is.finite(.data[[count_col]]),
      .data[[count_col]] >= min_count
    ) %>%
    transmute(
      chrN,
      bin_start = floor((window_pos_1 - 1) / 500000) * 500000 + 1
    ) %>%
    count(chrN, bin_start, name = "n_valid")
  wide_500 %>%
    transmute(
      chrN,
      window_pos_1,
      WindowMid = (window_pos_1 + window_pos_2) / 2,
      MeanY = .data[[value_col]],
      n_expected = (window_pos_2 - window_pos_1 + 1) / 10000
    ) %>%
    left_join(valid_10kb, by = c("chrN", "window_pos_1" = "bin_start")) %>%
    mutate(
      n_valid = coalesce(n_valid, 0L),
      keep = is.finite(MeanY) & n_expected > 0 & (n_valid / n_expected) >= 0.5
    ) %>%
    filter(keep) %>%
    arrange(chrN, WindowMid) %>%
    distinct(chrN, WindowMid, .keep_all = TRUE) %>%
    group_by(chrN) %>%
    mutate(
      step = WindowMid - lag(WindowMid),
      line_run = cumsum(replace_na(step, 0) > 550000)
    ) %>%
    ungroup() %>%
    select(chrN, WindowMid, MeanY, line_run)
}

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
pixy_folder_500kb <- paths$pixy_folder_500kb
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
  dir.create(file.path(project_root, "pixy"), recursive = TRUE, showWarnings = FALSE)
}

pixy_wide <- read_pixy_wide(pixy_folder, "10 kb")
pixy_wide_500kb <- read_pixy_wide(pixy_folder_500kb, "500 kb")

windowed_Fst_avg <- line_500kb(
  pixy_wide_500kb, pixy_wide, "avg_hudson_fst", "no_snps", min_snps_per_10kb
)
windowed_Fst_avg_list <- lapply(1:9, function(i) windowed_Fst_avg %>% filter(chrN == i))
names(windowed_Fst_avg_list) <- paste0("chr", 1:9)

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
trait_colors <- c(
  SD = "#0072B2",
  GR = "#E69F00",
  tough = "#8250C4",
  Nreabs = "#CC79A7",
  Chlreabs = "#009E73"
)

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
      aes(xintercept = peak_bp, color = trait, linetype = qtl_linetype),
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = trait_colors, name = NULL, drop = FALSE) +
    scale_color_manual(values = trait_colors, name = NULL, drop = FALSE) +
    scale_linetype_identity(guide = "none") +
    theme(legend.position = "none")
}

pixy_fst_10kb <- pixy_wide %>%
  filter(is.finite(avg_hudson_fst), is.finite(no_snps), no_snps >= min_snps_per_10kb)
cutoff99_fst <- {
  x <- pixy_fst_10kb$avg_hudson_fst
  x <- pmax(x, 0)
  quantile(x[is.finite(x)], 0.99)
}
cat("Genome-wide 99th-percentile Fst cutoff: ", signif(cutoff99_fst, 4), "\n", sep = "")

# Grey points first; QTL CI/peak next; 99th line + red outliers + 500 kb avg on top
plot_chr_fst_base <- function(raw_df, cutoff99) {
  ggplot(raw_df, aes(x = window_pos_1)) +
    geom_point_rast(
      data = raw_df %>% filter(pmax(avg_hudson_fst, 0) <= cutoff99),
      aes(y = pmax(avg_hudson_fst, 0)),
      size = 1,
      alpha = 0.35,
      color = "grey50",
      raster.dpi = fst_raster_dpi
    ) +
    labs(x = "Position (bp)", y = "Fst") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(ylim = c(0, 1))
}

add_fst_foreground <- function(p, raw_df, win_df, cutoff99) {
  p +
    geom_hline(
      yintercept = cutoff99,
      linetype = "dashed",
      color = "#D62728",
      linewidth = 0.35
    ) +
    # Vector points only: ggrastr washes out sparse outliers in PDF (pale red).
    # size 1.5 ≈ grey rast size 1 visually; alpha=1 keeps #D62728 fully opaque.
    geom_point(
      data = raw_df %>% filter(pmax(avg_hudson_fst, 0) > cutoff99),
      aes(x = window_pos_1, y = pmax(avg_hudson_fst, 0)),
      inherit.aes = FALSE,
      size = 1.5,
      alpha = 1,
      stroke = 0,
      shape = 16,
      color = "#D62728"
    ) +
    geom_line(
      data = win_df,
      aes(x = WindowMid, y = MeanY, group = line_run),
      inherit.aes = FALSE,
      color = "darkred",
      linewidth = 0.7
    )
}

chrom_list <- sort(unique(pixy_wide$chrN))
final_plots <- list()
for (i in seq_along(chrom_list)) {
  chr_num <- chrom_list[i]
  raw_df_i <- pixy_fst_10kb %>% filter(chrN == chr_num)
  win_df_i <- windowed_Fst_avg_list[[paste0("chr", chr_num)]]
  qtl_here <- all_QTL_CI %>% filter(chr == chr_num)
  plot_i <- plot_chr_fst_base(raw_df_i, cutoff99_fst)
  plot_i <- add_qtl_layers_fst(plot_i, qtl_here, trait_colors)
  plot_i <- add_fst_foreground(plot_i, raw_df_i, win_df_i, cutoff99_fst)
  final_plots[[i]] <- plot_i + ggtitle(paste("Chr", chr_num))
}
names(final_plots) <- paste0("Chrom", chrom_list)

# Shared subset-panel styling (2-panel and 3-panel)
fst_ylim_subset <- 0.75
fst_subset_title <- theme(
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 15),
  plot.title = element_text(hjust = 0.5, margin = margin(t = 2, b = 4))
)
fst_subset_no_x <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)
fst_subset_y <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.7))
fst_subset_top <- function(p) {
  p + coord_cartesian(ylim = c(0, fst_ylim_subset)) + fst_subset_y +
    labs(x = NULL) + fst_subset_title + fst_subset_no_x +
    theme(plot.margin = margin(2, 4, 0, 4))
}
fst_subset_mid <- function(p) {
  p + coord_cartesian(ylim = c(0, fst_ylim_subset)) + fst_subset_y +
    labs(x = NULL) + fst_subset_title + fst_subset_no_x +
    theme(plot.margin = margin(4, 4, 0, 4))
}
fst_subset_bot <- function(p) {
  p + coord_cartesian(ylim = c(0, fst_ylim_subset)) + fst_subset_y +
    fst_subset_title +
    theme(
      plot.margin = margin(4, 4, 4, 4),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

Fst_chr1_chr3_2panel <- fst_subset_top(final_plots$Chrom1) / fst_subset_bot(final_plots$Chrom3)
ggsave(file.path(project_root, "pixy", "Fst_chr1_chr3_2panel.png"), Fst_chr1_chr3_2panel,
       device = "png", width = 22, height = 12, units = "cm", dpi = 200)
ggsave(file.path(project_root, "pixy", "Fst_chr1_chr3_2panel.pdf"), Fst_chr1_chr3_2panel,
       device = "pdf", width = 22, height = 12, units = "cm")

Fst_chr1_chr6_2panel <- fst_subset_top(final_plots$Chrom1) / fst_subset_bot(final_plots$Chrom6)
ggsave(file.path(project_root, "pixy", "Fst_chr1_chr6_2panel.png"), Fst_chr1_chr6_2panel,
       device = "png", width = 22, height = 12, units = "cm", dpi = 200)
ggsave(file.path(project_root, "pixy", "Fst_chr1_chr6_2panel.pdf"), Fst_chr1_chr6_2panel,
       device = "pdf", width = 22, height = 12, units = "cm")

Fst_chr1_chr2_chr3_3panel <- (
  fst_subset_top(final_plots$Chrom1) /
  fst_subset_mid(final_plots$Chrom2) /
  fst_subset_bot(final_plots$Chrom3)
)
ggsave(
  file.path(project_root, "pixy", "Fst_chr1_chr2_chr3_3panel.png"),
  Fst_chr1_chr2_chr3_3panel,
  device = "png", width = 22, height = 14, units = "cm", dpi = 200
)
ggsave(
  file.path(project_root, "pixy", "Fst_chr1_chr2_chr3_3panel.pdf"),
  Fst_chr1_chr2_chr3_3panel,
  device = "pdf", width = 22, height = 14, units = "cm"
)
message("Wrote Fst subset 2-panel and 3-panel figures (ylim 0–", fst_ylim_subset, ")")

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
message("Wrote ", fst_pixy_png)

# Single-column stack: same 8 chromosomes, wider panels (cowplot avoids
# patchwork nesting issues when combining ncol=1 panels + legend)
fst_1col_w_cm <- 22
fst_1col_h_cm <- 48
fst_panels_1col <- plot_grid(
  final_plots$Chrom1, final_plots$Chrom2, final_plots$Chrom3, final_plots$Chrom5,
  final_plots$Chrom6, final_plots$Chrom7, final_plots$Chrom8, final_plots$Chrom9,
  ncol = 1, align = "v"
)
Fst_LASI_sep_chrs_w.QTL_1col <- plot_grid(
  fst_panels_1col, fst_legend_grob,
  ncol = 1, rel_heights = c(1, 0.04)
)
fst_1col_pdf <- file.path(
  project_root, "pixy",
  "Fst_10kb_LASIout_separate_chrs_w.500kb.avg_7QTL_narrowed.CIs_1col.pdf"
)
fst_1col_png <- sub("\\.pdf$", ".png", fst_1col_pdf)
ggsave(fst_1col_pdf, Fst_LASI_sep_chrs_w.QTL_1col, device = "pdf",
       width = fst_1col_w_cm, height = fst_1col_h_cm, units = "cm")
ggsave(fst_1col_png, Fst_LASI_sep_chrs_w.QTL_1col, device = "png",
       width = fst_1col_w_cm, height = fst_1col_h_cm, units = "cm", dpi = 200)
message("Wrote ", fst_1col_png)

if (dir.exists(final_figures_dir)) {
  ok <- TRUE
  tryCatch(
    {
      ggsave(
        file.path(final_figures_dir, "sup_genomewide_Fst_QTL_CIs.pdf"),
        Fst_LASI_sep_chrs_w.QTL,
        device = "pdf",
        width = fst_fig_w_cm,
        height = fst_fig_h_cm,
        units = "cm"
      )
      ggsave(
        file.path(final_figures_dir, "sup_genomewide_Fst_QTL_CIs.png"),
        Fst_LASI_sep_chrs_w.QTL,
        device = "png",
        width = fst_fig_w_cm,
        height = fst_fig_h_cm,
        units = "cm",
        dpi = 200
      )
      ggsave(
        file.path(final_figures_dir, "sup_genomewide_Fst_QTL_CIs_1col.pdf"),
        Fst_LASI_sep_chrs_w.QTL_1col,
        device = "pdf",
        width = fst_1col_w_cm,
        height = fst_1col_h_cm,
        units = "cm"
      )
      ggsave(
        file.path(final_figures_dir, "sup_genomewide_Fst_QTL_CIs_1col.png"),
        Fst_LASI_sep_chrs_w.QTL_1col,
        device = "png",
        width = fst_1col_w_cm,
        height = fst_1col_h_cm,
        units = "cm",
        dpi = 200
      )
    },
    error = function(e) {
      ok <<- FALSE
      message("Skipped Final_figures write: ", conditionMessage(e))
    }
  )
  if (ok) message("Wrote 2-col and 1-col Fst panels to ", final_figures_dir)
}
