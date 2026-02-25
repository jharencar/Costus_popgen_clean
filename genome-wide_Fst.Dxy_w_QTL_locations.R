# =============================================================================
# Genome-wide Fst and Dxy with QTL locations
# Extracted from Population_genetics_vill_alle.qmd (pixy section)
# =============================================================================
#
# Dependencies: tidyverse (readr, dplyr, tidyr), ggplot2, patchwork
# Paths: pixy_folder, all_QTL_CI CSV, and ggsave output dir "pixy/" may need
#        updating for your setup.
# =============================================================================

# ---- load dependencies ----
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

theme_set(theme_bw())

# ---- pixy: import data ----
# function from the pixy website (https://pixy.readthedocs.io/en/latest/plotting.html)
pixy_to_long <- function(pixy_files){

  pixy_df <- list()

  for(i in 1:length(pixy_files)){

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

# Import my data
pixy_folder <- "/Users/juliaharencar/Documents/Github/Population_genetics_vill_alle/pixy/output_251124_LASIref"
pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)
pixy_wide <- pivot_wider(pixy_df, names_from = statistic, values_from = value) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000)
# make chromosome number only col:
pixy_wide <- pixy_wide %>%
  mutate(chrN = as.numeric(str_extract(chromosome, "\\d+")))

# ---- windowed averages (500 kb windows, for red line) ----
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
plot_chr_combined <- function(raw_df, win_df) {
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
    labs(
      x = "Position (bp)",
      y = "Dxy"
    ) +
    theme_bw()
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

# QTL peaks and credible intervals (from rqtl script)
all_QTL_CI <- read.csv("/Users/juliaharencar/Documents/Github/alle_vill_QTL/all_QTL_CIs_bp_cm.251130.csv", header = TRUE)

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
Dxy_chr1_chr3_2panel <- final_plots$Chrom1 / final_plots$Chrom3
ggsave("pixy/Dxy_chr1_chr3_2panel.png", Dxy_chr1_chr3_2panel,
       device = "png", width = 22, height = 10, units = "cm")

Dxy_LASI_sep_chrs_w.QTL <- final_plots$Chrom1 + final_plots$Chrom2 + final_plots$Chrom3 +
  final_plots$Chrom5 + final_plots$Chrom6 + final_plots$Chrom7 + final_plots$Chrom8 + final_plots$Chrom9

ggsave("pixy/Dxy_10kb_LASIout_separate_chrs_w.500kb.avg_7QTL_narrowed.CIs.png",
       Dxy_LASI_sep_chrs_w.QTL, device = "png", width = 50, height = 20, units = "cm")

# ---- Fst: plot with QTL ----
plot_chr_combined <- function(raw_df, win_df) {
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
    labs(
      x = "Position (bp)",
      y = "Fst"
    ) +
    theme_bw()
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
Fst_chr1_chr3_2panel <- final_plots$Chrom1 / final_plots$Chrom3
ggsave("pixy/Fst_chr1_chr3_2panel.png", Fst_chr1_chr3_2panel,
       device = "png", width = 22, height = 10, units = "cm")

Fst_LASI_sep_chrs_w.QTL <- final_plots$Chrom1 + final_plots$Chrom2 + final_plots$Chrom3 +
  final_plots$Chrom5 + final_plots$Chrom6 + final_plots$Chrom7 + final_plots$Chrom8 + final_plots$Chrom9

ggsave("pixy/Fst_10kb_LASIout_separate_chrs_w.500kb.avg_7QTL_narrowed.CIs.png",
       Fst_LASI_sep_chrs_w.QTL, device = "png", width = 50, height = 20, units = "cm")
