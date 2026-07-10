# =============================================================================
# genome-wide_f_dM_dsuite_w_QTL_plots.R
# Along-genome f_dM (Dsuite w50s25) with 500 kb averages and QTL CIs — mirrors
# the Dxy multi-panel figure (Dxy_10kb_LASIout_separate_chrs_w.500kb.avg_7QTL_narrowed.CIs)
# but uses Dsuite local f_dM instead of pixy Dxy.
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

theme_set(theme_bw(base_size = 11))

# ---- Paths ----
project_root <- "/Users/admin-jgharenc/Documents/Github/Costus_popgen_clean"
dsuite_dir <- file.path(project_root, "Dsuite", "Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120")
qtl_ci_path <- file.path(project_root, "QTL", "all_QTL_CIs_bp_cm.251130.csv")

# ---- Dsuite: load window stats ----
file_list <- list.files(path = dsuite_dir, pattern = "*.txt", full.names = TRUE)
dsuite_wide <- data.frame()
for (file in file_list) {
  temp_df <- read.delim(file, sep = "\t")
  dsuite_wide <- rbind(dsuite_wide, temp_df)
}
dsuite_wide <- dsuite_wide %>%
  mutate(chrN = as.numeric(str_extract(chr, "\\d+")))

# ---- QTL and trait colors ----
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

# ---- 500 kb window averages (same logic as pixy Dxy in genome-wide_Fst.Dxy_w_QTL_plots.R) ----
window_size_bp <- 500000
min_frac_valid <- 0.25

stat_col <- "f_dM"

windowed_stat_avg <- dsuite_wide %>%
  mutate(bin_mid = (floor(windowStart / window_size_bp) + 0.5) * window_size_bp) %>%
  group_by(chrN, WindowMid = bin_mid) %>%
  summarise(
    n_total = n(),
    n_valid = sum(is.finite(.data[[stat_col]])),
    MeanY = mean(.data[[stat_col]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_total > 0, n_valid / n_total >= min_frac_valid) %>%
  select(chrN, WindowMid, MeanY)

windowed_stat_avg_list <- lapply(1:9, function(i) windowed_stat_avg %>% filter(chrN == i))
names(windowed_stat_avg_list) <- paste0("chr", 1:9)

# ---- Cutoff: 99.9th % of |stat| (symmetric tail highlight; raw signed y below) ----
cutoff99.9_abs_stat <- {
  x <- dsuite_wide[[stat_col]]
  x <- x[is.finite(x)]
  quantile(abs(x), 0.999)
}

plot_chr_combined <- function(raw_df, win_df, cutoff_abs, value_col = stat_col) {
  vc <- sym(value_col)
  raw_df <- raw_df %>% mutate(.abs_stat = abs(!!vc))
  ggplot(raw_df, aes(x = windowStart)) +
    geom_point(
      data = raw_df %>% filter(.abs_stat <= cutoff_abs),
      aes(y = !!vc),
      size = 1,
      alpha = 0.4,
      color = "grey40"
    ) +
    geom_point(
      data = raw_df %>% filter(.abs_stat > cutoff_abs),
      aes(y = !!vc),
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
      y = "f_dM"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)))
}

chroms <- 1:9

plots <- lapply(chroms, function(i) {
  raw_df_i <- dsuite_wide %>% filter(chrN == i, is.finite(.data[[stat_col]]))
  win_df_i <- windowed_stat_avg_list[[paste0("chr", i)]]
  plot_chr_combined(
    raw_df = raw_df_i,
    win_df = win_df_i,
    cutoff_abs = cutoff99.9_abs_stat
  )
})
names(plots) <- paste0("chr", chroms)

chrom_list <- sort(unique(dsuite_wide$chrN))

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

for (chr_num in chrom_list) {
  base_plot <- plots[[paste0("chr", chr_num)]]
  qtl_here <- all_QTL_CI %>% filter(chr == chr_num)
  plot_i <- add_qtl_layers(base_plot, qtl_here, trait_colors)
  plot_i <- plot_i + ggtitle(paste("Chromosome", chr_num))
  final_plots[[paste0("Chrom", chr_num)]] <- plot_i
}

f_dM_LASI_sep_chrs_w.QTL <- final_plots$Chrom1 + final_plots$Chrom2 + final_plots$Chrom3 +
  final_plots$Chrom5 + final_plots$Chrom6 + final_plots$Chrom7 + final_plots$Chrom8 + final_plots$Chrom9

if (!dir.exists(file.path(project_root, "Dsuite"))) dir.create(file.path(project_root, "Dsuite"))
ggsave(
  file.path(project_root, "Dsuite", "f_dM_w50s25_LASIout_separate_chrs_w.500kb.avg_7QTL_narrowed.CIs.png"),
  f_dM_LASI_sep_chrs_w.QTL,
  device = "png",
  width = 50,
  height = 20,
  units = "cm"
)
