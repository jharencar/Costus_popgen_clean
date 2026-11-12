# =============================================================================
# genome-wide_Fst.Dxy_w_QTL_plots.R
# Genome-wide Fst, Dxy, and f_dM with QTL locations and distribution histograms
# =============================================================================
#
# Single script: (1) imports Hudson pixy (10 kb and 500 kb) + Dsuite f_dM once,
# (2) along-genome Fst/Dxy plots (10 kb points with >=10 SNPs or >=1000 sites;
# 500 kb line is pixy's own window, kept when >=50% of the 10 kb slots pass
# that same filter)
# with QTL peaks/CIs (antagonistic peaks dashed), (3) 1 cM distributions (Fst = ratio of Hudson components
# over 10 kb windows with >=10 SNPs; Dxy is the mean of windows with >=1000 sites; f_dM stays an unfiltered mean)
# with QTL summaries, density curves, tail shading;
# (4) simplified 1 cM Fst density; (5) 10 kb Fst histogram with QTL lines at
# the max 10 kb window in each 95% CI; (6) 10 kb Fst, Dxy, and f_dM distributions
# with a density curve for the windows inside each QTL 95% CI (dashed if antagonistic).
# Those 10 kb Fst/Dxy distributions keep windows with >=10 SNPs (Fst) or >=1000
# sites (Dxy). (7) 0.1 cM distributions built from those filtered 10 kb windows
# (f_dM stays a mean of positive Dsuite windows; it has no extra filter).
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

# Ratio of averages across windows (Bhatia et al. 2013): sum of Hudson
# numerators over sum of Hudson denominators, not the mean of window Fst.
hudson_fst_ratio <- function(num, den) {
  num <- sum(num, na.rm = TRUE)
  den <- sum(den, na.rm = TRUE)
  if (!is.finite(den) || den <= 0) NA_real_ else num / den
}

# Every Fst and Dxy figure below uses these 10 kb floors.
min_snps_per_10kb <- 10
min_dxy_sites_per_10kb <- 1000

read_pixy_wide <- function(folder, label) {
  if (is.na(folder) || !nzchar(folder) || !dir.exists(folder)) {
    stop(
      label, " pixy folder not found: ", folder,
      "\nSet PIXY_FOLDER / PIXY_FOLDER_500KB, or place pixy_hudson_10kb and pixy_hudson_500kb under the archive."
    )
  }
  pixy_files <- list.files(
    folder,
    pattern = "_(pi|dxy|fst)\\.txt$",
    full.names = TRUE
  )
  if (length(pixy_files) == 0) {
    stop("No pixy *_pi.txt / *_dxy.txt / *_fst.txt files in ", folder)
  }
  wide <- pixy_to_long(pixy_files) %>%
    pivot_wider(names_from = statistic, values_from = value) %>%
    mutate(
      chr_position = ((window_pos_1 + window_pos_2) / 2) / 1000000,
      chrN = as.numeric(str_extract(chromosome, "\\d+"))
    )
  if (!"avg_hudson_fst" %in% names(wide)) {
    if ("avg_wc_fst" %in% names(wide)) {
      stop(
        label, " is Weir–Cockerham (avg_wc_fst). ",
        "Use the Hudson pixy run, not output_251124_LASIref."
      )
    }
    stop(label, " is missing avg_hudson_fst.")
  }
  num_cols <- intersect(
    c("avg_hudson_fst", "hudson_fst_num", "hudson_fst_den", "avg_dxy", "avg_pi", "no_snps", "no_sites"),
    names(wide)
  )
  wide %>% mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x))))
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
pixy_folder_500kb <- paths$pixy_folder_500kb
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
figures_tables_dir <- file.path(
  Sys.getenv("HOME"),
  "Library/CloudStorage/GoogleDrive-jgharencar@ucdavis.edu",
  "My Drive/Costus/PopGen/Figures&Tables"
)
if (nzchar(Sys.getenv("PLOT_OUT_ROOT", ""))) {
  project_root <- Sys.getenv("PLOT_OUT_ROOT")
  final_figures_dir <- project_root
  figures_tables_dir <- project_root
  dir.create(file.path(project_root, "Dsuite"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(project_root, "pixy"), recursive = TRUE, showWarnings = FALSE)
}

# ---- Pixy: 10 kb points, 500 kb line (each file is already a pixy window) ----
pixy_wide <- read_pixy_wide(pixy_folder, "10 kb")
if (!all(c("hudson_fst_num", "hudson_fst_den") %in% names(pixy_wide))) {
  stop(
    "10 kb pixy output is missing hudson_fst_num and hudson_fst_den. ",
    "Re-run pixy with --fst_type hudson --fst_components."
  )
}
pixy_wide_500kb <- read_pixy_wide(pixy_folder_500kb, "500 kb")

# Genome-wide Fst is the ratio of sums. The mean of window ratios is printed
# only so it can be compared with that value; it is not used in the figures.
mean_fst_windows <- mean(pixy_wide$avg_hudson_fst, na.rm = TRUE)
mean_fst_ratio <- hudson_fst_ratio(pixy_wide$hudson_fst_num, pixy_wide$hudson_fst_den)
cat(
  "Genome-wide Hudson Fst (alle vs vill):\n",
  "  ratio of sums across 10 kb windows = ", sprintf("%.6f", mean_fst_ratio), "\n",
  "  mean of 10 kb ratios (not used in figures) = ", sprintf("%.6f", mean_fst_windows), "\n",
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

# ---- Genome-wide f_dM figure (positive values only, 99th % cutoff) ----
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
cutoff99thPercentile <- quantile(abs(LASIout_w50s25_all9$f_dM[LASIout_w50s25_all9$f_dM < 0]), 0.99)

LASIout_positive <- LASIout_w50s25_all9 %>% filter(f_dM > 0)

fdM_xlim <- range(LASIout_positive$windowStart_cum, na.rm = TRUE)

LASI_99thper_w50s25_f_dM <- ggplot(data = LASIout_positive, aes(x = windowStart_cum, y = f_dM, color = as_factor(chr))) +
  geom_point(size = 3) +
  geom_hline(yintercept = cutoff99thPercentile, linetype = "dashed", color = "red") +
  scale_color_manual(values = rep(c("black", "darkgray"), ceiling(length(unique(axis_set_dsuite$chr)) / 2))) +
  guides(color = "none") +
  geom_point(
    data = subset(LASIout_positive, f_dM > cutoff99thPercentile),
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
ggsave(file.path(project_root, "Dsuite", "LASI_99thper_w50s25_f_dM_pos_only_all_chroms.png"),
       plot = LASI_99thper_w50s25_f_dM,
       width = 24, height = 3.8, units = "in", dpi = 300)

# ---- Genome-wide f_dM figure (all finite values; negatives included) ----
LASIout_fd <- LASIout_w50s25_all9 %>% filter(is.finite(f_dM))
fdM_xlim_all <- range(LASIout_fd$windowStart_cum, na.rm = TRUE)
cutoff99_abs_fdM <- quantile(abs(LASIout_fd$f_dM), 0.99, na.rm = TRUE)

LASI_f_dM_incl_negative <- ggplot(data = LASIout_fd, aes(x = windowStart_cum, y = f_dM, color = as_factor(chr))) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_point(size = 3) +
  scale_color_manual(values = rep(c("black", "darkgray"), ceiling(length(unique(axis_set_dsuite$chr)) / 2))) +
  guides(color = "none") +
  geom_point(
    data = LASIout_fd %>% filter(abs(f_dM) > cutoff99_abs_fdM),
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

ggsave(file.path(project_root, "Dsuite", "LASI_99thper_w50s25_f_dM_incl_negative_all_chroms.png"),
       plot = LASI_f_dM_incl_negative,
       width = 24, height = 5.5, units = "in", dpi = 300)

# ---- QTL and trait colors (used by along-genome plots and histograms) ----
# Opposite of expected species direction from raw genotype means (AA = vill, BB = allenii).
# Expected: higher SD, GR, Nreabs, and Chlreabs in vill; higher toughness in allenii.
# GR chromosome 5 is in this set and is already removed by the CI filter below.
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

# ---- Along-genome: 500 kb pixy windows (not a mean of the 10 kb rows) ----
# WindowMid is the midpoint of pixy's 500 kb window. MeanY keeps the column
# name the line geoms already use; the value is avg_dxy or avg_hudson_fst.
# A 500 kb window is dropped when half or more of the 10 kb slots in that span
# fail the 10 kb filter (Fst: >=10 SNPs; Dxy: >=1000 sites) or have no finite
# value. Absent 10 kb rows count as missing. The y value is still pixy's own
# 500 kb window. line_run is assigned after dropping failed windows, and a new
# segment starts only when the next kept window is more than one 500 kb step
# away, so a dropped window is not bridged. Pi rows in the wide table share a
# midpoint and must not split adjacent windows into their own segments.
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
windowed_Dxy_avg <- line_500kb(
  pixy_wide_500kb, pixy_wide, "avg_dxy", "no_sites", min_dxy_sites_per_10kb
)
windowed_Fst_avg <- line_500kb(
  pixy_wide_500kb, pixy_wide, "avg_hudson_fst", "no_snps", min_snps_per_10kb
)
cat(
  "500 kb Dxy line: kept ", nrow(windowed_Dxy_avg),
  " windows; dropped ",
  sum(is.finite(pixy_wide_500kb$avg_dxy)) - nrow(windowed_Dxy_avg),
  " with >=50% of 10 kb slots below 1000 sites or missing.\n",
  "500 kb Fst line: kept ", nrow(windowed_Fst_avg),
  " windows; dropped ",
  sum(is.finite(pixy_wide_500kb$avg_hudson_fst)) - nrow(windowed_Fst_avg),
  " with >=50% of 10 kb slots below 10 SNPs or missing.\n",
  sep = ""
)

# Split by chromosome for plotting (list with chr1, ..., chr9; empty tibble if no passing bins)
windowed_Dxy_avg_list <- lapply(1:9, function(i) windowed_Dxy_avg %>% filter(chrN == i))
names(windowed_Dxy_avg_list) <- paste0("chr", 1:9)
windowed_Fst_avg_list <- lapply(1:9, function(i) windowed_Fst_avg %>% filter(chrN == i))
names(windowed_Fst_avg_list) <- paste0("chr", 1:9)

# ---- Dxy: plot with QTL ----
pixy_dxy_10kb <- pixy_wide %>%
  filter(is.finite(avg_dxy), is.finite(no_sites), no_sites >= min_dxy_sites_per_10kb)
cutoff99_dxy <- {
  x <- pixy_dxy_10kb$avg_dxy
  x <- pmax(x, 0)
  quantile(x[is.finite(x)], 0.99)
}

plot_chr_combined <- function(raw_df, win_df, cutoff99) {
  ggplot(raw_df, aes(x = window_pos_1)) +
    geom_point(
      data = raw_df %>% filter(pmax(avg_dxy, 0) <= cutoff99),
      aes(y = pmax(avg_dxy, 0)),
      size = 1,
      alpha = 0.4,
      color = "grey40"
    ) +
    geom_point(
      data = raw_df %>% filter(pmax(avg_dxy, 0) > cutoff99),
      aes(y = pmax(avg_dxy, 0)),
      size = 2,
      alpha = 0.9,
      color = "red"
    ) +
    geom_line(
      data = win_df,
      aes(x = WindowMid, y = MeanY, group = line_run),
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

  raw_df_i <- pixy_dxy_10kb %>% filter(chrN == i)
  win_df_i <- windowed_Dxy_avg_list[[paste0("chr", i)]]

  plot_chr_combined(
    raw_df = raw_df_i,
    win_df = win_df_i,
    cutoff99 = cutoff99_dxy
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
      aes(xintercept = peak_bp, color = trait, linetype = qtl_linetype),
      linewidth = 1
    ) +
    scale_fill_manual(values = trait_colors) +
    scale_color_manual(values = trait_colors) +
    scale_linetype_identity(guide = "none")
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
      aes(xintercept = peak_bp, color = trait, linetype = qtl_linetype),
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = trait_colors, name = NULL, drop = FALSE) +
    scale_color_manual(values = trait_colors, name = NULL, drop = FALSE) +
    scale_linetype_identity(guide = "none") +
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
pixy_fst_10kb <- pixy_wide %>%
  filter(is.finite(avg_hudson_fst), is.finite(no_snps), no_snps >= min_snps_per_10kb)
cutoff99_fst <- {
  x <- pixy_fst_10kb$avg_hudson_fst
  x <- pmax(x, 0)
  quantile(x[is.finite(x)], 0.99)
}

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
    labs(
      x = "Position (bp)",
      y = "Fst"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
    # Display window only; do not drop boundary points with scale_y limits
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

# Reuse QTL data and helper (already in env)
final_plots <- list()

for (i in seq_along(chrom_list)) {

  chr_num <- chrom_list[i]
  raw_df_i <- pixy_fst_10kb %>% filter(chrN == chr_num)
  win_df_i <- windowed_Fst_avg_list[[paste0("chr", chr_num)]]
  qtl_here <- all_QTL_CI %>%
    filter(chr == chr_num)

  plot_i <- plot_chr_fst_base(raw_df_i, cutoff99_fst)
  plot_i <- add_qtl_layers_fst(plot_i, qtl_here, trait_colors)
  plot_i <- add_fst_foreground(plot_i, raw_df_i, win_df_i, cutoff99_fst)
  plot_i <- plot_i + ggtitle(paste("Chr", chr_num))
  final_plots[[i]] <- plot_i
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

# Single-column stack: same 8 chromosomes, wider panels (cowplot avoids
# patchwork nesting issues when combining ncol=1 panels + legend)
fst_1col_w_cm <- 22
fst_1col_h_cm <- 48
fst_panels_1col <- cowplot::plot_grid(
  final_plots$Chrom1, final_plots$Chrom2, final_plots$Chrom3, final_plots$Chrom5,
  final_plots$Chrom6, final_plots$Chrom7, final_plots$Chrom8, final_plots$Chrom9,
  ncol = 1, align = "v"
)
Fst_LASI_sep_chrs_w.QTL_1col <- cowplot::plot_grid(
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

if (dir.exists(final_figures_dir)) {
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
  message("Wrote 2-col and 1-col Fst panels to ", final_figures_dir)
}

# ---- Genome-wide distribution: 1 cM window histograms with QTL averages ----
trait_labels <- c(
  SD = "Seed dormancy", GR = "Growth rate", tough = "Leaf toughness",
  Nreabs = "Nitrogen reabsorption\nduring drought", Chlreabs = "Chlorophyll reabsorption\nduring drought"
)
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
# 5th and 95th percentile cutoffs (red tail shading)
p05_fdM <- quantile(f_dM_summary$avg_f_dM_per_cM, 0.05)
p95_fdM <- quantile(f_dM_summary$avg_f_dM_per_cM, 0.95)

f_dM_windows_hist_w.QTL.CIs <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = p05_fdM, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.12) +
  annotate("rect", xmin = p95_fdM, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.12) +
  geom_histogram(data = f_dM_summary, aes(x = avg_f_dM_per_cM, y = after_stat(density)), bins = 50, fill = "grey80", color = "black") +
  geom_density(data = f_dM_summary, aes(x = avg_f_dM_per_cM, y = after_stat(density)), color = "darkblue", linewidth = 0.8) +
  geom_vline(
    data = qtl_f_dM_bp,
    aes(xintercept = avg_f_dM_QTL, color = trait, linetype = qtl_linetype),
    linewidth = 1,
    show.legend = TRUE,
    key_glyph = "point"
  ) +
  scale_color_manual(values = trait_colors, labels = trait_labels) +
  scale_linetype_identity(guide = "none") +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5, linetype = 0))) +
  labs(x = expression(italic(f)[d * M] ~ "(1 cM windows)"), y = "probability density", color = "QTL Trait:") +
  theme(legend.position = "left", legend.justification = "center")

# Dxy: 1 cM bins, QTL window averages (same seven intervals as the other panels)
dxy <- pixy_dxy_10kb %>%
  filter(chromosome != "Chrom4") %>%
  mutate(window_center = chr_position * 10e5) %>%
  select(chrN, chromosome, window_pos_1, window_pos_2, avg_dxy, window_center) %>%
  rowwise() %>%
  mutate(cM_window_start = bp_to_cm_funs[[as.character(chrN)]](window_pos_1)) %>%
  ungroup()
dxy_summary <- dxy %>%
  mutate(cM_bin = floor(cM_window_start)) %>%
  group_by(chrN, cM_bin) %>%
  summarise(avg_dxy_per_cM = mean(avg_dxy, na.rm = TRUE), n_windows = n(), .groups = "drop")
qtl_dxy_bp <- all_QTL_CI %>%
  rowwise() %>%
  mutate(
    avg_dxy_QTL = mean(
      dxy$avg_dxy[dxy$chrN == chr & dxy$window_pos_1 <= upper95CI_bp & dxy$window_pos_2 >= lower95CI_bp],
      na.rm = TRUE
    )
  ) %>%
  ungroup()
# 5th and 95th percentile cutoffs (red tail shading)
p05_dxy <- quantile(dxy_summary$avg_dxy_per_cM, 0.05)
p95_dxy <- quantile(dxy_summary$avg_dxy_per_cM, 0.95)

dxy_windows_hist_w.QTL.CIs <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = p05_dxy, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.12) +
  annotate("rect", xmin = p95_dxy, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.12) +
  geom_histogram(data = dxy_summary, aes(x = avg_dxy_per_cM, y = after_stat(density)), bins = 50, fill = "grey80", color = "black") +
  geom_density(data = dxy_summary, aes(x = avg_dxy_per_cM, y = after_stat(density)), color = "darkblue", linewidth = 0.8) +
  geom_vline(
    data = qtl_dxy_bp,
    aes(xintercept = avg_dxy_QTL, color = trait, linetype = qtl_linetype),
    linewidth = 1,
    show.legend = TRUE,
    key_glyph = "point"
  ) +
  scale_color_manual(values = trait_colors, labels = trait_labels) +
  scale_linetype_identity(guide = "none") +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 5, linetype = 0))) +
  labs(x = expression(italic(D)[xy] ~ "(1 cM windows)"), y = "probability density", color = "QTL Trait:") +
  theme(legend.position = "left", legend.justification = "center")

# Fst: 1 cM bins. Pixy has no cM window, so each bin is the ratio of summed
# Hudson components from the 10 kb rows that fall in it (not the mean of those Fst values).
fst <- pixy_fst_10kb %>%
  filter(chromosome != "Chrom4") %>%
  mutate(window_center = chr_position * 10e5) %>%
  select(
    chrN, chromosome, window_pos_1, window_pos_2,
    avg_hudson_fst, hudson_fst_num, hudson_fst_den, window_center
  ) %>%
  rowwise() %>%
  mutate(cM_window_start = bp_to_cm_funs[[as.character(chrN)]](window_pos_1)) %>%
  ungroup()
fst_summary <- fst %>%
  mutate(cM_bin = floor(cM_window_start)) %>%
  group_by(chrN, cM_bin) %>%
  summarise(
    avg_fst_per_cM = hudson_fst_ratio(hudson_fst_num, hudson_fst_den),
    n_windows = n(),
    .groups = "drop"
  )
fst_cm <- fst %>%
  arrange(chrN, cM_window_start) %>%
  group_by(chrN) %>%
  mutate(cM_window_end = lead(cM_window_start)) %>%
  ungroup()
qtl_fst_cm <- all_QTL_CI %>%
  rowwise() %>%
  mutate(
    avg_fst_QTL = {
      hit <- fst_cm$chrN == chr &
        fst_cm$cM_window_start <= upper95CI_cm &
        fst_cm$cM_window_end >= lower95CI_cm
      hit[is.na(hit)] <- FALSE
      hudson_fst_ratio(fst_cm$hudson_fst_num[hit], fst_cm$hudson_fst_den[hit])
    }
  ) %>%
  ungroup()
# Cache 1 cM Fst objects so the simplified plot chunk can be re-run alone.
# New name so a leftover Weir–Cockerham Fst_1cM_window_cache.rds is not reused.
saveRDS(
  list(fst_summary = fst_summary, qtl_fst_cm = qtl_fst_cm, project_root = project_root),
  file.path(project_root, "Fst_1cM_window_cache_hudson.rds")
)
# 5th and 95th percentile cutoffs (red tail shading)
p05_fst <- quantile(fst_summary$avg_fst_per_cM, 0.05)
p95_fst <- quantile(fst_summary$avg_fst_per_cM, 0.95)

qtl_fst_cm <- qtl_fst_cm %>%
  mutate(
    effect_class = factor(
      ifelse(opposite_direction, "antagonistic", "non-antagonistic"),
      levels = c("non-antagonistic", "antagonistic")
    )
  )

fst_windows_hist_w.QTL.CIs <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = p05_fst, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.12) +
  annotate("rect", xmin = p95_fst, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.12) +
  geom_histogram(data = fst_summary, aes(x = avg_fst_per_cM, y = after_stat(density)), binwidth = 0.01, fill = "grey80", color = "black") +
  geom_density(data = fst_summary, aes(x = avg_fst_per_cM, y = after_stat(density)), color = "darkblue", linewidth = 0.8) +
  geom_vline(
    data = qtl_fst_cm,
    aes(xintercept = avg_fst_QTL, color = trait, linetype = effect_class),
    linewidth = 1,
    show.legend = TRUE
  ) +
  scale_color_manual(values = trait_colors, labels = trait_labels) +
  scale_linetype_manual(
    values = c("non-antagonistic" = "solid", "antagonistic" = "dashed")
  ) +
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(linetype = "solid", linewidth = 1.2)
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(color = "black", linewidth = 0.9)
    )
  ) +
  labs(
    x = "Fst (1 cM windows)",
    y = "probability density",
    color = "QTL Trait:",
    linetype = NULL
  ) +
  theme(legend.position = "left", legend.justification = "center") +
  scale_x_continuous(expand = c(-.01, -.01))

# Solo Fst distribution figure (includes antagonistic vs non-antagonistic linetype legend)
ggsave(file.path(project_root, "Fst_genomewide_distribution_w_QTL_windows.png"), fst_windows_hist_w.QTL.CIs,
  device = "png", width = 15, height = 5.5, units = "cm")

# Multi-panel reuse keeps the dashed QTL lines and hides the second legend
fst_windows_hist_w.QTL.CIs <- fst_windows_hist_w.QTL.CIs + guides(linetype = "none")

# First-pass simplified 1 cM Fst: density only, dark-grey QTL CI means, no histogram/percentiles
# This block is re-runnable alone if Fst_1cM_window_cache_hudson.rds exists (written above on a full run).
if (!exists("fst_summary") || !exists("qtl_fst_cm")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  library(ggplot2)
  cache_candidates <- c(
    if (exists("project_root")) file.path(project_root, "Fst_1cM_window_cache_hudson.rds") else NULL,
    "Fst_1cM_window_cache_hudson.rds",
    file.path("..", "Fst_1cM_window_cache_hudson.rds"),
    file.path("..", "..", "Fst_1cM_window_cache_hudson.rds")
  )
  cache_path <- cache_candidates[file.exists(cache_candidates)][1]
  if (is.na(cache_path) || !nzchar(cache_path)) {
    stop(
      "fst_summary/qtl_fst_cm not in memory and Fst_1cM_window_cache_hudson.rds not found. ",
      "Run the script once through the 1 cM Fst summary (needs the Hudson pixy tables, lepmap, and QTL CIs)."
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
    aes(xintercept = avg_fst_QTL, linetype = qtl_linetype),
    color = "grey60",
    linewidth = 1
  ) +
  scale_linetype_identity(guide = "none") +
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
        fst_df$avg_hudson_fst[
          fst_df$chrN == chr &
            fst_df$window_pos_1 <= upper95CI_bp &
            fst_df$window_pos_2 >= lower95CI_bp
        ],
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

ratio_fst_in_qtl_bp <- function(qtl_rows, fst_df) {
  qtl_rows %>%
    rowwise() %>%
    mutate(
      avg_fst_QTL = {
        hit <- fst_df$chrN == chr &
          fst_df$window_pos_1 <= upper95CI_bp &
          fst_df$window_pos_2 >= lower95CI_bp
        hit[is.na(hit)] <- FALSE
        hudson_fst_ratio(fst_df$hudson_fst_num[hit], fst_df$hudson_fst_den[hit])
      }
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
      aes(xintercept = max_fst_QTL, color = trait, linetype = qtl_linetype),
      linewidth = 1,
      show.legend = TRUE,
      key_glyph = "point"
    ) +
    scale_color_manual(values = trait_colors, labels = trait_labels) +
    scale_linetype_identity(guide = "none") +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 5, linetype = 0))) +
    labs(x = xlab, y = "probability density", color = "QTL Trait:") +
    theme(legend.position = "left", legend.justification = "center") +
    scale_x_continuous(expand = c(-0.01, -0.01))
  ggsave(outfile, p, device = "png", width = 15, height = 5.5, units = "cm")
  invisible(p)
}

# 10 kb: native pixy windows (histogram = per-window Fst; QTL = max window in CI)
fst_10kb_summary <- fst %>% select(avg_fst = avg_hudson_fst)
qtl_fst_10kb_max <- max_fst_in_qtl_bp(all_QTL_CI, fst)
fst_10kb_hist <- make_fst_bp_hist(
  fst_10kb_summary,
  "avg_fst",
  qtl_fst_10kb_max,
  "Fst (10 kb windows)",
  file.path(project_root, "Fst_genomewide_distribution_10kb_w_QTL_max_windows.png")
)

# 10 kb: mean Fst across windows in each QTL CI (no genome-wide histogram)
qtl_fst_10kb_avg <- ratio_fst_in_qtl_bp(all_QTL_CI, fst) %>%
  mutate(qtl_label = paste0(trait, " (chr", chr, ")"))
fst_10kb_qtl_mean_plot <- ggplot(qtl_fst_10kb_avg, aes(x = avg_fst_QTL, y = reorder(qtl_label, avg_fst_QTL), color = trait)) +
  geom_segment(aes(x = 0, xend = avg_fst_QTL, y = qtl_label, yend = qtl_label, linetype = qtl_linetype), linewidth = 0.9, color = "grey70") +
  scale_linetype_identity(guide = "none") +
  geom_point(size = 4) +
  scale_color_manual(values = trait_colors, labels = trait_labels) +
  labs(
    x = "Fst (ratio of sums, 10 kb windows in QTL 95% CI)",
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

# Two-panel supplement: Dxy and f_dM (1 cM window distributions vs QTL window means)
dxy_fdM_hist_w.QTL.CIs <- dxy_windows_hist_w.QTL.CIs /
  (f_dM_windows_hist_w.QTL.CIs + guides(color = "none")) +
  plot_layout(guides = "collect") &
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

dxy_fdM_fig_w_cm <- 17
dxy_fdM_fig_h_cm <- 12
dxy_fdM_repo_png <- file.path(project_root, "Dxy_fdM_genomewide_distribution_w_QTL_2panel.png")
ggsave(
  dxy_fdM_repo_png,
  dxy_fdM_hist_w.QTL.CIs,
  device = "png",
  width = dxy_fdM_fig_w_cm,
  height = dxy_fdM_fig_h_cm,
  units = "cm",
  dpi = 200
)

sup_dxy_fdM_name <- "sup_Dxy_fdM_QTL_window_distributions.pdf"
if (dir.exists(figures_tables_dir)) {
  ggsave(
    file.path(figures_tables_dir, sup_dxy_fdM_name),
    dxy_fdM_hist_w.QTL.CIs,
    device = "pdf",
    width = dxy_fdM_fig_w_cm,
    height = dxy_fdM_fig_h_cm,
    units = "cm"
  )
  message("Wrote ", sup_dxy_fdM_name, " to ", figures_tables_dir)
}

if (dir.exists(final_figures_dir)) {
  ggsave(
    file.path(final_figures_dir, "sup_Dxy_fdM_QTL_window_distributions.pdf"),
    dxy_fdM_hist_w.QTL.CIs,
    device = "pdf",
    width = dxy_fdM_fig_w_cm,
    height = dxy_fdM_fig_h_cm,
    units = "cm"
  )
  ggsave(
    file.path(final_figures_dir, "sup_Dxy_fdM_QTL_window_distributions.png"),
    dxy_fdM_hist_w.QTL.CIs,
    device = "png",
    width = dxy_fdM_fig_w_cm,
    height = dxy_fdM_fig_h_cm,
    units = "cm",
    dpi = 200
  )
  message("Wrote sup_Dxy_fdM_QTL_window_distributions to ", final_figures_dir)
}

# ---- 10 kb distributions: density of windows inside each QTL 95% CI ----
# Same seven QTL intervals and the same antagonistic dashes as the 1 cM figures.
# The genome-wide histogram and dark-blue density are the 10 kb windows (chr 4
# excluded). Each colored curve is the kernel density of the 10 kb windows that
# overlap that QTL's 95% CI, not the CI-wide mean.
# Fst and Dxy are the pixy 10 kb values. f_dM has no native 10 kb window, so each
# value is the mean of positive Dsuite windows (same positive-only set as the
# 1 cM panel) whose midpoint falls in that 10 kb bin.
# coord_cartesian limits are a display window only.
# Fst keeps 10 kb windows with at least 10 SNPs. Dxy keeps 10 kb windows with at
# least 1000 sites (10% of the window). f_dM is unchanged: each 10 kb value is
# the mean of positive Dsuite windows whose midpoint falls in that bin.
# On the Fst panel, a point on each QTL curve marks the 10 kb window that contains
# that QTL's peak (peak_bp), if that window passes the SNP filter.
windows_in_qtl_cis <- function(windows, value_col, qtl = all_QTL_CI) {
  parts <- lapply(seq_len(nrow(qtl)), function(i) {
    q <- qtl[i, ]
    hit <- windows$chrN == q$chr &
      windows$window_pos_1 <= q$upper95CI_bp &
      windows$window_pos_2 >= q$lower95CI_bp
    hit[is.na(hit)] <- FALSE
    tibble(
      qtl_id = paste(q$trait, q$chr, sep = "_chr"),
      trait = q$trait,
      effect_class = if (isTRUE(q$opposite_direction)) "antagonistic" else "non-antagonistic",
      value = windows[[value_col]][hit]
    )
  })
  bind_rows(parts) %>%
    filter(is.finite(value)) %>%
    mutate(
      trait = factor(trait, levels = names(trait_colors)),
      effect_class = factor(effect_class, levels = c("non-antagonistic", "antagonistic"))
    ) %>%
    arrange(effect_class, qtl_id) %>%
    mutate(qtl_id = factor(qtl_id, levels = unique(qtl_id)))
}

dist_10kb_qtl_ci <- function(gw, xcol, qtl_long, xlab, xlim, binwidth) {
  x <- gw[[xcol]]
  p05 <- quantile(x, 0.05, na.rm = TRUE)
  p95 <- quantile(x, 0.95, na.rm = TRUE)
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = p05, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.12) +
    annotate("rect", xmin = p95, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.12) +
    geom_histogram(
      data = gw,
      aes(x = .data[[xcol]], y = after_stat(density)),
      binwidth = binwidth,
      fill = "grey80",
      color = "grey35",
      linewidth = 0.15
    ) +
    geom_density(
      data = gw,
      aes(x = .data[[xcol]]),
      color = "darkblue",
      linewidth = 0.8
    ) +
    geom_density(
      data = qtl_long,
      aes(x = value, color = trait, linetype = effect_class, group = qtl_id),
      linewidth = 0.95,
      key_glyph = "path"
    ) +
    scale_color_manual(values = trait_colors, labels = trait_labels, drop = FALSE) +
    scale_linetype_manual(
      values = c("non-antagonistic" = "solid", "antagonistic" = "dashed"),
      drop = FALSE
    ) +
    guides(
      color = guide_legend(
        order = 1,
        override.aes = list(linetype = "solid", linewidth = 1.3)
      ),
      linetype = guide_legend(
        order = 2,
        override.aes = list(color = "black", linewidth = 0.9)
      )
    ) +
    labs(x = xlab, y = "probability density", color = "QTL Trait:", linetype = NULL) +
    coord_cartesian(xlim = xlim) +
    theme(legend.position = "right", legend.justification = "center")
}

fst_10kb_gw <- pixy_wide %>%
  mutate(
    no_snps = suppressWarnings(as.numeric(no_snps)),
    hudson_fst_num = suppressWarnings(as.numeric(hudson_fst_num)),
    hudson_fst_den = suppressWarnings(as.numeric(hudson_fst_den))
  ) %>%
  filter(chrN != 4, is.finite(avg_hudson_fst), is.finite(no_snps), no_snps >= min_snps_per_10kb) %>%
  transmute(
    chrN, window_pos_1, window_pos_2, avg_hudson_fst,
    hudson_fst_num, hudson_fst_den, no_snps
  )
dxy_10kb_gw <- pixy_wide %>%
  mutate(no_sites = suppressWarnings(as.numeric(no_sites))) %>%
  filter(chrN != 4, is.finite(avg_dxy), is.finite(no_sites), no_sites >= min_dxy_sites_per_10kb) %>%
  transmute(chrN, window_pos_1, window_pos_2, avg_dxy, no_sites)
message(
  "10 kb filters: ", nrow(fst_10kb_gw), " Fst windows with >=10 SNPs; ",
  nrow(dxy_10kb_gw), " Dxy windows with >=1000 sites."
)
fdm_10kb_gw <- pos_f_dM_clean %>%
  filter(chrN != 4, is.finite(f_dM)) %>%
  mutate(
    mid = (windowStart + windowEnd) / 2,
    window_pos_1 = floor((mid - 1) / 10000) * 10000 + 1,
    window_pos_2 = window_pos_1 + 9999
  ) %>%
  group_by(chrN, window_pos_1, window_pos_2) %>%
  summarise(avg_f_dM = mean(f_dM), .groups = "drop")

qtl_fst_10kb_ci <- windows_in_qtl_cis(fst_10kb_gw, "avg_hudson_fst")
qtl_dxy_10kb_ci <- windows_in_qtl_cis(dxy_10kb_gw, "avg_dxy")
qtl_fdm_10kb_ci <- windows_in_qtl_cis(fdm_10kb_gw, "avg_f_dM")
message(
  "10 kb Fst windows inside QTL 95% CIs:\n",
  paste(
    format(qtl_fst_10kb_ci %>% count(qtl_id, effect_class), n = Inf),
    collapse = "\n"
  )
)

# Height of a kernel density at x, using the same nrd0 bandwidth as geom_density
# so the peak marker sits on the drawn curve.
density_at <- function(x, at) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || !is.finite(at)) return(NA_real_)
  d <- stats::density(x, bw = stats::bw.nrd0(x), n = 512)
  stats::approx(d$x, d$y, xout = at, rule = 2)$y
}

qtl_fst_peak_pts <- all_QTL_CI %>%
  mutate(
    qtl_id = paste(trait, chr, sep = "_chr"),
    trait = factor(trait, levels = names(trait_colors))
  ) %>%
  rowwise() %>%
  mutate(
    peak_fst = {
      hit <- fst_10kb_gw$chrN == chr &
        fst_10kb_gw$window_pos_1 <= peak_bp &
        fst_10kb_gw$window_pos_2 >= peak_bp
      hit[is.na(hit)] <- FALSE
      vals <- fst_10kb_gw$avg_hudson_fst[hit]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) NA_real_ else vals[[1]]
    },
    peak_density = density_at(
      qtl_fst_10kb_ci$value[as.character(qtl_fst_10kb_ci$qtl_id) == qtl_id],
      peak_fst
    )
  ) %>%
  ungroup()
if (any(!is.finite(qtl_fst_peak_pts$peak_fst))) {
  missing_peaks <- qtl_fst_peak_pts %>% filter(!is.finite(peak_fst))
  warning(
    "No 10 kb Fst window covers the QTL peak for: ",
    paste(missing_peaks$qtl_id, collapse = ", ")
  )
}
message(
  "Fst of the 10 kb window at each QTL peak:\n",
  paste(
    format(
      qtl_fst_peak_pts %>% select(qtl_id, peak_bp, peak_fst, peak_density),
      n = Inf
    ),
    collapse = "\n"
  )
)

fst_10kb_ci_dist <- dist_10kb_qtl_ci(
  fst_10kb_gw, "avg_hudson_fst", qtl_fst_10kb_ci,
  "Fst (10 kb windows, >=10 SNPs)", c(-0.25, 0.65), 0.01
) +
  geom_point(
    data = qtl_fst_peak_pts %>% filter(is.finite(peak_fst), is.finite(peak_density)),
    aes(x = peak_fst, y = peak_density, fill = trait),
    shape = 21,
    color = "black",
    size = 2.8,
    stroke = 0.45,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = trait_colors, guide = "none")
dxy_10kb_ci_dist <- dist_10kb_qtl_ci(
  dxy_10kb_gw, "avg_dxy", qtl_dxy_10kb_ci,
  expression(italic(D)[xy] ~ "(10 kb windows, >= 1000 sites)"), c(0, 0.04), 0.001
)
fdm_10kb_ci_dist <- dist_10kb_qtl_ci(
  fdm_10kb_gw, "avg_f_dM", qtl_fdm_10kb_ci,
  expression(italic(f)[d * M] ~ "(10 kb windows)"), c(0, 0.55), 0.01
)

fst_dxy_fdm_10kb_ci_3panel <- fst_10kb_ci_dist /
  (dxy_10kb_ci_dist + guides(color = "none", linetype = "none")) /
  (fdm_10kb_ci_dist + guides(color = "none", linetype = "none")) +
  plot_layout(guides = "collect")
ggsave(
  file.path(project_root, "Fst_Dxy_fdM_genomewide_distribution_10kb_w_QTL_CI_3panel.png"),
  fst_dxy_fdm_10kb_ci_3panel,
  device = "png",
  width = 18.5,
  height = 18,
  units = "cm",
  dpi = 200
)
ggsave(
  file.path(project_root, "Fst_genomewide_distribution_10kb_w_QTL_CI.png"),
  fst_10kb_ci_dist,
  device = "png",
  width = 16,
  height = 6.5,
  units = "cm",
  dpi = 200
)
dxy_fdm_10kb_ci_2panel <- dxy_10kb_ci_dist /
  (fdm_10kb_ci_dist + guides(color = "none", linetype = "none")) +
  plot_layout(guides = "collect")
ggsave(
  file.path(project_root, "Dxy_fdM_genomewide_distribution_10kb_w_QTL_CI_2panel.png"),
  dxy_fdm_10kb_ci_2panel,
  device = "png",
  width = 18.5,
  height = 12,
  units = "cm",
  dpi = 200
)

# ---- 0.1 cM distributions, built from the filtered 10 kb windows ----
# Fst bins are the ratio of Hudson components across 10 kb windows with >=10 SNPs.
# Dxy bins are the mean of 10 kb windows with >=1000 sites. f_dM bins are the mean
# of positive Dsuite windows (no 10 kb site/SNP filter). A bin overlaps a QTL CI
# when its interval overlaps the CI on this script's linkage map (bp through
# bp_to_cm_funs, not the raw cM columns in the QTL table). The Fst point is the
# 0.1 cM bin that contains peak_bp on that map.
add_window_cm <- function(windows) {
  if (!exists("bp_to_cm_funs")) {
    stop("bp_to_cm_funs is missing. The 0.1 cM windows need the Lep-MAP3 map built above.")
  }
  windows %>%
    filter(as.character(chrN) %in% names(bp_to_cm_funs)) %>%
    group_by(chrN) %>%
    group_modify(function(d, key) {
      fun <- bp_to_cm_funs[[as.character(key$chrN)]]
      d %>% mutate(cM_window_start = fun(window_pos_1))
    }) %>%
    ungroup() %>%
    filter(is.finite(cM_window_start))
}

# CI bounds in the same cM coordinates as add_window_cm().
qtl_cm_on_map <- function(qtl = all_QTL_CI) {
  qtl %>%
    rowwise() %>%
    mutate(
      cm_a = bp_to_cm_funs[[as.character(chr)]](lower95CI_bp),
      cm_b = bp_to_cm_funs[[as.character(chr)]](upper95CI_bp),
      cm_peak = bp_to_cm_funs[[as.character(chr)]](peak_bp),
      cm_lo = min(cm_a, cm_b, na.rm = TRUE),
      cm_hi = max(cm_a, cm_b, na.rm = TRUE)
    ) %>%
    ungroup()
}

windows_in_qtl_cm <- function(bins, value_col, qtl = qtl_cm_on_map()) {
  parts <- lapply(seq_len(nrow(qtl)), function(i) {
    q <- qtl[i, ]
    hit <- bins$chrN == q$chr &
      bins$bin_lo < q$cm_hi &
      bins$bin_hi > q$cm_lo
    hit[is.na(hit)] <- FALSE
    tibble(
      qtl_id = paste(q$trait, q$chr, sep = "_chr"),
      trait = q$trait,
      effect_class = if (isTRUE(q$opposite_direction)) "antagonistic" else "non-antagonistic",
      value = bins[[value_col]][hit]
    )
  })
  bind_rows(parts) %>%
    filter(is.finite(value)) %>%
    mutate(
      trait = factor(trait, levels = names(trait_colors)),
      effect_class = factor(effect_class, levels = c("non-antagonistic", "antagonistic"))
    ) %>%
    arrange(effect_class, qtl_id) %>%
    mutate(qtl_id = factor(qtl_id, levels = unique(qtl_id)))
}

display_xlim <- function(x, probs = c(0.005, 0.995), pad = 0.06) {
  qs <- quantile(x, probs, na.rm = TRUE)
  span <- diff(as.numeric(qs))
  c(qs[[1]] - pad * span, qs[[2]] + pad * span)
}

fst_01cm <- add_window_cm(fst_10kb_gw) %>%
  mutate(cM_bin = floor(cM_window_start / 0.1)) %>%
  group_by(chrN, cM_bin) %>%
  summarise(
    avg_fst = hudson_fst_ratio(hudson_fst_num, hudson_fst_den),
    .groups = "drop"
  ) %>%
  filter(is.finite(avg_fst)) %>%
  mutate(bin_lo = cM_bin * 0.1, bin_hi = bin_lo + 0.1)
dxy_01cm <- add_window_cm(dxy_10kb_gw) %>%
  mutate(cM_bin = floor(cM_window_start / 0.1)) %>%
  group_by(chrN, cM_bin) %>%
  summarise(avg_dxy = mean(avg_dxy, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(avg_dxy)) %>%
  mutate(bin_lo = cM_bin * 0.1, bin_hi = bin_lo + 0.1)
fdm_01cm <- add_window_cm(
  pos_f_dM_clean %>%
    filter(chrN != 4, is.finite(f_dM)) %>%
    transmute(chrN, window_pos_1 = windowStart, f_dM)
) %>%
  mutate(cM_bin = floor(cM_window_start / 0.1)) %>%
  group_by(chrN, cM_bin) %>%
  summarise(avg_f_dM = mean(f_dM, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(avg_f_dM)) %>%
  mutate(bin_lo = cM_bin * 0.1, bin_hi = bin_lo + 0.1)

qtl_fst_01cm_ci <- windows_in_qtl_cm(fst_01cm, "avg_fst")
qtl_dxy_01cm_ci <- windows_in_qtl_cm(dxy_01cm, "avg_dxy")
qtl_fdm_01cm_ci <- windows_in_qtl_cm(fdm_01cm, "avg_f_dM")
message(
  "0.1 cM bins: Fst ", nrow(fst_01cm),
  "; Dxy ", nrow(dxy_01cm),
  "; f_dM ", nrow(fdm_01cm), "."
)
message(
  "0.1 cM bins inside QTL 95% CIs (Fst):\n",
  paste(format(qtl_fst_01cm_ci %>% count(qtl_id, effect_class), n = Inf), collapse = "\n")
)

qtl_fst_peak_01 <- qtl_cm_on_map() %>%
  mutate(
    qtl_id = paste(trait, chr, sep = "_chr"),
    trait = factor(trait, levels = names(trait_colors)),
    chr = as.integer(chr),
    cM_bin = floor(cm_peak / 0.1)
  ) %>%
  left_join(
    fst_01cm %>% mutate(chrN = as.integer(chrN)) %>% select(chrN, cM_bin, peak_fst = avg_fst),
    by = c("chr" = "chrN", "cM_bin")
  ) %>%
  rowwise() %>%
  mutate(
    peak_density = density_at(
      qtl_fst_01cm_ci$value[as.character(qtl_fst_01cm_ci$qtl_id) == qtl_id],
      peak_fst
    )
  ) %>%
  ungroup()
if (any(!is.finite(qtl_fst_peak_01$peak_fst))) {
  missing_peaks <- qtl_fst_peak_01 %>% filter(!is.finite(peak_fst))
  warning(
    "No 0.1 cM Fst bin covers the QTL peak for: ",
    paste(missing_peaks$qtl_id, collapse = ", ")
  )
}

fst_01cm_dist <- dist_10kb_qtl_ci(
  fst_01cm, "avg_fst", qtl_fst_01cm_ci,
  "Fst (0.1 cM windows)", display_xlim(fst_01cm$avg_fst), 0.005
) +
  geom_point(
    data = qtl_fst_peak_01 %>% filter(is.finite(peak_fst), is.finite(peak_density)),
    aes(x = peak_fst, y = peak_density, fill = trait),
    shape = 21,
    color = "black",
    size = 2.8,
    stroke = 0.45,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = trait_colors, guide = "none")
dxy_01cm_dist <- dist_10kb_qtl_ci(
  dxy_01cm, "avg_dxy", qtl_dxy_01cm_ci,
  expression(italic(D)[xy] ~ "(0.1 cM windows)"), display_xlim(dxy_01cm$avg_dxy), 0.0005
)
fdm_01cm_dist <- dist_10kb_qtl_ci(
  fdm_01cm, "avg_f_dM", qtl_fdm_01cm_ci,
  expression(italic(f)[d * M] ~ "(0.1 cM windows)"), display_xlim(fdm_01cm$avg_f_dM), 0.005
)

fst_dxy_fdm_01cm_3panel <- fst_01cm_dist /
  (dxy_01cm_dist + guides(color = "none", linetype = "none")) /
  (fdm_01cm_dist + guides(color = "none", linetype = "none")) +
  plot_layout(guides = "collect")
ggsave(
  file.path(project_root, "Fst_Dxy_fdM_genomewide_distribution_0.1cM_w_QTL_CI_3panel.png"),
  fst_dxy_fdm_01cm_3panel,
  device = "png",
  width = 18.5,
  height = 18,
  units = "cm",
  dpi = 200
)
ggsave(
  file.path(project_root, "Fst_genomewide_distribution_0.1cM_w_QTL_CI.png"),
  fst_01cm_dist,
  device = "png",
  width = 16,
  height = 6.5,
  units = "cm",
  dpi = 200
)
dxy_fdm_01cm_2panel <- dxy_01cm_dist /
  (fdm_01cm_dist + guides(color = "none", linetype = "none")) +
  plot_layout(guides = "collect")
ggsave(
  file.path(project_root, "Dxy_fdM_genomewide_distribution_0.1cM_w_QTL_CI_2panel.png"),
  dxy_fdm_01cm_2panel,
  device = "png",
  width = 18.5,
  height = 12,
  units = "cm",
  dpi = 200
)
