# =============================================================================
# Compare within-chromosome introgression landscapes: full vs noNearHybrid 1 km
#
# Dsuite w50s25 windows are variant-count windows, so coordinates shift when the
# sample set changes. Compare on a shared physical grid: assign each window's
# midpoint to 100 kb bins (and also summarize at 500 kb for tracks).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

out <- "04_introgression_vs_QTL/noNearHybrid_1km"
full_dir <- "Dsuite/Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120"
sub_dir <- "Dsuite/dsuite_noNearHybrid_1km_w50s25"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

bin_fine <- 100000L
bin_coarse <- 500000L

load_dsuite <- function(d, lab) {
  files <- list.files(d, pattern = "\\.txt$", full.names = TRUE)
  bind_rows(lapply(files, read.delim, sep = "\t", stringsAsFactors = FALSE)) %>%
    mutate(
      chrN = as.integer(str_extract(chr, "\\d+")),
      run = lab,
      across(any_of(c("D", "f_d", "f_dM", "d_f", "windowStart", "windowEnd")),
             ~ suppressWarnings(as.numeric(.x))),
      mid = (windowStart + windowEnd) / 2
    ) %>%
    filter(is.finite(mid), is.finite(chrN))
}

message("Loading...")
full <- load_dsuite(full_dir, "full")
sub <- load_dsuite(sub_dir, "subset")
message("Windows: full=", nrow(full), " subset=", nrow(sub))

bin_means <- function(df, bin_bp) {
  df %>%
    mutate(bin_start = floor(mid / bin_bp) * bin_bp) %>%
    group_by(chrN, bin_start, run) %>%
    summarise(
      n_windows = n(),
      mid = mean(mid),
      across(c(D, f_d, f_dM, d_f), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    )
}

join_bins <- function(bin_bp) {
  bf <- bin_means(full, bin_bp)
  bs <- bin_means(sub, bin_bp)
  inner_join(
    bf %>% select(-run),
    bs %>% select(-run),
    by = c("chrN", "bin_start"),
    suffix = c("_full", "_sub")
  ) %>%
    mutate(
      bin_mid = bin_start + bin_bp / 2,
      bin_bp = bin_bp,
      f_dM_delta = f_dM_sub - f_dM_full,
      d_f_delta = d_f_sub - d_f_full,
      D_delta = D_sub - D_full
    )
}

fine <- join_bins(bin_fine)
coarse <- join_bins(bin_coarse)
message("Shared 100 kb bins: ", nrow(fine), "; shared 500 kb bins: ", nrow(coarse))

stats <- c("D", "f_d", "f_dM", "d_f")

summarise_concordance <- function(joined, label) {
  map_dfr(stats, function(st) {
    xf <- paste0(st, "_full")
    xs <- paste0(st, "_sub")
    chr <- joined %>%
      group_by(chrN) %>%
      summarise(
        grid = label,
        statistic = st,
        n_bins = n(),
        pearson_r = cor(.data[[xf]], .data[[xs]], use = "complete.obs"),
        spearman_rho = cor(.data[[xf]], .data[[xs]], method = "spearman", use = "complete.obs"),
        mean_full = mean(.data[[xf]], na.rm = TRUE),
        mean_sub = mean(.data[[xs]], na.rm = TRUE),
        mean_delta = mean(.data[[xs]] - .data[[xf]], na.rm = TRUE),
        mae = mean(abs(.data[[xs]] - .data[[xf]]), na.rm = TRUE),
        rmse = sqrt(mean((.data[[xs]] - .data[[xf]])^2, na.rm = TRUE)),
        .groups = "drop"
      )
    gen <- tibble(
      chrN = NA_integer_,
      grid = label,
      statistic = st,
      n_bins = nrow(joined),
      pearson_r = cor(joined[[xf]], joined[[xs]], use = "complete.obs"),
      spearman_rho = cor(joined[[xf]], joined[[xs]], method = "spearman", use = "complete.obs"),
      mean_full = mean(joined[[xf]], na.rm = TRUE),
      mean_sub = mean(joined[[xs]], na.rm = TRUE),
      mean_delta = mean(joined[[xs]] - joined[[xf]], na.rm = TRUE),
      mae = mean(abs(joined[[xs]] - joined[[xf]]), na.rm = TRUE),
      rmse = sqrt(mean((joined[[xs]] - joined[[xf]])^2, na.rm = TRUE))
    )
    bind_rows(chr, gen)
  })
}

conc <- bind_rows(
  summarise_concordance(fine, "100kb"),
  summarise_concordance(coarse, "500kb")
)
write_csv(conc, file.path(out, "landscape_concordance_by_chr.csv"))

# ---- Outlier concordance on 100 kb f_dM means ----
# Use each run's genome-wide 99th percentile of positive 100 kb means among bins
# with mean f_dM > 0, or stick to same rule as plots: 99th of |neg| on raw windows
# projected via bin means exceeding each run's raw-window threshold.
thr_full <- as.numeric(quantile(abs(full$f_dM[is.finite(full$f_dM) & full$f_dM < 0]), 0.99))
thr_sub <- as.numeric(quantile(abs(sub$f_dM[is.finite(sub$f_dM) & sub$f_dM < 0]), 0.99))

fine <- fine %>%
  mutate(
    out_full = is.finite(f_dM_full) & f_dM_full > thr_full,
    out_sub = is.finite(f_dM_sub) & f_dM_sub > thr_sub,
    out_either = out_full | out_sub,
    out_both = out_full & out_sub
  )

outlier_chr <- fine %>%
  group_by(chrN) %>%
  summarise(
    n_bins = n(),
    n_out_full = sum(out_full),
    n_out_sub = sum(out_sub),
    n_out_both = sum(out_both),
    n_out_full_only = sum(out_full & !out_sub),
    n_out_sub_only = sum(out_sub & !out_full),
    jaccard = ifelse(sum(out_either) == 0, NA_real_, sum(out_both) / sum(out_either)),
    .groups = "drop"
  )
outlier_genome <- tibble(
  chrN = NA_integer_,
  n_bins = nrow(fine),
  n_out_full = sum(fine$out_full),
  n_out_sub = sum(fine$out_sub),
  n_out_both = sum(fine$out_both),
  n_out_full_only = sum(fine$out_full & !fine$out_sub),
  n_out_sub_only = sum(fine$out_sub & !fine$out_full),
  jaccard = sum(fine$out_both) / max(1, sum(fine$out_either))
)
outlier_tbl <- bind_rows(outlier_chr, outlier_genome) %>%
  mutate(thr_full = thr_full, thr_sub = thr_sub, grid = "100kb")
write_csv(outlier_tbl, file.path(out, "landscape_fdM_outlier_concordance.csv"))

# Also: top 1% of |f_dM| within each run on the shared 100 kb grid (rank-based)
fine <- fine %>%
  mutate(
    top1_full = is.finite(f_dM_full) & f_dM_full >= quantile(f_dM_full, 0.99, na.rm = TRUE),
    top1_sub = is.finite(f_dM_sub) & f_dM_sub >= quantile(f_dM_sub, 0.99, na.rm = TRUE)
  )
top1_jaccard <- sum(fine$top1_full & fine$top1_sub) / sum(fine$top1_full | fine$top1_sub)

theme_set(theme_bw(base_size = 11))
run_cols <- c(Full = "#444444", `No near hybrid` = "#D55E00")

# ---- Plots ----

# 1) Spearman by chromosome (100 kb grid)
rho_plot <- conc %>%
  filter(grid == "100kb", !is.na(chrN)) %>%
  mutate(statistic = factor(statistic, levels = stats)) %>%
  ggplot(aes(x = factor(chrN), y = spearman_rho, fill = statistic)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(yintercept = c(0.5, 0.8), linetype = "dotted", color = "grey50") +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Chromosome", y = "Spearman rho (100 kb bin means)", fill = NULL,
    title = "Within-chromosome concordance of introgression landscapes",
    subtitle = sprintf(
      "Shared 100 kb bins (n=%d). Window coords differ between runs; compare via physical bins.",
      nrow(fine)
    )
  ) +
  theme(legend.position = "bottom")
ggsave(file.path(out, "landscape_spearman_by_chr.png"), rho_plot, width = 8.5, height = 4.8, dpi = 150)

# 2) Genome hex-like density via geom_bin2d (no hexbin package needed)
dens <- ggplot(fine, aes(f_dM_full, f_dM_sub)) +
  geom_bin2d(bins = 50) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10", name = "n bins") +
  labs(
    x = expression(f[dM] ~ "(full, 100 kb mean)"),
    y = expression(f[dM] ~ "(no near hybrid, 100 kb mean)"),
    title = "Shared 100 kb bins: f_dM full vs subset",
    subtitle = {
      g <- conc %>% filter(grid == "100kb", is.na(chrN), statistic == "f_dM")
      sprintf("Pearson r=%.3f; Spearman rho=%.3f; mean delta=%.4f",
              g$pearson_r, g$spearman_rho, g$mean_delta)
    }
  )
ggsave(file.path(out, "landscape_fdM_bin2d_genome.png"), dens, width = 6.5, height = 5.5, dpi = 150)

# 3) Per-chr scatter of 100 kb f_dM
sc <- ggplot(fine, aes(f_dM_full, f_dM_sub)) +
  geom_point(alpha = 0.25, size = 0.7, color = "grey25") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~chrN, ncol = 3, labeller = labeller(chrN = \(x) paste("chr", x))) +
  labs(
    x = expression(f[dM] ~ "(full)"), y = expression(f[dM] ~ "(no near hybrid)"),
    title = "100 kb bin means of f_dM by chromosome"
  )
ggsave(file.path(out, "landscape_scatter_fdM_by_chr.png"), sc, width = 9, height = 9, dpi = 150)

# 4) Along-chromosome 500 kb tracks
plot_chr_track <- function(chr) {
  d <- coarse %>% filter(chrN == chr)
  long <- d %>%
    select(bin_mid, f_dM_full, f_dM_sub) %>%
    pivot_longer(-bin_mid, names_to = "run", values_to = "f_dM") %>%
    mutate(run = recode(run, f_dM_full = "Full", f_dM_sub = "No near hybrid"))
  p_top <- ggplot(long, aes(bin_mid / 1e6, f_dM, color = run)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 0.55) +
    scale_color_manual(values = run_cols) +
    labs(x = NULL, y = expression(f[dM]), color = NULL, title = paste("chr", chr)) +
    theme(legend.position = "none", plot.title = element_text(size = 10))
  p_bot <- ggplot(d, aes(bin_mid / 1e6, f_dM_delta)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey40") +
    geom_line(color = "#0072B2", linewidth = 0.5) +
    labs(x = "Mb", y = expression(Delta))
  p_top / p_bot + plot_layout(heights = c(2.2, 1))
}
track_all <- wrap_plots(lapply(1:9, plot_chr_track), ncol = 3) +
  plot_annotation(
    title = "f_dM landscape stability along chromosomes (500 kb means)",
    subtitle = "Top: full vs no-near-hybrid; bottom: subset − full",
    caption = "Grey = full; orange = no near hybrid (1 km)"
  )
ggsave(file.path(out, "landscape_fdM_tracks_by_chr.png"), track_all, width = 12, height = 14, dpi = 140)

# 5) Outlier Jaccard
jac <- ggplot(outlier_chr, aes(factor(chrN), jaccard)) +
  geom_col(fill = "#8250C4", width = 0.7) +
  geom_hline(yintercept = outlier_genome$jaccard, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Chromosome", y = "Jaccard of high-f_dM 100 kb bins",
    title = "Where do extreme introgression bins agree?",
    subtitle = sprintf(
      "Bin mean > own-run window 99|neg| threshold (full=%.3f, subset=%.3f). Genome Jaccard=%.2f (dashed). Top-1%% rank Jaccard=%.2f",
      thr_full, thr_sub, outlier_genome$jaccard, top1_jaccard
    )
  )
ggsave(file.path(out, "landscape_fdM_outlier_jaccard.png"), jac, width = 8, height = 4.5, dpi = 150)

# 6) Delta distribution by chromosome
delta_plot <- ggplot(fine, aes(x = factor(chrN), y = f_dM_delta)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_violin(fill = "#56B4E9", alpha = 0.5, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, alpha = 0.8) +
  labs(
    x = "Chromosome", y = expression(f[dM] * " subset − full (100 kb)"),
    title = "Distribution of local f_dM change by chromosome"
  )
ggsave(file.path(out, "landscape_fdM_delta_violin.png"), delta_plot, width = 8, height = 4.5, dpi = 150)

# ---- Summary text ----
g100 <- conc %>% filter(grid == "100kb", is.na(chrN))
g500 <- conc %>% filter(grid == "500kb", is.na(chrN))
chr_fdm <- conc %>% filter(grid == "100kb", !is.na(chrN), statistic == "f_dM") %>% arrange(chrN)

summary_txt <- c(
  "Introgression landscape comparison: full vs noNearHybrid 1 km",
  "Method: mean Dsuite stats in shared physical bins (window midpoints).",
  sprintf("Raw windows: full=%d, subset=%d", nrow(full), nrow(sub)),
  sprintf("Shared bins: 100 kb n=%d; 500 kb n=%d", nrow(fine), nrow(coarse)),
  "",
  "Genome-wide concordance (100 kb):",
  sprintf("  f_dM: r=%.3f, rho=%.3f, mean_delta=%.4f, MAE=%.4f",
          g100$pearson_r[g100$statistic == "f_dM"],
          g100$spearman_rho[g100$statistic == "f_dM"],
          g100$mean_delta[g100$statistic == "f_dM"],
          g100$mae[g100$statistic == "f_dM"]),
  sprintf("  d_f:  r=%.3f, rho=%.3f, mean_delta=%.4f, MAE=%.4f",
          g100$pearson_r[g100$statistic == "d_f"],
          g100$spearman_rho[g100$statistic == "d_f"],
          g100$mean_delta[g100$statistic == "d_f"],
          g100$mae[g100$statistic == "d_f"]),
  "",
  "Genome-wide concordance (500 kb smoothed):",
  sprintf("  f_dM: r=%.3f, rho=%.3f",
          g500$pearson_r[g500$statistic == "f_dM"],
          g500$spearman_rho[g500$statistic == "f_dM"]),
  sprintf("  d_f:  r=%.3f, rho=%.3f",
          g500$pearson_r[g500$statistic == "d_f"],
          g500$spearman_rho[g500$statistic == "d_f"]),
  "",
  sprintf("High-f_dM 100 kb outliers (bin mean > own-run 99|neg| thr): full=%d, subset=%d, both=%d, Jaccard=%.3f",
          outlier_genome$n_out_full, outlier_genome$n_out_sub,
          outlier_genome$n_out_both, outlier_genome$jaccard),
  sprintf("Top 1%% f_dM bins (rank within each run) Jaccard=%.3f", top1_jaccard),
  "",
  "Per-chromosome Spearman rho for f_dM (100 kb):",
  paste(sprintf("  chr%d: %.3f (n=%d)", chr_fdm$chrN, chr_fdm$spearman_rho, chr_fdm$n_bins),
        collapse = "\n"),
  "",
  "Interpretation sketch:",
  "  - Broad along-chromosome patterns are moderately to strongly conserved (rho often ~0.6-0.9).",
  "  - Absolute levels drop in the subset (negative mean delta).",
  "  - Extreme outlier bins are less shared than the overall landscape (lower Jaccard),",
  "    so peaks move/weaken more than the background pattern."
)
writeLines(summary_txt, file.path(out, "landscape_comparison_summary.txt"))
message(paste(summary_txt, collapse = "\n"))
message("\nWrote figures + CSVs to ", out)
