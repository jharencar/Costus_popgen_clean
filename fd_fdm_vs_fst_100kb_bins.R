# =============================================================================
# fd_fdm_vs_fst_100kb_bins.R
#
# Phase 1: Shared 100 kb bins — pixy W&C F_ST + Dsuite f_d, f_dM, d_f (w50 step 25).
# Dsuite: stride-2 (default) or greedy; overlap weights w_i = O_i / len_i.
#
# value_mode:
#   all — all finite bins
#   positive — F_ST > 0 and introgression stat > 0
#   positive_no_outliers — positive then drop |studentized residual| > 3 from lm(y ~ F_ST)
#
# Spearman rho; cor.test p is heuristic under spatial autocorrelation.
#
# CLI: Rscript fd_fdm_vs_fst_100kb_bins.R [stride2|greedy] [all|positive|positive_no_outliers]
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# ---- Paths ----
project_root <- "/Users/admin-jgharenc/Documents/Github/Costus_popgen_clean"
dsuite_dir <- file.path(
  project_root,
  "Dsuite",
  "Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120"
)
pixy_folder <- file.path(project_root, "pixy", "output_251124_LASIref")

window_size_bp <- 100000L
min_frac_valid_fst <- 0.25
student_thresh <- 3

argv <- commandArgs(trailingOnly = TRUE)
thin_method <- if (length(argv) >= 1L && argv[1] %in% c("greedy", "stride2")) {
  argv[1]
} else {
  "stride2"
}
value_mode <- if (
  length(argv) >= 2L &&
    argv[2] %in% c("all", "positive", "positive_no_outliers")
) {
  argv[2]
} else {
  "positive"
}

# ---- Pixy ----
pixy_to_long <- function(pixy_files) {
  pixy_df <- list()
  for (i in seq_along(pixy_files)) {
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    if (stat_file_type == "pi") {
      df <- read_delim(pixy_files[i], delim = "\t", show_col_types = FALSE)
      df <- df %>%
        gather(
          -pop, -window_pos_1, -window_pos_2, -chromosome,
          key = "statistic", value = "value"
        ) %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      pixy_df[[i]] <- df
    } else {
      df <- read_delim(pixy_files[i], delim = "\t", show_col_types = FALSE)
      df <- df %>%
        gather(
          -pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
          key = "statistic", value = "value"
        )
      pixy_df[[i]] <- df
    }
  }
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
}

assert_pixy_intervals_nonoverlapping <- function(pixy_fst_df) {
  bad <- pixy_fst_df %>%
    arrange(chromosome, window_pos_1) %>%
    group_by(chromosome) %>%
    mutate(
      prev_end = lag(window_pos_2),
      overlaps_prev = !is.na(prev_end) & window_pos_1 <= prev_end
    ) %>%
    filter(overlaps_prev)
  if (nrow(bad) > 0) {
    warning(
      "Pixy F_ST: found ",
      nrow(bad),
      " row(s) overlapping previous interval on same chromosome."
    )
  }
}

overlap_len_inclusive <- function(ws, we, bs, be) {
  as.integer(pmax(0L, pmin(we, be) - pmax(ws, bs) + 1L))
}

greedy_nonoverlap_pick <- function(df) {
  if (nrow(df) == 0) {
    return(df)
  }
  ord <- order(df$windowEnd, df$windowStart)
  x <- df[ord, , drop = FALSE]
  keep <- rep(FALSE, nrow(x))
  last_end <- -Inf
  for (i in seq_len(nrow(x))) {
    if (x$windowStart[i] > last_end) {
      keep[i] <- TRUE
      last_end <- x$windowEnd[i]
    }
  }
  x[keep, , drop = FALSE]
}

remove_lm_outliers <- function(df, y_col, x_col = "avg_wc_fst_wmean", thresh = student_thresh) {
  if (nrow(df) < 4L) {
    warning("Too few rows (", nrow(df), ") for studentized residuals; no outliers removed.")
    return(df)
  }
  f <- as.formula(paste(y_col, "~", x_col))
  fit <- tryCatch(lm(f, data = df), error = function(e) NULL)
  if (is.null(fit)) {
    warning("lm failed; no outliers removed.")
    return(df)
  }
  rs <- abs(stats::rstudent(fit))
  ok <- !is.na(rs) & rs <= thresh
  df[ok, , drop = FALSE]
}

pixy_files <- list.files(pixy_folder, full.names = TRUE)
pixy_df <- pixy_to_long(pixy_files)
pixy_wide <- pivot_wider(pixy_df, names_from = statistic, values_from = value)

pixy_fst <- pixy_wide %>%
  filter(pop1 == "alle", pop2 == "vill") %>%
  filter(!is.na(window_pos_1), !is.na(window_pos_2))

assert_pixy_intervals_nonoverlapping(pixy_fst)

pixy_fst_frag <- pixy_fst %>%
  mutate(
    chrN = as.numeric(str_extract(chromosome, "\\d+")),
    tile_len = as.integer(window_pos_2 - window_pos_1 + 1L),
    k_lo = as.integer(floor(window_pos_1 / window_size_bp)),
    k_hi = as.integer(floor(window_pos_2 / window_size_bp))
  ) %>%
  rowwise() %>%
  mutate(bin_k = list(seq(k_lo, k_hi))) %>%
  ungroup() %>%
  unnest(bin_k) %>%
  mutate(
    bin_start_bp = as.integer(bin_k) * window_size_bp + 1L,
    bin_end_bp = as.integer(bin_k + 1L) * window_size_bp,
    overlap_bp_pixy = overlap_len_inclusive(
      as.integer(window_pos_1),
      as.integer(window_pos_2),
      bin_start_bp,
      bin_end_bp
    ),
    w_pixy = ifelse(
      is.finite(no_snps) & no_snps > 0 & is.finite(avg_wc_fst),
      no_snps * overlap_bp_pixy / tile_len,
      NA_real_
    )
  )

pixy_bins <- pixy_fst_frag %>%
  group_by(chrN, chromosome, bin_k) %>%
  summarise(
    n_pixy_rows = n(),
    n_valid_fst = sum(is.finite(avg_wc_fst)),
    sum_no_snps_pixy = sum(no_snps[is.finite(no_snps) & no_snps > 0], na.rm = TRUE),
    avg_wc_fst_wmean = {
      ok <- is.finite(w_pixy) & w_pixy > 0 & is.finite(avg_wc_fst)
      if (!any(ok)) {
        NA_real_
      } else {
        weighted.mean(avg_wc_fst[ok], w = w_pixy[ok])
      }
    },
    .groups = "drop"
  ) %>%
  mutate(frac_valid = n_valid_fst / pmax(n_pixy_rows, 1L)) %>%
  filter(frac_valid >= min_frac_valid_fst) %>%
  mutate(
    bin_start_bp = as.integer(bin_k) * window_size_bp + 1L,
    bin_end_bp = as.integer(bin_k + 1L) * window_size_bp
  )

# ---- Dsuite ----
file_list <- list.files(path = dsuite_dir, pattern = "*.txt", full.names = TRUE)
if (length(file_list) == 0) {
  stop("No Dsuite *.txt files in ", dsuite_dir)
}

LASIout <- map_dfr(file_list, ~ read.delim(.x, sep = "\t", stringsAsFactors = FALSE))

LASIout_std <- LASIout %>%
  mutate(
    chrN = as.numeric(str_extract(chr, "\\d+")),
    windowStart = as.integer(windowStart),
    windowEnd = as.integer(windowEnd)
  )

expand_dsuite_to_bins <- function(df_rows) {
  df_rows %>%
    mutate(
      len_bp = windowEnd - windowStart + 1L,
      k_lo = as.integer(floor(windowStart / window_size_bp)),
      k_hi = as.integer(floor(windowEnd / window_size_bp))
    ) %>%
    rowwise() %>%
    mutate(bin_k = list(seq(k_lo, k_hi))) %>%
    ungroup() %>%
    unnest(bin_k) %>%
    mutate(
      bin_start_bp = bin_k * window_size_bp + 1L,
      bin_end_bp = (bin_k + 1L) * window_size_bp,
      overlap_bp = overlap_len_inclusive(
        windowStart, windowEnd, bin_start_bp, bin_end_bp
      ),
      w_frac = overlap_bp / len_bp
    ) %>%
    filter(
      overlap_bp > 0L,
      is.finite(f_d),
      is.finite(f_dM),
      is.finite(d_f),
      w_frac > 0
    )
}

if (thin_method == "stride2") {
  dsuite_thin <- LASIout_std %>%
    group_by(chr) %>%
    arrange(windowStart, windowEnd, .by_group = TRUE) %>%
    mutate(slide_i = row_number()) %>%
    ungroup() %>%
    filter(slide_i %% 2L == 1L)

  dsuite_expand <- expand_dsuite_to_bins(dsuite_thin)
  thinning_method_label <- "stride2_chr_order"
} else if (thin_method == "greedy") {
  dsuite_expand_raw <- expand_dsuite_to_bins(LASIout_std)
  dsuite_expand <- dsuite_expand_raw %>%
    group_by(chrN, bin_k) %>%
    group_modify(function(.x, .y) {
      bk <- as.integer(.y$bin_k[1])
      bs <- bk * window_size_bp + 1L
      be <- (bk + 1L) * window_size_bp
      g <- greedy_nonoverlap_pick(.x)
      if (nrow(g) == 0) {
        return(tibble())
      }
      g %>%
        mutate(
          bin_start_bp = bs,
          bin_end_bp = be,
          len_bp = windowEnd - windowStart + 1L,
          overlap_bp = overlap_len_inclusive(windowStart, windowEnd, bs, be),
          w_frac = overlap_bp / len_bp
        ) %>%
        filter(overlap_bp > 0L, w_frac > 0)
    }) %>%
    ungroup()
  thinning_method_label <- "greedy_bp_nonoverlap_per_bin"
} else {
  stop("thin_method must be \"stride2\" or \"greedy\", got: ", thin_method)
}

dsuite_bins <- dsuite_expand %>%
  group_by(chrN, bin_k) %>%
  summarise(
    n_dsuite_windows_used = n(),
    sum_w_frac = sum(w_frac),
    f_d_wmean = sum(w_frac * f_d) / sum(w_frac),
    f_dM_wmean = sum(w_frac * f_dM) / sum(w_frac),
    d_f_wmean = sum(w_frac * d_f) / sum(w_frac),
    dsuite_span_kb_mean = mean((windowEnd - windowStart + 1L) / 1000),
    dsuite_span_kb_median = median((windowEnd - windowStart + 1L) / 1000),
    dsuite_span_kb_min = min((windowEnd - windowStart + 1L) / 1000),
    dsuite_span_kb_max = max((windowEnd - windowStart + 1L) / 1000),
    sum_O_i_used = sum(overlap_bp),
    .groups = "drop"
  )

dsuite_expand_raw <- LASIout_std %>%
  mutate(
    len_bp = windowEnd - windowStart + 1L,
    k_lo = as.integer(floor(windowStart / window_size_bp)),
    k_hi = as.integer(floor(windowEnd / window_size_bp))
  ) %>%
  rowwise() %>%
  mutate(bin_k = list(seq(k_lo, k_hi))) %>%
  ungroup() %>%
  unnest(bin_k) %>%
  mutate(
    bin_start_bp = bin_k * window_size_bp + 1L,
    bin_end_bp = (bin_k + 1L) * window_size_bp,
    overlap_bp = overlap_len_inclusive(
      windowStart, windowEnd, bin_start_bp, bin_end_bp
    )
  ) %>%
  filter(overlap_bp > 0L)

n_raw_per_bin <- dsuite_expand_raw %>%
  count(chrN, bin_k, name = "n_dsuite_windows_raw")

dsuite_bins <- dsuite_bins %>%
  left_join(n_raw_per_bin, by = c("chrN", "bin_k"))

paired_union <- dsuite_expand %>%
  group_by(chrN, bin_k) %>%
  summarise(
    frac_bin_bp_covered = min(sum(overlap_bp) / window_size_bp, 1),
    .groups = "drop"
  )

paired <- pixy_bins %>%
  inner_join(dsuite_bins, by = c("chrN", "bin_k")) %>%
  left_join(paired_union, by = c("chrN", "bin_k")) %>%
  filter(
    is.finite(avg_wc_fst_wmean),
    is.finite(f_d_wmean),
    is.finite(f_dM_wmean),
    is.finite(d_f_wmean)
  ) %>%
  mutate(
    thinning_method = thinning_method_label,
    thin_method_param = thin_method
  )

# ---- Filter by value_mode (per stat column) ----
filter_positive <- function(df, ycol) {
  df %>%
    filter(.data[["avg_wc_fst_wmean"]] > 0, .data[[ycol]] > 0)
}

paired_fd <- if (value_mode == "all") {
  paired
} else if (value_mode == "positive") {
  filter_positive(paired, "f_d_wmean")
} else {
  remove_lm_outliers(filter_positive(paired, "f_d_wmean"), "f_d_wmean")
}

paired_fdm <- if (value_mode == "all") {
  paired
} else if (value_mode == "positive") {
  filter_positive(paired, "f_dM_wmean")
} else {
  remove_lm_outliers(filter_positive(paired, "f_dM_wmean"), "f_dM_wmean")
}

paired_df <- if (value_mode == "all") {
  paired
} else if (value_mode == "positive") {
  filter_positive(paired, "d_f_wmean")
} else {
  remove_lm_outliers(filter_positive(paired, "d_f_wmean"), "d_f_wmean")
}

spearman_block <- function(x, y) {
  rho <- cor(x, y, method = "spearman", use = "complete.obs")
  ct <- cor.test(x, y, method = "spearman", exact = FALSE)
  list(rho = rho, p = ct$p.value, n = sum(complete.cases(x, y)))
}

s_fd <- spearman_block(paired_fd$avg_wc_fst_wmean, paired_fd$f_d_wmean)
s_fdm <- spearman_block(paired_fdm$avg_wc_fst_wmean, paired_fdm$f_dM_wmean)
s_df <- spearman_block(paired_df$avg_wc_fst_wmean, paired_df$d_f_wmean)

message(
  "Spearman rho (100 kb bins; mode=", value_mode, "; thin=", thin_method, "):\n",
  "  f_d  vs F_ST: rho = ", signif(s_fd$rho, 5), "  (n = ", s_fd$n, ")\n",
  "  cor.test approximate p = ", signif(s_fd$p, 5), "\n",
  "  f_dM vs F_ST: rho = ", signif(s_fdm$rho, 5), "  (n = ", s_fdm$n, ")\n",
  "  cor.test approximate p = ", signif(s_fdm$p, 5), "\n",
  "  d_f  vs F_ST: rho = ", signif(s_df$rho, 5), "  (n = ", s_df$n, ")\n",
  "  cor.test approximate p = ", signif(s_df$p, 5), "\n",
  "Note: p-values assume independent bins; genome ordering induces positive autocorrelation."
)

mode_suffix <- value_mode

readr::write_tsv(paired_fd, file.path(project_root, paste0("fd_fst_100kb_bins_", mode_suffix, ".tsv")))
readr::write_tsv(paired_fdm, file.path(project_root, paste0("fdm_fst_100kb_bins_", mode_suffix, ".tsv")))
readr::write_tsv(paired_df, file.path(project_root, paste0("df_fst_100kb_bins_", mode_suffix, ".tsv")))
message("Wrote TSVs for mode: ", mode_suffix)

thin_lab <- if (thin_method == "stride2") {
  "Dsuite stride-2 + overlap weights"
} else {
  "Dsuite greedy non-overlap + overlap weights"
}

mode_subtitle <- switch(
  value_mode,
  all = "",
  positive = "Positive F_ST and stat only",
  positive_no_outliers = paste0(
    "Positive only; outliers removed (|studentized residual| > ",
    student_thresh,
    " from lm)"
  )
)

theme_scatter <- theme_bw(base_size = 12) +
  theme(aspect.ratio = 1)

plot_one <- function(dat, y, ylab_expr, title_expr, rho, n) {
  lbl <- paste0("Spearman rho = ", signif(rho, 4), " (n = ", n, ")")
  p <- ggplot(dat, aes(x = avg_wc_fst_wmean, y = .data[[y]])) +
    geom_point(alpha = 0.35, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.9) +
    labs(
      x = expression("F"[ST] ~ "(SNP-weighted mean, pixy 10 kb to 100 kb)"),
      y = ylab_expr,
      title = title_expr,
      subtitle = lbl
    ) +
    theme_scatter
  if (nzchar(mode_subtitle)) {
    p <- p + labs(caption = mode_subtitle)
  }
  p
}

p_fd <- plot_one(
  paired_fd,
  "f_d_wmean",
  bquote(italic(f)[d] ~ "(" * .(thin_lab) * ")"),
  expression(italic(f)[d] ~ " vs " ~ "F"[ST] ~ " (100 kb bins)"),
  s_fd$rho,
  s_fd$n
)
p_fdm <- plot_one(
  paired_fdm,
  "f_dM_wmean",
  bquote(italic(f)[d * M] ~ "(" * .(thin_lab) * ")"),
  expression(italic(f)[d * M] ~ " vs " ~ "F"[ST] ~ " (100 kb bins)"),
  s_fdm$rho,
  s_fdm$n
)
p_df <- plot_one(
  paired_df,
  "d_f_wmean",
  bquote(italic(d)[f] ~ "(" * .(thin_lab) * ")"),
  expression(italic(d)[f] ~ " vs " ~ "F"[ST] ~ " (100 kb bins)"),
  s_df$rho,
  s_df$n
)

ggsave(
  file.path(project_root, paste0("fd_vs_fst_100kb_bins_", mode_suffix, ".png")),
  p_fd,
  width = 6,
  height = 6.4,
  dpi = 150,
  bg = "white"
)
ggsave(
  file.path(project_root, paste0("fdm_vs_fst_100kb_bins_", mode_suffix, ".png")),
  p_fdm,
  width = 6,
  height = 6.4,
  dpi = 150,
  bg = "white"
)
ggsave(
  file.path(project_root, paste0("df_vs_fst_100kb_bins_", mode_suffix, ".png")),
  p_df,
  width = 6,
  height = 6.4,
  dpi = 150,
  bg = "white"
)
ggsave(
  file.path(project_root, paste0("fd_vs_fst_100kb_bins_", mode_suffix, ".pdf")),
  p_fd,
  width = 6,
  height = 6.4,
  bg = "white"
)
ggsave(
  file.path(project_root, paste0("fdm_vs_fst_100kb_bins_", mode_suffix, ".pdf")),
  p_fdm,
  width = 6,
  height = 6.4,
  bg = "white"
)
ggsave(
  file.path(project_root, paste0("df_vs_fst_100kb_bins_", mode_suffix, ".pdf")),
  p_df,
  width = 6,
  height = 6.4,
  bg = "white"
)

message("Wrote PNG/PDF plots for mode: ", mode_suffix)
