# =============================================================================
# Circular bp shift of QTL 95% CIs vs outlier peaks.
#
# Two one-sided tests, same null. Each CI keeps its observed base-pair length.
# One shared shift per chromosome, wrapping at the chromosome end:
#   1. High Hudson Fst peaks overlap QTL CIs more often than chance (enrichment).
#   2. High f_dM peaks overlap QTL CIs less often than chance (depletion).
#   3. High Dxy peaks overlap QTL CIs more often than chance (enrichment).
#
# High f_dM is the positive tail. Dsuite files are brac_alle_vill, so P1 =
# bracteatus, P2 = allenii, P3 = villosissimus. Positive f_dM is excess sharing
# between allenii and villosissimus (Malinsky et al. 2021). The negative tail
# is not used.
#
# QTL sets, in order: all retained CIs (huge Nreabs chr2 and GR chr5 dropped),
# then non-antagonistic only. Chromosome 4 is excluded. The four p-values share
# one placement null and are not independent.
#
# Two outlier definitions share the same bp shifts. N_PERM defaults to 10000.
#   strict (unsuffixed files): Hudson Fst and Dxy >= genome-wide 99th percentile;
#     high f_dM >= the 99th percentile of |f_dM| among negative windows.
#     That equals -quantile(negative f_dM, 0.01), the sign-flipped extreme tail.
#   relaxed (_p95 files): Hudson Fst and Dxy >= genome-wide 95th percentile;
#     high f_dM >= the 95th percentile of |f_dM| among negative windows
#     = -quantile(negative f_dM, 0.05).
# Outputs:
#   circular_shift_outlier_qtl_overlap.rds
#   circular_shift_outlier_qtl_overlap_summary.csv
#   circular_shift_outlier_qtl_overlap_per_ci.csv
#   the same three names with _p95 before the extension,
#   and _p95_noSD (95th thresholds, seed-dormancy CIs removed from both QTL sets).
# =============================================================================

suppressPackageStartupMessages(library(tidyverse))

for (p in c(
  file.path("scripts", "project_paths.R"),
  file.path("..", "scripts", "project_paths.R")
)) {
  if (file.exists(p)) {
    source(p)
    break
  }
}
if (!exists("project_paths")) {
  stop("Could not source scripts/project_paths.R. Run from the repo root or from 04_introgression_vs_QTL/.")
}

paths <- project_paths()
project_root <- paths$project_root
pixy_folder <- paths$pixy_folder
dsuite_dir <- paths$dsuite_dir
qtl_ci_path <- paths$qtl_ci_path

if (is.na(qtl_ci_path) || !nzchar(qtl_ci_path) || !file.exists(qtl_ci_path)) {
  alt <- file.path(project_root, "QTL", "all_QTL_CIs_bp_cm.251130.csv")
  if (file.exists(alt)) qtl_ci_path <- alt
}

stop_if_missing <- function(path, label) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    stop(label, " not found: ", path)
  }
}
stop_if_missing(pixy_folder, "Hudson 10 kb pixy folder")
stop_if_missing(dsuite_dir, "Dsuite window folder")
stop_if_missing(qtl_ci_path, "QTL CI table")

n_perm <- as.integer(Sys.getenv("N_PERM", unset = "10000"))
if (!is.finite(n_perm) || n_perm < 1L) stop("N_PERM must be a positive integer.")
set.seed(251130L)

out_dir <- if (nzchar(Sys.getenv("PLOT_OUT_ROOT", ""))) {
  Sys.getenv("PLOT_OUT_ROOT")
} else {
  file.path(project_root, "04_introgression_vs_QTL")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- QTL sets ----
opposite_direction_qtl <- tibble::tribble(
  ~trait, ~chr,
  "SD", 6L,
  "GR", 5L,
  "tough", 3L,
  "tough", 8L
)
all_QTL_CI <- read.csv(qtl_ci_path, header = TRUE, stringsAsFactors = FALSE) %>%
  filter(!(chr == 2 & trait == "Nreabs")) %>%
  filter(!(chr == 5 & trait == "GR")) %>%
  filter(as.integer(chr) != 4L) %>%
  mutate(
    opposite_direction = paste(trait, as.integer(chr)) %in%
      paste(opposite_direction_qtl$trait, opposite_direction_qtl$chr)
  )

qtl_sets <- list(
  all_retained = all_QTL_CI,
  non_antagonistic = all_QTL_CI %>% filter(!opposite_direction)
)

# ---- Outlier peaks ----
# Genome-wide thresholds use every chromosome, including chromosome 4, so they
# match the published cutoffs. Peaks on chromosome 4 are then dropped.
fst_files <- list.files(pixy_folder, pattern = "_fst\\.txt$", full.names = TRUE)
if (length(fst_files) == 0) stop("No *_fst.txt files in ", pixy_folder)
fst_win <- bind_rows(lapply(fst_files, read.delim, sep = "\t", stringsAsFactors = FALSE))
if (!"avg_hudson_fst" %in% names(fst_win)) {
  stop("Pixy Fst file is missing avg_hudson_fst. Use the Hudson run, not output_251124_LASIref.")
}
fst_win <- fst_win %>%
  mutate(
    chrN = as.integer(str_extract(chromosome, "\\d+")),
    avg_hudson_fst = suppressWarnings(as.numeric(avg_hudson_fst)),
    window_pos_1 = as.numeric(window_pos_1),
    window_pos_2 = as.numeric(window_pos_2)
  )
fst_for_peaks <- fst_win %>%
  filter(chrN != 4L, is.finite(avg_hudson_fst)) %>%
  arrange(chrN, window_pos_1)

dxy_files <- list.files(pixy_folder, pattern = "_dxy\\.txt$", full.names = TRUE)
if (length(dxy_files) == 0) stop("No *_dxy.txt files in ", pixy_folder)
dxy_win <- bind_rows(lapply(dxy_files, read.delim, sep = "\t", stringsAsFactors = FALSE))
if (!"avg_dxy" %in% names(dxy_win)) stop("Pixy Dxy file is missing avg_dxy.")
dxy_win <- dxy_win %>%
  mutate(
    chrN = as.integer(str_extract(chromosome, "\\d+")),
    avg_dxy = suppressWarnings(as.numeric(avg_dxy)),
    window_pos_1 = as.numeric(window_pos_1),
    window_pos_2 = as.numeric(window_pos_2)
  )
dxy_for_peaks <- dxy_win %>%
  filter(chrN != 4L, is.finite(avg_dxy)) %>%
  arrange(chrN, window_pos_1)

merge_consecutive_runs <- function(df) {
  peaks <- vector("list", length(unique(df$chrN)))
  k <- 0L
  for (ch in unique(df$chrN)) {
    sub <- df %>% filter(chrN == ch)
    runs <- rle(sub$is_out)
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    for (i in which(runs$values)) {
      rows <- seq.int(starts[i], ends[i])
      k <- k + 1L
      peaks[[k]] <- tibble(
        chr = ch,
        start = min(sub$window_pos_1[rows]),
        end = max(sub$window_pos_2[rows])
      )
    }
  }
  if (k == 0L) return(tibble(chr = integer(), start = numeric(), end = numeric()))
  bind_rows(peaks[seq_len(k)])
}

dsuite_files <- list.files(dsuite_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(dsuite_files) == 0) stop("No Dsuite txt files in ", dsuite_dir)
fdm_win <- bind_rows(lapply(dsuite_files, read.delim, sep = "\t", stringsAsFactors = FALSE)) %>%
  mutate(
    chrN = as.integer(str_extract(chr, "\\d+")),
    f_dM = suppressWarnings(as.numeric(f_dM)),
    windowStart = as.numeric(windowStart),
    windowEnd = as.numeric(windowEnd)
  )
neg_fdm <- fdm_win$f_dM[is.finite(fdm_win$f_dM) & fdm_win$f_dM < 0]
fdm_for_peaks <- fdm_win %>%
  filter(chrN != 4L, is.finite(f_dM))

merge_overlapping_outliers <- function(df) {
  peaks <- list()
  for (ch in unique(df$chrN)) {
    sub <- df %>%
      filter(chrN == ch, is_out) %>%
      arrange(windowStart, windowEnd)
    if (nrow(sub) == 0) next
    cs <- sub$windowStart[1]
    ce <- sub$windowEnd[1]
    if (nrow(sub) > 1) {
      for (i in seq.int(2, nrow(sub))) {
        if (sub$windowStart[i] <= ce) {
          ce <- max(ce, sub$windowEnd[i])
        } else {
          peaks[[length(peaks) + 1L]] <- tibble(chr = ch, start = cs, end = ce)
          cs <- sub$windowStart[i]
          ce <- sub$windowEnd[i]
        }
      }
    }
    peaks[[length(peaks) + 1L]] <- tibble(chr = ch, start = cs, end = ce)
  }
  if (length(peaks) == 0) return(tibble(chr = integer(), start = numeric(), end = numeric()))
  bind_rows(peaks)
}

# Positive f_dM cutoff from the negative tail. quantile(|f_dM|, p) among
# negative windows equals -quantile(f_dM, 1 - p) on those same windows.
fdm_positive_cutoff <- function(tail_prob) {
  as.numeric(quantile(abs(neg_fdm), tail_prob, na.rm = TRUE))
}

# Chromosome length is the furthest base in the windows or a retained CI.
# The circle is [0, L]. A coordinate that lands on an exact multiple of L is
# the chromosome end, so an interval that finishes there is not dropped.
chr_end <- bind_rows(
  fst_win %>% filter(chrN != 4L) %>% transmute(chr = chrN, end = window_pos_2),
  dxy_win %>% filter(chrN != 4L) %>% transmute(chr = chrN, end = window_pos_2),
  fdm_win %>% filter(chrN != 4L) %>% transmute(chr = chrN, end = windowEnd),
  all_QTL_CI %>% transmute(chr = as.integer(chr), end = upper95CI_bp)
) %>%
  group_by(chr) %>%
  summarise(L = max(end, na.rm = TRUE), .groups = "drop")
L_bp <- setNames(chr_end$L, as.character(chr_end$chr))

too_long <- all_QTL_CI %>%
  mutate(width = upper95CI_bp - lower95CI_bp, L = L_bp[as.character(chr)]) %>%
  filter(width >= L)
if (nrow(too_long) > 0) {
  stop("A retained CI is at least as long as its chromosome in bp.")
}

# ---- Shift and overlap ----
shift_ci_bp <- function(lo, hi, shift, chr_num) {
  L <- L_bp[[as.character(chr_num)]]
  if (!is.finite(L) || L <= 0) stop("Chromosome ", chr_num, " has no bp length.")
  width <- hi - lo
  if (width >= L) return(list(c(0, L)))
  a <- (lo + shift) %% L
  b <- (hi + shift) %% L
  if (b == 0 && a > 0) {
    list(c(a, L))
  } else if (a <= b) {
    list(c(a, b))
  } else {
    list(c(a, L), c(0, b))
  }
}

peak_hits_segs <- function(peak_start, peak_end, segs) {
  if (length(peak_start) == 0 || length(segs) == 0) {
    return(rep(FALSE, length(peak_start)))
  }
  hit <- rep(FALSE, length(peak_start))
  for (seg in segs) {
    hit <- hit | (peak_start <= seg[2] & peak_end >= seg[1])
  }
  hit
}

merge_bp_segs <- function(segs) {
  if (length(segs) == 0) return(segs)
  m <- do.call(rbind, segs)
  m <- m[order(m[, 1]), , drop = FALSE]
  out <- list()
  cs <- m[1, 1]
  ce <- m[1, 2]
  if (nrow(m) > 1) {
    for (i in seq.int(2, nrow(m))) {
      if (m[i, 1] <= ce) {
        ce <- max(ce, m[i, 2])
      } else {
        out[[length(out) + 1L]] <- c(cs, ce)
        cs <- m[i, 1]
        ce <- m[i, 2]
      }
    }
  }
  out[[length(out) + 1L]] <- c(cs, ce)
  out
}

overlap_bp <- function(peak_start, peak_end, segs) {
  segs <- merge_bp_segs(segs)
  if (length(segs) == 0 || length(peak_start) == 0) return(0)
  total <- 0
  for (i in seq_along(peak_start)) {
    for (seg in segs) {
      lo <- max(peak_start[i], seg[1])
      hi <- min(peak_end[i], seg[2])
      if (hi > lo) total <- total + (hi - lo)
    }
  }
  total
}

ci_original_segs <- function(qtl_row) {
  list(c(qtl_row$lower95CI_bp, qtl_row$upper95CI_bp))
}

# One bp shift per chromosome, shared by both QTL sets and both statistics.
# Each draw is Uniform(0, L) and the CI keeps its observed bp length.
chr_with_qtl <- as.character(sort(unique(all_QTL_CI$chr)))
shifts <- lapply(chr_with_qtl, function(ch) runif(n_perm, min = 0, max = L_bp[[ch]]))
names(shifts) <- chr_with_qtl

peaks_by_chr <- function(peak_df) split(peak_df, peak_df$chr)

score_placement <- function(qtl, peak_df, seg_fun) {
  by_chr <- peaks_by_chr(peak_df)
  n_peak <- 0L
  n_ci <- 0L
  bp <- 0
  per_ci <- logical(nrow(qtl))
  all_segs <- list()
  for (i in seq_len(nrow(qtl))) {
    ch <- as.character(qtl$chr[i])
    segs <- seg_fun(qtl[i, , drop = FALSE])
    all_segs[[ch]] <- c(all_segs[[ch]], segs)
    pk <- by_chr[[ch]]
    if (is.null(pk) || length(segs) == 0) {
      per_ci[i] <- FALSE
      next
    }
    per_ci[i] <- any(peak_hits_segs(pk$start, pk$end, segs))
  }
  for (ch in names(all_segs)) {
    pk <- by_chr[[ch]]
    if (is.null(pk)) next
    hit <- peak_hits_segs(pk$start, pk$end, all_segs[[ch]])
    n_peak <- n_peak + sum(hit)
    bp <- bp + overlap_bp(pk$start[hit], pk$end[hit], all_segs[[ch]])
  }
  n_ci <- sum(per_ci)
  list(n_peaks = n_peak, n_cis = n_ci, overlap_bp = bp, per_ci = per_ci)
}

observed_scores <- function(qtl, peak_df) {
  score_placement(qtl, peak_df, ci_original_segs)
}

null_scores <- function(qtl, peak_df) {
  n_peaks <- integer(n_perm)
  n_cis <- integer(n_perm)
  overlap_bp_null <- numeric(n_perm)
  qtl_chr <- split(seq_len(nrow(qtl)), qtl$chr)
  peak_chr <- peaks_by_chr(peak_df)
  for (p in seq_len(n_perm)) {
    n_peak <- 0L
    n_ci <- 0L
    bp <- 0
    for (ch in names(qtl_chr)) {
      pk <- peak_chr[[ch]]
      rows <- qtl_chr[[ch]]
      segs <- list()
      ci_hit <- 0L
      s <- shifts[[ch]][p]
      for (i in rows) {
        seg_i <- shift_ci_bp(qtl$lower95CI_bp[i], qtl$upper95CI_bp[i], s, qtl$chr[i])
        segs <- c(segs, seg_i)
        if (!is.null(pk) && any(peak_hits_segs(pk$start, pk$end, seg_i))) {
          ci_hit <- ci_hit + 1L
        }
      }
      if (!is.null(pk)) {
        hit <- peak_hits_segs(pk$start, pk$end, segs)
        n_peak <- n_peak + sum(hit)
        bp <- bp + overlap_bp(pk$start[hit], pk$end[hit], segs)
      }
      n_ci <- n_ci + ci_hit
    }
    n_peaks[p] <- n_peak
    n_cis[p] <- n_ci
    overlap_bp_null[p] <- bp
  }
  list(n_peaks = n_peaks, n_cis = n_cis, overlap_bp = overlap_bp_null)
}

one_sided_p <- function(null, obs, alternative) {
  if (alternative == "enrichment") {
    (1 + sum(null >= obs, na.rm = TRUE)) / (1 + length(null))
  } else if (alternative == "depletion") {
    (1 + sum(null <= obs, na.rm = TRUE)) / (1 + length(null))
  } else {
    stop("alternative must be enrichment or depletion")
  }
}

outlier_levels <- list(
  strict = list(
    label = "strict_99",
    fst_prob = 0.99,
    fdm_tail_prob = 0.99,
    file_suffix = "",
    drop_traits = character(),
    fst_threshold_name = "hudson_fst_99th",
    fdm_threshold_name = "fdm_abs_negative_99th",
    dxy_threshold_name = "dxy_99th"
  ),
  p95 = list(
    label = "p95",
    fst_prob = 0.95,
    fdm_tail_prob = 0.95,
    file_suffix = "_p95",
    drop_traits = character(),
    fst_threshold_name = "hudson_fst_95th",
    fdm_threshold_name = "fdm_abs_negative_95th",
    dxy_threshold_name = "dxy_95th"
  ),
  p95_noSD = list(
    label = "p95_noSD",
    fst_prob = 0.95,
    fdm_tail_prob = 0.95,
    file_suffix = "_p95_noSD",
    drop_traits = "SD",
    fst_threshold_name = "hudson_fst_95th",
    fdm_threshold_name = "fdm_abs_negative_95th",
    dxy_threshold_name = "dxy_95th"
  )
)

run_outlier_level <- function(level) {
  fst_threshold <- as.numeric(quantile(fst_win$avg_hudson_fst, level$fst_prob, na.rm = TRUE))
  fdm_threshold <- fdm_positive_cutoff(level$fdm_tail_prob)
  dxy_threshold <- as.numeric(quantile(dxy_win$avg_dxy, level$fst_prob, na.rm = TRUE))
  level_peaks <- list(
    fst = fst_for_peaks %>%
      mutate(is_out = avg_hudson_fst >= fst_threshold) %>%
      merge_consecutive_runs(),
    f_dM = fdm_for_peaks %>%
      mutate(is_out = f_dM >= fdm_threshold) %>%
      merge_overlapping_outliers(),
    dxy = dxy_for_peaks %>%
      mutate(is_out = avg_dxy >= dxy_threshold) %>%
      merge_consecutive_runs()
  )
  stat_meta <- list(
    fst = list(
      alternative = "enrichment",
      threshold = fst_threshold,
      threshold_name = level$fst_threshold_name,
      tail_prob = level$fst_prob
    ),
    f_dM = list(
      alternative = "depletion",
      threshold = fdm_threshold,
      threshold_name = level$fdm_threshold_name,
      tail_prob = level$fdm_tail_prob
    ),
    dxy = list(
      alternative = "enrichment",
      threshold = dxy_threshold,
      threshold_name = level$dxy_threshold_name,
      tail_prob = level$fst_prob
    )
  )
  sets <- lapply(qtl_sets, function(d) d %>% filter(!trait %in% level$drop_traits))
  if (any(vapply(sets, nrow, integer(1)) == 0L)) {
    stop("Dropping traits ", paste(level$drop_traits, collapse = ", "), " removed every CI from a set.")
  }
  results <- list()
  summary_rows <- list()
  per_ci_rows <- list()
  for (set_name in names(sets)) {
    qtl <- sets[[set_name]]
    for (stat_name in names(level_peaks)) {
      cat(
        "Running ", level$fst_threshold_name, " / ", set_name, " / ", stat_name,
        " (", n_perm, " shifts)\n",
        sep = ""
      )
      obs <- observed_scores(qtl, level_peaks[[stat_name]])
      null <- null_scores(qtl, level_peaks[[stat_name]])
      alt <- stat_meta[[stat_name]]$alternative
      results[[set_name]][[stat_name]] <- list(
        alternative = alt,
        threshold = stat_meta[[stat_name]]$threshold,
        threshold_name = stat_meta[[stat_name]]$threshold_name,
        tail_prob = stat_meta[[stat_name]]$tail_prob,
        n_peaks_genome = nrow(level_peaks[[stat_name]]),
        observed_n_peaks = obs$n_peaks,
        observed_n_cis = obs$n_cis,
        observed_overlap_bp = obs$overlap_bp,
        null_n_peaks = null$n_peaks,
        null_n_cis = null$n_cis,
        null_overlap_bp = null$overlap_bp,
        p_n_peaks = one_sided_p(null$n_peaks, obs$n_peaks, alt),
        p_n_cis = one_sided_p(null$n_cis, obs$n_cis, alt),
        p_overlap_bp = one_sided_p(null$overlap_bp, obs$overlap_bp, alt)
      )
      summary_rows[[length(summary_rows) + 1L]] <- tibble(
        outlier_level = level$label,
        qtl_set = set_name,
        statistic = stat_name,
        alternative = alt,
        tail_prob = stat_meta[[stat_name]]$tail_prob,
        n_qtl = nrow(qtl),
        n_peaks_genome = nrow(level_peaks[[stat_name]]),
        threshold = stat_meta[[stat_name]]$threshold,
        observed_n_peaks = obs$n_peaks,
        null_mean_n_peaks = mean(null$n_peaks),
        p_n_peaks = results[[set_name]][[stat_name]]$p_n_peaks,
        observed_n_cis = obs$n_cis,
        p_n_cis = results[[set_name]][[stat_name]]$p_n_cis,
        observed_overlap_bp = obs$overlap_bp,
        p_overlap_bp = results[[set_name]][[stat_name]]$p_overlap_bp
      )
      per_ci_rows[[length(per_ci_rows) + 1L]] <- qtl %>%
        mutate(
          outlier_level = level$label,
          qtl_set = set_name,
          statistic = stat_name,
          ci_contains_peak = obs$per_ci
        ) %>%
        select(
          outlier_level, qtl_set, statistic, trait, chr, opposite_direction,
          lower95CI_bp, upper95CI_bp, ci_contains_peak
        )
    }
  }
  summary_tbl <- bind_rows(summary_rows)
  per_ci_tbl <- bind_rows(per_ci_rows)
  out <- list(
    n_perm = n_perm,
    seed = 251130L,
    chr_skipped = 4L,
    outlier_level = level$label,
    dropped_traits = level$drop_traits,
    fst_tail_prob = level$fst_prob,
    fdm_negative_tail_prob = level$fdm_tail_prob,
    fst_threshold = fst_threshold,
    fdm_positive_threshold = fdm_threshold,
    dxy_threshold = dxy_threshold,
    shift = "circular bp; each CI keeps its observed base-pair length",
    chr_length_bp = L_bp,
    fdm_trio = "brac_alle_vill (P1 bracteatus, P2 allenii, P3 villosissimus); positive f_dM is alle–vill sharing",
    note = paste(
      "Primary p-values are p_n_peaks. Fst and high Dxy are enrichment; high f_dM is depletion.",
      "High f_dM uses quantile(|f_dM|, tail_prob) among negative windows,",
      "equal to -quantile(f_dM, 1 - tail_prob) on that negative set.",
      "The four tests share one bp placement null and are not independent."
    ),
    qtl_sets = lapply(sets, function(d) d %>% select(trait, chr, opposite_direction, lower95CI_bp, upper95CI_bp)),
    peaks = level_peaks,
    summary = summary_tbl,
    per_ci = per_ci_tbl,
    tests = results
  )
  suffix <- level$file_suffix
  rds_path <- file.path(out_dir, paste0("circular_shift_outlier_qtl_overlap", suffix, ".rds"))
  summary_path <- file.path(out_dir, paste0("circular_shift_outlier_qtl_overlap_summary", suffix, ".csv"))
  per_ci_path <- file.path(out_dir, paste0("circular_shift_outlier_qtl_overlap_per_ci", suffix, ".csv"))
  saveRDS(out, rds_path)
  write.csv(summary_tbl, summary_path, row.names = FALSE)
  write.csv(per_ci_tbl, per_ci_path, row.names = FALSE)
  cat("\n", out$outlier_level, ": Fst threshold ", signif(fst_threshold, 4),
      ", positive f_dM threshold ", signif(fdm_threshold, 4),
      ", Dxy threshold ", signif(dxy_threshold, 4), "\n", sep = "")
  print(summary_tbl %>% select(qtl_set, statistic, alternative, observed_n_peaks, null_mean_n_peaks, p_n_peaks))
  cat("Saved:\n  ", rds_path, "\n  ", summary_path, "\n  ", per_ci_path, "\n", sep = "")
  invisible(summary_tbl)
}

cat("Primary statistic: number of outlier peaks that intersect a QTL CI.\n")
cat("Fst and Dxy alternative: enrichment (null >= observed). f_dM alternative: depletion (null <= observed).\n")
cat("Both outlier definitions use the same bp shifts.\n")
for (level in outlier_levels) run_outlier_level(level)
