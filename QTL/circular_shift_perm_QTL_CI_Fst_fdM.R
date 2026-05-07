# =============================================================================
# Circular (shift) permutation of QTL 95% CIs in cM along each chromosome,
# wrapping at chromosome ends. Builds null distributions for mean Fst and mean
# f_dM averaged over QTL windows (cM overlap; see Fst/f_dM notes in script body).
# Chromosome 4 is excluded from map data, pixy/Dsuite windows, and QTL rows.
#
# Requires: lepmap marker files (same pattern as genome-wide_Fst.Dxy_w_QTL_plots.R).
# Set LEPMAP_DIR if the default sibling path does not exist.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

project_root <- "/Users/admin-jgharenc/Documents/Github/Costus_popgen_clean"
qtl_ci_path <- file.path(project_root, "QTL", "all_QTL_CIs_bp_cm.251130.csv")
pixy_folder <- file.path(project_root, "pixy", "output_251124_LASIref")
dsuite_dir <- file.path(
  project_root,
  "Dsuite",
  "Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120"
)

lepmap_dir <- Sys.getenv(
  "LEPMAP_DIR",
  unset = normalizePath(
    file.path(project_root, "..", "alle_vill_QTL", "lepmap3", "best_iter_ordermarkers"),
    mustWork = FALSE
  )
)

n_perm <- as.integer(Sys.getenv("N_PERM", unset = "10000"))
if (!is.finite(n_perm) || n_perm < 1L) {
  stop("N_PERM must be a positive integer (default 10000).")
}
set.seed(251130L)

# ---- QTL filters (match genome-wide_Fst.Dxy_w_QTL_plots.R) ----
all_QTL_CI <- read.csv(qtl_ci_path, header = TRUE, stringsAsFactors = FALSE) %>%
  filter(!chr == 2 | !trait == "Nreabs") %>%
  filter(!chr == 5 | !trait == "GR") %>%
  filter(as.integer(chr) != 4L)

all_QTL_CI_fst <- all_QTL_CI %>% filter(!chr == 6 | !trait == "SD")

# ---- Lepmap: bp <-> cM (chr 4 dropped) ----
lepmap_files <- list.files(
  lepmap_dir,
  pattern = "*_clean_LodLim11_theta0.001_allF2s_0.5covCO_withRef.txt",
  full.names = TRUE
)
if (length(lepmap_files) == 0) {
  stop(
    "No lepmap marker files found in lepmap_dir = ", lepmap_dir,
    ". Set env LEPMAP_DIR to the folder containing *_clean_LodLim11_...txt files."
  )
}

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

chr_len_cm <- bp_to_cm_iso %>%
  group_by(chrom) %>%
  summarise(L_cm = max(cM_iso, na.rm = TRUE), .groups = "drop")

bp_to_cm_funs <- bp_to_cm_iso %>%
  group_by(chrom) %>%
  group_map(~ approxfun(.x$bp_pos, .x$cM_iso, rule = 2))
names(bp_to_cm_funs) <- as.character(sort(unique(bp_to_cm_iso$chrom)))

# Inverse map (unique cM positions; avoids approx duplicate-x warnings)
cm_to_bp_funs <- bp_to_cm_iso %>%
  group_by(chrom) %>%
  group_map(~ {
    o <- order(.x$cM_iso, .x$bp_pos)
    cx <- .x$cM_iso[o]
    bx <- .x$bp_pos[o]
    keep <- !duplicated(cx)
    approxfun(cx[keep], bx[keep], rule = 2)
  })
names(cm_to_bp_funs) <- as.character(sort(unique(bp_to_cm_iso$chrom)))

# ---- Pixy -> Fst windows with cM (exclude chr 4) ----
pixy_to_long <- function(pixy_files) {
  pixy_df <- list()
  for (i in seq_along(pixy_files)) {
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
    if (stat_file_type == "pi") {
      df <- readr::read_delim(pixy_files[i], delim = "\t", show_col_types = FALSE)
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
          key = "statistic", value = "value"
        ) %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      pixy_df[[i]] <- df
    } else {
      df <- readr::read_delim(pixy_files[i], delim = "\t", show_col_types = FALSE)
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
pixy_df <- pixy_to_long(pixy_files)
pixy_wide <- pivot_wider(pixy_df, names_from = statistic, values_from = value) %>%
  mutate(chrN = as.numeric(str_extract(chromosome, "\\d+")))

fst <- pixy_wide %>%
  filter(is.finite(avg_wc_fst), chromosome != "Chrom4") %>%
  mutate(
    cM_window_start = map2_dbl(chrN, window_pos_1, function(cn, bp) {
      bp_to_cm_funs[[as.character(cn)]](bp)
    })
  )

fst_cm <- fst %>%
  arrange(chrN, cM_window_start) %>%
  group_by(chrN) %>%
  mutate(cM_window_end = lead(cM_window_start)) %>%
  ungroup()

# ---- Dsuite f_dM (exclude chr 4); overlap via bp after mapping shifted CI ----
file_list <- list.files(path = dsuite_dir, pattern = "*.txt", full.names = TRUE)
if (length(file_list) == 0) {
  stop("No Dsuite *.txt files in ", dsuite_dir)
}
LASIout <- data.frame()
for (file in file_list) {
  LASIout <- rbind(LASIout, read.delim(file, sep = "\t"))
}
pos_f_dM_clean <- LASIout %>%
  select(chr, windowStart, windowEnd, f_dM) %>%
  mutate(f_dM_pos = ifelse(f_dM > 0, f_dM, NA_real_)) %>%
  filter(is.finite(f_dM_pos)) %>%
  mutate(chrN = as.numeric(str_extract(chr, "\\d+"))) %>%
  filter(chrN != 4)

segments_cm_to_bp_ranges <- function(chr_num, segments_cm) {
  inv <- cm_to_bp_funs[[as.character(chr_num)]]
  lapply(segments_cm, function(seg) {
    p1 <- inv(seg[1])
    p2 <- inv(seg[2])
    range(c(p1, p2))
  })
}

mean_fdm_over_bp_ranges <- function(chr_num, bp_ranges) {
  sub <- pos_f_dM_clean %>% filter(chrN == chr_num)
  if (nrow(sub) == 0) {
    return(NA_real_)
  }
  ws <- sub$windowStart
  we <- sub$windowEnd
  hits <- rep(FALSE, nrow(sub))
  for (rg in bp_ranges) {
    hits <- hits | (ws <= rg[2] & we >= rg[1])
  }
  mean(sub$f_dM_pos[hits], na.rm = TRUE)
}

# ---- Shift CI in cM; return list of 1 or 2 intervals [lo, hi] in [0, L) ----
shift_ci_cm <- function(lo, hi, shift, L) {
  if (hi <= lo) {
    stop("QTL CI must have lower95CI_cm < upper95CI_cm")
  }
  if (hi - lo > L) {
    stop("CI length exceeds chromosome length in cM")
  }
  a <- (lo + shift) %% L
  b <- (hi + shift) %% L
  if (a <= b) {
    list(c(a, b))
  } else {
    list(c(a, L), c(0, b))
  }
}

mean_fst_over_ci_segments <- function(chr_num, segments_cm) {
  df <- fst_cm %>% filter(chrN == chr_num)
  if (nrow(df) == 0) {
    return(NA_real_)
  }
  w1 <- df$cM_window_start
  w2 <- df$cM_window_end
  hits <- rep(FALSE, nrow(df))
  for (seg in segments_cm) {
    hits <- hits | (w1 <= seg[2] & w2 >= seg[1] & !is.na(w2))
  }
  mean(df$avg_wc_fst[hits], na.rm = TRUE)
}

# ---- Observed statistics: Fst in cM (matches permutation); f_dM in bp (genome-wide script) ----
# Helper functions use chr_num (not chr) so rowwise() mutate can see column chr without masking.
obs_fst <- all_QTL_CI_fst %>%
  rowwise() %>%
  mutate(
    avg_fst_QTL = mean_fst_over_ci_segments(
      as.integer(chr),
      list(c(lower95CI_cm, upper95CI_cm))
    )
  ) %>%
  ungroup()

obs_fdm <- all_QTL_CI_fst %>%
  rowwise() %>%
  mutate(
    avg_f_dM_QTL = mean_fdm_over_bp_ranges(
      as.integer(chr),
      list(c(lower95CI_bp, upper95CI_bp))
    )
  ) %>%
  ungroup()

obs_mean_fst <- mean(obs_fst$avg_fst_QTL, na.rm = TRUE)
obs_mean_fdm <- mean(obs_fdm$avg_f_dM_QTL, na.rm = TRUE)

# Precompute per-chromosome shift draws: independent Uniform(0, L) per chr that has QTLs
chr_with_qtl <- sort(unique(as.integer(all_QTL_CI_fst$chr)))
L_vec <- chr_len_cm %>%
  filter(chrom %in% chr_with_qtl) %>%
  { stats::setNames(.$L_cm, as.character(.$chrom)) }

# ---- Permutation ----
null_mean_fst <- numeric(n_perm)
null_mean_fdm <- numeric(n_perm)

for (p in seq_len(n_perm)) {
  shifts <- map_dbl(names(L_vec), function(ch) runif(1, min = 0, max = L_vec[[ch]]))
  names(shifts) <- names(L_vec)

  fst_vals <- map2_dbl(
    all_QTL_CI_fst$chr,
    seq_len(nrow(all_QTL_CI_fst)),
    function(chr_val, i) {
      row <- all_QTL_CI_fst[i, , drop = FALSE]
      ch <- as.integer(chr_val)
      L <- L_vec[[as.character(ch)]]
      if (is.null(L) || is.na(L) || L <= 0) {
        return(NA_real_)
      }
      s <- shifts[[as.character(ch)]]
      segs <- shift_ci_cm(row$lower95CI_cm, row$upper95CI_cm, s, L)
      mean_fst_over_ci_segments(ch, segs)
    }
  )

  fdm_vals <- map2_dbl(
    all_QTL_CI_fst$chr,
    seq_len(nrow(all_QTL_CI_fst)),
    function(chr_val, i) {
      row <- all_QTL_CI_fst[i, , drop = FALSE]
      ch <- as.integer(chr_val)
      L <- L_vec[[as.character(ch)]]
      if (is.null(L) || is.na(L) || L <= 0) {
        return(NA_real_)
      }
      s <- shifts[[as.character(ch)]]
      segs <- shift_ci_cm(row$lower95CI_cm, row$upper95CI_cm, s, L)
      bp_rng <- segments_cm_to_bp_ranges(ch, segs)
      mean_fdm_over_bp_ranges(ch, bp_rng)
    }
  )

  null_mean_fst[p] <- mean(fst_vals, na.rm = TRUE)
  null_mean_fdm[p] <- mean(fdm_vals, na.rm = TRUE)
}

p_greater_fst <- (1 + sum(null_mean_fst >= obs_mean_fst, na.rm = TRUE)) / (1 + n_perm)
p_greater_fdm <- (1 + sum(null_mean_fdm >= obs_mean_fdm, na.rm = TRUE)) / (1 + n_perm)

out <- list(
  n_perm = n_perm,
  chr_skipped = 4L,
  obs_mean_fst = obs_mean_fst,
  obs_mean_fdm = obs_mean_fdm,
  null_mean_fst = null_mean_fst,
  null_mean_fdm = null_mean_fdm,
  p_one_sided_mean_fst_ge_obs = p_greater_fst,
  p_one_sided_mean_fdm_ge_obs = p_greater_fdm,
  obs_per_qtl = tibble(
    trait = all_QTL_CI_fst$trait,
    chr = all_QTL_CI_fst$chr,
    avg_fst_QTL = obs_fst$avg_fst_QTL,
    avg_f_dM_QTL = obs_fdm$avg_f_dM_QTL
  ),
  chr_lengths_cm = enframe(L_vec, name = "chr", value = "L_cm")
)

saveRDS(out, file.path(project_root, "QTL", "circular_shift_perm_QTL_CI_Fst_fdM.rds"))

cat(
  "Circular shift permutation (cM), chr 4 excluded.\n",
  "  n_perm = ", n_perm, "\n",
  "  mean Fst across QTL windows (obs): ", sprintf("%.6f", obs_mean_fst),
  "  p(>= obs) = ", sprintf("%.5f", p_greater_fst), "\n",
  "  mean f_dM across QTL windows (obs): ", sprintf("%.6f", obs_mean_fdm),
  "  p(>= obs) = ", sprintf("%.5f", p_greater_fdm), "\n",
  "  Saved: ", file.path(project_root, "QTL", "circular_shift_perm_QTL_CI_Fst_fdM.rds"), "\n",
  sep = ""
)
