# =============================================================================
# Supplementary noNearHybrid 1 km: Dsuite plots + full-vs-subset summary CSV
#
# Writes into this directory (does not overwrite main-repo Dsuite figures).
# Run from repo root:
#   Rscript 04_introgression_vs_QTL/noNearHybrid_1km/make_plots_and_summary.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
this_dir <- if (!is.na(file_arg) && nzchar(file_arg)) {
  dirname(normalizePath(file_arg))
} else {
  normalizePath(getwd())
}
project_root <- normalizePath(file.path(this_dir, "..", ".."), winslash = "/", mustWork = TRUE)
out_root <- normalizePath(this_dir, winslash = "/", mustWork = TRUE)

source(file.path(project_root, "scripts", "project_paths.R"))
paths <- project_paths()

full_dir <- file.path(
  project_root, "Dsuite",
  "Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120"
)
subset_dir <- file.path(project_root, "Dsuite", "dsuite_noNearHybrid_1km_w50s25")
if (!dir.exists(full_dir)) stop("Full-run Dsuite dir missing: ", full_dir)
if (!dir.exists(subset_dir)) stop("Subset Dsuite dir missing: ", subset_dir)

dir.create(file.path(out_root, "Dsuite"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "pixy"), recursive = TRUE, showWarnings = FALSE)

load_dsuite <- function(dir) {
  files <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)
  if (!length(files)) stop("No .txt in ", dir)
  dplyr::bind_rows(lapply(files, read.delim, sep = "\t", stringsAsFactors = FALSE)) %>%
    dplyr::mutate(
      chrN = as.integer(stringr::str_extract(chr, "\\d+")),
      dplyr::across(
        dplyr::any_of(c("D", "f_d", "f_dM", "d_f", "windowStart", "windowEnd")),
        ~ suppressWarnings(as.numeric(.x))
      )
    )
}

stat_summary_block <- function(x) {
  x <- x[is.finite(x)]
  neg <- x[x < 0]
  tibble::tibble(
    n = length(x),
    mean = if (length(x)) mean(x) else NA_real_,
    median = if (length(x)) median(x) else NA_real_,
    sd = if (length(x) > 1) sd(x) else NA_real_,
    mean_positive = {
      pos <- x[x > 0]
      if (length(pos)) mean(pos) else NA_real_
    },
    frac_positive = if (length(x)) mean(x > 0) else NA_real_,
    p99_abs_among_negative = if (length(neg)) as.numeric(quantile(abs(neg), 0.99)) else NA_real_
  )
}

summarise_run <- function(df, run_label) {
  stats <- c("D", "f_d", "f_dM", "d_f")
  scopes <- c(as.list(sort(unique(df$chrN[is.finite(df$chrN)]))), list("genome"))
  rows <- list()
  for (sc in scopes) {
    sub <- if (identical(sc, "genome")) df else df %>% dplyr::filter(chrN == sc)
    scope_name <- if (identical(sc, "genome")) "genome" else paste0("chr", sc)
    for (st in stats) {
      sm <- stat_summary_block(sub[[st]])
      rows[[length(rows) + 1]] <- dplyr::bind_cols(
        tibble::tibble(scope = scope_name, run = run_label, statistic = st),
        sm
      )
    }
  }
  dplyr::bind_rows(rows)
}

message("Loading full and subset Dsuite windows...")
full_df <- load_dsuite(full_dir)
subset_df <- load_dsuite(subset_dir)

long_sum <- dplyr::bind_rows(
  summarise_run(full_df, "full"),
  summarise_run(subset_df, "noNearHybrid_1km")
)

wide <- long_sum %>%
  tidyr::pivot_wider(
    id_cols = c(scope, statistic),
    names_from = run,
    values_from = c(n, mean, median, sd, mean_positive, frac_positive, p99_abs_among_negative),
    names_sep = "__"
  )

need <- c(
  "mean__full", "mean__noNearHybrid_1km",
  "mean_positive__full", "mean_positive__noNearHybrid_1km",
  "frac_positive__full", "frac_positive__noNearHybrid_1km"
)
missing <- setdiff(need, names(wide))
if (length(missing)) {
  stop(
    "Summary pivot missing columns: ", paste(missing, collapse = ", "),
    "\nHave: ", paste(names(wide), collapse = ", ")
  )
}

wide <- wide %>%
  dplyr::mutate(
    mean__delta_subset_minus_full = mean__noNearHybrid_1km - mean__full,
    mean_positive__delta_subset_minus_full =
      mean_positive__noNearHybrid_1km - mean_positive__full,
    frac_positive__delta_subset_minus_full =
      frac_positive__noNearHybrid_1km - frac_positive__full,
    scope_ord = ifelse(scope == "genome", 99L, as.integer(sub("^chr", "", scope))),
    stat_ord = match(statistic, c("D", "f_d", "f_dM", "d_f"))
  ) %>%
  dplyr::arrange(scope_ord, stat_ord) %>%
  dplyr::select(-scope_ord, -stat_ord)

sum_csv <- file.path(out_root, "dsuite_full_vs_noNearHybrid_1km_summary.csv")
long_csv <- file.path(out_root, "dsuite_full_vs_noNearHybrid_1km_summary_long.csv")
readr::write_csv(wide, sum_csv)
readr::write_csv(long_sum, long_csv)
message("Wrote ", sum_csv)
message("Wrote ", long_csv)

# ---- Rerun plot scripts with subset Dsuite, writing into out_root ----
pixy_folder <- if (!is.na(paths$pixy_folder) && dir.exists(paths$pixy_folder)) {
  paths$pixy_folder
} else {
  file.path(project_root, "pixy", "pixy_hudson_10kb")
}
pixy_500 <- if (!is.na(paths$pixy_folder_500kb) && dir.exists(paths$pixy_folder_500kb)) {
  paths$pixy_folder_500kb
} else {
  file.path(project_root, "pixy", "pixy_hudson_500kb")
}

Sys.setenv(
  PROJECT_ROOT = project_root,
  DSUITE_DIR = subset_dir,
  PLOT_OUT_ROOT = out_root,
  PIXY_FOLDER = pixy_folder,
  PIXY_FOLDER_500KB = pixy_500
)
if (!nzchar(Sys.getenv("N_PERM", ""))) Sys.setenv(N_PERM = "10000")

run_rscript <- function(rel) {
  path <- file.path(project_root, rel)
  if (!file.exists(path)) stop("Missing script: ", path)
  message("\n===== ", rel, " =====")
  status <- system2("Rscript", args = path, stdout = "", stderr = "")
  if (!identical(as.integer(status), 0L)) stop("Script failed (", status, "): ", rel)
  message("===== done =====")
}

run_rscript("04_introgression_vs_QTL/genome-wide_f_dM_dsuite_w_QTL_plots.R")
run_rscript("04_introgression_vs_QTL/genome-wide_d_f_dsuite_w_QTL_plots.R")
run_rscript("04_introgression_vs_QTL/combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R")
run_rscript("04_introgression_vs_QTL/genome-wide_Fst.Dxy_w_QTL_plots.R")
run_rscript("04_introgression_vs_QTL/circular_shift_outlier_qtl_overlap.R")

# fd_fdm_vs_fst_100kb_bins.R still uses WC pixy paths hardcoded; patch a temp copy.
message("\n===== fd_fdm_vs_fst_100kb_bins.R (subset Dsuite) =====")
wc_pixy <- file.path(project_root, "pixy", "output_251124_LASIref")
if (dir.exists(wc_pixy)) {
  src <- readLines(file.path(project_root, "fd_fdm_vs_fst_100kb_bins.R"))
  src <- sub(
    '^project_root <- ".*"$',
    paste0('project_root <- "', out_root, '"'),
    src
  )
  src <- sub(
    "dsuite_dir <- file\\.path\\(\\s*project_root,\\s*\"Dsuite\",\\s*\"Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120\"\\s*\\)",
    paste0('dsuite_dir <- "', subset_dir, '"'),
    src
  )
  src <- sub(
    'pixy_folder <- file\\.path\\(project_root, "pixy", "output_251124_LASIref"\\)',
    paste0('pixy_folder <- "', wc_pixy, '"'),
    src
  )
  tmp <- file.path(out_root, "_tmp_fd_fdm_vs_fst_100kb_bins.R")
  writeLines(src, tmp)
  status <- system2("Rscript", args = c(tmp, "stride2", "positive"), stdout = "", stderr = "")
  unlink(tmp)
  if (!identical(status, 0L)) {
    warning("fd_fdm_vs_fst_100kb_bins failed with status ", status)
  } else {
    message("===== done =====")
  }
} else {
  message("Skipping fd_fdm_vs_fst: WC pixy folder missing at ", wc_pixy)
}

message("\nAll done. Outputs in:\n  ", out_root)
message("Summary CSV:\n  ", sum_csv)
