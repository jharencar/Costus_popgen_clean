# =============================================================================
# Fully subset-based combo 4-panel + circular-shift (noNearHybrid 1 km)
#
# Requires local copies of:
#   Dsuite/dsuite_noNearHybrid_1km_w50s25/
#   pixy/pixy_hudson_10kb_noNearHybrid_1km/
#   pixy/pixy_hudson_500kb_noNearHybrid_1km/   (optional for some plots)
#
# Run from repo root:
#   Rscript 04_introgression_vs_QTL/noNearHybrid_1km/remake_combo_and_circular_shift.R
# =============================================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE)[1])
this_dir <- if (!is.na(file_arg) && nzchar(file_arg)) {
  dirname(normalizePath(file_arg))
} else {
  normalizePath(getwd())
}
project_root <- normalizePath(file.path(this_dir, "..", ".."), winslash = "/", mustWork = TRUE)
out_root <- normalizePath(this_dir, winslash = "/", mustWork = TRUE)

subset_dsuite <- file.path(project_root, "Dsuite", "dsuite_noNearHybrid_1km_w50s25")
subset_pixy_10kb <- file.path(project_root, "pixy", "pixy_hudson_10kb_noNearHybrid_1km")
subset_pixy_500kb <- file.path(project_root, "pixy", "pixy_hudson_500kb_noNearHybrid_1km")

stopifnot(dir.exists(subset_dsuite))
stopifnot(dir.exists(subset_pixy_10kb))
if (!dir.exists(subset_pixy_500kb)) {
  warning("500 kb subset pixy folder missing: ", subset_pixy_500kb)
}

Sys.setenv(
  PROJECT_ROOT = project_root,
  DSUITE_DIR = subset_dsuite,
  PIXY_FOLDER = subset_pixy_10kb,
  PIXY_FOLDER_500KB = if (dir.exists(subset_pixy_500kb)) subset_pixy_500kb else "",
  PLOT_OUT_ROOT = out_root,
  N_PERM = Sys.getenv("N_PERM", "10000")
)

run_rscript <- function(rel) {
  path <- file.path(project_root, rel)
  if (!file.exists(path)) stop("Missing script: ", path)
  message("\n===== ", rel, " =====")
  status <- system2("Rscript", args = path, stdout = "", stderr = "")
  if (!identical(as.integer(status), 0L)) stop("Script failed (", status, "): ", rel)
  message("===== done =====")
}

# Combo 4-panel with subset Dsuite + subset pixy
run_rscript("04_introgression_vs_QTL/combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R")

# Circular-shift QTL overlap with subset Fst/Dxy/fdM
run_rscript("04_introgression_vs_QTL/circular_shift_outlier_qtl_overlap.R")

message("\nWrote subset combo + circular-shift outputs under:\n  ", out_root)
