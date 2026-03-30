# QTL credible intervals in bp and cM — print table
#
# Source of the numbers: Costus_QTL.qmd — function extract_bayesint_bp_cm() applies
# rqtl::bayesint(..., expandtomarkers = TRUE): peak bp from marker names, peak/cI cM from
# bi$pos; CI bounds in bp are the flanking markers (rows 1 and 3 of the 3-row bayesint
# output). The combined table is written as all_QTL_CIs_bp_cm.251130.csv (see that chunk).

suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))

csv_path <- file.path("QTL", "all_QTL_CIs_bp_cm.251130.csv")
if (!file.exists(csv_path)) {
  stop("File not found: ", normalizePath(csv_path, mustWork = FALSE),
       "\nRun this script from the repository root (Population_genetics_vill_alle.Rproj).",
       call. = FALSE)
}

qtl_ci <- read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)

fmt_bp <- function(x) format(round(as.numeric(x)), big.mark = ",")
fmt_cm <- function(x) sprintf("%.3f", as.numeric(x))

bp_tbl <- qtl_ci %>%
  transmute(
    trait,
    chr = as.integer(chr),
    peak_bp = fmt_bp(peak_bp),
    lower95CI_bp = fmt_bp(lower95CI_bp),
    upper95CI_bp = fmt_bp(upper95CI_bp),
    CI_size_bp = fmt_bp(CI_size_bp)
  )

cm_tbl <- qtl_ci %>%
  transmute(
    trait,
    chr = as.integer(chr),
    peak_cm = fmt_cm(peak_cm),
    lower95CI_cm = fmt_cm(lower95CI_cm),
    upper95CI_cm = fmt_cm(upper95CI_cm),
    CI_size_cm = fmt_cm(CI_size_cm)
  )

cat("\nQTL peaks and 95% credible intervals\n")
cat("Source CSV:", normalizePath(csv_path), "\n\n")

cat("--- Base pairs ---\n")
print(as.data.frame(bp_tbl), row.names = FALSE, right = FALSE)

cat("\n--- Centimorgans (LepMap / rqtl map positions) ---\n")
print(as.data.frame(cm_tbl), row.names = FALSE, right = FALSE)
