# Final Waters mass-alignment chemistry QTL mapping.
# Run from project root or with working directory = QTL/ (script will switch to QTL/ if needed).
# Loads canonical most_dups genoprob cross, injects scaled chemistry phenotypes,
# per-trait AIC covariate selection, scanone (HK) + permutations for 10% FDR.

library(dplyr)
library(qtl)

for (p in c("log_all_warnings.R", file.path("QTL", "log_all_warnings.R"))) {
  if (file.exists(p)) {
    source(p, local = FALSE)
    break
  }
}
if (!exists("start_warning_log")) {
  stop("log_all_warnings.R not found in working directory or QTL/", call. = FALSE)
}
start_warning_log(default_basename = "WatersMass_chemistry_warnings.txt")

# ---- Config (override via environment) ----
default_genoprob <- "/Users/admin-jgharenc/Documents/Github/alle_vill_QTL/aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds"
genoprob_path <- Sys.getenv("GENOPROB_RDS", unset = default_genoprob)
n_perm <- as.integer(Sys.getenv("N_PERM", unset = "100"))
run_multi_qtl <- identical(Sys.getenv("RUN_MULTI_QTL", unset = "0"), "1")

phe_merged <- "F3phe_onehotcovar_newHI_with_WatersMass_chemistry.csv"
meta_file <- "WatersMass_feature_metadata.csv"

if (!file.exists(phe_merged) && file.exists(file.path("QTL", phe_merged))) {
  setwd("QTL")
  message("Set working directory to QTL/")
}
stopifnot(file.exists(phe_merged), file.exists(meta_file))

if (!file.exists(genoprob_path)) {
  stop(
    "Canonical genoprob not found: ", genoprob_path,
    "\nSet GENOPROB_RDS to the path of most_dups_genoprob.rds (see NOTES.md).",
    call. = FALSE
  )
}

message("GENOPROB_RDS: ", genoprob_path)
message("N_PERM: ", n_perm, if (run_multi_qtl) " | RUN_MULTI_QTL=1" else "")

# ---- Load cross and inject chemistry phenotypes ----
Costus_prob <- readRDS(genoprob_path)
phe <- read.csv(phe_merged, check.names = FALSE)
feature_meta <- read.csv(meta_file, check.names = FALSE)

#' Resolve individual IDs on cross$pheno (rownames are often 1..n; Lab IDs live in ID column).
resolve_cross_ids <- function(cross, phe_ids) {
  pheno <- cross$pheno
  phe_ids <- as.character(phe_ids)
  candidates <- list()
  if (!is.null(rownames(pheno)) && length(rownames(pheno)) == nrow(pheno)) {
    candidates$rownames <- as.character(rownames(pheno))
  }
  if ("ID" %in% colnames(pheno)) {
    candidates$ID <- as.character(pheno$ID)
  }
  if (length(candidates) == 0) {
    stop("Cannot determine individual IDs: no rownames or ID column on cross$pheno.", call. = FALSE)
  }
  rates <- vapply(candidates, function(ids) mean(ids %in% phe_ids), numeric(1))
  src <- names(which.max(rates))
  ids <- candidates[[src]]
  message(
    "Cross ID source: ", src, " (", round(100 * rates[[src]], 1), "% match to phenotype file)"
  )
  if (rates[[src]] < 0.5) {
    warning(
      "Low phenotype ID match rate (", round(100 * rates[[src]], 1), "%). ",
      "Check that F3phe ID column matches cross individual IDs.",
      call. = FALSE
    )
  }
  ids
}

chem_trait_names <- intersect(
  feature_meta$trait,
  names(phe)[grepl("^row[0-9]+$", names(phe))]
)
if (length(chem_trait_names) == 0) {
  stop("No row### chemistry traits found in ", phe_merged, ". Run prepare_WatersMass_chemistry.R first.")
}

cross_ids <- resolve_cross_ids(Costus_prob, phe$ID)
phe_idx <- match(cross_ids, as.character(phe$ID))
n_unmatched <- sum(is.na(phe_idx))
if (n_unmatched > 0) {
  warning("Cross individuals without merged phenotype row: ", n_unmatched)
}

for (tr in chem_trait_names) {
  vals <- rep(NA_real_, length(cross_ids))
  vals[!is.na(phe_idx)] <- as.numeric(phe[[tr]][phe_idx[!is.na(phe_idx)]])
  Costus_prob$pheno[, tr] <- vals
}
n_markers_raw <- tryCatch(nmar(Costus_prob), error = function(e) NULL)
n_markers <- if (is.null(n_markers_raw)) {
  sum(vapply(Costus_prob$geno, function(g) ncol(g$data), integer(1)))
} else if (length(n_markers_raw) == 1L) {
  as.integer(n_markers_raw)
} else {
  sum(as.integer(n_markers_raw))
}
message(
  "Loaded genoprob: ", nind(Costus_prob), " individuals, ",
  n_markers, " markers; ", length(chem_trait_names), " chemistry traits in phenotype file."
)

# ---- Mappable traits (variance > 0, >= 20 non-NA among cross individuals) ----
min_complete <- 20
mappable_traits <- character(0)
for (tr in chem_trait_names) {
  y <- as.numeric(Costus_prob$pheno[, tr])
  if (sum(!is.na(y)) >= min_complete && var(y, na.rm = TRUE) > 0) {
    mappable_traits <- c(mappable_traits, tr)
  }
}
message("Mappable chemistry traits: ", length(mappable_traits), " of ", length(chem_trait_names))
if (length(mappable_traits) == 0) stop("No mappable chemistry traits after QC filters.")

# ---- Covariate helpers (same logic as prelim_UV_chemistry_QTL.R / Costus_QTL.qmd) ----
min_covar_n <- 10
pn <- names(Costus_prob$pheno)
fam_onehot <- c("F21_238", "F21_259", "F21_241", "F21_284", "F21_232")
coh_cols <- pn[grepl("^coh[0-9]+$", pn)]
coh_n <- colSums(Costus_prob$pheno[, coh_cols, drop = FALSE], na.rm = TRUE)
coh_ok <- coh_cols[coh_n >= min_covar_n]
candidate_covars <- c("hybrid_index", fam_onehot[fam_onehot %in% pn], coh_ok)
candidate_covars <- candidate_covars[candidate_covars %in% pn]
message(
  "Candidate covariates (n=", length(candidate_covars), "): ",
  paste(head(candidate_covars, 8), collapse = ", "),
  if (length(candidate_covars) > 8) " ..." else ""
)

fam_onehot_all <- c("F21_238", "F21_259", "F21_241", "F21_284", "F21_232")
ensure_full_rank_covariate_names <- function(names) {
  out <- names
  fam_in <- fam_onehot_all[fam_onehot_all %in% out]
  if (length(fam_in) >= 5) out <- setdiff(out, fam_onehot_all[5])
  coh_in <- out[grepl("^coh[0-9]+$", out)]
  if (length(coh_in) >= 2) out <- setdiff(out, coh_in[length(coh_in)])
  out
}

select_covariates_for_trait <- function(cross, pheno.col, candidate_names, min_complete = 20) {
  y <- as.numeric(qtl::pull.pheno(cross, pheno.col))
  X <- as.data.frame(qtl::pull.pheno(cross, candidate_names))
  ok <- complete.cases(y, X)
  if (sum(ok) < min_complete) return(NULL)
  dat <- cbind(y = y[ok], X[ok, , drop = FALSE])
  fit <- lm(y ~ ., data = dat)
  fit_step <- step(fit, trace = 0)
  sel <- setdiff(names(coef(fit_step)), "(Intercept)")
  if (length(sel) == 0) return(NULL)
  sel <- ensure_full_rank_covariate_names(sel)
  if (length(sel) == 0) return(NULL)
  covar <- qtl::pull.pheno(cross, sel)
  if (length(sel) == 1) covar <- as.matrix(covar, ncol = 1, dimnames = list(NULL, sel))
  covar
}

# ---- scanone per trait ----
scanone_chemistry <- list()
covar_by_trait <- list()
failed_traits <- character(0)
n_ind <- nind(Costus_prob)

for (trait in mappable_traits) {
  covar <- select_covariates_for_trait(Costus_prob, trait, candidate_covars)
  covar_by_trait[[trait]] <- if (!is.null(covar)) colnames(covar) else character(0)
  addcovar <- covar
  if (!is.null(addcovar) && ncol(addcovar) >= n_ind - 5) addcovar <- NULL
  so <- tryCatch(
    scanone(Costus_prob, method = "hk", pheno.col = trait, addcovar = addcovar),
    error = function(e) {
      covar_by_trait[[trait]] <<- character(0)
      tryCatch(
        scanone(Costus_prob, method = "hk", pheno.col = trait, addcovar = NULL),
        error = function(e2) {
          failed_traits <<- c(failed_traits, trait)
          NULL
        }
      )
    }
  )
  if (!is.null(so)) scanone_chemistry[[trait]] <- so
}
if (length(failed_traits) > 0) {
  warning(
    "scanone failed for ", length(failed_traits), " trait(s): ",
    paste(head(failed_traits, 5), collapse = ", "),
    if (length(failed_traits) > 5) " ..." else ""
  )
}
message(
  "Completed scanone for ", length(scanone_chemistry), " traits.",
  if (length(failed_traits) > 0) paste0(" Skipped ", length(failed_traits), ".") else ""
)
if (length(scanone_chemistry) == 0) stop("All scanone runs failed.")

summ <- do.call(rbind, lapply(names(scanone_chemistry), function(tr) {
  s <- scanone_chemistry[[tr]]
  lod_col <- setdiff(names(s), c("chr", "pos"))
  lod_vals <- s[[lod_col]]
  w <- which.max(lod_vals)
  data.frame(
    trait = tr,
    chr = s$chr[w],
    pos = s$pos[w],
    lod = lod_vals[w],
    marker = rownames(s)[w],
    n_covar = length(covar_by_trait[[tr]]),
    stringsAsFactors = FALSE
  )
}))
summ <- summ[order(-summ$lod), ]
rownames(summ) <- NULL
summ <- dplyr::left_join(summ, feature_meta, by = "trait")
message("Top 10 peaks (by LOD):")
print(head(summ, 10))

saveRDS(scanone_chemistry, "WatersMass_chemistry_scanone_by_trait.rds")
saveRDS(covar_by_trait, "WatersMass_chemistry_covariates_by_trait.rds")
write.csv(summ, "WatersMass_chemistry_scanone_summary.csv", row.names = FALSE)

# ---- Permutations for 10% and 50% LOD thresholds ----
perm_cutoff_by_trait <- setNames(rep(NA_real_, length(scanone_chemistry)), names(scanone_chemistry))
perm_cutoff_50_by_trait <- setNames(rep(NA_real_, length(scanone_chemistry)), names(scanone_chemistry))

if (n_perm > 0) {
  message("Running ", n_perm, " permutations per trait for LOD thresholds (10% and 50%)...")
  for (tr in names(scanone_chemistry)) {
    addcovar <- if (length(covar_by_trait[[tr]]) > 0) {
      qtl::pull.pheno(Costus_prob, covar_by_trait[[tr]])
    } else {
      NULL
    }
    perm_tr <- tryCatch(
      scanone(Costus_prob, method = "hk", pheno.col = tr, addcovar = addcovar, n.perm = n_perm),
      error = function(e) NULL
    )
    if (!is.null(perm_tr)) {
      thresh10 <- summary(perm_tr, alpha = 0.1)
      perm_cutoff_by_trait[[tr]] <- if (is.data.frame(thresh10)) thresh10[1, 1] else as.numeric(thresh10)[1]
      thresh50 <- summary(perm_tr, alpha = 0.5)
      perm_cutoff_50_by_trait[[tr]] <- if (is.data.frame(thresh50)) thresh50[1, 1] else as.numeric(thresh50)[1]
    }
  }
  message("Permutation cutoffs computed.")
}

significant_10 <- summ
significant_10$cutoff_10 <- perm_cutoff_by_trait[significant_10$trait]
if (n_perm > 0) {
  significant_10 <- significant_10[
    significant_10$lod > significant_10$cutoff_10 & !is.na(significant_10$cutoff_10),
  ]
  significant_10 <- significant_10[order(-significant_10$lod), ]
} else {
  significant_10 <- significant_10[0, , drop = FALSE]
  message("n_perm=0: skipping significant_10pct filtering (no permutation thresholds).")
}
write.csv(significant_10, "WatersMass_chemistry_significant_10pct.csv", row.names = FALSE)
message("Significant traits at 10% FDR: ", nrow(significant_10))

# ---- Shared peaks above 50% LOD threshold ----
peak_above_50 <- summ
peak_above_50$cutoff_50 <- perm_cutoff_50_by_trait[peak_above_50$trait]
peak_above_50 <- peak_above_50[
  peak_above_50$lod > peak_above_50$cutoff_50 & !is.na(peak_above_50$cutoff_50),
]
peak_above_50 <- peak_above_50[order(peak_above_50$chr, peak_above_50$pos), ]
cM_window <- 15

if (nrow(peak_above_50) > 0) {
  peak_above_50$group_id <- NA_integer_
  grp <- 0
  for (i in seq_len(nrow(peak_above_50))) {
    if (!is.na(peak_above_50$group_id[i])) next
    grp <- grp + 1
    peak_above_50$group_id[i] <- grp
    chr_i <- peak_above_50$chr[i]
    pos_i <- peak_above_50$pos[i]
    same_region <- peak_above_50$chr == chr_i & abs(peak_above_50$pos - pos_i) <= cM_window
    peak_above_50$group_id[same_region] <- grp
  }
  repeat {
    merged <- FALSE
    for (g in unique(peak_above_50$group_id)) {
      w <- which(peak_above_50$group_id == g)
      chrs <- peak_above_50$chr[w]
      pos <- peak_above_50$pos[w]
      for (j in seq_len(nrow(peak_above_50))) {
        if (peak_above_50$group_id[j] == g) next
        near <- any(peak_above_50$chr[j] == chrs & abs(peak_above_50$pos[j] - pos) <= cM_window)
        if (near) {
          old_g <- peak_above_50$group_id[j]
          peak_above_50$group_id[peak_above_50$group_id == old_g] <- g
          merged <- TRUE
        }
      }
    }
    if (!merged) break
  }
  shared_peaks <- aggregate(
    list(traits = peak_above_50$trait),
    by = list(chr = peak_above_50$chr, group_id = peak_above_50$group_id),
    FUN = function(x) paste(x, collapse = ", ")
  )
  shared_peaks$n_traits <- lengths(strsplit(shared_peaks$traits, ", "))
  pos_med <- aggregate(pos ~ chr + group_id, data = peak_above_50, FUN = median)
  shared_peaks$pos_median_cM <- pos_med$pos[
    match(paste(shared_peaks$chr, shared_peaks$group_id), paste(pos_med$chr, pos_med$group_id))
  ]
  shared_peaks <- shared_peaks[order(-shared_peaks$n_traits), c("chr", "pos_median_cM", "n_traits", "traits")]
  write.csv(shared_peaks, "WatersMass_chemistry_shared_peaks_above50pct.csv", row.names = FALSE)
  message(
    "Peaks above 50% threshold: ", nrow(peak_above_50), " traits in ",
    nrow(shared_peaks), " shared region(s)."
  )
} else {
  write.csv(
    data.frame(chr = integer(), pos_median_cM = numeric(), n_traits = integer(), traits = character()),
    "WatersMass_chemistry_shared_peaks_above50pct.csv",
    row.names = FALSE
  )
  message("No peaks above 50% LOD threshold; wrote empty shared-peaks file.")
}

# ---- Plots: significant traits only ----
pdf_scanone <- "WatersMass_chemistry_scanone_plots_significant.pdf"
sig_traits <- significant_10$trait
if (length(sig_traits) > 0) {
  pdf(pdf_scanone, width = 14, height = 5)
  par(mar = c(3, 3, 2, 1), mfrow = c(1, 1))
  for (tr in sig_traits) {
    meta_row <- feature_meta[feature_meta$trait == tr, , drop = FALSE]
    main_lbl <- if (nrow(meta_row)) {
      sprintf("%s (m/z=%s, RT=%s)", tr, meta_row$mz[1], meta_row$rt[1])
    } else {
      tr
    }
    plot(scanone_chemistry[[tr]], lodcolumn = 1, ylab = "LOD", main = main_lbl)
    cutoff <- perm_cutoff_by_trait[[tr]]
    if (!is.na(cutoff)) abline(h = cutoff, col = "red", lwd = 3)
  }
  dev.off()
  message("Saved: ", pdf_scanone, " (", length(sig_traits), " traits)")
} else {
  message("No significant traits at 10% FDR; skipping scanone PDF.")
}

# ---- Optional multi-QTL (10% FDR traits only; set RUN_MULTI_QTL=1) ----
if (run_multi_qtl && nrow(significant_10) > 0) {
  message("RUN_MULTI_QTL: stepwiseqtl / fitqtl / bayesint for ", nrow(significant_10), " trait(s)...")
  multi_summary <- list()
  for (tr in significant_10$trait) {
    addcovar <- if (length(covar_by_trait[[tr]]) > 0) {
      qtl::pull.pheno(Costus_prob, covar_by_trait[[tr]])
    } else {
      NULL
    }
    sc2_perm <- max(100L, n_perm)
    sc2 <- tryCatch(
      scantwo(Costus_prob, pheno.col = tr, method = "hk", addcovar = addcovar, n.perm = sc2_perm),
      error = function(e) NULL
    )
    if (is.null(sc2)) {
      warning("scantwo failed for ", tr, "; skipping multi-QTL.")
      next
    }
    pen <- calc.penalties(sc2, alpha = 0.05)
    step_fit <- tryCatch(
      stepwiseqtl(
        Costus_prob, pheno.col = tr, method = "hk", model = "normal",
        covar = addcovar, penalties = pen, max.qtl = 3,
        keeplodprofile = TRUE, verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(step_fit) || nqr(step_fit) == 0) next
    rqtl <- tryCatch(
      refineqtl(Costus_prob, pheno.col = tr, qtl = step_fit, covar = addcovar, method = "hk", maxit = 10),
      error = function(e) NULL
    )
    if (is.null(rqtl)) next
    fit_formula <- attr(step_fit, "formula")
    if (is.null(fit_formula)) {
      nq <- nqr(rqtl)
      covar_terms <- if (!is.null(addcovar)) colnames(addcovar) else character(0)
      fit_formula <- stats::as.formula(paste(
        "y ~",
        paste(c(covar_terms, paste0("Q", seq_len(nq))), collapse = " + ")
      ))
    }
    fit_out <- tryCatch(
      fitqtl(
        Costus_prob, pheno.col = tr, method = "hk", model = "normal",
        covar = addcovar, get.ests = TRUE, qtl = rqtl, formula = fit_formula, dropone = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(fit_out)) next
    bi <- tryCatch(bayesint(rqtl, expandtomarkers = TRUE), error = function(e) NULL)
    multi_summary[[tr]] <- list(step = step_fit, refined = rqtl, fit = fit_out, bayesint = bi)
  }
  if (length(multi_summary) > 0) {
    saveRDS(multi_summary, "WatersMass_chemistry_multi_QTL_models.rds")
    message("Saved WatersMass_chemistry_multi_QTL_models.rds (", length(multi_summary), " traits).")
  }
}

message("Waters mass-alignment chemistry QTL mapping complete.")
