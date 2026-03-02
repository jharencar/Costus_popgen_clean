# Preliminary QTL setup for UV chemistry data
# Run from project root or with working directory = QTL/ (script will switch to QTL/ if needed).
# Setup mirrors first 5 chunks of Costus_QTL.qmd, without the secondary
# duplicate-loci removal (no scantwo for this preliminary analysis).
# Lab_ID in prelim_UV_chemistry_transposed_w_LabID.csv matches ID in
# F3phe_onehotcovar_newHI_20250925.csv and the genotype file.

library(dplyr)
library(qtl)
if (requireNamespace("qtl2ggplot", quietly = TRUE)) library(qtl2ggplot)  # optional; script uses base plot()

# Permutations per chemistry trait for LOD threshold (10% alpha); set to 0 to skip perms and threshold line
n_perm <- 100

# Use QTL/ as working directory if data files are there (allows running from project root)
phe_file <- "F3phe_onehotcovar_newHI_20250925.csv"
if (!file.exists(phe_file) && file.exists(file.path("QTL", phe_file))) {
  setwd("QTL")
  message("Set working directory to QTL/")
}

# ---- Load cross (skip time-intensive setup: merge, filter, est.rf) ----
# Uncomment and run sections 1-5 below once to create aHMM_prelim_UV_chemistry_est.rf.rds.
out_rds <- "aHMM_prelim_UV_chemistry_est.rf.rds"
if (file.exists(out_rds)) {
  rds_path <- out_rds
} else if (file.exists(file.path("QTL", out_rds))) {
  setwd("QTL")
  message("Set working directory to QTL/")
  rds_path <- out_rds
} else {
  stop("Run the commented setup (sections 1-5) once to create ", out_rds, " in QTL/")
}
Costus <- readRDS(rds_path)
chem_trait_names <- names(Costus$pheno)[grepl("^V[0-9]+$", names(Costus$pheno))]
message("Loaded ", out_rds, "; ", nind(Costus), " individuals, ", length(chem_trait_names), " chemistry traits.")

# ## ---- 1. Merge chemistry into phenotype ----
# ## Keep all F3phe individuals; add chemistry columns (NA where no chemistry).
# phe <- read.csv(phe_file, check.names = FALSE)
# chem <- read.csv("prelim_UV_chemistry_transposed_w_LabID.csv", check.names = FALSE)
# chem_only <- chem %>% select(-field_ID)
# phe_with_chem <- phe %>% left_join(chem_only, by = c("ID" = "Lab_ID"))
# phefile_merged <- "F3phe_onehotcovar_newHI_with_prelim_UV_chemistry.csv"
# write.csv(phe_with_chem, phefile_merged, row.names = FALSE)
#
# ## ---- 2. Import cross ----
# genfile <- "f3gen_250924_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001.csv"
# OG.cross <- read.cross("csvs", ".", genfile = genfile, phefile = phefile_merged, estimate.map = FALSE)
# OG.cross <- jittermap(OG.cross)
#
# ## ---- 3. Filter ----
# n_mar <- nmar(OG.cross)
# OG.cross <- subset(OG.cross, ind = (ntyped(OG.cross) > 0.8 * n_mar))
# n_ind <- nind(OG.cross)
# nt.bymar <- ntyped(OG.cross, "mar")
# todrop <- names(nt.bymar[nt.bymar < 0.8 * n_ind])
# filtered.cross <- drop.markers(OG.cross, todrop)
# dup <- findDupMarkers(filtered.cross, exact.only = TRUE, adjacent.only = TRUE)
# Costus_w.o_dups <- drop.markers(filtered.cross, unlist(dup))
#
# ## ---- 4. Duplicate individuals ----
# cg <- comparegeno(Costus_w.o_dups)
# wh <- which(cg > 0.95, arr.ind = TRUE)
# wh <- wh[wh[, 1] < wh[, 2], ]
# g1 <- pull.geno(Costus_w.o_dups, chr = 1)
# g2 <- pull.geno(Costus_w.o_dups, chr = 2)
# g3 <- pull.geno(Costus_w.o_dups, chr = 3)
# g4 <- pull.geno(Costus_w.o_dups, chr = 4)
# g5 <- pull.geno(Costus_w.o_dups, chr = 5)
# g6 <- pull.geno(Costus_w.o_dups, chr = 6)
# g7 <- pull.geno(Costus_w.o_dups, chr = 7)
# g8 <- pull.geno(Costus_w.o_dups, chr = 8)
# g9 <- pull.geno(Costus_w.o_dups, chr = 9)
# if (nrow(wh) > 0) {
#   for (i in seq_len(nrow(wh))) {
#     tozero <- !is.na(g1[wh[i, 1], ]) & !is.na(g1[wh[i, 2], ]) & g1[wh[i, 1], ] != g1[wh[i, 2], ]
#     Costus_w.o_dups$geno[[1]]$data[wh[i, 1], tozero] <- NA
#   }
#   for (i in seq_len(nrow(wh))) { tozero <- !is.na(g2[wh[i, 1], ]) & !is.na(g2[wh[i, 2], ]) & g2[wh[i, 1], ] != g2[wh[i, 2], ]; Costus_w.o_dups$geno[[2]]$data[wh[i, 1], tozero] <- NA }
#   for (i in seq_len(nrow(wh))) { tozero <- !is.na(g3[wh[i, 1], ]) & !is.na(g3[wh[i, 2], ]) & g3[wh[i, 1], ] != g3[wh[i, 2], ]; Costus_w.o_dups$geno[[3]]$data[wh[i, 1], tozero] <- NA }
#   for (i in seq_len(nrow(wh))) { tozero <- !is.na(g4[wh[i, 1], ]) & !is.na(g4[wh[i, 2], ]) & g4[wh[i, 1], ] != g4[wh[i, 2], ]; Costus_w.o_dups$geno[[4]]$data[wh[i, 1], tozero] <- NA }
#   for (i in seq_len(nrow(wh))) { tozero <- !is.na(g5[wh[i, 1], ]) & !is.na(g5[wh[i, 2], ]) & g5[wh[i, 1], ] != g5[wh[i, 2], ]; Costus_w.o_dups$geno[[5]]$data[wh[i, 1], tozero] <- NA }
#   for (i in seq_len(nrow(wh))) { tozero <- !is.na(g6[wh[i, 1], ]) & !is.na(g6[wh[i, 2], ]) & g6[wh[i, 1], ] != g6[wh[i, 2], ]; Costus_w.o_dups$geno[[6]]$data[wh[i, 1], tozero] <- NA }
#   for (i in seq_len(nrow(wh))) { tozero <- !is.na(g7[wh[i, 1], ]) & !is.na(g7[wh[i, 2], ]) & g7[wh[i, 1], ] != g7[wh[i, 2], ]; Costus_w.o_dups$geno[[7]]$data[wh[i, 1], tozero] <- NA }
#   for (i in seq_len(nrow(wh))) { tozero <- !is.na(g8[wh[i, 1], ]) & !is.na(g8[wh[i, 2], ]) & g8[wh[i, 1], ] != g8[wh[i, 2], ]; Costus_w.o_dups$geno[[8]]$data[wh[i, 1], tozero] <- NA }
#   for (i in seq_len(nrow(wh))) { tozero <- !is.na(g9[wh[i, 1], ]) & !is.na(g9[wh[i, 2], ]) & g9[wh[i, 1], ] != g9[wh[i, 2], ]; Costus_w.o_dups$geno[[9]]$data[wh[i, 1], tozero] <- NA }
#   Costus.clean <- subset(Costus_w.o_dups, ind = -wh[, 2])
# } else { Costus.clean <- Costus_w.o_dups }
#
# ## ---- 5. est.rf and save ----
# Costus <- est.rf(Costus.clean)
# saveRDS(Costus, file = out_rds)

# ---- 6. Preliminary QTL: model selection and scanone ----
# Same model-selection logic as Costus_QTL.qmd (AIC stepwise) but streamlined.
# Covariates for each chemistry peak are chosen via step(lm(y ~ candidates)).

Costus_prob <- calc.genoprob(Costus)

# Candidate covariates: hybrid_index, family one-hots, growth-cohort one-hots with >= min_n
min_covar_n <- 10
pn <- names(Costus$pheno)
fam_onehot <- c("F21_238", "F21_259", "F21_241", "F21_284", "F21_232")
coh_cols <- pn[grepl("^coh[0-9]+$", pn)]
coh_n <- colSums(Costus$pheno[, coh_cols, drop = FALSE], na.rm = TRUE)
coh_ok <- coh_cols[coh_n >= min_covar_n]
candidate_covars <- c("hybrid_index", fam_onehot[fam_onehot %in% pn], coh_ok)
candidate_covars <- candidate_covars[candidate_covars %in% pn]
message("Candidate covariates (n=", length(candidate_covars), "): ", paste(head(candidate_covars, 8), collapse = ", "), " ...")
message("Chemistry traits to map: ", length(chem_trait_names))

#' Run AIC stepwise model selection for one trait; return addcovar matrix for scanone.
#' Drops rows with NA in y or covariates; uses same candidate set as main qmd.
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

#' Drop one reference from each one-hot group so addcovar is full rank (avoids DGELS error in scanone).
fam_onehot_all <- c("F21_238", "F21_259", "F21_241", "F21_284", "F21_232")
ensure_full_rank_covariate_names <- function(names) {
  out <- names
  fam_in <- fam_onehot_all[fam_onehot_all %in% out]
  if (length(fam_in) >= 5) out <- setdiff(out, fam_onehot_all[5])  # drop F21_232 as reference
  coh_in <- out[grepl("^coh[0-9]+$", out)]
  if (length(coh_in) >= 2) out <- setdiff(out, coh_in[length(coh_in)])  # drop one cohort as reference
  out
}

# Run scanone for each chemistry trait with its own selected covariates
scanone_chemistry <- list()
covar_by_trait <- list()
failed_traits <- character(0)
n_ind <- nind(Costus_prob)
for (trait in chem_trait_names) {
  covar <- select_covariates_for_trait(Costus, trait, candidate_covars)
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
  warning("scanone failed for ", length(failed_traits), " trait(s) (e.g. too few non-NA or singular design): ", paste(head(failed_traits, 5), collapse = ", "), if (length(failed_traits) > 5) " ..." else "")
}
message("Completed scanone for ", length(scanone_chemistry), " chemistry traits.", if (length(failed_traits) > 0) paste0(" Skipped ", length(failed_traits), ".") else "")

# Summary: max LOD and position per trait (scanone stores LOD in column named after phenotype)
if (length(scanone_chemistry) == 0) {
  stop("No chemistry traits could be mapped (all scanone runs failed). Check phenotype NAs or try fewer covariates.")
}
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
message("Top 10 peaks (by LOD):")
print(head(summ, 10))

# Save results
saveRDS(scanone_chemistry, "prelim_UV_chemistry_scanone_by_trait.rds")
saveRDS(covar_by_trait, "prelim_UV_chemistry_covariates_by_trait.rds")
write.csv(summ, "prelim_UV_chemistry_scanone_summary.csv", row.names = FALSE)
message("Saved: prelim_UV_chemistry_scanone_by_trait.rds, prelim_UV_chemistry_covariates_by_trait.rds, prelim_UV_chemistry_scanone_summary.csv")

# Permutations per trait for 10% and 50% LOD thresholds (as in main qmd)
perm_cutoff_by_trait <- setNames(rep(NA_real_, length(scanone_chemistry)), names(scanone_chemistry))
perm_cutoff_50_by_trait <- setNames(rep(NA_real_, length(scanone_chemistry)), names(scanone_chemistry))
if (n_perm > 0) {
  message("Running ", n_perm, " permutations per trait for LOD thresholds (10% and 50%)...")
  for (tr in names(scanone_chemistry)) {
    addcovar <- if (length(covar_by_trait[[tr]]) > 0) qtl::pull.pheno(Costus, covar_by_trait[[tr]]) else NULL
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
  saveRDS(perm_cutoff_by_trait, "prelim_UV_chemistry_perm_cutoffs_10pct.rds")
  saveRDS(perm_cutoff_50_by_trait, "prelim_UV_chemistry_perm_cutoffs_50pct.rds")
  message("Permutation cutoffs saved (10% and 50%).")
}

# Identify chemicals with shared or close peaks above 50% LOD threshold
peak_above_50 <- summ
peak_above_50$cutoff_50 <- perm_cutoff_50_by_trait[peak_above_50$trait]
peak_above_50 <- peak_above_50[peak_above_50$lod > peak_above_50$cutoff_50 & !is.na(peak_above_50$cutoff_50), ]
peak_above_50 <- peak_above_50[order(peak_above_50$chr, peak_above_50$pos), ]
cM_window <- 15  # traits within this many cM on same chr considered "shared" peak
if (nrow(peak_above_50) > 0) {
  peak_above_50$group_id <- NA_integer_
  grp <- 0
  for (i in seq_len(nrow(peak_above_50))) {
    if (!is.na(peak_above_50$group_id[i])) next
    grp <- grp + 1
    peak_above_50$group_id[i] <- grp
    chr_i <- peak_above_50$chr[i]
    pos_i <- peak_above_50$pos[i]
    # all on same chr within cM_window get same group (then merge overlapping in next pass)
    same_region <- peak_above_50$chr == chr_i & abs(peak_above_50$pos - pos_i) <= cM_window
    peak_above_50$group_id[same_region] <- grp
  }
  # merge groups: if any member of group A is within cM_window of any member of group B, same group
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
  shared_peaks$pos_median_cM <- pos_med$pos[match(paste(shared_peaks$chr, shared_peaks$group_id), paste(pos_med$chr, pos_med$group_id))]
  shared_peaks <- shared_peaks[order(-shared_peaks$n_traits), c("chr", "pos_median_cM", "n_traits", "traits")]
  write.csv(shared_peaks, "prelim_UV_chemistry_shared_peaks_above50pct.csv", row.names = FALSE)
  message("Chemicals with peak above 50% threshold: ", nrow(peak_above_50), " traits in ", nrow(shared_peaks), " shared/close peak region(s).")
  message("Shared peaks (above 50% LOD, within ", cM_window, " cM) written to prelim_UV_chemistry_shared_peaks_above50pct.csv")
  if (nrow(shared_peaks) > 0) print(shared_peaks)
} else {
  message("No peaks above 50% LOD threshold (or no permutation run); skipping shared-peak summary.")
}

# Preliminary scanone plots: single column (one plot per page) for wider view of QTL along chromosomes
pdf_scanone <- "prelim_UV_chemistry_scanone_plots.pdf"
pdf(pdf_scanone, width = 14, height = 5)
par(mar = c(3, 3, 2, 1), mfrow = c(1, 1))
traits_ordered <- summ$trait
for (tr in traits_ordered) {
  plot(scanone_chemistry[[tr]], lodcolumn = 1, ylab = "LOD", main = tr)
  cutoff <- perm_cutoff_by_trait[[tr]]
  if (!is.na(cutoff)) abline(h = cutoff, col = "red", lwd = 3)
}
dev.off()
message("Saved: ", pdf_scanone)

message("Preliminary UV chemistry QTL setup and scanone complete.")
