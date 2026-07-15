# Subset scanone panel: SD, growth rate, chlorophyll reabsorption on chr 1/2/3/5/6.
# Run from repo root: Rscript QTL/plot_subset_scanone_chr12356.R

if (basename(getwd()) == "QTL") setwd("..")
repo_rlibs <- file.path(getwd(), "R_libs")
if (dir.exists(repo_rlibs)) .libPaths(c(repo_rlibs, .libPaths()))

library(qtl)

genoprob_path <- "/Users/admin-jgharenc/Documents/Github/alle_vill_QTL/aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds"
if (!file.exists(genoprob_path)) {
  stop("Genoprob RDS not found: ", genoprob_path)
}

message("Loading genoprob cross...")
Costus_prob <- readRDS(genoprob_path)

covariates_SD_select <- pull.pheno(Costus_prob, c("s.coh4", "s.coh7", "s.coh8", "s.coh9", "s.coh14", "s.coh15"))
covariates_GR_select <- pull.pheno(Costus_prob, c(
  "hybrid_index", "F21_259", "F21_238", "F21_284", "F21_232",
  "coh7", "coh10", "coh11", "coh13", "coh18", "coh19", "coh24", "coh25",
  "coh27", "coh28", "coh29", "coh30", "coh31", "coh32", "coh34"
))
covariates_Chlreabs_select <- pull.pheno(Costus_prob, c(
  "hybrid_index", "F21_238", "F21_284", "F21_232",
  "coh13", "coh16", "coh19", "coh21", "coh22", "coh24", "coh25",
  "coh27", "coh28", "coh30", "coh31", "coh32"
))

# 10% LOD FDR cutoffs (1000 permutations; Costus_QTL.qmd)
SD_10cutoff_lod <- 3.26
GR_10cutoff_lod <- 3.43
Chlreabs_10cutoff_lod <- 3.38

message("Running scanone (seed dormancy)...")
seed.dormancy.scanone.hk <- scanone(
  Costus_prob, method = "hk", pheno.col = "log_days_to_radicle", addcovar = covariates_SD_select
)
message("Running scanone (growth rate)...")
growth.rate.scanone.hk <- scanone(
  Costus_prob, method = "hk", pheno.col = "GR_cm.dy", addcovar = covariates_GR_select
)
message("Running scanone (chlorophyll reabsorption)...")
Chlreabs.scanone.hk <- scanone(
  Costus_prob, method = "hk", pheno.col = "chloro_reabsorp", addcovar = covariates_Chlreabs_select
)

subset_chr <- c(1, 2, 3, 5, 6)
subset_scanone <- function(scanone_obj, chr_keep) {
  scanone_obj[scanone_obj$chr %in% chr_keep, ]
}

scanone_subset_plots <- list(
  list(obj = subset_scanone(seed.dormancy.scanone.hk, subset_chr), cutoff = SD_10cutoff_lod, title = "seed dormancy"),
  list(obj = subset_scanone(growth.rate.scanone.hk, subset_chr), cutoff = GR_10cutoff_lod, title = "growth rate"),
  list(obj = subset_scanone(Chlreabs.scanone.hk, subset_chr), cutoff = Chlreabs_10cutoff_lod, title = "% chlorophyll reabsorption")
)

dir.create("QTL_plots", showWarnings = FALSE)
out_pdf <- "QTL_plots/sup_subset_QTL_scanone_chr12356.pdf"
out_png <- "QTL_plots/sup_subset_QTL_scanone_chr12356.png"

draw_subset <- function() {
  # Large margins + cex => compact LOD panels, prominent labels
  par(mfrow = c(3, 1), mar = c(7.5, 7, 6, 2), mgp = c(3.2, 0.9, 0))
  for (p in scanone_subset_plots) {
    plot(
      p$obj,
      lodcolumn = 1,
      ylab = "LOD score",
      main = paste0(p$title, " (chr ", paste(subset_chr, collapse = ","), ")"),
      cex.lab = 2.0,
      cex.main = 1.9,
      cex.axis = 1.65
    )
    abline(h = p$cutoff, col = "red", lwd = 2.5)
  }
}

pdf(out_pdf, width = 10, height = 10)
draw_subset()
dev.off()

png(out_png, width = 2000, height = 2000, res = 300)
draw_subset()
dev.off()

final_figures_dir <- file.path(
  Sys.getenv("HOME"),
  "Library/CloudStorage/GoogleDrive-jgharencar@ucdavis.edu",
  "My Drive/Costus/PopGen/Figures&Tables/Final_figures"
)
if (dir.exists(final_figures_dir)) {
  file.copy(out_pdf, file.path(final_figures_dir, basename(out_pdf)), overwrite = TRUE)
  file.copy(out_png, file.path(final_figures_dir, basename(out_png)), overwrite = TRUE)
  message("Copied to ", final_figures_dir)
}

message("Wrote ", out_pdf)
message("Wrote ", out_png)
