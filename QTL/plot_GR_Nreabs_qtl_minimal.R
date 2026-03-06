# Minimal script to produce Growth rate and Nreabs QTL plots.
# Run from repo root (Costus_popgen_clean). No permutations; uses saved LOD cutoffs.

library(qtl)

# If run from QTL/, switch to repo root
if (basename(getwd()) == "QTL") setwd("..")

# Genotype-probability cross (same object as in Costus_QTL.qmd)
Costus_prob <- readRDS(file = "/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds")

# Covariates (same as in Costus_QTL.qmd)
covariates_GR_select <- pull.pheno(Costus_prob, c(
  "hybrid_index", "F21_259", "F21_238", "F21_284", "F21_232",
  "coh7", "coh10", "coh11", "coh13", "coh18", "coh19", "coh24", "coh25",
  "coh27", "coh28", "coh29", "coh30", "coh31", "coh32", "coh34"
))
covariates_Nreabs_select <- pull.pheno(Costus_prob, c(
  "hybrid_index", "F21_238", "F21_259", "F21_232", "F21_284",
  "coh10", "coh16", "coh18", "coh21", "coh22", "coh24", "coh25",
  "coh27", "coh29", "coh30", "coh32", "coh34"
))

# Scanone (no permutations)
growth.rate.scanone.hk <- scanone(Costus_prob, method = "hk", pheno.col = "GR_cm.dy",
                                  addcovar = covariates_GR_select)
Nreabs.scanone.hk <- scanone(Costus_prob, method = "hk", pheno.col = "N_reabsorp",
                             addcovar = covariates_Nreabs_select)

# Saved 10% permutation cutoffs (from Costus_QTL.qmd)
GR_10cutoff_lod    <- 3.43
Nreabs_10cutoff_lod <- 3.42

# Output directory (same as QMD)
dir.create("./QTL_plots", showWarnings = FALSE)

# --- Growth rate plot ---
pdf("./QTL_plots/growth.rate-qtl_251016.pdf", 20, 5)
par(mar = c(5, 5, 4, 2), cex.axis = 2)
plot(growth.rate.scanone.hk, lodcolumn = 1,
     ylab = "LOD score",
     main = "Growth rate",
     cex.lab = 2.5,
     cex.main = 2.5)
abline(h = GR_10cutoff_lod, col = "red", lwd = 3)
dev.off()

# --- Nreabs plot ---
pdf("./QTL_plots/Nreabs-qtl.pdf", 20, 5)
par(mar = c(5, 5, 4, 2), cex.axis = 2)
plot(Nreabs.scanone.hk, lodcolumn = 1,
     ylab = "LOD score",
     main = "% leaf nitrogen reabsorption during drought",
     cex.lab = 2.5,
     cex.main = 2.5)
abline(h = Nreabs_10cutoff_lod, col = "red", lwd = 3)
dev.off()
