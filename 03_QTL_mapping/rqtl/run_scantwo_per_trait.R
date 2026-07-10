#!/usr/bin/env Rscript
# scantwo permutations for stepwiseqtl penalty calculation (Farm or local).
# Usage: TRAIT=SD Rscript run_scantwo_per_trait.R

suppressPackageStartupMessages(library(qtl))

trait_config <- list(
  SD = list(pheno = "log_days_to_radicle",
    covar_cols = c("s.coh4", "s.coh7", "s.coh8", "s.coh9", "s.coh14", "s.coh15"),
    outfile = "SD.scantwo.hk.perm.rds"),
  GR = list(pheno = "GR_cm.dy",
    covar_cols = c("hybrid_index", "F21_259", "F21_238", "F21_284", "F21_232",
      "coh7", "coh10", "coh11", "coh13", "coh18", "coh19", "coh24",
      "coh25", "coh27", "coh28", "coh29", "coh30", "coh31", "coh32", "coh34"),
    outfile = "GR.scantwo.hk.perm.rds"),
  tough = list(pheno = "ave_tough",
    covar_cols = c("hybrid_index", "F21_232", "F21_238", "F21_259", "F21_284",
      "coh10", "coh11", "coh16", "coh19", "coh22", "coh24", "coh27",
      "coh28", "coh29", "coh30", "coh32"),
    outfile = "tough.scantwo.hk.perm.rds"),
  Nreabs = list(pheno = "N_reabsorp",
    covar_cols = c("hybrid_index", "F21_238", "F21_259", "F21_232", "F21_284",
      "coh10", "coh16", "coh18", "coh21", "coh22", "coh24", "coh25",
      "coh27", "coh29", "coh30", "coh32", "coh34"),
    outfile = "Nreabs.scantwo.hk.perm.rds"),
  Chlreabs = list(pheno = "chloro_reabsorp",
    covar_cols = c("hybrid_index", "F21_238", "F21_284", "F21_232",
      "coh13", "coh16", "coh19", "coh21", "coh22", "coh24", "coh25",
      "coh27", "coh28", "coh30", "coh31", "coh32"),
    outfile = "Chlreabs.scantwo.hk.perm.rds")
)

trait <- Sys.getenv("TRAIT", "")
if (!nzchar(trait) || !trait %in% names(trait_config)) {
  stop("Set TRAIT to one of: ", paste(names(trait_config), collapse = ", "))
}
cfg <- trait_config[[trait]]
for (p in c(
  file.path("scripts", "project_paths.R"),
  file.path("..", "..", "scripts", "project_paths.R"),
  file.path("03_QTL_mapping", "rqtl", "..", "..", "scripts", "project_paths.R")
)) {
  if (file.exists(p)) { source(p); break }
}
paths <- if (exists("project_paths")) project_paths() else list()
workdir <- Sys.getenv("WORKDIR", if (!is.null(paths$rqtl_dir)) paths$rqtl_dir else getwd())
setwd(workdir)

cross_rds <- Sys.getenv("CROSS_RDS",
  if (!is.null(paths$genoprob_rds) && nzchar(paths$genoprob_rds)) {
    paths$genoprob_rds
  } else {
    "aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_least_dups_est.rf.rds"
  })
n_perm <- as.integer(Sys.getenv("N_PERM", "1000"))
n_cores <- as.integer(Sys.getenv("N_CORES", "8"))
outdir <- Sys.getenv("OUTDIR", if (!is.null(paths$scantwo_perm_dir)) paths$scantwo_perm_dir else "scantwo.perms")

for (p in c("log_all_warnings.R", file.path("QTL", "log_all_warnings.R"),
            file.path("03_QTL_mapping/rqtl/log_all_warnings.R"))) {
  if (file.exists(p)) { source(p, local = FALSE); start_warning_log(paste0(trait, "_scantwo_warnings.txt")); break }
}
if (!file.exists(cross_rds)) stop("CROSS_RDS not found: ", cross_rds)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

Costus <- readRDS(cross_rds)
covariates <- pull.pheno(Costus, cfg$covar_cols)
Costus <- calc.genoprob(Costus, step = 5, error.prob = 0.01)
perm <- scantwo(Costus, method = "hk", pheno.col = cfg$pheno, n.perm = n_perm,
  addcovar = covariates, n.cluster = n_cores, verbose = TRUE)
outfile <- file.path(outdir, cfg$outfile)
saveRDS(perm, outfile)
message("Saved: ", outfile)
print(calc.penalties(perm, alpha = 0.05))
