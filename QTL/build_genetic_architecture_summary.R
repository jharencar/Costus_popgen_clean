# build_genetic_architecture_summary.R
# Merge REML h2, fitqtl PVE, and fixed-effect partitions.
# Run after genomic_heritability_QTL_traits.R and genetic_variance_partition.R.

suppressPackageStartupMessages(library(dplyr))

repo_root <- Sys.getenv("REPO_ROOT", unset = normalizePath("."))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) repo_root <- normalizePath(args[1L])
qtl_dir <- file.path(repo_root, "QTL")

h2_path <- file.path(qtl_dir, "genomic_heritability_QTL_traits.csv")
part_path <- file.path(qtl_dir, "genetic_variance_partition_summary.csv")
sens_path <- file.path(qtl_dir, "genomic_heritability_GRM_sensitivity.csv")

for (p in c(h2_path, part_path)) {
  if (!file.exists(p)) stop("Missing required file: ", p, "\nRun genomic_heritability_QTL_traits.R and genetic_variance_partition.R first.")
}

h2 <- read.csv(h2_path, stringsAsFactors = FALSE)
part <- read.csv(part_path, stringsAsFactors = FALSE)
sens <- if (file.exists(sens_path)) read.csv(sens_path, stringsAsFactors = FALSE) else NULL

merged <- part %>%
  select(
    trait, n_obs, pct_QTL, pct_hybrid_index, pct_F21, pct_environment,
    v_QTL, v_genetic_fixed, h2_fixed_approx
  ) %>%
  left_join(
    h2 %>% select(trait, n = n, sigma_g, sigma_e, vp_cond, h2 = h2, n_markers_grm, family_icc),
    by = "trait"
  )

if (!is.null(sens)) {
  merged <- merged %>% left_join(sens %>% select(trait, h2_dosage, h2_hard), by = "trait")
}

merged <- merged %>%
  mutate(
    frac_genetic_QTL = ifelse(is.finite(h2) & h2 > 0, pct_QTL / (h2 * 100), NA_real_),
    frac_genetic_QTL_note = "PVE_QTL (%) / h2_REML; h2 = sigma_g^2/(sigma_g^2+sigma_e^2) conditional on QTL covariates"
  )

out_path <- file.path(qtl_dir, "genetic_architecture_summary.csv")
meta_path <- file.path(qtl_dir, "genetic_architecture_summary_column_metadata.csv")
write.csv(merged, out_path, row.names = FALSE)
message("Saved: ", out_path)
if (file.exists(meta_path)) {
  message("Column metadata: ", meta_path)
  message("Data dictionary: ", file.path(qtl_dir, "genetic_architecture_summary_DATA_DICTIONARY.md"))
}
print(merged)
