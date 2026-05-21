# genetic_variance_partition.R
# Supplementary fixed-effect variance partition from fitqtl drop-one tables.
# Run from repo root: Rscript QTL/genetic_variance_partition.R

suppressPackageStartupMessages(library(dplyr))

repo_root <- Sys.getenv("REPO_ROOT", unset = normalizePath("."))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) repo_root <- normalizePath(args[1L])
qtl_dir <- file.path(repo_root, "QTL")

source(file.path(qtl_dir, "fitqtl_summary_utils.R"), local = TRUE)

model_files <- list.files(qtl_dir, pattern = glob2rx("*_QTL_model_output_251203"), full.names = TRUE)
if (length(model_files) == 0) stop("No model output files in ", qtl_dir)

rows <- list()
drop_rows <- list()

for (path in model_files) {
  trait <- sub("_QTL_model_output_251203$", "", basename(path))
  p <- parse_fitqtl_summary_file(path)
  vp_total <- if (is.finite(p$ss_total) && p$n_obs > 1) p$ss_total / (p$n_obs - 1) else NA_real_

  do <- p$drop_one
  do$component <- vapply(do$term, classify_fitqtl_term, character(1))
  do$trait <- trait
  do$v_term <- (do$pctvar / 100) * p$ss_total
  drop_rows[[trait]] <- do

  sums <- do %>%
    filter(!is.na(component)) %>%
    group_by(component) %>%
    summarise(
      pct_var = sum(pctvar, na.rm = TRUE),
      v_term = sum(v_term, na.rm = TRUE),
      .groups = "drop"
    )

  pct_qtl <- sum(do$pctvar[do$component == "QTL"], na.rm = TRUE)
  pct_hi <- sum(do$pctvar[do$term == "hybrid_index"], na.rm = TRUE)
  pct_f21 <- sum(do$pctvar[do$component == "other_genotypic" & do$term != "hybrid_index"], na.rm = TRUE)
  pct_env <- sum(do$pctvar[do$component == "environment"], na.rm = TRUE)

  v_qtl <- sum(do$v_term[do$component == "QTL"], na.rm = TRUE)
  v_hi <- sum(do$v_term[do$term == "hybrid_index"], na.rm = TRUE)
  v_f21 <- sum(do$v_term[do$component == "other_genotypic" & do$term != "hybrid_index"], na.rm = TRUE)
  v_env <- sum(do$v_term[do$component == "environment"], na.rm = TRUE)
  v_genetic_fixed <- v_qtl + v_hi + v_f21

  rows[[trait]] <- data.frame(
    trait = trait,
    n_obs = p$n_obs,
    ss_total = p$ss_total,
    vp_total = vp_total,
    full_model_pctvar = p$full_model_pctvar,
    pct_unexplained = 100 - p$full_model_pctvar,
    pct_QTL = pct_qtl,
    pct_hybrid_index = pct_hi,
    pct_F21 = pct_f21,
    pct_environment = pct_env,
    v_QTL = v_qtl,
    v_hybrid_index = v_hi,
    v_F21 = v_f21,
    v_environment = v_env,
    v_genetic_fixed = v_genetic_fixed,
    h2_fixed_approx = if (is.finite(vp_total) && vp_total > 0) v_genetic_fixed / vp_total else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(rows)
drop_df <- bind_rows(drop_rows)

out_sum <- file.path(qtl_dir, "genetic_variance_partition_summary.csv")
out_drop <- file.path(qtl_dir, "genetic_variance_partition_dropone.csv")
write.csv(summary_df, out_sum, row.names = FALSE)
write.csv(drop_df, out_drop, row.names = FALSE)
message("Saved: ", out_sum)
message("Saved: ", out_drop)
print(summary_df)
