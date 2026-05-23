# SD_secondary_hybrid_index_analysis.R
# Secondary seed-dormancy (SD) analysis adding hybrid_index to covariates.
# Compares to primary model in SD_QTL_model_output_251203 (cohorts only).
# Run from repo root: Rscript QTL/SD_secondary_hybrid_index_analysis.R

local({
  rl <- file.path(Sys.getenv("REPO_ROOT", getwd()), "R_libs")
  if (dir.exists(rl)) .libPaths(c(rl, .libPaths()))
})

suppressPackageStartupMessages({
  library(qtl)
  library(rrBLUP)
  library(sommer)
})

repo_root <- Sys.getenv("REPO_ROOT", unset = normalizePath("."))
qtl_dir <- file.path(repo_root, "QTL")
GENOPROB_RDS <- Sys.getenv(
  "GENOPROB_RDS",
  unset = "/Users/admin-jgharenc/Documents/Github/alle_vill_QTL/aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds"
)

source(file.path(qtl_dir, "fitqtl_summary_utils.R"), local = TRUE)

MAF_MIN <- 0.05
MISSING_MAX <- 0.20
G_DIAG <- 1e-6

extract_dosage_matrix <- function(cross, type = "dosage") {
  mats <- lapply(cross$geno, function(g) {
    if (type == "hard") {
      gdata <- g$data
      gdata[gdata == 1] <- 0
      gdata[gdata == 2] <- 1
      gdata[gdata == 3] <- 2
      return(gdata)
    }
    apply(g$prob, c(1, 2), function(p) p[3] + 0.5 * p[2])
  })
  do.call(cbind, mats)
}

filter_markers <- function(Z, maf_min = 0.05, miss_max = 0.20) {
  p <- colMeans(Z, na.rm = TRUE) / 2
  maf <- pmin(p, 1 - p)
  miss <- colMeans(is.na(Z))
  Z[, maf >= maf_min & miss <= miss_max, drop = FALSE]
}

build_grm <- function(Z, diag_eps = 1e-6) {
  Z <- Z[, apply(Z, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  for (j in seq_len(ncol(Z))) {
    na <- is.na(Z[, j])
    if (any(na)) Z[na, j] <- mean(Z[, j], na.rm = TRUE)
  }
  G <- A.mat(scale(Z, center = TRUE, scale = FALSE))
  G + diag(diag_eps, nrow(G))
}

extract_sigma_g_e <- function(fit) {
  sg <- fit$sigma$`u:id`[1, 1]
  se <- fit$sigma$units[1, 1]
  vc <- tryCatch(summary(fit)$varcomp, error = function(e) NULL)
  sg_se <- se_se <- NA_real_
  if (!is.null(vc) && nrow(vc) >= 2) {
    rn <- rownames(vc)
    sg_se <- vc[grepl("^u:id", rn), "VarCompSE"][1]
    se_se <- vc[grepl("^units", rn), "VarCompSE"][1]
  }
  list(sigma_g = as.numeric(sg), sigma_e = as.numeric(se), sigma_g_se = as.numeric(sg_se), sigma_e_se = as.numeric(se_se))
}

fit_grm_reml <- function(y, X, id, G) {
  dat <- data.frame(y = y, id = factor(id))
  if (!is.null(X) && ncol(X) > 0) dat <- cbind(dat, X)
  cc <- complete.cases(dat)
  dat <- dat[cc, , drop = FALSE]
  G <- G[cc, cc, drop = FALSE]
  rownames(G) <- colnames(G) <- as.character(dat$id)
  f <- if (is.null(X) || ncol(X) == 0) y ~ 1 else as.formula(paste("y ~", paste(colnames(X), collapse = " + ")))
  fit <- mmer(fixed = f, random = ~vsr(id, Gu = G), data = dat, verbose = FALSE)
  sig <- extract_sigma_g_e(fit)
  vp <- sig$sigma_g + sig$sigma_e
  h2 <- if (is.finite(vp) && vp > 0) sig$sigma_g / vp else NA_real_
  list(n = nrow(dat), sigma_g = sig$sigma_g, sigma_e = sig$sigma_e, vp_cond = vp, h2 = h2)
}

if (!file.exists(GENOPROB_RDS)) stop("GENOPROB_RDS not found: ", GENOPROB_RDS)
message("Loading: ", GENOPROB_RDS)
cross <- readRDS(GENOPROB_RDS)

coh_vars <- c("s.coh4", "s.coh7", "s.coh8", "s.coh9", "s.coh14", "s.coh15")
covar_primary <- pull.pheno(cross, coh_vars)
covar_hi <- pull.pheno(cross, c("hybrid_index", coh_vars))

# QTL positions from final SD model (2@68.6, 6@55.5)
SD.qtl <- makeqtl(cross, chr = c("2", "6"), pos = c(68.577962, 55.490804), what = "prob")
SD.rqtl <- refineqtl(
  cross, pheno.col = "log_days_to_radicle", qtl = SD.qtl,
  covar = covar_hi, method = "hk", maxit = 10, maxit.fitqtl = 1000
)

sd_formula_hi <- y ~ Q1 + Q2 + hybrid_index + s.coh4 + s.coh7 + s.coh8 + s.coh9 + s.coh14 + s.coh15
SD.fit.hi <- fitqtl(
  cross, pheno.col = "log_days_to_radicle", qtl = SD.rqtl, method = "hk",
  get.ests = TRUE, covar = covar_hi, formula = sd_formula_hi, dropone = TRUE
)

out_fitqtl <- file.path(qtl_dir, "SD_QTL_model_output_251203_with_hybrid_index")
sink(out_fitqtl)
cat("SD secondary analysis: hybrid_index added to covariates\n")
cat("Primary reference: SD_QTL_model_output_251203 (cohorts only)\n\n")
print(summary(SD.fit.hi))
sink()

message("Saved fitqtl summary: ", out_fitqtl)

# Parse for partition metrics
parsed_hi <- parse_fitqtl_summary_file(out_fitqtl)
do <- parsed_hi$drop_one
do$component <- vapply(do$term, classify_fitqtl_term, character(1))
pct_qtl <- sum(do$pctvar[do$component == "QTL"], na.rm = TRUE)
pct_hi <- sum(do$pctvar[do$term == "hybrid_index"], na.rm = TRUE)
pct_env <- sum(do$pctvar[do$component == "environment"], na.rm = TRUE)
v_qtl <- sum(do$ss_type3[do$component == "QTL"], na.rm = TRUE)
v_hi <- sum(do$ss_type3[do$term == "hybrid_index"], na.rm = TRUE)
vp_total <- parsed_hi$ss_total / (parsed_hi$n_obs - 1)

# REML h2 with hybrid_index in fixed effects
Z <- filter_markers(extract_dosage_matrix(cross, "dosage"), MAF_MIN, MISSING_MAX)
G <- build_grm(Z, G_DIAG)
phe <- as.data.frame(cross$pheno)
ids <- as.character(pull.pheno(cross, "ID"))
y <- as.numeric(phe$log_days_to_radicle)
X_hi <- as.data.frame(phe[, c("hybrid_index", coh_vars), drop = FALSE])
for (nm in coh_vars) X_hi[[nm]] <- as.numeric(X_hi[[nm]])
cc <- complete.cases(cbind(y, X_hi))
res_reml <- fit_grm_reml(y[cc], X_hi[cc, , drop = FALSE], ids[cc], G[cc, cc])

# Primary SD row for comparison
primary_path <- file.path(qtl_dir, "SD_QTL_model_output_251203")
parsed_pri <- parse_fitqtl_summary_file(primary_path)
do_pri <- parsed_pri$drop_one
do_pri$component <- vapply(do_pri$term, classify_fitqtl_term, character(1))
pct_qtl_pri <- sum(do_pri$pctvar[do_pri$component == "QTL"], na.rm = TRUE)

# AIC comparison (phenotype-only models, no QTL)
dat_lm <- phe[cc, c("log_days_to_radicle", "hybrid_index", coh_vars), drop = FALSE]
dat_lm$log_days_to_radicle <- as.numeric(dat_lm$log_days_to_radicle)
f_coh <- as.formula(paste("log_days_to_radicle ~", paste(coh_vars, collapse = " + ")))
f_hi <- as.formula(paste("log_days_to_radicle ~ hybrid_index +", paste(coh_vars, collapse = " + ")))
aic_coh <- AIC(lm(f_coh, data = dat_lm))
aic_hi <- AIC(lm(f_hi, data = dat_lm))
lr <- anova(lm(f_coh, data = dat_lm), lm(f_hi, data = dat_lm))

comparison <- data.frame(
  model = c("primary_cohorts_only", "secondary_with_hybrid_index"),
  covariates = c(
    paste(coh_vars, collapse = " + "),
    paste(c("hybrid_index", coh_vars), collapse = " + ")
  ),
  n = c(parsed_pri$n_obs, parsed_hi$n_obs),
  pct_QTL = c(pct_qtl_pri, pct_qtl),
  pct_hybrid_index = c(0, pct_hi),
  pct_environment = c(
    sum(do_pri$pctvar[do_pri$component == "environment"], na.rm = TRUE),
    pct_env
  ),
  full_model_pctvar = c(parsed_pri$full_model_pctvar, parsed_hi$full_model_pctvar),
  h2_REML = c(NA_real_, res_reml$h2),
  sigma_g = c(NA_real_, res_reml$sigma_g),
  sigma_e = c(NA_real_, res_reml$sigma_e),
  vp_cond = c(NA_real_, res_reml$vp_cond),
  frac_genetic_QTL = c(
    NA_real_,
    if (is.finite(res_reml$h2) && res_reml$h2 > 0) pct_qtl / (res_reml$h2 * 100) else NA_real_
  ),
  AIC_lm_cohorts_only = c(aic_coh, aic_coh),
  AIC_lm_with_HI = c(NA_real_, aic_hi),
  stringsAsFactors = FALSE
)

# Fill primary h2 from existing summary if present
sum_path <- file.path(qtl_dir, "genomic_heritability_QTL_traits.csv")
if (file.exists(sum_path)) {
  prev <- read.csv(sum_path, stringsAsFactors = FALSE)
  pri_h2 <- prev$h2[prev$trait == "SD"]
  if (length(pri_h2) == 1) {
    comparison$h2_REML[1] <- pri_h2
    comparison$frac_genetic_QTL[1] <- if (pri_h2 > 0) pct_qtl_pri / (pri_h2 * 100) else NA_real_
  }
}

out_cmp <- file.path(qtl_dir, "SD_hybrid_index_secondary_comparison.csv")
write.csv(comparison, out_cmp, row.names = FALSE)

# Side-by-side row for main architecture table
arch_hi <- data.frame(
  trait = "SD_with_hybrid_index",
  analysis = "secondary_sensitivity",
  n_obs = parsed_hi$n_obs,
  pct_QTL = pct_qtl,
  pct_hybrid_index = pct_hi,
  pct_F21 = 0,
  pct_environment = pct_env,
  h2_REML = res_reml$h2,
  Vg = res_reml$sigma_g,
  Ve = res_reml$sigma_e,
  vp_cond = res_reml$vp_cond,
  frac_genetic_QTL = comparison$frac_genetic_QTL[2],
  fitqtl_output_file = basename(out_fitqtl),
  stringsAsFactors = FALSE
)
out_arch <- file.path(qtl_dir, "genetic_architecture_summary_SD_with_HI.csv")
write.csv(arch_hi, out_arch, row.names = FALSE)

message("Saved: ", out_cmp)
message("Saved: ", out_arch)
print(comparison)
print(lr)
