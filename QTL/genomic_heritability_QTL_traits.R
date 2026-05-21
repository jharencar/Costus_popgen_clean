# genomic_heritability_QTL_traits.R
# Narrow-sense genomic heritability (GRM + REML via sommer) for QTL mapping traits.
# Run from repo root: Rscript QTL/genomic_heritability_QTL_traits.R

local({
  rp <- Sys.getenv("R_LIBS_PROJECT", "")
  if (nzchar(rp) && dir.exists(rp)) .libPaths(c(rp, .libPaths()))
  rl <- file.path(
    Sys.getenv("REPO_ROOT", getwd()),
    "R_libs"
  )
  if (dir.exists(rl)) .libPaths(c(rl, .libPaths()))
})

suppressPackageStartupMessages({
  library(qtl)
  library(rrBLUP)
  library(sommer)
  library(lme4)
})

# ---- paths (edit GENOPROB_RDS on another machine) ----
args <- commandArgs(trailingOnly = TRUE)
repo_root <- Sys.getenv("REPO_ROOT", unset = normalizePath("."))
if (length(args) >= 1L) repo_root <- normalizePath(args[1L])

qtl_dir <- file.path(repo_root, "QTL")
GENOPROB_RDS <- Sys.getenv(
  "GENOPROB_RDS",
  unset = "/Users/admin-jgharenc/Documents/Github/alle_vill_QTL/aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds"
)
MODEL_PATTERN <- "*_QTL_model_output_251203"
MAF_MIN <- 0.05
MISSING_MAX <- 0.20
G_DIAG <- 1e-6

source(file.path(qtl_dir, "fitqtl_summary_utils.R"), local = TRUE)

if (!file.exists(GENOPROB_RDS)) {
  stop("GENOPROB_RDS not found: ", GENOPROB_RDS)
}

message("Loading genoprob cross: ", GENOPROB_RDS)
cross <- readRDS(GENOPROB_RDS)

#' Dosage matrix from genotype probabilities (AA=1, AB=2, BB=3 -> 0,1,2).
extract_dosage_matrix <- function(cross, type = c("dosage", "hard")) {
  type <- match.arg(type)
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
  keep <- maf >= maf_min & miss <= miss_max
  Z[, keep, drop = FALSE]
}

build_grm <- function(Z, diag_eps = 1e-6) {
  Z <- Z[, apply(Z, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  for (j in seq_len(ncol(Z))) {
    na <- is.na(Z[, j])
    if (any(na)) Z[na, j] <- mean(Z[, j], na.rm = TRUE)
  }
  # Z: individuals x markers (rrBLUP::A.mat expects n x p)
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
  list(
    sigma_g = as.numeric(sg),
    sigma_e = as.numeric(se),
    sigma_g_se = as.numeric(sg_se),
    sigma_e_se = as.numeric(se_se)
  )
}

fit_grm_reml <- function(y, X, id, G) {
  dat <- data.frame(y = y, id = factor(id))
  if (!is.null(X) && ncol(X) > 0) dat <- cbind(dat, X)
  cc <- complete.cases(dat)
  dat <- dat[cc, , drop = FALSE]
  G <- G[cc, cc, drop = FALSE]
  rownames(G) <- colnames(G) <- as.character(dat$id)
  if (is.null(X) || ncol(X) == 0) {
    f <- y ~ 1
  } else {
    f <- as.formula(paste("y ~", paste(colnames(X), collapse = " + ")))
  }
  fit <- mmer(
    fixed = f,
    random = ~vsr(id, Gu = G),
    data = dat,
    verbose = FALSE
  )
  sig <- extract_sigma_g_e(fit)
  vp <- sig$sigma_g + sig$sigma_e
  h2 <- if (is.finite(vp) && vp > 0) sig$sigma_g / vp else NA_real_
  # Delta-method SE for h2 = sg / (sg + se)
  h2_se <- NA_real_
  if (is.finite(h2) && is.finite(sig$sigma_g_se) && is.finite(sig$sigma_e_se)) {
    dh2_dg <- sig$sigma_e / (vp^2)
    dh2_de <- -sig$sigma_g / (vp^2)
    h2_se <- sqrt((dh2_dg * sig$sigma_g_se)^2 + (dh2_de * sig$sigma_e_se)^2)
  }
  list(
    fit = fit,
    n = nrow(dat),
    sigma_g = sig$sigma_g,
    sigma_e = sig$sigma_e,
    sigma_g_se = sig$sigma_g_se,
    sigma_e_se = sig$sigma_e_se,
    vp_cond = vp,
    h2 = h2,
    h2_se = h2_se
  )
}

fit_family_icc <- function(y, X, family) {
  dat <- data.frame(y = y, family = factor(family))
  if (!is.null(X) && ncol(X) > 0) dat <- cbind(dat, X)
  cc <- complete.cases(dat)
  dat <- dat[cc, , drop = FALSE]
  f <- if (ncol(X) > 0) {
    as.formula(paste("y ~", paste(colnames(X), collapse = " + "), "+ (1|family)"))
  } else {
    y ~ (1 | family)
  }
  fit <- tryCatch(lmer(f, data = dat, REML = TRUE), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(n = nrow(dat), var_family = NA_real_, var_resid = NA_real_, icc = NA_real_, note = fit$message))
  }
  vc <- as.data.frame(VarCorr(fit))
  vf <- vc$vcov[vc$grp == "family"]
  vr <- vc$vcov[vc$grp == "Residual"]
  icc <- vf / (vf + vr)
  list(n = nrow(dat), var_family = vf, var_resid = vr, icc = icc, note = NA_character_)
}

# ---- GRM (dosage and hard-call) ----
message("Building dosage GRM (", MAF_MIN, " MAF, ", MISSING_MAX * 100, "% max missing)...")
Z_dos <- filter_markers(extract_dosage_matrix(cross, "dosage"), MAF_MIN, MISSING_MAX)
G_dos <- build_grm(Z_dos, G_DIAG)
message("Markers in GRM: ", ncol(Z_dos))

message("Building hard-call GRM for sensitivity...")
Z_hard <- filter_markers(extract_dosage_matrix(cross, "hard"), MAF_MIN, MISSING_MAX)
G_hard <- build_grm(Z_hard, G_DIAG)

phe <- as.data.frame(cross$pheno)
ids <- as.character(pull.pheno(cross, "ID"))
family <- as.character(pull.pheno(cross, "family"))

model_files <- list.files(qtl_dir, pattern = glob2rx(MODEL_PATTERN), full.names = TRUE)
if (length(model_files) == 0) stop("No model output files in ", qtl_dir)

audit_rows <- list()
result_rows <- list()
sens_rows <- list()

for (path in model_files) {
  trait <- sub("_QTL_model_output_251203$", "", basename(path))
  cfg <- trait_pheno_config(trait)
  parsed <- parse_fitqtl_summary_file(path)
  covs <- covariates_from_formula_rhs(parsed$formula_rhs)

  y <- as.numeric(phe[[cfg$pheno]])
  X <- as.data.frame(phe[, covs, drop = FALSE], stringsAsFactors = FALSE)
  for (nm in names(X)) {
    if (grepl("^coh|^s\\.coh", nm)) X[[nm]] <- as.numeric(X[[nm]])
  }

  cc <- complete.cases(cbind(y, X))
  n_cc <- sum(cc)
  audit_rows[[trait]] <- data.frame(
    trait = trait,
    pheno = cfg$pheno,
    n_fitqtl = parsed$n_obs,
    n_complete = n_cc,
    n_match = n_cc == parsed$n_obs,
    stringsAsFactors = FALSE
  )
  if (!n_cc == parsed$n_obs) {
    warning(trait, ": complete cases (", n_cc, ") != fitqtl n (", parsed$n_obs, ")")
  }

  y_sub <- y[cc]
  X_sub <- X[cc, , drop = FALSE]
  id_sub <- ids[cc]
  fam_sub <- family[cc]

  res_dos <- fit_grm_reml(y_sub, X_sub, id_sub, G_dos[cc, cc])
  res_hard <- fit_grm_reml(y_sub, X_sub, id_sub, G_hard[cc, cc])
  fam_chk <- fit_family_icc(y_sub, X_sub, fam_sub)

  result_rows[[trait]] <- data.frame(
    trait = trait,
    pheno = cfg$pheno,
    grm_type = "dosage",
    n = res_dos$n,
    n_markers_grm = ncol(Z_dos),
    sigma_g = res_dos$sigma_g,
    sigma_e = res_dos$sigma_e,
    sigma_g_se = res_dos$sigma_g_se,
    sigma_e_se = res_dos$sigma_e_se,
    vp_cond = res_dos$vp_cond,
    h2 = res_dos$h2,
    h2_se = res_dos$h2_se,
    family_icc = fam_chk$icc,
    family_icc_note = fam_chk$note,
    stringsAsFactors = FALSE
  )
  sens_rows[[trait]] <- data.frame(
    trait = trait,
    h2_dosage = res_dos$h2,
    h2_hard = res_hard$h2,
    stringsAsFactors = FALSE
  )
}

audit_df <- do.call(rbind, audit_rows)
out_h2 <- do.call(rbind, result_rows)
out_sens <- do.call(rbind, sens_rows)

audit_path <- file.path(qtl_dir, "QTL_trait_sample_audit.csv")
h2_path <- file.path(qtl_dir, "genomic_heritability_QTL_traits.csv")
sens_path <- file.path(qtl_dir, "genomic_heritability_GRM_sensitivity.csv")

write.csv(audit_df, audit_path, row.names = FALSE)
write.csv(out_h2, h2_path, row.names = FALSE)
write.csv(out_sens, sens_path, row.names = FALSE)

message("Saved: ", audit_path)
message("Saved: ", h2_path)
message("Saved: ", sens_path)
print(audit_df)
print(out_h2)
