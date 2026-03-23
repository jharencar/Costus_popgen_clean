# genomic_heritability_QTL_traits.R
# Narrow-sense genomic heritability (h^2) via VanRaden GRM + REML (sommer::mmer),
# aligned with QTL trait scales and covariate sets from Costus_QTL.qmd.
# Optional: lme4::lmer with (1|family) as a coarse F2-line check (n = 5 groups).
#
# Usage (from repo root):  Rscript QTL/genomic_heritability_QTL_traits.R
# Or from QTL/:            Rscript genomic_heritability_QTL_traits.R
#
# Requires: sommer, Matrix, lme4; optional: data.table (faster I/O)

suppressPackageStartupMessages({
  if (!requireNamespace("sommer", quietly = TRUE)) {
    stop("Install sommer: install.packages(\"sommer\")")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Install Matrix")
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Install lme4: install.packages(\"lme4\")")
  }
  library(sommer)
  library(Matrix)
  library(lme4)
})

has_dt <- requireNamespace("data.table", quietly = TRUE)
if (has_dt) {
  library(data.table)
}

# ---- paths: repo root = parent of QTL/ if script lives in QTL/ ----
args <- commandArgs(trailingOnly = TRUE)
cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}
root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], winslash = "/")
} else if (basename(script_dir) == "QTL") {
  normalizePath(file.path(script_dir, ".."), winslash = "/")
} else {
  script_dir
}

phe_path <- file.path(root, "QTL", "F3phe_onehotcovar_newHI_20250925.csv")
gen_path <- file.path(root, "QTL", "f3gen_250924_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001.csv")

if (!file.exists(phe_path)) {
  phe_path <- file.path(root, "F3phe_onehotcovar_newHI_20250925.csv")
}
if (!file.exists(gen_path)) {
  gen_path <- file.path(root, "f3gen_250924_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001.csv")
}

# ---- QC for SNPs used in GRM ----
maf_min <- 0.05
max_snp_missing_rate <- 0.2

# ---- Specification (matches Costus_QTL.qmd + methods note) ----
# Phenotype scale: SD on log scale (log_days_to_radicle); other traits on original scale.
# h^2 model: y = X*beta + g + e, g ~ N(0, sigma_g^2 G), e ~ N(0, sigma_e^2 I), with the same
# fixed-effect covariates as the final scanone / fitqtl models per trait (not residualized y).
trait_cfgs <- list(
  SD = list(
    label = "SD (seed dormancy)",
    pheno_col = "days_to_radicle",
    transform = "log",
    covariates = c("s.coh4", "s.coh7", "s.coh8", "s.coh9", "s.coh14", "s.coh15")
  ),
  GR = list(
    label = "GR (growth rate)",
    pheno_col = "GR_cm.dy",
    transform = "identity",
    covariates = c(
      "hybrid_index", "F21_259", "F21_238", "F21_284", "F21_232",
      "coh7", "coh10", "coh11", "coh13", "coh18", "coh19", "coh24", "coh25",
      "coh27", "coh28", "coh29", "coh30", "coh31", "coh32", "coh34"
    )
  ),
  tough = list(
    label = "tough (ave_tough)",
    pheno_col = "ave_tough",
    transform = "identity",
    covariates = c(
      "hybrid_index", "F21_232", "F21_238", "F21_259", "F21_284",
      "coh10", "coh11", "coh16", "coh19", "coh22", "coh24", "coh27",
      "coh28", "coh29", "coh30", "coh32"
    )
  ),
  Nreabs = list(
    label = "Nreabs",
    pheno_col = "N_reabsorp",
    transform = "identity",
    covariates = c(
      "hybrid_index", "F21_238", "F21_259", "F21_232", "F21_284",
      "coh10", "coh16", "coh18", "coh21", "coh22", "coh24", "coh25",
      "coh27", "coh29", "coh30", "coh32", "coh34"
    )
  ),
  Chlreabs = list(
    label = "Chlreabs",
    pheno_col = "chloro_reabsorp",
    transform = "identity",
    covariates = c(
      "hybrid_index", "F21_238", "F21_284", "F21_232",
      "coh13", "coh16", "coh19", "coh21", "coh22", "coh24", "coh25",
      "coh27", "coh28", "coh30", "coh31", "coh32"
    )
  )
)

read_geno_matrix <- function(path) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  hdr <- readLines(con, n = 1L)
  cn <- strsplit(hdr, ",", fixed = TRUE)[[1L]]
  if (has_dt) {
    geno <- fread(path, skip = 3L, header = FALSE, data.table = FALSE)
  } else {
    geno <- utils::read.csv(path, skip = 3L, header = FALSE, check.names = FALSE)
  }
  colnames(geno) <- cn
  id <- as.character(geno[[1L]])
  M <- as.matrix(geno[, -1L, drop = FALSE])
  rownames(M) <- id
  list(id = id, M = M, marker_names = colnames(geno)[-1L])
}

#' VanRaden method 1 GRM from SNP dosages in {0,1,2}.
#' Accepts R/qtl 1/2/3 codes or allele codes A/H/B (A=0, H=1, B=2).
make_grm <- function(M, ids, maf_min = 0.05, max_miss = 0.2) {
  if (is.character(M) || is.factor(M)) {
    M <- matrix(as.character(M), nrow = nrow(M), ncol = ncol(M))
    key <- c(A = 0, H = 1, B = 2, a = 0, h = 1, b = 2)
    M <- matrix(key[M], nrow = nrow(M), ncol = ncol(M))
  }
  storage.mode(M) <- "double"
  rng <- suppressWarnings(range(M, na.rm = TRUE))
  if (all(is.finite(rng)) && rng[1L] >= 1 && rng[2L] <= 3) {
    M <- M - 1
  }
  miss <- is.na(M)
  p <- colMeans(M, na.rm = TRUE) / 2
  p[!is.finite(p)] <- 0.5
  snp_miss <- colMeans(miss)
  keep <- which(p >= maf_min & p <= (1 - maf_min) & snp_miss <= max_miss)
  if (length(keep) < 10L) {
    stop("Too few SNPs after MAF/missing filters")
  }
  M <- M[, keep, drop = FALSE]
  miss <- is.na(M)
  for (j in seq_len(ncol(M))) {
    v <- M[, j]
    mu <- mean(v, na.rm = TRUE)
    if (!is.finite(mu)) {
      mu <- 2 * p[j]
    }
    v[miss[, j]] <- mu
    M[, j] <- v
  }
  p <- colMeans(M) / 2
  Z <- sweep(M, 2L, 2 * p, "-")
  denom <- 2 * sum(p * (1 - p))
  if (denom <= 0) {
    stop("Invalid denominator for VanRaden GRM")
  }
  G <- tcrossprod(Z) / denom
  dimnames(G) <- list(ids, ids)
  list(G = G, n_snps = ncol(M))
}

prepare_trait_pheno <- function(phe, cfg) {
  y_raw <- phe[[cfg$pheno_col]]
  if (is.null(y_raw)) {
    stop("Missing phenotype column: ", cfg$pheno_col)
  }
  y <- if (identical(cfg$transform, "log")) {
    log(as.numeric(y_raw))
  } else {
    as.numeric(y_raw)
  }
  cov_cols <- cfg$covariates
  miss_cov <- cov_cols[!cov_cols %in% names(phe)]
  if (length(miss_cov)) {
    stop("Missing covariate columns: ", paste(miss_cov, collapse = ", "))
  }
  X <- as.data.frame(phe[, cov_cols, drop = FALSE])
  for (nm in names(X)) {
    X[[nm]] <- as.numeric(X[[nm]])
  }
  data.frame(
    ID = as.character(phe$ID),
    family = as.character(phe$family),
    y = y,
    X,
    stringsAsFactors = FALSE
  )
}

quote_term <- function(nm) {
  if (make.names(nm) == nm) {
    nm
  } else {
    paste0("`", gsub("`", "\\\\`", nm, fixed = TRUE), "`")
  }
}

fit_greml <- function(dat, G) {
  ids <- dat$ID
  if (!all(ids %in% rownames(G))) {
    stop("ID/G mismatch")
  }
  Gsub <- G[ids, ids, drop = FALSE]
  ev <- eigen(Gsub, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 1e-8) {
    Gsub <- as.matrix(Matrix::nearPD(Gsub)$mat)
    dimnames(Gsub) <- list(ids, ids)
  }
  dat$IDf <- factor(ids, levels = ids)
  covnms <- setdiff(names(dat), c("ID", "family", "y", "IDf"))
  rhs <- paste(vapply(covnms, quote_term, character(1L)), collapse = " + ")
  f <- if (nzchar(rhs)) {
    as.formula(paste0("y ~ ", rhs))
  } else {
    y ~ 1
  }
  mmer(
    f,
    random = ~ vsr(IDf, Gu = Gsub),
    rcov = ~ units,
    data = dat,
    verbose = FALSE
  )
}

h2_from_mmer <- function(fit) {
  s <- summary(fit)$varcomp
  if (is.null(s) || !nrow(s)) {
    return(list(h2 = NA_real_, se = NA_real_, sigma_g2 = NA_real_, sigma_e2 = NA_real_))
  }
  rn <- rownames(s)
  ig <- grep("^u:IDf", rn)
  ie <- grep("^units", rn)
  if (!length(ig) || !length(ie)) {
    ig <- grep("IDf", rn)
    ie <- grep("^units", rn)
  }
  if (!length(ig) || !length(ie)) {
    return(list(h2 = NA_real_, se = NA_real_, sigma_g2 = NA_real_, sigma_e2 = NA_real_))
  }
  sg2 <- s[ig[1L], "VarComp"]
  se2 <- s[ie[1L], "VarComp"]
  h2 <- as.numeric(sg2 / (sg2 + se2))
  se_h2 <- NA_real_
  s_se <- s[ig[1L], "VarCompSE"]
  e_se <- s[ie[1L], "VarCompSE"]
  if (is.finite(s_se) && is.finite(e_se)) {
    dhds <- se2 / (sg2 + se2)^2
    dhde <- -sg2 / (sg2 + se2)^2
    se_h2 <- sqrt(dhds^2 * s_se^2 + dhde^2 * e_se^2)
  }
  list(h2 = h2, se = se_h2, sigma_g2 = sg2, sigma_e2 = se2)
}

# Coarse F2-line check: Var(family)/(Var(family)+Var(residual)) with the same fixed
# covariates as GBLUP. Five F2 families only—wide uncertainty; not equivalent to h^2.
fit_family_lmer <- function(dat) {
  covnms <- setdiff(names(dat), c("ID", "family", "y", "IDf"))
  rhs <- paste(vapply(covnms, quote_term, character(1L)), collapse = " + ")
  if (nzchar(rhs)) {
    f <- as.formula(paste0("y ~ ", rhs, " + (1|family)"))
  } else {
    f <- y ~ (1 | family)
  }
  fit <- tryCatch(
    lmer(f, data = dat, REML = TRUE, control = lmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(icc_family = NA_real_, var_family = NA_real_, var_resid = NA_real_))
  }
  vc <- as.data.frame(VarCorr(fit))
  vf <- vc$vcov[vc$grp == "family"]
  vr <- vc$vcov[vc$grp == "Residual"]
  if (!length(vf) || !length(vr)) {
    return(list(icc_family = NA_real_, var_family = NA_real_, var_resid = NA_real_))
  }
  icc <- vf / (vf + vr)
  list(icc_family = as.numeric(icc), var_family = as.numeric(vf), var_resid = as.numeric(vr))
}

# ---- load phenotypes ----
phe <- if (has_dt) {
  fread(phe_path, data.table = FALSE)
} else {
  read.csv(phe_path, stringsAsFactors = FALSE)
}

message("Reading genotypes from ", gen_path, " ...")
geno <- read_geno_matrix(gen_path)
message("Building GRM (MAF >= ", maf_min, ") ...")
grm <- make_grm(geno$M, ids = geno$id, maf_min = maf_min, max_miss = max_snp_missing_rate)
G <- grm$G
message("GRM uses ", grm$n_snps, " SNPs.")

results <- list()

for (tn in names(trait_cfgs)) {
  cfg <- trait_cfgs[[tn]]
  df <- prepare_trait_pheno(phe, cfg)
  ok <- stats::complete.cases(df[c("y", names(df)[names(df) %in% cfg$covariates])])
  df <- df[ok, , drop = FALSE]
  # intersect with genotyped IDs
  df <- df[df$ID %in% rownames(G), , drop = FALSE]
  message(
    tn, ": n = ", nrow(df), " (complete cases + genotyped)"
  )
  if (nrow(df) < 30L) {
    warning(tn, ": too few individuals; skipping")
    next
  }
  df$family <- factor(df$family)
  fit <- tryCatch(
    fit_greml(df, G),
    error = function(e) {
      warning(tn, ": mmer failed: ", conditionMessage(e))
      NULL
    }
  )
  h2s <- if (!is.null(fit)) h2_from_mmer(fit) else list(h2 = NA_real_, se = NA_real_, sigma_g2 = NA_real_, sigma_e2 = NA_real_)
  fam <- fit_family_lmer(df)
  results[[length(results) + 1L]] <- data.frame(
    trait = tn,
    label = cfg$label,
    n = nrow(df),
    n_snps_grm = grm$n_snps,
    transform = cfg$transform,
    h2_g = h2s$h2,
    se_h2_g = h2s$se,
    sigma_g2 = h2s$sigma_g2,
    sigma_e2 = h2s$sigma_e2,
    icc_family_lmer = fam$icc_family,
    var_family = fam$var_family,
    var_resid_lmer = fam$var_resid,
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, results)
row.names(out) <- NULL

out_path <- file.path(root, "QTL", "genomic_heritability_QTL_traits.csv")
if (!dir.exists(dirname(out_path))) {
  out_path <- file.path(root, "genomic_heritability_QTL_traits.csv")
}
write.csv(out, out_path, row.names = FALSE)
message("Wrote ", out_path)

print(out)
