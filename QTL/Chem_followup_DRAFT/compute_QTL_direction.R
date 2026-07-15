# Peak-marker allele effects for chemistry follow-up traits.
# Allele coding (confirmed vs hybrid_index = villosissimus ancestry):
#   AA = villosissimus, BB = allenii

.libPaths(c("/Users/admin-jgharenc/Documents/Github/Costus_popgen_clean/R_libs", .libPaths()))
library(qtl)

root <- "/Users/admin-jgharenc/Documents/Github/Costus_popgen_clean"
setwd(file.path(root, "QTL"))

genoprob_path <- "/Users/admin-jgharenc/Documents/Github/alle_vill_QTL/aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds"
Costus_prob <- readRDS(genoprob_path)
phe <- read.csv("F3phe_onehotcovar_newHI_with_WatersMass_chemistry.csv", check.names = FALSE)
follow <- read.csv("Chem_followup_DRAFT/Chem_QTL_followup_traits_DRAFT_for_review.csv", check.names = FALSE)
summ <- read.csv("Chem_260703/WatersMass_chemistry_scanone_summary.csv", check.names = FALSE)

resolve_cross_ids <- function(cross, phe_ids) {
  pheno <- cross$pheno
  phe_ids <- as.character(phe_ids)
  candidates <- list()
  if (!is.null(rownames(pheno))) candidates$rownames <- as.character(rownames(pheno))
  if ("ID" %in% colnames(pheno)) candidates$ID <- as.character(pheno$ID)
  rates <- vapply(candidates, function(ids) mean(ids %in% phe_ids), numeric(1))
  candidates[[names(which.max(rates))]]
}

cross_ids <- resolve_cross_ids(Costus_prob, phe$ID)
phe_idx <- match(cross_ids, as.character(phe$ID))
traits <- unique(as.character(follow$trait))
for (tr in traits) {
  vals <- rep(NA_real_, length(cross_ids))
  vals[!is.na(phe_idx)] <- as.numeric(phe[[tr]][phe_idx[!is.na(phe_idx)]])
  Costus_prob$pheno[, tr] <- vals
}

hi <- as.numeric(pull.pheno(Costus_prob, "hybrid_index"))
doseB <- rowMeans(pull.geno(Costus_prob) - 1, na.rm = TRUE)
message(sprintf(
  "Allele coding check: cor(hybrid_index=villos ancestry, B dosage)=%.3f => AA=villos, BB=allenii",
  cor(hi, doseB, use = "pairwise")
))

min_covar_n <- 10
pn <- names(Costus_prob$pheno)
fam_onehot <- c("F21_238", "F21_259", "F21_241", "F21_284", "F21_232")
coh_cols <- pn[grepl("^coh[0-9]+$", pn)]
coh_n <- colSums(Costus_prob$pheno[, coh_cols, drop = FALSE], na.rm = TRUE)
coh_ok <- coh_cols[coh_n >= min_covar_n]
candidate_covars <- c("hybrid_index", fam_onehot[fam_onehot %in% pn], coh_ok)
candidate_covars <- candidate_covars[candidate_covars %in% pn]
fam_onehot_all <- fam_onehot

ensure_full_rank_covariate_names <- function(nms) {
  out <- nms
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
  if (sum(ok) < min_complete) return(character(0))
  dat <- cbind(y = y[ok], X[ok, , drop = FALSE])
  fit <- lm(y ~ ., data = dat)
  fit_step <- step(fit, trace = 0)
  sel <- setdiff(names(coef(fit_step)), "(Intercept)")
  ensure_full_rank_covariate_names(sel)
}

rows <- vector("list", nrow(follow))
for (i in seq_len(nrow(follow))) {
  tr <- as.character(follow$trait[i])
  sr <- summ[summ$trait == tr, , drop = FALSE][1, ]
  chr <- as.character(sr$chr)
  pos <- as.numeric(sr$pos)
  mk <- find.marker(Costus_prob, chr = chr, pos = pos)
  geno_mat <- pull.geno(Costus_prob, chr = chr)
  g <- geno_mat[, mk]
  y <- as.numeric(pull.pheno(Costus_prob, tr))
  cov_names <- select_covariates_for_trait(Costus_prob, tr, candidate_covars)

  raw_means <- tapply(y, g, mean, na.rm = TRUE)
  raw_n <- tapply(y, g, function(z) sum(!is.na(z)))
  mean_AA_raw <- unname(raw_means["1"])
  mean_H_raw <- unname(raw_means["2"])
  mean_BB_raw <- unname(raw_means["3"])

  g_f <- factor(g, levels = c(1, 2, 3), labels = c("AA", "H", "BB"))
  dat <- data.frame(y = y, g = g_f, stringsAsFactors = FALSE)
  if (length(cov_names)) {
    for (cn in cov_names) dat[[cn]] <- as.numeric(Costus_prob$pheno[[cn]])
  }
  ok <- complete.cases(dat)
  dat <- dat[ok, , drop = FALSE]

  if (nrow(dat) < 20 || nlevels(droplevels(dat$g)) < 2) {
    high_parent <- NA_character_
    delta_BB_AA <- NA_real_
    mean_AA_adj <- mean_AA_raw
    mean_H_adj <- mean_H_raw
    mean_BB_adj <- mean_BB_raw
    cor_add <- NA_real_
  } else {
    form <- if (length(cov_names)) {
      as.formula(paste("y ~ g +", paste(cov_names, collapse = " + ")))
    } else {
      y ~ g
    }
    fit <- lm(form, data = dat)
    newdata <- data.frame(g = factor(c("AA", "H", "BB"), levels = c("AA", "H", "BB")))
    if (length(cov_names)) {
      for (cn in cov_names) newdata[[cn]] <- mean(dat[[cn]], na.rm = TRUE)
    }
    preds <- as.numeric(predict(fit, newdata = newdata))
    mean_AA_adj <- preds[1]
    mean_H_adj <- preds[2]
    mean_BB_adj <- preds[3]
    delta_BB_AA <- mean_BB_adj - mean_AA_adj
    if (is.na(delta_BB_AA) || abs(delta_BB_AA) < 1e-12) {
      high_parent <- "unclear"
    } else if (delta_BB_AA > 0) {
      high_parent <- "allenii"
    } else {
      high_parent <- "villosissimus"
    }

    dat2 <- data.frame(y = y, doseB = as.numeric(g) - 1)
    if (length(cov_names)) {
      for (cn in cov_names) dat2[[cn]] <- as.numeric(Costus_prob$pheno[[cn]])
    }
    ok2 <- complete.cases(dat2)
    dat2 <- dat2[ok2, , drop = FALSE]
    if (length(cov_names)) {
      y_res <- resid(lm(as.formula(paste("y ~", paste(cov_names, collapse = "+"))), data = dat2))
    } else {
      y_res <- dat2$y - mean(dat2$y)
    }
    cor_add <- cor(y_res, dat2$doseB, use = "pairwise")
  }

  rows[[i]] <- data.frame(
    group = follow$group[i],
    trait = tr,
    row_ID = follow$row_ID[i],
    chr = follow$chr[i],
    pos = follow$pos[i],
    lod = follow$lod[i],
    marker = mk,
    chem_class = follow$chem_class[i],
    n_AA = unname(raw_n["1"]),
    n_H = unname(raw_n["2"]),
    n_BB = unname(raw_n["3"]),
    mean_AA_raw = mean_AA_raw,
    mean_H_raw = mean_H_raw,
    mean_BB_raw = mean_BB_raw,
    mean_AA_adj = mean_AA_adj,
    mean_H_adj = mean_H_adj,
    mean_BB_adj = mean_BB_adj,
    delta_BB_minus_AA_adj = delta_BB_AA,
    cor_Bdosage_resid = cor_add,
    high_expression_parent = high_parent,
    covars = paste(cov_names, collapse = ";"),
    scantwo_region = follow$scantwo_region[i],
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, rows)
out_path <- "Chem_followup_DRAFT/Chem_QTL_followup_QTL_direction.csv"
write.csv(out, out_path, row.names = FALSE)
message("Wrote ", out_path)
print(table(out$group, out$high_expression_parent))

show <- out[, c(
  "group", "trait", "chr", "pos", "lod", "chem_class",
  "mean_AA_adj", "mean_H_adj", "mean_BB_adj",
  "delta_BB_minus_AA_adj", "high_expression_parent", "scantwo_region"
)]
show$pos <- round(show$pos, 1)
show$lod <- round(show$lod, 2)
show$mean_AA_adj <- round(show$mean_AA_adj, 3)
show$mean_H_adj <- round(show$mean_H_adj, 3)
show$mean_BB_adj <- round(show$mean_BB_adj, 3)
show$delta_BB_minus_AA_adj <- round(show$delta_BB_minus_AA_adj, 3)
print(show, row.names = FALSE)
