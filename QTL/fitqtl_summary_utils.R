# Shared parsers for fitqtl summary text files (*_QTL_model_output_251203).
# Source from genomic_heritability_QTL_traits.R and genetic_variance_partition.R.

#' Parse a single fitqtl summary text file.
#' @return list(n_obs, formula_rhs, full_model_pctvar, ss_total, ss_model, ss_error,
#'   drop_one data.frame: term, pctvar, ss_type3)
parse_fitqtl_summary_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  n_obs <- NA_integer_
  formula_rhs <- NA_character_
  full_pctvar <- NA_real_
  ss_total <- ss_model <- ss_error <- NA_real_

  for (i in seq_along(lines)) {
    if (grepl("^Number of observations", lines[i])) {
      n_obs <- as.integer(sub(".*:\\s*", "", lines[i]))
    }
    if (grepl("^Model formula:", lines[i])) {
      formula_rhs <- trimws(sub("^Model formula:\\s*y\\s*~\\s*", "", lines[i]))
      j <- i + 1L
      while (j <= length(lines) && !grepl("^\\s+df\\s+", lines[j]) && !grepl("^Model\\s+", lines[j])) {
        if (grepl("~", lines[j]) || grepl("^\\s+[A-Za-z0-9_.]+", lines[j])) {
          extra <- trimws(lines[j])
          if (nzchar(extra) && !grepl("^df\\s", extra)) {
            formula_rhs <- paste(formula_rhs, extra)
          }
        }
        j <- j + 1L
      }
      formula_rhs <- gsub("\\s+", " ", formula_rhs)
    }
  }

  in_full <- FALSE
  for (i in seq_along(lines)) {
    if (grepl("Full model result", lines[i], fixed = TRUE)) in_full <- TRUE
    if (!in_full) next
    if (grepl("^Total\\s+", lines[i])) {
      parts <- strsplit(trimws(lines[i]), "\\s+")[[1]]
      if (length(parts) >= 3) ss_total <- as.numeric(parts[3])
    }
    if (grepl("^Error\\s+", lines[i])) {
      parts <- strsplit(trimws(lines[i]), "\\s+")[[1]]
      if (length(parts) >= 3) ss_error <- as.numeric(parts[3])
    }
    if (grepl("^Model\\s+", lines[i])) {
      parts <- strsplit(trimws(lines[i]), "\\s+")[[1]]
      if (length(parts) >= 6 && !is.na(suppressWarnings(as.numeric(parts[2])))) {
        ss_model <- as.numeric(parts[3])
        full_pctvar <- as.numeric(parts[6])
      }
    }
    if (grepl("^Drop one QTL", lines[i], fixed = TRUE)) break
  }

  drop_one <- data.frame(
    term = character(),
    pctvar = numeric(),
    ss_type3 = numeric(),
    stringsAsFactors = FALSE
  )
  in_drop <- FALSE
  for (i in seq_along(lines)) {
    if (grepl("Drop one QTL at a time", lines[i], fixed = TRUE)) {
      in_drop <- TRUE
      next
    }
    if (in_drop && trimws(lines[i]) == "---") break
    if (in_drop && nzchar(trimws(lines[i]))) {
      parts <- strsplit(trimws(lines[i]), "\\s+")[[1]]
      if (identical(parts[1], "df") || identical(parts[1], "Model")) next
      if (length(parts) >= 5 && !is.na(suppressWarnings(as.numeric(parts[2])))) {
        drop_one <- rbind(
          drop_one,
          data.frame(
            term = parts[1],
            pctvar = as.numeric(parts[5]),
            ss_type3 = as.numeric(parts[3]),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  list(
    n_obs = n_obs,
    formula_rhs = formula_rhs,
    full_model_pctvar = full_pctvar,
    ss_total = ss_total,
    ss_model = ss_model,
    ss_error = ss_error,
    drop_one = drop_one
  )
}

#' Classify drop-one term into component.
classify_fitqtl_term <- function(term) {
  term <- trimws(term)
  if (grepl("^[0-9]+@[0-9.]+$", term)) return("QTL")
  if (term == "hybrid_index") return("other_genotypic")
  if (grepl("^F21_", term)) return("other_genotypic")
  if (grepl("^coh[0-9]+$", term) || grepl("^s\\.coh", term)) return("environment")
  NA_character_
}

#' Covariate names for REML (exclude QTL terms chr@pos).
covariates_from_formula_rhs <- function(formula_rhs) {
  terms <- strsplit(formula_rhs, "\\s*\\+\\s*")[[1]]
  terms <- trimws(terms)
  terms[!grepl("^[0-9]+@[0-9.]+$", terms) & !grepl("^Q[0-9]+$", terms)]
}

#' Trait phenotype column and optional transform.
trait_pheno_config <- function(trait) {
  switch(trait,
    SD = list(pheno = "log_days_to_radicle", label = "seed dormancy (log)"),
    GR = list(pheno = "GR_cm.dy", label = "growth rate"),
    tough = list(pheno = "ave_tough", label = "leaf toughness"),
    Nreabs = list(pheno = "N_reabsorp", label = "N reabsorption"),
    Chlreabs = list(pheno = "chloro_reabsorp", label = "chlorophyll reabsorption"),
    stop("Unknown trait: ", trait)
  )
}
