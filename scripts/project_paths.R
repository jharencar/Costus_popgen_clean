# Shared path resolution for publication R scripts.
# Usage at top of a script:
#   source(file.path("scripts", "project_paths.R"))
#   paths <- project_paths()
#
# Large data live on the dissert archive (not in git). Override with environment:
#   COSTUS_ARCHIVE  root archive dir (default: /Volumes/dissert/Costus_archive on macOS)
#   PIXY_FOLDER, DSUITE_DIR, LEPMAP_DIR, GENOPROB_RDS, POSTERIOR_DIR, WILD_AHMM_DIR

.first_existing_dir <- function(candidates) {
  for (c in candidates) {
    if (nzchar(c) && dir.exists(c)) return(normalizePath(c, winslash = "/", mustWork = FALSE))
  }
  NA_character_
}

.first_existing_file <- function(candidates) {
  for (c in candidates) {
    if (nzchar(c) && file.exists(c)) return(normalizePath(c, winslash = "/", mustWork = FALSE))
  }
  NA_character_
}

project_paths <- function() {
  root <- Sys.getenv("PROJECT_ROOT", "")
  if (!nzchar(root)) {
    cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    candidates <- c(
      file.path(cwd, "03_QTL_mapping", "rqtl", "Costus_QTL.qmd"),
      file.path(cwd, "README.md"),
      file.path(dirname(cwd), "README.md")
    )
    for (marker in candidates) {
      if (file.exists(marker)) {
        root <- if (basename(cwd) %in% c("rqtl", "aHMM", "pixy", "dsuite", "lepmap3")) {
          normalizePath(file.path(cwd, "..", ".."), winslash = "/")
        } else if (file.exists(file.path(cwd, "03_QTL_mapping"))) {
          cwd
        } else {
          dirname(cwd)
        }
        break
      }
    }
    if (!nzchar(root)) root <- cwd
  }
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)

  archive <- normalizePath(
    Sys.getenv("COSTUS_ARCHIVE", "/Volumes/dissert/Costus_archive"),
    winslash = "/",
    mustWork = FALSE
  )

  pixy_folder <- .first_existing_dir(c(
    Sys.getenv("PIXY_FOLDER", ""),
    file.path(archive, "06_pixy", "output_251124_LASIref"),
    file.path(archive, "06_pixy", "from_hummingbird", "pixy", "output_latest"),
    file.path(root, "pixy", "output_251124_LASIref")
  ))

  dsuite_dir <- .first_existing_dir(c(
    Sys.getenv("DSUITE_DIR", ""),
    file.path(archive, "07_dsuite", "Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120"),
    file.path(root, "Dsuite", "Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120")
  ))

  lepmap_dir <- .first_existing_dir(c(
    Sys.getenv("LEPMAP_DIR", ""),
    file.path(archive, "04_lepmap3", "lepmap3", "best_iter_ordermarkers"),
    file.path(root, "03_QTL_mapping", "lepmap3", "best_iter_ordermarkers")
  ))

  genoprob_rds <- .first_existing_file(c(
    Sys.getenv("GENOPROB_RDS", ""),
    file.path(archive, "05_rqtl",
      "aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds"),
    file.path(root, "03_QTL_mapping", "rqtl",
      "aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds")
  ))

  posterior_dir <- .first_existing_dir(c(
    Sys.getenv("POSTERIOR_DIR", ""),
    file.path(archive, "03_aHMM_posteriors", "QTL_grouped_run",
      "grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO",
      "no_thin_cated_posteriors_minDif0.7_grp_allF2s_CO0.05")
  ))

  wild_ahmm_dir <- .first_existing_dir(c(
    Sys.getenv("WILD_AHMM_DIR", ""),
    file.path(archive, "03_aHMM_posteriors", "wild_AHV_LASIref2", "AHV_LASIref2_260210"),
    file.path(root, "aHMM", "AHV_LASIref2_260210")
  ))

  rqtl_dir <- Sys.getenv("RQTL_DIR", file.path(archive, "05_rqtl"))
  scantwo_perm_dir <- Sys.getenv(
    "SCANTWO_PERM_DIR",
    .first_existing_dir(c(
      file.path(rqtl_dir, "scantwo.perms", "from_alle_vill_QTL"),
      file.path(rqtl_dir, "scantwo.perms"),
      file.path(root, "03_QTL_mapping", "rqtl", "scantwo.perms")
    ))
  )

  list(
    project_root = root,
    costus_archive = archive,
    pixy_folder = pixy_folder,
    dsuite_dir = dsuite_dir,
    lepmap_dir = lepmap_dir,
    genoprob_rds = genoprob_rds,
    posterior_dir = posterior_dir,
    wild_ahmm_dir = wild_ahmm_dir,
    vcf_dir = Sys.getenv("VCF_DIR", file.path(archive, "02_VCFs")),
    rqtl_dir = rqtl_dir,
    scantwo_perm_dir = scantwo_perm_dir,
    qtl_ci_path = Sys.getenv(
      "QTL_CI_PATH",
      file.path(root, "03_QTL_mapping", "data", "all_QTL_CIs_bp_cm.251130.csv")
    )
  )
}
