# Prepare Waters mass-alignment chemistry for QTL mapping.
# Run from project root or with working directory = QTL/ (script will switch to QTL/ if needed).
# Reads WatersMassAlingment.xlsx (sheet "in"), per-sample max-scales intensities,
# matches field_ID -> Lab_ID via Kay_lab_database.csv, merges into F3phe phenotype file.

library(dplyr)
library(tidyr)
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' required. Install with install.packages('readxl').")
}
library(readxl)

for (p in c("log_all_warnings.R", file.path("QTL", "log_all_warnings.R"))) {
  if (file.exists(p)) {
    source(p, local = FALSE)
    break
  }
}
if (!exists("start_warning_log")) {
  stop("log_all_warnings.R not found in working directory or QTL/", call. = FALSE)
}
start_warning_log(default_basename = "prepare_WatersMass_chemistry_warnings.txt")

excel_path <- "WatersMassAlingment.xlsx"
kay_path <- "Kay_lab_database.csv"
phe_file <- "F3phe_onehotcovar_newHI_20250925.csv"
meta_out <- "WatersMass_feature_metadata.csv"
scaled_out <- "WatersMass_chemistry_transposed_w_LabID.csv"
unscaled_out <- "WatersMass_chemistry_transposed_w_LabID_unscaled.csv"
phe_merged_out <- "F3phe_onehotcovar_newHI_with_WatersMass_chemistry.csv"

if (!file.exists(excel_path) && file.exists(file.path("QTL", excel_path))) {
  setwd("QTL")
  message("Set working directory to QTL/")
}
stopifnot(file.exists(excel_path), file.exists(kay_path), file.exists(phe_file))

# ---- Load aligned chemistry (features x samples) ----
raw <- read_excel(excel_path, sheet = "in")
meta_cols <- c("row ID", "row m/z", "row retention time")
sample_cols <- setdiff(names(raw), meta_cols)
stopifnot(length(meta_cols) == 3, length(sample_cols) > 0)

feature_meta <- raw %>%
  transmute(
    trait = paste0("row", as.integer(`row ID`)),
    row_ID = as.integer(`row ID`),
    mz = as.numeric(`row m/z`),
    rt = as.numeric(`row retention time`)
  )
write.csv(feature_meta, meta_out, row.names = FALSE)
message("Wrote ", meta_out, " (", nrow(feature_meta), " features)")

trait_cols <- feature_meta$trait

chem_long <- raw %>%
  mutate(trait = paste0("row", as.integer(`row ID`))) %>%
  select(trait, all_of(sample_cols)) %>%
  pivot_longer(-trait, names_to = "field_ID", values_to = "value")

chem_wide_unscaled <- chem_long %>%
  pivot_wider(names_from = trait, values_from = value)

# Per-sample max scaling (top peak = 1; skip samples where max = 0)
chem_mat <- as.matrix(chem_wide_unscaled[, trait_cols, drop = FALSE])
storage.mode(chem_mat) <- "numeric"
for (i in seq_len(nrow(chem_mat))) {
  smax <- max(chem_mat[i, ], na.rm = TRUE)
  if (is.finite(smax) && smax > 0) chem_mat[i, ] <- chem_mat[i, ] / smax
}
chem_wide_scaled <- chem_wide_unscaled
chem_wide_scaled[, trait_cols] <- as.data.frame(chem_mat)

# ---- Kay lab database: details (field_ID) -> Lab_ID ----
kay <- read.csv(kay_path, check.names = FALSE)
if (names(kay)[1] == "") names(kay)[1] <- "row"
kay <- kay %>%
  mutate(
    details = trimws(as.character(details)),
    Lab_ID = trimws(as.character(`Lab ID`))
  ) %>%
  filter(details != "", !is.na(details)) %>%
  mutate(Lab_ID = gsub("-", "_", Lab_ID, fixed = TRUE))

lookup <- kay %>%
  select(details, Lab_ID) %>%
  distinct(details, .keep_all = TRUE)

add_lab_id <- function(df) {
  out <- df %>%
    left_join(lookup, by = c("field_ID" = "details")) %>%
    relocate(Lab_ID, .after = field_ID)
  out
}

chem_unscaled <- add_lab_id(chem_wide_unscaled)
chem_scaled <- add_lab_id(chem_wide_scaled)

write.csv(chem_unscaled, unscaled_out, row.names = FALSE)
write.csv(chem_scaled, scaled_out, row.names = FALSE)
message("Wrote ", scaled_out, " and ", unscaled_out)

# ---- Match report ----
n_total <- nrow(chem_scaled)
n_kay <- sum(!is.na(chem_scaled$Lab_ID))
unmatched <- chem_scaled$field_ID[is.na(chem_scaled$Lab_ID)]

phe <- read.csv(phe_file, check.names = FALSE)
f3_ids <- unique(trimws(as.character(phe$ID)))
n_f3phe <- sum(chem_scaled$Lab_ID %in% f3_ids, na.rm = TRUE)

message(
  "Kay match: ", n_kay, "/", n_total, " chemistry samples -> Lab_ID",
  " (expected ~104/130)."
)
message(
  "F3phe overlap: ", n_f3phe, "/", n_total,
  " chemistry samples with Lab_ID in ", phe_file,
  " (expected ~83/130)."
)
if (length(unmatched) > 0) {
  message("Unmatched field_IDs (n=", length(unmatched), "): ", paste(unmatched, collapse = ", "))
}

# ---- Merge scaled chemistry onto full F3 phenotype table ----
chem_only <- chem_scaled %>% select(-field_ID)
phe_with_chem <- phe %>%
  left_join(chem_only, by = c("ID" = "Lab_ID"))
write.csv(phe_with_chem, phe_merged_out, row.names = FALSE)
message("Wrote ", phe_merged_out)
