# Match prelim_UV_chemistry sample IDs (field IDs) to genetics/phenotype IDs (Lab ID)
# via Kay_lab_database. Output: transposed chemistry with field_ID and Lab_ID for
# use with f3gen and F3phe_onehotcovar (which use Lab ID format).

library(tidyverse)

# Paths (relative to project root or set as needed)
chem_path   <- "QTL/prelim_UV_chemistry.csv"
kay_path    <- "QTL/Kay_lab_database.csv"
output_path <- "QTL/prelim_UV_chemistry_transposed_w_LabID.csv"

# ---- Load chemistry ----
chem <- read.csv(chem_path, check.names = FALSE)

# First 7 columns: (empty/peak), rt, start, end, sd, FWHM, r.squared
# Column 8 onward: sample IDs (field_ID) as column names
n_meta <- 7
if (names(chem)[1] == "") names(chem)[1] <- "peak"
peak_col <- names(chem)[1]

# Transpose: one row per sample (field_ID), columns = peak IDs (V4, V6, ...)
chem_long <- chem %>%
  select(1, (n_meta + 1):ncol(chem)) %>%
  pivot_longer(-1, names_to = "field_ID", values_to = "value")

chem_wide <- chem_long %>%
  pivot_wider(names_from = all_of(peak_col), values_from = value)

# field_ID is first column; columns 2+ are peak IDs
stopifnot(names(chem_wide)[1] == "field_ID")

# ---- Load Kay lab database: details -> Lab ID ----
kay <- read.csv(kay_path, check.names = FALSE)
if (names(kay)[1] == "") names(kay)[1] <- "row"
kay <- kay %>%
  mutate(
    details = trimws(as.character(details)),
    Lab_ID  = trimws(as.character(`Lab ID`))
  ) %>%
  filter(details != "", !is.na(details))

# Normalize Lab ID to underscore (e.g. 22-195 -> 22_195) to match F3phe/f3gen
kay <- kay %>%
  mutate(Lab_ID = gsub("-", "_", Lab_ID, fixed = TRUE))

# One row per details (field ID); keep first Lab ID if duplicate details
lookup <- kay %>%
  select(details, Lab_ID) %>%
  distinct(details, .keep_all = TRUE)

# ---- Add Lab_ID to chemistry ----
chem_out <- chem_wide %>%
  left_join(lookup, by = c("field_ID" = "details"))

# Put Lab_ID as second column (field_ID, Lab_ID, then peaks)
lab_id_col <- chem_out %>% select(Lab_ID)
chem_out <- chem_out %>% select(-Lab_ID)
chem_out <- chem_out %>%
  bind_cols(lab_id_col) %>%
  select(field_ID, Lab_ID, everything())

# Report matches
n_matched <- sum(!is.na(chem_out$Lab_ID))
n_total   <- nrow(chem_out)
message(
  "Matched ", n_matched, " of ", n_total, " chemistry samples to Lab ID. ",
  "Unmatched field_IDs: ",
  paste(chem_out$field_ID[is.na(chem_out$Lab_ID)], collapse = ", ")
)

# Save
write.csv(chem_out, output_path, row.names = FALSE)
message("Wrote ", output_path)
