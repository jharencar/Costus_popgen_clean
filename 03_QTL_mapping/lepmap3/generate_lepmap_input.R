#!/usr/bin/env Rscript
# aHMM posteriors -> Lep-MAP3 .call genotype file.
# Original: generate_lepmap_input_files_forJulia.R (verbatim in _originals/)
#
# Pedigree construction remains in the original script (hardcoded F2 family structure).
#
# Usage:
#   POSTERIOR_DIR=/path/to/posteriors \
#   OUT_CALL=lepmap3/genotypes.call.txt \
#   Rscript generate_lepmap_input.R

suppressPackageStartupMessages(library(dplyr))

posterior_dir <- Sys.getenv("POSTERIOR_DIR", "")
out_call <- Sys.getenv("OUT_CALL", "lepmap3_genotypes.call.txt")
prob_threshold <- as.numeric(Sys.getenv("PROB_THRESHOLD", "0.9"))

if (!nzchar(posterior_dir) || !dir.exists(posterior_dir)) {
  stop("Set POSTERIOR_DIR to directory of concatenated *.posterior files.")
}

setwd(posterior_dir)
files <- list.files(".", pattern = "*posterior")
if (!length(files)) stop("No posterior files found in ", posterior_dir)

df <- NULL
for (i in files) {
  id <- strsplit(i, ".", fixed = TRUE)[[1]][1]
  sample_id <- substr(id, 1, nchar(id) - 3)
  post <- read.table(i, sep = "\t", header = FALSE)
  names(post) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
  post <- post %>%
    mutate(gencode = case_when(
      ref_homo > prob_threshold ~ "AA",
      het > prob_threshold ~ "AB",
      alt_homo > prob_threshold ~ "BB",
      TRUE ~ "NoCall"
    ))
  if (is.null(df)) {
    df <- rbind(
      c("ID", paste0(post$chrom, "_", post$position)),
      c("", substr(post$chrom, nchar(post$chrom), nchar(post$chrom))),
      c("", post$position)
    )
  }
  df <- rbind(df, c(sample_id, post$gencode))
}

colnames(df) <- df[1, ]
df <- df[-1, ]

df.t <- t(df)
df.t[1, 1] <- "CHROM"
df.t[1, 2] <- "POS"

dir.create(dirname(out_call), showWarnings = FALSE, recursive = TRUE)
write.table(df.t, out_call, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
message("Wrote Lep-MAP3 call file: ", out_call, " (", nrow(df.t), " rows)")
