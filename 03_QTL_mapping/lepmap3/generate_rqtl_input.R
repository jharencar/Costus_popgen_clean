#!/usr/bin/env Rscript
# aHMM posteriors -> rQTL genotype CSV. Original: generate_rQTL_input_files.R
suppressPackageStartupMessages(library(dplyr))
posterior_dir <- Sys.getenv("POSTERIOR_DIR", "")
lepmap_dir <- Sys.getenv("LEPMAP_DIR", "")
out_csv <- Sys.getenv("OUT_CSV", "f3gen_rqtl.csv")
if (!nzchar(posterior_dir) || !dir.exists(posterior_dir)) stop("Set POSTERIOR_DIR")
setwd(posterior_dir)
files <- list.files(".", pattern = "*posterior")
df <- NULL
for (i in files) {
  id <- strsplit(i, ".", fixed = TRUE)[[1]][1]
  sample_id <- substr(id, 1, nchar(id) - 3)
  post <- read.table(i, sep = "\t", header = FALSE)
  names(post) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
  post <- post %>% mutate(gencode = case_when(
    ref_homo > 0.9 ~ "A", het > 0.9 ~ "H", alt_homo > 0.9 ~ "B", TRUE ~ "-"),
    dist_cM = position / 1e6)
  if (is.null(df)) {
    df <- rbind(c("ID", paste0(post$chrom, "_", post$position)),
                c("", substr(post$chrom, nchar(post$chrom), nchar(post$chrom))),
                c("", post$dist_cM))
  }
  df <- rbind(df, c(sample_id, post$gencode))
}
colnames(df) <- df[1, ]; df <- df[-1, ]
if (nzchar(lepmap_dir) && dir.exists(lepmap_dir)) {
  orders <- do.call(rbind, lapply(Sys.glob(file.path(lepmap_dir, "order*_manMap2.mapped")), read.table))
  for (j in 2:ncol(df)) {
    cM <- orders[orders[,1] %in% colnames(df)[j], 3]
    if (length(cM) == 1) df[2, j] <- cM
  }
}
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(df, out_csv, row.names = FALSE, quote = FALSE)
message("Wrote ", out_csv)
