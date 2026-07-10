### Consolidating aHMM results to generate R/qtl input data
### originally by K. Uckele; Feb. 5, 2024
### Modified by J. Harencar last on Sept. 24, 2025

## clear workspace
rm(list = ls())

library(dplyr)
library(tidyr)

#### import posterior files ######
#setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/DP1.3_minDif0.2/refV2_concatenated_files_thinned_50000_minDP1.3_minDif0.2")
#setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/LASI_ref_single_samps/minGTP3_minDif0.7_dist5k_minDP0.5_pulse0.001/concatenated_F3_posteriors")
#setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/grouped_run_minDif0.7_LASIref_F1parent_F3s/concatenated_F3_posteriors")
setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO/no_thin_cated_posteriors_minDif0.7_grp_allF2s_CO0.05")

files <- list.files(".", pattern = "*posterior")
#files <- files[grep("F3", files)]

# initiate empty list
all_samples_long <- list()

# Code for if all posteriors have the same sites ####
df <- c()

# loop to import posterior data from all individuals into one df
for (i in files) {
  id <- strsplit(i, ".", fixed = TRUE)[[1]][1]
  sample_id <- substr(id, 1, nchar(id) - 3)sample_id <- substr(id, 1, nchar(id) - 3)
  #sample_id <- strsplit(i, "_", fixed = TRUE)[[1]][1]
  #pop <- substr(id, nchar(id) - 3, nchar(id))
  post <- read.table(i, sep = "\t", header = FALSE)
  # rename cols
  names(post) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
  # assign genetic codes
  post <- post %>%
    mutate(gencode = case_when(
    ref_homo > 0.9  ~ "A" ,
    het > 0.9  ~ "H" ,
    alt_homo > 0.9  ~ "B",
    ref_homo <= 0.9 & het <= 0.9 & alt_homo < 0.9 ~ "-"),
    dist_cM = position / 1e6
  )

  if(exists("head1") == "FALSE") {
    head1 <- c("ID", paste0(post$chrom, "_", post$position))
    df <- rbind(df, head1)
  }
  #if(exists("head2") == "FALSE") {
  #  head2 <- c("", rep(1, nrow(post)))
  #  df <- rbind(df, head2)
  #}
  if(exists("head2") == "FALSE") {
    head2 <- c("", substr(post$chrom, nchar(post$chrom), nchar(post$chrom)))
    df <- rbind(df, head2)
  }
  if(exists("head3") == "FALSE") {
    head3 <- c("", post$dist_cM)
    df <- rbind(df, head3)
  }
  df <- rbind(df, c(sample_id, post$gencode))
}

rownames(df) <- NULL

dim(df)
#KAU   502 16082
#JGH  394 9534 (lost one F3 somehow, but more sites bc last was trying fixed diffs: 390 393 23823)
# JGH (07/22/25): 371 42471 (368 individuals)

colnames(df) <- df[1,]
df <- df[-1,]

View(df)
write.csv(df, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/rQTL_genetic_mapping/f3gen_250722_minGTP3_minDif0.7_dist5k_minDP0.5_pulse0.0001.csv", row.names = FALSE, quote = FALSE)


#### replace third line with cM estimates from lep map ####
setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL")
order1 <- read.table("./lepmap3/order1_manMap2.mapped")
order2 <- read.table("./lepmap3/order2_manMap2.mapped")
order3 <- read.table("./lepmap3/order3_manMap2.mapped")
order4 <- read.table("./lepmap3/order4_manMap2.mapped")
order5 <- read.table("./lepmap3/order5_manMap2.mapped")
order6 <- read.table("./lepmap3/order6_manMap2.mapped")
order7 <- read.table("./lepmap3/order7_manMap2.mapped")
order8 <- read.table("./lepmap3/order8_manMap2.mapped")
order9 <- read.table("./lepmap3/order9_manMap2.mapped")

orders <- rbind(order1, order2, order3, order4, order5, order6, order7, order8, order9)

old_cM <- df[2,1:100]

for (i in 2:ncol(df)) {
  cM <- orders[orders[,1] %in% colnames(df)[i],3]
  df[2,i] <- cM
}

View(df)


#write.csv(df, "~/Documents/R:QTL/data.csv", row.names = FALSE, quote = FALSE)
write.csv(df, "F3gen_lemMap2_cm.csv", row.names = FALSE, quote = FALSE)

#write.csv(df[,1:100], "~/Documents/R:QTL/F2gen_subset.csv", row.names = FALSE, quote = FALSE)

## Test to see if there are any errors when loading into R/qtl
library(qtl)
cross <- read.cross("csv", dir="~/Documents/R:QTL", genfile = "F2gen.csv", phefile = "F2phe.csv", estimate.map=FALSE)


cross <- read.cross("csv", dir="~/Downloads", genfile = "F3gen.csv", phefile = "F3phe.csv", estimate.map=FALSE)

# # Version for singly run hybs:
# for (i in files) {
#   sample_id <- strsplit(i, "_", fixed = TRUE)[[1]][1]
#   
#   post <- read.table(i, sep = "\t", header = FALSE)
#   names(post) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
#   
#   post <- post %>%
#     mutate(
#       gencode = case_when(
#         ref_homo > 0.9 ~ "A",
#         het > 0.9 ~ "H",
#         alt_homo > 0.9 ~ "B",
#         TRUE ~ "-"
#       ),
#       site = paste0(chrom, "_", position),
#       sample_id=sample_id
#     ) %>%
#     select(sample_id, site, gencode)
#   all_samples_long[[sample_id]] <- post
# }
# 
# # Combine all samples into long format
# long_df <- bind_rows(all_samples_long)
# 
# # check dimensions
# #table(long_df$sample_id)  # should have 1 row count per site × sample; looks good
# 
# #check that all files were included:
# #length(unique(long_df$sample_id)) == length(files) # TRUE! yay
# 
# # replace all NAs from rbind with the null value preffered by rqtl ("-")
# long_df[is.na(long_df)] <- "-"
# 
# # Pivot to wide format
# wide_df <- pivot_wider(long_df, names_from = site, values_from = gencode, values_fill = "-")
# 
# 
# #TEST_wide_df <- head(wide_df[,1:5])
# # # Remove column names and sample_ID column
# # colnames(wide_df) <- wide_df[1,]
# # wide_df <- wide_df[-1,-1]
# 
# # rename ID column
# wide_df <- wide_df %>% rename(ID = sample_id)
# 
# # Separate site names into chrom and position
# site_names <- names(wide_df)[-1]  # drop the ID column
# chrom_parts <- sub("_.*", "", site_names)       # e.g., Chrom1
# chrom_num <- sub("Chrom", "", chrom_parts)      # e.g., "1"
# position_bp <- as.numeric(sub(".*_", "", site_names))  # e.g., 193991
# position_mbp <- position_bp / 1e6
# 
# # Create header rows
# header_chrom <- c("", chrom_num)
# header_pos   <- c("", sprintf("%.6f", position_mbp))  # round nicely
# 
# # Convert wide_df to character so we can rbind these rows
# wide_chr <- wide_df %>%
#   mutate(across(everything(), as.character))
# 
# # Prepend headers
# final_wide_df <- rbind(header_chrom, header_pos, wide_chr)
# #final_df[c(1,2),1] <- ""
# 
# # check that it looks right
# head(final_wide_df[,1:5])
# 
# # Write to CSV
# write.csv(final_wide_df, "f3_genotype_matrix_mindif0.7_dist5k.csv", row.names = FALSE)
