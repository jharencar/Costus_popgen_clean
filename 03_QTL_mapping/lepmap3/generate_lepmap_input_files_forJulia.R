### Generating lep-map-3 input files
### K. Uckele
### March 12, 2024 ####
#### modified by J. Harencar 2025

## clear workspace
rm(list = ls())

#library(readr)
library(dplyr)

# setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/shared_homo_sites/minGTP2_minDif0.1_minDP0.2/F2_F3s/fixed_diffs_minGTP2_minDif0.1_minDP0.2")
setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO/no_thin_cated_posteriors_minDif0.7_grp_allF2s_CO0.05")

files <- list.files(".", pattern = "*posterior")
#files <- c(files[grep("HYBR", files)], "19_117_LASI.posterior", "19_211_BRAC.posterior")

df <- c()

for (i in files) {
  id <- strsplit(i, ".", fixed = TRUE)[[1]][1]
  sample_id <- substr(id, 1, nchar(id) - 3)
  #pop <- substr(id, nchar(id) - 3, nchar(id))
  post <- read.table(i, sep = "\t", header = FALSE)
  # rename cols 
  names(post) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
  # assign genetic codes
  post <- post %>% 
    mutate(gencode = case_when(
      ref_homo > 0.9  ~ "AA" ,
      het > 0.9  ~ "AB" ,
      alt_homo > 0.9  ~ "BB",
      ref_homo <= 0.9 & het <= 0.9 & alt_homo < 0.9 ~ "NoCall"),
      #dist_cM = position / 1e6
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
    head3 <- c("", post$position) # change back to post$dist_cM and un comment the calc above to get rough cM if needed
    df <- rbind(df, head3)
  }
  df <- rbind(df, c(sample_id, post$gencode))
}

rownames(df) <- NULL

dim(df)
#KAU   502 16082
#JGH   415  5304
#JGH after removing 24 F3s with less than 0.2x ave coverage and only including 5 parental F2s, not all 21. 
# 377 10718 # thinned to 50k distance 
# might want to consider more intentional thinning... like deep the marker with the best data in a certain window??

#JGH w/o thinning and now also including non-parental F2s (16 additional F2s)
# there should be 43200 markers based on raw aHMM output - seems right
# I think there should only be 377+16 F2s individuals - 393 indivs - great
# 393 43201

#w/o thinning and with 23 0.05-0.2 coverage individuals added back in
# 414 40917

 # PICKUP! for some reason, not getting the marker names as col names like before?> 
colnames(df) <- df[1,]
df <- df[-1,]

View(df)

write.csv(df, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/f3gen_250731_allF2s_minGTP3_minDif0.7_dist5k_minDP0.5_pulse0.0001.csv", row.names = FALSE, quote = FALSE)

#write.csv(df[,1:100], "~/Documents/R:QTL/F2gen_subset.csv", row.names = FALSE, quote = FALSE)

## Test to see if there are any errors when loading into R/qtl
#library(qtl)
#cross <- read.cross("csv", dir="/Users/juliaharencar/Documents/Github/alle_vill_QTL", genfile = "f3gen_fxDif.csv", phefile = "F2phe.csv", estimate.map=FALSE)



################################################################################
## create input file for lep map 3 
################################################################################
# 
#df <- df[-1,]
df.t <- t(df)
#df.t[,1] <- paste0("Chrom",df.t[,1])

# name chrom and pos cols as first row val
df.t[1,1] <- "CHROM"
df.t[1,2] <- "POS"


#View(df.t)

#rownames(df.t)[1] <- "SNP"


write.table(df.t, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/affx_minDif0.7_genotypes_250819_no_thin_allF2s_0.05covCO.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# write.table(df.t, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/affx_minDif0.7_genotypes_250724.txt", quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
# write.table(df.t, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/affx_fxDif_genotypes_241028.txt", quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
df.t <- read.table("/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/affx_fxDif_genotypes_thinned_50k_241025.txt")
#################################################################################
## create pedigree file for lep map 3
#################################################################################
## dimensions should be: 6 rows, 499+2 columns

#CHR POS F   F   F   F   F   F
#CHR POS female  male    progeny_1   progeny_2   progeny_3   progeny_4
#CHR POS 0   0   male    male    male    male
#CHR POS 0   0   female  female  female  female
#CHR POS 2   1   0   0   0   0
#CHR POS 0   0   0   0   0   0

#First line is the family name
# second line is individual name
# third and fourth are the father and mother
# Line 5 contains the sex of each individual (1 male, 2 female, 0 unknown)
# sixth line is the phenotype (can be 0 for all individuals, this is not currently used). 

firstline <- c("CHR", "POS", rep("F", 420))
secondline <- c("CHR", "POS", df.t[1,-1])
thirdline <- c("CHR", "POS","0","0", "19_467","19_467", rep("19_615M", 26), "21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_232M","21_238M","21_238M","21_232M","21_238M","21_238M","21_238M","21_232M","21_238M","21_238M","21_238M","21_238M","21_232M","21_238M","21_238M","21_232M","21_238M","21_232M","21_232M","21_238M","21_232M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_238M","21_238M","21_259M","21_259M","21_284M","21_284M","21_259M","21_284M","21_259M","21_259M","21_259M","21_259M","21_259M","21_284M","21_284M","21_284M","21_259M","21_259M","21_284M","21_259M","21_259M","21_259M","21_259M","21_284M","21_259M","21_284M","21_284M","21_259M","21_284M","21_284M","21_284M","21_284M","21_259M","21_259M","21_284M","21_284M","21_259M","21_284M","21_259M","21_259M","21_259M","21_259M","21_284M","21_259M","21_259M","21_259M","21_284M","21_284M","21_259M","21_284M","21_284M","21_284M","21_284M","21_284M","21_284M","21_232M","21_232M","21_232M","21_238M","21_232M","21_238M","21_238M","21_238M","21_238M","21_232M","21_238M","21_238M","21_232M","21_238M","21_232M","21_238M","21_238M","21_238M","21_259M","21_238M","21_259M","21_238M","21_238M","21_238M","21_238M","21_238M","21_284M","21_232M","21_232M","21_232M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_284M","21_284M","21_259M","21_259M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_284M","21_284M","21_232M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_232M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_284M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_232M","21_238M","21_238M","21_238M","21_238M","21_232M","21_241M","21_241M","21_284M","21_241M","21_232M","21_241M","21_241M","21_241M","21_238M","21_238M","21_241M","21_238M","21_232M","21_238M","21_238M","21_241M","21_241M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_232M","21_232M","21_232M","21_238M","21_232M","21_238M","21_238M","21_238M","21_232M","21_238M","21_238M","21_238M","21_238M","21_238M","21_238M","21_232M","21_238M","21_241M","21_259M","21_238M","21_232M","21_232M","21_238M","21_238M","21_259M","21_232M","21_232M","21_232M","21_232M","21_238M","21_232M","21_284M") # not actually sure if this is what I should do here
fourthline <- c("CHR", "POS", "0", "0", "19_500", "19_500", rep("19_615F", 26),"21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_232F","21_238F","21_238F","21_232F","21_238F","21_238F","21_238F","21_232F","21_238F","21_238F","21_238F","21_238F","21_232F","21_238F","21_238F","21_232F","21_238F","21_232F","21_232F","21_238F","21_232F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_238F","21_238F","21_259F","21_259F","21_284F","21_284F","21_259F","21_284F","21_259F","21_259F","21_259F","21_259F","21_259F","21_284F","21_284F","21_284F","21_259F","21_259F","21_284F","21_259F","21_259F","21_259F","21_259F","21_284F","21_259F","21_284F","21_284F","21_259F","21_284F","21_284F","21_284F","21_284F","21_259F","21_259F","21_284F","21_284F","21_259F","21_284F","21_259F","21_259F","21_259F","21_259F","21_284F","21_259F","21_259F","21_259F","21_284F","21_284F","21_259F","21_284F","21_284F","21_284F","21_284F","21_284F","21_284F","21_232F","21_232F","21_232F","21_238F","21_232F","21_238F","21_238F","21_238F","21_238F","21_232F","21_238F","21_238F","21_232F","21_238F","21_232F","21_238F","21_238F","21_238F","21_259F","21_238F","21_259F","21_238F","21_238F","21_238F","21_238F","21_238F","21_284F","21_232F","21_232F","21_232F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_284F","21_284F","21_259F","21_259F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_284F","21_284F","21_232F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_232F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_284F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_232F","21_238F","21_238F","21_238F","21_238F","21_232F","21_241F","21_241F","21_284F","21_241F","21_232F","21_241F","21_241F","21_241F","21_238F","21_238F","21_241F","21_238F","21_232F","21_238F","21_238F","21_241F","21_241F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_232F","21_232F","21_232F","21_238F","21_232F","21_238F","21_238F","21_238F","21_232F","21_238F","21_238F","21_238F","21_238F","21_238F","21_238F","21_232F","21_238F","21_241F","21_259F","21_238F","21_232F","21_232F","21_238F","21_238F","21_259F","21_232F","21_232F","21_232F","21_232F","21_238F","21_232F","21_284F") 
fifthline <- c("CHR", "POS", "1", "2", "1", "2", rep("0", 3), "1", "2", rep("0", 3), "1", "2", rep("0", 2), "1", "2", "0", "1","2", rep("0", 7), "1", "2", rep("0", 390)) # Is this correct? 1 = vill; 2 = alle (best guess is that alle was mom to wild, selfed F1)
sixthline <- c("CHR", "POS", rep("0", 420))

ped <- rbind(firstline, secondline, thirdline, fourthline, fifthline, sixthline)
View(ped)

#write.table(ped, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/pedigree.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


#### make a manual map ### 
genotypes <- read.table("/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/affx_minDif0.7_genotypes_250724_linux.txt", header = TRUE, row.names = NULL)
View(genotypes)
nrow(genotypes)
ch1 <- sum(substr(genotypes$row.names, 6, 6) == 1)
#[1] 2195 JGH 2167, 1187
ch2 <- sum(substr(genotypes$SNP, 6, 6) == 2)
#[1] 2381 JGH 2318, 1298
ch3 <- sum(substr(genotypes$SNP, 6, 6) == 3)
#[1] 1766 JGH 1736, 912
ch4 <- sum(substr(genotypes$SNP, 6, 6) == 4)
#[1] 1588 JGH 1509, 820
ch5 <- sum(substr(genotypes$SNP, 6, 6) == 5)
# 2080 JGH 2113, 1140
ch6 <- sum(substr(genotypes$SNP, 6, 6) == 6)
# 1702 JGH 1699, 901
ch7 <- sum(substr(genotypes$SNP, 6, 6) == 7)
# 1419 JGH 424, 786
# 1564 JGH 1567, 840
ch9 <- sum(substr(genotypes$SNP, 6, 6) == 9)
# 1586 JGH 1548, 798

# not sure what above lines are doing but pretty sure this is the only line needed for making the manual map
linkage_groups <- substr(genotypes$row.names, 6, 6)

write.table(linkage_groups, "/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/minDif0.7_grouped_run_manual_map.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

map1 <- read.table("/Users/juliaharencar/Documents/Github/alle_vill_QTL/lepmap3/map1.txt")
head(map1)
unique(map1$V1)
