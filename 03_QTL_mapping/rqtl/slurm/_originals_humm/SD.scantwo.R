library(qtl)

#Costus <- read.cross("csvs", ".", genfile = "F3gen.csv", phefile = "F3phe_onehotcovar.csv", estimate.map=FALSE)
#Costus <- subset(Costus, ind=(ntyped(Costus)>11607))
## identify names of markers to drop, based on genotypes-per-marker plot
# drop all markers with less than ~80% genotyped individuals (391*0.8 = 313)
#nt.bymar <- ntyped(Costus, "mar")
#todrop <- names(nt.bymar[nt.bymar < 313]) 
## drop those markers
#Costus <- drop.markers(Costus, todrop)

#dup <- findDupMarkers(Costus, exact.only=FALSE, adjacent.only=TRUE)
#Costus <- drop.markers(Costus, unlist(dup))

#estimate the recombination fraction 
#Costus <- est.rf(Costus)
#saveRDS(Costus, file="./Costus_AV_qtl_data_est.rfdup_dropped.dup.markers.rds")

## calculate conditional genotype probabilities given multipoint marker data
#Costus_prob <- calc.genoprob(Costus)

## load est.rf results after running on cluster/serv <- readRDS(file = "./Costus_AV_qtl_data_est.rf.rds")
Costus <- readRDS(file = "./aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_est.rf.rds")

covariates_SD_select <- pull.pheno(Costus, c("s.coh4", "s.coh7", "s.coh8", "s.coh9", "s.coh14", "s.coh15"))

# Pre-calculate with larger step
Costus <- calc.genoprob(Costus)

library(parallel)
SD.scantwo.hk.perm <- scantwo(
  Costus,
  method = "hk",
  pheno.col = "log_days_to_radicle",
  n.perm = 1000,
  addcovar = covariates_SD_select,
  n.cluster = 67,
  verbose = TRUE
)

#SD.scantwo.hk.perm_full_set <- scantwo(Costus, method="hk", pheno.col = "log_days_to_radicle", n.cluster = 30, n.perm=990, addcovar=covariates_SD_select, verbose = TRUE)
summary(SD.scantwo.hk.perm)
saveRDS(SD.scantwo.hk.perm, "SD.scantwo.hk.perm.rds")

#calculate penalties to be used in stepwiseqtl based on scantwo permutations
SD_pen <- calc.penalties(SD.scantwo.hk.perm, alpha = 0.05)

#this step took forever, so writing the results to csv files
write.csv(SD_pen, file = "seed.dormancy_penalties.csv")

#SD_pen <- read.csv(file = "seed.dormancy_penalties.csv", header = TRUE, row.names = 1)

#stepwiseqtl: Stepwise selection for multiple QTL
#seed.dormancy
SD.qtl <- stepwiseqtl(Costus, pheno.col = "log_days_to_radicle", method = "hk", model = "normal", covar=covariates_SD_select, penalties = SD_pen, max.qtl = 8, keeplodprofile = TRUE, verbose = TRUE, refine.locations = TRUE)
summary(SD.qtl)

#fitqtl: Fit a multiple QTL model
seed.dormancy.fit <- fitqtl(Costus, pheno.col = "log_days_to_radicle", method = "hk", model = "normal", covar=covariates_SD_select, get.ests = TRUE, qtl = SD.qtl, formula = attr(SD.qtl,"formula"), dropone = TRUE)
summary(seed.dormancy.fit)
