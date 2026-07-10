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
Costus <- readRDS(file = "./aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_least_dups_est.rf.rds")

covariates_GR_select <- pull.pheno(Costus,c("hybrid_index", "F21_259", "F21_238", "F21_284", "F21_232","coh7","coh10","coh11", "coh13","coh18","coh19","coh24","coh25","coh27","coh28","coh29","coh30","coh31","coh32","coh34"))

# Pre-calculate with larger step
Costus <- calc.genoprob(Costus, step = 5, error.prob = 0.01)

library(parallel)
GR.scantwo.hk.perm <- scantwo(
  Costus,
  method = "hk",
  pheno.col = "GR_cm.dy",
  n.perm = 1000,
  addcovar = covariates_GR_select,
  n.cluster = 72,
  verbose = TRUE
)

#SD.scantwo.hk.perm_full_set <- scantwo(Costus, method="hk", pheno.col = "log_days_to_radicle", n.cluster = 30, n.perm=990, addcovar=covariates_SD_select, verbose = TRUE)
summary(GR.scantwo.hk.perm)
saveRDS(GR.scantwo.hk.perm, "GR.scantwo.hk.perm.rds")

#calculate penalties to be used in stepwiseqtl based on scantwo permutations
GR_pen <- calc.penalties(GR.scantwo.hk.perm, alpha = 0.05)

#this step took forever, so writing the results to csv files
write.csv(GR_pen, file = "growth.rate_penalties.csv")

#SD_pen <- read.csv(file = "seed.dormancy_penalties.csv", header = TRUE, row.names = 1)

#stepwiseqtl: Stepwise selection for multiple QTL
#seed.dormancy
GR.qtl <- stepwiseqtl(Costus, pheno.col = "GR_cm.dy", method = "hk", model = "normal", covar=covariates_GR_select, penalties = as.vector(GR_pen["GR_cm.dy",1:3]), max.qtl = 8, keeplodprofile = TRUE, verbose = TRUE, refine.locations = TRUE)
summary(GR.qtl)

#fitqtl: Fit a multiple QTL model
GR.fit <- fitqtl(Costus, pheno.col = "GR_cm.dy", method = "hk", model = "normal", covar=covariates_GR_select, get.ests = TRUE, qtl = GR.qtl, formula = attr(GR.qtl,"formula"), dropone = TRUE)
summary(GR.fit)
