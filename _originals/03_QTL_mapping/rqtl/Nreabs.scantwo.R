library(qtl)

## load cross object after filtering and est.rf ; note this cross has fewer marks (second round of more stringent dup removal) than in final fitqtl models - should check and see if those look different with the exact same geno probs and markers as here in tems of which peaks and interactions are significant
Costus <- readRDS(file = "./aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_least_dups_est.rf.rds")

# set covariates based on previous model selection
covariates_Nreabs_select <- pull.pheno(Costus, c("hybrid_index", "F21_238","F21_259","F21_232","F21_284","coh10","coh16","coh18","coh21","coh22","coh24","coh25","coh27","coh29","coh30","coh32","coh34"))

# Pre-calculate with larger step
Costus <- calc.genoprob(Costus, step = 5, error.prob = 0.01)

library(parallel)
Nreabs.scantwo.hk.perm <- scantwo(
  Costus,
  method = "hk",
  pheno.col = "N_reabsorp",
  n.perm = 1008,
  addcovar = covariates_Nreabs_select,
  n.cluster = 72,
  verbose = TRUE
)

summary(Nreabs.scantwo.hk.perm)
saveRDS(Nreabs.scantwo.hk.perm, "Nreabs.scantwo.hk.perm.rds")

#calculate penalties to be used in stepwiseqtl based on scantwo permutations
Nreabs_pen <- calc.penalties(Nreabs.scantwo.hk.perm, alpha = 0.05)

#this step took forever, so writing the results to csv files
write.csv(Nreabs_pen, file = "seed.dormancy_penalties.csv")

#Nreabs_pen <- read.csv(file = "seed.dormancy_penalties.csv", header = TRUE, row.names = 1)

#stepwiseqtl: Stepwise selection for multiple QTL
#seed.dormancy
Nreabs.qtl <- stepwiseqtl(Costus, pheno.col = "log_days_to_radicle", method = "hk", model = "normal", covar=covariates_Nreabs_select, penalties = Nreabs_pen, max.qtl = 8, keeplodprofile = TRUE, verbose = TRUE, refine.locations = TRUE)
summary(Nreabs.qtl)

#fitqtl: Fit a multiple QTL model
Nreabs.fit <- fitqtl(Costus, pheno.col = "log_days_to_radicle", method = "hk", model = "normal", covar=covariates_Nreabs_select, get.ests = TRUE, qtl = Nreabs.qtl, formula = attr(Nreabs.qtl,"formula"), dropone = TRUE)
summary(Nreabs.fit)
