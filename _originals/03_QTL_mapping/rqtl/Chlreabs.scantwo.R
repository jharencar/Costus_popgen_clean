library(qtl)

## load cross object after filtering and est.rf ; note this cross has fewer marks (second round of more stringent dup removal) than in final fitqtl models - should check and see if those look different with the exact same geno probs and markers as here in tems of which peaks and interactions are significant
Costus <- readRDS(file = "./aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_least_dups_est.rf.rds")

# set covariates based on previous model selection
covariates_Chlreabs_select <- pull.pheno(Costus, c("hybrid_index", "F21_238", "F21_284", "F21_232", "coh13", "coh16", "coh19", "coh21", "coh22", "coh24", "coh25", "coh27", "coh28", "coh30", "coh31", "coh32"))

# Pre-calculate with larger step
Costus <- calc.genoprob(Costus, step = 5, error.prob = 0.01)

library(parallel)
Chlreabs.scantwo.hk.perm <- scantwo(
  Costus,
  method = "hk",
  pheno.col = "chloro_reabsorp",
  n.perm = 1008,
  addcovar = covariates_Chlreabs_select,
  n.cluster = 72,
  verbose = TRUE
)

summary(Chlreabs.scantwo.hk.perm)
saveRDS(Chlreabs.scantwo.hk.perm, "Chlreabs.scantwo.hk.perm.rds")

#calculate penalties to be used in stepwiseqtl based on scantwo permutations
Chlreabs_pen <- calc.penalties(Chlreabs.scantwo.hk.perm, alpha = 0.05)

#this step took forever, so writing the results to csv files
write.csv(Chlreabs_pen, file = "seed.dormancy_penalties.csv")

#Chlreabs_pen <- read.csv(file = "seed.dormancy_penalties.csv", header = TRUE, row.names = 1)

#stepwiseqtl: Stepwise selection for multiple QTL
#seed.dormancy
Chlreabs.qtl <- stepwiseqtl(Costus, pheno.col = "log_days_to_radicle", method = "hk", model = "normal", covar=covariates_Chlreabs_select, penalties = Chlreabs_pen, max.qtl = 8, keeplodprofile = TRUE, verbose = TRUE, refine.locations = TRUE)
summary(Chlreabs.qtl)

#fitqtl: Fit a multiple QTL model
Chlreabs.fit <- fitqtl(Costus, pheno.col = "log_days_to_radicle", method = "hk", model = "normal", covar=covariates_Chlreabs_select, get.ests = TRUE, qtl = Chlreabs.qtl, formula = attr(Chlreabs.qtl,"formula"), dropone = TRUE)
summary(Chlreabs.fit)
