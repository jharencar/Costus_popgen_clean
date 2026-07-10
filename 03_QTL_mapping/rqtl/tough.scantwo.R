library(qtl)

## load cross object after filtering and est.rf ; note this cross has fewer marks (second round of more stringent dup removal) than in final fitqtl models - should check and see if those look different with the exact same geno probs and markers as here in tems of which peaks and interactions are significant
Costus <- readRDS(file = "./aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_least_dups_est.rf.rds")

# set covariates based on previous model selection
covariates_tough_select <- pull.pheno(Costus, c("hybrid_index", "F21_232","F21_238","F21_259","F21_284","coh10","coh11","coh16","coh19","coh22","coh24","coh27","coh28","coh29","coh30","coh32"))

# Pre-calculate with larger step
Costus <- calc.genoprob(Costus, step = 5, error.prob = 0.01)

library(parallel)
tough.scantwo.hk.perm <- scantwo(
  Costus,
  method = "hk",
  pheno.col = "ave_tough",
  n.perm = 1008,
  addcovar = covariates_tough_select,
  n.cluster = 72,
  verbose = TRUE
)

summary(tough.scantwo.hk.perm)
saveRDS(tough.scantwo.hk.perm, "tough.scantwo.hk.perm.rds")

#calculate penalties to be used in stepwiseqtl based on scantwo permutations
tough_pen <- calc.penalties(tough.scantwo.hk.perm, alpha = 0.05)

#this step took forever, so writing the results to csv files
write.csv(tough_pen, file = "seed.dormancy_penalties.csv")

#tough_pen <- read.csv(file = "seed.dormancy_penalties.csv", header = TRUE, row.names = 1)

#stepwiseqtl: Stepwise selection for multiple QTL
#seed.dormancy
tough.qtl <- stepwiseqtl(Costus, pheno.col = "log_days_to_radicle", method = "hk", model = "normal", covar=covariates_tough_select, penalties = tough_pen, max.qtl = 8, keeplodprofile = TRUE, verbose = TRUE, refine.locations = TRUE)
summary(tough.qtl)

#fitqtl: Fit a multiple QTL model
tough.fit <- fitqtl(Costus, pheno.col = "log_days_to_radicle", method = "hk", model = "normal", covar=covariates_tough_select, get.ests = TRUE, qtl = tough.qtl, formula = attr(tough.qtl,"formula"), dropone = TRUE)
summary(tough.fit)
