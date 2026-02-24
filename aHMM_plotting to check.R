#Code for Generating a plot of Ancestry_HMM output
#Julia Harenčár
#

library(readr)
library("ggplot2")
library(dplyr)
theme_set(theme_bw())
library(cowplot)

## clear workspace
dev.off()
rm(list = ls())

# popgen wd
# setwd("/Users/juliaharencar/Documents/Github/Population_genetics_vill_alle/aHMM/final_params_0.02er_0.0001pulse")

# QTL wd
#setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/full_final_posteriors/")
#setwd("/Users/juliaharencar/Documents/Github/alle_vill_QTL/aHMM/smlBatchTest_DP1_minDif0.2")
setwd("/Users/juliaharencar/Documents/Github/Population_genetics_vill_alle/aHMM/pre_filter_DP0.5_minDif0.3")

directory <- "./chr3/"
chromosome <- "chr3"

files <- list.files(directory, pattern = "*posterior")
#files <- files[grep("HYB", files)]
IDs <- substr(files, 1, 9)

plot_list = list()
for (i in 1:length(files)) {
  post <- read_tsv(paste0(directory, files[i]))
  # rename cols 
  names(post) <- c("chrom", "position", "ref_homo", "het", "alt_homo")
  # Filter out posteriors with a 0.9 threshold
  post <- post %>% filter(ref_homo > 0.9 | alt_homo > 0.9 | het > 0.9) 
  # Create allele frequency cols for plotting 
  post <- post %>% mutate (ancestry = ref_homo+het*0.5)
  p <- ggplot(post, aes(position, ancestry)) + 
    ggtitle(paste(IDs[i], chromosome)) +
    geom_point() +
    ylim(-0.01,1.01)
  plot_list[[i]] = p
}

pdf("chr3.popgen_F1s_DP0.5_minDif0.3.pdf", height = 10, width = 12)
 # looking at just 15 from mapping pop
plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
          plot_list[[5]],
          nrow = 3,
          ncol = 2)
dev.off()

