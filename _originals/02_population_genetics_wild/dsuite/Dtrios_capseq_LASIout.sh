#!/bin/bash
#SBATCH --partition 128x24
#SBATCH --job-name=capseq_Dtrios
#SBATCH --output=capseq_Dtrios.out
#SBATCH -e capseq_Dtrios.err
#SBATCH --mail-type=BEGAN,END,FAIL
#SBATCH --mail-user=jharenca@ucsc.edu

# move into working directory
cd /hb/groups/kay_lab/Julia_popgen/Dsuite

# run Dsuite 
# -k, --JKnum       (default=20) the number of Jackknife blocks to divide the dataset into; should be at least 20 for the whole dataset
# -j, --JKwindow    (default=NA) Jackknife block size in number of informative SNPs (as used in v0.2) when specified, this is used in place of the --JKnum option
# play with/research and think about -k val or using -j
# starting by running defaults

#-t, --tree=TREE_FILE.nwk  (optional) a file with a tree in the newick format specifying the relationships between populations/species
#                          D and f4-ratio values for trios arranged according to the tree will be output in a file with _tree.txt suffix
# I think I probably want to do this - include a tree with the known relationships - look at example tutorial and maybe oother pprs. 

# -l NUMLINES           (optional) the number of lines (SNPs) in the VCF
#                       input(s) - speeds up operation if known. If N
#                       INPUT_FILE.vcf files are provided, there must be N
#                       comma-separated integers provided without whitespace
#                       between them.
# Should definitely figure this out and do this eventually, especially after filtering (might take to long without any filters and spoecifying... )
# 

Dsuite Dtrios -k 20 -t tree.nwk --out-prefix=capseq_LASIout --run-name TEST no19113_cat_allclusters_supercontig.vcf.biallelic.maaft.trim.spur.ts.rn4.vcf capseq_LASIout_SETS.txt 
