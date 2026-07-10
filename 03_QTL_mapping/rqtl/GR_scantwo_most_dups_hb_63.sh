#!/bin/sh
#SBATCH --partition 512x64-ib 
#SBATCH --qos=ib
#SBATCH --account=ib
#SBATCH --exclusive=user
#SBATCH --job-name=most_dups_rqtl_scantwo_GR63
#SBATCH --output=most_dups_rqtl_scantwo_GR63.out
#SBATCH -e most_dups_rqtl_scantwo_GR63.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --mem=511G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=63 # max is 63

# load R QTL
module load miniconda3
conda activate /hb/groups/kay_lab/R_QTL

cd /hb/groups/kay_lab/Julia_popgen/rqtl

Rscript ./GR.scantwo_most_dups_63.R
