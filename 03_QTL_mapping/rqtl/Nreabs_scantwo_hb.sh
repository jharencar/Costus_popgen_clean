#!/bin/sh
#SBATCH --partition 512x64-ib 
#SBATCH --qos=ib
#SBATCH --account=ib
#SBATCH --exclusive=user
#SBATCH --job-name=rqtl_scantwo_Nreabs
#SBATCH --output=rqtl_scantwo_Nreabs.out
#SBATCH -e rqtl_scantwo_Nreabs.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --mem=1000G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72 # max is 63

# load R QTL
module load miniconda3
conda activate /hb/groups/kay_lab/R_QTL

cd /hb/groups/kay_lab/Julia_popgen/rqtl

Rscript ./Nreabs.scantwo.R
