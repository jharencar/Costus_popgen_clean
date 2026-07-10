#!/bin/bash
#SBATCH --partition 128x24 # this means each node (there are 19 nodes in this partition) has 128 GB RAM and 24 cores/cpus 
#SBATCH --job-name=chrs_abba_baba2 # Job name
#SBATCH --time=72:00:00 # adjust once I know how long it takes
#SBATCH --output=chrs_abba_baba2.out
#SBATCH -e chrs_abba_baba2.err
#SBATCH --mail-type=END, FAIL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jharenca@ucsc.edu  # Where to send mail
#SBATCH --mem=6G 
#SBATCH --cpus-per-task=12
#SBATCH --array=1-9 # set to number of chromosomes

# load modules
module load angsd
module load R/R-3.6.1 miniconda3.9 # extra R packages loaded in miniconda JGH env activated below
source activate JGH_conda_env

# change directories 
cd /hb/groups/kay_lab/Julia_popgen/ANGSD/abba_baba 

# run doabbababba2
angsd -doAbbababa2 1 -bam BAV_clean_ordered_0.15co_240205.bam.filelist -sizeFile sizeFile.size -doCounts 1 -out lasiOUT_chr${SLURM_ARRAY_TASK_ID}_240219 -anc /hb/groups/kay_lab/genomes/C.lasius_genome_NCBI.fasta -minQ 30 -minMapQ 30 -only_proper_pairs 1 -baq 1 -ref /hb/groups/kay_lab/genomes/C.lasius_genome_NCBI.fasta -blockSize 5000000 -nThreads 96 -r Chrom${SLURM_ARRAY_TASK_ID}

# calculate the Dstat with R script
#Rscript /hb/software/apps/angsd/gnu-937/angsd/R/estAvgError.R angsdFile="lasiOUT_chr${SLURM_ARRAY_TASK_ID}_240219" out="lasiOUT_chr${SLURM_ARRAY_TASK_ID}_240219" sizeFile=sizeFile.size nameFile=popNames.name
