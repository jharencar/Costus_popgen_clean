#!/bin/bash
#SBATCH --partition 128x24 # this means each node (there are 19 nodes in this partition) has 128 GB RAM and 24 cores/cpus 
#SBATCH --job-name=filter_samps # Job name
#SBATCH --time=72:00:00 # adjust once I know how long it takes
#SBATCH --output=filter_samps.out
#SBATCH -e filter_samps.err
#SBATCH --mail-type=END, FAIL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jharenca@ucsc.edu  # Where to send mail
#SBATCH --mem=6G 
#SBATCH --cpus-per-task=12
#SBATCH --array=1-9 # set to number of chromosomes

module load bcftools

bcftools view -S AV_only.samps.file AHVF1_allSites_filtered_chr${SLURM_ARRAY_TASK_ID}.vcf.gz > AV_only_filtered_chr${SLURM_ARRAY_TASK_ID}.vcf

bgzip -c AV_only_filtered_chr${SLURM_ARRAY_TASK_ID}.vcf > AV_only_filtered_chr${SLURM_ARRAY_TASK_ID}.vcf.gz 

bcftools index AV_only_filtered_chr${SLURM_ARRAY_TASK_ID}.vcf.gz


pixy --stats pi fst dxy \
--vcf AV_only_filtered_chr${SLURM_ARRAY_TASK_ID}.vcf.gz \
--populations AV_sampleIDs_popfile.txt \
--window_size 1000 \
--n_cores 12 \
--chromosomes '${SLURM_ARRAY_TASK_ID}'