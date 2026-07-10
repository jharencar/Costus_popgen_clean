#!/bin/bash
#SBATCH --partition 512x64-ib 
#SBATCH --qos=ib
#SBATCH --account=ib
#SBATCH --job-name=aHMM
#SBATCH --output=filter_.%A_%a.out
#SBATCH -e aHMM_.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-9

# load modules
module load bcftools

# First ID shared heterozygous sites
# Take vcf for a chromosome, subset for 1 parent, then the other, filter for min allele frequency of 0.001 (sites must have at least)
bcftools view ALLE_only_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz -q .001 -Q 0.999 -H | cut -f 1-2 >VILL.heterozygous.txt
bcftools view VILL_only_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz -q .001 -Q 0.999 -H | cut -f 1-2 >ALLE.heterozygous.txt

# Compare sets of heterozygous sites from each parent and keep those shared between both
comm -12 <(sort VILL.heterozygous.txt) <(sort VILL.heterozygous.txt) | sort -nk 2 >shared.heterozygous.txt

# Next ID shared homozygous sites - will include those with the same allele but those will be removed later (still don't fully get this code)
# subset for parents, set max allele freq to 0.01, so will include all where allele 1 (max) is 0 (homo for min allele)
# next get the opposite homozygous sites, set min allele freq to 0.99
bcftools view VILL_only_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz -Q 0.01 -H | cut -f 1-2 >VILL.homozygous.txt
bcftools view VILL_only_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz -q 0.99 -H | cut -f 1-2 >>VILL.homozygous.txt
bcftools view ALLE_only_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz -Q 0.01 -H | cut -f 1-2 >ALLE.homozygous.txt
bcftools view ALLE_only_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz -q 0.99 -H | cut -f 1-2 >>ALLE.homozygous.txt

comm -12 <(sort VILL.homozygous.txt) <(sort VILL.homozygous.txt) | sort -nk 2  > shared.homozygous.txt

# combining the above sites we basically just remove sites homozygous in one and heterozygous in the other
cat shared.homozygous.txt shared.heterozygous.txt | sort -nk 2 >shared.sites.txt

