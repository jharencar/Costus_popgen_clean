#!/bin/bash
#SBATCH --partition 128x24 # this means each node (there are 19 nodes in this partition) has 128 GB RAM and 24 cores/cpus 
#SBATCH --job-name=ANGSD1_beagles# Job name
#SBATCH --time=24:00:00
#SBATCH --output=ANGSD1_beagles_output.txt
#SBATCH --mail-type=END, FAIL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jharenca@ucsc.edu  # Where to send mail
#SBATCH --nodes=1 

## load ANGSD
module load angsd/gnu-930 samtools/samtools-1.10

## change into the directory with .bams and script
cd /hb/scratch/jharenca/vill_alle_plate1/newdata/whole_plate_mapping/ANGSD

## call the code to run

# -GL - Estimate genotype likelihoods. Options - 1: SAMtools 2: GATK 3: SOAPsnp 4: SYK
# helpful summary of some of the main diffs between GATK and samtools: https://www.biostars.org/p/57149/
# Sounds like samtools is what I want as GATK is optimized for human genomics and I will want to manually set filters (I don't have piles of data to train GATK with to infer filters)
# -doGlf: glf is a simple format for genotype likelihoods. Options - 1: binary glf (10 log likes) .glf.gz 2: beagle likelihood file  .beagle.gz 3: binary 3 times likelihood	.glf.gz 4: text version (10 log likes)	.glf.gz
# -doMajorMinor: infer the major and minor alleles. Options - 1: Infer major and minor from GL (THIS is what I want) 2: Infer major and minor from allele counts 3: use major and minor from a file (requires -sites file.txt) 4: Use reference allele as major (requires -ref) 5: Use ancestral allele as major (requires -anc)
# -SNP_pval: p value to be used for snp calling, 1e-6 seems to be a good default
# -doMaf: estimate allele frequencices. Options - 0: (Calculate persite frequencies '.mafs.gz') 1: Frequency (fixed major and minor based on maj and minor inferred from -doMajorMinor 1) 2: Frequency (fixed major unknown minor) 4: Frequency from genotype probabilities 8: AlleleCounts based method (known major minor)

# to generate beagle files for NGSadmix and PCAngsd
angsd -GL 1 -out AVH_only.genolike -nThreads 24 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam AVH_only.bam.filelist