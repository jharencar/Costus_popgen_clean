#!/bin/bash
#SBATCH --partition 128x24
#SBATCH --job-name=vcf_unfiltered
#SBATCH --output=vcf_unfiltered.out
#SBATCH -e vcf_unfiltered.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

#load bcftools and samtools
module load miniconda3.9
source activate ET-conda-env
module load samtools

# change into the directory with .bams
cd /hb/scratch/YOUR/PATH/HERE

# set variables
BAM_LIST=$(ls *.bam) # make list of file names to iterate through

# generate index .bai files for each .bam in your directory
for BAM_NAME in $BAM_LIST; do
	samtools index ${BAM_NAME}
done

## Code to generate initial vcf in Nicolas Alexandre's snp filtering 3 step code; code originally written by Tyler Linderoth
# modified by Julia Harenčár 211030

# change into the directory with bam.filelist (explicit paths in list, of in same directory as bams)
cd /hb/scratch/YOUR/PATH/HERE

BAM_LIST="YOUR.BAM.LIST.bam.filelist"
REF="/hb/groups/kay_lab/genomes/C.lasius_genome_NCBI.fasta"
THREADS=48

# step 1: for each of 21 sets of scaffolds
for i in {1..21}; do
        bcftools mpileup -b ${BAM_LIST} -R scaffolds_${i}.list -f ${REF} --threads ${THREADS} -d 100000000 -Q 0 -a AD,DP,SP,INFO/ADF,INFO/ADR -L 100000000 | bcftools call --threads ${THREADS} -c -A -f GQ - -O z > unfiltered_scaffold_set${i}.vcf.gz
done

for i in {1..21}; do
        echo "bcftools mpileup -b ${BAM_LIST} -R scaffolds_${i}.list -f ${REF} --threads ${THREADS} -d 100000000 -Q 0 -a AD,DP,SP,INFO/ADF,INFO/ADR -L 100000000 | bcftools call --threads ${THREADS} -c -A -f GQ - -O z > unfiltered_scaffold_set${i}.vcf.gz"
done > job.list



