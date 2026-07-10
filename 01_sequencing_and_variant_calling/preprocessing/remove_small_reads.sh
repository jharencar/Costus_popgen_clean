#!/bin/bash
#SBATCH --partition 128x24 # this means each node (there are 19 nodes in this partition) has 128 GB RAM and 24 cores/cpus 
#SBATCH --job-name=rm_sml# Job name
#SBATCH --time=24:00:00 
#SBATCH --output=rm_sml_.%A_%a.out
#SBATCH -e rm_sml_.%A_%a.err
#SBATCH --mem=105G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --array=[1-156]%3  # []%3 limits it so only 3 run in parallel at once
#submit with SBATCH array, submit as normal or to specify certain array numbers: sbatch --array=x-y job_script.sbatch

# Code to remove all reads less than 100bp in length with samtools view and awk
#by: Julia Harenčár, 210817

## load samtools
module load samtools/samtools-1.10

## change into the directory with .bams
cd /hb/scratch/jharenca/pop_gen/combinedP1P2/merged_bams

# set variables
BAM_FILE=$(ls *q20.bam | sed -n ${SLURM_ARRAY_TASK_ID}p) # make list of file names to iterate through
SAMPLE_ID=$(echo $BAM_FILE |cut  -d "_" -f 1,2)

### filtering out small reads ###
# -h in samtools view: export SAM file headers
# In awk, the substr function is used to keep header lines, and the rest two condition specify forward and reverse reads with the desired insert sizes, respectively
# Last but not least, samtools view -b is called to write the filtered reads into new a BAM file.

echo " samtools view -h ${BAM_FILE} | \
  awk 'substr($0,1,1)=="@" || ($9>=100) || ($9<=-100)' | \
  samtools view -b > ${SAMPLE_ID}_no_sml.bam "

samtools view -h ${BAM_FILE} | \
  awk 'substr($0,1,1)=="@" || ($9>=100) || ($9<=-100)' | \
  samtools view -b > ${SAMPLE_ID}_no_sml.bam
