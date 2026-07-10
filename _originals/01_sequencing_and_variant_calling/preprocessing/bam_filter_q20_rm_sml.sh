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
#SBATCH --array=1-161  # []%3 limits it so only 3 run in parallel at once
#submit with SBATCH array, submit as normal or to specify certain array numbers: sbatch --array=x-y job_script.sbatch

# Code to remove all reads less than 100bp in length with samtools view and awk
#by: Julia Harenčár, 210817

## load samtools
module load samtools/samtools-1.10

## change into the directory with .bams
cd /hb/scratch/YOUR/PATH/HERE

# set variables
BAM_FILE=$(ls *merged.bam | sed -n ${SLURM_ARRAY_TASK_ID}p) # make list of file names to iterate through
SAMPLE_ID=$(echo $BAM_FILE |cut  -d "_" -f 1,2)

# indexing/sorting bams
samtools index ${BAM_FILE}
echo " samtools index ${BAM_FILE}"

### filtering out small reads ###
# -h in samtools view: export SAM file headers
# In awk, the substr function is used to keep header lines, and the other two conditions specify forward and reverse reads with the desired insert sizes, respectively
# -q 20 removes anything with mapping quality below 20
# Last but not least, samtools view -b is called to write the filtered reads into new a BAM file.

## NOTE: this is a still untested combo of removing small and filtering for lower than q20 - I am not 100% sure the filtering and writing will work correctly 
samtools view -h ${BAM_FILE} | \
  awk 'substr($0,1,1)=="@" || ($9>=35) || ($9<=-35)' | \
  samtools view -q 20 -b > ${SAMPLE_ID}_q20_no_sml.bam
