#!/bin/bash
#SBATCH --partition 128x24   
#SBATCH --job-name=bwa_mapping # Job name
#SBATCH --time=7-24:00:00
#SBATCH --output=bwa_mapping.%A_%a.out
#SBATCH -e bwa_mapping.%A_%a.err
#SBATCH --mail-type=END,FAIL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jharenca@ucsc.edu  # Where to send mail
#SBATCH --nodes=1                    
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --array=[1-96]%3  # []%6 limits it so only 6 run in parallel at once
#submit with SBATCH array, submit as normal or to specify certain array numbers: sbatch --array=x-y job_script.sbatch

##ARRAY NOTES: Tasks can mean multiple things, depending on context.
##(more info on arrays: https://support.ceci-hpc.be/doc/_contents/SubmittingJobs/SlurmFAQ.html)
###Here, array over 3 nodes, 1 operation per node taking 24 cpus each (24 threads)

##get bwa and samtools modules
module load bwakit/bwa-0.7.15 samtools

##move into proper folder with all the .fq.gz files
cd /hb/scratch/jharenca/vill_alle_plate1/newdata/AVH_plate2/trimmed

# variable setting
NUMTHREADS=24
LIBRARY_NAME="plate2"
REF="C_lasius_ref.fa"

# Set up array variables so that the same SLURM_ARRAY_TASK_ID references both input files for a given individual (F and R reads)
R1=$(ls *_F_paired.fq.gz| sed -n ${SLURM_ARRAY_TASK_ID}p)
R2=$(ls *_R_paired.fq.gz| sed -n ${SLURM_ARRAY_TASK_ID}p)

# This is what I want: "19.572_LAEV" from 19.572_LAEV_F_paired.fq.gz 
SAMPLE_NAME=$(echo $R1 |cut  -d "_" -f 1,2 )

### NOTE: index reference before submitting SLURM: bwa index $REF 
# bwa mapping; -t number of treads; -M is a verboseLevel ; -R readGroupHeader ; 
echo "bwa mem -t ${NUMTHREADS} -M -R "@RG\tID:${LIBRARY_NAME}\tSM:${SAMPLE_NAME}\tLB:${LIBRARY_NAME}\tPL:ILLUMINA\tPU:none" ${REF} ${R1} ${R2} | samtools view -hb -@ ${NUMTHREADS} - | samtools sort -@ ${NUMTHREADS} - -o ${SAMPLE_NAME}_output.bam"
bwa mem -t ${NUMTHREADS} -M -R "@RG\tID:${LIBRARY_NAME}\tSM:${SAMPLE_NAME}\tLB:${LIBRARY_NAME}\tPL:ILLUMINA\tPU:none" ${REF} ${R1} ${R2} | samtools view -hb -@ ${NUMTHREADS} - | samtools sort -@ ${NUMTHREADS} - -o ${SAMPLE_NAME}_output.bam


