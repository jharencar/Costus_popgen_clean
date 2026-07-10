#!/bin/bash
#SBATCH --partition 128x24 
#SBATCH --job-name=aHMM
#SBATCH -D /hb/scratch/jharenca/QTL_mapping/aHMM/run_all_single_samp_full_ALLEref_VCF
#SBATCH --output=/hb/scratch/jharenca/QTL_mapping/aHMM/run_all_single_samp_full_ALLEref_VCF/slurm_log/aHMM_.%A_%a.out
#SBATCH -e /hb/scratch/jharenca/QTL_mapping/aHMM/run_all_single_samp_full_ALLEref_VCF/slurm_log/aHMM_.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-425

#load modules
module load gnu12/12.2.0 bcftools
module load py3-numpy

# Get the sample name for each array job from the sample list
SAMPLE_LIST="hyb_samps.list"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

# Create a directory for the sample
mkdir -p "${SAMPLE}" || { echo "Failed to create directory for $SAMPLE"; exit 1; }

# Change into the sample directory
cd "${SAMPLE}" || { echo "Failed to change directory to $SAMPLE"; exit 1; }

# Loop through all 9 chromosomes
for CHR in {1..9}; do
    # Create a subdirectory for each chromosome
    mkdir -p chr${CHR} || { echo "Failed to create directory for chr${CHR}"; exit 1; }

    # Change into the chromosome directory
    cd chr${CHR} || { echo "Failed to change directory to chr${CHR}"; exit 1; }

    # Convert VCF file to aHMM input
    python3 ../../../new_vcf2ahmm.py --vcf ../../renamed_aHMMsamps_ALLEref_filtered_chr${CHR}.vcf --pop ../../${SAMPLE}.pop --rate 0.7 --minGTP 2 --minDif 0.1 --dist 1000 --minDP 1 1>chr${CHR}.panel || { echo "python3 script failed for $SAMPLE, chr${CHR}"; exit 1; }

    echo "Panel generation done for $SAMPLE, chr${CHR}"

    # Run ancestry_hmm
    /hb/home/jharenca/Software/Ancestry_HMM/src/ancestry_hmm -i chr${CHR}.panel -s ../../${SAMPLE}.ploidy -a 2 0.5 0.5 -p 0 100000 0.5 -p 1 0.0000001 0.5 -e 0.02 || { echo "ancestry_hmm failed for $SAMPLE, chr${CHR}"; exit 1; }

    echo "Ancestry HMM done for $SAMPLE, chr${CHR}"

    # Move back to the sample directory
    cd ..
done

echo "All chromosomes processed for $SAMPLE."