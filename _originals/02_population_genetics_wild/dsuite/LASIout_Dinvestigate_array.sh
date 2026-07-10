#!/bin/bash
#SBATCH --partition 512x64-ib 
#SBATCH --qos=ib
#SBATCH --account=ib
#SBATCH --job-name=Dinvestigate_w50s25_LASIout
#SBATCH --output=Dinvestigate_w50s25_LASIout.%A_%a.out
#SBATCH -e Dinvestigate_w50s25_LASIout.%A_%a.err
#SBATCH --mail-type=BEGAN,END,FAIL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --array=1-9

# move into working directory
cd /hb/groups/kay_lab/Julia_popgen/Dsuite

# run Dsuite Dinvestigate to get f_dM
# -w SIZE,STEP --window=SIZE,STEP         (required) D, f_D, f_dM, and d_f statistics for windows containing SIZE useable SNPs, moving by STEP (default: 50,25) 
# Maybe I should try multiple window sizes? What are resonable sizes?
# -n, --run-name                          run-name will be included in the output file name

echo "Dsuite Dinvestigate -w 50,25 -n LASIout_chr${SLURM_ARRAY_TASK_ID} /hb/groups/kay_lab/Julia_popgen/VCFs/BAVLasi_LASIref_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz LASIout_SETS.txt test_trios.txt"
Dsuite Dinvestigate -w 50,25 -n LASIout_chr${SLURM_ARRAY_TASK_ID} /hb/groups/kay_lab/Julia_popgen/VCFs/BAVLasi_LASIref_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz LASIout_SETS.txt test_trios.txt
echo "scaffold set ${SLURM_ARRAY_TASK_ID} done"
