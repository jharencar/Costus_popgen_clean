#!/bin/bash
#SBATCH --partition 512x64-ib 
#SBATCH --qos=ib
#SBATCH --account=ib
#SBATCH -D /hb/scratch/jharenca/QTL_mapping/aHMM/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO
#SBATCH --output=/hb/scratch/jharenca/QTL_mapping/aHMM/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO/slurm_log/aHMM_.%A_%a.out
#SBATCH -e /hb/scratch/jharenca/QTL_mapping/aHMM/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO/slurm_log/aHMM_.%A_%a.err
#SBATCH --job-name=aHMM
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-9

#load modules
module load gnu12/12.2.0 bcftools
module load py3-numpy

# make new output directory
mkdir chr${SLURM_ARRAY_TASK_ID}/ || { echo "Failed to create directory"; exit 1; }

# change into new output directory
cd chr${SLURM_ARRAY_TASK_ID}/ || { echo "Failed to change directory"; exit 1; }

# Convert VCF file to aHMM input   
python3 ../../new_vcf2ahmm.py --vcf ../../LASIref_renamed_QTL_w.F1.F2_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf --pop ../F3s0.05covCO_allF2s_F1parent.pop --rate 0.7 --minGTP 3 --minDif 0.7 --dist 5000 --minDP 0.5 1>chr${SLURM_ARRAY_TASK_ID}.panel || { echo "python3 script failed"; exit 1; }
echo "done with panel gen"

## --minGTP minimum number of genotype calls for each parent population (default = 20)
## --minDP minimum rate of depth>1 for admixed population (if no -geno flag)

# Run ancestry_hmm
/hb/home/jharenca/Software/Ancestry_HMM/src/ancestry_hmm -i chr${SLURM_ARRAY_TASK_ID}.panel -s ../F3s0.05covCO_allF2s_F1parent.ploidy -a 2 0.5 0.5 -p 0 100000 0.5 -p 1 0.0001 0.5 -e 0.02 || { echo "ancestry_hmm failed"; exit 1; }
echo "done"
