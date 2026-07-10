#!/bin/bash
#SBATCH --partition 512x64-ib 
#SBATCH --qos=512x64
#SBATCH --account=512x64
#SBATCH --job-name=aHMM_240904
#SBATCH --output=aHMM_240904_%A_%a.out
#SBATCH -e aHMM_240304_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --array=1-9

# make new output directory
#mkdir ALLEsub_final_params_refV2_chr${SLURM_ARRAY_TASK_ID}/

# change into new output directory
cd /hb/groups/kay_lab/Julia_popgen/aHMM/ALLEsub_final_params_refV2_chr${SLURM_ARRAY_TASK_ID}/

# Rename samples in VCF file   
#module load bcftools
#bcftools reheader --samples ../AHV_vcf_rename.txt -o ../renamed_AHV_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz ../AHV_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Unzip VCF file
#gunzip ../renamed_AHV_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz
   
# Convert VCF file to aHMM input   
module load python-3.6.2
python3 ../vcf2aHMM.py --vcf ../renamed_AHV_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf --pop ../aHMM_pop_subALLE_parents_oldstyle_240409.pop --rate 0.7 --minGT 1 --minDif 0.3 --dist 1000 --minDP 1.5 1>chr${SLURM_ARRAY_TASK_ID}.panel
  
# Run ancestry_hmm   
module load miniconda3.9
source activate JGH_conda_env
/hb/home/jharenca/Software/Ancestry_HMM/src/ancestry_hmm -i chr${SLURM_ARRAY_TASK_ID}.panel -s ../AHV_240409_ahmm.ploidy -a 2 0.5 0.5 -p 0 100000 0.5 -p 1 0.0001 0.5 -e 0.02
#/hb/home/jharenca/Software/Ancestry_HMM/src/ancestry_hmm -i chr${SLURM_ARRAY_TASK_ID}final_params240304_filtered.panel -s ../fixed_ahmm.ploidy -a 2 0.5 0.5 -p 0 100000 0.5 -p 1 0.0001 0.5 -e 0.02

# ancestry_hmm parameters:
# --minGT : minimum mean rate of genotype calls each parent population (1)
# --minDP : minimum mean depth for admixed population (1)
# --rate : mean recombination rate estimate in cM/Mb (1)
# --dist : distance between sites 1000 a pretty good default/min (1000)
# --minDif : allele freq difference between parentals (approximates Fst) (0.1)
