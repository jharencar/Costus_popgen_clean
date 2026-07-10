#!/bin/bash
#SBATCH --partition 128x24
#SBATCH --job-name=aHMM
#SBATCH --output=aHMM_scaf${SLURM_ARRAY_TASK_ID}.out
#SBATCH -e aHMM_scaf${SLURM_ARRAY_TASK_ID}.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kuckele@ucsc.edu
#SBATCH --mem=124G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5,7-11

# make new output directory
mkdir chr${SLURM_ARRAY_TASK_ID}/

# change into new output directory
cd /hb/scratch/kuckele/Davis_data/all_fastq/deduped/trimmomatic/ancestryInfer6/chr${SLURM_ARRAY_TASK_ID}/

# Rename samples in VCF file   
module load bcftools
bcftools reheader --samples ../rename.samples.list -o renamed_LBGM_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz ../../LBGM_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz

# Unzip VCF file
gunzip renamed_LBGM_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz
   
# Convert VCF file to aHMM input   
module load python-3.6.2
python3 ../vcf2aHMM.py --vcf renamed_LBGM_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf --pop ../ahmm.pop --rate 1 --minGT 1 --minDif 0.1 --dist 1000 --minDP 1 1>chr${SLURM_ARRAY_TASK_ID}.panel
  
# Run ancestry_hmm   
module load miniconda3.9
source activate ancestry_hmm
ancestry_hmm -i chr${SLURM_ARRAY_TASK_ID}.panel -s ../ahmm_sample.info -a 2 0.5 0.5 -p 0 100000 0.5 -p 1 1 0.5 -e 1e-3
	



# ancestry_hmm parameters:
# --minGT : minimum mean rate of genotype calls each parent population (1)
# --minDP : minimum mean depth for admixed population (1)
# --rate : mean recombination rate estimate in cM/Mb (1)
# --dist : distance between sites 1000 a pretty good default/min (1000)
# --minDif : allele freq difference between parentals (approximates Fst) (0.1)


