#!/bin/bash
#SBATCH --partition 128x24 
#SBATCH --job-name=bcftools_filter
#SBATCH --output=bcftools_filter.%A_%a.out
#SBATCH --error=bcftools_filter.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=1-9
#SBATCH --mem=4G

# submit with SBATCH array, submit as normal or to specify certain array numbers: sbatch --array=x-y job_script.sbatch

# load bcftools and vcftools
module load vcftools gnu12/12.2.0 bcftools
module load py3-numpy

# move into directory with VCFs needing filtering
cd /hb/groups/kay_lab/Julia_popgen/VCFs

# set variables
VCF_FILE="popgen_ALLEref_unfiltered_scaffold_set${SLURM_ARRAY_TASK_ID}.vcf.gz"
SCAFFOLD_SET=${SLURM_ARRAY_TASK_ID}
OUT="popgen_ALLEref_filtered_scaf_${SLURM_ARRAY_TASK_ID}"
OUTPUT_VCF="popgen_ALLEref_filtered_scaf_${SLURM_ARRAY_TASK_ID}.vcf.gz"

# generate data file with per site mean depths
vcftools --gzvcf "$VCF_FILE" --site-mean-depth --out "$OUT"

# Create a temporary file to store mean depth values
tmp_mean_depth_file=$(mktemp)

# Extract mean depth values using awk and save to the temporary file
awk 'NR>1 {print $3}' "${OUT}.ldepth.mean" > "$tmp_mean_depth_file"

# Calculate the 99th quantile using Python
quantile_99=$(python -c "
import numpy as np
depths = np.loadtxt('${tmp_mean_depth_file}')
quantile_99 = np.quantile(depths, 0.99)
print(quantile_99)
")

# Clean up temporary file
rm "$tmp_mean_depth_file"

# Calculate max depth threshold
num_indivs=$(bcftools query -l "$VCF_FILE" | wc -l)
max_depth=$(awk "BEGIN {print ${quantile_99} * ${num_indivs}}")
echo "max depth is ${max_depth}"

# some filtering/cleaning
condition="INFO/DP<=$max_depth"

echo "bcftools filter --SnpGap 3 -i 'QUAL>19 && INFO/MQ>19' --threads 2 -O z $VCF_FILE | \
bcftools view -T ^/hb/groups/kay_lab/genomes/genmap/ALLE_non_callable_sites.bed -m2 -M2 -i \"$condition\" -v snps -U --threads 2 -q .05 -Q .95 -o ${OUTPUT_VCF} -O z"

bcftools filter --SnpGap 3 -i 'QUAL>19 && INFO/MQ>19' --threads 2 -O z "$VCF_FILE" | \
bcftools view -T ^/hb/groups/kay_lab/genomes/genmap/ALLE_non_callable_sites.bed -m2 -M2 -i "$condition" -v snps -U --threads 2 -q .05 -Q .95 -o "$OUTPUT_VCF" -O z
