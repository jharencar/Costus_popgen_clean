#!/bin/bash
#SBATCH --partition 128x24
#SBATCH --job-name=pixy
#SBATCH --output=pixy.%A_%a.out
#SBATCH -e pixy.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --array=1-9

module load miniconda3
conda activate JGH_conda_env|| { echo "Failed to activate Conda environment"; exit 1; }

cd /hb/groups/kay_lab/Julia_popgen/pixy

#run pixy
pixy --stats pi fst dxy \
--vcf /hb/groups/kay_lab/Julia_popgen/VCFs/renamed_AVonly_all_sites_mappability_maxDepth_filtered_chr${SLURM_ARRAY_TASK_ID}.vcf.gz \
--populations AV_sampleIDs_popfile.txt \
--window_size 10000 \
--n_cores 11 \
--output_folder output_map_filt_${SLURM_ARRAY_TASK_ID} \
--output_prefix chr${SLURM_ARRAY_TASK_ID}
