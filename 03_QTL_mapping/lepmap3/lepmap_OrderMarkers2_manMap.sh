#!/bin/bash
#SBATCH --partition 512x64-ib 
#SBATCH --qos=512x64
#SBATCH --account=512x64
#SBATCH --job-name=lepmap_OrderMarkers2
#SBATCH --output=lepmap_OrderMarkers2_%A_%a.out
#SBATCH -e lepmap_OrderMarkers2_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --mem=80
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=1-9

module load java
cd /hb/groups/kay_lab/Julia_popgen/lepmap

java -cp /hb/groups/kay_lab/Software/Lep-MAP3/bin OrderMarkers2 map=manual_map3.txt data=refV2_thinned_50k_data.call scale=2M/N sexAveraged=1 outputPhasedData=0 chromosome=${SLURM_ARRAY_TASK_ID} > order${SLURM_ARRAY_TASK_ID}_LodLim4.txt
