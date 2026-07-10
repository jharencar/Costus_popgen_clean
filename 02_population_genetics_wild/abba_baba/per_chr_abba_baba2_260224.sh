#!/bin/bash
#SBATCH --partition 512x64
#SBATCH --qos 512x64
#SBATCH --account 512x64
#SBATCH --job-name=AHV_ANGSD_abbababa2
#SBATCH -e AHV_ANGSD_abbababa2.%A_%a.err
#SBATCH --output=AHV_ANGSD_abbababa2.%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-9

# load modules
module load angsd samtools

# change directories
cd /hb/groups/kay_lab/Julia_popgen/ANGSD/abba_baba

# set variables
BAMLIST="v2_ref_BAVLasi_clean_0.15co_240405.bam.filelist"
SITES="filtered_BAVLasi_LASI2ref_chr${SLURM_ARRAY_TASK_ID}.sites"
CHR=$(cut -f1 "$SITES" | head -n 1)
OUT="lasiOUT_${CHR}_260224"

# index sites file
angsd sites index "$SITES"

# run doAbbababa2 with site filtering
angsd \
    -doAbbababa2 1 \
    -bam "$BAMLIST" \
    -sizeFile sizeFile.size \
    -doCounts 1 \
    -out "$OUT" \
    -anc /hb/groups/kay_lab/genomes/C.lasius_genome_NCBI_v2.fasta \
    -ref /hb/groups/kay_lab/genomes/C.lasius_genome_NCBI_v2.fasta \
    -minQ 30 \
    -minMapQ 30 \
    -only_proper_pairs 1 \
    -baq 1 \
    -blockSize 5000000 \
    -nThreads 8 \
    -sites "$SITES" \
    -rf <(echo "$CHR")\
    -useLast 1
