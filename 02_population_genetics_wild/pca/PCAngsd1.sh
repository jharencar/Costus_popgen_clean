#!/bin/bash
#SBATCH --partition 128x24 # this means each node (there are 19 nodes in this partition) has 128 GB RAM and 24 cores/cpus 
#SBATCH --job-name=PCAngsd1 # Job name
#SBATCH -e PCAngsd1.err
#SBATCH --output=PCAngsd1.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR.EMAIL
#SBATCH --time=24:00:00  # ran in roughly 9 hours for all 96 merged bams of plate 1
#SBATCH --mem=80G # this is pretty arbitrary/set when I really had no idea what I was doing, not sure if it is worth optimizing now...
#SBATCH --nodes=2 # this is pretty arbitrary/set when I really had no idea what I was doing, not sure if it is worth optimizing now...

## load ANGSD
module load angsd samtools

## change into the directory with .bams
cd /hb/YOUR/PATH/HERE

## set directories
IN_DIR=/hb/YOUR/PATH/HERE/ANGSD/angsd_in
OUT_DIR=/hb/YOUR/PATH/HERE/ANGSD/PCAngsd_out

# call pcangsd
# options (and defaults): http://www.popgen.dk/software/index.php/PCAngsd
# default params used in these runs 
#(-iter 0, means without estimating individual allele freqs. this seems preferable given low individual coverage, 
# but the plots more clearly demonstrate the pattern when we do the estimate. See PCAngsd_angsd1.Rmd for comparison of plots)

#run on LGVHA
echo "python /hb/groups/kay_lab/Software/pcangsd/pcangsd.py -beagle $IN_DIR/YOUR_FILE_NAME.genolike.beagle.gz -o $OUT_DIR/YOUR_FILE_NAME.PCAngsd"
python /hb/groups/kay_lab/Software/pcangsd/pcangsd.py -beagle $IN_DIR/YOUR_FILE_NAME.genolike.beagle.gz -o $OUT_DIR/YOUR_FILE_NAME.PCAngsd