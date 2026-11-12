#!/bin/bash
# RETIRED. 50 kb Weir–Cockerham pixy job.
# Do not submit this. The Hudson rerun is run_pixy_array.slurm
# (WINDOW_SIZE=10000 and WINDOW_SIZE=500000).
# Full script, as it was: _originals/wc_fst_original_submission/02_population_genetics_wild/pixy/run_pixy_avonly.sh

#SBATCH --job-name=pixy_AVonly_retired
#SBATCH -A jrigrp
#SBATCH -p low
#SBATCH --time=00:05:00

echo "run_pixy_avonly.sh is retired (50 kb Weir-Cockerham)." >&2
echo "Submit run_pixy_array.slurm with WINDOW_SIZE=10000 or WINDOW_SIZE=500000." >&2
exit 1
