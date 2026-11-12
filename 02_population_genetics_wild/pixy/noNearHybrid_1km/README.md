# Supplementary pixy: no individuals within 1 km of hybrids

Same parent subset as `dsuite/noNearHybrid_1km/` (9 *allenii* + 12 *villosissimus*
present in the AVonly VCFs). Original Hudson pixy runs are unchanged.

## Farm

Scripts sync to `~/Costus_scripts/pixy_noNearHybrid_1km/`. Reuses
`run_pixy_array.slurm` from `~/Costus_scripts/pixy_whole_genome/` with
environment overrides.

```bash
cd ~/Costus_scripts/pixy_noNearHybrid_1km

POPFILE=$PWD/avonly_populations_noNearHybrid_1km.txt
SLURM=~/Costus_scripts/pixy_whole_genome/run_pixy_array.slurm

sbatch --export=ALL,WINDOW_SIZE=10000,POPFILE=$POPFILE,OUTDIR=/quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_10kb_noNearHybrid_1km $SLURM
sbatch --export=ALL,WINDOW_SIZE=500000,POPFILE=$POPFILE,OUTDIR=/quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_500kb_noNearHybrid_1km $SLURM
```

## Outputs

- `/quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_10kb_noNearHybrid_1km/`
- `/quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_500kb_noNearHybrid_1km/`

Plotting uses these with the subset Dsuite folder via `PLOT_OUT_ROOT` /
`PIXY_FOLDER` / `DSUITE_DIR` (see `04_introgression_vs_QTL/noNearHybrid_1km/`).
