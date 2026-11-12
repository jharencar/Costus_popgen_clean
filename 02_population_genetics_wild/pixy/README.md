# PIXY genome-wide differentiation

Hudson Fst (Bhatia et al. 2013), pi, and Dxy. Each window's Fst is pixy's ratio of summed Hudson components (`avg_hudson_fst`). `--fst_components` also writes `hudson_fst_num` and `hudson_fst_den` so a multi-window summary can be `sum(num) / sum(den)` rather than a mean of window Fst values.

| Script | Purpose |
|--------|---------|
| `run_pixy_array.slurm` | Farm job. Submit twice: `WINDOW_SIZE=10000` and `WINDOW_SIZE=500000` |
| `aggregate_pixy_genomewide.py` | Genome-wide pi, Dxy, and Hudson Fst from the component columns |
| `plot_pixy.R` | Per-chromosome plots of the 10 kb Hudson output |

The 500 kb along-genome line in `04_introgression_vs_QTL/genome-wide_Fst.Dxy_w_QTL_plots.R` is the 500 kb pixy file, not a mean of the 10 kb rows.

Do not submit `run_pixy_avonly.sh`. That 50 kb Weir–Cockerham job is retired. The original-submission Weir–Cockerham scripts, including the full `run_pixy_avonly.sh`, are in [`_originals/wc_fst_original_submission/`](../../_originals/wc_fst_original_submission/). Hummingbird jobs `run_pixy.sh` and `512x64_run_pixy.sh` are historical and are not the Farm path.

```bash
sbatch --export=ALL,WINDOW_SIZE=10000  run_pixy_array.slurm
sbatch --export=ALL,WINDOW_SIZE=500000 run_pixy_array.slurm
# after the 10 kb array finishes:
python3 aggregate_pixy_genomewide.py /quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_10kb
```

Defaults: all-sites per-chromosome VCFs and `avonly_populations.txt` on Farm, conda env `pixy` (needs pixy >= 2.0). Output:

- `/quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_10kb`
- `/quobyte/jrigrp/juliah/Costus_data_tmp/pixy_hudson_500kb`

Supplementary parent subset (same 9 alle + 12 vill as Dsuite noNearHybrid 1 km):
`noNearHybrid_1km/` → Farm `~/Costus_scripts/pixy_noNearHybrid_1km/`, outputs
`pixy_hudson_{10kb,500kb}_noNearHybrid_1km/`.
