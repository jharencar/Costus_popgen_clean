# Original-submission Weir–Cockerham pixy workflow

Copies of the scripts that produced and plotted the November 2025 10 kb
Weir–Cockerham `avg_wc_fst` tables (`output_251124_LASIref`) for the original
submission, taken before the Hudson rerun. The Hudson scripts in the live tree
are the ones that will replace this version. The July 2026 verbatim tree
elsewhere under `_originals/` is unchanged.

Live scripts run Hudson Fst (`--fst_type hudson --fst_components`) on Farm
at 10 kb and 500 kb. The 500 kb along-genome line is pixy's own window, not a
mean of 10 kb rows.

| Archived file | Was |
|---------------|-----|
| `02_population_genetics_wild/pixy/run_pixy_array.slurm` | Farm 10 kb job, default Weir–Cockerham, no components |
| `02_population_genetics_wild/pixy/run_pixy_avonly.sh` | Farm 50 kb Weir–Cockerham job with `--fst_components` |
| `02_population_genetics_wild/pixy/aggregate_pixy_genomewide.py` | Genome-wide Fst from `wc_fst_a/b/c` only |
| `02_population_genetics_wild/pixy/plot_pixy.R` | Per-chromosome plots of `avg_wc_fst` |
| `04_introgression_vs_QTL/genome-wide_Fst.Dxy_w_QTL_plots.R` | 500 kb line = mean of 10 kb windows |
| `04_introgression_vs_QTL/combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R` | Four-panel figure using `avg_wc_fst` |
