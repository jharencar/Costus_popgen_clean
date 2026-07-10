# PIXY genome-wide differentiation

| Script | Purpose |
|--------|---------|
| `run_pixy_array.slurm` | Farm-ready: pi, Fst, Dxy per chromosome (10 kb windows) |
| `run_pixy_avonly.sh` | Farm aggregate workflow (`humm_backup` copy) |
| `aggregate_pixy_genomewide.py` | Combine per-chr pixy outputs |
| `plot_pixy.R` | Per-chromosome plotting |
| `VCF_sample_subset.sh` | Subset VCF before pixy |
| `run_pixy.sh`, `512x64_run_pixy.sh` | Original Hummingbird jobs (`_originals/`) |
