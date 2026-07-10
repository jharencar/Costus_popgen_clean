# scantwo permutations (R/qtl)

**Canonical script for new Farm runs:** [`run_scantwo_per_trait.R`](run_scantwo_per_trait.R)  
**Farm array job:** [`slurm/run_scantwo_all_traits.slurm`](slurm/run_scantwo_all_traits.slurm)

Set `TRAIT` to `SD`, `GR`, `tough`, `Nreabs`, or `Chlreabs`. Outputs go to `scantwo.perms/*.scantwo.hk.perm.rds`.

## Archives (Hummingbird originals)

Per-trait scripts copied from Hummingbird are kept for reference:

- `SD.scantwo.R`, `GR.scantwo.R`, etc. — R scripts used on Hummingbird
- `*_scantwo_hb.sh` — Hummingbird slurm wrappers
- `slurm/_originals_humm/` — duplicate backup of Hummingbird versions

Use `run_scantwo_per_trait.R` for reproducible Farm runs unless you need to match an exact historical cross object (e.g. `most_dups` vs `least_dups`).
