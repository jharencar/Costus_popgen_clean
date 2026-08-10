# QTL tabular inputs (Dryad only)

The phenotype, genotype, and QTL CI CSV files formerly kept here are archived **only in
the Dryad data deposit** (CC0), under `data_qtl/small_inputs/`:

| File | Role |
|------|------|
| `F3phe_onehotcovar_newHI_20250925.csv` | F3 phenotypes and covariates |
| `f3gen_250924_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001.csv` | Genotype table for R/qtl |
| `all_QTL_CIs_bp_cm.251130.csv` | Published QTL peaks and 95% Bayes CIs |

They are intentionally omitted from this code tree and from the Zenodo software archive to
avoid licensing ambiguity with Dryad (CC0 data vs MIT code). Download the Dryad package and
point scripts at those paths (see Dryad README / `scripts/project_paths.R`).
