# Changelog (publication repository cleanup)

Original scripts are preserved verbatim under `_originals/` before any cleaned copies are added.

| Date | File | Change |
|------|------|--------|
| 2026-07-08 | Repository structure | Initial methods-aligned layout; scripts copied from Costus_popgen_clean, Population_genetics_vill_alle, alle_vill_QTL, seq_analysis, humm_backup |
| 2026-07-08 | `sample_manifest/build_combined_manifest.py` | New: merge wild + QTL NCBI metadata |
| 2026-07-08 | `03_QTL_mapping/rqtl/run_scantwo_per_trait.R` | New: unified Farm scantwo script |
| 2026-07-08 | `03_QTL_mapping/lepmap3/generate_rqtl_input.R` | Clean copy of `generate_rQTL_input_files.R` (env-var paths, fixed line-29 typo) |
| 2026-07-08 | `sample_manifest/combined_sample_manifest.tsv` | Manual revision + dedup to 524 unique `sample_id` |
| 2026-07-08 | Hummingbird scripts added | Dsuite, Lep-MAP3 OrderMarkers2, pixy; backed up to `_originals/` |
| 2026-07-08 | `run_Dinvestigate_fdm_array.slurm`, `run_Dtrios_capseq.slurm`, `run_pixy_array.slurm`, `run_OrderMarkers2_array.slurm` | Farm-ready wrappers with env-var paths |
| 2026-07-08 | `scripts/project_paths.R` | Shared path resolution for introgression plotting scripts |
| 2026-07-08 | `04_introgression_vs_QTL/*.R` | Use `project_paths()` instead of hardcoded user paths |
| 2026-07-08 | `03_QTL_mapping/lepmap3/generate_lepmap_input.R` | Clean Lep-MAP3 .call builder |
| 2026-07-08 | `sample_manifest/deduplicate_manifest.py` | Utility: one row per `sample_id` |
