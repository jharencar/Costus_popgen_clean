# Lep-MAP3 linkage mapping

| Script | Purpose |
|--------|---------|
| `generate_lepmap_input.R` | aHMM posteriors → `.call` file (env-var paths) |
| `generate_lepmap_input_files_forJulia.R` | Original full script incl. pedigree (`_originals/`) |
| `generate_rqtl_input.R` | aHMM posteriors → rQTL CSV |
| `generate_rQTL_input_files.R` | Original rQTL input script (`_originals/`) |
| `run_OrderMarkers2_array.slurm` | OrderMarkers2 per chr (Farm-ready) |
| `lepmap_OrderMarkers2_manMap.sh` | Original Hummingbird job (`_originals/`) |
| `ordermarker_plots.R` | Genetic map QC plots |
| `manual_map3.txt` | Manual map for chr4 join / centromere fixes |

**Methods parameters:** LOD 11, theta 0.001, distortionLod 11; hard calls >90%; ParentCall2 with pedigree.

**Note:** Pedigree rows are hardcoded in `generate_lepmap_input_files_forJulia.R`. For a full rerun, use that original or export pedigree to a separate TSV.
