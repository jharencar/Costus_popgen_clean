# Lep-MAP3 linkage mapping

| Script | Purpose |
|--------|---------|
| `generate_lepmap_input.R` | aHMM posteriors → `.call` file (env-var paths) |
| `generate_rqtl_input.R` | aHMM posteriors → rQTL CSV |
| `run_OrderMarkers2_array.slurm` | OrderMarkers2 per chr (Farm-ready) |
| `lepmap_OrderMarkers2_manMap.sh` | Original Hummingbird job (`_originals/`) |
| `ordermarker_plots.R` | Genetic map QC plots |
| `manual_map3.txt` | Manual map for chr4 join / centromere fixes |

**Methods parameters:** LOD 11, theta 0.001, distortionLod 11; hard calls >90%; ParentCall2 with pedigree.

**Note:** `generate_lepmap_input.R` omits the hardcoded F2 pedigree; the full pedigree-building version is preserved in the repository under `_originals/03_QTL_mapping/lepmap3/generate_lepmap_input_files_forJulia.R` (and in git history). For a full rerun, restore that script or export the pedigree to a separate TSV.
