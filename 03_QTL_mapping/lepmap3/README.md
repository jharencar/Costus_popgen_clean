# Lep-MAP3 linkage mapping

| Script | Purpose |
|--------|---------|
| `generate_lepmap_input.R` | aHMM posteriors → `.call` file (env-var paths) |
| `generate_rqtl_input.R` | aHMM posteriors → rQTL CSV |
| `pedigree.txt` | F2 mapping-population pedigree (ParentCall2 input) |
| `run_OrderMarkers2_array.slurm` | OrderMarkers2 per chr (Farm-ready) |
| `lepmap_OrderMarkers2_manMap.sh` | Original Hummingbird job (`_originals/`) |
| `ordermarker_plots.R` | Genetic map QC plots |
| `manual_map3.txt` | Manual map for chr4 join / centromere fixes |

**Methods parameters:** LOD 11, theta 0.001, distortionLod 11; hard calls >90%; ParentCall2 with pedigree.

**Pedigree:** `pedigree.txt` is the tab-delimited Lep-MAP3 pedigree used for the manuscript run (6 header rows — family, individual, father, mother, sex, phenotype — over 419 individuals; the "all F2s, 0.05× coverage cutoff" version). Prepend it to the genotype `.call` file from `generate_lepmap_input.R` to build the ParentCall2 input. The full pedigree-*building* script (`generate_lepmap_input_files_forJulia.R`) is preserved in the repository under `_originals/` and in git history.
