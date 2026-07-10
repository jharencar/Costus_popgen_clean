# Costus allenii × C. villosissimus — population genomics and QTL code

Code accompanying the manuscript on ecogeographic isolation, introgression, and the genetic basis of habitat-associated traits in *Costus allenii* and *C. villosissimus*.

Methods prose: [`popgen_manu_methods.txt`](popgen_manu_methods.txt). Internal dev notes: [`NOTES.md`](NOTES.md).

## Repository map

| Folder | Methods section |
|--------|-----------------|
| [`01_sequencing_and_variant_calling/`](01_sequencing_and_variant_calling/) | DNA extraction, sequencing, trimming, alignment, variant filters |
| [`02_population_genetics_wild/`](02_population_genetics_wild/) | PCA, aHMM, ABBA-BABA, Dsuite, PIXY |
| [`03_QTL_mapping/`](03_QTL_mapping/) | F3 aHMM, Lep-MAP3, R/qtl, chemistry QTL |
| [`04_introgression_vs_QTL/`](04_introgression_vs_QTL/) | Fst vs QTL window comparisons |
| [`sample_manifest/`](sample_manifest/) | Combined wild + mapping population sample table |

Original unmodified scripts: [`_originals/`](_originals/). Cleaned copies use env-var paths (see [`CHANGELOG.md`](CHANGELOG.md)). Shared R paths: [`scripts/project_paths.R`](scripts/project_paths.R) (`COSTUS_ARCHIVE`, `PIXY_FOLDER`, `DSUITE_DIR`, etc.).

**HPC shell scripts** (`*_hb.sh`, `LASIout_*.sh`, `01_sequencing_and_variant_calling/preprocessing/`): document Hummingbird cluster paths for historical reruns. Farm/local wrappers (`run_*_array.slurm`) use environment variables instead. Set `JULIA_POPGEN` when adapting cluster scripts to a new filesystem.

## Working on this machine

- **Publication cleanup** uses only the local git repo (`Costus_popgen_clean`). Avoid browsing `/Volumes/dissert` in Finder if the volume does not mount in the sidebar — it can stress Finder. Use Terminal for external-drive paths when needed.
- **Sample manifest:** after manual edits, deduplicate with `python3 sample_manifest/deduplicate_manifest.py` (one row per `sample_id`).

## Software

| Tool | Version (methods) |
|------|-------------------|
| Trimmomatic | 0.39 |
| BBTools dedupe | 39.01 |
| BWA mem | 0.7.15 |
| genmap | k=15, merge 100 bp |
| ANGSD / PCAngsd | — |
| Ancestry HMM | uniform error 0.02 |
| Dsuite | Dtrios, Dinvestigate |
| PIXY | 10 kb windows |
| Lep-MAP3 | LOD 11 |
| R/qtl | Haley–Knott, 1000 permutations |

## Quick start

1. **Sample manifest:** `python3 sample_manifest/build_combined_manifest.py`
2. **QTL mapping:** open [`03_QTL_mapping/rqtl/Costus_QTL.qmd`](03_QTL_mapping/rqtl/Costus_QTL.qmd) or run chunks sequentially. Inputs in [`03_QTL_mapping/data/`](03_QTL_mapping/data/).
3. **Farm scantwo perms:** `bash 03_QTL_mapping/rqtl/setup_rqtl_conda_farm.sh` then `sbatch 03_QTL_mapping/rqtl/slurm/run_scantwo_all_traits.slurm`
4. **Chemistry QTL:** `sbatch 03_QTL_mapping/rqtl/chemistry/run_WatersMass_chemistry_QTL.slurm` (swap in final chemistry file when ready)

## External data (not in repository)

Large analysis files are archived on the **dissert** drive at `Costus_archive/` (not committed to git). Set `COSTUS_ARCHIVE` if the mount path differs (e.g. `/run/media/julia/dissert/Costus_archive` on Linux).

| Used by (methods code) | Archive path (under `$COSTUS_ARCHIVE`) |
|------------------------|----------------------------------------|
| `04_introgression_vs_QTL/*.R`, `02_population_genetics_wild/pixy/plot_pixy.R` | `06_pixy/output_251124_LASIref/` (or `06_pixy/from_hummingbird/pixy/output_latest/`) |
| `04_introgression_vs_QTL/*.R` | `07_dsuite/Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120/` |
| `02_population_genetics_wild/aHMM/plot_AHV_LASIref2_all_individuals_genomewide_stripes.R` | `03_aHMM_posteriors/wild_AHV_LASIref2/AHV_LASIref2_260210/` |
| `03_QTL_mapping/lepmap3/generate_*_input.R`, `ordermarker_plots.R` | `03_aHMM_posteriors/QTL_grouped_run/.../no_thin_cated_posteriors_.../` (posteriors); `04_lepmap3/lepmap3/` (order files) |
| `03_QTL_mapping/rqtl/Costus_QTL.qmd`, `run_scantwo_per_trait.R`, chemistry QTL | `05_rqtl/` (`*_genoprob.rds`, `scantwo.perms/`) |
| Dsuite / pixy HPC reruns (`run_*_array.slurm`) | `02_VCFs/BAVLasi_LASIref_filtered_scaf/`, `02_VCFs/renamed_AVonly_pixy/` |
| ABBA-BABA (`02_population_genetics_wild/abba_baba/`) | `08_ANGSD/abba_baba/from_hummingbird/` |

R scripts resolve paths via [`scripts/project_paths.R`](scripts/project_paths.R) (`PIXY_FOLDER`, `DSUITE_DIR`, `GENOPROB_RDS`, `POSTERIOR_DIR`, etc.). Small CSV inputs ship in [`03_QTL_mapping/data/`](03_QTL_mapping/data/).

**Archive layout and rsync history:** [`sample_manifest/HUMMINGBIRD_DATA_ARCHIVE.md`](sample_manifest/HUMMINGBIRD_DATA_ARCHIVE.md).

**QTL notebook:** render `Costus_QTL.qmd` with working directory `$COSTUS_ARCHIVE/05_rqtl` (or set `GENOPROB_RDS` / symlink `.rds` files into `03_QTL_mapping/rqtl/`).

## HPC scripts

Hummingbird originals are in `_originals/` and beside Farm-ready `run_*.slurm` wrappers in each subfolder. See README in `02_population_genetics_wild/dsuite/`, `03_QTL_mapping/lepmap3/`, and `03_QTL_mapping/rqtl/`.

## Figure reproduction

| Script | Figure type |
|--------|-------------|
| `02_population_genetics_wild/aHMM/plot_AHV_LASIref2_all_individuals_genomewide_stripes.R` | Wild hybrid ancestry stripes |
| `04_introgression_vs_QTL/combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R` | Genome-wide introgression + differentiation |
| `04_introgression_vs_QTL/genome-wide_Fst.Dxy_w_QTL_plots.R` | Fst vs QTL CIs |
| `03_QTL_mapping/rqtl/Costus_QTL.qmd` | QTL scanone / model figures |

## Citation

*(Add citation when manuscript is published.)*
