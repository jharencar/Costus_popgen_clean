# Costus allenii × C. villosissimus — population genomics and QTL code

Code accompanying the manuscript on ecogeographic isolation, introgression, and the genetic basis of habitat-associated traits in *Costus allenii* and *C. villosissimus*.

Methods prose: [`popgen_manu_methods.txt`](popgen_manu_methods.txt). Developer notes: [`NOTES.md`](NOTES.md).

## Repository map

| Folder | Methods section |
|--------|-----------------|
| [`01_sequencing_and_variant_calling/`](01_sequencing_and_variant_calling/) | Read processing, alignment, variant filters |
| [`02_population_genetics_wild/`](02_population_genetics_wild/) | PCA, aHMM, ABBA-BABA, Dsuite, PIXY |
| [`03_QTL_mapping/`](03_QTL_mapping/) | F3 aHMM, Lep-MAP3, R/qtl |
| [`04_introgression_vs_QTL/`](04_introgression_vs_QTL/) | Fst / introgression vs QTL window comparisons |
| [`sample_manifest/`](sample_manifest/) | Combined wild + mapping population sample table |
| [`scripts/project_paths.R`](scripts/project_paths.R) | Shared path resolution (`COSTUS_ARCHIVE`, `PIXY_FOLDER`, `DSUITE_DIR`, `GENOPROB_RDS`, …) |

Verbatim pre-cleanup copies: [`_originals/`](_originals/). Farm-ready `run_*.slurm` wrappers use environment variables; Hummingbird shell scripts keep historical `/hb/` paths.

## Sequences and reference (NCBI)

| Resource | Accession / location |
|----------|----------------------|
| Reference genome (*C. lasius*) | [GCA_027563935.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027563935.1/) |
| Low-coverage WGS FASTQs (wild + mapping pop) | NCBI SRA BioProject *(accession TBD — update when deposited)* |

Do **not** re-upload FASTQs to Dryad. Reconstruct BAMs / VCFs / aHMM posteriors from NCBI FASTQs using the scripts in `01_`–`03_` if a full pipeline rerun is required.

## What is regenerable vs provided elsewhere

| Regenerable with this repo + NCBI FASTQs | Provided for convenience (Dryad, not GitHub) |
|------------------------------------------|-----------------------------------------------|
| Trimmed reads, BAMs, filtered VCFs | Canonical QTL `*_genoprob.rds`, scantwo permutation `.rds` |
| aHMM `*.posterior`, Lep-MAP3 intermediates | Small phenotype / genotype CSVs (also in `03_QTL_mapping/data/`) |
| Dsuite / PIXY genome-wide runs from VCFs | PIXY + Dsuite **plot-input tables** used by figure scripts |

GitHub holds **code + small CSVs only**. Large objects ship in the Dryad package (see below) or a local `COSTUS_ARCHIVE` mirror.

## Dryad package

Dryad DOI: *(placeholder — add after deposit)*.

Minimal Dryad zip contents (see `DRYAD_README.md` inside the package):

- This publication code tree (`01_`–`04_`, `scripts/`, `sample_manifest/`, README, methods text)
- QTL objects: `most_dups_genoprob.rds` (+ related crosses / scantwo perms)
- Small inputs: phenotype CSV, QTL CI CSV, genotype CSV, sample manifest
- Figure-input tables: PIXY `output_251124_LASIref/`, Dsuite `Fstats_w50s25_.../`

Excluded from Dryad: FASTQs, BAMs, VCFs, aHMM posteriors, unpublished analyses.

## Quick start

1. Clone this repository; install R packages used by the figure / QTL scripts (`tidyverse`, `qtl`, `patchwork`, `ggrastr`, …).
2. Download the Dryad package (or set `COSTUS_ARCHIVE` to a local archive with the same layout).
3. `source("scripts/project_paths.R")` — or export env vars from [`sample_manifest/archive_paths.env.example`](sample_manifest/archive_paths.env.example).
4. QTL notebook: [`03_QTL_mapping/rqtl/Costus_QTL.qmd`](03_QTL_mapping/rqtl/Costus_QTL.qmd) with working directory pointing at Dryad/`05_rqtl` (or set `GENOPROB_RDS`).
5. Figures: scripts listed below under Figure reproduction.

## Software (methods)

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

## Figure reproduction

| Script | Figure type |
|--------|-------------|
| `02_population_genetics_wild/aHMM/plot_AHV_LASIref2_all_individuals_genomewide_stripes.R` | Wild hybrid ancestry stripes |
| `04_introgression_vs_QTL/combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R` | Genome-wide introgression + differentiation |
| `04_introgression_vs_QTL/genome-wide_Fst.Dxy_w_QTL_plots.R` | Fst vs QTL CIs |
| `03_QTL_mapping/rqtl/Costus_QTL.qmd` | QTL scanone / model figures |

## Citation

*(Add citation when manuscript is published.)*
