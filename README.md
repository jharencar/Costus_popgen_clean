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

## What is regenerable vs provided elsewhere

| Regenerable with this repo + NCBI FASTQs | Provided for convenience (Dryad, not GitHub) |
|------------------------------------------|-----------------------------------------------|
| Trimmed reads, BAMs, filtered VCFs | Canonical QTL `*_genoprob.rds`, scantwo permutation `.rds` |
| aHMM `*.posterior`, Lep-MAP3 intermediates | Small phenotype / genotype CSVs (also in `03_QTL_mapping/data/`) |
| Dsuite / PIXY genome-wide runs from VCFs | PIXY + Dsuite **plot-input tables** used by figure scripts |

## Dryad package

Dryad DOI: *(placeholder — add after deposit)*.

Minimal Dryad zip contents (see `DRYAD_README.md` inside the package):

- This publication code tree (`01_`–`04_`, `scripts/`, `sample_manifest/`, README, methods text)
- QTL objects: `most_dups_genoprob.rds` (+ related crosses / scantwo perms)
- Small inputs: phenotype CSV, QTL CI CSV, genotype CSV, sample manifest
- Figure-input tables: PIXY `output_251124_LASIref/`, Dsuite `Fstats_w50s25_.../`

## Citation

*(Add citation when manuscript is published.)*
