# Costus allenii × C. villosissimus — population genomics and QTL code

Code accompanying the manuscript on ecogeographic isolation, introgression, and the genetic basis of habitat-associated traits in *Costus allenii* and *C. villosissimus*.

Methods: see the associated manuscript, "Divergent habitat selection across many loci maintains species boundaries during introgression." Developer notes: [`NOTES.md`](NOTES.md).

## Repository map

| Folder | Methods section |
|--------|-----------------|
| [`01_sequencing_and_variant_calling/`](01_sequencing_and_variant_calling/) | Read processing, alignment, variant filters |
| [`02_population_genetics_wild/`](02_population_genetics_wild/) | PCA, aHMM, ABBA-BABA, Dsuite, PIXY |
| [`03_QTL_mapping/`](03_QTL_mapping/) | F3 aHMM, Lep-MAP3, R/qtl |
| [`04_introgression_vs_QTL/`](04_introgression_vs_QTL/) | Fst / introgression vs QTL window comparisons |
| [`sample_manifest/`](sample_manifest/) | Pointer to Dryad sample table (TSV not in this code archive) |
| [`scripts/project_paths.R`](scripts/project_paths.R) | Shared path resolution (`COSTUS_ARCHIVE`, `PIXY_FOLDER`, `DSUITE_DIR`, `GENOPROB_RDS`, …) |

Verbatim pre-cleanup copies: [`_originals/`](_originals/). Farm-ready `run_*.slurm` wrappers use environment variables; Hummingbird shell scripts keep historical `/hb/` paths.

## Sequences and reference (NCBI)

| Resource | Accession / location |
|----------|----------------------|
| Reference genome (*C. lasius*) | [GCA_027563935.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027563935.1/) |
| Low-coverage WGS FASTQs (wild + mapping pop) | NCBI SRA BioProject [PRJNA1479894](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1479894) |

## How Dryad (data) and Zenodo (code) fit together

Code and data are archived separately so licenses stay unambiguous:

- **Zenodo** (this repository snapshot, [MIT License](LICENSE)): analysis scripts and docs only. Zenodo DOI: *(placeholder — add after deposit)*.
- **Dryad** (CC0): all tabular / figure-input / permutation **data**. Dryad DOI: *(placeholder — add after deposit)*.

**No data files are duplicated in Zenodo.** Download both deposits, unpack them side by side, and point the code at the Dryad paths via `COSTUS_ARCHIVE` (see Dryad README). Full reuse instructions live in the **single Dryad README** (prefer Dryad's in-browser README field).

| Regenerable with this code + NCBI FASTQs | Provided on Dryad only |
|------------------------------------------|------------------------|
| Trimmed reads, BAMs, filtered VCFs | Sample manifest TSV |
| aHMM `*.posterior`, Lep-MAP3 intermediates | Phenotype / genotype / QTL CI CSVs |
| Dsuite / PIXY genome-wide runs from VCFs | PIXY + Dsuite plot-input tables |
| QTL `*_genoprob.rds` / `est.rf.rds` (from Dryad CSVs) | scantwo permutation `.rds` |

## License

Code in this repository is released under the [MIT License](LICENSE). Data deposited
on Dryad is released under CC0 per Dryad policy and is not redistributed with this code.

## Citation

*(Add citation when manuscript is published.)*
