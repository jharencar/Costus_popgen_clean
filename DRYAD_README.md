# Costus allenii × C. villosissimus — population genomics and QTL

[Access this dataset on Dryad](DOI_PLACEHOLDER)

This repository holds the data and code used for analyses in the manuscript, "Divergent habitat selection across many loci maintains species boundaries during introgression." The analyses investigate ecogeographic isolation, introgression, and QTL for habitat-associated traits in Neotropical *Costus allenii* and *C. villosissimus*. It contains (1) analysis and figure code, (2) expensive-to-regenerate QTL R objects and permutation results, (3) small phenotype/genotype/CI tables, and (4) PIXY and Dsuite window tables used by manuscript figure scripts.

**Not included here** (available elsewhere or regenerable): raw FASTQs, BAMs, VCFs, Ancestry_HMM posterior files, and other large intermediates. See **Sharing/Access information**.

Package build date: `2026-07-22`. Zip checksum (SHA-256): `5ce59113f8688082ad141585bad4773bb9302d964d1f40cdd4ebb861df4976f2` (`manifest/`).

---



## Description of the data and file structure



### Study system and experimental summary

- **Species:** *Costus allenii*, *C. villosissimus*, wild hybrids; *C. bracteatus* (ingroup for ABBA-BABA / Dsuite); *C. lasius* (outgroup and reference genome).
- **Wild WGS:** Central Panama collections (plus *C. bracteatus* from Costa Rica; *C. lasius* from Panama and Peru). 
- **Mapping population:** F3 progeny of an *allenii* × *villosissimus* cross (plus F1 parent and F2s used in linkage mapping). QTL mapping retains **385** genotyped individuals after QC.
- **Traits scored in the F3 population (habitat-associated):** seed dormancy (days to radicle), growth rate (cm/day), leaf toughness and thickness, percent nitrogen and chlorophyll reabsorption during drought-induced senescence, and leaf nitrogen content.
- **Analyses covered by this package + NCBI sequences:** low-coverage WGS processing; wild population genetics (PCA, Ancestry_HMM, ABBA-BABA, Dsuite, PIXY); F3 Ancestry_HMM → Lep-MAP3 → R/qtl; comparison of differentiation / introgression statistics with QTL intervals.

Full methods prose is in the manuscript, "Divergent habitat selection across many loci maintains species boundaries during introgression."

### Package layout

```
Costus_Dryad_package/
  DRYAD_README.md                 # this file
  code/                           # publication code tree (methods scripts + docs)
  data_qtl/
    05_rqtl/                      # R/qtl cross objects (genoprob, est.rf)
    scantwo.perms/                # scantwo permutation RDS by trait
    small_inputs/                 # phenotype / genotype / QTL CI CSVs
  data_figure_inputs/
    pixy_output_251124_LASIref/   # PIXY π, Dxy, Fst (10 kb windows)
    Dsuite_Fstats_w50s25_LASIref_251120/  # Dsuite local F-statistics (w50 s25)
  manifest/                       # build date + zip checksum
```

### How the pieces relate

1. Use `code/` with NCBI FASTQs (and the reference genome) to regenerate pipelines from scratch, **or**
2. Point `code/scripts/project_paths.R` at this repo (see **Code/Software**) to reload QTL objects and figure-input tables without re-running cluster jobs.
3. `data_qtl/small_inputs/` are the tabular phenotype/genotype/CI inputs used to build / annotate the R/qtl objects.
4. `code/sample_manifest/combined_sample_manifest.tsv` is the sample metadata table (wild + mapping population).
5. `data_figure_inputs/` are the windowed summary tables read by introgression-vs-QTL and genome-wide figure scripts in `code/04_introgression_vs_QTL/` (and related plotting scripts documented in `code/NOTES.md`).



### File descriptions

#### `code/`

Publication GitHub tree for this project (`https://github.com/jharencar/Costus_popgen_clean`), including:


| Path                                               | Contents                                                 |
| -------------------------------------------------- | -------------------------------------------------------- |
| `01_sequencing_and_variant_calling/`               | Read processing, alignment, variant filters              |
| `02_population_genetics_wild/`                     | PCA, aHMM, ABBA-BABA, Dsuite, PIXY workflows             |
| `03_QTL_mapping/`                                  | F3 aHMM, Lep-MAP3, R/qtl scripts                         |
| `04_introgression_vs_QTL/`                         | Fst / introgression vs QTL comparisons and figures       |
| `sample_manifest/`                                 | Final sample table only (`combined_sample_manifest.tsv`) |
| `scripts/project_paths.R`                          | Shared path resolution via environment variables         |
| `README.md`, `NOTES.md` | Documentation                                            |




#### `code/sample_manifest/combined_sample_manifest.tsv`

Final sample metadata for wild collections and the QTL mapping population (**524** rows; one row per unique `sample_id`). This is the only sample table in the package (build scripts and NCBI source tables are not included).


| Column                  | Definition                            |
| ----------------------- | ------------------------------------- |
| `dataset`               | `wild_collection` or `qtl_mapping`    |
| `sample_id`             | Lab / analysis sample ID              |
| `organism`              | Species or hybrid designation         |
| `collection_date`       | Collection date (when applicable)     |
| `geographic_location`   | Locality string                       |
| `field_collection_id`   | Field collection ID                   |
| `latitude`, `longitude` | Decimal degrees (WGS84); may be blank |
| `elevation_m`           | Elevation in meters; may be blank     |


**Missing data:** empty cells for optional locality fields when not recorded.

#### `data_qtl/05_rqtl/` (R binary objects)

All are R serialized objects (`.rds`) for use with **R/qtl**. Shared filename stem:

`aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_<dup_level>_<stage>.rds`


| File                        | Description                                                                                                                                                                                    |
| --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `…_most_dups_genoprob.rds`  | **Canonical** R/qtl `cross` after `calc.genoprob()` on the moderate duplicate-marker thinning (~14,889 markers; 385 individuals). Use for published `fitqtl` / `refineqtl` / multi-QTL models. |
| `…_least_dups_genoprob.rds` | More aggressively thinned genoprob cross (~4,013 markers); used for some scans/permutations, not the primary published multi-QTL fits.                                                         |
| `…_least_dups_est.rf.rds`   | Cross after `est.rf()` (recombination fractions; before or without the genoprob stage used above).                                                                                             |


**Allele coding (locked):** AA = *C. villosissimus*; BB = *C. allenii*; AB = heterozygote. Genome-wide `hybrid_index` is high when ancestry is more *villosissimus*.

#### `data_qtl/scantwo.perms/`

R/qtl `scantwo` permutation objects (Haley–Knott) used to obtain penalties for `stepwiseqtl` via `calc.penalties`.


| File                           | Trait abbreviation | Phenotype column in R/qtl                 |
| ------------------------------ | ------------------ | ----------------------------------------- |
| `SD.scantwo.hk.perm.rds`       | SD                 | `log_days_to_radicle` (log seed dormancy) |
| `GR.scantwo.hk.perm.rds`       | GR                 | `GR_cm.dy`                                |
| `tough.scantwo.hk.perm.rds`    | tough              | `ave_tough`                               |
| `Nreabs.scantwo.hk.perm.rds`   | Nreabs             | `N_reabsorp`                              |
| `Chlreabs.scantwo.hk.perm.rds` | Chlreabs           | `chloro_reabsorp`                         |




#### `data_qtl/small_inputs/`



##### `F3phe_onehotcovar_newHI_20250925.csv`

Phenotypes and covariates for the F3 mapping population (**385** rows after QTL QC; one row per individual).


| Column(s)               | Definition                                                               | Units / codes                       |
| ----------------------- | ------------------------------------------------------------------------ | ----------------------------------- |
| `ID`                    | Individual identifier (matches genotype file)                            | character                           |
| `family`                | F2 family label                                                          | character / factor                  |
| `hybrid_index`          | Genome-wide *C. villosissimus* ancestry proportion                       | 0–1                                 |
| `seed.cohort`, `cohort` | Seed-sowing cohort / measurement cohort                                  | categorical                         |
| `F21_*`                 | One-hot indicators for F2 families retained as covariates                | 0/1                                 |
| `s.coh*`, `coh*`        | One-hot indicators for seed / measurement cohorts retained as covariates | 0/1                                 |
| `days_to_radicle`       | Seed dormancy: days from sowing to radicle emergence                     | days (QTL models use log transform) |
| `GR_cm.dy`              | Growth rate                                                              | cm / day                            |
| `ave_thick`             | Mean leaf thickness (3 points)                                           | mm                                  |
| `ave_tough`             | Mean leaf toughness (penetrometer)                                       | g force                             |
| `N_reabsorp`            | Leaf nitrogen reabsorption during drought senescence                     | percent                             |
| `chloro_reabsorp`       | Chlorophyll reabsorption (SPAD-based) during drought senescence          | percent                             |
| `N_wght_p`              | Leaf nitrogen content                                                    | percent by weight                   |


**Missing data:** empty cells (read as `NA` in R). Trait-specific complete-case *n* varies (some individuals lack leaf or dormancy measures).

##### `f3gen_250924_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001.csv`

Wide genotype table for the mapping population (~391 rows including header alignment with IDs; ~36,783 marker columns before R/qtl QC thinning).


| Column(s)   | Definition                                             |
| ----------- | ------------------------------------------------------ |
| `ID`        | Individual identifier                                  |
| `ChromN_bp` | Marker columns named `Chrom<chromosome>_<position_bp>` |


**Genotype codes:** `A` = homozygous *villosissimus* (AA); `B` = homozygous *allenii* (BB); `H` = heterozygote; `-` = missing. (R/qtl may also represent these as AA/AB/BB depending on import.)

##### `all_QTL_CIs_bp_cm.251130.csv`

Published QTL peaks and 95% Bayes credible intervals (**9** rows).


| Column                                    | Definition                                            | Units                       |
| ----------------------------------------- | ----------------------------------------------------- | --------------------------- |
| `trait`                                   | Trait code: `SD`, `GR`, `tough`, `Nreabs`, `Chlreabs` | —                           |
| `chr`                                     | Chromosome / linkage group                            | 1–9                         |
| `peak_bp`, `lower95CI_bp`, `upper95CI_bp` | Peak and CI in physical coordinates                   | bp on *C. lasius* reference |
| `peak_cm`, `lower95CI_cm`, `upper95CI_cm` | Peak and CI on the linkage map                        | cM                          |
| `CI_size_bp`, `CI_size_cm`                | Interval widths                                       | bp, cM                      |




#### `data_figure_inputs/pixy_output_251124_LASIref/`

PIXY output for *allenii* vs *villosissimus* in **10 kb** non-overlapping windows on the *C. lasius* reference. Files are named `chr{1–9}_{pi,dxy,fst}.txt` (tab-delimited).

`*_pi.txt` **columns:** `pop`, `chromosome`, `window_pos_1`, `window_pos_2`, `avg_pi`, `no_sites`, `count_diffs`, `count_comparisons`, `count_missing`

`*_dxy.txt` **columns:** `pop1`, `pop2`, `chromosome`, `window_pos_1`, `window_pos_2`, `avg_dxy`, `no_sites`, `count_diffs`, `count_comparisons`, `count_missing`

`*_fst.txt` **columns:** `pop1`, `pop2`, `chromosome`, `window_pos_1`, `window_pos_2`, `avg_wc_fst`, `no_snps`

Population codes: `alle` = *C. allenii*, `vill` = *C. villosissimus*. **Missing windows:** `NA` (or empty numeric fields) when `no_sites` / `no_snps` is 0.

#### `data_figure_inputs/Dsuite_Fstats_w50s25_LASIref_251120/`

Dsuite `Dinvestigate` local F-statistics for trio *bracteatus*–*allenii*–*villosissimus* with *C. lasius* outgroup. Windows: **50** variants, step **25** variants. One file per chromosome:

`brac_alle_vill_localFstats_LASIout_chr{1–9}_50_25.txt`


| Column                     | Definition                               |
| -------------------------- | ---------------------------------------- |
| `chr`                      | Chromosome name (`Chrom1` …)             |
| `windowStart`, `windowEnd` | Physical span of the variant window (bp) |
| `D`                        | Patterson’s *D*                          |
| `f_d`                      | *f*_d introgression statistic            |
| `f_dM`                     | *f*_dM                                   |
| `d_f`                      | *d*_f                                    |




#### `manifest/`


| File                                  | Contents                                   |
| ------------------------------------- | ------------------------------------------ |
| `BUILD_DATE.txt`                      | ISO timestamp when this repo was assembled |
| `Costus_popgen_Dryad_20260722.sha256` | SHA-256 of the deposit zip                 |




### Abbreviations


| Abbreviation                    | Meaning                                                                         |
| ------------------------------- | ------------------------------------------------------------------------------- |
| aHMM                            | Ancestry_HMM local ancestry inference                                           |
| bp / cM                         | Base pairs / centiMorgans                                                       |
| CI                              | Credible interval (Bayes; R/qtl `bayesint`)                                     |
| Fst / F~ST~                   | Weir–Cockerham fixation index (PIXY `avg_wc_fst`)                               |
| Dxy / D~xy~                   | Absolute divergence (PIXY)                                                      |
| π / pi                          | Nucleotide diversity (PIXY)                                                     |
| fdM / *f*_dM                    | Dsuite introgression statistic                                                  |
| SD, GR, tough, Nreabs, Chlreabs | Seed dormancy, growth rate, toughness, N reabsorption, chlorophyll reabsorption |
| HK                              | Haley–Knott regression (R/qtl)                                                  |
| LASI / LASIref                  | *C. lasius* reference genome                                                    |
| alle / vill / brac              | *allenii* / *villosissimus* / *bracteatus* population codes                     |




### Details that affect reuse

- Chromosome **4** is often excluded from cM-scaled QTL vs Fst summaries because of poor Lep-MAP3 joining (physical-window plots may still include it).
- QTL models use trait-specific covariates (hybrid index, selected F2 families, cohorts); see the manuscript and R/qtl scripts.
- Seed dormancy enters QTL models as `log_days_to_radicle`, not raw days.
- This package is intentionally **minimal**: regenerating BAMs/VCFs/posteriors from SRA is expected for full pipeline reproduction.

---



## Sharing/Access information

Links to other publicly accessible locations of the data:

- Publication code (GitHub): [https://github.com/jharencar/Costus_popgen_clean](https://github.com/jharencar/Costus_popgen_clean)
- Reference genome (*C. lasius*): NCBI [GCA_027563935.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027563935.1/)
- Low-coverage WGS FASTQs (wild + mapping population): NCBI SRA BioProject [PRJNA1479894](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1479894)
- Dryad DOI: *(add after deposit)*
- Related manuscript: *(add citation / DOI when available)*

Data were derived from the following sources:

- Newly generated low-coverage whole-genome sequencing of wild *Costus* and the F3 mapping population (SRA BioProject above), aligned to *C. lasius* GCA_027563935.1
- Phenotypes measured on the greenhouse F3 mapping population (methods in the associated manuscript)
- Downstream summary statistics (PIXY, Dsuite) computed from filtered genotype data produced by the pipelines in `code/`

**License note:** Dryad releases data under Creative Commons Zero (CC0). Prefer linking to the GitHub repository for code under whatever license is declared there rather than treating scripts as independently relicensed beyond Dryad’s terms for the deposited files.

---



## Code/Software



### Recommended environment variables

After unpacking this repo:

```bash
export COSTUS_ARCHIVE="/path/to/Costus_Dryad_package"
export PROJECT_ROOT="${COSTUS_ARCHIVE}/code"
export GENOPROB_RDS="${COSTUS_ARCHIVE}/data_qtl/05_rqtl/aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds"
export PIXY_FOLDER="${COSTUS_ARCHIVE}/data_figure_inputs/pixy_output_251124_LASIref"
export DSUITE_DIR="${COSTUS_ARCHIVE}/data_figure_inputs/Dsuite_Fstats_w50s25_LASIref_251120"
export SCANTWO_PERM_DIR="${COSTUS_ARCHIVE}/data_qtl/scantwo.perms"
```

Then in R from `code/`:

```r
source("scripts/project_paths.R")
paths <- project_paths()
```



### Primary tools


| Stage               | Software (representative)                                                   |
| ------------------- | --------------------------------------------------------------------------- |
| Read QC / trim      | Trimmomatic (v0.39), BBTools dedupe                                         |
| Align               | BWA-MEM                                                                     |
| Variants / popgen   | bcftools; ANGSD; PCAngsd; Ancestry_HMM (+ `vcf2aHMM` helper); Dsuite; PIXY  |
| Linkage map         | Lep-MAP3                                                                    |
| QTL                 | R and **R/qtl** (`scanone`, `scantwo`, `stepwiseqtl`, `fitqtl`, `bayesint`) |
| Figures / wrangling | R (`tidyverse` / `ggplot2` stack as used in plotting scripts)               |


---



## Citation

*(Add Dryad DOI and manuscript citation after deposit / acceptance.)*