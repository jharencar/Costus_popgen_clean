# Costus allenii × C. villosissimus — population genomics and QTL

[Access this dataset on Dryad](DOI_PLACEHOLDER)

This is the **data** deposit for analyses in the manuscript, "Divergent habitat selection across many loci maintains species boundaries during introgression." The analyses investigate ecogeographic isolation, introgression, and QTL for habitat-associated traits in Neotropical *Costus allenii* and *C. villosissimus*.

## How to use the Dryad and Zenodo deposits together

This manuscript is archived as **two linked deposits** so data and code can carry appropriate licenses:

| Deposit | Contents | License |
| ------- | -------- | ------- |
| **This Dryad deposit** | Data tables and expensive-to-recompute QTL permutation objects (listed below) | **CC0** (Dryad policy) |
| **Linked Zenodo software deposit** | Analysis scripts and documentation only (`01_`–`04_`, `scripts/`, docs, `LICENSE`) | **MIT** |

**No data files are duplicated between the two deposits.** Phenotype / genotype / QTL-CI CSVs, the sample manifest, PIXY and Dsuite window tables, and `scantwo` permutation `.rds` files live **only on Dryad**. Zenodo contains **code only**.

### Recommended workflow

1. Download **this Dryad data package** and unpack it (e.g. to `Costus_Dryad_data/`).
2. Download the **Zenodo code archive** (linked from this Dryad record) and unpack it (e.g. to `Costus_popgen_clean/`).
3. Point the code at the data with environment variables (see **Code/Software**), then run scripts from the code tree. For QTL reuse: use the Dryad CSVs in `data_qtl/small_inputs/` with the R/qtl scripts in the Zenodo archive; use Dryad `data_qtl/scantwo.perms/` for stepwise QTL penalties; use Dryad `data_figure_inputs/` for genome-wide figure scripts in `04_introgression_vs_QTL/`.
4. For full pipeline regeneration from raw sequence, also download NCBI SRA BioProject [PRJNA1479894](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1479894) and the *C. lasius* reference [GCA_027563935.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027563935.1/).

Development copies of the code also live on GitHub: [https://github.com/jharencar/Costus_popgen_clean](https://github.com/jharencar/Costus_popgen_clean). Prefer the Zenodo archive for a citable, frozen snapshot.

**Not included in either deposit** (regenerable): large R/qtl cross objects (`*_genoprob.rds`, `est.rf.rds`); raw FASTQs, BAMs, VCFs, Ancestry_HMM posteriors, and other large intermediates.

Package build date: `2026-07-27`. A SHA-256 checksum of the data zip is provided alongside the download (`Costus_popgen_Dryad_data_20260727.zip.sha256`).

---



## Description of the data and file structure



### Study system and experimental summary

- **Species:** *Costus allenii*, *C. villosissimus*, wild hybrids; *C. bracteatus* (ingroup for ABBA-BABA / Dsuite); *C. lasius* (outgroup and reference genome).
- **Wild WGS:** Central Panama collections (plus *C. bracteatus* from Costa Rica; *C. lasius* from Panama and Peru). 
- **Mapping population:** F3 progeny of an *allenii* × *villosissimus* cross (plus F1 parent and F2s used in linkage mapping). QTL mapping retains **385** genotyped individuals after QC.
- **Traits scored in the F3 population (habitat-associated):** seed dormancy (days to radicle), growth rate (cm/day), leaf toughness and thickness, percent nitrogen and chlorophyll reabsorption during drought-induced senescence, and leaf nitrogen content.
- **Analyses covered by this deposit + code archive + NCBI sequences:** low-coverage WGS processing; wild population genetics (PCA, Ancestry_HMM, ABBA-BABA, Dsuite, PIXY); F3 Ancestry_HMM → Lep-MAP3 → R/qtl; comparison of differentiation / introgression statistics with QTL intervals.

Full methods prose is in the manuscript, "Divergent habitat selection across many loci maintains species boundaries during introgression."

### Package layout

```
Costus_Dryad_data/
  sample_manifest/
    combined_sample_manifest.tsv  # sample metadata (wild + mapping population)
  data_qtl/
    scantwo.perms/                # scantwo permutation RDS by trait
    small_inputs/                 # phenotype / genotype / QTL CI CSVs
  data_figure_inputs/
    pixy_output_251124_LASIref/   # PIXY π, Dxy, Fst (10 kb windows)
    Dsuite_Fstats_w50s25_LASIref_251120/  # Dsuite local F-statistics (w50 s25)
  manifest/                       # build date
```

This README is the single documentation file for the paired deposits (entered in Dryad's
README field; not duplicated inside the data zip or the Zenodo code zip). Analysis **code**
is only on Zenodo/GitHub. Large R/qtl cross objects (`*_genoprob.rds`, `est.rf.rds`) are
omitted and are regenerated from the CSVs in `data_qtl/small_inputs/` by the R/qtl scripts.

### How the pieces relate

1. Use the **Zenodo code archive** with NCBI FASTQs (and the reference genome) to regenerate pipelines from scratch, **or**
2. Point the code's `scripts/project_paths.R` at **this Dryad data deposit** (see **Code/Software**) to reload the figure-input tables, QTL CSVs, and permutation objects without re-running cluster jobs.
3. `data_qtl/small_inputs/` are the tabular phenotype/genotype/CI inputs used to build / annotate the R/qtl cross objects (the cross objects themselves are regenerated, not deposited).
4. `sample_manifest/combined_sample_manifest.tsv` is the sample metadata table (wild + mapping population); **Dryad only**.
5. `data_figure_inputs/` are the windowed summary tables read by the introgression-vs-QTL and genome-wide figure scripts in the Zenodo code archive (`04_introgression_vs_QTL/`).

### File descriptions

#### Code (Zenodo only — not in this zip)

The analysis code is archived on **Zenodo** under the **MIT license** (DOI in **Sharing/Access information**) and developed openly at [https://github.com/jharencar/Costus_popgen_clean](https://github.com/jharencar/Costus_popgen_clean). The code tree contains read-processing, population-genetics, QTL-mapping, and introgression-vs-QTL scripts (`01_`–`04_`), `scripts/project_paths.R`, and documentation (`README.md`, `NOTES.md`, `LICENSE`). It does **not** include any of the data files listed below.

#### `sample_manifest/combined_sample_manifest.tsv`

Final sample metadata for wild collections and the QTL mapping population (**524** rows; one row per unique `sample_id`). This is the only sample table in the deposit (build scripts and NCBI source tables are not included).


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

#### R/qtl cross objects (`*_genoprob.rds`, `est.rf.rds`) — regenerable, not deposited

The R/qtl `cross` objects are **not** included in this deposit because they are large and
fully regenerable from the deposited inputs. Rebuild them with the R/qtl scripts in the code
archive from `data_qtl/small_inputs/` (genotype + phenotype CSVs): import the cross, then
`est.rf()` and `calc.genoprob()`. The canonical published object corresponds to the moderate
duplicate-marker thinning (`most_dups`, ≈14,889 markers; 385 individuals) used for `fitqtl` /
`refineqtl` / multi-QTL models; a more aggressively thinned `least_dups` cross (≈4,013 markers)
is used for some scans.

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

Wide genotype table for the mapping population (≈391 rows including header alignment with IDs; ≈36,783 marker columns before R/qtl QC thinning).


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

This deposit is the **Weir–Cockerham** 10 kb run (`avg_wc_fst`). The current plotting scripts do not read it. They expect a later Hudson rerun (`--fst_type hudson --fst_components`) at 10 kb and 500 kb, with `avg_hudson_fst`, `hudson_fst_num`, and `hudson_fst_den`. Those tables are not in this deposit. Point `PIXY_FOLDER` and `PIXY_FOLDER_500KB` at them when they exist.

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


| File             | Contents                                 |
| ---------------- | ---------------------------------------- |
| `BUILD_DATE.txt` | ISO date when this deposit was assembled |




### Abbreviations


| Abbreviation                    | Meaning                                                                         |
| ------------------------------- | ------------------------------------------------------------------------------- |
| aHMM                            | Ancestry_HMM local ancestry inference                                           |
| bp / cM                         | Base pairs / centiMorgans                                                       |
| CI                              | Credible interval (Bayes; R/qtl `bayesint`)                                     |
| Fst / F_ST                      | This deposit: Weir–Cockerham (`avg_wc_fst`). Current code: Hudson (`avg_hudson_fst`) |
| Dxy / D_xy                      | Absolute divergence (PIXY)                                                      |
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
- This deposit is intentionally **minimal**: the analysis code (Zenodo/GitHub), the large R/qtl cross objects, and BAMs/VCFs/posteriors are regenerated (from SRA and the deposited CSVs) rather than stored here.

---



## Sharing/Access information

Links to other publicly accessible locations of the data:

- Analysis code (Zenodo archive, MIT license): *(add Zenodo DOI after deposit)*
- Analysis code (development repository, GitHub): [https://github.com/jharencar/Costus_popgen_clean](https://github.com/jharencar/Costus_popgen_clean)
- Reference genome (*C. lasius*): NCBI [GCA_027563935.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027563935.1/)
- Low-coverage WGS FASTQs (wild + mapping population): NCBI SRA BioProject [PRJNA1479894](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1479894)
- Dryad DOI (this data deposit): *(add after deposit)*
- Related manuscript: *(add citation / DOI when available)*

Data were derived from the following sources:

- Newly generated low-coverage whole-genome sequencing of wild *Costus* and the F3 mapping population (SRA BioProject above), aligned to *C. lasius* GCA_027563935.1
- Phenotypes measured on the greenhouse F3 mapping population (methods in the associated manuscript)
- Downstream summary statistics (PIXY, Dsuite) computed from filtered genotype data produced by the pipelines in the code archive

**License note:** This Dryad data deposit is released under Creative Commons Zero (**CC0**), per Dryad policy. The analysis **code** is licensed separately under the **MIT License** and archived on Zenodo (with a linked DOI) and on GitHub; refer to those deposits for code reuse terms rather than treating the scripts as CC0.

---



## Code/Software



### Recommended environment variables

After unpacking **this Dryad data deposit** and the **Zenodo code archive**:

```bash
export COSTUS_ARCHIVE="/path/to/Costus_Dryad_data"       # this Dryad data deposit
export PROJECT_ROOT="/path/to/Costus_popgen_clean"        # Zenodo code archive
# Phenotype / genotype / CI CSVs live only on Dryad:
export QTL_CI_PATH="${COSTUS_ARCHIVE}/data_qtl/small_inputs/all_QTL_CIs_bp_cm.251130.csv"
export SAMPLE_MANIFEST="${COSTUS_ARCHIVE}/sample_manifest/combined_sample_manifest.tsv"
# Cross object is not deposited; point this at one you regenerate from the Dryad CSVs:
export GENOPROB_RDS="/path/to/your_regenerated_most_dups_genoprob.rds"
# Deposited 10 kb tables are Weir–Cockerham and are refused by the current scripts:
# export PIXY_FOLDER="${COSTUS_ARCHIVE}/data_figure_inputs/pixy_output_251124_LASIref"
# Current scripts (Hudson rerun; not in this deposit):
export PIXY_FOLDER="/path/to/pixy_hudson_10kb"
export PIXY_FOLDER_500KB="/path/to/pixy_hudson_500kb"
export DSUITE_DIR="${COSTUS_ARCHIVE}/data_figure_inputs/Dsuite_Fstats_w50s25_LASIref_251120"
export SCANTWO_PERM_DIR="${COSTUS_ARCHIVE}/data_qtl/scantwo.perms"
```

Copy the Dryad CSVs from `data_qtl/small_inputs/` into a working directory (or set paths
explicitly) when running `03_QTL_mapping/rqtl/Costus_QTL.qmd`, which historically expected
`F3phe_onehotcovar_newHI_20250925.csv` and `f3gen_250924_...csv` in the working directory.

Then in R from the code archive:

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