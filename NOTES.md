# Project notes (Costus_popgen_clean)

This file captures decisions and work from development sessions so future tooling and agents stay aligned with local conventions. **If you are an agent, read this file when working in this repository.**

## Methods summary (manuscript)

Full methods prose is in the associated manuscript. Below: data types, sample sets, and filters agents need when writing or modifying analysis/plotting scripts.

### Study system and samples

- **Focal species:** *Costus allenii*, *C. villosissimus*, wild hybrids between them.
- **Tests / outgroups:** *C. bracteatus* (ingroup for ABBA-BABA and Dsuite); *C. lasius* (outgroup; also **reference genome**).
- **Wild population genomics (collected):** 38 *allenii*, 23 *villosissimus*, 12 hybrids (Central Panama); plus 21 *bracteatus* (Costa Rica), 2 *lasius* (Panama + Peru).
- **Mapping population:** 390 F3 (+ wild F1 parent, 21 F2s including five F2 parents of F3s). F3 minimizes maternal/environmental noise vs field; **cohort** tracks shared greenhouse conditions.
- **Coverage removals (all downstream):** wild samples with **<0.15×** mean coverage or **<10%** breadth removed (**10** *allenii*, **3** *villosissimus*); one F3 with **<0.05×** removed → post-filter wild totals **28** *allenii*, **20** *villosissimus* (before PCA-hybrid exclusion or analysis-specific subsampling).

### Reference, sequencing, and read processing

| Step | Setting |
|------|---------|
| Reference | *C. lasius* — NCBI **GCA_027563935.1** (avoids reference bias vs focal clade) |
| Sequencing | Low-coverage WGS; **150 bp** paired-end (wild: Fulgent; mapping pop + 6 wild: UC Davis) |
| Trim / dedupe | Trimmomatic v0.39; BBTools dedupe |
| Align | BWA mem → *C. lasius* |
| Read filters | Unmapped removed; **min length 35 bp**; **MAPQ ≥ 20** |

### Core variant filters (wild popgen pipeline)

Apply unless an analysis block below overrides.

- **Mappability:** genmap, **k = 150**, merge distance **100 bp**, keep regions **>= 0.7**.
- **SNP calling:** base + map qual **20**; **MAF 0.05**; sites **within 3 bp of indel** removed; depth **> 99th percentile** removed.

### Analysis-specific filters and data types

#### PCA (recent hybridization)

- Extra SNP filters: **p-value 1e-6**; base/map qual **30**; **>= 1×** in **>= 80%** of samples per species **and** hybrids.
- **ANGSD** genotype likelihoods (major/minor inferred) → **PCAngsd** PCA.
- Hybrid individuals in PCA are used for visualization; **excluded** from ABBA-BABA and Dsuite *f*~dM~ (see below).

#### aHMM — wild hybrids (local ancestry)

- **Parent panels:** **18** *villosissimus* + **24** *allenii* (wild references; not necessarily identical to post-coverage wild counts).
- **vcf2aHMM** (`new_vcf2aHMM.py`): mean recomb **0.7**; min **3** genotype calls per parent pop; min mean depth admixed **0.5**; **5000 bp** between sites; min AF difference **0.7**.
- **Run:** uniform error **0.02**; single bidirectional pulse (**50%** each parent) at time **0.0001** (disfavors tiny blocks).
- **Repo plotting convention (not in methods prose):** posterior TSV → keep sites with **ref_homo > 0.9 | alt_homo > 0.9 | het > 0.9**; ancestry score = `ref_homo + het * 0.5` (see aHMM session below).

#### ABBA-BABA — ancient introgression (wild)

- **Hybrids removed** (PCA hybrids) → **28** *allenii* + **18** *villosissimus*.
- **ANGSD** `-doabbababa2`: **5 Mb** blocks; qual **30**; ingroup **21** *bracteatus*; outgroup *C. lasius* reference; per-chromosome **FDR**.

#### Dsuite — Ecuador capture test

- **756-gene capture**, hard genotypes; **Dtrios**; **25,978** variants; **20** jackknife blocks (~1297 variants/block).
- Pops: 1 *bracteatus* (CR), 1 *allenii* (Panama), 1 *villosissimus* (Ecuador), 2 *lasius* (Panama + Peru).

#### PIXY — π, D~xy~, and Hudson F~ST~

- **All-sites VCF** (invariant sites included): bcftools mpileup + core filters above + **<80% missing per species**.
- Windows: **10 kb** and **500 kb**, both non-overlapping. π and D~xy~ are pixy's ratio of `count_diffs / count_comparisons`.
- **F~ST~:** Hudson estimator (`--fst_type hudson`), which is the two-population estimator recommended by Bhatia et al. 2013. Each window's `avg_hudson_fst` is already the ratio of summed numerator and denominator. `--fst_components` keeps `hudson_fst_num` and `hudson_fst_den`.
- **500 kb line** on the along-genome plots is the 500 kb pixy file, not a mean of 10 kb windows. A 500 kb window is drawn when at least half of its 10 kb slots pass the filter below. Adjacent kept windows are connected; a dropped window breaks the line.
- **10 kb filter (all current F~ST~ and D~xy~ plots):** F~ST~ windows need **≥10 SNPs**. D~xy~ windows need **≥1000 sites** (10% of a 10 kb window). The 1 cM and 0.1 cM summaries are built from those passing 10 kb rows. Pi is not filtered.
- **1 cM F~ST~** and QTL-interval F~ST~ are `sum(hudson_fst_num) / sum(hudson_fst_den)` over the passing 10 kb rows in the bin. Pixy cannot emit centimorgan windows. 1 cM D~xy~ is the mean of the passing 10 kb `avg_dxy` values.
- **Original-submission Weir–Cockerham tables** (`output_251124_LASIref`, column `avg_wc_fst`) are the November 2025 run. Current scripts refuse that column. Snapshot of those scripts: `_originals/wc_fst_original_submission/`. The Hudson scripts are the version that replaces that submission.
- **Repo:** `02_population_genetics_wild/pixy/run_pixy_array.slurm` (Farm; submit at both window sizes), `genome-wide_Fst.Dxy_w_QTL_plots.R`, pixy chromosome plots; combo figures stack pixy with Dsuite (shared cumulative x).

#### Dsuite — *f*~dM~ (introgression sliding windows)

- **50 variant** windows, **25 variant** step (`w50s25` in paths/filenames); **Dinvestigate**; **biallelic SNPs** only; **>3 bp** from indels; AF **0.05–0.95**; same qual/coverage/mappability as core pipeline.
- **Samples:** **21** *bracteatus*, **21** *allenii*, **18** *villosissimus*, **2** *lasius* + *lasius* reference; **no PCA hybrids**.
- **Repo only:** highlight extreme |*f*~dM~|, *d*~f~, D~xy~, and F~ST~ at the **99th percentile** (two-tailed stats use the 1st and 99th percentiles). The older 99.9th-percentile highlight is retired.

#### F~ST~ vs introgression / QTL (repo extensions)

- Manuscript QTL comparison: **10 kb** Hudson F~ST~ windows + **1 cM** ratios of Hudson components; **chr4 excluded** from cM mapping (poor Lep-MAP3 join). The along-genome line is a separate **500 kb** pixy run.
- **Repo:** `fd_fdm_vs_fst_100kb_bins.R` — **100 kb** bins, pixy SNP-weighted F~ST~ from 10 kb tiles + Dsuite *f*~d~/*f*~dM~/*d*~f~ with **stride-2** (or greedy) thinning on w50s25 windows (not in manuscript).

### QTL / mapping population (shorter)

#### aHMM (F3)

- Same parameters as wild hybrids; F3 + pedigree used for Lep-MAP3 input.

#### Allele coding (AA / BB / hybrid_index / aHMM plots) — locked 2026-07-14/15

**Resolved.** Mapping genotypes used Ancestry_HMM with  
`--pop ../F3s0.05covCO_allF2s_F1parent.pop`  
([`03_QTL_mapping/aHMM_F3/512x64_LASIref_ahmm_pipeline_F3s_F2parents_F1parent.sh`](03_QTL_mapping/aHMM_F3/512x64_LASIref_ahmm_pipeline_F3s_F2parents_F1parent.sh)).  
User confirmed that `.pop` has **parent1 = *C. villosissimus***, **parent2 = *C. allenii***.

**Pipeline → symbol map**

| Source | Villosissimus | Allenii |
|--------|---------------|---------|
| `.pop` | parent1 | parent2 |
| aHMM panel / HMM state | state 0 | state 1 |
| Posterior column | `ref_homo` | `alt_homo` |
| Hard call (Lep-MAP / R/qtl) | **AA** | **BB** |
| `hybrid_index` | high (≈ genome-wide vill) | low |
| Strip / posterior **plot y** (`ref_homo + 0.5×het`) | **≈ 1** | **≈ 0** |

Hard-call rule: `ref_homo > 0.9 → AA`, `alt_homo > 0.9 → BB`.  
HI tracks AA dosage (*r*≈0.92). `Costus_QTL.qmd` “HI = vill ancestry” is correct; ignore older “HI = allenii” comments.

**Plot y ≠ HMM state ID.** Internal state **0** is parent1 (vill); plotting maps P(homozygous state 0) to **y ≈ 1**. So “1 on the plot = homozygous vill” and “parent1 = vill” agree.

**Sanity check (AHV_LASIref2_260210):** labeled `19.587.726_VILL` and high-vill PCA hybrids **19.589** / **19.727** have high mean plotted ancestry and many `ref_homo` (plot ≈1) sites; allenii-leaning hybrids have more `alt_homo` (plot ≈0). Confirms plot 1 / `ref_homo` / parent1 = vill.

**R/qtl `fitqtl`:** **a = (μ_AA − μ_BB)/2**, **d = μ_AB − mid-parent**.  
Positive **a** = vill higher than allenii *at that QTL after covariates*. Manuscript additive caption must say substituting *allenii* → *vill* (not the reverse).  
For reabsorption, genome-wide **HI coef is positive** (vill↑ → reabs↑) while local QTL **a can be negative** (conditional on HI); use HI/raw means for the species check, not table `a` alone. Raw AA/BB means are in [`QTL/QTL_peaks_table.tex`](QTL/QTL_peaks_table.tex) (also under Google Drive `Figures&Tables/`).

Full investigation + F2-family/HI exploration + draft LaTeX QTL table (not used in MS):  
[`_unpublished_local/allele_coding_HI_background/`](_unpublished_local/allele_coding_HI_background/) (gitignored).

#### Lep-MAP3 linkage map

- Data: wild F1, **21** F2s, **389** F3s (one low-coverage F3 dropped vs 390 collected).
- aHMM genotypes → hard calls **>90%** prob; **LOD 11**, theta **0.001**, `distortionLod` **11**.
- Manual fixes: **chr4** three LGs joined; small centromeric translocations removed on chr2/chr5 → **9** linkage groups.
- **OrderMarkers2** ×100 per chr; best likelihood run kept.

#### R/qtl

- Drop individual **<80%** markers; markers in **<80%** individuals; dup thinning → **14,889** markers; remove **3** duplicate-genotype individuals → **385** individuals.
- **scanone** (Haley–Knott) + **10% LOD FDR** (1000 permutations); **stepwiseqtl** with penalties from **scantwo** perms; **bayesint** 95% CIs.
- Covariates (AIC stepwise): **hybrid index**, **F2 family**, **cohort** (one-hot where kept).
- **stepwiseqtl** tractability: extra dup removal; genoprob at **5 cM** steps, error **0.01**; final **fitqtl** on full marker genoprob with **refineqtl** peaks.
- Trait notes: log(seed dormancy); growth rate / chlorophyll reabsorption have manuscript exceptions to strict stepwise-only reporting.

#### Relating introgression to QTL

- F~ST~ in **10 kb** windows; QTL **95% Bayes CIs** overlaid; compare genome-wide **1 cM** mean F~ST~ to QTL windows **<15 cM** (~10 Mbp) — excludes very large CIs (e.g. N reabsorption chr2, growth rate chr5).

#### QTL / aHMM `.rds` cross objects (reference)

These are **rQTL `cross` objects** (or scanone/scantwo perm outputs), built from aHMM hard calls in `QTL/f3gen_250924_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001.csv` + phenotypes in `QTL/F3phe_onehotcovar_newHI_20250925.csv`. Documented in [`Costus_QTL.qmd`](Costus_QTL.qmd). **Most live outside this repo** (e.g. `alle_vill_QTL/` on the author machine); paths below are filenames—adjust `readRDS()` paths locally.

**Shared filename stem:** `aHMM_<date>_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_<dup_level>_<stage>.rds`

| Suffix | Stage | What it contains |
|--------|--------|------------------|
| `*_est.rf.rds` | After `est.rf()` | Hard genotypes + recombination fractions estimated; **no** genotype probabilities yet |
| `*_genoprob.rds` | After `calc.genoprob()` | Genotype **probabilities** at each locus (for Haley–Knott); use as `Costus_prob` in mapping |

**Duplicate-removal level (`<dup_level>`):**

| Tag | Pipeline step | ~Markers (order of magnitude) | Typical use |
|-----|----------------|--------------------------------|-------------|
| **`more_dups`** | First dup pass: `findDupMarkers(exact.only=TRUE, adjacent.only=TRUE)` on filtered cross; then individual dedup → `est.rf` | ~**14,889** (per qmd after first thinning) | Parent of genoprob; saved as `*_more_dups_est.rf.rds` |
| **`most_dups`** | **`calc.genoprob(Costus)`** on the `more_dups` / moderate-clean cross | **Many** (full moderate-clean map; same analysis set as final QTL) | **Canonical for final QTL:** `refineqtl`, `fitqtl`, `stepwiseqtl`, GRM/h² (planned). Loaded as `Costus_prob` in qmd line ~1004 |
| **`least_dups`** | Second pass on genoprob: `findDupMarkers(exact.only=FALSE, adjacent.only=FALSE)` + drop | ~**4,013** | Smaller cross for some **scanone/scantwo perms** and stepwise exploration; **not** the primary object for published multi-QTL `fitqtl` models |
| **`exact.F_adjacent.T`** | Alternate dup pass on `est.rf` object (adjacent-only, after Oct 2024 dup logic) | Intermediate | Saved in qmd; superseded by `least_dups` for some scans |

**Naming quirk:** In `Costus_QTL.qmd`, one chunk **saves** `cleaned_more_dups_genoprob.rds` but **reads** `cleaned_most_dups_genoprob.rds` for `Costus_prob`. Treat **`most_dups_genoprob`** as canonical; confirm on disk whether `more_dups` and `most_dups` genoprob files are duplicates or only one exists.

**Canonical genoprob path (author machine):**

`/Users/admin-jgharenc/Documents/Github/alle_vill_QTL/aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001_cleaned_most_dups_genoprob.rds`

Also referenced in [`QTL/plot_GR_Nreabs_qtl_minimal.R`](QTL/plot_GR_Nreabs_qtl_minimal.R). Embedded phenotypes include `hybrid_index`, `F21_*`, cohorts, traits, and `log_days_to_radicle`.

**Building pipeline (short):**

1. `read.cross` → drop sparse inds/markers → first dup thinning → drop 3 near-duplicate individuals (`Costus.clean`).
2. `est.rf` → save `*_more_dups_est.rf.rds` (or `*_251006_*` dated variant).
3. `calc.genoprob` → **`most_dups_genoprob.rds`** ← **use this for mapping + heritability**.
4. Optional: aggressive dup drop → `least_dups_genoprob.rds` (~4013 markers).

**Other QTL-related `.rds` files:**

| File (location) | Contents |
|-----------------|----------|
| `log.seed.dormancy.scanone-hk_1kperm_251122.rds`, `growth.rate.scanone-hk_1kperm_251120.rds`, etc. | **scanone** permutation results (1000 perms) for 10% LOD FDR cutoffs |
| `scantwo.perms/*.scantwo.hk.perm.rds` | **scantwo** permutations for `calc.penalties` / `stepwiseqtl` |
| `QTL/circular_shift_perm_QTL_CI_Fst_fdM.rds` | Circular-shift null for QTL CI vs Fst/fdM overlap (in repo) |

**Agent defaults:** For new QTL or genomic-heritability scripts, read **`most_dups_genoprob.rds`** unless explicitly reproducing an old scan that used **`least_dups`**. Do not assume `.rds` files are in `Costus_popgen_clean`—set a configurable path at the top of scripts.

---

## Session: genomic heritability and genetic variance (QTL traits)

### Scripts (repo `QTL/`)

- **`fitqtl_summary_utils.R`** — shared parser for `*_QTL_model_output_251203` (used by variance plot, partition, heritability audit).
- **`QTL_model_variances_plot.R`** — `% phenotypic variance` by component; ggsave **5 × 1.8 in** (shortened from 2.5).
- **`genomic_heritability_QTL_traits.R`** — GRM + `sommer` REML; env `GENOPROB_RDS`; optional project **`R_libs/`** for dependencies.
- **`genetic_variance_partition.R`** — SS-based supplementary partition from fitqtl drop-one tables.
- **`build_genetic_architecture_summary.R`** — merges REML + partition → **`genetic_architecture_summary.csv`**.
- **`genetic_architecture_summary_column_metadata.csv`** + **`genetic_architecture_summary_DATA_DICTIONARY.md`** — column definitions and formulas for the summary table.

### Conventions

- **Primary** narrow-sense \( \hat{h}^2 = \hat{\sigma}^2_g / (\hat{\sigma}^2_g + \hat{\sigma}^2_e) \) conditional on the **same fixed covariates** as each trait’s final `fitqtl` model (not QTL genotypes in the REML fixed effect).
- **Supplementary** `h2_fixed_approx` from summed Type-III `%var` (QTL + hybrid_index + F21) — can exceed 1; do not treat as \( h^2 \).
- **`frac_genetic_QTL`** = `pct_QTL` / (`h2` × 100).
- Sample audit: `QTL_trait_sample_audit.csv` (complete-case `n` must match fitqtl `n` for all five traits).

**SD + hybrid_index (secondary):** `QTL/SD_secondary_hybrid_index_analysis.R` refits SD QTL with `hybrid_index` added to seed-cohort covariates. Outputs: `SD_QTL_model_output_251203_with_hybrid_index`, `SD_hybrid_index_secondary_comparison.csv`, `genetic_architecture_summary_SD_with_HI.csv`. In this run, HI was **not** significant after cohorts in nested `lm` (p ≈ 0.55) and **increased** AIC; QTL `%var` and REML `h2` changed only slightly vs primary.

### Manuscript vs repo — quick flags

| Topic | Manuscript | Repo / NOTES |
|-------|------------|--------------|
| Wild *villosissimus* after coverage filter | 20 (38−10−3 allenii only listed) | ABBA-BABA / *f*~dM~ use **18** (hybrid exclusion + analysis set) |
| F3 count | 390 collected; **389** in Lep-MAP | **385** in final QTL cross |
| aHMM posterior filter | Not stated | **>0.9** on ref/alt/het (plotting) |
| Outlier highlight | **99th** percentile | F~ST~, D~xy~, *f*~dM~, and *d*~f~ figures all use the **99th** percentile |
| F~ST~ vs *f*~dM~ correlation | — | **100 kb** bins + stride-2 (`fd_fdm_vs_fst_100kb_bins.R`) |
| PIXY / F~ST~ | **10 kb** Weir–Cockerham in the 251124 deposit | Current code: Hudson **10 kb** points, Hudson **500 kb** line, 1 cM F~ST~ = ratio of components |

---

## Session: aHMM genome-wide strip plot (AHV LASIref2 260210)

### Deliverable

- **`aHMM/plot_AHV_LASIref2_all_individuals_genomewide_stripes.R`** — standalone script (also sourced from `aHMM/popgen_aHMM.qmd`). Writes **`aHMM/AHV_LASIref2_260210_all_individuals_genomewide.pdf`** (path adjusts if cwd is repo root vs `aHMM/`).

### Posterior handling (matches older qmd convention)

- Read TSV posteriors; rename columns to `ref_homo`, `het`, `alt_homo`.
- Filter: `ref_homo > 0.9 | alt_homo > 0.9 | het > 0.9`.
- **Ancestry** for plotting: `ref_homo + het * 0.5`.

### Plot design

- **Layout:** One row per individual (`ggh4x::facet_grid2` + `strip_themed`); cumulative x = concatenated chromosomes (chr1–chr9).
- **Strip labels:** File stems like `19.516_HYB` are shown as **numeric only** (`gsub("[^0-9.]", "", ...)`) to narrow the strip column.
- **Row background:** Alternating **grey95** / **white** for strip + panel (`geom_rect` per row + matching strip backgrounds). Point colors by chr alternate **grey35** / **black** (not full-panel stripes; earlier version used grey/black panel fills and was replaced).
- **x-axis:** Chromosome tick labels at chr midpoints (1–9); **vertical lines** at chr boundaries; **y-axis** breaks `0`, `0.5`, `1`; **`coord_cartesian(ylim = c(0, 1))`** without `scale_y` `limits` (avoids dropping points at boundaries).
- **Size:** PDF width **17.8 cm**; height scales with number of individuals (`min(26, max(14, n * 1.9))` cm).
- **Dependency:** **`ggh4x`** (CRAN) for per-strip styling.

### Qmd

- Section in `aHMM/popgen_aHMM.qmd` documents the figure and sources the script with `source("plot_AHV_LASIref2_all_individuals_genomewide_stripes.R")` (cwd when rendering is typically `aHMM/`).

---

## Session: combo Dsuite + pixy four-panel figure

### File

- **`combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R`** — panels: **f_dM**, **d_f**, **Dxy**, **Fst** (no pi). `project_root` is set at the top of the script (adjust for your machine).

### Shared x-axis (aligned with other genome-wide scripts)

- After both Dsuite and pixy cumulative positions exist:
  - `dsuite_xlim <- range(LASIout_w50s25_all9$windowStart_cum, na.rm = TRUE)`
  - `pixy_xlim <- range(pixy_wide$windowStart_cum, na.rm = TRUE)`
  - **`combo_xlim`** = min of the two mins, max of the two maxes (pixy can extend slightly past Dsuite at the right).
- All four panels use the same **`scale_x_continuous(..., limits = combo_xlim, expand = x_expand)`** so patchwork stacks line up.
- **Removed** hardcoded **`coord_cartesian(xlim = c(24398, 1120430367))`** on **d_f** (it clipped real data). **d_f** x labels use **`chr_label`** like **f_dM**.

### Horizontal “side gap”

- Controlled by **`x_expand <- expansion(mult = ..., add = 0)`**. Started from genome-wide convention **0.018**; tightened to **0.012** for slightly less margin at chr1 / chr9 edges (tune this one value if needed).

### Panel heights (2026-11-11)

- `plot_layout(heights = c(0.9, 0.9, 1.2, 1.45, 0.08))` for **f_dM**, **d_f**, **Dxy**, **Fst**, legend.
  Was `c(0.55, 0.55, 1.45, 1.45, 0.08)` — introgression panels taller, Dxy a bit shorter.

### Warnings

- Large **`geom_point` “Removed … rows”** messages are often **missing y** (e.g. `avg_dxy` / `avg_wc_fst` NA), not necessarily x limits.

---

## Session: Dsuite introgression vs F_ST (100 kb bins, Spearman)

### Deliverable

- **`fd_fdm_vs_fst_100kb_bins.R`** — Shared **100 kb** grid: pixy **SNP-weighted** `avg_wc_fst` (10 kb tiles split at bin edges); Dsuite **`f_d`**, **`f_dM`**, **`d_f`** with **stride-2** (or **`greedy`**) + overlap weights. Modes: **`all`**, **`positive`**, **`positive_no_outliers`** (positive then drop **`|rstudent| > 3`** from per-statistic `lm(stat ~ F_ST)`). Summary: **`fd_fdm_fst_100kb_report.md`**.

### CLI

- `Rscript fd_fdm_vs_fst_100kb_bins.R [stride2|greedy] [all|positive|positive_no_outliers]` — second arg defaults to **`positive`** if omitted.

### Outputs (repo root)

- Per mode suffix **`_all`**, **`_positive`**, **`_positive_no_outliers`**: `fd_*`, `fdm_*`, `df_*` TSVs and `*_vs_fst_100kb_bins_<suffix>.{png,pdf}`.

---

## Session: Fst distribution histogram at 10 kb

### Script

- **`genome-wide_Fst.Dxy_w_QTL_plots.R`** — after the existing 1 cM Fst histogram (`Fst_genomewide_distribution_w_QTL_windows.png`):
  - **`Fst_genomewide_distribution_10kb_w_QTL_max_windows.png`** — genome-wide 10 kb histogram + dark-blue **`geom_density()`** curve; colored vertical lines = **max** Fst in each QTL 95% CI (bp). Left spike at Fst ≈ 0 = windows with zero/near-zero differentiation (often few SNPs), stacked in the first `binwidth = 0.01` bin.
  - **`Fst_10kb_QTL_CI_mean_only.png`** — **mean** Fst across 10 kb windows in each QTL CI only (no genome-wide background).

1 cM plot unchanged (mean Fst in cM windows + mean Fst in QTL CI in cM).

---

## Session: chemistry QTL (not in manuscript)

Chemistry QTL is **out of scope** for the public repo and Dryad. Materials were removed from the active tree (git history retained) and backed up on dissert at `Costus_local_only_not_Dryad/` (and locally under gitignored `_unpublished_local/` if present).

---

## Convention: log all R warnings (batch / HPC scripts)

Long-running or HPC R scripts should record **every** warning to a text file, not only the first 50 echoed to Slurm stdout.

**Utility:** [`QTL/log_all_warnings.R`](QTL/log_all_warnings.R)

At the top of a script (after `library()` calls):

```r
for (p in c("log_all_warnings.R", file.path("QTL", "log_all_warnings.R"))) {
  if (file.exists(p)) { source(p, local = FALSE); break }
}
start_warning_log(default_basename = "my_script_warnings.txt")
```

- Sets `options(warn = 1)` so each warning prints when it occurs.
- `sink(..., type = "message", split = TRUE)` writes warnings to a file **and** stderr (Slurm log still gets them).
- `on.exit()` closes the log when the script finishes or errors.

**Override path:** `WARNING_LOG=/path/to/warnings.txt Rscript my_script.R`

**On Farm:** copy `log_all_warnings.R` alongside other scripts in `~/Costus_scripts/`.

---

## How to use this file

- **`AGENTS.md`** and **`CLAUDE.md`** in the repo root tell humans and agents to **read `NOTES.md`** for ongoing project context before making substantive edits.

---

## Session: publication repository cleanup (2026-07-08)

### Layout

Methods-aligned folders `01_`–`04_` under repo root; verbatim pre-cleanup scripts in `_originals/`. See root [`README.md`](README.md) and [`CHANGELOG.md`](CHANGELOG.md).

### Conventions added

- **`scripts/project_paths.R`** — `project_paths()` resolves `PROJECT_ROOT`, pixy/Dsuite/QTL paths from env or repo layout.
- **Farm slurm** — `run_*.slurm` wrappers with env vars; Hummingbird originals kept alongside.
- **Sample manifest** — `sample_manifest/combined_sample_manifest.tsv` only (524 unique `sample_id`). Build scripts / NCBI source tables live in gitignored `_unpublished_local/sample_manifest_build/`.
- **Avoid Finder on `/Volumes/dissert`** if the volume only appears in Disk Utility/Terminal (can crash Finder). Do publication file work in the git repo; use Terminal for dissert paths.

### QTL scantwo on Farm

Canonical: `03_QTL_mapping/rqtl/run_scantwo_per_trait.R` + `slurm/run_scantwo_all_traits.slurm`. Per-trait Hummingbird `.R`/`.sh` files kept as archives in `03_QTL_mapping/rqtl/`.

## Session: QTL line style + histogram percentile shading (2026-07-15)

**Opposite-of-expected QTL peaks → dashed** (solid = matches species expectation). Direction from raw genotype means (AA=vill, BB=allenii), not conditional fitqtl `a`:

| QTL | Direction | Style on plots |
|-----|-----------|----------------|
| SD chr2, GR chr3, Nreabs chr1, Chlreabs chr1 | matches expected | solid |
| SD chr6, tough chr3, tough chr8 | opposite | dashed |
| GR chr5 | opposite | dashed if shown; already excluded by CI filter |

Scripts: `genome-wide_Fst.Dxy_w_QTL_plots.R`, `combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R`.  
1 cM Fst/Dxy/fdM histograms: 5th/95th percentiles as **red tail shading** (not dashed red lines).  
`sup_Dxy_fdM_QTL_window_distributions.pdf` also written to `Figures&Tables/` (parent of `Final_figures/`).

---

## Session: simplified 1 cM Fst density (2026-08-27)

First-pass companion to `Fst_genomewide_distribution_w_QTL_windows.png` in `04_introgression_vs_QTL/genome-wide_Fst.Dxy_w_QTL_plots.R`. Same 1 cM window averages and QTL CI means; **no histogram**, **no 5th/95th percentile shading**, **no trait colors**.

- Black `geom_density()` of genome-wide 1 cM mean Fst
- Light-grey (`grey60`) solid `geom_vline()` at mean Fst within each QTL 95% CI
- White panel, no grid; x-axis label **Distribution of Fst (1cM windows)**
- Output: `Fst_genomewide_distribution_w_QTL_windows_simplified.png` (project root)
- After a full run, `Fst_1cM_window_cache.rds` (gitignored) lets the simplified `ggsave` chunk be re-run alone to tweak size without reloading pixy.

The original colored histogram figure is unchanged.

---

## Session: Hudson Fst rerun (2026-10-13)

Reviewer point (Bhatia et al. 2013): combine SNPs as a ratio of averages, and for two populations prefer the Hudson estimator over Weir–Cockerham.

- Farm job [`02_population_genetics_wild/pixy/run_pixy_array.slurm`](02_population_genetics_wild/pixy/run_pixy_array.slurm): `--fst_type hudson --fst_components`, submitted at **10 kb** and **500 kb**. Output dirs `pixy_hudson_10kb` and `pixy_hudson_500kb` on quobyte. `run_pixy_avonly.sh` (50 kb Weir–Cockerham) exits immediately; the full original-submission script is under `_originals/wc_fst_original_submission/`.
- Along-genome **500 kb line** is `avg_hudson_fst` / `avg_dxy` from the 500 kb pixy files. The R mean of 10 kb windows was removed.
- **1 cM** and QTL-interval F~ST~ stay in `genome-wide_Fst.Dxy_w_QTL_plots.R` because the bins are centimorgans. They use `sum(hudson_fst_num) / sum(hudson_fst_den)`. The simplified-density cache is `Fst_1cM_window_cache_hudson.rds`; do not reload `Fst_1cM_window_cache.rds`.
- 1 cM D~xy~ is still `mean(avg_dxy)`. Dsuite `f_dM` / `d_f` 500 kb lines are unchanged (pixy does not calculate them).
- Plotting scripts stop if the folder has `avg_wc_fst` and no `avg_hudson_fst`, so `output_251124_LASIref` cannot be plotted by the current code.
- Figures are regenerated only after both Farm tables are on disk (`PIXY_FOLDER`, `PIXY_FOLDER_500KB`).

---

## Session: 10 kb QTL-CI distributions (2026-10-28)

Additional versions of the 1 cM distribution figures, in `04_introgression_vs_QTL/genome-wide_Fst.Dxy_w_QTL_plots.R`. The 1 cM plots are unchanged.

- Genome-wide background: histogram plus dark-blue density of **10 kb** windows (chromosome 4 excluded). Red bands mark the genome-wide 5th and 95th percentiles.
- Each QTL is a **density of the 10 kb windows that overlap its 95% CI**, colored by trait. Antagonistic QTL (SD chr6, toughness chr3, toughness chr8) are dashed. Both toughness QTL share that color and dash.
- **Fst** keeps 10 kb windows with **≥10 SNPs**. **Dxy** keeps 10 kb windows with **≥1000 sites** (10% of the window). **f_dM** has no extra filter: each 10 kb value is the mean of positive Dsuite w50s25 windows whose midpoint falls in that bin.
- Axis limits are a display window only (`coord_cartesian`): Fst −0.25 to 0.65, Dxy 0 to 0.04, f_dM 0 to 0.55. Densities use the filtered windows.
- Files, project root: `Fst_Dxy_fdM_genomewide_distribution_10kb_w_QTL_CI_3panel.png`, `Fst_genomewide_distribution_10kb_w_QTL_CI.png`, `Dxy_fdM_genomewide_distribution_10kb_w_QTL_CI_2panel.png`.
- Fst panels only: a filled point on each within-CI curve at the Fst of the 10 kb pixy window that contains `peak_bp`, if that window passes the SNP filter.

### 0.1 cM distributions (2026-10-29)

Same figure design, in the same script. The 1 cM figures are unchanged.

- **Fst** 0.1 cM bins are `sum(hudson_fst_num) / sum(hudson_fst_den)` over 10 kb windows with ≥10 SNPs whose start falls in the bin.
- **Dxy** 0.1 cM bins are the mean of 10 kb `avg_dxy` from windows with ≥1000 sites.
- **f_dM** 0.1 cM bins are the mean of positive Dsuite windows in the bin (no 10 kb filter).
- A bin is inside a QTL CI when its interval overlaps the CI on the script's linkage map (`bp_to_cm_funs`: chromosome flips, then isotonic cM). The QTL table's raw cM columns are not used for this overlap. The Fst point is the bin that contains `peak_bp` on that map.
- Files: `Fst_Dxy_fdM_genomewide_distribution_0.1cM_w_QTL_CI_3panel.png`, `Fst_genomewide_distribution_0.1cM_w_QTL_CI.png`, `Dxy_fdM_genomewide_distribution_0.1cM_w_QTL_CI_2panel.png`.

The same floors now apply to the other F~ST~ and D~xy~ figures: along-genome 10 kb points, the 1 cM histograms (including the simplified F~ST~ density and the max-in-CI 10 kb histogram), the four-panel combo D~xy~ and F~ST~ panels, and `plot_pixy.R`. Pi is not filtered. The 500 kb line is still pixy's own 500 kb estimate, drawn when at least half of the 10 kb slots in that span pass the matching floor. f~dM~ is not filtered. The circular-shift outlier test is not a plot and still uses every finite 10 kb window.

## Session: supplementary Dsuite noNearHybrid 1 km (2026-11-08)

Exclude *C. allenii* / *C. villosissimus* within **1 km** of a PCA hybrid from a
supplementary Dinvestigate run. Original `Fstats_w50s25_..._251120` unchanged.

- **Kept parents:** geographic include list was 12 allenii + 12 vill; BAVLasi VCF
  intersection is **9 alle + 12 vill** (`19.497`, `19.508`, `19.553` absent from VCF).
  Map outlines and `alle_vill_include_IDs.txt` match that VCF set (21 parents).
- **Map:** `Figures&Tables/OG_sample_map_Dsuite_noNearHybrid_1km.pdf`
- **Farm package:** `~/Costus_scripts/dsuite/noNearHybrid_1km/` (also in repo
  `02_population_genetics_wild/dsuite/noNearHybrid_1km/`).
- **Staged VCFs:** dissert `Costus_archive/02_VCFs/BAVLasi_LASIref_filtered_scaf/` →
  `/quobyte/jrigrp/juliah/Costus_data_tmp/VCFs/BAVLasi_LASIref_filtered_scaf/`.
- **SETS:** `LASIout_SETS_noNearHybrid_1km.txt` (brac=20 alle=9 vill=12 Outgroup=2 xxx=20).
- **Submitted:** array job `40549436` (after fixing `WORKDIR` to use `SLURM_SUBMIT_DIR`;
  first attempt `40549387` failed because Slurm spool path broke SETS lookup).
- **Outputs:** `/quobyte/jrigrp/juliah/Costus_data_tmp/dsuite_noNearHybrid_1km_w50s25/`
- **Dsuite binary:** built on Farm into `~/bin/Dsuite` (conda `dsuite` env for libstdc++).
- **Subset pixy (Farm, same 21 parents):** jobs `40563512` (10 kb) and `40563513`
  (500 kb) **COMPLETED**; popfile `~/Costus_scripts/pixy_noNearHybrid_1km/avonly_populations_noNearHybrid_1km.txt`;
  Farm outputs `pixy_hudson_{10kb,500kb}_noNearHybrid_1km/`; local copies under `pixy/`.
  Fully subset combo 4-panel + circular-shift remade in
  `04_introgression_vs_QTL/noNearHybrid_1km/` (subset Dsuite + subset pixy).

## Session: Fst chromosome-panel point styling (2026-11-11)

Along-genome Fst panels in `04_introgression_vs_QTL/genome-wide_Fst.Dxy_w_QTL_plots.R`
(`Fst_10kb_LASIout_separate_chrs_w.500kb.avg_7QTL_narrowed.CIs`, 2-panel chr
figures, and `sup_genomewide_Fst_QTL_CIs`):

- Background 10 kb points: `size = 1`, `alpha = 0.35`, `grey50`
- 99th-percentile outliers: vector `geom_point` (not rast — PDF sparse-dot washout), `size = 1.5`
  (≈ grey rast `size = 1`), `alpha = 1`, `shape = 16`, `#D62728`
- Thin dashed `#D62728` `geom_hline` at the genome-wide 99th percentile
- All panels: `coord_cartesian(ylim = c(0, 1))` (display window only)
- 500 kb darkred line unchanged
- Fast replot while iterating (pixy + QTL only): `Rscript 04_introgression_vs_QTL/replot_fst_chr_panels.R`
- Also writes a **single-column** 8-panel stack (`*_1col.pdf/png`, 22 × 48 cm) so each chromosome track is wider than in the 2×4 grid.
- Layer order: grey 10 kb points → QTL CI/peak → dashed 99th line + red outliers + darkred 500 kb average (so outliers/avg sit on top of QTL).
- Subset panels share styling: `ylim = c(0, 0.75)`, y breaks `0/0.25/0.5/0.7`, no x ticks/labels,
  "Position (bp)" only on bottom panel, 2-panels at 22 × 12 cm.
- Three-panel **chr 1 / 2 / 3**: `pixy/Fst_chr1_chr2_chr3_3panel.{png,pdf}` (22 × 14 cm).
- Combo four-panel Fst: `coord_cartesian(ylim = c(0, NA))` (clip negative Hudson Fst in the display).
