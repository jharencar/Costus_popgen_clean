# Project notes (Costus_popgen_clean)

This file captures decisions and work from development sessions so future tooling and agents stay aligned with local conventions. **If you are an agent, read this file when working in this repository.**

## Methods summary (manuscript)

Full prose in `popgen_manu_methods.txt`. Below: data types, sample sets, and filters agents need when writing or modifying analysis/plotting scripts.

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

- **Mappability:** genmap, **k = 15**, merge distance **100 bp**, keep regions **>= 0.7**.
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

#### PIXY — π and D~xy~ (genome-wide differentiation)

- **All-sites VCF** (invariant sites included): bcftools mpileup + core filters above + **<80% missing per species**.
- Windows: **10 kb non-overlapping**; PIXY → average π within species, D~xy~ between species.
- **Repo:** `genome-wide_Fst.Dxy_w_QTL_plots.R`, pixy two-panel chr plots; combo figures stack pixy with Dsuite (shared cumulative x).

#### Dsuite — *f*~dM~ (introgression sliding windows)

- **50 variant** windows, **25 variant** step (`w50s25` in paths/filenames); **Dinvestigate**; **biallelic SNPs** only; **>3 bp** from indels; AF **0.05–0.95**; same qual/coverage/mappability as core pipeline.
- **Samples:** **21** *bracteatus*, **21** *allenii*, **18** *villosissimus*, **2** *lasius* + *lasius* reference; **no PCA hybrids**.
- **Repo only:** highlight extreme |*f*~dM~| at **99.9th percentile** among negative windows (qmd/combo convention); manuscript F~ST~ figures use **99th percentile** outlier highlighting (different stat/plot).

#### F~ST~ vs introgression / QTL (repo extensions)

- Manuscript QTL comparison: **10 kb** F~ST~ windows + **1 cM** windowed averages; **chr4 excluded** from cM mapping (poor Lep-MAP3 join).
- **Repo:** `fd_fdm_vs_fst_100kb_bins.R` — **100 kb** bins, pixy SNP-weighted F~ST~ from 10 kb tiles + Dsuite *f*~d~/*f*~dM~/*d*~f~ with **stride-2** (or greedy) thinning on w50s25 windows (not in manuscript).

### QTL / mapping population (shorter)

#### aHMM (F3)

- Same parameters as wild hybrids; F3 + pedigree used for Lep-MAP3 input.

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
| `QTL/aHMM_prelim_UV_chemistry_est.rf.rds` | Preliminary UV chemistry QTL cross (`est.rf` only; in repo) |
| `QTL/prelim_UV_chemistry_scanone_by_trait.rds`, `*_covariates_by_trait.rds`, `*_perm_cutoffs_*.rds` | Prelim chemistry scanone outputs (in repo) |
| `QTL/circular_shift_perm_QTL_CI_Fst_fdM.rds` | Circular-shift null for QTL CI vs Fst/fdM overlap (in repo) |

**Agent defaults:** For new QTL or genomic-heritability scripts, read **`most_dups_genoprob.rds`** unless explicitly reproducing an old scan that used **`least_dups`**. Do not assume `.rds` files are in `Costus_popgen_clean`—set a configurable path at the top of scripts.

---

## Session: genomic heritability and genetic variance (QTL traits)

### Scripts (repo `QTL/`)

- **`fitqtl_summary_utils.R`** — shared parser for `*_QTL_model_output_251203` (used by variance plot, partition, heritability audit).
- **`genomic_heritability_QTL_traits.R`** — GRM + `sommer` REML; env `GENOPROB_RDS`; optional project **`R_libs/`** for dependencies.
- **`genetic_variance_partition.R`** — SS-based supplementary partition from fitqtl drop-one tables.
- **`build_genetic_architecture_summary.R`** — merges REML + partition → **`genetic_architecture_summary.csv`**.

### Conventions

- **Primary** narrow-sense \( \hat{h}^2 = \hat{\sigma}^2_g / (\hat{\sigma}^2_g + \hat{\sigma}^2_e) \) conditional on the **same fixed covariates** as each trait’s final `fitqtl` model (not QTL genotypes in the REML fixed effect).
- **Supplementary** `h2_fixed_approx` from summed Type-III `%var` (QTL + hybrid_index + F21) — can exceed 1; do not treat as \( h^2 \).
- **`frac_genetic_QTL`** = `pct_QTL` / (`h2` × 100).
- Sample audit: `QTL_trait_sample_audit.csv` (complete-case `n` must match fitqtl `n` for all five traits).

### Manuscript vs repo — quick flags

| Topic | Manuscript | Repo / NOTES |
|-------|------------|--------------|
| Wild *villosissimus* after coverage filter | 20 (38−10−3 allenii only listed) | ABBA-BABA / *f*~dM~ use **18** (hybrid exclusion + analysis set) |
| F3 count | 390 collected; **389** in Lep-MAP | **385** in final QTL cross |
| aHMM posterior filter | Not stated | **>0.9** on ref/alt/het (plotting) |
| F~ST~ outlier highlight | **99th** percentile (figures) | *f*~dM~ plots often **99.9th** percentile |
| F~ST~ vs *f*~dM~ correlation | — | **100 kb** bins + stride-2 (`fd_fdm_vs_fst_100kb_bins.R`) |
| PIXY / F~ST~ window size | **10 kb** non-overlapping | Same in repo histogram/QTL scripts |

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

## How to use this file

- **`AGENTS.md`** and **`CLAUDE.md`** in the repo root tell humans and agents to **read `NOTES.md`** for ongoing project context before making substantive edits.
