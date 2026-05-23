# Data dictionary: `genetic_architecture_summary.csv`

One row per mapped trait (Chlreabs, GR, Nreabs, SD, tough). Combines **Tier A** (GRM + REML narrow-sense *h*²) and **Tier B** (fixed-effect variance partition from `fitqtl` drop-one tables).

**Build order:** `genomic_heritability_QTL_traits.R` → `genetic_variance_partition.R` → `build_genetic_architecture_summary.R`

**Inputs:**

| Input | Role |
|-------|------|
| `QTL/*_QTL_model_output_251203` | Final Haley–Knott `fitqtl` summaries (Dec 2025 models) |
| `aHMM_..._cleaned_most_dups_genoprob.rds` | Genotype probabilities for GRM (see `NOTES.md`) |

Machine-readable column metadata: [`genetic_architecture_summary_column_metadata.csv`](genetic_architecture_summary_column_metadata.csv)

---

## Column reference

| Column | What it is | How it is calculated |
|--------|------------|----------------------|
| **trait** | Trait shorthand | Label only |
| **n_obs** | Sample size in QTL model | `Number of observations` from `fitqtl` output |
| **pct_QTL** | % phenotypic variance from QTL | Sum of drop-one **%var** for QTL terms (`chr@pos`) |
| **pct_hybrid_index** | % phenotypic variance from genome-wide ancestry | Drop-one **%var** for `hybrid_index` (0 for SD primary) |
| **pct_F21** | % phenotypic variance from F2 line | Sum of drop-one **%var** for `F21_*` dummies |
| **pct_environment** | % phenotypic variance from batch covariates | Sum of drop-one **%var** for `coh*` or `s.coh*` |
| **v_QTL** | QTL variance on trait scale (approx.) | Σ (`pctvar`/100 × `SS_total`) for QTL terms |
| **v_genetic_fixed** | Summed “genetic” fixed-effect variance | `v_QTL` + hybrid_index + F21 increments |
| **h2_fixed_approx** | Supplementary ratio only | `v_genetic_fixed` ÷ (`SS_total` / (`n_obs` − 1)) — **not** narrow-sense *h*² |
| **n** | REML sample size | Same individuals as `n_obs` (complete cases) |
| **sigma_g** | **Vg** (additive genetic variance) | REML variance of genomic random effect *g* |
| **sigma_e** | Residual variance **Ve** | REML residual variance |
| **vp_cond** | Conditional phenotypic variance | `sigma_g` + `sigma_e` |
| **h2** | **Primary narrow-sense genomic *h*²** | `sigma_g` / (`sigma_g` + `sigma_e`) |
| **n_markers_grm** | Markers in GRM | After MAF ≥ 0.05 and ≤ 20% missing |
| **family_icc** | F2-line ICC check | `lmer(..., (1|family))` — not primary *h*² |
| **h2_dosage** | Sensitivity *h*² | Same as **h2** (dosage GRM is primary) |
| **h2_hard** | Sensitivity *h*² | REML with hard-call GRM |
| **frac_genetic_QTL** | QTL share of genetic variance | `pct_QTL` / (`h2` × 100) |
| **frac_genetic_QTL_note** | Formula reminder | Text constant per row |

---

## Primary vs supplementary metrics

| Use in text | Columns |
|-------------|---------|
| **Primary heritability** | `h2`, `sigma_g`, `vp_cond` |
| **QTL phenotypic contribution** | `pct_QTL` |
| **QTL as % of genetic variance** | `frac_genetic_QTL` |
| **Supplementary / do not call *h*²** | `h2_fixed_approx`, `pct_*` sums, `v_genetic_fixed` |

---

## Per-trait covariate structure (REML and fitqtl)

| Trait | Phenotype | Fixed covariates in model |
|-------|-----------|---------------------------|
| SD | `log_days_to_radicle` | `s.coh4`, `s.coh7`, `s.coh8`, `s.coh9`, `s.coh14`, `s.coh15` + QTL |
| GR, tough, Nreabs, Chlreabs | See `Costus_QTL.qmd` | `hybrid_index` + selected `F21_*` + selected `coh*` + QTL |

QTL loci are **not** included as fixed effects in the REML model (only in `fitqtl`).

---

## Related files

- `genetic_variance_partition_summary.csv` — full partition including `ss_total`, `vp_total`, `full_model_pctvar`
- `genomic_heritability_QTL_traits.csv` — REML details including `h2_se`
- `SD_hybrid_index_secondary_comparison.csv` — SD sensitivity with `hybrid_index` added (not in main summary table)
