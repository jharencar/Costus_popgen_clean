# Allele / hybrid_index coding investigation (2026-07-13)

## Bottom line

**Not yet 100% locked to species.** What *is* locked:

| Quantity | Relationship | Evidence |
|----------|--------------|----------|
| `hybrid_index` | ≈ genome-wide **ref** ancestry = `mean(ref_homo + 0.5×het)` | *r* ≈ **0.936** across 385 F3s with QTL-run posteriors |
| R/qtl **AA** dosage | tracks `hybrid_index` | *r*(HI, doseA) ≈ **0.922** |
| R/qtl **BB** dosage | opposite of HI | *r*(HI, doseB) ≈ **−0.922** |

So: **HI high = more AA = more aHMM reference (parent1 / ancestry state 0)**.  
Species assignment **requires** the QTL `.pop` used for `new_vcf2ahmm.py`.

## QTL aHMM run that built the mapping genotypes

Pipeline that matches the genoprob stem  
`aHMM_251111_minGTP3_minDif0.7_dist5k_minDP0.05_pulse0.0001`:

[`03_QTL_mapping/aHMM_F3/512x64_LASIref_ahmm_pipeline_F3s_F2parents_F1parent.sh`](../../03_QTL_mapping/aHMM_F3/512x64_LASIref_ahmm_pipeline_F3s_F2parents_F1parent.sh)

```text
--pop ../F3s0.05covCO_allF2s_F1parent.pop
--minGTP 3 --minDif 0.7 --dist 5000 --minDP 0.5
-p 0 100000 0.5 -p 1 0.0001 0.5
```

Working dir on Farm:  
`/hb/scratch/jharenca/QTL_mapping/aHMM/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO`

**MISSING LOCALLY and on `/Volumes/dissert` (searched 2026-07-13):**  
`F3s0.05covCO_allF2s_F1parent.pop` / `.ploidy` were **not** archived with the posteriors.

- QTL posteriors on dissert: `Costus_archive/03_aHMM_posteriors/QTL_grouped_run/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO/` (chr `.panel` + `.posterior` only).
- Wild posteriors: `Costus_archive/03_aHMM_posteriors/wild_AHV_LASIref2/AHV_LASIref2_260210/` (also no `.pop`).
- **Zero** `*.pop` / `*.ploidy` anywhere under `/Volumes/dissert`.
- Legacy scripts on dissert still use `--pop ../F3s0.05covCO_allF2s_F1parent.pop`.

Next places on Farm (from hummingbird inventory / scratch path):
- `/hb/groups/kay_lab/Julia_popgen/aHMM/old_.pops`
- `/hb/scratch/jharenca/QTL_mapping/aHMM/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO/`

## How parent1 becomes AA

From `vcf2aHMM` / `new_vcf2ahmm`:

1. Panel columns after chrom/pos: **parent1 counts, then parent2 counts**, then admixed.
2. Ancestry_HMM posterior columns (as used in this project): `ref_homo`, `het`, `alt_homo`  
   → `ref_homo` = homozygous ancestry **0** = **parent1**.
3. Hard-call encoding (lepmap / rQTL input):  
   `ref_homo > 0.9 → AA`, `alt_homo > 0.9 → BB`.

Therefore: **AA = parent1**, **BB = parent2**. Period.  
Only the **labels inside the `.pop` file** assign parent1/2 to vill vs allenii.

## Local / historical `.pop` files (parent order)

| File | parent1 | parent2 | Notes |
|------|---------|---------|-------|
| `alle_vill_QTL/aHMM/QTL_aHMM_parent_subset.pop` | VILL | ALLE | Older QTL subset run (different params) |
| `alle_vill_QTL/aHMM/aHMM_pop_subALLE_parents_oldstyle_240409.pop` | VILL | ALLE | Same parents; admixed list differs |
| `Population_genetics_vill_alle/aHMM/aHMM_pop_parent_subset_old_style.pop` | VILL | ALLE | Wild-oriented |
| `generate_pop_files.sh` / `untitled.sh` | VILL | ALLE | Always writes vill first |

**Convention in every recoverable file: parent1 = villosissimus, parent2 = allenii.**  
That *suggests* AA = vill for the QTL run **if** `F3s0.05covCO_allF2s_F1parent.pop` follows the same convention — but that file must be confirmed on dissert/Farm.

Wild AHV LASIref2 (`AHV_LASIref2_260210`) should use the **same parent1/parent2 species order** for strip plots and QTL genotype labels to mean the same thing. Confirm its `.pop` when dissert is mounted (`COSTUS_ARCHIVE/.../wild_AHV_LASIref2/` and QTL grouped_run).

## Conflicting `hybrid_index` comments (do not trust blindly)

| Location | Comment |
|----------|---------|
| `03_QTL_mapping/rqtl/Costus_QTL.qmd` L227 | `# C.villosissimus ancestry proportion` |
| `alle_vill_QTL/alle_vill_qtl_mapping_fresh.R` | villos |
| `alle_vill_QTL/qtl_mapping.R` L320 | `# C.allenii ancestry proportion` |

Comments disagree. Empirically HI = **ref/parent1** ancestry, not “whatever the comment says.”

## Phenotype cross-checks (user priors)

User ground truth used here:

- Allenii = tougher leaves  
- Vill = higher Chl & N reabsorption  
- Vill = higher flavonoids (independent chem)  
- Seed dormancy & growth rate: user reports **opposing QTL** (one peak each direction), not a single species direction

Data (genotype class means / correlations with HI):

| Trait | Direction of AA / high-HI | If AA=vill (pop convention) | Matches user? |
|-------|---------------------------|-----------------------------|---------------|
| Chloro reabs | AA higher; cor(HI,~0.20) | Vill high | **Yes** |
| N reabs | AA higher; cor(HI,~0.27) | Vill high | **Yes** |
| Leaf toughness | AA higher; cor(HI,~0.10) | Vill high | **No** (user: allenii tougher) |
| Days to radicle | AA slightly higher (genomewide weak); SD QTLs bidirectional | Ambiguous | Partial |
| Chem flav QTL cluster (chr1) | BB higher (prior direction analysis) | Allenii high | **No** (user: vill higher flav) |

**Reabsorption strongly supports AA = vill if parents were labeled as usual.**  
Toughness and flav independently argue the opposite *or* those QTL directions ≠ species mean differences.

→ Phenotype priors alone are **not** sufficient to flip the coding; they currently conflict with each other under a single AA↔species map.  
→ Best resolver = open `F3s0.05covCO_allF2s_F1parent.pop` (+ wild AHV `.pop`).

## When dissert is remounted — checklist

```bash
# Default archive root
export COSTUS_ARCHIVE=/Volumes/dissert/Costus_archive

# QTL pop used for mapping genotypes (highest priority)
find "$COSTUS_ARCHIVE" /Volumes/dissert -name 'F3s0.05covCO_allF2s_F1parent.pop' 2>/dev/null
head -2 /path/to/F3s0.05covCO_allF2s_F1parent.pop

# Wild AHV pop (compare parent1/parent2 order)
find "$COSTUS_ARCHIVE" -path '*AHV_LASIref2*' -name '*.pop' 2>/dev/null
# Also Farm if available:
# /hb/scratch/jharenca/QTL_mapping/aHMM/grouped_run_minDif0.7_LASIref_w.F2s_0.05covCO/F3s0.05covCO_allF2s_F1parent.pop
```

Expected if consistent with local convention:

```text
parent1 ..._VILL,...
parent2 ..._ALLE,...
```

Then: **AA = villosissimus, BB = allenii, HI = vill ancestry** — and reabs supports that; toughness/flav become biological surprises (or need separate checking), not evidence to flip labels.

If parent1 is ALLE in that file, **all interpretations above flip**.

## Interim action for chemistry directions

Until the QTL `.pop` is verified:

1. Prefer reporting directions as **AA vs BB** / **high-HI vs low-HI**, not species names.  
2. Or report provisional species labels with explicit caveat: “assuming parent1=VILL as in all local `.pop` files.”  
3. Fix/flag conflicting HI comments after pop verification (do not “fix” species wording in docs until then).
