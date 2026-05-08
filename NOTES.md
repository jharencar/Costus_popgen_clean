# Project notes (Costus_popgen_clean)

This file captures decisions and work from development sessions so future tooling and agents stay aligned with local conventions. **If you are an agent, read this file when working in this repository.**

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

## How to use this file

- **`AGENTS.md`** and **`CLAUDE.md`** in the repo root tell humans and agents to **read `NOTES.md`** for ongoing project context before making substantive edits.
