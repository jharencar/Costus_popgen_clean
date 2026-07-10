# Hummingbird + dissert data archive

Archive target: **`/Volumes/dissert/Costus_archive/`** on macOS (`/run/media/julia/dissert/Costus_archive/` on Linux). Reorganized **2026-07-10**; Hummingbird gap-fill completed **2026-07-10**.

**Verified on dissert (2026-07-10):** `07_dsuite/from_hummingbird/` (2.4 GB), `06_pixy/from_hummingbird/` (215 MB), `05_rqtl/scantwo.perms/from_hummingbird/`, `08_ANGSD/abba_baba/from_hummingbird/` (167 MB), plus Mac copies: VCFs chr1–9, QTL posteriors (9232 files), lepmap3, genoprob.

Cross-checked against `humm_inventory_20260710.txt` and `dissert_inventory_20260710.txt`.

**Use Terminal for dissert paths** — avoid Finder on this volume (see `NOTES.md`).

---

## Workflow scripts (`sample_manifest/`)

| Script | Purpose |
|--------|---------|
| `inventory_dissert.sh` | Re-scan dissert; save to `dissert_inventory_YYYYMMDD.txt` |
| `reorganize_dissert_archive.sh` | Move `humm_backup` → `Costus_archive/`; copy QTL data from Mac |
| `rsync_from_hummingbird_missing.sh` | Rsync **gaps only** from Hummingbird |
| `cleanup_dissert_redundant.sh` | Remove broken symlinks / optional space recovery |
| `cleanup_unwanted_humm_rsync_vcfs.sh` | Delete OG / legacy VCFs pulled by mistake into `from_hummingbird/` |

```bash
# Already done (2026-07-10):
bash sample_manifest/reorganize_dissert_archive.sh --execute

# Next: pull gaps from Hummingbird (dry-run, then execute)
bash sample_manifest/rsync_from_hummingbird_missing.sh
bash sample_manifest/rsync_from_hummingbird_missing.sh --execute

# Optional: free ~20 GB on dissert
bash sample_manifest/cleanup_dissert_redundant.sh
bash sample_manifest/cleanup_dissert_redundant.sh --execute
```

---

## Costus_archive layout (post-reorganization)

```
/run/media/julia/dissert/Costus_archive/          (~180 GB; Linux mount)
  02_VCFs/                                29 GB
    BAVLasi_LASIref_filtered_scaf/        3.7 GB  (chr1–9 + .tbi)
    renamed_AVonly_pixy/                    4.2 GB  (chr1–9 + .tbi)
    other_filtered/                       2.2 GB
    legacy_uncompressed_AHV/               19 GB  ← optional delete
  03_aHMM_posteriors/                     10 GB
    QTL_grouped_run/                      9.7 GB  (from alle_vill_QTL)
    wild_AHV_LASIref2/                    391 MB  (from repo)
  04_lepmap3/                             2.9 GB  (from alle_vill_QTL)
  05_rqtl/                                4.3 GB  (genoprob + scantwo perms)
  06_pixy/                                 70 MB  (partial; pull full from Hummingbird)
  07_dsuite/                               27 MB  (Fstats only; pull full ~6.6 GB)
  08_ANGSD/                               133 GB  (abba redo BAMs + abba_baba outputs)
  09_raw_sequencing/                      symlinks → top-level FASTQ/BAM dirs
  10_scripts_legacy/                      old humm_backup scripts
  manifest/
```

Large raw data stays at `/Volumes/dissert/` top level; `09_raw_sequencing/` holds symlinks only.

| Top-level dir | Size | Role |
|---------------|------|------|
| `Costus_mapping_pop_raw_for_NCBI/` | 317 GB | Mapping pop FASTQs (canonical) |
| `All_wild_samps_raw_fastqs/` | 92 GB | Wild collection FASTQs |
| `OG_Combined_cleaned_bams/` | 44 GB | Wild cleaned BAMs |
| `humm_backup/` | ~119 MB | Residual lepmap order files only |

---

## FASTQ redundancy investigation (2026-07-10)

### Mapping population (`Costus_mapping_pop_raw_for_NCBI/`, 317 GB)

| Location | Files | Size | Verdict |
|----------|-------|------|---------|
| **Top level** (Illumina names, e.g. `21_228_S70_L002_R1_001.fastq.gz`) | **818** | **~312 GB** | **KEEP — canonical copy** (409 samples × R1/R2) |
| `fastqs_for_upload/` | 818 symlinks | ~0 | **DELETE — all 818 symlinks broken** (point to `/Volumes/Radishes/...`) |
| `_extras/duplicate_lanes/` | 2 | 212 MB | Optional delete — alternate lane for `22_033` (not in main set) |
| `_extras/not_mapping_pop/` | 4 | 1.9 GB | **KEEP** — samples excluded from mapping pop (`21_234`, `22_148`) |

**No duplicate 312 GB copy found.** The apparent “double storage” is broken symlinks in `fastqs_for_upload/`, not extra FASTQ data. Inflated file counts (~3284) include macOS `._` resource forks and broken symlinks.

### Wild collection (`All_wild_samps_raw_fastqs/`, 92 GB)

540 FASTQs with lab-style names (`18.005_amazonicusVar.gual_A_R1.fastq.gz`). **No basename overlap** with mapping pop. **KEEP.**

### Other locations

| Location | FASTQs | Notes |
|----------|--------|-------|
| `humm_backup/` | 0 | Analysis outputs only (post-reorg) |
| `_extras/` (top) | 0 real `.fastq.gz` | Metadata/scripts; `_extras/_extras/` has 1.9 GB clumped `.fq.gz` (unrelated legacy) |
| `greenhouse_costus/` | 14 | 6.9 GB — separate greenhouse experiment |

### Space recovery summary

| Action | Frees | Risk |
|--------|-------|------|
| Delete `fastqs_for_upload/` | ~0 (broken symlinks) | None — cleanup only |
| Delete `02_VCFs/legacy_uncompressed_AHV/` | **~19 GB** | Low after verifying gz VCFs in `02_VCFs/` |
| Delete `_extras/_extras/` | **~1.9 GB** | Low — unrelated clumped files |
| Delete `_extras/duplicate_lanes/` | **~212 MB** | Only if alternate `22_033` lane not needed |
| Delete mapping FASTQs | 312 GB | **Do not** — only copy |

---

## Still needed from Hummingbird (~7 GB)

**Status: gap-fill completed 2026-07-10** (see verified sizes above). Re-run `rsync_from_hummingbird_missing.sh` only if files are missing.

After reorganization, rsync **only these gaps** (script excludes all `*.vcf*`, `*OG*`, `*unfiltered*`, `trash/`, `old_*`):

| Source (Hummingbird) | Destination (dissert) | Size |
|----------------------|----------------------|------|
| `Julia_popgen/Dsuite/` | `Costus_archive/07_dsuite/from_hummingbird/` | ~6.6 GB → smaller without VCFs |
| `Julia_popgen/pixy/` | `Costus_archive/06_pixy/from_hummingbird/` | ~62 MB |
| `Julia_popgen/rqtl/scantwo.perms.done/` | `Costus_archive/05_rqtl/scantwo.perms/from_hummingbird/` | `.rds` only |
| `Julia_popgen/ANGSD/abba_baba/` | `Costus_archive/08_ANGSD/abba_baba/from_hummingbird/` | merge |

**If rsync is pulling unwanted VCFs:** press `Ctrl+C` in the terminal, pull the updated script, remove partial junk, restart:

```bash
# 1. Stop running rsync: Ctrl+C

# 2. Preview/delete OG / any VCFs already copied into from_hummingbird/
bash sample_manifest/cleanup_unwanted_humm_rsync_vcfs.sh
bash sample_manifest/cleanup_unwanted_humm_rsync_vcfs.sh --execute

# 3. Dry-run then re-transfer (skips files already present; won't re-pull excluded names)
bash sample_manifest/rsync_from_hummingbird_missing.sh
bash sample_manifest/rsync_from_hummingbird_missing.sh --execute
```

Canonical publication VCFs stay in `Costus_archive/02_VCFs/` only — **never** rsync `Julia_popgen/VCFs/` from Hummingbird.

**Skip entirely:**

- Filtered VCFs (now in `Costus_archive/02_VCFs/`)
- Davis / mapping FASTQs (317 GB on dissert)
- QTL posteriors, lepmap3, genoprob (copied from `alle_vill_QTL`)
- Wild aHMM LASIref2 (copied from repo)
- Whole `Julia_popgen` (2 TB), `QTL_mapping` scratch (2.1 TB)

### Manual rsync (if not using script)

```bash
DEST="/Volumes/dissert/Costus_archive"
HUMM="jharenca@hb.ucsc.edu"
X='--exclude=trash/ --exclude=old_* --exclude=*OG* --exclude=*unfiltered* --exclude=*.vcf --exclude=*.vcf.gz --exclude=*.vcf.gz.tbi'

# Dry-run each; remove -n to transfer
rsync -avP --partial -n $X "${HUMM}:/hb/groups/kay_lab/Julia_popgen/Dsuite/" \
  "${DEST}/07_dsuite/from_hummingbird/"

rsync -avP --partial -n $X "${HUMM}:/hb/groups/kay_lab/Julia_popgen/pixy/" \
  "${DEST}/06_pixy/from_hummingbird/"

rsync -avP --partial -n --include='*.rds' --include='*/' --exclude='*' \
  "${HUMM}:/hb/groups/kay_lab/Julia_popgen/rqtl/scantwo.perms.done/" \
  "${DEST}/05_rqtl/scantwo.perms/from_hummingbird/"

rsync -avP --partial -n $X "${HUMM}:/hb/groups/kay_lab/Julia_popgen/ANGSD/abba_baba/" \
  "${DEST}/08_ANGSD/abba_baba/from_hummingbird/"
```

Log transfers: `tee "${DEST}/manifest/rsync_logs/dsuite_$(date +%Y%m%d).log"`

**rsync exit code 23:** partial transfer. On exFAT dissert this is often harmless (can't set Unix permissions). Check the log:

```bash
grep -iE 'error|denied|vanished|No space' /Volumes/dissert/Costus_archive/manifest/rsync_logs/*.log | tail -20
```

If only `failed to set permissions` / `chown` / `utime` messages, data likely transferred OK — re-run `--execute` to fill gaps. Updated script uses `--no-perms --no-owner --no-group` to avoid this.

---

## Hummingbird inventory summary (2026-07-10)

| Path | Size | Pull? |
|------|------|-------|
| `/hb/groups/kay_lab/Julia_popgen` | 2.0 TB | Subdirs only (see gaps above) |
| `/hb/scratch/jharenca/QTL_mapping` | 2.1 TB | Skip (posteriors on Mac) |
| `/hb/groups/kay_lab/Davis_data/plates_8.9.10/Data` | 194 GB | Skip (on dissert) |
| `/hb/home/jharenca/pop_gen` | 21 GB | Skip unless wild pipeline needed |
| `Julia_popgen/Dsuite/` | 6.6 GB | **Yes** |
| `Julia_popgen/pixy/` | 62 MB | **Yes** |
| `Julia_popgen/VCFs/` | 412 GB | Skip (on dissert) |
| `Julia_popgen/rqtl/scantwo.perms.done/` | part of 2.2 GB | Diff / pull |
| Concatenated QTL posteriors on scratch | 0 files | Skip (on Mac in `alle_vill_QTL`) |

---

## After Hummingbird pull

1. Re-run `bash sample_manifest/inventory_dissert.sh > sample_manifest/dissert_inventory_$(date +%Y%m%d).txt`
2. `shasum -a 256` on key VCFs and genoprob for Dryad
3. Update root `README.md` external data paths
4. Do **not** commit large files to git

---

*Last updated 2026-07-10 after dissert reorganization and FASTQ redundancy audit.*
