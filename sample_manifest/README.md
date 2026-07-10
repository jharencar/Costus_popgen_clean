# Sample manifest

| File | Description |
|------|-------------|
| `combined_sample_manifest.tsv` | **524 rows** — one row per `sample_id` (wild + QTL mapping) |
| `combined_sample_manifest_before_dedup.tsv` | Pre-deduplication backup (590 rows; 66 redundant duplicates removed) |
| `wild_sample_metadata.tsv` | Source: `/Volumes/dissert/_extras/` (NCBI prep) |
| `qtl_mapping_sample_metadata.tsv` | Source: `/Volumes/dissert/Costus_mapping_pop_raw_for_NCBI/` |
| `Kay_lab_database_full.csv` | Lab ID crosswalk |
| `build_combined_manifest.py` | Regenerate combined manifest from NCBI source TSVs |

**Columns:** `dataset`, `sample_id`, `organism`, `collection_date`, `geographic_location`, `field_collection_id`, `latitude`, `longitude`, `elevation_m`

**Deduplication:** Redundant rows sharing the same `sample_id` were removed (kept first occurrence). All removed duplicates were identical across columns.

## Data archive (dissert + Hummingbird)

External analysis data lives on **`/Volumes/dissert/Costus_archive/`** (reorganized 2026-07-10). On Linux set `COSTUS_DISS=/run/media/julia/dissert`. See:

| File | Purpose |
|------|---------|
| [`HUMMINGBIRD_DATA_ARCHIVE.md`](HUMMINGBIRD_DATA_ARCHIVE.md) | What's on dissert, FASTQ redundancy audit, Hummingbird gap rsync |
| [`archive_paths.env.example`](archive_paths.env.example) | `COSTUS_ARCHIVE` env vars for publication R scripts |
| `dissert_inventory_20260710.txt` | Dissert scan (post-reorg) |
| `humm_inventory_20260710.txt` | Hummingbird cluster scan |
| `inventory_dissert.sh` / `inventory_hummingbird.sh` | Re-run inventories |
| `reorganize_dissert_archive.sh` | Move scattered backup → `Costus_archive/` |
| `rsync_from_hummingbird_missing.sh` | Pull ~7 GB still missing from Hummingbird |
| `cleanup_dissert_redundant.sh` | Remove broken FASTQ symlinks; optional ~21 GB recovery |

```bash
# Pull remaining Hummingbird gaps (after reorganize)
bash sample_manifest/rsync_from_hummingbird_missing.sh --execute
```
