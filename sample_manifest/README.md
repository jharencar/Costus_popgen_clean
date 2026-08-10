# Sample manifest

The sample metadata table is archived **only in the Dryad data deposit** (CC0), not in
this code repository / Zenodo software archive:

`sample_manifest/combined_sample_manifest.tsv` — one row per `sample_id` (524 unique IDs; wild collections + QTL mapping population).

| Column | Meaning |
|--------|---------|
| `dataset` | `wild_collection` or `qtl_mapping` |
| `sample_id` | Lab / analysis sample ID |
| `organism` | Species or hybrid designation |
| `collection_date` | Collection date (when applicable) |
| `geographic_location` | Locality string |
| `field_collection_id` | Field collection ID |
| `latitude`, `longitude` | Decimal degrees (WGS84); may be blank |
| `elevation_m` | Elevation in meters; may be blank |

Build scripts and NCBI source tables used to assemble this file are kept locally (not in the public tree). See the Dryad README for how to use the Dryad and Zenodo deposits together.
