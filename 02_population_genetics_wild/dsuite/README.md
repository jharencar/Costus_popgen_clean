# Dsuite analyses

| Script | Purpose | HPC |
|--------|---------|-----|
| `run_Dinvestigate_fdm_array.slurm` | Genome-wide f_dM (w50s25) per chromosome | Farm-ready |
| `run_Dtrios_capseq.slurm` | Ecuador capture ABBA-BABA (Dtrios) | Farm-ready |
| `noNearHybrid_1km/` | Supplementary f_dM excluding samples within 1 km of hybrids | Farm (`~/Costus_scripts/dsuite/`) |
| `LASIout_Dinvestigate_array.sh` | Original Hummingbird Dinvestigate job | `_originals/` |
| `Dtrios_capseq_LASIout.sh` | Original Hummingbird Dtrios job | `_originals/` |

**Required inputs (not in repo):** filtered VCFs, `LASIout_SETS.txt`, `test_trios.txt` (f_dM); capture VCF + `capseq_LASIout_SETS.txt` + `tree.nwk` (Dtrios).

**Plotting:** `04_introgression_vs_QTL/combo_Dsuite_pixy_4panel_fdm_df_dxy_fst.R`

### Supplementary: noNearHybrid_1km

See [`noNearHybrid_1km/README.md`](noNearHybrid_1km/README.md). Farm copy:
`~/Costus_scripts/dsuite/noNearHybrid_1km/`. Needs BAVLasi VCFs staged under
`/quobyte/jrigrp/juliah/Costus_data_tmp/VCFs/BAVLasi_LASIref_filtered_scaf/`
(from dissert `Costus_archive/02_VCFs/BAVLasi_LASIref_filtered_scaf/`).
