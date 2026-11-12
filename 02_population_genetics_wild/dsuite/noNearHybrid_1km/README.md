# Supplementary Dsuite: no individuals within 1 km of hybrids

Supplementary `f_dM` / `d_f` run. The original
`Fstats_w50s25_Mappability.Depth.Transl.Filtered.LASIref_251120` analysis is
**not** modified.

## Sample filter

From the Panama Dsuite/BAV sample set, drop any *C. allenii* or
*C. villosissimus* within **1 km** of a PCA hybrid. Hybrids remain excluded
(as in the original). All 21 *C. bracteatus* and both *C. lasius* stay.

| Set | allenii | vill | Total |
|-----|---------|------|-------|
| Original BAV Dsuite parents | 28 | 18 | 46 |
| Excluded by 1 km | 16 | 6 | 22 |
| **This run** | **12** | **12** | **24** |

IDs: `alle_vill_include_IDs.txt`. Distance table:
`dsuite_no_near_hybrid_samples.csv`.

## Farm paths

Scripts live in:

```text
~/Costus_scripts/dsuite/noNearHybrid_1km/
```

Expected VCFs (stage once from dissert):

```text
/quobyte/jrigrp/juliah/Costus_data_tmp/VCFs/BAVLasi_LASIref_filtered_scaf/
  BAVLasi_LASIref_filtered_scaf_{1..9}.vcf.gz
  BAVLasi_LASIref_filtered_scaf_{1..9}.vcf.gz.tbi
```

Source on the dissert archive:

```text
/Volumes/dissert/Costus_archive/02_VCFs/BAVLasi_LASIref_filtered_scaf/
```

Outputs (separate from the original):

```text
/quobyte/jrigrp/juliah/Costus_data_tmp/dsuite_noNearHybrid_1km_w50s25/
```

## One-time setup on Farm

```bash
# 1) Install Dsuite from GitHub (not on bioconda)
cd ~/Costus_scripts/dsuite/noNearHybrid_1km
bash install_dsuite_farm.sh
export PATH="$HOME/bin:$PATH"

# 2) Stage BAVLasi VCFs (from a machine with dissert mounted):
# rsync -avP /Volumes/dissert/Costus_archive/02_VCFs/BAVLasi_LASIref_filtered_scaf/ \
#   farm:/quobyte/jrigrp/juliah/Costus_data_tmp/VCFs/BAVLasi_LASIref_filtered_scaf/

# 3) Build SETS from the VCF sample list
bash build_LASIout_SETS_noNearHybrid_1km.sh \
  /quobyte/jrigrp/juliah/Costus_data_tmp/VCFs/BAVLasi_LASIref_filtered_scaf/BAVLasi_LASIref_filtered_scaf_1.vcf.gz

# 4) Submit
mkdir -p ~/Costus_scripts/slurm_log
sbatch --export=ALL,PATH="$HOME/bin:$PATH" run_Dinvestigate_noNearHybrid_1km.slurm
```

Trio file `test_trios.txt` matches the original `brac_alle_vill` naming.
