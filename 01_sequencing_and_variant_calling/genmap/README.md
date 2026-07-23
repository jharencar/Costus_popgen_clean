# Mappability (genmap)

Manuscript filter: **genmap** with **k = 150**, merge distance **100 bp**, exclude regions with mappability **&lt; 0.7** (*C. lasius* reference GCA_027563935.1).

| File | Purpose |
|------|---------|
| `run_genmap_mask.sh` | Index genome → genmap bedGraph → `non_callable_sites.bed` |

Downstream VCF filters exclude the mask with `bcftools view -T ^non_callable_sites.bed` (see `../vcf_filtering/`). On Hummingbird the production BED was stored as e.g. `/hb/groups/kay_lab/genomes/genmap/ALLE_non_callable_sites.bed` (or LASIref equivalent).
