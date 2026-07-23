# Sequencing and variant calling

Methods section: *DNA extraction, sequencing, and initial data cleaning*.

| Subfolder | Contents |
|-----------|----------|
| `preprocessing/` | Trimmomatic, BWA, dedupe, BAM filters |
| `vcf_filtering/` | bcftools VCF generation and LASIref filters |
| `genmap/` | genmap mappability index/map → `non_callable_sites.bed` (k=150, merge 100 bp, mask &lt;0.7) |

HPC paths in scripts are from UCSC Hummingbird / UC Davis Farm; retarget via environment variables before running (`REF_FASTA`, `GENMAP_OUTDIR`, etc.).
