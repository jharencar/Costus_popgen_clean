# Sequencing and variant calling

Methods section: *DNA extraction, sequencing, and initial data cleaning*.

| Subfolder | Contents |
|-----------|----------|
| `preprocessing/` | Trimmomatic, BWA, dedupe, BAM filters (`seq_analysis`) |
| `vcf_filtering/` | bcftools VCF generation and LASIref filters (`humm_backup`) |
| `genmap/` | *(add genmap mask scripts when located on Hummingbird)* |

HPC paths in scripts are from UCSC Hummingbird / UC Davis Farm; retarget via environment variables before running.
