#!/bin/bash
#SBATCH --partition 128x24
#SBATCH --job-name=preANGSD_q20_array
#SBATCH --output=preANGSD_q20_array.%A_%a.out
#SBATCH -e preANGSD_q20_array.%A_%a.err
#SBATCH --mail-type=BEGAN,END,FAIL
#SBATCH --mail-user=jharenca@ucsc.edu
#SBATCH --mem=6G 
#SBATCH --cpus-per-task=12
#SBATCH --array=1-11 # set to number of scaffold.list files

# (consider adding #SBATCH --time=24:00:00 back when I have an idea of time to run the script)

# Code to add fields to initial vcf (Nicolas' snp filtering 3 step code - steps 2 and 3); code originally written by Tyler Linderoth
# Julia Harenčár 230719

# load bcftools
module load bcftools

# step 2: adding more fields to initial vcfs from mpileup 
for i in AHVF1_unfiltered*.vcf.gz; do 
       echo "bcftools +fill-tags $i -Oz -o ${i/.vcf.gz/_2.vcf.gz} -- -t 'AN,AC,MAF,AC_Hom,AC_Het,ExcHet'"
done > AHVF1_add_tags_job.list

for i in AHVF1_unfiltered*.vcf.gz; do 
       bcftools +fill-tags $i -Oz -o ${i/.vcf.gz/_2.vcf.gz} -- -t 'AN,AC,MAF,AC_Hom,AC_Het,ExcHet'
done

# open permissions on newly created VCFs
chmod u=rwx *vcf*

# step 3: filters for snp present in one and deleted in other
#exclude -T unless have preselected/filtered out other sites already (optional)
#Seperate N_PASS (will have one && and two N_PASS one for ea spp) Use NPASS but think about coverage across 

# him: first 510, bcftools view -i '%QUAL>=20' calls.bcf starting with 0, in bam. filelist, need at least 1 x cov in at least 408 (80% of my bams for that species). first gives initial vcf then moves into a second vcf, basically only change genome after -f
#will get .sites file of sites for all individuals. 

# my notes: will need to change values of # included
# view -m2 -M2 -v snps : keep only biallelic
# -i '%QUAL>=20' : include only sites with quality score over 20
# -i 'N_PASS(FMT/DP[0-39] >= 1) >= 32 && [ETC]' : keep only sites present in at least 80% (30/37 and 22/27) of samples for each species (must have at least 1x coverage at the site for it to be included) 
# -e INFO/MQ < 20 : exclude anything with mapping quality less than 20
# INFO/DP > 1000 : Exclude sites with greater than 1000 depth (FMT/DP = 'genotype read depth' while INFO/DP = 'mean depth per variant site')
# INFO/FS < 1e-100 : exclude sites with a Fisher's strand bias p-value below 1e-100
# INFO/BQBZ < 1e-100 : exclude sites with a Base Quality bias p-value below 1e-100
# INFO/MQBZ < 1e-100 : exclude sites with a Mapping Quality bias p-value below 1e-100
# INFO/RPBZ < 1e-100 : exclude sites with a Read Position bias p-value below 1e-100
# INFO/ExcHet < 1e-2' ; tests for excess heterozygosity (0 is bad, 1 is good)
# -M 2 :  exclude sites with more than 2 variants (exclude multiallelic snps)
# bcftools norm --atom-overlaps '*' -f C_lasius_ref.fa -m +any : what chatGPT thinks: "Overall, this part of the command is ensuring that the variants in the VCF file are properly normalized, taking into account any overlapping multiallelic variants"
# indels not explicitely excluded here because excluded from initial vcf gen
# side note: || means 'or'

#(CONSIDER: quality score filter, ie bcftools view -i '%QUAL>=20' calls.bcf; not sure if I have already accounted for quality or if there are other reasons not to do this...)

#bcftools view -i 'N_PASS(FMT/DP[0-21] >= 1) >= 17 && N_PASS(FMT/DP[22-34] >= 1) >= 10 && N_PASS(FMT/DP[35-52] >= 1) >= 14 && N_PASS(FMT/DP[53-59] >= 1) >= 5' AHVF1_unfiltered_scaffold_set${SLURM_ARRAY_TASK_ID}_2.vcf.gz | bcftools view -e 'QUAL<20 || INFO/MQ < 20 || INFO/DP > 1000 || INFO/FS < 1e-100 || INFO/BQBZ < 1e-100 || INFO/MQBZ < 1e-100 || INFO/RPBZ < 1e-100 || INFO/ExcHet < 1e-2' -M 2 -| bcftools norm --atom-overlaps '*' -f /hb/groups/kay_lab/genomes/C.lasius_genome_NCBI.fasta -m +any | bcftools query -f '%CHROM\t%POS\n' > AHVF1_filtered_q20_scaf${SLURM_ARRAY_TASK_ID}.sites

#echo "bcftools view -i 'N_PASS(FMT/DP[0-22] >= 1) >= 18 && N_PASS(FMT/DP[23-35] >= 1) >= 10 && N_PASS(FMT/DP[36-53] >= 1) >= 14 && N_PASS(FMT/DP[54-60] >= 1) >= 5' AHVF1_unfiltered_scaffold_set${SLURM_ARRAY_TASK_ID}_2.vcf.gz | bcftools view -e 'QUAL<20 || INFO/MQ < 20 || INFO/DP > 1000 || INFO/FS < 1e-100 || INFO/BQBZ < 1e-100 || INFO/MQBZ < 1e-100 || INFO/RPBZ < 1e-100 || INFO/ExcHet < 1e-2' -M 2 -| bcftools norm --atom-overlaps '*' -f /hb/groups/kay_lab/genomes/C.lasius_genome_NCBI.fasta -m +any | bcftools query -f '%CHROM\t%POS\n' > AHVF1_filtered_q20_scaf${SLURM_ARRAY_TASK_ID}.sites"


# next
#for i in {1..11}; do
#        cat AHVF1_filtered_q20_scaf${i}.sites
#done >> AHVF1_filtered.sites

#for i in {1..11}; do
#        echo "cat Filtered_scaf${i}.sites"
#done

echo "done!  └(^o^ )┘"
