# filter VCF by mapability (uniqueness)

mappability_k: 150
mappability_merge: 100
mappability_min: 1


genmap index -F assembly.fasta \
        -I assembly.genmap_index 
genmap map -K $mappability_k -E 0 -I assembly.genmap_index/ \
    -O genmap/assembly.genmap150_out \
        -bg -T 10 -v &> genmap/assembly.genmap.150.log 
sort -k1,1 -k2,2n genmap/assembly.genmap150_out.bedgraph \
    genmap/assembly.genmap150_out.s.bg 
awk 'BEGIN{OFS="\t";FS="\t"} { if($4>='$mappability_min'} print $1,$2,$3 }' \
    assembly.genmap150_out.s.bg > temporary_map.bed
bedtools sort -i temporary_map.bed | \
    bedtools merge -d $mappability_merge -i - > callable_sites.bed

cat assembly.fasta.fai | cut -f1,2  | sort -k1,1  > sorted.chrom.sizes
bedtools complement -i callable_sites.bed -g sorted.chrom.sizes > non_callable_sites.bed
bcftools view -R {input.bed} -m2 -M2 \
            {input.vcf} -O u -o {output.vcf}


# can add -S 20 to use less memory