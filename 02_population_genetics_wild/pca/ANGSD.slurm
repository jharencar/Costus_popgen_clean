#!/bin/bash
#SBATCH --partition 128x24 # this means each node (there are 19 nodes in this partition) has 128 GB RAM and 24 cores/cpus 
#SBATCH --job-name=ANGSD_beagles# Job name
#SBATCH --time=10:00:00  # ran in roughly 9 hours for all 96 merged bams of plate 1
#SBATCH --output=ANGSD_beagles_output.txt
#SBATCH --mail-type=END, FAIL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jharenca@ucsc.edu  # Where to send mail
#SBATCH --nodes=3 

## load ANGSD
module load angsd/gnu-930 samtools/samtools-1.10

## change into the directory with .bams
cd /hb/scratch/jharenca/vill_alle_plate1/newdata/whole_plate_mapping/ANGSD

## call the code to run


# helpful summary of some of the main diffs between GATK and samtools: https://www.biostars.org/p/57149/
# (I am very unsure about the SNP_pval param. ask Merly or Molly or...) 
# to generate beagle files for NGSadmix and PCAngsd. 
echo "angsd -GL 1 -out genolike -sites -nThreads 24 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam bam.filelist"
angsd -GL 1 \ # -GL: Estimate genotype likelihoods. Options - 1: SAMtools 2: GATK 3: SOAPsnp 4: SYK
-out genolike \ # file output type
-nThreads 144 \
-doGlf 2 \ # -doGlf: glf is a simple format for genotype likelihoods. Options-	1: binary glf (10 log likes) .glf.gz 2: beagle likelihood file  .beagle.gz 3: binary 3 times likelihood	.glf.gz 4: text version (10 log likes)	.glf.gz
-doMajorMinor 1 \ # -doMajorMinor: infer the major and minor alleles. Options - 1: Infer major and minor from GL (THIS is what I want) 2: Infer major and minor from allele counts 3: use major and minor from a file (requires -sites file.txt) 4: Use reference allele as major (requires -ref) 5: Use ancestral allele as major (requires -anc)
-SNP_pval 1e-6 \ # -SNP_pval: p value to be used for snp calling, 1e-6 seems to be a good default (talk to Merly?)
-minMapQ 30 \ # -minMapQ minimum MapQ score, which in bams created by bwa is equivalent to Phred score, 99.9% accuracy is minMapQ 30 (1 in 1000 incorect probability) Start with at least 20 (99%). 
-doMaf 1 \ # -minMaf remove sites with MAF below x (is MAF major or minor freq? must be minor... in which case it seems like a good idea to remove below .05 right? those are likely erroneous bp calls? WOuld be good to know how often this is even an issue in my data... )
-bam bam.filelist

#to run NGSadmix (located in the mics folder in angsd)
# -K: number of clusters (i chose 17 because I have that many species in full plate, nice not to have to use Evanno...)
# -P: number of threads
# -minMaf: Minimum minor allele frequency (5% is default, seems fine??)
#NGSadmix -likes input.gz -K 17 -P 10 -o plate_1_NGSadmix_out -minMaf 0.05 

#ANGSD documentation: http://www.popgen.dk/angsd/index.php/ANGSD

#NGSAdmix documentation: http://www.popgen.dk/software/index.php/NgsAdmixTutorial
# full tutorial (NGSAdmix) to unpack: $ANGSD -bam all.files -GL 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 2e-6 -minMapQ 30 -minQ 20 -minInd 25 -minMaf 0.05 -doGlf 2 -out $IN_DIR/Demo1input.gz -P 5


### #generate pop info file example

#generate bam filelist (angsd input to gen beagle file)
ls *ALLE* *VILL* *HYB* > AVH_only.bam.filelist

#generate pop info file
paste -d " " <( cut -f 2 -d "_" AVH_only.bam.filelist) <(cut -f 1 -d"_" AVH_only.bam.filelist) > AVH_only.pop.info

#To view a summary of the population information file, cut the first column, sort and count:
cut -f 1 -d " " AVH_only.pop.info | sort | uniq -c

# make a population label file and place it in the output directory
cut -f1 -d" " AVH_only.pop.info > AVH_only.poplabel


#generate bam filelist (angsd input to gen beagle file)
ls *ALLE* *VILL* *HYB* *GUAN* > GAVH.bam.filelist

#generate pop info file
paste -d " " <( cut -f 2 -d "_" GAVH.bam.filelist) <(cut -f 1 -d"_" GAVH.bam.filelist) > GAVH.pop.info

#To view a summary of the population information file, cut the first column, sort and count:
cut -f 1 -d " " GAVH.pop.info | sort | uniq -c

# make a population label file and place it in the output directory
cut -f1 -d" " GAVH.pop.info > GAVH_only.poplabel


GAVH.genolike.beagle.gz





