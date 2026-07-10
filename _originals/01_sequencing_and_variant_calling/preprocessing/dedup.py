# Script to remove duplicates from all samples in a folder
# By: Julia Harencar
# example of functioning one line (one sample R1 and R2) version:
## clumpify.sh in1=SE5869_SA49746_S1_L002_R1_001.fastq.gz in2=SE5869_SA49746_S1_L002_R2_001.fastq.gz out1=SA49746_S1_clumped1.fq.gz out2=SA49746_S1_clumped2.fq.gz dedupe subs=0 reorder

import os, sys
import glob
import subprocess

#get a list of prefixes for the sampling
files_for_prefix = glob.glob("*_R1_001.fastq.gz")
samples = map(lambda each:each.split("_L001")[0], files_for_prefix)

for sample in samples:
    pe1 = str(sample) + "_L001_R1_001.fastq.gz"
    pe2 = str(sample) + "_L001_R2_001.fastq.gz"
    out1 = str(samples) + "_clumped1.fq.gz"
    out2 = str(samples) + "_clumped2.fq.gz"
    null = "null"
    cmd = ["clumpify.sh"," in1=",pe1," in2=",pe2," out1=",out1," out2=",out2," dedupe subs=0 reorder"]
    cmd = "".join(cmd)
    print (cmd)
    process = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
    output, error = process.communicate()
        
print("done └(^o^ )┘")