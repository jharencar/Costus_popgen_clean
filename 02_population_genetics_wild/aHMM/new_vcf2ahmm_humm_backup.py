#!/usr/bin/env python3
# Maximilian Genetti (maxgenetti@gmail.com)

import os
import sys
import argparse
import numpy as np
from collections import defaultdict

'''
Reads in a vcf filtered for LD, indels and biallelic snps.
Reads in a file containing parent1, parent2, admixed followed by comma seperated names or indexes starting at 1.
Reads in a recombination map if supplied, or requires a default value in cM/mb.
Converts to aHMM format where it has reads for each allele or genotypes if -g flag is used.
'''

class VCF2aHMM :
    '''
    Args:
        Reads a filtered vcf, population assignment file, and optionally a recombination map file
        If no recombination map is given, it takes a mean recombination rate estimate.

    Returns:
        Ancestry HMM input file
    '''

    def __init__ (self) :
        '''constructor: saves attributes'''
        args = self.get_args() #reads command line
        self.vName = args.vcf
        self.mName = args.map
        self.pName = args.pop
        self.rate = args.rate *0.00000001 #convert cM/Mb to morgan/bp 
        self.minFreq = args.minDif
        self.minGTP = args.minGTP
        self.minGTA = args.minGTA
        self.minDP = args.minDP
        self.dist = args.dist
        self.distM = args.distM/100 #convert cM to morgans 
        self.geno = args.geno
        self.recMap = {0:0} 
        self.ancestry = {}
        self.table = []
        self.chrom = ''
        self.header = self.readVCFheader()
        self.ahmm = []
        self.current = 0
        self.readSamples()
        self.readRecMap()
        self.siteFreq = []

    def get_args(self):
        '''
        Notes:
            reads command line arguments
        '''
        parser = argparse.ArgumentParser(description='Read in files and paremeters')
        parser.add_argument('--vcf', type=str, help='name of the vcf')
        parser.add_argument('--map', type=str, default='check rate', help='name of the recombination map file (units in cM/Mb)')
        parser.add_argument('--pop', type=str, help='name of the population identity file')
        parser.add_argument('--rate', type=float, default=0, help='mean recombination rate estimate in cM/Mb')
        parser.add_argument('--dist', type=int, default=1000, help='minimum distance between sites in bp')
        parser.add_argument('--distM', type=float, default=0.0001, help='minimum distance between sites in Cm')
        parser.add_argument('--minDif', type=float, default=0.1, help='minimum frequency difference in parental genotypes')
        parser.add_argument('--minGTP', type=int, default=20, help='minimum number of genotype calls for each parent population')
        parser.add_argument('--minGTA', type=float, default=0.5, help='minimum rate of genotype calls for admixed population with -geno flag')
        parser.add_argument('--minDP', type=float, default=0.5, help='minimum rate of depth>1 for admixed population (if no -geno flag)')
        parser.add_argument('-geno', action='store_true', default=False, help='use admixed genotypes instead of read counts')
        return parser.parse_args()

    def readVCFheader (self) :
        '''
        Args:
            vcf file name
        Notes:
            stores header
        '''
        with open(self.vName) as fileV:
            header = []
            data = []
            # skip to VCF header
            line = fileV.readline()
            while not line.startswith('#CHROM') :
                line = fileV.readline()
            self.chrom = fileV.readline().split()[0] #reads chrom from first position in VCF

            fileV.close()
            return line.split()

    def readSamples (self) :
        '''
        Args:
            reads in the population assignment file
        Notes:
            stores dictionary of population indexes
        '''
        i = self.header.index('FORMAT') #for name listed samples

        with open(self.pName) as fileP:
            self.ancestry['parent1'] = []
            self.ancestry['parent2'] = []
            self.ancestry['admixed'] = []
            key = ''
            samples = []
            for line in fileP :
                try :
                    key = line.split()[0]
                    samples = line.split()[1].split(',')
                except :
                    return
                for sample in samples:
                    try :
                        self.ancestry[key].append(int(sample)) #read list of numbers
                    except :
                        try :
                            self.ancestry[key] += [item for item in range(int(sample.split("-")[0]), int(sample.split("-")[1])+1)] #read range of numbers
                        except :
                            self.ancestry[key] += [self.header.index(sample)-i] #read sample names

    def readVCF (self) :
        '''
        Args:
            vcf file name
        Notes:
            stores array of data and header
        '''

        with open(self.vName) as fileV:

            header = []
            data = []
            # skip to first fasta header
            line = fileV.readline()
            while not line.startswith('#CHROM') :
                line = fileV.readline()
            #print("Starting", file = sys.stderr)
            for line in fileV :
                self.table = line.split()
                self.convert()

    def readRecMap (self) :
        '''
        Args:
            reads in the recombination map file if supplied (format based on corbett-detig et al. 2014)
        Notes:
            Stores recombination map file as dictionary of all positions on chromosome. self.map[pos]= cM
            Can be very large.
        '''
        if self.mName != "check rate" :
            print('recombination rate read from %s' %(self.mName), file = sys.stderr)
            with open(self.mName) as fileM :
                """
                7 column : Species, CHROM, block, start cM, middle cM, end cM, rate
                outputs library of all positions on chromsome.
                """
                prevPos = 1
                prevDist = 0
                for line in fileM :
                    l = line.split()
                    if len(l) == 7 :
                        chrom = l[1]
                        if chrom == self.chrom :
                            pos = int(l[2].split("_")[1])+1
                            dist = [float(d) for d in l[3:6]]
                            p1 = pos-50000
                            p2 = pos+50000
                            d1 = dist[0]
                            d2 = dist[2]
                            m = max((d2-d1)/(p2-p1), 2e-10) #rate in M/bp
                            for i in range(p1,p2+1) :
                                self.recMap[i] = m*(i-p1) + d1
                            if p1 != prevPos :
                                pp1 = prevPos
                                pp2 = p1
                                dd1 = prevDist
                                dd2 = d1
                                mm = max((dd2-dd1)/(pp2-pp1), 2e-10)
                                for i in range(pp1,pp2+1) :
                                    self.recMap[i] = mm*(i-pp1) + dd1
                            prevDist = d2
                            prevPos = p2
                        else :
                            continue
                    elif len(l) == 3 :
                        chrom = l[0]
                        if chrom == self.chrom :
                            pos = int(l[1])
                            p1 = prevPos
                            d1 = prevDist
                            p2 = int(l[1])
                            d2 = float(l[2])
                            m = (d2-d1)/(p2-p1)
                            for i in range(p1,p2+1) :
                                self.recMap[i] = m*(i-p1) + d1
                            prevDist = d2
                            prevPos = p2
            self.rate = 2e-10

        elif self.rate == 0 :
            print('Set the --rate parameter or include a recombination map', file = sys.stderr)
            exit

    def convert (self) :
        '''
        Converts vcf into np array with values for ahmm input file
        '''
        newData = []
        i = self.header.index('FORMAT') #Find Format column
        if not self.geno :
            try:
                ad = self.table[i].split(':').index('AD') #Find allele depth in format
            except:
                print(self.table[i], file = sys.stderr)

        gt = self.table[i].split(':').index('GT') #Find genotype call in format
        prev = self.current #assume first position is zero

        site = self.table
        chrom = site[0]

        #calculate morgans
        pos = int(site[1])
        if pos - prev <= self.dist :
            return
        try :
            morgans = self.recMap[pos]-self.recMap[prev]
        except :
            if self.rate == 2e-10 :
                return
            else :
                recRate = self.rate #uses input rate or 1e-12 as default if over map bounds
                morgans = (pos - prev)*recRate
                #morgans = "%.12f" % morgans
        if morgans < self.distM :
            return

        parent1, parent2 = [0, 0], [0, 0]   #default zero counts for allele frequencies
        admixed = [] #empty list for allele frequencies

        for sample in self.ancestry['parent1'] :
            for s in site[i + sample].split(':')[gt].split('/') :
                if s == "0" or s == "1" :   #check all alleles
                    parent1[int(s)] += 1

        for sample in self.ancestry['parent2'] :
            for s in site[i + sample].split(':')[gt].split('/') :
                if s == "0" or s == "1" :   #check all alleles
                    parent2[int(s)] += 1

        if self.geno == True : #if -g flag for admixed genotypes
            for sample in self.ancestry['admixed'] :
                alleles = [0,0]
                for s in site[i + sample].split(':')[gt].split('/') :
                    if s == "0" or s == "1" :   #check all alleles
                        alleles[int(s)] += 1
                admixed += alleles

        else: #if no -geno flag count reads
            for sample in self.ancestry['admixed'] :
                alleles = [0,0]
                c = site[i + sample].split(':')[ad].split(",")
                try :
                    alleles[0] = int(c[0])
                    alleles[1] = int(c[1])
                except :
                    alleles = [0,0]
                admixed += alleles

        if self.geno == True and sum(admixed)/len(self.ancestry['admixed']) <= self.minGTA : #check admixed genotype calls
            return
        elif self.geno == False and np.count_nonzero(admixed)/len(self.ancestry['admixed']) <= self.minDP : #check admixed mean depth
            return
        elif sum(parent1) >= self.minGTP  and sum(parent2) >= self.minGTP : #check parental genotype calls
            if -self.minFreq <= parent1[0]/(parent1[0]+parent1[1]) - parent2[0]/(parent2[0]+parent2[1]) <= self.minFreq : #genotype frequency difference check
                return
            else :
                a = [chrom] + [str(pos)] + parent1 + parent2 + [morgans] + admixed
                self.siteFreq += [abs(parent1[0]/(parent1[0]+parent1[1]) - parent2[0]/(parent2[0]+parent2[1]))]
                self.current = int(site[1]) #update current last position
                print('\t'.join([str(x) for x in a]), file = sys.stdout)
        else :
            return

def main() :
    aHMM = VCF2aHMM()
    aHMM.readVCF()
    print("Sites in Panel             : " + str(len(aHMM.siteFreq)), file = sys.stderr)
    print("Mean Frequency Difference  : " + str(np.mean(aHMM.siteFreq)), file = sys.stderr)
    print("Frequency Difference stDev : " + str(np.std(aHMM.siteFreq)), file = sys.stderr)
if __name__ == "__main__":
    main()
