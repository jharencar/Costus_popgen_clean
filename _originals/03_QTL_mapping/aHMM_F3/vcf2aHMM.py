#!/usr/bin/env python3
# Max Genetti (mgenetti)

import os
import sys
import argparse
import numpy as np
from collections import defaultdict
#import allel

'''
Reads in a vcf filtered for LD, indels and biallelic snps.
Reads in a tsv file containing parent1, parent2, admixed followed by names or indexes starting at 1.
Reads in a recombination map if supplied, or requires a default value.
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
        self.minGT = args.minGT
        self.minDP = args.minDP
        self.dist = args.dist
        self.geno = args.geno
        self.recMap = {}
        self.ancestry = {}
        self.table = []
        self.header = ''
        self.ahmm = []

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
        parser.add_argument('--minDif', type=float, default=0.2, help='minimum frequency difference in parental genotypes')
        parser.add_argument('--minGT', type=float, default=0.5, help='minimum mean rate of genotype calls each parent population (and admixed with -geno flag')
        parser.add_argument('--minDP', type=float, default=0.5, help='minimum mean depth for admixed population (if no -g flag)')
        parser.add_argument('-geno', action='store_true', default=False, help='use admixed genotypes instead of read counts')
        return parser.parse_args()

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

            self.header = line.split()

            for line in fileV :
                self.table.append(line.split())

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
            for line in fileP :
                key = line.split()[0]
                samples = line.split()[1].split(',')
                for sample in samples:
                    try :
                        self.ancestry[key].append(int(sample))
                    except :
                        try :
                            self.ancestry[key] += [item for item in range(int(sample.split("-")[0]), int(sample.split("-")[1])+1)]
                        except :
                            self.ancestry[key] += [self.header.index(sample)-i]

    def readRecMap (self) :
        '''
        Args:
            reads in the recombination map file if supplied
        Notes:
            stores recombination map file as dictionary of all positions. Can be very large.
        '''
        if self.mName != "check rate" :
            with open(self.mName) as fileM :
                prevPos = 0
                prevRate = 0
                for line in fileM :
                    l = line.split()
                    pos = int(l[0])
                    rate = float(l[1]) * 0.00000001 #convert cM/Mb to morgan/bp
                    self.recMap[pos] = rate / 1000000
                    m = (rate-prevRate)/(pos-prevPos) 
                    b = (prevPos*rate-pos*prevRate)/(prevPos-pos)
                    if pos-prevPos > 1:
                        for i in range(prevPos+1,pos) :
                            self.recMap[i] = m*i + b
                    prevRate = rate
                    prevPos = pos
                self.rate == 1e-12
        elif self.rate == 0 :
            print('Set the --rate parameter or include a recombination map')
            exit

    def convert (self) :
        '''
        Converts vcf into np array with values for ahmm input file
        '''
        newData = []
        
        i = self.header.index('FORMAT') #Find Format column
        dp = self.table[0][i].split(':').index('DP') #Find read depth in format
        ad = self.table[0][i].split(':').index('AD') #Find read depth in format
        gt = self.table[0][i].split(':').index('GT') #Find genotype call in format
        prev = 0 #assume first position is zero

        for site in self.table : #Check all sites

            chrom = site[0]

            #calculate morgans
            pos = int(site[1])
            recRate = self.recMap.setdefault(pos, self.rate) #uses input rate or 1e-12 as default if over map bounds

            morgans = (pos - prev)*recRate
            morgans = "%.12f" % morgans
            
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
            
            else: #if no -g flag count reads
                for sample in self.ancestry['admixed'] :
                    alleles = [0,0]
                    c = site[i + sample].split(':')[ad].split(",")
                    alleles[0] = int(c[0])
                    alleles[1] = int(c[1])

                    admixed += alleles

            if self.geno == True and sum(admixed)/len(self.ancestry['admixed']) < self.minGT : #check admixed genotype calls
                continue
            elif self.geno == False and sum(admixed)/len(self.ancestry['admixed']) < self.minDP : #check admixed mean depth
                continue
            if sum(parent1)/len(self.ancestry['parent1']) >= self.minGT  and sum(parent2)/len(self.ancestry['parent2']) >= self.minGT : #check parental genotype calls

                if -self.minFreq <= parent1[0]/(parent1[0]+parent1[1]) - parent2[0]/(parent2[0]+parent2[1]) <= self.minFreq : #genotype frequency difference check
                    continue
                elif -self.minFreq <= parent1[1]/(parent1[0]+parent1[1]) - parent2[1]/(parent2[0]+parent2[1]) <= self.minFreq : #genotype frequency difference check
                    continue
                if pos - prev >= self.dist : #minimum distance
                    a = [chrom] + [str(pos)] + parent1 + parent2 + [morgans] + admixed
                    newData.append(a)
                    prev = pos #set current to previous
                else :
                    continue
            else :
                continue

        self.ahmm = np.asarray(newData)


def main() :

    aHMM = VCF2aHMM()
    aHMM.readVCF()
    aHMM.readRecMap()
    aHMM.readSamples()
    aHMM.convert()
    
    for line in aHMM.ahmm :
        print('\t'.join([str(x) for x in line]))


if __name__ == "__main__":
    main()
