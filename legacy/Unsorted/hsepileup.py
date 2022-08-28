#!/usr/bin/env python
import numpy as np
import sys
import os.path
import pandas as pd
import os

n=1000
def count_range(list1, l, r):
    finallist = []
    for x in list1:
        if l <= x <= r:
            finallist.append(x)
    if len(finallist) > 0:
        return(finallist)
    else:
        return([0])
def brian_counts(p,n_limit,sliding_i, size_mod):
    l=0
    for j in range(0,len(p)):
        l+=(n_limit+sliding_i-size_mod > p[j] > sliding_i)
    return(l)
def count_greater(p,lim):
    l=0
    for j in range(0,len(p)):
        l+=(p[j] > lim)
    return(l)

outfolder = "/home/bvtsu/data/i5k_out/"

organism_dir = sys.argv[1]
try:
    assert(os.path.isfile(outfolder+organism_dir+".compiled.chrom.sizes"))
    assert(os.path.isfile(outfolder+organism_dir+"_extendedFIMO.txt"))
    if os.path.isfile(outfolder+organism_dir+".HSE.counts.txt") == True:
        os.remove(outfolder+organism_dir+".HSE.counts.txt")
    chrom_sizes = pd.read_csv(outfolder+organism_dir+".compiled.chrom.sizes", sep = "\t");
    genome_size = chrom_sizes[["size"]].values.sum()
    limitless_fimo = pd.read_csv(outfolder+organism_dir+"_extendedFIMO.txt")[["sequence_name","start"]];
    for chromosome_ids in limitless_fimo["sequence_name"].unique():
        chr_subsetter_HSE =  [x for x in limitless_fimo[limitless_fimo["sequence_name"] == chromosome_ids]["start"]];
        length_of_chr = chrom_sizes[chrom_sizes["chrom"] == chromosome_ids]["size"].iloc[0];
        with open(outfolder+organism_dir+".HSE.counts.txt", "a+") as hse_counts:
            hse_counts.writelines([str(brian_counts(chr_subsetter_HSE,n,i,14))+"\n" for i in range(1,length_of_chr,n)])
        chr_subsetter_HSE = pd.DataFrame()
    limitless_fimo = pd.DataFrame()
    chrom_sizes =  pd.DataFrame()
    m = 0
    standard_dev = 0
    with open(outfolder+organism_dir+".HSE.counts.txt", "r") as hse_counts:
        HSE_counts = [int(x) for x in hse_counts.readlines()];
    m = np.mean(HSE_counts)
    standard_dev = np.std(HSE_counts)
    indices_above5SD = count_greater(HSE_counts,(m+5*standard_dev))
    indices_above10SD = count_greater(HSE_counts,(m+10*standard_dev))
    indices_above10 = count_greater(HSE_counts,10)
    indices_above25 = count_greater(HSE_counts,25)
    indices_above50 = count_greater(HSE_counts,50)
    indices_above75 = count_greater(HSE_counts,75)
    indices_above100 = count_greater(HSE_counts,100)
    HSE_counts = []
    with open(outfolder+"i5k.HSE.pileup.txt", "a+") as pileup:
        pileup.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(organism_dir,str(genome_size),str(m),str(standard_dev),str(indices_above5SD),str(indices_above10SD),str(indices_above10),str(indices_above25),str(indices_above50),str(indices_above75),str(indices_above100)))
except:
    print("Nothing")
