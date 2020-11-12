#!/usr/bin/env python

import sys

import glob, os
import pandas as pd 
import numpy as np

import gzip
import os.path
import urllib
from pathlib import Path

from itertools import tee
import statistics
from statistics import mean
import csv

organism=sys.argv[2]
outfolder = 'refseq_out/'

def merge_overlaps(sorted_by_lower_bound):
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return(merged)

def HSE_Helitron_counter(mainfolder_out,organism_dir):
    limitless_fimo = pd.read_csv('refseq_motif_files/{}_genomic.fna.motif.tsv'.format(sys.argv[1]), sep='\t').sort_values(['sequence_name', 'start'], ascending=[True, True]);
    limitless_fimo = limitless_fimo.drop_duplicates(subset = ['sequence_name','start','stop'])
    limitless_fimo.to_csv('{}{}_extendedFIMO.txt'.format(mainfolder_out,organism_dir))
    with open('rmsk.out.refseq/{}_genomic.fna.ori.out'.format(sys.argv[1])) as f:
        helitron_input = f.readlines()
    helitron_file = pd.DataFrame([i.split()[0:10] for i in helitron_input])[[4,5,6,8,9]]
    helitron_file.columns = ['genoName','genoStart','genoEnd','strand','repName']
    helitron_file.repName.str.lower()
    hit_indices = [i for i, s in enumerate(helitron_file.repName.str.lower()) if 'heli' in s]
    helitron_file = helitron_file.iloc[hit_indices,:]
    helitron_file[['genoStart','genoEnd']] = helitron_file[['genoStart','genoEnd']].astype('int')
    num_HSEs = []
    limitless_fimo_start_stop = limitless_fimo[['sequence_name','start','stop']];
    helitron_file_start_stop = helitron_file[['genoStart','genoEnd','genoName']]
    for index,row_vals in helitron_file_start_stop.iterrows():
        chr_subsetter =  limitless_fimo_start_stop.loc[limitless_fimo_start_stop['sequence_name'] == row_vals['genoName']];
        num_HSEs.append(len(chr_subsetter.loc[(chr_subsetter['start'] >= row_vals['genoStart']) & (chr_subsetter['stop'] <= row_vals['genoEnd']), 'start'].values))
    helitron_file['Num_HSEs'] = num_HSEs
    helitron_file.to_csv('{}{}_helitron_HSEs.txt'.format(mainfolder_out,organism_dir))
    num_HSEs = []
    summary_df = pd.DataFrame(columns = ['#name','number.w.HSEs','number.total','fraction.w.HSEs','grouped.number.total.HSEs'])
    unique_list = helitron_file['repName'].unique();
    for heli_ind in range(0,len(unique_list)):
        helitron_types = unique_list[heli_ind]
        subsetter_heli= helitron_file[helitron_file['repName']==helitron_types]
        subsetter_heli_withHSE = subsetter_heli['Num_HSEs'].astype(bool).sum(axis=0)
        len_subsetter_heli = len(subsetter_heli['repName'])
        subsetter_heli_HSE_sum = subsetter_heli['Num_HSEs'].sum()
        summary_df.loc[heli_ind] = [helitron_types,subsetter_heli_withHSE,len_subsetter_heli,subsetter_heli_withHSE/len_subsetter_heli,subsetter_heli_HSE_sum]
    summary_df.sort_values('#name').to_csv('{}{}_Helitrons.HSEs.summary.txt'.format(mainfolder_out,organism_dir))

def summarize_coverage(mainfolder_out,organism_dir):
    helitron_file = pd.read_csv('{}{}_helitron_HSEs.txt'.format(mainfolder,organism_dir))[['genoName','genoStart','genoEnd']];
    limitless_fimo = pd.read_csv('{}{}_extendedFIMO.txt'.format(mainfolder,organism_dir))[['sequence_name','start','stop']];
    try:
        num_HSEs_in_heli = []
        merged_coordinates_heli = []
        merged_coordinates_hse = []
        for chromosome_ids in helitron_file['genoName'].unique():
            limitless_fimo['start'] = limitless_fimo['start'] #bedtools is zero-index, fimo is 1
            limitless_fimo['stop'] = limitless_fimo['stop'] #unchanged, bedtools is 1-index for stop, so is fimo
            chr_subsetter_HSE =  limitless_fimo.loc[limitless_fimo['sequence_name'] == chromosome_ids];
            coordinate_tuple_hse=list(zip(chr_subsetter_HSE['start'], chr_subsetter_HSE['stop']));
            sorted_coordinate_lower_bound_hse = sorted(coordinate_tuple_hse, key=lambda tup: tup[0]);
            merged_coordinates_hse_chr=merge_overlaps(sorted_coordinate_lower_bound_hse);
            merged_coordinates_hse += merged_coordinates_hse_chr#merged_coordinates_hse.extend(merged_coordinates_hse_chr)
                
            chr_subsetter_heli = helitron_file.loc[helitron_file['genoName'] == chromosome_ids];
            coordinate_tuple_heli=list(zip(chr_subsetter_heli['genoStart'], chr_subsetter_heli['genoEnd']));
            sorted_coordinate_lower_bound_heli = sorted(coordinate_tuple_heli, key=lambda tup: tup[0]);
            merged_coordinates_heli_chr=merge_overlaps(sorted_coordinate_lower_bound_heli);
            merged_coordinates_heli += merged_coordinates_heli_chr#merged_coordinates_heli.extend(merged_coordinates_heli_chr)
            for i in merged_coordinates_heli_chr:
                num_HSEs_in_heli.append(len(chr_subsetter_HSE.loc[(chr_subsetter_HSE['start'] >= i[0]) & (chr_subsetter_HSE['stop'] <= i[1]), 'start'].values))
            
        genome_coverage_hse = sum(([i[1]-i[0]+1 for i in merged_coordinates_hse]))
        #print("Bases covered by HSEs (merged overlaps):", (genome_coverage_hse))
            
        genome_coverage_heli = sum([i[1]-i[0] for i in merged_coordinates_heli])
        #print("Bases covered by Helitrons (merged overlaps):", (genome_coverage_heli))
            
        sum_num_HSEs_in_heli = sum(num_HSEs_in_heli)
        #print("Number of HSEs in Helitrons (merged overlaps):", (sum_num_HSEs_in_heli))
    except:
        merged_coordinates_heli = 0
        merged_coordinates_hse = 0
        genome_coverage_heli = 0
        genome_coverage_hse = 0
        sum_num_HSEs_in_heli = 0
    return(len(helitron_file),merged_coordinates_heli,len(limitless_fimo),merged_coordinates_hse,genome_coverage_heli,genome_coverage_hse,sum_num_HSEs_in_heli)

def full_organism_summary(mainfolder_out,summary_cov_path,organism_dir):
    my_futurefile = Path('{}{}_Helitrons.HSEs.summary.txt'.format(mainfolder_out,organism_dir))
    if my_futurefile.is_file():
        summary_final = pd.read_csv(my_futurefile)[['#name','number.w.HSEs','number.total','fraction.w.HSEs','grouped.number.total.HSEs']]
        if len(summary_final) >= 1:
            total_helis = summary_final['number.total'].sum()
            try:
                summary_final['grouped.total.to.number.w.HSEs'] = summary_final['grouped.number.total.HSEs']/summary_final['number.w.HSEs'];
                summary_final.to_csv('{}{}_Helitrons.HSEs.summary.txt'.format(mainfolder_out,organism_dir));
                limitless_fimo_final = pd.read_csv('{}{}_extendedFIMO.txt'.format(mainfolder_out,organism_dir), sep=',');
                #print("test")
                summary_sizes_heli_HSE = pd.read_csv('{}{}'.format(mainfolder_out,summary_cov_path),sep='\t')
                HSEs_in_merged_heli = summary_sizes_heli_HSE[summary_sizes_heli_HSE['Assembly'] == organism_dir]['num.HSEs.in.merged_overlapping.heli'].values[0]
                total_merged_heli = summary_sizes_heli_HSE[summary_sizes_heli_HSE['Assembly'] == organism_dir]['Total.Helitrons.merged_overlap'].values[0]
                total_heli_bases = summary_sizes_heli_HSE[summary_sizes_heli_HSE['Assembly'] == organism_dir]['Helitrons.merged_overlaps.bases'].values[0]
                total_hse_bases = summary_sizes_heli_HSE[summary_sizes_heli_HSE['Assembly'] == organism_dir]['HSEs.merged_overlaps.bases'].values[0]
                #chrom_sizes = pd.read_csv('{}{}.chrom.sizes'.format(mainfolder,organism_dir,organism_dir), sep='\t', header = None);
                #print("test")
                genome_size = int(open('{}{}.fa.repeatmasker.tbl'.format(mainfolder,organism_dir)).readlines()[3].split()[2])#chrom_sizes[[1]].values.sum()
                #print(genome_size)
                org_name = organism_dir.split('.')[0]#conversion_table.loc[conversion_table['Assembly']==organism_dir, 'Common'].values[0];
                org_assembly = organism_dir.split('.')[1]
                summary_final_sum = summary_final['grouped.number.total.HSEs'].sum()
                hses_notin_mergedheli = len(limitless_fimo_final)-HSEs_in_merged_heli
                total_hses = len(limitless_fimo_final)
                #print("test")
                ratio_hses_in_mergedheli_to_total = HSEs_in_merged_heli/len(limitless_fimo_final)
                ratio_hses_outside_to_total = (len(limitless_fimo_final)-HSEs_in_merged_heli)/len(limitless_fimo_final)
                ratio_hse_genome = total_hse_bases/genome_size
                ratio_heli_genome = total_heli_bases/genome_size
            #print("test")
            except:
                summary_final_sum = 0
                HSEs_in_merged_heli = 0
                hses_notin_mergedheli = 0
                total_hses = 0
                ratio_hses_in_mergedheli_to_total = 0
                ratio_hses_outside_to_total = 0
                
                total_merged_heli = 0
                total_hse_bases= 0
                total_heli_bases= 0
                genome_size=0
                ratio_hse_genome=0
                ratio_heli_genome=0
            return(org_name,org_assembly,summary_final_sum,HSEs_in_merged_heli,hses_notin_mergedheli,total_hses,ratio_hses_in_mergedheli_to_total,ratio_hses_outside_to_total,total_helis,total_merged_heli,total_hse_bases,total_heli_bases,genome_size,ratio_hse_genome,ratio_heli_genome)

HSE_Helitron_counter(outfolder,organism)
