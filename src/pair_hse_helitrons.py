#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd 

from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(prog = 'pair_hse_helitrons.py', conflict_handler = 'resolve')
    parser.add_argument('-genome', type = str, required = True, help = '=> organism genome id')
    parser.add_argument('-o', type = str, required = True, help = '=> path/to/outfolder')
    return(parser.parse_args())

def HSE_Helitron_counter(mainfolder_out,organism_dir):
    limitless_fimo = pd.read_csv('{}/{}.motif.tsv'.format(mainfolder_out,organism_dir), sep='\t'
                        ).sort_values(['sequence_name', 'start'], ascending=[True, True]);
    limitless_fimo = limitless_fimo.drop_duplicates(subset = ['sequence_name','start','stop'])
    limitless_fimo.to_csv('{}/{}_extendedFIMO.txt'.format(mainfolder_out,organism_dir))
    with open('{}/rmsk.txt'.format(mainfolder_out,organism_dir)) as f:
        helitron_input = f.readlines()
    helitron_file = pd.DataFrame([i.split()[0:11] for i in helitron_input])[[5,6,7,9,10]]
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
    helitron_file.to_csv('{}/{}_helitron_HSEs.txt'.format(mainfolder_out,organism_dir))
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
    summary_df.sort_values('#name').to_csv('{}/{}_Helitrons.HSEs.summary.txt'.format(mainfolder_out,organism_dir))

def main():
    args = parse_args()
    organism= args.genome
    outfolder = args.o
    HSE_Helitron_counter(outfolder,organism)

if __name__ == '__main__':
    main()