#!/usr/bin/env python

import argparse
import pandas as pd 

def parse_args():
    parser = argparse.ArgumentParser(prog = 'pair_hse_helitrons.py', conflict_handler = 'resolve')
    parser.add_argument('-genomes', type = str, required = True, help = '=> organism genome ids')
    parser.add_argument('-o', type = str, required = True, help = '=> path/to/outfolder')
    return(parser.parse_args())

def merge_overlaps(sorted_by_lower_bound):
    """
    if HSE coords overlap, merge the coordinates
    """
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return(merged)

def summarize_coverage(mainfolder_out,organism_dir):
    helitron_file = pd.read_csv('{}/{}/{}_helitron_HSEs.txt'.format(mainfolder_out,organism_dir,organism_dir))[['genoName','genoStart','genoEnd']];
    limitless_fimo = pd.read_csv('{}/{}/{}_extendedFIMO.txt'.format(mainfolder_out,organism_dir,organism_dir))[['sequence_name','start','stop']];
    chrom_sizes = pd.read_csv('{}/{}/{}.compiled.chrom.sizes'.format(mainfolder_out,organism_dir,organism_dir), sep='\t');
    genome_size = chrom_sizes['size'].sum()
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
            merged_coordinates_hse += merged_coordinates_hse_chr #merged_coordinates_hse.extend(merged_coordinates_hse_chr)
                
            chr_subsetter_heli = helitron_file.loc[helitron_file['genoName'] == chromosome_ids];
            coordinate_tuple_heli=list(zip(chr_subsetter_heli['genoStart'], chr_subsetter_heli['genoEnd']));
            sorted_coordinate_lower_bound_heli = sorted(coordinate_tuple_heli, key=lambda tup: tup[0]);
            merged_coordinates_heli_chr=merge_overlaps(sorted_coordinate_lower_bound_heli);
            merged_coordinates_heli += merged_coordinates_heli_chr #merged_coordinates_heli.extend(merged_coordinates_heli_chr)
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
    return(len(limitless_fimo),genome_coverage_heli,genome_coverage_hse,sum_num_HSEs_in_heli,genome_size)

def main():
    args = parse_args()
    organisms = args.genomes
    outfolder = args.o
    with open("{}".format(organisms), 'r') as organism_list:
        with open('{}/HSE_in_out_summary.txt'.format(outfolder), 'w') as HSE_in_out_summary:
            HSE_in_out_summary.write('Org_Assembly\tHSEs.in.merged_overlap.heli\tHSEs.outside\tHSEs.total\tratio.HSEs.in.heli\tratios.HSEs.outside\tnum.HSEs.bases\tnum.heli.bases\tgenome.Size.bases\tHSE.perc.of.genome\tHelitron.perc.of.genome\n')
            for organism in organism_list:
                org = organism.strip()
                (total_hses, heli_bases, hse_bases, num_HSEs_in_helis,genome_bases) = summarize_coverage(outfolder,org)
                HSE_in_out_summary.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                                                    org,
                                                                    num_HSEs_in_helis,
                                                                    total_hses-num_HSEs_in_helis,
                                                                    total_hses,
                                                                    num_HSEs_in_helis/total_hses,
                                                                    1-(num_HSEs_in_helis/total_hses),
                                                                    hse_bases,
                                                                    heli_bases,
                                                                    genome_bases,
                                                                    hse_bases/genome_bases,
                                                                    heli_bases/genome_bases))

if __name__ == '__main__':
    main()