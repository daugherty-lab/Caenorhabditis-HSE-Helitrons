#size_summary.write('Assembly\tTotal.Helitrons\tTotal.Helitrons.merged_overlap\tTotal.HSEs\tTotal.HSEs.merged_overlap\tHelitrons.merged_overlaps.bases\tHSEs.merged_overlaps.bases\tnum.HSEs.in.merged_overlapping.heli\n')

conversion_table =  pd.read_csv("/Users/briantsu/Public/notebooks/HSE_project/assembly_conversion_list.txt");
organism_list = ['cb3']
with open('{}HSE_in_out_summary.txt'.format(mainfolder), 'w') as HSE_in_out_summary:
    HSE_in_out_summary.write('Org_name\tOrg_Assembly\tHSEs.in.UNMERGED.heli\tHSEs.in.merged_overlap.heli\tHSEs.outside\tHSEs.total\tratio.HSEs.in.heli\tratios.HSEs.outside\tnumber.total.all.Helitrons\tnumber.total.merged.Helitrons\tnum.HSEs.bases\tnum.heli.bases\tgenome.Size.bases\tHSE.perc.of.genome\tHelitron.perc.of.genome\n')
    for organism_dir in organism_list:#total_org_list:#organism_list:
        my_futurefile = Path('{}{}/{}_Helitrons.HSEs.summary.txt'.format(mainfolder,organism_dir,organism_dir))
        if my_futurefile.is_file():
            summary_brian = pd.read_csv(my_futurefile)[['#name','number.w.HSEs','number.total','fraction.w.HSEs','grouped.number.total.HSEs','grouped.total.to.number.w.HSEs']]#'{}{}/{}_Helitrons.HSEs.summary.txt'.format(mainfolder,organism_dir,organism_dir));
            if len(summary_brian) >= 1:
                #summary_brian['grouped.total.to.number.w.HSEs'] = summary_brian['grouped.number.total.HSEs']/summary_brian['number.w.HSEs'];
                #summary_brian.to_csv('{}{}/{}_Helitrons.HSEs.summary.txt'.format(mainfolder,organism_dir,organism_dir));
                limitless_fimo_brian = pd.read_csv('{}{}/{}_extendedFIMO.txt'.format(mainfolder,organism_dir,organism_dir), sep=',');
                summary_sizes_heli_HSE = pd.read_csv('{}summary_sizes.txt'.format(mainfolder),sep='\t')
                HSEs_in_merged_heli = summary_sizes_heli_HSE[summary_sizes_heli_HSE['Assembly'] == organism_dir]['num.HSEs.in.merged_overlapping.heli'][0]
                total_merged_heli = summary_sizes_heli_HSE[summary_sizes_heli_HSE['Assembly'] == organism_dir]['Total.Helitrons.merged_overlap'][0]
                total_heli_bases = summary_sizes_heli_HSE[summary_sizes_heli_HSE['Assembly'] == organism_dir]['Helitrons.merged_overlaps.bases'][0]
                total_hse_bases = summary_sizes_heli_HSE[summary_sizes_heli_HSE['Assembly'] == organism_dir]['HSEs.merged_overlaps.bases'][0]
                chrom_sizes = pd.read_csv('{}{}/{}.chrom.sizes'.format(mainfolder,organism_dir,organism_dir), sep='\t', header = None);
                genome_size = chrom_sizes[[1]].values.sum()
                org_name = conversion_table.loc[conversion_table['Assembly']==organism_dir, 'Common'].values[0];
                HSE_in_out_summary.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(org_name,
                                                                 organism_dir,
                                                                 summary_brian['grouped.number.total.HSEs'].sum(),
                                                                 HSEs_in_merged_heli,
                                                                 len(limitless_fimo_brian)-HSEs_in_merged_heli,
                                                                 len(limitless_fimo_brian),
                                                                 HSEs_in_merged_heli/len(limitless_fimo_brian),
                                                                 (len(limitless_fimo_brian)-HSEs_in_merged_heli)/len(limitless_fimo_brian),
                                                                 summary_brian['number.total'].sum(),
                                                                 total_merged_heli,
                                                                 total_hse_bases,
                                                                 total_heli_bases,
                                                                 total_hse_bases/genome_size,
                                                                 total_heli_bases/genome_size))
#summary_df