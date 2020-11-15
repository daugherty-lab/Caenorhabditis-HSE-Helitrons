using StatPlots
using DataFrames
using Gadfly
using RDatasets
using CSV

mainfolder = "/Users/briantsu/Public/notebooks/HSE_project/"
HSE_in_out_summary = sort(CSV.read("/Users/briantsu/Public/notebooks/HSE_project/HSE_in_out_summary.txt", delim = '\t')[[2,end-2]],2)[:Org_Assembly];

n=2000
theme1 = Theme(default_color=RGBA(1, 0.5, 0, 0.2))
orgs_cb3 = String["PS312","ce11"]
statistics_df = DataFrame(Org_Assembly = String[], Genome_Size = Int64[], Mean = Float64[], SD= Float64[], Num_above5_SD = Int64[], Num_above10_SD = Int64[], Num_above10 = Int64[], Num_above25 = Int64[], Num_above50 = Int64[])
@parallel for organism_dir in orgs_cb3#HSE_in_out_summary#orgs_with_flagged_heli
    tic()
    print("Starting "*organism_dir*"\n")
    chrom_sizes = CSV.read(mainfolder*organism_dir*"/"*organism_dir*".chrom.sizes", datarow=1, delim = '\t')
    limitless_fimo = CSV.read(mainfolder*organism_dir*"/"*organism_dir*"_extendedFIMO.txt")[4:6]
    open(mainfolder*organism_dir*"/chromosome_scanning/"*organism_dir*".HSE.counts.txt", "w") do io
    end
    genome_size = 0
    for chromosome_ids in unique(limitless_fimo[:sequence_name])
        print(chromosome_ids,"\n")
        chr_subsetter_HSE =  limitless_fimo[limitless_fimo[:sequence_name] .== chromosome_ids,:][:start];
        length_of_chr = chrom_sizes[chrom_sizes[:Column1] .== chromosome_ids,:][:Column2][1];
        genome_size = genome_size+length_of_chr
        open(mainfolder*organism_dir*"/chromosome_scanning/"*organism_dir*".HSE.counts.txt", "a") do io
            writedlm(io, [(count(x->(i<x<(i+n-14)),chr_subsetter_HSE)) for i in 1:n:length_of_chr][1:end-1])
        end
        chr_subsetter_HSE = DataFrame()
    end
    limitless_fimo = DataFrame()
    chrom_sizes =  DataFrame()
    
    max_HSE_value = 0
    m = 0
    standard_dev = 0
    HSE_counts = readdlm(mainfolder*organism_dir*"/chromosome_scanning/"*organism_dir*".HSE.counts.txt");
    try
        max_HSE_value = max(HSE_counts...)
        m = mean(HSE_counts)
        standard_dev = stdm(HSE_counts, m; corrected=true)
    end
    indices_above5SD = find( x->(x > m+5*standard_dev), HSE_counts)
    indices_above10SD = find( x->(x > m+10*standard_dev), HSE_counts)
    indices_above10 = find( x->(x > 10), HSE_counts)
    indices_above25 = find( x->(x > 25), HSE_counts)
    indices_above50 = find( x->(x > 50), HSE_counts)
    push!(statistics_df, [organism_dir,genome_size,m,standard_dev,length(indices_above5SD),length(indices_above10SD),length(indices_above10),length(indices_above25),length(indices_above50)])
	toc()
end
CSV.write(mainfolder*"74_Org_Assemblys.HSE.cluster.stats.csv",statistics_df; delim = ',')
statistics_df = DataFrame()