#!/bin/bash
#PBS -q hotel
#PBS -N ucsc
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00
#PBS -o ucsc.sh.out
#PBS -e ucsc.sh.err
#PBS -V
#PBS -M bvtsu@ucsd.edu

start=$SECONDS

# 1st phase, make main directory if it doesn't already exist
if [ -e $outfolder ]
then
    echo "Main output directory exists."
else
    mkdir $outfolder
    echo "Main output directory created."
fi

# make species-specific directory if it doesn't already exist
if [ -e ${outfolder}/${INDEX}/ ]
then
    echo "Species directory exists."
else
    mkdir ${outfolder}/${INDEX}/
    echo "Species directory created."
fi

# 2nd phase, search HSF-1 binding motif if not already done
if [ -e ${outfolder}/${INDEX}/${INDEX}.motif.tsv -a -e ${outfolder}/${INDEX}/${INDEX}.compiled.chrom.sizes ]
then
    echo "Already processed."

# if not done, check if genome.2bit exists
else
    if [ -e ${outfolder}/${INDEX}/${INDEX}.2bit ]
    then
        echo "ok"

    # if genome.2bit doesn't exist, wget the genome.2bit
    else
        echo "nok"
        wget -w 2 -P ${outfolder}/${INDEX}/ ftp://hgdownload.cse.ucsc.edu/goldenPath/${INDEX}/bigZips/${INDEX}.2bit
    fi
    # Convert .2bit to Fa sequence format
    twobittofa ${outfolder}/${INDEX}/${INDEX}.2bit ${outfolder}/${INDEX}/${INDEX}.fa
    
    # Generate chromsizes for each chromosome -- this is important for matching sites
    python src/chrom_size.py ${outfolder}/${INDEX}/${INDEX}.fa ${outfolder}/${INDEX}/${INDEX}.compiled.chrom.sizes

    # now search for the 15-mer HSF-1 binding motif, which we will now call HSEs, within the genome Fa file
    # Using regex
    python src/motif_search.py ${outfolder}/${INDEX}/${INDEX}.fa

    # Alternatively, can use FIMO here, but need to get background occurrence of nucs in seq
    #

    echo Job exited with status $?
    duration=$(( SECONDS - start ))
    echo "Completed in $duration seconds"
    echo "Cleaning up directory"
    rm ${outfolder}/${INDEX}/${INDEX}.fa
    rm ${outfolder}/${INDEX}/${INDEX}.2bit
fi

# 3rd phase, UCSC already has rmsk.txt files. 
# Download these to side-step using computationally demanding RepeatMasker search
if [ -e ${outfolder}/${INDEX}/rmsk.txt.gz ]
then
    echo "${INDEX} rmsk.txt.gz exists."

# if not done, check if genome.2bit exists
else
    wget -P ${outfolder}/${INDEX}/ ftp://hgdownload.cse.ucsc.edu/goldenPath/${INDEX}/database/rmsk.txt.gz
fi

# decompress
gzip -d ${outfolder}/${INDEX}/rmsk.txt.gz

# 4th phase
python src/pair_hse_helitrons.py -genome ${INDEX} -o ${outfolder}/${INDEX}

# recompress
gzip ${outfolder}/${INDEX}/rmsk.txt