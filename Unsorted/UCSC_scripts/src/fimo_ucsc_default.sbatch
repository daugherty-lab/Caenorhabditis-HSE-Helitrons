#!/bin/bash
#PBS -q hotel
#PBS -N fimo
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00
#PBS -o fimo.sh.out
#PBS -e fimo.sh.err
#PBS -V
#PBS -M bvtsu@ucsd.edu

#1st phase, make main directory if it doesn't already exist
if [ -e /home/bvtsu/data/UCSC/ ]
then
    echo "Directory exists."
else
    mkdir /home/bvtsu/data/UCSC/
    echo "Directory created."
fi

#make species-specific directory if it doesn't already exist
if [ -e /home/bvtsu/data/UCSC/${INDEX}/ ]
then
    echo "Species directory exists."
else
    mkdir /home/bvtsu/data/UCSC/${INDEX}/
    echo "Species directory created."
fi

#2nd phase, UCSC already has rmsk.txt files. Download these to side-step using computationally demanding RepeatMasker search
wget -P /home/bvtsu/data/UCSC/${INDEX}/ ftp://hgdownload.cse.ucsc.edu/goldenPath/${INDEX}/database/rmsk.txt.gz

#3rd phase, check if HSF-1 binding motif search was already done
if [ -e /home/bvtsu/data/UCSC/${INDEX}_fimo.tsv ]
then
    echo "Already processed."
#if not done, check if genome.2bit exists
else
    if [ -e /home/bvtsu/ref_genome/ucsc_genomes/${INDEX}.2bit ]
    then
        echo "ok"
#if genome.2bit doesn't exist, wget the genome.2bit
    else
        echo "nok"
        wget -w 2 -P /home/bvtsu/ref_genome/ucsc_genomes/ ftp://hgdownload.cse.ucsc.edu/goldenPath/${INDEX}/bigZips/${INDEX}.2bit
    fi
#Convert .2bit to Fa sequence format and get background presence of nucleotides across the genome
    /home/bvtsu/./twoBitToFa /home/bvtsu/ref_genome/ucsc_genomes/${INDEX}.2bit /home/bvtsu/ref_genome/ucsc_genomes/${INDEX}.fa
    start=$SECONDS
    fasta-get-markov /home/bvtsu/ref_genome/ucsc_genomes/${INDEX}.fa /home/bvtsu/ref_genome/ucsc_genomes/${INDEX}_bg.txt &
    wait $!

#now search for the 15-mer HSF-1 binding motif, which we will now call HSEs, within the genome Fa file
#Using the MEME file
    fimo --oc /home/bvtsu/data/UCSC/${INDEX}/ --bfile /home/bvtsu/ref_genome/ucsc_genomes/${INDEX}_bg.txt --verbosity 1 --max-stored-scores 80000000 --motif RGAADNTTCTAGAAR /home/bvtsu/code/helitron.hses.012019/HSF1meme_final.txt /home/bvtsu/ref_genome/ucsc_genomes/${INDEX}.fa &

#Using regex
#
    wait $!

    echo Job exited with status $?
    duration=$(( SECONDS - start ))
    echo "Completed in $duration seconds"
    rm /home/bvtsu/ref_genome/ucsc_genomes/${INDEX}.fa
    rm /home/bvtsu/ref_genome/ucsc_genomes/${INDEX}.2bit
fi