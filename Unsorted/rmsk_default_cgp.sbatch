#!/bin/bash
#PBS -q hotel
#PBS -N rmsk
#PBS -l nodes=1:ppn=8
#PBS -l walltime=4:00:00
#PBS -o rmsk.sh.out
#PBS -e rmsk.sh.err
#PBS -V
#PBS -M bvtsu@ucsd.edu

motif_dataDIR=/home/bvtsu/data/fimo_cgp_default/
rmsk_dataDIR=/home/bvtsu/data/rmsk.out.cgp/
genomeDIR=/home/bvtsu/genomes/
caeHelitronDIR=/home/bvtsu/Caenorhabditis-HSE-Helitrons/Caenorhabditis_genomes_project/

if [ -e ${motif_dataDIR}${INDEX}_fimo.tsv ]
then
    echo "Already processed."
else
    gzip -d ${genomeDIR}${INDEX}.fa.gz
    start=$SECONDS
    
    fasta-get-markov ${genomeDIR}${INDEX}.fa ${genomeDIR}${INDEX}_bg.txt &
    wait $! 
    
    fimo --oc ${motif_dataDIR}${INDEX}/ --bfile ${genomeDIR}${INDEX}_bg.txt --verbosity 1 --max-stored-scores 80000000 --motif RGAADNTTCTAGAAR ${caeHelitronDIR}HSF1meme_final.txt ${genomeDIR}${INDEX}.fa &

    wait $!

    echo Job exited with status $?
    duration=$(( SECONDS - start ))
    echo "Completed in $duration seconds"
    mkdir ${motif_dataDIR}${INDEX}/ 
    mv ${motif_dataDIR}${INDEX}/fimo.tsv ${motif_dataDIR}${INDEX}/${INDEX}_fimo.tsv
fi


RepeatMasker -s -pa 8 -lib ${caeHelitronDIR}caeHelitrons_final.fa -nolow -norna -gff -no_is -xsmall -dir ${rmsk_dataDIR} ${genomeDIR}${INDEX}.fa

gzip ${genomeDIR}${INDEX}.fa