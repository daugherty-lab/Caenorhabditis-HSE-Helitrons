#!/bin/bash
#PBS -q hotel
#PBS -N rmsk
#PBS -l nodes=1:ppn=8
#PBS -l walltime=4:00:00
#PBS -o rmsk.sh.out
#PBS -e rmsk.sh.err
#PBS -V
#PBS -M bvtsu@ucsd.edu


if [ -e /home/bvtsu/data/fimo_cgp_default/${INDEX}_fimo.tsv ]
then
    echo "Already processed."
else
    gzip -d /home/bvtsu/genomes/${INDEX}.fa.gz
    start=$SECONDS
    
    fasta-get-markov /home/bvtsu/genomes/${INDEX}.fa /home/bvtsu/genomes/${INDEX}_bg.txt &
    wait $! 
    
    fimo --oc /home/bvtsu/data/fimo_cgp_default/${INDEX}/ --bfile /home/bvtsu/genomes/${INDEX}_bg.txt --verbosity 1 --max-stored-scores 80000000 --motif RGAADNTTCTAGAAR /home/bvtsu/subproject/HSF1meme_final.txt /home/bvtsu/genomes/${INDEX}.fa &

    wait $!

    echo Job exited with status $?
    duration=$(( SECONDS - start ))
    echo "Completed in $duration seconds"
    mv /home/bvtsu/data/fimo_cgp_default/${INDEX}/fimo.tsv /home/bvtsu/data/fimo_cgp_default/${INDEX}/${INDEX}_fimo.tsv
fi


RepeatMasker -s -pa 8 -lib /home/bvtsu/subproject/caeHelitrons_final.fa -nolow -norna -gff -no_is -xsmall -dir /home/bvtsu/data/rmsk.out.cgp/ /home/bvtsu/genomes/${INDEX}.fa

gzip /home/bvtsu/genomes/${INDEX}.fa