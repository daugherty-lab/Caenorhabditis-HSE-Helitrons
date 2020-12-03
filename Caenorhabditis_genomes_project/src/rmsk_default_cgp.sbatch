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
    gzip -d /home/bvtsu/ref_genome/CGP/${INDEX}.fa.gz
    start=$SECONDS
    
    fasta-get-markov /home/bvtsu/ref_genome/CGP/${INDEX}.fa /home/bvtsu/ref_genome/CGP/${INDEX}_bg.txt &
    wait $! 
    
    fimo --oc /home/bvtsu/data/fimo_cgp_default/${INDEX}/ --bfile /home/bvtsu/ref_genome/CGP/${INDEX}_bg.txt --verbosity 1 --max-stored-scores 80000000 --motif RGAADNTTCTAGAAR /home/bvtsu/code/helitron.hses.012019/HSF1meme_final.txt /home/bvtsu/ref_genome/CGP/${INDEX}.fa &

    wait $!

    echo Job exited with status $?
    duration=$(( SECONDS - start ))
    echo "Completed in $duration seconds"
    mv /home/bvtsu/data/fimo_cgp_default/${INDEX}/fimo.tsv /home/bvtsu/data/fimo_cgp_default/${INDEX}/${INDEX}_fimo.tsv
fi


#RepeatMasker -s -pa 8 -lib /home/bvtsu/code/helitron.hses/caeHelitrons_final.fa -nolow -norna -gff -no_is -xsmall -dir /home/bvtsu/data/rmsk.out.cgp/ /home/bvtsu/ref_genome/CGP/${INDEX}.fa

gzip /home/bvtsu/ref_genome/CGP/${INDEX}.fa