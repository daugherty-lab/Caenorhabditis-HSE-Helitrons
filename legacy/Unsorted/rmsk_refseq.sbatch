#!/bin/bash
#PBS -q hotel
#PBS -N "${INDEX}"
#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -o rmsk-"${INDEX}".sh.out
#PBS -e rmsk-"${INDEX}".sh.err
#PBS -V
#PBS -M bvtsu@ucsd.edu

wget_url=$(grep "${INDEX}" /home/bvtsu/ref_genome/assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
genome_accession=${wget_url##*/}
if [ -e /home/bvtsu/data/rmsk.out.refseq/${genome_accession}_genomic.fna.ori.out ]
then
    echo "Already processed."
else
    if [ -e /home/bvtsu/ref_genome/refseq/${genome_accession}_genomic.fna.gz ]
    then
        echo "ok"
    else
        echo "nok"
        wget -w 2 -P /home/bvtsu/ref_genome/refseq/ ${wget_url}/${genome_accession}_genomic.fna.gz
        gzip -d /home/bvtsu/ref_genome/refseq/${genome_accession}_genomic.fna.gz
    fi
    start=$SECONDS
    RepeatMasker -s -pa 8 -lib /home/bvtsu/code/helitron.hses/caeHelitrons_final.fa -nolow -norna -gff -no_is -xsmall -dir /home/bvtsu/data/rmsk.out.refseq/ /home/bvtsu/ref_genome/refseq/${genome_accession}_genomic.fna &

    wait $!

    echo Job exited with status $?
    duration=$(( SECONDS - start ))
    echo "Completed in $duration seconds"
    rm /home/bvtsu/data/rmsk.out.refseq/${genome_accession}_genomic.fna.masked
    rm /home/bvtsu/ref_genome/refseq/${genome_accession}_genomic.fna
fi
