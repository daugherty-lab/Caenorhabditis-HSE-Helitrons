#!/bin/bash
#PBS -q hotel
#PBS -N "${INDEX}"
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00
#PBS -o motifsearch-"${INDEX}".sh.out
#PBS -e motifsearch-"${INDEX}".sh.err
#PBS -V
#PBS -M bvtsu@ucsd.edu

wget_url=$(grep "${INDEX}" /home/bvtsu/ref_genome/assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
genome_accession=${wget_url##*/}
if [ -e /home/bvtsu/data/refseq_motif_files/${genome_accession}_genomic.fna.motif.tsv ]
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
    python /home/bvtsu/code/helitron.hses/src/motif_search.py /home/bvtsu/ref_genome/refseq/${genome_accession}_genomic.fna
    wait $!

    echo Job exited with status $?
    duration=$(( SECONDS - start ))
    echo "Completed in $duration seconds"
    rm /home/bvtsu/ref_genome/refseq/${genome_accession}_genomic.fna
fi
