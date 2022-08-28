#!/bin/bash
wget_url=$(grep "${INDEX}" assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
genome_accession=${wget_url##*/}
if [ -e refseq_out/"${INDEX//[[:blank:]]/}".compiled.chrom.sizes ]
then
    echo "Already processed."
else
    #touch refseq_out/"${INDEX//[[:blank:]]/}".compiled.chrom.sizes
    if [ -e sequence/${genome_accession}_genomic.fna.gz ]
    then
        echo "ok"
    else
        echo "nok"
        wget -w 2 -P sequence/ ${wget_url}/${genome_accession}_genomic.fna.gz
        gzip -d sequence/${genome_accession}_genomic.fna.gz
    fi
    start=$SECONDS
    python chrom_size.py sequence/${genome_accession}_genomic.fna refseq_out/"${INDEX//[[:blank:]]/}".compiled.chrom.sizes
    wait $!

    echo Job exited with status $?
    duration=$(( SECONDS - start ))
    echo "Completed in $duration seconds"
    rm sequence/${genome_accession}_genomic.fna
fi
