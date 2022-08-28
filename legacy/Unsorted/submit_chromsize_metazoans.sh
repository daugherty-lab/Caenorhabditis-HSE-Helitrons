#!/bin/bash

# This script submits jobs to the cluster to run RepeatMasker on each input line from PLANFILE

#export PLANFILE=/home/bvtsu/code/helitron.hses/refpiece1.txt 
export PLANFILE=refpiece1.txt

sed 1d $PLANFILE |
while read l; do
    export INDEX="${l}"
    echo ${INDEX}
    wget_url=$(grep "${INDEX}" assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
    genome_accession=${wget_url##*/}
    if [ -e refseq_out/"${INDEX//[[:blank:]]/}".compiled.chrom.sizes ]
    then
        echo "Already processed."
    else
        echo "nok"
        ./chrom_size_refseq.sbatch || echo "missing" &
        wait $!
    fi
done
