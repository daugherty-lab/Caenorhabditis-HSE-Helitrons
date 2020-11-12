#!/bin/bash

# This script submits jobs to the cluster to run RepeatMasker on each input line from PLANFILE

#export PLANFILE=/home/bvtsu/code/helitron.hses/refseq_metazoans_planfile.txt 
export PLANFILE=refseq_metazoans_planfile.txt

sed 1d $PLANFILE |
while read l; do
    export INDEX="${l}"
    echo ${INDEX}
    wget_url=$(grep "${INDEX}" assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
    genome_accession=${wget_url##*/}
    if [ -e refseq_out/"${INDEX//[[:blank:]]/}"_extendedFIMO.txt ]
    then
        echo "already processed"
    else
        echo "nok"
        ./hse_helitron.sbatch || echo "missing" &
        wait $!
    fi
done
