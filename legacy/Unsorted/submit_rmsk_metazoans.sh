#!/bin/bash

# This script submits jobs to the cluster to run RepeatMasker on each input line from PLANFILE

export PLANFILE=/home/bvtsu/code/helitron.hses/refseq_metazoans_planfile.txt 

sed 1d $PLANFILE |
while read l; do
    export INDEX="${l}"
    echo ${INDEX}
    wget_url=$(grep "${INDEX}" /home/bvtsu/ref_genome/assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
    genome_accession=${wget_url##*/}
    if [ -e /home/bvtsu/data/rmsk.out.refseq/${genome_accession}_genomic.fna.ori.out ]
    then
        echo "Already processed."
    else
        echo "nok"
        qsub  -V -N ${INDEX//[[:blank:]]/} \
            -o log/rmsk-"${INDEX//[[:blank:]]/}".out \
            -e log/rmsk-"${INDEX//[[:blank:]]/}".err \
            /home/bvtsu/code/helitron.hses/src/rmsk_refseq.sbatch
    fi
done
