#!/bin/bash

# This script submits jobs to the cluster to run RepeatMasker on each input line from PLANFILE

export PLANFILE=/home/bvtsu/code/helitron.hses/refseq_metazoans_planfile.txt 

sed 1d $PLANFILE |
while read l; do
    export INDEX="${l}"
    echo ${INDEX}
    wget_url=$(grep "${INDEX}" /home/bvtsu/ref_genome/assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
    genome_accession=${wget_url##*/}
    if [ -e /home/bvtsu/data/refseq_motif_files/${genome_accession}_genomic.fna.motif.tsv ]
    then
        echo "Already processed."
    else
        echo "nok"
        qsub  -V -N ${INDEX//[[:blank:]]/} \
            -o log/motifsearch-"${INDEX//[[:blank:]]/}".out \
            -e log/motifsearch-"${INDEX//[[:blank:]]/}".err \
            /home/bvtsu/code/helitron.hses/src/motif_search_metazoans.sbatch
    fi
done
