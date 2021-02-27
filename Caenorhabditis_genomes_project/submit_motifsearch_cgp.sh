#!/bin/bash

# This script submits jobs to the cluster to run motif_search.py on each input line from PLANFILE
#export PLANFILE=example/CGP_planfile.txt #Note, there must be a newline at the end for this to read each line
export PLANFILE=example/CGP_planfile_full.txt
output_folder=/home/bvtsu/data/regex_cgp_default/

if [ -e ${output_folder} ] #check specified data folder location for existing outputs
then
    echo "outfolder exists" #don't need to re-create existing data
else
    echo "nok"
    mkdir ${output_folder}
fi

while read l; do #start a while loop to read the line (l) as a string
    export INDEX=${l} #store string into variable INDEX
    echo ${INDEX} #print the stored string
    if [ -e ${output_folder}${INDEX}.fa.motif.tsv ]
    then
        echo "Already processed."
    else #queues the affiliated command to TSCC, make log for stdout and stderr
        echo "Identifying motifs"
        qsub  -V -N ${INDEX} \
            -o log/motifsearch-${INDEX}.out \
            -e log/motifsearch-${INDEX}.err \
            src/motif_search_cgp.sbatch
    fi
done < ${PLANFILE}
