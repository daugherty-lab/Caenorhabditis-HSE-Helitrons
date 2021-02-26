#!/bin/bash

# This script submits jobs to the cluster to run pair_hse_helitrons.py on each input line from PLANFILE

#Assign your planfile to the PLANFILE variable
export PLANFILE=example/CGP_planfile.txt #planfile exists in example folder, but there must be a newline at the end for this to read each line
outfolder=/home/bvtsu/data/CGP_out/

while read l; do #start a while loop to read the line (l) as a string
    export INDEX=${l} #store string into variable INDEX
    echo ${INDEX} #print the stored string
    #wget_url=$(grep "${INDEX}" assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
    #genome_accession=${wget_url##*/}
    if [ -e ${outfolder}${INDEX}_extendedFIMO.txt ] #check specified data folder location for existing outputs
    then
    #    echo "already processed" #don't need to re-create existing data
    #else
        echo "nok"
        #queues the affiliated command to TSCC, make log for stdout and stderr
        qsub  -V -N ${INDEX} \
                -o log/hse_heli-${INDEX}.out \
                -e log/hse_heli-${INDEX}.err \
                src/hse_helitron.sbatch #call src script to analyze hse-helitron pairs
    fi
done < $PLANFILE #This calls the defined planfile above to be used for the entire while loop