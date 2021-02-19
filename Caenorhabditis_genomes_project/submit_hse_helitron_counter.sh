#!/bin/bash

# This script submits jobs to the cluster to run pair_hse_helitrons.py on each input line from PLANFILE

#Assign your planfile to the PLANFILE variable
export PLANFILE=example/CGP_planfile.txt #planfile exists in example folder

while read l; do #start a while loop to read the line (l) as a string
    export INDEX="${l}" #store string into variable INDEX
    echo ${INDEX} #print the stored string
    #wget_url=$(grep "${INDEX}" assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
    #genome_accession=${wget_url##*/}
    if [ -e ../../data/fimo_cgp_default/"${INDEX//[[:blank:]]/}"_extendedFIMO.txt ] #check specified data folder location for existing outputs
    then
        echo "already processed" #don't need to re-create existing data
    else
        echo "nok"
        qsub  -V -N ${INDEX//[[:blank:]]/} \ #qsub queues the affiliated command to TSCC
            -o log/hsehelitron-"${INDEX//[[:blank:]]/}".out \ #make log for terminal stdout
            -e log/hsehelitron-"${INDEX//[[:blank:]]/}".err \ #make log for terminal stderr
            src/hse_helitron.sbatch #call src script to analyze hse-helitron pairs
    fi
done < $PLANFILE #This calls the defined planfile above to be used for the entire while loop