#!/bin/bash

# This script submits jobs to the cluster to run RepeatMasker on each input line from PLANFILE
# Output files are in /home/bvtsu/data/rmsk.out.cgp/

export PLANFILE=example/CGP_planfile.txt #Note, there must be a newline at the end for this to read each line

while read l; do #start a while loop to read the line (l) as a string
    export INDEX=${l} #store string into variable INDEX
    echo ${INDEX} #print the stored string
    #queues the affiliated command to TSCC, make log for stdout and stderr
    qsub  -V -N ${INDEX} \
            -o log/rmsk-${INDEX}.out \
            -e log/rmsk-${INDEX}.err \
            src/rmsk_default_cgp.sbatch #call src script to generate fimo.tsv (HSE motifs) and rmsk files (helitrons)
done < $PLANFILE  #This calls the defined planfile above to be used for the entire while loop