#!/bin/bash

# This script submits jobs to the cluster to run RepeatMasker on each input line from PLANFILE
# Output files are in /home/bvtsu/data/rmsk.out.cgp/

export PLANFILE=/home/bvtsu/subproject/CGP_planfile.txt

while read l; do
    export INDEX=${l}
    echo ${INDEX}
    qsub  -V -N ${INDEX} \
            -o log/rmsk-${INDEX}.out \
            -e log/rmsk-${INDEX}.err \
            /home/bvtsu/subproject/src/rmsk_default_cgp.sbatch
done < $PLANFILE