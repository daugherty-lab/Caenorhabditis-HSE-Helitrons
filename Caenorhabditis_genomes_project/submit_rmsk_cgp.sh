#!/bin/bash

# This script submits jobs to the cluster to run RepeatMasker on each input line from PLANFILE
# Output files are in /home/bvtsu/data/rmsk.out

export PLANFILE=/home/bvtsu/code/helitron.hses/CGP_planfile.txt 

sed 1d $PLANFILE |
while read l; do
    export INDEX=${l}
    echo ${INDEX}
    qsub  -V -N ${INDEX} \
            -o log/rmsk-${INDEX}.out \
            -e log/rmsk-${INDEX}.err \
            /home/bvtsu/code/helitron.hses/src/rmsk_default_cgp.sbatch
done