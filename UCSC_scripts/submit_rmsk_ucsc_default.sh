#!/bin/bash

# This script submits jobs to the cluster to download rmsk files, .2bit files from UCSC genome browser for subsequent HSE motif searching using species assembly names in each input line from PLANFILE
# Output files are in /home/bvtsu/data/UCSC/
#num="11"
#export PLANFILE=/home/bvtsu/code/helitron.hses.012019/planfile_rmsk.split.txt00${num}
export PLANFILE=/home/bvtsu/code/helitron.hses.012019/UCSC_planfile.txt

sed 1d $PLANFILE |
while read l; do
    export INDEX=${l}
    echo ${INDEX}

    qsub  -V -N ${INDEX} \
            -o log/fimo-${INDEX}.out \
            -e log/fimo-${INDEX}.err \
            /home/bvtsu/code/helitron.hses.012019/src/fimo_ucsc_default.sbatch
done