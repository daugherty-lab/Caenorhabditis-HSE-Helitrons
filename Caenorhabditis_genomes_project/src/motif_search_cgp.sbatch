#!/bin/bash
#PBS -q hotel
#PBS -N "${INDEX}"
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00
#PBS -o motifsearch-"${INDEX}".sh.out
#PBS -e motifsearch-"${INDEX}".sh.err
#PBS -V
#PBS -M bvtsu@ucsd.edu

motif_dataDIR=/home/bvtsu/data/regex_cgp_default/
genomeDIR=/home/bvtsu/genomes/
motifpyDIR=/home/bvtsu/Caenorhabditis-HSE-Helitrons/Caenorhabditis_genomes_project/src/

gzip -d ${genomeDIR}${INDEX}.fa.gz
start=$SECONDS

python ${motifpyDIR}motif_search.py ${genomeDIR}${INDEX}.fa ${motif_dataDIR}
wait $!

echo Job exited with status $?
duration=$(( SECONDS - start ))
echo "Completed in $duration seconds"
gzip ${genomeDIR}${INDEX}.fa
