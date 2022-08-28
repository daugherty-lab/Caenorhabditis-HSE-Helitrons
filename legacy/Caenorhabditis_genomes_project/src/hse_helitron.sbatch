#!/bin/bash
#PBS -q hotel
#PBS -N rmsk
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00
#PBS -o hse-heli.sh.out
#PBS -e hse-heli.sh.err
#PBS -V
#PBS -M bvtsu@ucsd.edu

infolder=/home/bvtsu/data/
outfolder=/home/bvtsu/data/CGP_out/
pypath=/home/bvtsu/Caenorhabditis-HSE-Helitrons/Caenorhabditis_genomes_project/src/
rmskDIR=/home/bvtsu/data/rmsk.out.cgp/download.caenorhabditis.org/v1/repeatmasker/
###What is this doing here? wget_url=$(grep "${INDEX}" assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
#genome_accession=${wget_url##*/}
#echo ${INDEX//[[:blank:]]/}
echo ${INDEX}
#echo ${genome_accession} 

gzip -d ${rmskDIR}${INDEX}.fa.repeatmasker.out.gz
gzip -d ${rmskDIR}${INDEX}.fa.repeatmasker.tbl.gz

python ${pypath}pair_hse_helitrons.py ${infolder} ${INDEX} ${outfolder} &
wait $!
echo Job exited with status $?
gzip ${rmskDIR}${INDEX}.fa.repeatmasker.out
gzip ${rmskDIR}${INDEX}.fa.repeatmasker.tbl