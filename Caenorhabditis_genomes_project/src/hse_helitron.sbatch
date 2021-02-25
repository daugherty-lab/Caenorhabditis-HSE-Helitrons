#!/bin/bash
#PBS -q hotel
#PBS -N rmsk
#PBS -l nodes=1:ppn=8
#PBS -l walltime=4:00:00
#PBS -V
#PBS -M bvtsu@ucsd.edu

infolder=/home/bvtsu/data/
outfolder=/home/bvtsu/data/CGP_out/

###What is this doing here? wget_url=$(grep "${INDEX}" assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
#genome_accession=${wget_url##*/}
#echo ${INDEX//[[:blank:]]/}
echo ${INDEX}
#echo ${genome_accession} 
python pair_hse_helitrons.py ${infolder} ${INDEX} ${outfolder} &
wait $!
echo Job exited with status $?
