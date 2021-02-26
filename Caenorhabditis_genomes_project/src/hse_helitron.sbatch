#!/bin/bash
#PBS -q hotel
#PBS -N rmsk
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00
#PBS -o hse-heli.sh.out
#PBS -e hse-heli.sh.err
#PBS -V
#PBS -M bvtsu@ucsd.edu

infolder=/home/bvtsu/data/
outfolder=/home/bvtsu/data/CGP_out/

###What is this doing here? wget_url=$(grep "${INDEX}" assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
#genome_accession=${wget_url##*/}
#echo ${INDEX//[[:blank:]]/}
echo ${INDEX}
#echo ${genome_accession} 
if [ -e ${outfolder} ] #check specified data folder location for existing outputs
then
    echo "outfolder exists" #don't need to re-create existing data
else
    echo "nok"
    mkdir ${outfolder}
fi
python /home/bvtsu/Caenorhabditis-HSE-Helitrons/Caenorhabditis_genomes_project/src/pair_hse_helitrons.py ${infolder} ${INDEX} ${outfolder} &
wait $!
echo Job exited with status $?
