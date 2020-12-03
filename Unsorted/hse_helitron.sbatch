#!/bin/bash
wget_url=$(grep "${INDEX}" assembly_summary_refseq.txt | grep -v virus | sort -rnk13 | head -n 1 | cut -f20)
genome_accession=${wget_url##*/}
echo ${INDEX//[[:blank:]]/}
echo ${genome_accession} 
python pair_hse_helitrons.py ${genome_accession} ${INDEX//[[:blank:]]/} &
wait $!
echo Job exited with status $?
