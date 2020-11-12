#!/bin/bash
#PBS -q hotel
#PBS -N "${INDEX}"
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00
#PBS -o hsepileup-"${INDEX}".sh.out
#PBS -e hsepileup-"${INDEX}".sh.err
#PBS -V
#PBS -M bvtsu@ucsd.edu

start=$SECONDS
python /home/bvtsu/code/helitron.hses/src/hsepileup.py "${INDEX//[[:blank:]]/}"
wait $!

echo Job exited with status $?
duration=$(( SECONDS - start ))
echo "Completed in $duration seconds"
