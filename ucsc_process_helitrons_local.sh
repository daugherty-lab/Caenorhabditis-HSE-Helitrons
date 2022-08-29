#!/bin/bash

# This script submits jobs to the cluster to download 
# UCSC genomes using PLANFILE lines as input.
# Output files are in:

# "$#" = value of the total number of command line args passed
# As long as "$#" is greater than (-gt) 0 args, keep while loop alive
while [[ "$#" -gt 0 ]]; do
    # Check each case (options/flags) until match is found
    case $1 in
        # get input following arg option, then shift to next str
        -i|--inputlist) inputlist="$2"; shift ;;
        -o|--oc) outputdir="$2"; shift ;;
        
        # if extra, unmatched options show up, exit
        # Exit code 0 - Success
        # Exit code 1 - General/misc errors, such as "divide by zero" and other impermissible operations
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    
    # end case search (case spelled backwards)
    esac
    shift # to the next str, if any, then loop
done

PLANFILE=$inputlist


while read line; # start a while loop to read the line (l) as a string
do 
    export INDEX=$line # store string into variable INDEX
    export outfolder=$outputdir
    echo $INDEX # print the line
    echo $outfolder # print the output destination
    src/./ucsc_local.sbatch &
done < $PLANFILE  #This calls the defined planfile above to be used for the entire while loop