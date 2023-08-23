#!/bin/bash

#USAGE
#copy the script into the folder
#./IterateIQtree.sh

#read all files matching pattern .fasta in the current folder
alignfiles=(*.fasta)
echo "identified ${#alignfiles[@]} alignments"

#go through each alignment, identify tree and date file corresponding and run all through iqtree (NOTE needs lots of mem)
for i in "${alignfiles[@]}"; do
    echo "initiating iqtree for $i"
    treefile=$(echo "$i" | sed 's/.fasta/.tre/g')
    datefile=$(echo "$i" | sed 's/.fasta/_dates.txt/g')
    iqtree2 -s "$i" --date "$datefile" -te "$treefile" --dating LSD --clock-sd 0.2 -nt 64 -mem 70G -m GTR+I+R3 --date-ci 50 -o "EPI_ISL_402125/2019-12-31"

    
done
