#!/bin/bash

#USAGE
#./ IterateFasttree.sh

#read all files matching pattern .fasta in the current folder
files=(*.fasta)
echo "identified ${#files[@]} alignments"

#go through each alignment to build a fasttree and export as same file name with .tre
for i in "${files[@]}"; do
    echo "making a tree for $i"
    echo "$i"
    outfile=$(echo "$i" | sed 's,.fasta,.tre,g')
    Fasttree -gtr -nt -fastest "$i" > "$outfile"
done
