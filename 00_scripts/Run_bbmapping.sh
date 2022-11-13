#!/bin/bash

infile=""  # 09b_Extracted flanking genome seqs.
outbase="" # Filename base for outfiles
genome1="" # path to genome.fasta or genome.fasta.gz
genome2="" 
maxindel="50" # set higher e.g. 400 if using denovo stacks

echo "bbmapping to a bunch of genomes\n"
echo "1) Genome 1\n"

bbmap.sh ordered=t maxindel=$maxindel build=1 ref=$genome1 in=$infile out="$outbase-all.sam" outm="$outbase-matched.sam" 

echo "\n2) Genome 2\n"

bbmap.sh ordered=t maxindel=$maxindel build=1 ref=$genome1 in=$infile out="$outbase-2-all.sam" outm="$outbase-2-matched.sam" 

# etcetera
