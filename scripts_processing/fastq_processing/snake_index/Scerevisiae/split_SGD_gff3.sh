#!/bin/bash

# Set input file name
input_file="saccharomyces_cerevisiae_R64-3-1_20210421.gff"

# Set output file names
gff_file="saccharomyces_cerevisiae_R64-3-1_20210421_nofasta.gff"
fasta_file="saccharomyces_cerevisiae_R64-3-1_20210421_allchrom.fasta"

# Use awk to split the file based on the ##FASTA delimiter
awk '/^##FASTA/{flag=1; next} {print > (flag ? "'"$fasta_file"'" : "'"$gff_file"'")}' $input_file
