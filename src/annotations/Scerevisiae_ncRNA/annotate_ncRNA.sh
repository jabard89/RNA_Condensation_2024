#!/bin/bash
curl -C - -o SGD_230511_rna_coding.fasta.gz http://sgd-archive.yeastgenome.org/sequence/S288C_reference/rna/rna_coding.fasta.gz
curl -C - -o 230511_SGD_ncRNA_xref.txt http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_ncRNA_xref.txt
gzip -dfc SGD_230511_rna_coding.fasta.gz > SGD_230511_rna_coding.fasta
bash extract_fasta_lengths.sh *rna_coding.fasta > SGD_230511_rna_coding_lengths.tsv