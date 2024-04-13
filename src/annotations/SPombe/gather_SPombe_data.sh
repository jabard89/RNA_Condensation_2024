#!/bin/bash

curl -C - -o Schizosaccharomyces_pombe_all_chromosomes.fa.gz https://www.pombase.org/data/genome_sequence_and_features/genome_sequence/Schizosaccharomyces_pombe_all_chromosomes.fa.gz
gzip -dfc Schizosaccharomyces_pombe_all_chromosomes.fa.gz > Schizosaccharomyces_pombe_all_chromosomes.fa

curl -C - -o Schizosaccharomyces_pombe_all_chromosomes.gff3.gz https://www.pombase.org/data/genome_sequence_and_features/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz
gzip -dfc Schizosaccharomyces_pombe_all_chromosomes.gff3.gz > Schizosaccharomyces_pombe_all_chromosomes.gff3

curl -C - -o Schizosaccharomyces_pombe_cds.fa.gz https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/cds.fa.gz
gzip -dfc Schizosaccharomyces_pombe_cds.fa.gz > Schizosaccharomyces_pombe_cds.fa

curl -C - -o gene_IDs_names_products.tsvz https://www.pombase.org/data/names_and_identifiers/gene_IDs_names_products.tsv
gzip -dfc gene_IDs_names_products.tsvz > Schizosaccharomyces_pombe_230503_gene_IDs_names_products.tsv

curl -C - -o 5UTR.fa.gz https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/UTR/5UTR.fa.gz
gzip -dfc 5UTR.fa.gz > Schizosaccharomyces_pombe_230503_5UTR.fa
curl -C - -o 3UTR.fa.gz https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/UTR/3UTR.fa.gz
gzip -dfc 3UTR.fa.gz > Schizosaccharomyces_pombe_230503_3UTR.fa


#python gffutils_add_geneid.py Schizosaccharomyces_pombe_all_chromosomes.gff3  Schizosaccharomyces_pombe_all_chromosomes_geneid.gff3
#python count_fasta_length.py Schizosaccharomyces_pombe_cds.fa  Schizosaccharomyces_pombe_cds_length.tsv
python count_fasta_length.py Schizosaccharomyces_pombe_230503_5UTR.fa  Schizosaccharomyces_pombe_230503_5UTR_length.tsv
python count_fasta_length.py Schizosaccharomyces_pombe_230503_3UTR.fa  Schizosaccharomyces_pombe_230503_3UTR_length.tsv


