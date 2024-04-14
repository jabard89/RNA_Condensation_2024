#!/bin/bash
DIR=.
curl -C - -o $DIR/S288C_reference_genome_R64-3-1_20210421.tgz http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-3-1_20210421.tgz
tar -xvzf $DIR/S288C_reference_genome_R64-3-1_20210421.tgz

gzip -dfc $DIR/S288C_reference_genome_R64-3-1_20210421/saccharomyces_cerevisiae_R64-3-1_20210421.gff.gz > $DIR/saccharomyces_cerevisiae_R64-3-1_20210421.gff
gzip -dfc $DIR/S288C_reference_genome_R64-3-1_20210421/orf_coding_all_R64-3-1_20210421.fasta.gz > $DIR/orf_coding_all_R64-3-1_20210421.fasta
gzip -dfc $DIR/S288C_reference_genome_R64-3-1_20210421/rna_coding_R64-3-1_20210421.fasta.gz > $DIR/rna_coding_R64-3-1_20210421.fasta

curl -C - -o $DIR/Schizosaccharomyces_pombe_all_chromosomes.fa.gz https://www.pombase.org/data/genome_sequence_and_features/genome_sequence/Schizosaccharomyces_pombe_all_chromosomes.fa.gz
gzip -dfc $DIR/Schizosaccharomyces_pombe_all_chromosomes.fa.gz > $DIR/Schizosaccharomyces_pombe_all_chromosomes.fa

curl -C - -o $DIR/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz https://www.pombase.org/data/genome_sequence_and_features/gff3/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz
gzip -dfc $DIR/Schizosaccharomyces_pombe_all_chromosomes.gff3.gz > $DIR/Schizosaccharomyces_pombe_all_chromosomes.gff3

curl -C - -o $DIR/Schizosaccharomyces_pombe_cds.fa.gz https://www.pombase.org/data/genome_sequence_and_features/feature_sequences/cds.fa.gz
gzip -dfc $DIR/Schizosaccharomyces_pombe_cds.fa.gz > $DIR/Schizosaccharomyces_pombe_cds.fa

bash relabel_chroms_gff3.sh $DIR/Schizosaccharomyces_pombe_all_chromosomes.gff3 > $DIR/Schizosaccharomyces_pombe_all_chromosomes_relabelled.gff3
bash relabel_fasta_pombe.sh $DIR/Schizosaccharomyces_pombe_all_chromosomes.fa > $DIR/Schizosaccharomyces_pombe_all_chromosomes_relabelled.fa
bash split_SGD_gff3.sh

python gffutils_add_geneid.py $DIR/Schizosaccharomyces_pombe_all_chromosomes_relabelled.gff3  $DIR/Schizosaccharomyces_pombe_all_chromosomes_relabelled_geneid.gff3
python gffutils_add_geneid.py $DIR/saccharomyces_cerevisiae_R64-3-1_20210421_nofasta.gff  $DIR/saccharomyces_cerevisiae_R64-3-1_20210421_nofasta_geneid.gff

cat $DIR/saccharomyces_cerevisiae_R64-3-1_20210421_nofasta_geneid.gff $DIR/Schizosaccharomyces_pombe_all_chromosomes_relabelled_geneid.gff3 > \
$DIR/saccharomyces_cerevisiae_R64-3-1_20210421_Schizosaccharomyces_pombe_all_chromosomes_relabelled_geneid.gff3
cat $DIR/saccharomyces_cerevisiae_R64-3-1_20210421_allchrom.fasta $DIR/Schizosaccharomyces_pombe_all_chromosomes_relabelled.fa > \
$DIR/saccharomyces_cerevisiae_R64-3-1_20210421_Schizosaccharomyces_pombe_all_chromosomes_relabelled.fasta
cat $DIR/orf_coding_all_R64-3-1_20210421.fasta $DIR/rna_coding_R64-3-1_20210421.fasta $DIR/Schizosaccharomyces_pombe_cds.fa > \
$DIR/Scerevisiae_orf_coding_all_Scerevisiae_rna_coding_Spombe_cds.fasta

rm -rf $DIR/S288C*
rm -f $DIR/Schizosaccharomyces_pombe_*
rm -f $DIR/saccharomyces_cerevisiae_R64-3-1_20210421.gff
rm -f $DIR/saccharomyces_cerevisiae_R64-3-1_20210421_allchrom.fasta
rm -f $DIR/saccharomyces_cerevisiae_R64-3-1_20210421_nofasta.gff
rm -f $DIR/orf_coding_all_R64-3-1_20210421.fasta
rm -f $DIR/rna_coding_R64-3-1_20210421.fasta
