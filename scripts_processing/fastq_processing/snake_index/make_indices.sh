#!/bin/bash
# genome=Scerevisiae_spike/spike_saccharomyces_cerevisiae_R64-3-1_20210421.fasta
# gtf=Scerevisiae_spike/spike_saccharomyces_cerevisiae_R64-3-1_20210421_geneid.gff3
# transcriptome=Scerevisiae_spike/spike_Scerevisiae_orf_coding_all_Scerevisiae_rna_coding.fasta
# snakemake --config genome=${genome} \
# gtf=${gtf} \
# transcriptome=${transcriptome} \
# --cores 1 -p --use-conda


# genome=Scerevisiae_Spombe_spike/spike_saccharomyces_cerevisiae_R64-3-1_20210421_Schizosaccharomyces_pombe_all_chromosomes_relabelled.fasta
# gtf=Scerevisiae_Spombe_spike/spike_saccharomyces_cerevisiae_R64-3-1_20210421_Schizosaccharomyces_pombe_all_chromosomes_relabelled_geneid.gff3
# transcriptome=Scerevisiae_Spombe_spike/spike_Scerevisiae_orf_coding_all_Scerevisiae_rna_coding_Spombe_cds.fasta
# snakemake --config genome=${genome} \
# gtf=${gtf} \
# transcriptome=${transcriptome} \
# --cores 1 -p --use-conda

# genome=Scerevisiae_Spombe/saccharomyces_cerevisiae_R64-3-1_20210421_Schizosaccharomyces_pombe_all_chromosomes_relabelled.fasta
# gtf=Scerevisiae_Spombe/saccharomyces_cerevisiae_R64-3-1_20210421_Schizosaccharomyces_pombe_all_chromosomes_relabelled_geneid.gff3
# transcriptome=Scerevisiae_Spombe/Scerevisiae_orf_coding_all_Scerevisiae_rna_coding_Spombe_cds.fasta
# snakemake --config genome=${genome} \
# gtf=${gtf} \
# transcriptome=${transcriptome} \
# --cores 1 -p --use-conda

# genome=Scerevisiae/saccharomyces_cerevisiae_R64-3-1_20210421_allchrom.fasta
# gtf=Scerevisiae/saccharomyces_cerevisiae_R64-3-1_20210421_nofasta_geneid.gff
# transcriptome=Scerevisiae/Scerevisiae_orf_coding_all_Scerevisiae_rna_coding.fasta
# snakemake --config genome=${genome} \
# gtf=${gtf} \
# transcriptome=${transcriptome} \
# --cores 1 -p --use-conda

genome=Scerevisiae_NL_RL/NL_RL_saccharomyces_cerevisiae_R64-3-1_20210421.fasta
gtf=Scerevisiae_NL_RL/NL_RL_saccharomyces_cerevisiae_R64-3-1_20210421_geneid.gff3
transcriptome=Scerevisiae_NL_RL/NL_RL_Scerevisiae_orf_coding_all_Scerevisiae_rna_coding.fasta
snakemake --config genome=${genome} \
gtf=${gtf} \
transcriptome=${transcriptome} \
--cores 1 -p --use-conda
