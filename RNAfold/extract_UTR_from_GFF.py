# -*- coding: utf-8 -*-
"""
Created on 9/27/2021
@author: Jared Bard
Given an annotated GFF file (annotated by add_UTR_GFF.py), pulls the 5' UTR, 3' UTR,
and the 5' UTR + x residues into a tab formatted file

Remember python is 0-based and half-open intervals. Luckily both SeqIO and HTseq are too, so they match well
"""

import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
import HTSeq
import itertools
import re
# requires bedops package

def gff_todict(gff_filename):
    # reads a gff file and organizes all the features by gene
    # traces parents from feature -> mRNA -> gene
    gene_dict = {}
    transcript_dict = {}  # going to store transcripts to make it easy to find their parent genes
    gff_file = HTSeq.GFF_Reader(gff_filename,end_included=True)
    print_flag = False
    for feature in gff_file:
        if feature.type == "gene":
            try:
                if 'gene_id' in feature.attr:  # this is right for cerevisiae
                    gene_id = feature.attr['gene_id']
                    if not print_flag:
                        print("Names from the 'gene_id' attribute")  # print this to the output once
                        print_flag = True
                else:
                    gene_id = feature.attr['Name']  # this seems better for all other provided gffs in yeasTSS
                    if not print_flag:
                        print("Name from the 'Name' attribute")
                        print_flag = True
            except:
                print("No name found for {}".format(feature.name))

            # now add the gene to the gene dictionary
            gene_dict[feature.name] = { 'gene_id':gene_id,'feature':feature,
                                        'transcripts':[],'five_prime_UTR':[],'three_prime_UTR':[],
                                        'exons':[],'CDSs':[] }

    # scan again for mRNAs associated with each gene
    for feature in gff_file:
        if feature.type == "mRNA":
            if not feature.name in transcript_dict:
                transcript_dict[feature.name] = feature.attr['Parent']
            try:
                gene_dict[feature.attr['Parent']]['transcripts'].append(feature)
            except:
                print("Cannot find parent {} of transcript {}".format(feature.attr['Parent'],feature.name))

    # scan for features
    def find_parent(feature):
        parent_gene = None
        try:
            parent_gene = transcript_dict[feature.attr['Parent']]
        except:
            print("Cannot find parent of {}".format(feature.name))
        return(parent_gene)

    for feature in gff_file:
        if feature.type not in ['exon','CDS','five_prime_UTR','three_prime_UTR']:
            continue
        if not find_parent(feature):
            continue
        parent_gene = find_parent(feature)
        if feature.type == 'exon':
            gene_dict[parent_gene]['exons'].append(feature)
        elif feature.type == 'CDS':
            gene_dict[parent_gene]['CDSs'].append(feature)
        elif feature.type == 'five_prime_UTR':
            gene_dict[parent_gene]['five_prime_UTR'].append(feature)
        elif feature.type == 'three_prime_UTR':
            gene_dict[parent_gene]['three_prime_UTR'].append(feature)

    # output a dictionary but with gene_id as keys
    out_dict = {}
    for gene in gene_dict:
        out_dict[gene_dict[gene]['gene_id']] = {
            'feature':gene_dict[gene]['feature'],
            'transcripts':gene_dict[gene]['transcripts'],
            'five_prime_UTR':gene_dict[gene]['five_prime_UTR'],
            'three_prime_UTR':gene_dict[gene]['three_prime_UTR'],
            'exons':gene_dict[gene]['exons'],
            'CDSs':gene_dict[gene]['CDSs']
        }
    return(out_dict)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Generate options")
    # Required arguments
    parser.add_argument(dest="gff_file",default=None,type=str,help="input GFF sequence")
    parser.add_argument(dest="fasta",default=None,type=str,help="Genome fasta file")
    parser.add_argument(dest="output_file",default=None,type=str,help="output file")
    parser.add_argument("--add_nts",dest="add_nts",default=0,type=int,help="how many nts of orf to add to UTR5")
    args = parser.parse_args()
    # pyscript_dir="/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools"
    # args = parser.parse_args(['--add_nts=20',
    #                     pyscript_dir+'/UTRs/Saccharomyces_cerevisiae.R64-1-1.109_yeasTSS_Pelechano2013.gff3',
    #                     pyscript_dir+'/UTRs/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa',
    #                     pyscript_dir+'/test.csv'])
    
    # Write out parameters
    with open(args.output_file,'w') as f:
        f.write("# Run started {}\n".format(datetime.now()))
        f.write("# Command: {}\n".format(' '.join(sys.argv)))
        f.write("# Parameters:\n")
        optdict = vars(args)
        for (k,v) in optdict.items():
            f.write("#\t{k}: {v}\n".format(k=k,v=v))

    gff_dict = gff_todict(args.gff_file)
    chrom_dict = SeqIO.to_dict(SeqIO.parse(args.fasta,"fasta"))


    # create a dictionary referenced by gene_ID for every gene
    def extract_seq_fromiv(iv,chrom_dict):
        chrom = iv.chrom
        strand = iv.strand
        if strand == "+":
            seq = str(chrom_dict[chrom][iv.start:iv.end].seq)
        elif strand == "-":
            seq = str(chrom_dict[chrom][iv.start:iv.end].reverse_complement().seq)
        return (seq)


    out_dict = { }
    for gene in gff_dict:
        out_dict[gene] = { 'ORF':gene,'UTR5':'','UTR5_plus' + str(args.add_nts):'','UTR3':'' }
        if gff_dict[gene]['five_prime_UTR']:
            out_dict[gene]['UTR5'] = extract_seq_fromiv(gff_dict[gene]['five_prime_UTR'][0].iv,chrom_dict)
            ivplus = gff_dict[gene]['five_prime_UTR'][0].iv.copy()
            ivplus.length = ivplus.length + args.add_nts
            out_dict[gene]['UTR5_plus' + str(args.add_nts)] = extract_seq_fromiv(ivplus,chrom_dict)
        if gff_dict[gene]['three_prime_UTR']:
            out_dict[gene]['UTR3'] = extract_seq_fromiv(gff_dict[gene]['three_prime_UTR'][0].iv,chrom_dict)

    # output
    out_frame = pd.DataFrame.from_dict(out_dict,orient='index')
    out_frame.to_csv(args.output_file,mode='a',sep=',',header=True,index=False)