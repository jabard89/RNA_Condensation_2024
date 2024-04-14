"""
Created on 02/21/2023
@author: Jared Bard
For use with STAR --quantMode, add a gene_id attribute to every feature in a group
Uses gffutils
"""

import sys, os, math, random, argparse, csv
from datetime import datetime
import pandas as pd
import numpy as np
from Bio import SeqIO
import gffutils

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate options")
	# Required arguments
	parser.add_argument(dest="gff_file",default=None,type=str,help="input GFF3 file")
	parser.add_argument(dest="output_file",default=None,type=str,help="file to export gff3 to")
	args = parser.parse_args()
	# args = parser.parse_args(["/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/src/"
	# 							"annotations/Schizosaccharomyces_pombe_all_chromosomes.gff3",                          
	#                           	"/home/jabard89/Dropbox/code_JB/repos/rnaseq_tools/src/"
	# 							"annotations/Spombe_test.gff3"])
	
	for file in [args.gff_file]:
		if not os.path.isfile(file):
			raise IOError("# Error: file {} does not exist".format(file))	
	
	for file in [args.output_file]:
		if os.path.isfile(file):
			raise IOError("# Error: file {} already exists, please remove it".format(file))
	
	# Write out parameters
	with open(args.output_file,'w') as f:
		f.write("# Run started {}\n".format(datetime.now()))
		f.write("# Command: {}\n".format(' '.join(sys.argv)))
		f.write("# Parameters:\n")
		optdict = vars(args)
		for (k,v) in optdict.items():
			f.write("#\t{k}: {v}\n".format(k=k,v=v))
		
	db = gffutils.create_db(args.gff_file, ":memory:")

	with open(args.output_file, 'a') as fout:
		for d in db.directives:
			fout.write('## {0}\n'.format(d))

		for feature in db.all_features():
			if feature.featuretype == 'gene':
				feature['gene_id'] = feature.id
			else:
				gene_parent = list(db.parents(feature,featuretype='gene'))
				if gene_parent:
					if len(gene_parent) > 1:
						print("feature {} has multiple gene parents".format(feature.id))
					else:
						feature['gene_id'] = gene_parent[0].id
			fout.write(str(feature) + '\n')
		

	