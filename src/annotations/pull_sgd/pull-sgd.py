#! python

import sys, os, math, random, argparse, collections, gzip
from urllib.parse import unquote
import wget, ssl
import util, biofile, translate, buildgene

"""
The purpose of this script is to extract gene and protein sequences in FASTA format and features in tab-delimited format from the Saccharomyces Genome Database
weekly build.

The Saccharomyces Genome Database GFF file, updated weekly, is at https://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff

This GFF file contains both a proper GFF and the genome sequence to which it refers.

Example usage:

# This writes out features to the console and puts nucleotide and protein files into scer-nt.fasta and scer-aa.fasta
python3 pull-sgd.py 

# This writes the features to a file
python3 pull-sgd.py --out scer-features.txt

# This specifies all the file locations
python3 pull-sgd.py --out ../genome/scer-features.txt --nucleotide-out ../genome/scer-nt.fasta --protein-out ../genome/scer-aa.fasta 

# This allows you to use a local copy of the GFF feature file (can be combined with other options)
python3 pull-sgd.py --gff saccharomyces_cerevisiae.gff

# Use local feature file, output proteins, CDS, transcripts, and features to separate files
python3 pull-sgd.py --gff saccharomyces_cerevisiae.gff.gz --protein-out scer-aa.fa --nucleotide-out scer-cds.fa --transcript-out scer-rna.fa --out scer-features.txt

# Include UTR sequences and lengths
python3 pull-sgd.py --gff saccharomyces_cerevisiae.gff.gz --utr-file ../RNA/Pelechano2013_SuppData3_medianTranscripts.txt --protein-out scer-aa.fa --nucleotide-out scer-cds.fa --transcript-out scer-rna.fa --out scer-features.txt
python3 pull-sgd.py --utr-file ../RNA/Pelechano2013_SuppData3_medianTranscripts.txt --protein-out scer-aa.fa --nucleotide-out scer-cds.fa --transcript-out scer-rna.fa --out scer-features.txt
"""

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Generate genes and features from the Saccharomyces Genome Database weekly build")
	# Optional arguments
	parser.add_argument("--gff", dest="gff", default="http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz", help="GFF file")
	parser.add_argument("--mapping", dest="chr_mapping", default=None, help="chromosome mapping file")
	parser.add_argument("--attr-parser", dest="attribute_parser", default=None, help="parser for attributes ('Default','JGI')")
	parser.add_argument("--id-field", dest="id_field", default='ID', help="name of GFF attribute storing identifiers")
	parser.add_argument("--gene-name-field", dest="gene_name_field", default='gene', help="name of GFF attribute storing gene name")
	#parser.add_argument("--id-field", dest="id_field", default='ID', help="name of GFF attribute storing gene identifiers")
	parser.add_argument("--gene-prefix", dest="gene_prefix", default='', help="prefix to put before gene name")
	# Add support for UTRs. Lump all together for starters: single argument
	# --utrs 
	parser.add_argument("--utr-file", dest="utr_filename", type=str, default=None, help="tab-delimited UTR length file")
	parser.add_argument("--split-utrs", dest="split_utrs", default=False, action="store_true", help="split 5' and 3' UTRs separately in CDS FASTA output?")
	parser.add_argument("--out", dest="out_fname", default=None, help="output filename")
	parser.add_argument("--transcript-out", dest="rna_fasta_out_fname", default="scer-rna.fasta", help="output FASTA filename for mRNAs")
	parser.add_argument("--nucleotide-out", dest="nt_fasta_out_fname", default="scer-nt.fasta", help="output FASTA filename for nucleotides")
	parser.add_argument("--protein-out", dest="aa_fasta_out_fname", default="scer-aa.fasta", help="output FASTA filename for proteins")
	options = parser.parse_args()

	info_outs = util.OutStreams(sys.stdout)
	data_outs = util.OutStreams()

	# Start up output
	if not options.out_fname is None:
		outf = open(options.out_fname,'w')
		data_outs.addStream(outf)
	else:
		# By default, write to stdout
		data_outs.addStream(sys.stdout)

	attr_parser = biofile.parseAttributesDefault
	if options.attribute_parser == 'JGI':
		attr_parser = biofile.parseAttributesJGI1

	# This SSL context is necessary on Mac because certificates are not available by default
	# to the wget URL loader
	# cf. https://stackoverflow.com/questions/35569042/ssl-certificate-verify-failed-with-python3
	ssl._create_default_https_context = ssl._create_unverified_context
	# Download the GFF file
	if options.gff.startswith("http"): # assume this is a URL and pull from it
		gff_fname = wget.download(options.gff)
	else:
		gff_fname = options.gff
	# Split the GFF file
	gff_filename = os.path.expanduser("./gff.tmp")
	fasta_filename = os.path.expanduser("./fasta.tmp")
	gff_only_outf = open(gff_filename, 'w',encoding='utf8')
	fasta_only_outf = open(fasta_filename, 'w',encoding='utf8')
	cur_outf = gff_only_outf

	if gff_fname.endswith(".gz"):
		inf = gzip.open(gff_fname,'rt',encoding='utf8')
	else:
		inf = open(gff_fname,'r',encoding='utf8')
	for line in inf:
		if line.startswith('##FASTA'):
			cur_outf = fasta_only_outf
		else:
			cur_outf.write(line)
	inf.close()
	gff_only_outf.close()
	fasta_only_outf.close()

	# "lambda x: x" here reflects that the chromosome name is the first (and only) identifier for each FASTA entry for SGD
	chr_dict = biofile.readFASTADict(fasta_filename, lambda x: x.split()[0])
	# Load sequence
	chromcol = buildgene.ChromosomeCollection()
	for chr_name in chr_dict:
		chromosome = buildgene.Chromosome(chr_name, chr_dict[chr_name])
		chromcol.addChromosome(chromosome)
		#print(chr_name)
	# Load genes
	with open(gff_filename,'r') as inf:
		chromcol.loadRegions(inf) #, buildgene.nameMapperFactory({"chrmt":"chrMito"}))
	# Load UTRs
	have_utrs = False
	if not options.utr_filename is None:
		with open(os.path.expanduser(options.utr_filename), 'r') as inf:
			chromcol.loadUTRLengths(inf)
			have_utrs = True

	# Write out parameters
	data_outs.write("# Run started {}\n".format(util.timestamp()))
	data_outs.write("# Command: {}\n".format(' '.join(sys.argv)))
	data_outs.write("# Parameters:\n")
	optdict = vars(options)
	for (k,v) in optdict.items():
		data_outs.write("#\t{k}: {v}\n".format(k=k, v=v))

	gene_names = sorted([gene.name for gene in chromcol.genes])

	prot_seqs = {}
	cds_seqs = {}
	rna_seqs = {}
	header_dict = {}
	for gene in chromcol.genes:
		gene_name = gene.name
		common_name = gene.attr("gene")
		if common_name is None:
			common_name = gene_name
		attr_string = ';'.join(["{}={}".format(k, v) for (k,v) in gene.attributes.items()])
		#print(gene_name, common_name, attr_string)
		header = "{:s} {:s} {:s}".format(gene.name, common_name, attr_string)
		header_dict[gene_name] = header
		rna = gene.getSplicedTranscript()
		rna_seqs[gene_name] = rna
		if gene.isCodingSequence():
			seq = gene.getCodingSequence()
			cds_seqs[gene_name] = seq
			prot = translate.translateRaw(seq)
			if not prot is None:
				# Not an empty sequence
				if prot.count('*')>1:
					# Multiple stop codons encountered. Maybe wrong genetic code.
					# Try mitochondrial code
					prot_mito = translate.translateRaw(seq, translate._scer_mito_code)
					if prot_mito.count('*') <= 1:
						# Guess that this is the correct translation
						prot = prot_mito
				if prot.endswith("*"):
					prot = prot[:-1]
				prot_seqs[gene_name] = prot

	if not options.nt_fasta_out_fname is None:
		nt_outf = open(options.nt_fasta_out_fname,'w')
		# Write the CDS
		keys = sorted(cds_seqs.keys())
		biofile.writeFASTA([cds_seqs[k] for k in keys], nt_outf, headers=[header_dict[k] for k in keys])

	if not options.rna_fasta_out_fname is None:
		rna_outf = open(options.rna_fasta_out_fname,'w')	
		# Write the transcripts
		keys = sorted(rna_seqs.keys())
		biofile.writeFASTA([rna_seqs[k] for k in keys], rna_outf, headers=[header_dict[k] for k in keys])

	if not options.aa_fasta_out_fname is None:
		aa_outf = open(options.aa_fasta_out_fname,'w')
		# Write the proteins
		keys = sorted(prot_seqs.keys())
		biofile.writeFASTA([prot_seqs[k] for k in keys], aa_outf, headers=[header_dict[k] for k in keys])

	# Write output
	dout = util.DelimitedOutput()
	dout.addHeader('ORF','Systematic gene name','s')
	dout.addHeader('gene','Alternative or common gene name','s')
	dout.addHeader('classification','classification of gene','s')
	dout.addHeader('type','type of gene','s')
	dout.addHeader('length.cds','Length of coding sequence in nucleotides','d')
	dout.addHeader('length.transcript','Length of spliced RNA in nucleotides','d')
	dout.addHeader('length.protein','Length of protein in amino acids','d')
	if have_utrs:
		dout.addHeader('utr5.length', 'Length of 5\'UTR in nucleotides','d')
		dout.addHeader('utr3.length', 'Length of 3\'UTR in nucleotides','d')
	dout.addHeader('num.exons','Number of exons','d')
	dout.addHeader('desc','Description','s')
	# Write the header descriptions
	dout.describeHeader(data_outs)
	# Write the header fiels
	dout.writeHeader(data_outs)
	#format = dout.getFormat(named=True)
	n_written = 0
	for gene_name in gene_names:
		gene = chromcol.getGene(gene_name)
		is_transcribed = gene.isCodingSequence() or gene.isNoncodingSequence()
		if is_transcribed:
			result = dout.createResult(default=None)
			result['ORF'] = gene.name
			result['type'] = gene.type
			result['gene'] = gene.attr("gene") #gene_dict[gene].name
			result['classification'] = gene.attr("orf_classification") #gene_dict[gene].classification
			result['num.exons'] = gene.getNumberOfExons() #cds_dict[gene].numExons()
			cds_length = None
			prot_length = None
			transcript_length = len(rna_seqs[gene_name])
			utr5_length = None
			utr3_length = None
			if gene.isCodingSequence():
				try:
					cds_length = len(cds_seqs[gene_name])
					prot_length = len(prot_seqs[gene_name])
					utr5_length = len(gene.get5PrimeUTRSequence())
					utr3_length = len(gene.get3PrimeUTRSequence())
				except KeyError:
					pass
			result['length.protein'] = prot_length
			result['length.cds'] = cds_length
			result['length.transcript'] = transcript_length
			result['utr5.length'] = utr5_length
			result['utr3.length'] = utr3_length
			desc = None
			try:
				desc = unquote(gene.attributes["Note"])
			except KeyError:
				pass
			result['desc'] = desc
			# Parse the values, convert Nones to NA, etc.
			line = dout.formatLine(result)
			# A more manual approach:
			# line = format.format(column1="the answer is", column2=42)
			data_outs.write(line)
			n_written += 1
	# Write out stopping time
	data_outs.write("# Run finished {}\n".format(util.timestamp()))

	# Shut down output
	if not options.out_fname is None:
		info_outs.write("# Wrote {} lines to {}\n".format(n_written, options.out_fname))
		outf.close()
