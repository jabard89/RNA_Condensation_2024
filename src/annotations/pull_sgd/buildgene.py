import sys, os, math, gzip, collections
import util, biofile, translate
from urllib.parse import unquote
import scipy as sp
import numpy as np

class CountRecord(object):
	chromosome = None
	count = None
	strand = None
	coordinate = None

def parseSGRLine(line_list, strand):
	"""Parse a single SGR entry.
	Entries look like:
	chrIV	1339815	4
	"""
	rec = CountRecord()
	chrom = line_list[0]
	rec.chromosome = chrom
	rec.strand = strand
	rec.coordinate = int(line_list[1])
	rec.count = int(line_list[2])
	return(rec)

def getChromosomeNumber(first_line):
	"""
	When reading chromosomes from a fasta file we want to nice the chromosome I or IV etc. 
	The default FASTA parcer gives us the wrong name. So this fuction solves that.
	The first line starts with something like ">ref|NC_001224| " then continues with information in brackets [etc]
	If we split by space, then [S. cerevisiae] is inappropriately borken up. So split by open bracket, then strip
	off the close bracket later. Now information in brackets is of form [identifier = value]. So lets make a dictionary
	of values keyed by identifiers. Chromosomes have [Chromosome=I] so we just return the dictionary value. But the 
	mitochondria doesn't have a Chromosome field. It has a locaiton field though. So in the absence of chromosome, 
	we return location.
	God help us if something doesn't have location.
	Input: first line of a fast entry
	output: string indicating location of fasta entry such as chromosome number or mitochondrion 
	"""

	res = first_line.split(' [')
	header_dict={}
	for info in res[1:]:
		header_dict[info.strip(']').split('=')[0]] = info.strip(']').split('=')[1]
	if 'chromosome' in header_dict:
		#record objects from the GFF has chromosomes listed in the form 'chrIV' etc
		return 'chr' + header_dict['chromosome']
	elif 'location' in header_dict:
		if header_dict['location'] == 'mitochondrion':
			#record objects from the GFF have mito genes in the forms chrmt or chrMito
			return 'chr' + 'Mito'

def nameMapperFactory(mapdict):
	"""This function returns a function which maps names onto other names
	with a dictionary. If it has no mapping, it returns
	the passed-in value.

	Example call: mapper = nameMapperFactory({"chrmt":"chrMito"})
	Then:
		mapper("chrmt") = "chrMito"
		mapper("novalue") = "novalue"
	"""
	def fn(name):
		res = name
		try:
			res = mapdict[name]
		except KeyError:
			pass
		return res
	return fn

def noop(x):
	"""This does exactly nothing."""
	return x


###
#
#  Main object model:
#  Genes are made up of Regions.
#  Regions are named chunks of sequence corresponding to 5'UTR, 3'UTR, CDS, and introns.
#  Chromosomes are sequence and count repositories. Regions fetch data from Chromosomes.
#  ChromosomeCollection represents the genome. It provides methods for loading multi-chromosome
#  data and looking up genes without knowing their chromosome.
#
###

class Gene(object):
	def __init__(self, gene_name, chromosome):
		self.name = gene_name
		self.type = None
		self.chromosome = chromosome
		self.regions = []
		self.strand = None
		self._commonNames = []

	def addAttributes(self, attribute_dict):
		self._attr_dict = attribute_dict

	@property
	def systematicName(self):
		return self.name

	@property
	def commonName(self):
		if len(self._commonNames)>0:
			return self._commonNames[0]
		else:
			return self.name

	def getAttributes(self):
		return self._attr_dict

	@property
	def attributes(self):
		return self._attr_dict
	

	def attr(self, attr_name):
		res = None
		try:
			res = self._attr_dict[attr_name]
		except KeyError:
			pass
		return res

	def addRegion(self, feature, strand, start, end):
		region = Region(self.chromosome)
		region.initialize(feature, strand, start, end)	
		self.regions.append(region)
		if self.strand is None:
			self.strand = strand
		self.sort()

	def addRegionFromGFF(self, gff_record):
		region = Region(self.chromosome)
		region.createFromGFF(gff_record)
		if self.strand is None:
			self.strand = region.strand
		self.regions.append(region)
		self.sort()

	def get5PrimeEnd(self):
		"""Return the 5' coordinate"""
		res = None
		if len(self.regions)>0:
			if self.strand == '+':
				res = self.regions[0].start
			elif self.strand == '-':
				res = self.regions[-1].end
		return res

	def get3PrimeEnd(self):
		"""Return the 3' coordinate"""
		res = None
		if len(self.regions)>0:
			if self.strand == '+':
				res = self.regions[-1].end
			elif self.strand == '-':
				res = self.regions[0].start
		return res

	def sort(self):
		"""Sort regions by coordinate"""
		self.regions = sorted(self.regions, key=lambda x: x.start)

	def addRegions(self, gff_records):
		for rec in gff_records:
			self.addRegionFromGFF(rec)

	def getRegions(self, feature_types=None):
		if isinstance(feature_types, str):
			feature_types = [feature_types]
		if feature_types is None:
			res = [reg for reg in self.regions]
		else:
			res = [reg for reg in self.regions if reg.feature in feature_types]
		return res

	def getRegion(self, feature_type):
		res = [reg for reg in self.regions if reg.feature == feature_type]
		if len(res) > 0:
			return res[0]
		else:
			return None

	def add5PrimeUTR(self, utr_length):
		start = self.get5PrimeEnd()
		utr5_start = None
		utr5_end = None
		if self.strand == "+":
			utr5_start = start-1 - (utr_length - 1) # inclusive
			utr5_end = start-1
		else:
			utr5_start = start+1
			utr5_end = start+1 + utr_length - 1 # inclusive
		assert utr5_end-utr5_start + 1 == utr_length
		self.addRegion("UTR5", self.strand, utr5_start, utr5_end)

	def add3PrimeUTR(self, utr_length):
		start = self.get3PrimeEnd()
		utr3_start = None
		utr3_end = None
		if self.strand == "-":
			utr3_start = start-1 - (utr_length - 1)
			utr3_end = start-1
		else:
			utr3_start = start+1
			utr3_end = start+1 + utr_length - 1
		assert utr3_end-utr3_start + 1 == utr_length
		self.addRegion("UTR3", self.strand, utr3_start, utr3_end)

	def get5PrimeUTRSequence(self):
		seq = ''
		seqreg_list = [reg for reg in self.regions if reg.feature == 'UTR5']
		if self.strand == "-":
			seqreg_list = seqreg_list[::-1]
		for seqreg in seqreg_list:
			seq += seqreg.nucleotides
		return seq

	def get3PrimeUTRSequence(self):
		seq = ''
		seqreg_list = [reg for reg in self.regions if reg.feature == 'UTR3']
		if self.strand == "-":
			seqreg_list = seqreg_list[::-1]
		for seqreg in seqreg_list:
			seq += seqreg.nucleotides
		return seq

	def getCodingSequence(self):
		seq = ''
		cds_list = [reg for reg in self.regions if reg.feature == 'CDS']
		if self.strand == "-":
			cds_list = cds_list[::-1]
		for cds in cds_list:
			seq += cds.nucleotides
		return seq

	def getNumberOfExons(self):
		return len(self.regions)

	def getSplicedTranscript(self):
		"""
		Retrieve the spliced transcript: 5'UTR, coding sequence, and 3'UTR concatenated, in that order.
		Introns are eliminated.
		"""
		seq = ''
		if self.isCodingSequence():
			mature_regions = ['UTR5','CDS','UTR3']
			cds_list = [reg for reg in self.regions if reg.feature in mature_regions]
			if self.strand == "-":
				cds_list = cds_list[::-1]
			for cds in cds_list:
				seq += cds.nucleotides
		elif self.isNoncodingSequence():
			mature_regions = ['noncoding_exon']
			seq_list = [reg for reg in self.regions if reg.feature in mature_regions]
			if self.strand == "-":
				seq_list = seq_list[::-1]
			for seqregion in seq_list:
				seq += seqregion.nucleotides
		return seq

	def isCodingSequence(self):
		"""
		Is this a coding sequence?
		"""
		region_names = [reg.feature for reg in self.regions]
		return 'CDS' in region_names

	def isNoncodingSequence(self):
		"""
		Is this a noncoding sequence?
		"""
		region_names = [reg.feature for reg in self.regions]
		return 'noncoding_exon' in region_names

	def getBoundAndUnboundSequences(self, region_id, seq_length, threshold, n_samples):
		#Maybe later
		pass

class Region(object):
	"""
	A region of a chromosome.
	This object does not store sequence or count information. It retrieves these on demand
	through its .chromosome property.
	"""
	def __init__(self, chromosome):
		self.chromosome = chromosome
		self.feature = ''
		self.strand = ""
		self.start = -1
		self.end = -1
		#self.counts = []

	def initialize(self, feature, strand, start, end):
		self.feature = feature
		self.strand = strand
		self.start = start
		self.end = end
		#self.counts = []

	def initializeFromGFF(self, gff_record):
		self.initialize(gff_record.feature, gff_record.strand, gff_record.start, gff_record.end)

	def __len__(self):
		"""Return the length of the region."""
		return self.end - self.start + 1

	@property
	def nucleotides(self):
		"""Return the nucleotide sequence."""
		# CDK: I can't reason why start must be -1, but empirically, it works
		# DAD: Python sequence indices are zero-based and exclusive at the upper end: [0,3)
		# GFF indices are 1-based and inclusive at the upper end.
		# So given a sequence ACGC, the first three nucleotides are:
		# Python: sequence[0:3] = ACG
		# GFF: genome coordinate 1-3 = ACG
		nucs = self.chromosome.sequence[self.start-1: self.end]
		if self.strand == '-':
			nucs = translate.reverseComplement(nucs)
		return nucs

	@property
	def counts(self):
		"""Return the counts."""
		strand_counts = self.chromosome.strand_counts[self.strand]
		if self.strand == '+':
			counts = strand_counts[self.start-1: self.end]
		elif self.strand == '-':
			# Reverse the counts
			counts = strand_counts[self.start-1: self.end]
			counts = counts[::-1]
		return counts

	def __str__(self):
		return self.nucleotides

class Chromosome(object):
	"""
	A coherent collection of genes.

	Genes belong to their Chromosome object. Creation of a Gene requires a Chromosome, because
	the Chromosome object stores sequence and other nucleotide-level information. Genes are
	references to information in Chromosome objects.
	"""
	def __init__(self, name, dna_sequence):
		self.name = name
		self.gene_dict = {}
		self.sequence = dna_sequence
		self.strand_counts = {"+":[0]*len(dna_sequence), "-":[0]*len(dna_sequence)} # A sequence of integers, default = 0, as long as the DNA sequence
		pass

	def addCount(self, strand, coordinate_1, count):
		"""Add count information for this coordinate on this strand. coordinate_1 is 1-based coordinate."""
		try:
			self.strand_counts[strand][coordinate_1-1] = count
		except IndexError as ie:
			#print(self.name, strand, coordinate_1, len(self.sequence))
			if len(self.sequence) > coordinate_1 + 1:
				raise(ie)

	def addGene(self, gene_name):
		"""Add a new gene."""
		gene = Gene(gene_name, self)
		self.gene_dict[gene_name] = gene
		return gene

	def addRegion(self, gene_name, region_type, strand, start_coord_1i, end_coord_1i, record=None):
		"""
		Add a new region to the named gene. If the gene does not have a corresponding
		object, create a new one.
		"""
		gene = None
		try:
			# Retrieve the corresponding gene
			gene = self.gene_dict[gene_name]
		except:
			# Make a new gene if this one hasn't been seen before
			gene = self.addGene(gene_name)
			gene.type = region_type
			if not record is None: # Additional information in the record
				# Parse attributes and add them
				gene.addAttributes(record.attributes)
			gene.strand = strand
		# We've got the gene. Add the region.
		gene.addRegion(region_type, strand, start_coord_1i, end_coord_1i)

	def addRegionFromGFF(self, gene_name, gff_record):
		self.addRegion(gene_name, gff_record.feature, gff_record.strand, gff_record.start, gff_record.end, record=gff_record)

	def getGene(self, gene_name):
		"""Retrieve gene by name."""
		return self.gene_dict[gene_name]
	
	def __getitem__(self, gene_name):
		"""Shortcut for getting genes by ID"""
		return self.gene_dict[gene_name]
	
	@property
	def genes(self):
		for gene_name in self.gene_dict:
			yield self.gene_dict[gene_name]

	def getCounts(self, strand, start, end):
		return self.strand_counts[strand][(start-1):end]

class ChromosomeCollection(object):
	"""
	A collection of chromosomes -- a genome, in short.
	
	"""
	def __init__(self):
		self.chromosomes = {}
		self.gene_to_chromosome_mapping = {}

	def addChromosome(self, chromosome):
		self.chromosomes[chromosome.name] = chromosome

	def __getitem__(self, chromosome_name):
		return self.chromosomes[chromosome_name]

	def getChromosome(self, chromosome_name):
		return self.chromosomes[chromosome_name]

	def getGene(self, gene_name):
		chrom = self.gene_to_chromosome_mapping[gene_name]
		return chrom.getGene(gene_name)

	def loadRegions(self, gff_infile, nameMapper=noop):
		"""
		Returns a dictionary keyed by orf whose entires are lists of record objects.
		"""
		features_of_interest = ['gene','CDS','intron', 'ncRNA_gene','tRNA_gene','transposable_element_gene','snRNA_gene','rRNA_gene','snoRNA_gene','telomerase_RNA_gene',
			'external_transcribed_spacer_region', 'internal_transcribed_spacer_region',
			'blocked_reading_frame',
			'noncoding_exon']
		#orf_record_dict = collections.defaultdict(list)
		#orf_list = []
		gfr = biofile.GFFReader(gff_infile) #, unescape=True)
		for rec in gfr.entries:
			#id = rec['ID']
			chromosome_name = nameMapper(rec.seqname)
			chromosome = self.chromosomes[chromosome_name]
			if rec.feature in features_of_interest:
				feature_id = rec['Name'] # E.g., Name=YAL067C_CDS
				gene_name = feature_id.split('_')[0]
				#orf_list.append(orf)
				#orf_record_dict[orf].append(rec)
				chromosome.addRegionFromGFF(gene_name, rec)
				self.gene_to_chromosome_mapping[gene_name] = chromosome

		# Potentially multiple entries per gene, if there are introns
		#return (orf_record_dict, orf_list)

	def loadStrandCounts(self, strand, count_infile, nameMapper=noop):
		"""
		Load the SGR counts.
		These are structured as chromosome, coordinate, count tuples,
		for a specific strand.
		"""
		count_dlr = util.DelimitedLineReader(count_infile, header=False)
		cur_chromosome_name = None
		cur_chromosome = None
		for line in count_dlr.entries:
			# Parse the line
			chrom_name = line[0]
			coordinate = line[1]
			count = line[2]
			if chrom_name != cur_chromosome_name:
				# Cache chromosome so we don't need to look it up every time
				cur_chromosome_name = chrom_name
				cur_chromosome = self.chromosomes[nameMapper(chrom_name)]
				#print("# New chromosome", cur_chromosome_name)		
			cur_chromosome.addCount(strand, coordinate, count)

	def loadUTRLengths(self, utr_length_infile):
		"""
		This assumes the format of Pelechano et al. 2013:

		# S. cerevisiae 5' and 3' UTR estimates.
		# 
		# Extensive transcriptional heterogeneity revealed by isoform profiling
		# Vicent Pelechano Wu Wei, & Lars M. Steinmetz, Nature 2013
		# Supplementary Data 3, http://www.nature.com/nature/journal/v497/n7447/full/nature12121.html
		# header annotated by Edward Wallace
		# 
		# Column key
		#   chr: chromosome
		#   strand: DAN strand (determined from annotation)
		#   gene: systematic ORF name
		#   type: Verified - this file only includes verified ORFs
		#   mTIF_number: number of detected reads with median 5' UTR and median 3' UTR length
		#   median5: median 5' UTR length
		#   median3: median 3' UTR length
		#   sd5: standard deviation of 5' UTR length (interquartile range would be better)
		#   sd3: standard deviation of 3' UTR length
		#
		chr	strand	gene	type	mTIF_number	median5	median3	sd5	sd3
		1	-	YAL067C	Verified	1	33	222	NA	NA
		1	-	YAL063C-A	Uncharacterized	3	932	32	234.540685880581	180.133283987163
		1	-	YAL054C	Verified	1	51	163	NA	NA
		1	-	YAL049C	Verified	74	66	119	174.54993376684	28.0249003365666

		"""
		feat_dlr = util.DelimitedLineReader(utr_length_infile, header=True)
		for linedict in feat_dlr.dictentries:
			orf = linedict['gene']
			try:
				gene = self.getGene(orf)
				fiveprime_utr_length = int(np.ceil(linedict['median5']))
				threeprime_utr_length = int(np.ceil(linedict['median3']))
				gene.add5PrimeUTR(fiveprime_utr_length)
				gene.add3PrimeUTR(threeprime_utr_length)
			except KeyError as bad_orf:
				# YFL057C, for example, is a pseudogene in R64 2-1 release and will
				# not appear as a gene, yet has reported UTR lengths
				pass
	
	@property
	def genes(self):
		"""Iterates over genes"""
		for chr_name in self.chromosomes:
			chrom = self.chromosomes[chr_name]
			for gene in chrom.genes:
				yield gene



