import sys, os, math, string
# DAD: right now, translate is only for CodingSequence's call to reverseTranslate.
import translate, util
from Bio import SeqIO

class BioFileError(BaseException):
	"""Class for packaging errors reading biological data from files"""
	def __init__(self, s):
		return

class CodingExonRecord:
	def __init__( self ):
		self.gene_ID = ''
		self.transcript_ID = ''
		self.peptide_ID = ''
		self.exon_ID = ''
		self.strand = ''
		self.exon_start = 0
		self.exon_end = 0
		self.coding_start = 0
		self.coding_end = 0
		self.chromosome = ''
		self.descriptors = None

	def getKey( self ):
		return "{0}-{1}".format(self.peptide_ID, self.coding_start)

	def addDescriptors(self, desc_dict):
		if self.descriptors is None:
			self.descriptors = desc_dict
		else:
			self.descriptors = dict(self.descriptors.items() + desc_dict.items())

	def getDescriptor(self, key):
		return self.descriptors.get(key,None)

	def pullCodingSequence(self, raw_seq, n_bases_upstream=0, n_bases_downstream=0):
		# DAD: be careful with 0-based or 1-based indexing!
		# coding_start/end are 1-based, and inclusive
		start = max(0, self.coding_start-1-n_bases_upstream)
		end = min(len(raw_seq), self.coding_end+n_bases_downstream)
		return raw_seq[start:end]

	def extend(self, n_bases_upstream, n_bases_downstream):
		if self.strand == '-':
			self.coding_start = max(0,self.coding_start-n_bases_downstream)
			self.coding_end = self.coding_end+n_bases_upstream
		else:
			self.coding_start = max(0,self.coding_start-n_bases_upstream)
			self.coding_end = self.coding_end+n_bases_downstream

	def contains(self, pos):
		return self.coding_start <= pos and self.coding_end >= pos

	def getCodingLength(self):
		return self.coding_end - self.coding_start + 1

	def printAttributes( self ):
		print(self)

	def printRecord( self ):
		print(self)

	def __repr__(self):
		return "{} {} {} {} {} {} {} {} {} {}".format( self.chromosome, self.gene_ID, self.transcript_ID, self.peptide_ID,
												   self.exon_ID, self.strand, self.exon_start, self.exon_end, self.coding_start, self.coding_end )

class CodingSequence:
	def __init__(self):
		self.exons = []
		#self.touched_sequence = True
		self.sequence = None
		self.strand = None
		self.chromosome = None
		self.id = None
		self.gene = None

	def add(self, rec):
		self.exons.append(rec)
		#self.touched_sequence = True
		# Assume that CDS come from the same chromosome and strand always.
		self.strand = rec.strand
		self.chromosome = rec.chromosome
		gene = rec.getDescriptor('gene')
		if not gene is None:
			self.gene = gene

	def setID(self, id):
		self.id = id
		if self.gene is None:
			self.gene = id

	def getID(self):
		return self.id

	def getSequence(self, raw_seq, n_bases_upstream=0, n_bases_downstream=0):
		# Sort exon records
		self.exons.sort(key=lambda x: x.coding_start)
		# Figure out which direction "upstream" is
		if self.strand == '-':
			# Swap -- downstream is really upstream on the negative strand
			(n_bases_upstream, n_bases_downstream) = (n_bases_downstream, n_bases_upstream)
		seq = ''
		if len(self.exons) == 1:
			seq += self.exons[0].pullCodingSequence(raw_seq, n_bases_upstream, n_bases_downstream)
		else:
			for (i,x) in enumerate(self.exons):
				if i==0:
					seq += x.pullCodingSequence(raw_seq, n_bases_upstream, 0)
				elif i==(len(self.exons)-1):
					seq += x.pullCodingSequence(raw_seq, 0, n_bases_downstream)
				else:
					seq += x.pullCodingSequence(raw_seq, 0, 0)
		if self.strand == '-':
			seq = translate.reverseComplement(seq)
		return seq

	#def truncate(self, length):
	#	trunc_cds =
	#	assert False, "truncate not implemented"
	def extend(self, n_upstream, n_downstream):
		self.exons.sort(key=lambda x: x.coding_start)
		self.exons[0].extend(n_upstream,0)
		self.exons[-1].extend(0,n_downstream)

	def getStrand(self):
		return self.strand

	def getChromosome(self):
		return self.chromosome

	def getCodingStart(self):
		res = None
		self.exons.sort(key=lambda x: x.coding_start)
		if self.strand == '-':
			res = self.exons[-1].coding_end
		else:
			res = self.exons[0].coding_start
		return res

	def getCodingEnd(self):
		res = None
		self.exons.sort(key=lambda x: x.coding_start)
		if self.strand == '-':
			res = self.exons[0].coding_start
		else:
			res = self.exons[-1].coding_end
		return res

	def mapPosition(self, pos):
		# Get the index of pos in nucleotides from the coding start (from the A of the ATG), 1-based.
		# First, find the containing exon
		self.exons.sort(key=lambda x: x.coding_start)
		ind = None
		res = None
		coding_length = 0
		if self.strand == '-':
			# Count from the back
			for i in range(len(self.exons)-1,-1,-1):
				e = self.exons[i]
				if e.contains(pos):
					ind = i
				if ind is None:
					coding_length += e.getCodingLength()
			if not ind is None:
				# Get the total coding length of the trailing exons, plus
				# the displacement into the exon containing pos
				res = coding_length + (self.exons[ind].coding_end - pos) + 1
		else:
			# Count from the front
			for i in range(len(self.exons)):
				e = self.exons[i]
				if e.contains(pos):
					ind = i
				if ind is None:
					coding_length += e.getCodingLength()
			if not ind is None:
				# Get the total coding length of the previous exons, plus
				# the displacement into the exon containing pos
				res = coding_length + pos - self.exons[ind].coding_start + 1
		return res

	def getFASTAHeader(self):
		res = '{0} {1} {2} from '.format(self.id, self.gene, self.chromosome)
		ind_flds = ','.join(['{0}-{1}'.format(r.coding_start, r.coding_end) for r in self.exons])
		res += ind_flds + ' ({0}) '.format(self.strand)
		exclude_keys = ['note']
		if len(self.exons)>0:
			d = self.exons[0].descriptors
			# DAD: Hack!  Use URL encoder from somewhere for this...
			res += '"{0}" '.format(d['note'].replace("%20",' ').replace("%3B",';').replace("%2C",',').replace("%2F",'/').replace("%3A",':').replace("%2B",'+'))
			ks = d.keys()
			ks.sort()
			res += ';'.join(["{0}={1}".format(k,d[k]) for k in ks if not k in exclude_keys])
		return res

	def numExons(self):
		return len(self.exons)

'''
http://song.sourceforge.net/gff3.shtml

Undefined fields are replaced with the "." character, as described in
the original GFF spec.

Column 1: "seqid"

The ID of the landmark used to establish the coordinate system for the
current feature. IDs may contain any characters, but must escape any
characters not in the set [a-zA-Z0-9.:^*$@!+_?-|].  In particular, IDs
may not contain unescaped whitespace and must not begin with an
unescaped ">".

Column 2: "source"

The source is a free text qualifier intended to describe the algorithm
or operating procedure that generated this feature.  Typically this is
the name of a piece of software, such as "Genescan" or a database
name, such as "Genbank."  In effect, the source is used to extend the
feature ontology by adding a qualifier to the type creating a new
composite type that is a subclass of the type in the type column.

Column 3: "type"

The type of the feature (previously called the "method").  This is
constrained to be either: (a) a term from the "lite" sequence
ontology, SOFA; or (b) a SOFA accession number.  The latter
alternative is distinguished using the syntax SO:000000.

Columns 4 & 5: "start" and "end"

The start and end of the feature, in 1-based integer coordinates,
relative to the landmark given in column 1.  Start is always less than
or equal to end.

For zero-length features, such as insertion sites, start equals end
and the implied site is to the right of the indicated base in the
direction of the landmark.

Column 6: "score"

The score of the feature, a floating point number.  As in earlier
versions of the format, the semantics of the score are ill-defined.
It is strongly recommended that E-values be used for sequence
similarity features, and that P-values be used for ab initio gene
prediction features.

Column 7: "strand"

The strand of the feature.  + for positive strand (relative to the
landmark), - for minus strand, and . for features that are not
stranded.  In addition, ? can be used for features whose strandedness
is relevant, but unknown.

Column 8: "phase"

For features of type "CDS", the phase indicates where the feature
begins with reference to the reading frame.  The phase is one of the
integers 0, 1, or 2, indicating the number of bases that should be
removed from the beginning of this feature to reach the first base of
the next codon. In other words, a phase of "0" indicates that the next
codon begins at the first base of the region described by the current
line, a phase of "1" indicates that the next codon begins at the
second base of this region, and a phase of "2" indicates that the
codon begins at the third base of this region. This is NOT to be
confused with the frame, which is simply start modulo 3.

For forward strand features, phase is counted from the start
field. For reverse strand features, phase is counted from the end
field.

The phase is REQUIRED for all CDS features.

Column 9: "attributes"

A list of feature attributes in the format tag=value.  Multiple
tag=value pairs are separated by semicolons.  URL escaping rules are
used for tags or values containing the following characters: ",=;".
Spaces are allowed in this field, but tabs must be replaced with the
%09 URL escape.

These tags have predefined meanings:

    ID	   Indicates the name of the feature.  IDs must be unique
	   within the scope of the GFF file.

    Name   Display name for the feature.  This is the name to be
           displayed to the user.  Unlike IDs, there is no requirement
	   that the Name be unique within the file.

    Alias  A secondary name for the feature.  It is suggested that
	   this tag be used whenever a secondary identifier for the
	   feature is needed, such as locus names and
	   accession numbers.  Unlike ID, there is no requirement
	   that Alias be unique within the file.

    Parent Indicates the parent of the feature.  A parent ID can be
	   used to group exons into transcripts, transcripts into
	   genes, an so forth.  A feature may have multiple parents.
	   Parent can *only* be used to indicate a partof
	   relationship.

    Target Indicates the target of a nucleotide-to-nucleotide or
	   protein-to-nucleotide alignment.  The format of the
	   value is "target_id start end [strand]", where strand
	   is optional and may be "+" or "-".  If the target_id
	   contains spaces, they must be escaped as hex escape %20.

    Gap   The alignment of the feature to the target if the two are
          not collinear (e.g. contain gaps).  The alignment format is
	  taken from the CIGAR format described in the
	  Exonerate documentation.
	  (http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate
          ?cvsroot=Ensembl).  See "THE GAP ATTRIBUTE" for a description
	  of this format.

    Derives_from
          Used to disambiguate the relationship between one
          feature and another when the relationship is a temporal
          one rather than a purely structural "part of" one.  This
          is needed for polycistronic genes.  See "PATHOLOGICAL CASES"
	  for further discussion.

    Note   A free text note.

    Dbxref A database cross reference.  See the section
	   "Ontology Associations and Db Cross References" for
	   details on the format.

    Ontology_term  A cross reference to an ontology term.  See
           the section "Ontology Associations and Db Cross References"
	   for details.

Multiple attributes of the same type are indicated by separating the
values with the comma "," character, as in:

       Parent=AF2312,AB2812,abc-3

Note that attribute names are case sensitive.  "Parent" is not the
same as "parent".

All attributes that begin with an uppercase letter are reserved for
later use.  Attributes that begin with a lowercase letter can be used
freely by applications.

From http://genome.ucsc.edu/FAQ/FAQformat.html#format3
GFF (General Feature Format) lines are based on the GFF standard file format. GFF lines have nine required fields that must be tab-separated. If the fields are separated by spaces instead of tabs, the track will not display correctly. For more information on GFF format, refer to http://www.sanger.ac.uk/resources/software/gff.

Here is a brief description of the GFF fields:

seqname - The name of the sequence. Must be a chromosome or scaffold.
source - The program that generated this feature.
feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
start - The starting position of the feature in the sequence. The first base is numbered 1.
end - The ending position of the feature (inclusive).
score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
group - All lines with the same group are linked together into a single item.
Example:
Here's an example of a GFF-based track. Click here for a copy of this example that can be pasted into the browser without editing. NOTE: Paste operations on some operating systems will replace tabs with spaces, which will result in an error when the GFF track is uploaded. You can circumvent this problem by pasting the URL of the above example (http://genome.ucsc.edu/goldenPath/help/regulatory.txt) instead of the text itselfinto the custom annotation track text box.

track name=regulatory description="TeleGene(tm) Regulatory Regions"
chr22  TeleGene enhancer  1000000  1001000  500 +  .  touch1
chr22  TeleGene promoter  1010000  1010100  900 +  .  touch1
chr22  TeleGene promoter  1020000  1020000  800 -  .  touch2
'''
'''
import HTMLParser
html_parser = HTMLParser.HTMLParser()
unescaped = html_parser.unescape(my_string)
'''

def unquote(x):
	return x.replace('"','')

def parseAttributes(attr_string, value_sep):
	return dict([(x.strip().split(value_sep)[0], unquote(x.strip().split(value_sep)[1])) for x in attr_string.split(';')])

def parseAttributesDefault(attr_string, value_sep='='):
	return dict([(x.split(value_sep)[0], unquote(x.split(value_sep)[1])) for x in attr_string.split(';')])

def parseAttributesJGI1(attr_string):
	# name "e_gw1.1.1434.1"; proteinId 30297; exonNumber 2; product_name ""
	return parseAttributes(attr_string, value_sep=None)


class GFFRecord:
	def __init__( self, attribute_parser=parseAttributesDefault):
		self.seqname = ''
		self.source = ''
		self.feature = ''
		self.start = 0
		self.end = 0
		self.score = 0
		self.strand = '+'
		self.phase = '.'
		self._attributes_str = ''
		self._attr_parser = attribute_parser
		self.attributes = None

	def read(self, line):
		flds = line.strip().split('\t')
		self.readFromFields(flds)

	def readFromFields(self, flds):
		self.seqname = flds[0]
		self.source = flds[1]
		self.feature = flds[2]
		self.start = int(flds[3])
		self.end = int(flds[4])
		self.score = flds[5]
		self.strand = flds[6]
		self.phase = flds[7]
		if len(flds)>8:
			self._attributes_str = flds[8]
			self.attributes = self._attr_parser(self._attributes_str)

	def readFrom(self, flds):
		self.seqname = flds[0]
		self.source = flds[1]
		self.feature = flds[2]
		self.start = int(flds[3])
		self.end = int(flds[4])
		self.score = flds[5]
		self.strand = flds[6]
		self.phase = flds[7]
		self._attributes_str = flds[8]
		self.attributes = self._attr_parser(self._attributes_str)

	def copy(self):
		n = GFFRecord()
		n.seqname = self.seqname
		n.source = self.source
		n.feature = self.feature
		n.start = self.start
		n.end = self.end
		n.score = self.score
		n.strand = self.strand
		n.phase = self.phase
		n._attr_parser = self._attr_parser
		n._attributes_str = self._attributes_str
		n.attributes = dict(self.attributes.items())
		return n

	def toCodingExonRecord(self):
		cer = CodingExonRecord()
		cer.gene_ID = ''
		cer.transcript_ID = ''
		cer.peptide_ID = ''
		cer.exon_ID = ''
		cer.strand = self.strand
		cer.exon_start = self.start
		cer.exon_end = self.end
		cer.coding_start = self.start
		cer.coding_end = self.end
		cer.chromosome = self.seqname
		cer.addDescriptors(self.attributes)
		cer.gene_ID = cer.descriptors.get("name",'')
		return cer

	def getAttribute(self, key):
		return self.attributes.get(key,None)

	def __getitem__(self, key):
		return self.getAttribute(key)

	@property
	def attributes_string(self):
		return self._attributes_str

	def contains(self, index):
		# Does this record span this index, inclusive?
		return index >= self.start and index <= self.end

	def __str__(self):
		return "{x.seqname}\t{x.source}\t{x.feature}\t{x.start}\t{x.end}\t{x.score}\t{x.strand}\t{x.phase}\t{x._attributes_str}\n".format(x=self)

	def write(self, stream):
		stream.write(str(self))

class GFFRecordCollection:
	def __init__(self):
		self._collection = []

	def add(self, rec):
		self._collection.append(rec)

	def hits(self, index):
		hit_records = []
		for r in self._collection:
			if r.contains(index):
				hit_records.append(r)
		return hit_records

	def getCollection(self):
		return self._collection[:]

	collection = property(getCollection, None, None, None)

class GFFRecordTracker:
	def __init__(self, gff_rec):
		self.rec = gff_rec
		self.new_start = None
		self.new_end = None
		self._mutations = []

	def getOriginal(self):
		return self.rec

	def getUpdated(self):
		rec_copy = self.rec.copy()
		rec_copy.end = self.new_end
		rec_copy.start = self.new_start
		return rec_copy

	def atEnd(self, ref_site):
		return (self.rec.end-1 <= ref_site) and not (self.new_end is None) and not (self.new_start is None)

	def getMutations(self):
		return self._mutations[:]

	def update(self, ref_site, net_difference, mutation):
		if not mutation is None:
			self._mutations.append(mutation)
		# ref_site is 0-based index, start/end are 1-based indices delimiting ends of record
		# If ref_site is at or past the beginning of the record, and the start site hasn't been adjusted, adjust it
		if self.rec.start-1 <= ref_site and self.new_start is None:
			self.new_start = self.rec.start + net_difference
		elif self.rec.end-1 <= ref_site and self.new_end is None:
			# If ref_site is at or past the end of the record, and the end site hasn't been adjusted, adjust it
			# and return that we're at the end.
			self.new_end = self.rec.end + net_difference
		# Otherwise, do nothing -- changes to the net difference are meaningless.

'''
Col	Field	Description
1	CHROM	Chromosome name
2	POS	1-based position. For an indel, this is the position preceding the indel.
3	ID	Variant identifier. Usually the dbSNP rsID.
4	REF	Reference sequence at POS involved in the variant. For a SNP, it is a single base.
5	ALT	Comma delimited list of alternative seuqence(s).
6	QUAL	Phred-scaled probability of all samples being homozygous reference.
7	FILTER	Semicolon delimited list of filters that the variant fails to pass.
8	INFO	Semicolon delimited list of variant information.
9	FORMAT	Colon delimited list of the format of individual genotypes in the following fields.
10+	Sample(s)	Individual genotype information defined by FORMAT.

##fileformat=VCFv4.0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mdy748-raw-reads-cleaned-map.sorted.bam
ref|NC_001133|	19	.	CCCAC	CCACACCAC	99	.	INDEL;DP=1623;AF1=0.5;CI95=0.5,0.5;DP4=78,42,73,14;MQ=60;PV4=0.0026,1,1,0.011	PL:GT:GQ	255,0,255:0/1:99
ref|NC_001133|	24	.	ACACCACACCAC	ACCCACCACACCAC	11.8	.	INDEL;DP=1747;AF1=0.5;CI95=0.5,0.5;DP4=109,52,32,2;MQ=60;PV4=0.0012,1,1,1	PL:GT:GQ	49,0,255:0/1:52
ref|NC_001133|	25	.	CACCACA	CCA	99	.	INDEL;DP=1736;AF1=0.5;CI95=0.5,0.5;DP4=92,51,16,14;MQ=60;PV4=0.3,1,1,1	PL:GT:GQ	255,0,222:0/1:99
'''
class VCFRecord:
	def __init__( self ):
		self.chromosome = ''
		self.position = 0
		self.id = ''
		self.ref = 0
		self.alt = None
		self.qual = 0
		self.filter = '+'
		self.info = '.'
		self.format = ''
		self.other = None

	def read(self, line):
		flds = line.strip().split('\t')
		self.readFromFields(flds)

	def readFromFields(self, flds):
		self.chromosome = flds[0]
		self.position = int(flds[1])
		self.id = flds[2]
		self.ref = flds[3]
		self.alt = flds[4].split(',')
		self.qual = float(flds[5])
		self.filter = flds[6]
		self.info = flds[7]
		self.format = flds[8]
		if len(flds)>9:
			self.other = '/'.join(flds[9:])
	
	def getStrandRef(self, strand):
		res = self.ref
		if strand == '-':
			res = translate.reverseComplement(self.ref)
		return res

	def getStrandAlts(self, strand):
		res = self.alt
		if strand == '-':
			res = [translate.reverseComplement(a) for a in self.alt]
		return res

	def __str__(self):
		s = "{x.chromosome}\t{x.position}\t{x.id}\t{x.ref}\t{x.alt}\t{x.qual}\t{x.filter}\t{x.info}\t{x.format}".format(x=self)
		if not self.other is None:
			s += '\t{0}'.format(self.other)
		return s

	def write(self, stream):
		stream.write(str(self))

	def getReferenceStart(self):
		return self.position

	def getReferenceEnd(self):
		return self.position + len(self.ref)

	start = property(getReferenceStart, None, None, None)
	end = property(getReferenceEnd, None, None, None)
	
#-----------------------------------
# Multiple-alignment FASTA reader

def firstField(x):
	h = x.split()[0].strip()
	h = h.split('/')[0].strip()
	return h

def secondField(x):
	h = x.split()[1].strip()
	return h

def lastField(x):
	h = x.strip().split()[-1]
	return h

def secondOrFirstField(x):
	try:
		h = secondField(x)
	except IndexError:
		h = firstField(x)
	return h

class FASTAEntry(object):
	def __init__(self, header, sequence):
		self.header = header
		self.sequence = sequence

class UCSCExonHeader(object):
	def __init__(self, header):
		header_flds = header[1:].split()
		id_flds = header_flds[0].split('_')
		self.id = id_flds[0]
		self.species = id_flds[1]
		self.exon_total = id_flds[-1]
		self.exon_number = id_flds[-2]
		self.other_header_info = ' '.join(header_flds[2:])
		self.original_header = header

	def __str__(self):
		return self.exonHeader()
		
	def exonHeader(self):
		return self.original_header
	
	def CDSHeader(self):
		return "{h.id}_{h.species} {h.other_header_info}".format(h=self)

class NoExonFASTAHeader(object):
	def __init__(self, header, header_fxn=firstField):
		self.id = header_fxn(header[1:])
		self.exon_total = 1
		self.exon_number = 1
		self.other_header_info = header
		self.original_header = header
	
	def __str__(self):
		return self.exonHeader()

	def exonHeader(self):
		return self.original_header
	
	def CDSHeader(self):
		return self.original_header
		
class GeneSpeciesHeader(object):
	# Header of the style > gene_id species_id other_stuff
	def __init__(self, header):
		header_flds = header[1:].strip().split()
		self.id = header_flds[0]
		self.species = header_flds[1]
		self.exon_total = 1
		self.exon_number = 1
		self.other_header_info = ''
		if len(header_flds)>2:
			self.other_header_info = ' '.join(header_flds[2:])
		self.original_header = header.strip()[1:]
	
	def __str__(self):
		return self.exonHeader()

	def exonHeader(self):
		return self.original_header
	
	def CDSHeader(self):
		return self.original_header
		

class MultipleFASTAReader(object):
	"""Allows iteration over exons and CDSs from a multiple-alignment FASTA file.
	
	This implementation never loads the entire file into memory.
	"""
	def __init__(self, file_instance, header_class):
		self._file = file_instance
		self._header_class = header_class
		self._cache = util.LineCache(file_instance)
	
	def nextLine(self):
		line = self._cache.pop().strip()
		return line
		
	def atEnd(self):
		return self._cache.isEmpty()

	def exons(self):
		"""Iterates over a list of exon alignments"""
		id = None
		target_id = None
		while target_id == id and not self.atEnd():
			# Read header
			header = self.nextLine()
			while header == '' and not self.atEnd():
				header = self.nextLine()
			exon_list = []
			at_least_one_exon = False
			# A blank line separates exons; two blank lines separates CDS
			#while header == '':
			#	header = self.nextLine()
			while  target_id == id and header != '' and not self.atEnd():
				#print "${}^".format(header)
				head = self._header_class(header)
				id = head.id
				# Initialize target ID if we don't have one yet
				if target_id is None:
					target_id = head.id
				if id == target_id:
					seq = self.nextLine()
					if len(seq)>0:
						exon_list.append(FASTAEntry(head, seq))
						at_least_one_exon = True
					if not self.atEnd():
						header = self.nextLine()
				else:
					# Put back header
					#print "putting back", header
					self._cache.push(header)
					#print "end of cache", self._cache.cache[0]
					#print "end of cache pop", self._cache.pop().strip()
					#print "next line", self.nextLine()
			if at_least_one_exon and target_id == id:
				#if target_id != head.id:
				#	# Put back last line
				#	self._cache.push(header)
				yield exon_list
	
	def CDSs(self):
		"""Iterates over a list of coding sequences (CDSs)"""
		while not self.atEnd():
			cds = ''
			cds_list = []
			ex_alignment_lists = [x for x in self.exons()]
			if len(ex_alignment_lists)==0:
				continue
			#print len(ex_alignment_lists), len(ex_alignment_lists[0])
			sentinel_list = ex_alignment_lists[0]
			sentinel_entry = ex_alignment_lists[0][0]
			head = sentinel_entry.header
			#print head.id
			# First sort the CDS entries by exon number
			# DAD: assume that they are in sorted order for now
			# Then assemble CDS for each species by matching species ids
			species = [entry.header.species for entry in sentinel_list]
			species_headers = dict([(entry.header.species, entry.header) for entry in sentinel_list])
			species_cds = dict([(s,'') for s in species])
			for exon_list in ex_alignment_lists:
				for entry in exon_list:
					species_cds[entry.header.species] += entry.sequence
			cds_list = []
			for spec in species:
				cds = species_cds[spec]
				head = species_headers[spec]
				#new_header = '{h.id} {species} {h.other_header_info} {h.exon_total}'.format(h=head, species=spec)
				cds_list.append(FASTAEntry(head, cds))
			yield cds_list
	
	
		

#-----------------------------------

def writeFASTA(seq_list, outfile, headers=None):
	n_written = 0
	if headers is None:
		headers = ["{0}".format(j+1) for j in range(len(seq_list))]
	for (hdr,seq) in zip(headers,seq_list):
		line = ">{0}\n{1}\n".format(hdr, seq)
		outfile.write(line)
		n_written += 1
	return n_written

def writeFASTADict(seq_dict, outfile):
	headers = seq_dict.keys()
	seqs = [seq_dict[k] for k in headers]
	return writeFASTA(seqs, outfile, headers=headers)

#-----------------------------------------------------------------------------------
def readFASTA(infile):
	"""Reads the sequences and headers from a FASTA file.

	'infile' is a FASTA file.  Reads a list of all of the headers in the FASTA file
	and all of the sequences in the FASTA file and returns them as the tuple
	(headers, sequences) where headers[i] is the header corresponding to
	sequences[i].
	Removes the '>' from the headers.
	Returns None if there is a problem processing the file."""

	headers = []
	sequences = []
	for record in SeqIO.parse(infile, "fasta"):
		#headers.append("{:s} {:s}".format(record.id, record.description))
		headers.append("{:s}".format(record.description))
		sequences.append(str(record.seq))

	'''
	if isinstance(infile,str):
		infile_name = os.path.expanduser(infile)
		if not os.path.isfile(infile_name):
			raise BioFileError("Cannot find the FASTA file {}.".format(infile_name))
		else:
			infile = open(infile_name, 'r')
	seq = None
	headers = []
	sequences = []
	lines = infile.readlines()
	if len(lines)>0:
		for line in lines:  # read the lines from the file
			if line[0] == '#':
				# Skip comment
				continue
			if line[0] == '>':  # a new header
				if seq:
					frag = ''.join(seq)
					if frag[-1] == '*': # Remove trailing stop
						frag = frag[:-1]
					sequences.append(frag.upper())
				seq = []
				# Strip off leading '>'
				h = line[1:].strip()
				headers.append(h)
				#print h
				#print line
			else:
				frag = line.strip()
				if seq:
					seq.append(frag)
		frag = ''.join(seq)
		if frag[-1] == '*': # Remove trailing stop
			frag = frag[:-1]
		sequences.append(frag.upper())
	infile.close()
	assert len(headers) == len(sequences), "Error, headers and sequences have different lengths, {} != {}.".format(len(headers), len(sequences))
	'''
	return (headers, sequences)

def degapAlignment(seqs, gap='-'):
	"""
	Remove any column that contains a gap.
	"""
	ungapped_columns = []
	N = len(seqs)
	L = len(seqs[0])
	for i in range(L):
		col = [s[i] for s in seqs]
		if col.count(gap) < N:
			ungapped_columns.append(i)

	new_seqs = []
	for j in range(N):
		new_seqs.append(''.join([seqs[j][i] for i in ungapped_columns]))
	return new_seqs

def unalignAlignment(seqs, gap='-'):
	"""Remove gap characters in an existing alignment."""
	return [s.replace(gap,'') for s in seqs]

def getPeptideID(header):
	try:
		h = header.split("pep:")[1].split()[0].strip()
	except:
		h = firstField(header)
	return h

def getGeneID(header):
	try:
		h = header.split("gene:")[1].split()[0].strip()
	except:
		h = firstField(header)
	return h

def getTranscriptID(header):
	try:
		h = header.split("trans:")[1].split()[0].strip()
	except:
		h = firstField(header)
	return h

def getFlybaseGeneID(header):
	flds = header.split()
	id = flds[0]
	flybase_dict = dict([tuple(x[:-1].split('=')) for x in flds[1:]])
	#print flybase_dict
	gene_ids = [x for x in flybase_dict["parent"].split(',') if x.startswith('FBgn')]
	return gene_ids[0]

def getFlybaseTranscriptID(header):
	flds = header.split()
	id = flds[0]
	flybase_dict = dict([tuple(x[:-1].split('=')) for x in flds[1:]])
	gene_ids = [x for x in flybase_dict["parent"].split(',') if x.startswith('FBtr')]
	return gene_ids[0]

def getIDFunction(s):
	if s == 'pep' or s == 'peptide':
		return getPeptideID
	elif s == 'trans' or s == 'transcript':
		return getTranscriptID
	elif s == 'gene':
		return getGeneID
	elif s == 'flybase-gene':
		return getFlybaseGeneID
	elif s == 'flybase-transcript':
		return getFlybaseTranscriptID
	return firstField

def readFASTADict(infile_name, key_fxn = firstField):
	fdict = {}
	(headers, sequences) = readFASTA(infile_name)
	for (i,h) in enumerate(headers):
		h_key = key_fxn(h)
		s = sequences[i]
		fdict[h_key] = s
	return fdict

class GFFReader(object):
	def __init__(self, infile, attribute_parser=parseAttributesDefault):
		self._infile = infile
		self._attr_parser = attribute_parser
		if isinstance(infile, str):
			infile_name = os.path.expanduser(infile)
			if not os.path.isfile(infile_name):
				raise BioFileError("Cannot find the FASTA file {}.".format(infile_name))
			else:
				self._infile = open(infile_name, 'r')
		self._dlr = util.DelimitedLineReader(self._infile, header=False, sep='\t')

	def reader(self):
		return self._dlr

	@property
	def entries(self):
		# ##gff-version 3
		#<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
		gff_req_fields = ['seqname','source','feature','start','end','score','strand','frame']
		n_req_fields = len(gff_req_fields)
		gff_addl_fields = ['attributes','comments']
		for flds in self._dlr.entries:
			n = len(flds)
			if n >= n_req_fields:
				names = gff_req_fields + gff_addl_fields[:(n-n_req_fields)] # Only take additional fields if provided
				rec = GFFRecord(self._attr_parser)
				rec.readFromFields(flds)
				yield rec
