import os, random, string, sys, math, traceback, unittest
import time, datetime, collections
#import pytz
import na

class listdict(collections.abc.MutableMapping):
	def __init__(self, thedict=None):
		self._dict = {}
		if not thedict is None:
			for (k,v) in thedict.items():
				self._dict[k] = [v]
	
	def __getitem__(self, key):
		res = []
		try:
			res = self._dict[key]
		except KeyError:
			self._dict[key] = res
		return res
	
	def __setitem__(self, key, value):
		if isinstance(value,list):
			self._dict[key] = value
		else:
			self._dict[key] = [value]
	
	def __delitem__(self, key):
		del self._dict[key]
	
	def __len__(self):
		return len(self._dict)
	
	def __iter__(self):
		for it in self._dict:
			yield it

def shuffle_string(seq):
	rs = [x for x in seq]
	random.shuffle(rs)
	return ''.join(rs)

# Float equality.  On my system (WinXP, Python 2.6), smallest distinguishable float difference is 7.45e-9.
def feq(f1,f2,eps=1e-8):
	return abs(f1-f2)<eps

# Easy timestamp production
def timestamp(zone='US/Central', timeformat='%Y-%m-%d %H:%M:%S %Z%z'):
	#tz = pytz.timezone(zone)
	#localtime = tz.localize(datetime.datetime.fromtimestamp(time.time()))
	timest = datetime.datetime.fromtimestamp(time.time()).strftime(timeformat)
	#timest = localtime.strftime(timeformat)
	return timest

def printTiming(func):
	"""Use as follows. Given a function foo(arg1,arg2), define as:
	@util.printTiming
	def foo(arg1,arg2):
		<some stuff>
	When foo(a,b) is called, printTiming will execute.
	"""
	# Define a wrapper function for func which stores the time before and after calling func.
	def wrapper(*arg):
		# Store the current time, in microseconds
		t1 = time.time()
		# Call the function
		res = func(*arg)
		# Store the current time again, in microseconds
		t2 = time.time()
		print('{0} took {1:0.3f} ms'.format(func.__name__, (t2-t1)*1000.0)) #, file=sys.stdout)
		return res
	# Return the wrapper function
	return wrapper

class OutStreams:
	def __init__(self, stream_list=[]):
		if not isinstance(stream_list, list):
			stream_list = [stream_list]
		self.streams = stream_list[:]

	def write(self, line, flush=False):
		for outs in self.streams:
			#print(line, file=outs, flush=flush)
			outs.write(line)

	def addStream(self, stream):
		self.streams.append(stream)

	def removeStream(self, stream):
		self.streams.remove(stream)

	def flush(self):
		for outs in self.streams:
			outs.flush()

def isNA(x):
	sys.stderr.write("util.isNA() should be updated to na.isNA()")
	return na.isNA(x)

""" Class for applying formatting automatically.
"""
class FieldFormatter:
	def __init__(self, var, format, transform=lambda x: x):
		self.var = var
		self.format = format
		self.transform = transform

	def __str__(self):
		res = None
		if not na.isNA(self.var):
			try:
				trans_var = self.transform(self.var)
				res = self.format.format(trans_var)
			except ValueError:
				pass
			except TypeError:
				pass
		else:
			res = na.NA
		return res

def dictParser(x, entry_sep=';', key_sep='='):
	entries = x.split(entry_sep)
	return dict([entry.split(key_sep) for entry in entries])

#########################
# DelimitedLineReader
#########################

def isBlank(s):
	return s.strip() == ''

def isComment(s, comment_char='#'):
	return s.strip().startswith(comment_char)

def looseIntParser(x):
	v = None
	try:
		v = int(x)
	except ValueError:
		if not na.isNA(x):
			v = naFloatParser(x)
	return v

def naIntParser(x):
	v = None
	try:
		v = int(x)
	except ValueError as ve:
		if not na.isNA(x):
			raise ve
	return v

def naFloatParser(x):
	v = None
	try:
		v = float(x)
	except ValueError as ve:
		if not na.isNA(x):
			raise ve
	return v

def naSciParser(x):
	v = None
	try:
		v = float(x)
	except ValueError as ve:
		if not na.isNA(x):
			raise ve
	return v

def naStringParser(x):
	"""A parser that respects NA's."""
	v = None
	if not na.isNA(x):
		v = str(x)
	return v

def maxQuantHeader(header_flds):
	# Header line
	header_line = '\t'.join(header_flds)
	header_line = header_line.lower()
	header_line = header_line.replace("h/l", "hl")
	header_line = header_line.replace("h/m", "hm")
	header_line = header_line.replace("m/l", "ml")
	header_line = header_line.replace(" ", ".")
	header_line = header_line.replace("(", ".")
	header_line = header_line.replace(".[%]", "")
	header_line = header_line.replace("[", ".")
	header_line = header_line.replace("]", "")
	header_line = header_line.replace(")", "")
	header_line = header_line.replace("/", '.')
	header_line = header_line.replace("+", "plus")
	header_line = header_line.replace(".\t", '\t')
	header_line = header_line.replace(".\n", '\n')
	header_line = header_line.replace("..", '.')
	header_flds = header_line.strip().split('\t')
	return header_flds

def defaultHeader(header_flds):
	return header_flds


class LineCache:
	"""Class for caching lines."""
	def __init__(self, instream, comment_str='#', save_comments=False):
		self.cache = []
		self.instream = instream
		self.cache_size = 100 # lines
		self.comment_str = comment_str
		self.comment_cache = None
		self._save_comments = save_comments
		if save_comments:
			self.comment_cache = []
		self.refill()

	def add(self, line):
		self.cache.append(line)

	def push(self, line):
		self.cache = [line] + self.cache

	def pop(self):
		# Use standard Python queue pattern
		res = self.cache.pop(0)
		if len(self.cache) == 0:
			self.refill()
		return res

	def refill(self):
		# Refill the cache
		# Presently only guarantees reading self.cache_size lines -- if these are all comments, no data will be read.
		eof = False
		for li in range(self.cache_size):
			line = self.instream.readline()
			# If end of file, bail.
			if not line:
				eof = True
				break
			else:
				if not isComment(line, self.comment_str): #line.strip().startswith(self.comment_str):
					if  not line.strip()=='':
						self.add(line)
				else:
					if self._save_comments:
						self.comment_cache.append(line)
		return eof
	
	def popcomment(self):
		res = None
		if len(self.comment_cache) > 0:
			res = self.comment_cache.pop()
		return res

	def getLine(self, index):
		eof = False
		while not eof and (index >= len(self.cache)):
			eof = self.refill()
		if index >= len(self.cache):
			raise ReaderEOFError("Attempt to read past end of stream ({:d}, {:d})".format(index, len(self.cache)))
		return self.cache[index]

	def isEmpty(self):
		return len(self.cache) == 0

	def len(self):
		return len(self.cache)

class ReaderError(Exception):
	"""Exception class for Readers"""

class ReaderEOFError(ReaderError):
	"""Error when Reader goes past end of file"""

def basicHeaderFixer(flds):
	# Make headers unique
	res = flds
	if len(flds) > len(set(flds)):
		new_flds = ['']*len(flds)
		for fi in range(len(flds)):
			f = flds[fi]
			f_count = flds.count(f)
			if f_count == 1 or fi == 0:
				new_flds[fi] = f
			else:  # more than one instance
				i = 1
				new_fld_name = '{:s}.{:d}'.format(f,i)
				# only look at field names up to this one.
				# if this is the first instance of a repeated field name,
				# it will not be changed.
				while new_fld_name in set(new_flds + flds[0:fi-1]):
					i += 1
					new_fld_name = '{:s}.{:d}'.format(f,i)
				new_flds[fi] = new_fld_name
		# Ensure we've actually fixed the problem
		assert len(new_flds) == len(set(new_flds))
		res = new_flds
	return res

class DelimitedLineReader:
	"""
	Class which parses delimited files line-by-line.
	Typical usage:
	dlr = DelimitedLineReader(file(...,'r'), delim='\t')
	headers = dlr.getHeader()
	for fields in dlr.entries:
	    print fields[1]
	# or
	for fields in dlr.dictentries:
		print fields['foo'] + fields['bar']
	"""
	handler_dict = {"s":str, "f":naFloatParser, "d":naIntParser, "e":naSciParser}

	def __init__(self, in_file, header=True, field_defs=None, sep="\t", skip=0, strip=False, comment_str="#", save_comments=False, custom_handler_dict=None, header_name_processor=basicHeaderFixer):
		self.infile = in_file
		self.delim = sep
		self.strip = strip
		self.skip = skip
		self.comment_str = comment_str
		# Data
		self.cache = LineCache(self.infile, comment_str=self.comment_str)
		self.comment_cache = None
		self.cur_flds = None
		self.cur_line = None
		self.n_lines_read = 0
		self.has_header = header
		self.headers = None
		self.field_defs = field_defs
		self.handlers = []
		self.header_name_processor = header_name_processor
		
		# Move past lines to skip
		for nskip in range(0,self.skip):
			self.next(process=False)

		if self.has_header:
			# Position reader right after header
			self.getHeader(move_to_data=True)

		if field_defs is None:
			# Attempt to infer the handlers
			self.field_defs = self.inferHandlers()
		else:
			for h in self.field_defs:
				try:
					self.handlers.append(self.handler_dict[h])
				except KeyError: # unrecognized field definition
					if not custom_handler_dict is None:
						try:
							self.handlers.append(custom_handler_dict[h])
						except KeyError as ke:
							raise ReaderError("No custom handler provided for field-type {:s}".format(h))

	def next(self, process=True, apply_handlers=True):
		if not self.atEnd():
			self.cur_line = self.cache.pop()
			self.n_lines_read += 1
		else:
			raise ReaderEOFError("Attempt to read past end of stream")
		res = None
		if self.isValid():
			res = self.cur_line
			if process: # Split the line into fields, parse by types according to handlers if specified
				self.cur_flds = self.process(apply_handlers)
				res = self.cur_flds
		return res

	def nextDict(self, apply_handlers=True):
		## Return fields as a dictionary, keyed by the header names
		## If no headers, return dictionary keyed by column numbers
		headers = self.getHeader()
		res = None
		flds = self.next(True, apply_handlers)
		if not headers is None:
			res = dict(zip(headers, flds))
		else:
			res = dict(zip(range(len(flds)),flds))
			#raise ReaderError, "Attempt to return dictionary of fields, but header is empty"
		return res

	def process(self, apply_handlers=True):
		res = None
		if self.isValid():
			res = self.processLine(self.cur_line, apply_handlers)
		return res

	def processLine(self, line, apply_handlers=True):
		if self.strip:
			line = line.strip()
		if line != "":
			flds = line.split(self.delim)
			flds[-1] = flds[-1].strip() # Get rid of \n
			if apply_handlers:
				if len(flds) > len(self.handlers):
					raise ReaderError("Line {} has {} fields, but only {} expected".format(self.n_lines_read, len(flds), len(self.handlers)))
				assert len(flds) <= len(self.handlers)
				if not self.handlers:
					# Infer handlers are strings
					self.handlers = [str for i in range(len(flds))]
				done_processing = False
				while not done_processing:
					res = []
					for hi in range(len(flds)):
						try:
							# Apply the new handlers to get the data.
							res.append(self.handlers[hi](flds[hi]))
						except ValueError:
							# Adaptively update handlers if there was a value error.
							handler_key = self.inferHandlerKey(flds[hi])
							self.handlers[hi] = self.handler_dict[handler_key]
							res.append(self.handlers[hi](flds[hi]))
							# Reapply the new handlers to get the data.
							#res = [self.handlers[i](flds[i]) for i in range(len(flds))]
					done_processing = True
			else:
				res = flds
		else:
			res = None
		return res

	def isValid(self):
		return not self.cur_line is None

	def atEnd(self):
		return self.cache.isEmpty()

	@property
	def ncol(self):
		res = None
		if self.handlers:
			res = len(self.handlers)
		return res

	def getRawLine(self):
		return self.cur_line

	def get(self):
		return self.cur_flds

	def getField(self, field):
		if self.headers:
			res = self.cur_flds[self.headers.index(field)]
		else:
			res = self.cur_flds[field]
		return res

	def getNumRead(self):
		return self.n_lines_read

	def isComment(self):
		res = False
		if self.isValid():
			#print self.cur_line, self.cache.len()
			line = self.cur_line.strip()
			if len(line)>0:
				res = line[0] == self.comment_str
		return res

	def isBlank(self):
		res = False
		if self.isValid():
			#print self.cur_line, self.cache.len()
			line = self.cur_line.strip()
			res = isBlank(line)
		return res

	def getHeaders(self, move_to_data=True):
		return self.getHeader(move_to_data)
	
	def getHeader(self, move_to_data=True):
		res = self.headers
		if res is None and self.has_header:
			li = 0
			self.cur_line = self.cache.getLine(li)
			while self.isComment():
				li += 1
				self.cur_line = self.cache.getLine(li)
			if self.isValid() and not self.isComment():
				self.headers = self.process(apply_handlers=False)
				if move_to_data:
					# Position reader at line following header
					for pop_li in range(li+1):
						self.cache.pop()
						self.n_lines_read += 1
			else:
				raise ReaderError("Never encountered a header line")
			# Now process the headers -- e.g., de-dupe column names.
			self.headers = self.header_name_processor(self.headers)
			res = self.headers
		return res

	def inferHandlerKey(self, fld):
		# int -> float -> string
		found_key = None
		for handler_key in "dfse":  # DAD: should include custom handlers up front
			handler = self.handler_dict[handler_key]
			try:
				res = handler(fld)
				# If we get here without an exception, save the key
				found_key = handler_key
				break
			except ValueError:
				continue
		return found_key


	def inferHandlers(self, max_lines=100):
		# DAD: run through fields until we've seen at least one non-NA for each.
		handlers_identified = False
		li = 0
		self.cur_line = self.cache.getLine(li)
		self.handlers = None
		inferred_string = []
		while not handlers_identified and li < max_lines and self.isValid():
			if not self.isComment() and not self.isBlank():
				# Not a comment line -- parse it.
				if self.strip:
					self.cur_line = self.cur_line.strip()
				flds = self.cur_line.split(self.delim)
				flds[-1] = flds[-1].strip() # Get rid of \n
				# Initialize empty handler list if we haven't done so already
				if self.handlers is None:
					self.handlers = [None]*len(flds)
					inferred_string = ['X']*len(flds)
				#if len(flds) != len(self.handlers):
				#	print flds
				assert len(flds) == len(self.handlers), "Number of fields {} not equal to number of handlers {}".format(len(flds), len(self.handlers))
				for hi in range(len(self.handlers)):
					fld = flds[hi]
					if self.handlers[hi] is None:
						if not na.isNA(fld):
							handler_key = self.inferHandlerKey(fld)
							inferred_string[hi] = handler_key
							self.handlers[hi] = self.handler_dict[handler_key]
					else: # handler has already been found; just confirm, and upgrade if necessary
						try:
							val = self.handlers[hi](fld)
						except ValueError:
							#print "upgrading handler", inferred_string[hi],
							handler_key = self.inferHandlerKey(fld)
							inferred_string[hi] = handler_key
							self.handlers[hi] = self.handler_dict[handler_key]
							#print "to", handler_key

				# We're finished when all handlers are not None.
				handlers_identified = len([h for h in self.handlers if h is None]) == 0
			if not handlers_identified:
				li += 1
				try:
					self.cur_line = self.cache.getLine(li)
				except ReaderEOFError:
					# We've reached the end of the file with an inconclusive result -- some fields
					# still can't have types inferred.
					# Just assume everything's a string.
					for hi in range(len(self.handlers)):
						if self.handlers[hi] is None:
							self.handlers[hi] = self.handler_dict["s"]
					handlers_identified = True
		if not handlers_identified and li >= max_lines:
			# Went past the allowed number of lines to look ahead; set all unset handlers to strings
			for hi in range(len(self.handlers)):
				if self.handlers[hi] is None:
					self.handlers[hi] = self.handler_dict["s"]
		#print inferred_string
		inferred_string = ''.join(inferred_string)
		return inferred_string

	def setHandlerType(self, handler_index, type_string):
		try:
			self.handlers[handler_index] = self.handler_dict[type_string]
		except KeyError as ke:
			raise ReaderError("Unknown handler type {}".format(type_string))
		except IndexError:
			raise ReaderError("Bad handler index {}".format(handler_index))
	
	def setColumnType(self, column_name, type_string):
		headers = self.getHeader()
		if not headers is None:
			try:
				column_id = headers.index(column_name)
				self.setHandlerType(column_id, type_string)
			except ValueError:
				raise ReaderError("Column {} not in header".format(column_name))
		else:
			# Try to interpret as column index
			try:
				column_id = int(column_name)
				if column_id < len(self.handlers):
					self.setHandlerType(column_id, type_string)
			except ValueError:
				raise ReaderError("Column {} not in header".format(column_name))
	
	def setHeaderNames(self, names):
		headers = self.getHeader()
		assert(len(names)==len(headers))
		self.headers = names
	
	# Generators for iteration
	@property
	def entries(self):
		while not self.atEnd():
			yield self.next()

	@property
	def dictentries(self):
		while not self.atEnd():
			yield self.nextDict()


# Read

class LightDataFrame:
	def __init__(self, header, data):
		self._header = basicHeaderFixer(header)
		# Data must be a list of lists, each list is a column.
		self._data = data
		self._header_lookup = dict([(h, self._header.index(h)) for h in self._header])

	def row(self, row_index):
		return [d[row_index] for d in self._data]

	def rowDict(self, row_index):
		return dict(zip(self._header,self.row(row_index)))

	def col(self, column_key):
		"""Index by column, given by column_key"""
		col_index = self._header_lookup[column_key]
		return self.coli(col_index)
	
	def coli(self, column_index):
		return self._data[column_index]

	def __getitem__(self, column_key):
		"""Index by column, given by column_key"""
		col_index = self._header_lookup[column_key]
		return self._data[col_index]

	@property
	def header(self):
		return self._header

	@property
	def rows(self):
		for ri in range(self.nrows):
			yield self.row(ri)

	@property
	def dictrows(self):
		for ri in range(self.nrows):
			yield self.rowDict(ri)

	@property
	def nrows(self):
		return len(self._data[0])

	@property
	def ncols(self):
		return len(self._data)

	@property
	def headers(self):
		return self._headers[:]
	
	def __str__(self):
		s = "\t".join(self.header)+"\n"
		for ri in range(self.nrows):
			rd = self.rowDict(ri)
			for h in self.header:
				s += '\t' + str(rd[h])
			s += '\n'
		return s


def readTable(stream, header=True, sep='\t', skip=0, header_name_processor=defaultHeader):
	dlr = DelimitedLineReader(stream, header=header, sep=sep, skip=skip, header_name_processor=header_name_processor)
	data_dict = {}
	initialized = False
	for flds in dlr.dictentries:
		if not initialized:
			if header:
				header_flds = dlr.getHeader()
			else:
				header_flds = range(len(flds))
			for h in header_flds:
				data_dict[h] = [flds[h]]
			initialized = True
		else:
			for h in header_flds:
				data_dict[h].append(flds[h])
	return LightDataFrame(header_flds, [data_dict[h] for h in header_flds])


class Header(object):
	format_types = {'f':'float', 'd':'int', 's':'string', 'e':'float'}
	def __init__(self, name, description, format):
		self.name = name
		self.description = description
		self.format = format[-1].lower()
		try:
			self.type_desc = Header.format_types[self.format]
		except KeyError:
			self.type_desc = ''


class DelimitedOutput(object):
	def __init__(self, sep='\t'):
		self._sep = sep
		self._header_list = []
	
	def addHeader(self, header, description='', format='s'):
		self._header_list.append(Header(header, description, format))

	def addHeaders(self, headers, descriptions=None):
		if not descriptions is None:
			assert(len(descriptions) == len(headers))
			for (h,d) in zip(headers, descriptions):
				self.addHeader(h,d)
		else:
			for h in headers:
				self.addHeader(h)

	def writeHeader(self, stream):
		line = self._sep.join([h.name for h in self._header_list]) + '\n'
		stream.write(line)
	
	def describeHeader(self, stream, comment='#'):
		stream.write("{}\n".format(comment))
		for h in self._header_list:
			line = "{com}\t{h.name} [{h.type_desc}] = {h.description}\n".format(com=comment, h=h)
			stream.write(line)

	def _makeFormattable(self, str):
		for ch in ' .#$%\t()':
			str = str.replace(ch,'_')
		return str

	def getFormat(self, named=True, force_string=False):
		#print "{name:s}\:{fmt:s}".format(name='a',fmt='b')
		formatfxn = lambda x: x
		if force_string:
			formatfxn = lambda x: 's'
		if named:
			res = self._sep.join(["{{{name:s}:{fmt:s}}}".format(name=self._makeFormattable(h.name), fmt=formatfxn(h.format)) for h in self._header_list])+'\n'
		else:
			res = self._sep.join(["{{:{fmt:s}}}".format(fmt=formatfxn(h.format)) for h in self._header_list])+'\n'
		return res

	def formatLine(self, entry):
		""" Format the entry by headers. """
		line = ""
		tab = '\t'
		# If entry is a dictionary...
		for h in self._header_list:
			try:
				fmt = "{:"+h.format+'}'
				if line != '':
					line += tab
				line += fmt.format(entry[h.name])
			except ValueError as ve:
				if h.format=='s': # string
					line += fmt.format(naStringParser(entry[h.name]))
				else:
					raise ve
			except TypeError as te:
				# This will occur if entry[h.name] is None. Format as NA.
				val = entry[h.name]
				res = na.formatNA(val)
				fmt = "{:s}"
				line += fmt.format(res)
		return line + '\n'

	def createResult(self, default=None):
		""" Generate a new result dictionary with default values populated. """
		entry = dict([(h.name, default) for h in self._header_list])
		return entry

	@property
	def headers(self):
		for h in self._header_list:
			yield h

