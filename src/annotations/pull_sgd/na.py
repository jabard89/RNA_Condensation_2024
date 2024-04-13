import math

# Definition of the default NA string.
NA = 'NA'

def isNA(x):
	res = x is None
	if not res:
		if type(x) == str:
			x = x.upper()
			res = (x==NA or x=='' or x=='NAN')
		elif type(x) == float:
			res = math.isnan(x)
	return res

def strNA(x):
	if isNA(x):
		return NA
	else:
		return str(x)

def formatNA(x, format="{0}", sep=None):
	if isinstance(x, list):
		flds = format.split(sep)
		if sep is None:
			sep = '\t'
		assert len(flds) == len(x)
		res = []
		for i in range(len(x)):
			if isNA(x[i]):
				res.append(NA)
			else:
				res.append(flds[i].format(x[i]))
		return sep.join(res)
	else:
		if isNA(x):
			return NA
		else:
			return format.format(x)
