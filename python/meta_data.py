#******************************************************************************
#** Program : meta_data.py
#** Author  : Neil Massey
#** Date    : 12/11/13
#** Purpose : Load in the meta data from a pure binary file as a python
#**           dictionary
#******************************************************************************

from bin_file_utils import *
import struct

def read_meta_data(fh):
	# return the metadata as a dictionary
	meta_data = {}
	b = fh.read(4)
	# if the file contains metadata
	if b == "META":
		# read the number of metadata articles
		n_meta_items = read_int(fh)
		for i in range(0, n_meta_items):
			key = read_string(fh)
			value = read_string(fh)
			meta_data[key] = value
	else:
		# rewind to beginning of file
		fh.seek(0)
	return meta_data