#******************************************************************************
#** Program : bin_file_utils.py
#** Author  : Neil Massey
#** Date    : 26/07/13
#** Purpose : utilities to read in various data types from a binary file
#******************************************************************************

from vector_3D import vector_3D
import struct

#******************************************************************************

def read_vector(fh):
	# get the data
	b = fh.read(12)
	X = struct.unpack('f', b[0:4])[0]
	Y = struct.unpack('f', b[4:8])[0]
	Z = struct.unpack('f', b[8:12])[0]
	new_vec = vector_3D(X, Y, Z)
	return new_vec

#******************************************************************************

def write_vector(fh, V):
	# pack the data into bytes
	fh.write(struct.pack('f', V[0], V[1], V[2]))

#******************************************************************************
	
def read_string(fh):
	# get the string length first
	b = fh.read(4)
	strlen = struct.unpack('i', b)[0]
	b = fh.read(strlen)
	return b

#******************************************************************************

def write_string(fh, S):
	# write the string length first
	fh.write(struct.pack('B', strlen(S)))
	fh.write(S)

#******************************************************************************

def read_float(fh):	
	b = fh.read(4)
	return struct.unpack('f', b)[0]

#******************************************************************************

def write_float(fh, F):
	# write a floating point number out as binary
	fh.write(struct.pack('f', F))

#******************************************************************************

def read_int(fh):
	b = fh.read(4)
	return struct.unpack('i', b)[0]

#******************************************************************************

def write_int(fh, I):
	fh.write(struct.pack('i', I))

#******************************************************************************

def read_int_as_byte(fh):
	b = fh.read(1)
	return struct.unpack('b', b)[0]
	
#******************************************************************************

def write_int_as_byte(fh, I):
	fh.write(struct.pack('b', I))
	
#******************************************************************************

def read_label(fh):
	# get the data first - long int for the label, so 8 bytes
	L = fh.read(8)
	label = struct.unpack('l', L)
	# single byte for max level
	b = fh.read(1)
	max_level = struct.unpack('b', b)
	return [label[0], max_level[0]]
