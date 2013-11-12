#******************************************************************************
#** Program : data_store.py
#** Author  : Neil Massey
#** Date    : 31/07/13
#** Purpose : class to represent the data storage used in the triangular mesh
#**			  regridding
#******************************************************************************

from bin_file_utils import *
import numpy
from meta_data import *

class data_store:

	#******************************************************************************

	def __init__(self):
		self.ds = numpy.zeros([0], 'f')
		self.mv = 2e20
		self.meta_data = {}
		
	#******************************************************************************

	def __getitem__(self, I):
		I = list(I)
		for i in I:
			if i >= self.ds.shape[1]:
				I.remove(i)
		if len(I) < 2:
			return self.mv
		return self.ds[I[0], I[1]]
		
	#******************************************************************************

	def __setitem__(self, I, V):
		I = list(I)
		for i in I:
			if i >= self.ds.shape[1]:
				I.remove(i)
		if len(I) < 2:
			return

		self.ds[I[0], I[1]] = V
		
	#******************************************************************************

	def set_size(self, n_t_steps, n_idxs):
		self.ds.resize([n_t_steps, n_idxs])

	#******************************************************************************
	
	def get_missing_value(self):
		return self.mv

	#******************************************************************************

	def get_meta_data(self):
		return self.meta_data

	#******************************************************************************

	def get_n_t_steps(self):
		return self.ds.shape[0]

	#******************************************************************************

	def get_n_idxs(self):
		return self.ds.shape[1]

	#******************************************************************************

	def max(self, time_step=-1, idx_list=[]):
		# idx_list allows you to specify which indices to take the max over - so
		# you can take it over just a single level in the tri-grid (for example)
		if idx_list == []:
			idx_list = [i for i in range(0, len(self.ds))]
		for i in idx_list:
			if i >= self.ds.shape[1]:
				idx_list.remove(i)
		if time_step == -1:
			D = self.ds[:, idx_list]
			return numpy.max(D[numpy.abs(D) < 0.99*abs(self.mv)])
		else:
			T = self.ds[time_step]
			D = T[:, idx_list]
			# require some fault tolerance
			return numpy.max(D[numpy.abs(D) < 0.99*abs(self.mv)])

	#******************************************************************************
	
	def min(self, time_step=-1, idx_list=[]):
		# idx_list allows you to specify which indices to take the max over - so
		# you can take it over just a single level in the tri-grid (for example)
		if idx_list == []:
			idx_list = [i for i in range(0, len(self.ds))]
		for i in idx_list:
			if i >= self.ds.shape[1]:
				idx_list.remove(i)
		if time_step == -1:
			D = self.ds[:, idx_list]
			return numpy.min(D[numpy.abs(D) < 0.99*abs(self.mv)])
		else:
			T = self.ds[time_step]
			D = T[:, idx_list]
			return numpy.min(D[numpy.abs(D) < 0.99*abs(self.mv)])
	
	#******************************************************************************

	def load(self, fh):
		self.meta_data = read_meta_data(fh)
		n_t_steps = read_int(fh)
		assert(n_t_steps < 1e7)
		n_idxs = read_int(fh)
		assert(n_idxs < 1e7)		# fewer than 1 million pts - check that we are
									# reading the correct file
		self.mv = read_float(fh)
		self.ds = numpy.zeros([n_t_steps, n_idxs], 'f')
		# there must be a faster way of reading than two nested loops
		for t in range(0, n_t_steps):
			for i in range(0, n_idxs):
				self.ds[t, i] = read_float(fh)
				
	#******************************************************************************

	def save(self, fh):
		# just the reverse of the above procedure
		write_int(fh, self.get_n_t_steps())
		write_int(fh, self.get_n_idxs())
		write_float(fh, self.mv)
		for t in range(0, self.get_n_t_steps()):
			for i in range(0, self.get_n_idxs()):
				write_float(fh, self.ds[t, i])