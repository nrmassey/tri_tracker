#******************************************************************************
#** Program : indexed_tri_3D.py
#** Author  : Neil Massey
#** Date    : 26/07/13
#** Purpose : class to hold a 3D triangle with index into a rectangular grid
#**			  created by gen_grid
#******************************************************************************

import numpy
from bin_file_utils import *

#******************************************************************************

class indexed_tri_3D:

	#**************************************************************************
	
	def __init__(self, ipoint_cloud_instance):
		self.pc_idx = numpy.zeros(3, 'i')
		self.point_cloud_instance = ipoint_cloud_instance
		self.grid_indices = []
		self.adjacency = [[],[]]
		
	#**************************************************************************

	def get_label(self):
		return self.label
		
	#**************************************************************************

	def get_ds_index(self):
		return self.ds_index
		
	#**************************************************************************

	def get_grid_indices(self):
		return self.grid_indices
		
	#**************************************************************************

	def get_adjacency_list(self, type):
		return self.adjacency[type]
		
	#**************************************************************************

	def __getitem__(self, i):
		return self.point_cloud_instance[self.pc_idx[i]]

	#**************************************************************************

	def get_centroid(self):
		return self.centroid

	#**************************************************************************

	def __calculate_centroid(self):
		self.centroid = (self[0] + self[1] + self[2])/3.0

	#**************************************************************************

	def load(self, fh):
		try:
			# read the label
			self.label = read_label(fh)
			# read the indices into the point cloud
			for i in range(0,3):
				self.pc_idx[i] = read_int(fh)
			
			# calculate the centroid
			self.__calculate_centroid();
			# read the target index
			self.ds_index = read_int(fh)
			# read the length of index list in
			n_idx = read_int(fh)
			for i in range(0, n_idx):
				i = read_int(fh)
				j = read_int(fh)
				cart_coord = read_vector(fh)
				self.grid_indices.append((i,j,cart_coord))
			# read the point and adjacency list
			for a in range(0, 2):
				# get the size first
				s = read_int(fh)
				for i in range(0, s):
					self.adjacency[a].append(read_label(fh))
			return True
		except:
			return False
