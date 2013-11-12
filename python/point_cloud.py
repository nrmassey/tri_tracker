#******************************************************************************
#** Program : point_cloud.py
#** Author  : Neil Massey
#** Date    : 26/07/13
#** Purpose : class to represent a point cloud and to load a point cloud in
#**			  from the binary mesh file
#******************************************************************************

from vector_3D import vector_3D
from bin_file_utils import *

class point_cloud:
	
	#******************************************************************************

	def __init__(self):
		self.point_cloud_store = []
	
	#******************************************************************************

	def __getitem__(self, i):
		return self.point_cloud_store[i]

	#******************************************************************************

	def load(self, fh):
		# read the number of points in the point cloud
		npts = read_int(fh)
		assert (npts > 0 and npts < 1e10)
		for i in range(0, npts):
			self.point_cloud_store.append(read_vector(fh))
