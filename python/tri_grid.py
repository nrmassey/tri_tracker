#******************************************************************************
#** Program : tri_grid.py
#** Author  : Neil Massey
#** Date    : 26/07/13
#** Purpose : class to represent a triangular grid and to load the binary
#**           file format created by gen_grid
#******************************************************************************

from point_cloud import *
from indexed_tri_3D import *
from quad_tree import *
from StringIO import *
from meta_data import *

class tri_grid:

	#**************************************************************************

	def __init__(self):
		self.triangles = []
		self.meta_data = {}
	
	#**************************************************************************

	def get_node(self, label):
		# get the triangle node from the label.
		# the label is now completely numeric and has the following format
	
		# two least significant figures : parent triangle (0->19)
		# next sig fig : level 1 triangle (100,200,300,400)
		# next sig fig : level 2 triangle (1000,2000,3000,4000) etc.
		# Add these together to get the label: 2318, 1100, etc.
	
		# first two powers of ten are the first triangle
		tri_number = label[0] % 100;
		# loop through the rest of the label, taking the requisite child
		current = self.triangles[tri_number].get_root()
		# these two values keep track of which significant figure we are at
		mp1 = 1000;
		mp2 = 100;
		for i in range(0, label[1]):
			child_number = (label[0] % mp1)/mp2
			current = current.get_child(child_number-1)
			assert(not isinstance(current, NoneType))
			mp1 *= 10
			mp2 *= 10
		return current

	#**************************************************************************

	def get_triangle(self, label):
		return self.get_node(label).get_data()

	#**************************************************************************
	
	def load(self, fh):
		# read the grid in from the binary file format
		self.meta_data = read_meta_data(fh)
		# create and read the point cloud
		self.point_cloud_instance = point_cloud()
		self.point_cloud_instance.load(fh)
		# read the number of base triangles in
		n_tris = read_int(fh);
		assert(n_tris < 1e7)		# less than a million triangles - this is
									# just a check that we are reading the 
									# correct file
		# load the triangles
		current_tri = indexed_tri_3D(self.point_cloud_instance)
		current_tri.load(fh)
		for i in range(0, n_tris):
			self.triangles.append(quad_tree())
			current_node = self.triangles[i].get_root()
			current_tri  = self.__load_node(fh, current_node, current_tri)

	#**************************************************************************

	def get_triangles_at_level(self, level):
		# get the first triangle to create a list
		tri_node_list = []
		# now extend the list to contain the extra triangles
		for i in range(0, len(self.triangles)):
			this_tri_list = self.triangles[i].get_all_nodes_at_level(level)
			tri_node_list.extend(this_tri_list)
		return tri_node_list
		
	#**************************************************************************

	def get_ds_indices_at_level(self, level):
		# get the ds (datastore) indices for a level in the mesh
		tris = self.get_triangles_at_level(level)
		ds_idxs = []
		for tri in tris:
			ds_idxs.append(tri.get_data().get_ds_index())
		return ds_idxs

	#**************************************************************************
	
	def get_meta_data(self):
		return self.meta_data
	
	#**************************************************************************

	def __load_node(self, fh, current_node, current_tri):
		# assign the current data
		current_node.set_data(current_tri)
		# get the next triangle in the file
		next_tri = indexed_tri_3D(self.point_cloud_instance)
		# try to load
		if not next_tri.load(fh):
			return next_tri
		# label[1] contains the maximum level of the triangle
		# check the current triangle label against the previous triangle label
		if next_tri.get_label()[1] == current_tri.get_label()[1]:
			# occur at same level so do not create any children
			return next_tri
			
		if next_tri.get_label()[1] > current_tri.get_label()[1]:
			# do not occur at same level so add children and start the recursion
			for i in range(0,4):
				# add child and assign the next triangle to it
				current_node = current_node.add_child(indexed_tri_3D(self.point_cloud_instance))
				next_tri = self.__load_node(fh, current_node, next_tri)
				# need to go back up to the parent node so as to create the siblings correctly
				# otherwise we would be reproducing with our own siblings
				if not isinstance(current_node.get_parent(), NoneType):
					current_node = current_node.get_parent()
		# control will never reach here but will return next_tri anyway to avoid warnings
		return next_tri

	#**************************************************************************
