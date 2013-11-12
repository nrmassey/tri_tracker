#******************************************************************************
#** Program : quad_tree.py
#** Author  : Neil Massey
#** Date    : 26/07/13
#** Purpose : class to hold quad tree to hold the triangular mesh information
#******************************************************************************

from types import *

#******************************************************************************

class quad_tree_node:

	#******************************************************************************

	def __init__(self, parent=None, data=None):
		self.children = [None, None, None, None]
		self.nc = 0
		self.parent = parent
		self.data = data

	#******************************************************************************
	
	def add_child(self, data=None):
		if self.nc >= 4:
			raise Exception("Cannot add any more children to node")
		else:
			self.children[self.nc] = quad_tree_node(parent=self, data=data)
			self.nc += 1
		return self.children[self.nc-1]
			
	#******************************************************************************

	def get_children(self):
		return self.children	
	
	#******************************************************************************

	def get_child(self, child_number):
		return self.children[child_number]

	#******************************************************************************
		
	def get_parent(self):
		return self.parent
	
	#******************************************************************************
	
	def get_data(self):
		return self.data

	#******************************************************************************
	
	def set_data(self, data):
		self.data = data
		
	#******************************************************************************
	
	def is_leaf(self):
		isleaf = True
		for i in range(0, 4):
			isleaf = isleaf & isinstance(self.children[i], NoneType)
		return isleaf

	#******************************************************************************
		
	def get_level(self):
		# get the depth the node occurs at
		l = 0
		cn = self.parent;
		while (not isinstance(cn, NoneType)):
			l+=1
			cn = cn.parent
		return l
		
	#*****************************************************************************
		
	def get_max_level(self, node):
		# get the maximum number of levels under this node
		if isinstance(node, NoneType):
			return 0;
		else:
			# get the max depth of the children
			mD = [0,0,0,0]
			for i in range(0, 4):
				mD[i] = self.get_max_level(node.get_child(i))
			mD.sort()
			return mD[3]+1;

#******************************************************************************

class quad_tree:

	#******************************************************************************

	def __init__(self, data=None):
		self.root = quad_tree_node(data=data)

	#******************************************************************************

	def get_root(self):
		return self.root

	#******************************************************************************

	def get_max_level(self):
		return self.root.get_max_level(self.root)

	#******************************************************************************

	def	get_all_nodes_at_level(self, level):
		node_list = []
		self.__get_nodes_at_level(self.root, level, node_list)
		return node_list

	#******************************************************************************
				
	def __get_nodes_at_level(self, current_node, level, node_list):
		# check whether the current node equals the level
		if current_node.get_level() == level:
			node_list.append(current_node)

		# recursively get the child nodes
		for i in range(0,4):
			c_node = current_node.get_child(i)
			if not isinstance(c_node, NoneType):
				self.__get_nodes_at_level(c_node, level, node_list)
