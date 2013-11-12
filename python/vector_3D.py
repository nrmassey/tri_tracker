#******************************************************************************
#** Program : vector_3D.py
#** Author  : Neil Massey
#** Date    : 26/07/13
#** Purpose : class to hold a 3D vector
#******************************************************************************

import numpy

class vector_3D:

	#**************************************************************************

	def __init__(self, xx=0, yy=0, zz=0):
		self.V = numpy.zeros(3, 'f');
		self.V[0] = xx
		self.V[1] = yy
		self.V[2] = zz

	#**************************************************************************

	def __getitem__(self, i):
		return self.V[i]
		
	#**************************************************************************
	
	def __setitem__(self, i, A):
		self.V[i] = A
	
	#**************************************************************************

	def __add__(self, rhs):
		new_vec = vector_3D()
		if isinstance(rhs, vector_3D):
			for i in range(0,3):
				new_vec.V[i] = self.V[i] + rhs.V[i]
		elif isinstance(rhs, (int, float)):
			for i in range(0,3):
				new_vec.V[i] = self.V[i] + rhs
		else:
			raise Exception("Unsupported right hand side type in operand")
		return new_vec

	#**************************************************************************

	def __sub__(self, rhs):
		new_vec = vector_3D()
		if isinstance(rhs, vector_3D):
			for i in range(0,3):
				new_vec.V[i] = self.V[i] - rhs.V[i]
		elif isinstance(rhs, (int, float)):
			for i in range(0,3):
				new_vec.V[i] = self.V[i] - rhs
		else:
			raise Exception("Unsupported right hand side type in operand")
		return new_vec

	#**************************************************************************

	def __mul__(self, rhs):
		new_vec = vector_3D()
		if isinstance(rhs, vector_3D):
			for i in range(0,3):
				new_vec.V[i] = self.V[i] * rhs.V[i]
		elif isinstance(rhs, (int, float)):
			for i in range(0,3):
				new_vec.V[i] = self.V[i] * rhs
		else:
			raise Exception("Unsupported right hand side type in operand")
		return new_vec

	#**************************************************************************

	def __div__(self, rhs):
		new_vec = vector_3D()
		if isinstance(rhs, vector_3D):
			for i in range(0,3):
				new_vec.V[i] = self.V[i] / rhs.V[i]
		elif isinstance(rhs, (int, float)):
			for i in range(0,3):
				new_vec.V[i] = self.V[i] / rhs
		else:
			raise Exception("Unsupported right hand side type in operand")
		return new_vec

	#**************************************************************************

	def __eq__(self, rhs):
		if not (isinstance(rhs, vector_3D)):
			return False
		eq = True
		for i in range(0,3):
			eq &= (self.V[i] == rhs[i])
		return eq

	#**************************************************************************

	def __ne__(self, rhs):
		if not (isinstance(rhs, vector_3D)):
			return True
		eq = True
		for i in range(0,3):
			eq &= (self.V[i] != rhs[i])
		return eq

	#**************************************************************************

	def dp(self, rhs):
		# dot product
		return self.V[0]*rhs.V[0] + self.V[1]*rhs.V[1] + self.V[2]*rhs.V[2]

	#**************************************************************************

	def xp(self, rhs):
		return vector_3D(self.V[1]*rhs.V[2] - self.V[2]*rhs.V[1],
						 self.V[2]*rhs.V[0] - self.V[0]*rhs.V[2],
						 self.V[1]*rhs.V[0] - self.V[1]*rhs.V[0])

	#**************************************************************************

	def mag(self):
		return numpy.sqrt(self.V[0]*self.V[0] + self.V[1]*self.V[1] + self.V[2]*self.V[2])

	#**************************************************************************

	def sqmag(self):
		return (self.V[0]*self.V[0] + self.V[1]*self.V[1] + self.V[2]*self.V[2])

	#**************************************************************************
