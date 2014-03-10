#!/usr/bin/env python

###############################################################################
# Program : measure_angles_in_iso.py
# Author  : Neil Massey
# Date	  : 03/03/14
# Purpose : Program to measure the angles in the triangles of the original
#           icosahedron, to ensure that all 3 points are specified in a
#           clockwise direction
###############################################################################

import sys, os, getopt

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy
import numpy.linalg

from plot_tri_tracker import load_mesh_file
import math

if __name__ == "__main__":
	mesh_file   = ""
	opts, args = getopt.getopt(sys.argv[1:], "m:o:l:",
							   ['mesh_file='])
	l = 0
	for opt, val in opts:
		if opt in ['--mesh_file', '-m']:
			mesh_file = val
	tg = load_mesh_file(mesh_file)

	tri_list = tg.get_triangles_at_level(l)
	
	for t in range(0, len(tri_list)):
		M = numpy.zeros([3,3], 'f')
		TRI = tri_list[t].get_data()
		C = TRI.centroid
		# calculate surface normal
		N = (TRI[1]-TRI[0]).xp(TRI[2]-TRI[1])
		print t, N.dp(C)