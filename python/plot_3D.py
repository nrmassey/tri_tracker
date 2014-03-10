#!/usr/bin/env python

###############################################################################
# Program : plot_3D.py
# Author  : Neil Massey
# Date	  : 08/07/13
# Purpose : Program to plot the grid from the regional version of tri_tracker
#           as a 3D solid
###############################################################################

import sys, os, getopt

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from plot_tri_tracker import load_mesh_file

###############################################################################

def project_3D_to_2D(P):
	# project a 3D point to 2D
	d = (2 - P[2])
	x = P[0] / d
	y = P[1] / d
	return [x,y]

###############################################################################

def plot_mesh(sp, tg, l, tl=-1):
	print "Plotting mesh"
	if tl == -1:
		tl = l
	# just plot the mesh at the level from the tri_grid
	# get the triangles
	tri_list = tg.get_triangles_at_level(l)
	# loop over each triangle
	for t in tri_list:
		if t.is_leaf() or t.get_level() == tl:
			# get each corner of the triangle and convert to model coordinates (lat/lon)
			TRI = t.get_data()
			P = [project_3D_to_2D(TRI[0]), project_3D_to_2D(TRI[1]), project_3D_to_2D(TRI[2])]
			# don't draw if it goes around the date line
			if abs(P[0][0] - P[1][0]) > 180 or abs(P[1][0] - P[2][0]) > 180 or abs(P[2][0] - P[1][0]) > 180:
				continue
			# get centroid
			C = TRI.centroid
			fc = [1,1,1,1]
			ec = [0,0,0,1]
			sp.fill([P[0][0], P[1][0], P[2][0], P[0][0]], \
					[P[0][1], P[1][1], P[2][1], P[0][1]], \
					facecolor=fc, edgecolor=ec, lw=0.25, zorder=C[2])

###############################################################################

if __name__ == "__main__":
	mesh_file   = ""
	opts, args = getopt.getopt(sys.argv[1:], "m:o:l:",
							   ['mesh_file=', 'out_name=', 'level='])
	l = 0
	for opt, val in opts:
		if opt in ['--mesh_file', '-m']:
			mesh_file = val
		if opt in ['--out_name', '-o']:
			out_name = val
		if opt in ['--level', '-l']:
			l = int(val)

	sp = plt.subplot(111)			
	tg = load_mesh_file(mesh_file)
	plot_mesh(sp, tg, l)
	fig = plt.gcf()
	fig.set_size_inches(16,12)
	plt.savefig(out_name, dpi=100)	
#	plt.show()