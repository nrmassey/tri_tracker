#! /usr/bin/env python

#########################################################################################
#
#  Program   : draw_tri_mesh.py
#  Author    : Neil Massey
#  Date      : 07/02/14
#  Purpose   : Draw triangular grid onto Buckminster Fuller Dymaxion (tm) grid
#
#########################################################################################

import sys, getopt
import matplotlib.pyplot as plt
from read_continents import read_continents
from dymax_proj import *
sys.path.append("/Users/massey/Coding/tri_tracker/python")
from tri_grid import *
from geo_convert import *
from scipy.io.netcdf import *
from conv_rot_grid import *

#########################################################################################

def draw_map(sp):
	# draw the continents on the map projection
	map = read_continents()
	for c in map:
		# loop through continents
		px = []
		py = []
		for p in c:
			px0, py0, tri, lcd = ll_to_dx(p[0], p[1])
			px.append(px0)
			py.append(py0)
		sp.fill(px, py, ec='k', fc='#dddddd', lw=0.25)
	
#########################################################################################

def draw_triangle(sp, tri_qn, max_level):
	# draw a triangle with the quad node at tri_qn
	# check whether this triangle has any children
	draw = False
	L = tri_qn.get_level() 
	if tri_qn.is_leaf():
		draw = True
	if L == max_level:
		draw = True
	
	tri_pairs = [(0,1), (1,2), (2,0)]
	edge_adjacents = tri_qn.get_data().get_adjacency_list(1)
	c = 'k'
	# height of a single triangle in outer mesh
	th = math.sqrt(0.75)
	if draw:
		if L == 0:
			return
		C = tri_qn.get_data().get_centroid()
		for p in range(0, 3):
			# don't draw if adjacent triangles at a higher level								
			TP = tri_pairs[p]
			P0 = tri_qn.get_data()[TP[0]]
			px0, py0, tri0, lcd0 = cart_to_dx(P0.V)
			P1 = tri_qn.get_data()[TP[1]]
			px1, py1, tri1, lcd1 = cart_to_dx(P1.V)

			# check for crossing between triangle 5 and 4
			if (tri0 != tri1):
				# Force the projection to use the triangle at which
				# the centroid of the triangle falls
				C = tri_qn.get_data().get_centroid()
				pxc, pyc, tric, lcdc = cart_to_dx(C.V)
				if (tri0 != tric):
					px0, py0, tri0, lcd0 = cart_to_dx(P0.V, tric)
				if (tri1 != tric):
					px1, py1, tri1, lcd1 = cart_to_dx(P1.V, tric)
			
			# does this edge share an edge in the adjacent triangles?
			for a in range(0, len(edge_adjacents)):
				adj_tn  = tg.get_node(edge_adjacents[a])
				adj_tri = adj_tn.get_data()
				for q in range(0, 2):
					r = q + 1
					if (P0 == adj_tri[q] and P1 == adj_tri[r]) or\
					   (P0 == adj_tri[r] and P1 == adj_tri[q]):
						if not adj_tn.is_leaf():
							draw = False
			# always draw if leaf node
			if L == max_level:
				draw = True
			# snap to grid
			if (L>=6):
				gh5 = th/(2**5) # grid size
				gh6 = th/(2**6)
				# snap to the level 5 grid for triangles  that are close to the level 5 grid
				py00 = numpy.floor([py0/gh5+0.5])[0] * gh5
				py10 = numpy.floor([py1/gh5+0.5])[0] * gh5
				if (abs(py00 - py0) <= gh6*0.4):
					py0 = py00
				if (abs(py10 - py1) <= gh6*0.4):
					py1 = py10
			else:
				gh = th/(2**L) # grid size
				py0 = numpy.floor([py0/gh+0.5])[0] * gh
				py1 = numpy.floor([py1/gh+0.5])[0] * gh
				

			# check for distances
			d = math.sqrt((py1-py0)**2 + (px1-px0)**2)
			TRI_L = 1.0 / (2**L) * 1.5
			if d > TRI_L:
				draw = False		

			if draw:
				sp.plot([px0,px1], [py0,py1], color=c, lw=0.5)
	else:
		# otherwise get the children and plot them
		for c in range(0, 4):
			T = tri_qn.get_child(c)
			draw_triangle(sp, T, max_level)

#########################################################################################

def draw_mesh(sp, tg, max_level):
	for i in range(0, len(tg.triangles)):
		T = tg.triangles[i].get_root()
		draw_triangle(sp, T, max_level)

#########################################################################################

def draw_outer_mesh():
	# 2D coordinates of the outer triangle after transformation
	# each triangle side length is 1
	# length and height of plot:
	lp = 5.5
	hp = math.sqrt(7.75)
	# height of a single triangle
	th = math.sqrt(0.75)
	
	# see 1st diagram at:
	# http://www.rwgrayprojects.com/rbfnotes/maps/graymap6.html
	# for point numbering
	
	# triangle points:
	tri_points = [(),				  # 0
	 			  ( 0.5,  0.0),       # 1
	              ( 1.5,  0.0),       # 2
	              ( 2.5,  0.0),       # 3
	              ( 3.5,  0.0),       # 4
	              ( 4.5,  0.0),       # 5
	              ( 0.75, 0.5*th),    # 6
	              ( 1.5,  2.0/3.0*th),# 7
	              ( 2.0,  1.0/3.0*th),# 8
	              ( 0.0,  th),        # 9
	              ( 1.0,  th),		  # 10
	              ( 2.0,  th),		  # 11
	              ( 3.0,  th),        # 12
	              ( 4.0,  th),        # 13
	              ( 5.0,  th),        # 14
	              ( 5.5,  th),        # 15
	              
	              ( 0.5, 2*th),		  # 16
	              ( 1.5, 2*th),       # 17
	              ( 2.5, 2*th),       # 18
	              ( 3.5, 2*th),       # 19
	              ( 4.5, 2*th),       # 20
	              ( 5.5, 2*th),       # 21
	              
	              ( 1.0, 3*th),       # 22
	              ( 2.0, 3*th),       # 23
	              ( 3.0, 3*th),       # 24
	              ( 4.0, 3*th),       # 25
	              ( 5.0, 3*th),       # 26
                 ]
    
    # triangle index lists
	tri_lists  = [(23, 18, 17),		  #  1
                  (17, 18, 11),       #  2
                  (18, 12, 11),       #  3
                  (18, 19, 12),       #  4
                  (23, 24, 18),       #  5
                  (22, 23, 17),       #  6
                  (16, 17, 10),       #  7
                  (10, 17, 11),       #  8
                  (11,  3,  8),       #  9a
                  (10, 11,  7, 2),    #  9b
                  (11, 12,  3),       # 10
                  (12, 13,  4),       # 11
                  (19, 13, 19),       # 12
                  (19, 20, 13),       # 13
                  (25, 20, 19),       # 14
                  (26, 21, 20),       # 15
                  (21, 15, 14),       # 16a
                  ( 9, 10,  6),       # 16b
                  ( 1, 10,  2),       # 17
                  ( 4, 13,  5),       # 18
                  (13, 20, 14),       # 19
                  (20, 21, 14),       # 20
                 ]
	for tri in tri_lists:
		for p in range(0, len(tri)):
			I0 = tri[p]
			P0 = tri_points[I0]
			px0, py0 = P0[0], P0[1]
	
			s = p+1
			if s >= len(tri):
				s = 0		
			I1 = tri[s]
			P1 = tri_points[I1]
			px1, py1 = P1[0], P1[1]
			sp.plot([px0,px1], [py0,py1], color='k', lw=1.5)
			
#	for p in range (1, len(tri_points)):
#		P = tri_points[p]
#		px0, py0 = P[0], P[1]
#		sp.text(px0, py0, str(p))
			
#########################################################################################

def draw_original_grid(sp, orig_file, orig_var):
    nc_file = netcdf_file(orig_file, "r")
    nc_var = nc_file.variables[orig_var]
    # get rotated pole attribute and coordinates
    rot_pole_var_name = nc_var._attributes["grid_mapping"]
    pole_lon = nc_file.variables[rot_pole_var_name]._attributes["grid_north_pole_longitude"]
    pole_lat = nc_file.variables[rot_pole_var_name]._attributes["grid_north_pole_latitude"]
    
    # get the rotated coordinates
    lat_name = nc_var.dimensions[2]		# CF compliance should always have lat / lon
    lon_name = nc_var.dimensions[3]	    # in this order
    lats = nc_file.variables[lat_name][:]
    lons = nc_file.variables[lon_name][:]
    
    # loop through
    for rlat in lats[::2]:
    	for rlon in lons[::2]:
    		# convert from rotated grid to lat lon grid and then to dymaxion grid
    		glon, glat = rot2glob(rlon, rlat, pole_lon, pole_lat)
    		dx, dy, tri, lcd = ll_to_dx(glon, glat)
    		sp.plot(dx, dy, 'r.', ms=2.0, zorder=3)
    nc_file.close()

#########################################################################################

if __name__ == "__main__":
	opts, args = getopt.getopt(sys.argv[1:], "m:o:l:x:y:a:g:v:",
							   ['mesh_file=','out_name=','level=','x_limits=', 
							    'y_limits=','max_level=','orig_file=','orig_var='])
	mesh_file = ""
	l = 0
	out_name = ""
	x_limits = [0,5.5]
	y_limits = [-0.1,math.sqrt(7.75)]
	for opt, val in opts:
		if opt in ['--mesh_file', '-m']:
			mesh_file = val
		if opt in ['--out_name',  '-o']:
			out_name = val
		if opt in ['--level',     '-l']:
			l = int(val)
		if opt in ['--max_level', '-a']:
			max_level = int(val)
		if opt in ['--y_limits',  '-y']:
			y = val.strip("[]").split(",")
			y_limits = [float(y[0]), float(y[1])]
		if opt in ['--x_limits',  '-x']:
			x = val.strip("[]").split(",")
			x_limits = [float(x[0]), float(x[1])]
		if opt in ['--orig_file', '-g']:
			orig_file = val
		if opt in ['--orig_var',  '-v']:
			orig_var = val	

	
	# initialise the Dymaxion map
	init()
 	# load the triangular grid in
 	tg = tri_grid()
 	fh = open(mesh_file, 'rb')
 	tg.load(fh)
 	fh.close()
 	# create the plotting surface
 	sp = plt.subplot(111)
 
 	# draw the triangles
 	draw_mesh(sp, tg, max_level)
 	# draw original twelve triangles
 	draw_outer_mesh()
 	# draw original grid
	draw_original_grid(sp, orig_file, orig_var)
	# draw the continents
	draw_map(sp)
	# remove the axis
	sp.axis('off')
	# set the map size
	sp.set_xlim(x_limits)
	sp.set_ylim(y_limits)
	# set the aspect to 1
	sp.set_aspect(1.0)
	plt.tight_layout()
	plt.savefig(out_name)
