#! /usr/bin/env python

#########################################################################################
#
#  Program   : dymax_proj.py
#  Author    : Neil Massey
#  Date      : 06/02/14
#  Purpose   : Project triangular grid to Buckminster Fuller's Dymaxion (tm) grid
#  Reference : Code adapted from Robert W Gray's website:
#              http://www.rwgrayprojects.com/rbfnotes/maps/graymap1.html
#
#              Robert W Gray's papers on the Dymaxion Map:
#              Gray, Robert W., Fuller's DymaxionTM Map, Cartography and Geographic Information Systems, 21(4): 243-246, 1994.
# 			   Gray, Robert W., Exact Transformation Equations For Fuller's World Map, Cartographica, 32(3): 17-25, 1995.
#
#########################################################################################

import numpy, math

# global variables
V = numpy.zeros([13, 3], 'f')
C = numpy.zeros([21, 3], 'f')
garc = 0.0
gt   = 0.0
gdve = 0.0
gel  = 0.0

#########################################################################################

def init():
	# initializes the global variables which includes the 
	# vertex coordinates and mid-face coordinates.        

	# global variables
	global V
	global C
	global garc
	global gt
	global gdve
	global gel 

	# Cartesian coordinates for the 12 vertices of icosahedron

	V[1,:]  = [ 0.420152426708710003,  0.078145249402782959,  0.904082550615019298]
	V[2,:]  = [ 0.995009439436241649, -0.091347795276427931,  0.040147175877166645]
	V[3,:]  = [ 0.518836730327364437,  0.835420380378235850,  0.181331837557262454]
	V[4,:]  = [-0.414682225320335218,  0.655962405434800777,  0.630675807891475371]
	V[5,:]  = [-0.515455959944041808, -0.381716898287133011,  0.767200992517747538]
	V[6,:]  = [ 0.355781402532944713, -0.843580002466178147,  0.402234226602925571]
	V[7,:]  = [ 0.414682225320335218, -0.655962405434800777, -0.630675807891475371]
	V[8,:]  = [ 0.515455959944041808,  0.381716898287133011, -0.767200992517747538]
	V[9,:]  = [-0.355781402532944713,  0.843580002466178147, -0.402234226602925571]
	V[10,:] = [-0.995009439436241649,  0.091347795276427931, -0.040147175877166645]
	V[11,:] = [-0.518836730327364437, -0.835420380378235850, -0.181331837557262454]
	V[12,:] = [-0.420152426708710003, -0.078145249402782959, -0.904082550615019298]

	## now calculate mid face coordinates

	H = (V[1] + V[2] + V[3]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[1,:] = H / magn

	H = (V[1] + V[3] + V[4]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[2,:] = H / magn

	H = (V[1] + V[4] + V[5]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[3,:] = H / magn

	H = (V[1] + V[5] + V[6]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[4,:] = H / magn

	H = (V[1] + V[2] + V[6]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[5,:] = H / magn

	H = (V[2] + V[3] + V[8]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[6,:] = H / magn

	H = (V[8] + V[3] + V[9]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[7,:] = H / magn

	H = (V[9] + V[3] + V[4]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[8,:] = H / magn

	H = (V[10] + V[9] + V[4]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[9,:] = H / magn

	H = (V[5] + V[10] + V[4]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[10,:] = H / magn

	H = (V[5] + V[11] + V[10]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[11,:] = H / magn

	H = (V[5] + V[6] + V[11]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[12,:] = H / magn

	H = (V[11] + V[6] + V[7]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[13,:] = H / magn

	H = (V[7] + V[6] + V[2]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[14,:] = H / magn

	H = (V[8] + V[7] + V[2]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[15,:] = H / magn

	H = (V[12] + V[9] + V[8]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[16,:] = H / magn

	H = (V[12] + V[9] + V[10]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[17,:] = H / magn

	H = (V[12] + V[11] + V[10]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[18,:] = H / magn

	H = (V[12] + V[11] + V[7]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[19,:] = H / magn

	H = (V[12] + V[8] + V[7]) / 3.0
	magn = math.sqrt(H[0]**2 + H[1]**2 + H[2]**2)
	C[20,:] = H / magn

	garc = 2.0 * math.asin( math.sqrt( 5 - math.sqrt(5)) / math.sqrt(10) )
	gt = garc / 2.0;

	gdve = math.sqrt( 3 + math.sqrt(5) ) / math.sqrt( 5 + math.sqrt(5) )
	gel = math.sqrt(8) / math.sqrt(5 + math.sqrt(5))

#########################################################################################

def s_tri_info(P, tri=-1):
	# Get which of the twenty triangles the Cartesian coordinate (x,y,z) is in
	# P = numpy array of (x,y,z)
	
	# variables required
	h_dist = numpy.zeros([3], 'f')
	h = numpy.zeros([3], 'f')

	h_tri = 0;
	h_dist[0] = 9999.0;

	# Which triangle face center is the closest to the given point
	# is the triangle in which the given point is in.

	if tri != -1:
		h_tri = tri
	else:
		for i in range(1, 20):	
			h = C[i] - P;
			h_dist[1] = math.sqrt(h[0]**2 + h[1]**2 + h[2]**2)
			if (h_dist[1] < h_dist[0]):
				h_tri = i
				h_dist[0] = h_dist[1]

	# Now the LCD triangle is determined.
	if h_tri == 1:
		v1 =  1; v2 =  3; v3 =  2
	elif h_tri == 2:
		v1 =  1; v2 =  4; v3 =  3
	elif h_tri == 3:
		v1 =  1; v2 =  5; v3 =  4 
	elif h_tri == 4:
		v1 =  1; v2 =  6; v3 =  5
	elif h_tri == 5:
		v1 =  1; v2 =  2; v3 =  6 
	elif h_tri == 6:
		v1 =  2; v2 =  3; v3 =  8 
	elif h_tri == 7:
		v1 =  3; v2 =  9; v3 =  8 
	elif h_tri == 8:
		v1 =  3; v2 =  4; v3 =  9 
	elif h_tri == 9:
		v1 =  4; v2 = 10; v3 =  9 
	elif h_tri == 10:
		v1 =  4; v2 =  5; v3 = 10 
	elif h_tri == 11:
		v1 =  5; v2 = 11; v3 = 10 
	elif h_tri == 12:
		v1 =  5; v2 =  6; v3 = 11 
	elif h_tri == 13:
		v1 =  6; v2 =  7; v3 = 11 
	elif h_tri == 14:
		v1 =  2; v2 =  7; v3 =  6 
	elif h_tri == 15:
		v1 =  2; v2 =  8; v3 =  7 
	elif h_tri == 16:
		v1 =  8; v2 =  9; v3 = 12 
	elif h_tri == 17:
		v1 =  9; v2 = 10; v3 = 12 
	elif h_tri == 18:
		v1 = 10; v2 = 11; v3 = 12 
	elif h_tri == 19:
		v1 = 11; v2 =  7; v3 = 12 
	elif h_tri == 20:
		v1 =  8; v2 = 12; v3 =  7 

	h = P - V[v1]
	h_dist[0] = math.sqrt(h[0]**2 + h[1]**2 + h[2]**2)

	h = P - V[v2];
	h_dist[1] = math.sqrt(h[0]**2 + h[1]**2 + h[2]**2)

	h = P - V[v3];
	h_dist[2] = math.sqrt(h[0]**2 + h[1]**2 + h[2]**2)

	if ( (h_dist[0] <= h_dist[1]) and (h_dist[1] <= h_dist[2]) ):
		h_lcd = 1
	if ( (h_dist[0] <= h_dist[2]) and (h_dist[2] <= h_dist[1]) ):
		h_lcd = 6
	if ( (h_dist[1] <= h_dist[0]) and (h_dist[0] <= h_dist[2]) ):
		h_lcd = 2
	if ( (h_dist[1] <= h_dist[2]) and (h_dist[2] <= h_dist[0]) ):
		h_lcd = 3
	if ( (h_dist[2] <= h_dist[0]) and (h_dist[0] <= h_dist[1]) ):
		h_lcd = 5
	if ( (h_dist[2] <= h_dist[1]) and (h_dist[1] <= h_dist[0]) ):
		h_lcd = 4

	return h_tri, h_lcd

#########################################################################################

def conv_ll_t_sc(lng, lat):
	# convert (long., lat.) point into spherical polar coordinates
	# with r=radius=1.  Angles are given in radians.              

	h_theta = 90.0 - lat ;
	h_phi = lng;
	if (lng < 0.0):
		h_phi = lng + 360.0
	theta = math.radians(h_theta)
	phi = math.radians(h_phi)
	return theta, phi

#########################################################################################

def s_to_c(theta, phi):
	# Covert spherical polar coordinates to cartesian coordinates.
	# The angles are given in radians.

	v = numpy.zeros([3], 'f')

	v[0] = math.sin(theta) * math.cos(phi)
	v[1] = math.sin(theta) * math.sin(phi)
	v[2] = math.cos(theta)

	return v

#########################################################################################

def c_to_s(P):
	# convert cartesian coordinates into spherical polar coordinates.
	# The angles are given in radians.                                

	if (P[0]>0.0 and P[1]>0.0):
		a = math.radians(0.0)
	if (P[0]<0.0 and P[1]>0.0):
		a = math.radians(180.0)
	if (P[0]<0.0 and P[1]<0.0):
		a = math.radians(180.0)
	if (P[0]>0.0 and P[1]<0.0):
		a = math.radians(360.0)
		
	lat = math.acos(P[2]);
	
	if (P[0]==0.0 and P[1]>0.0):
		lng = math.radians(90.0)
	if (P[0]==0.0 and P[1]<0.0):
		lng = math.radians(270.0)
	if (P[0]>0.0  and P[1]==0.0):
		lng = math.radians(0.0)
	if (P[0]<0.0  and P[1]==0.0):
		lng = math.radians(180.0)
	if (P[0]!=0.0 and P[1]!=0.0):
		lng = math.atan(P[1]/P[0]) + a
		
	return lng, lat

#########################################################################################

def r2(axis, alpha, P):
	# Rotate a 3-D point about the specified axis.

	Q = numpy.zeros([3], 'f')
	if (axis == 1):
		Q[0] = P[0]
		Q[1] = P[1] * math.cos(alpha) + P[2] * math.sin(alpha)
		Q[2] = P[2] * math.cos(alpha) - P[1] * math.sin(alpha)

	if (axis == 2):
		Q[0] = P[0] * math.cos(alpha) - P[2] * math.sin(alpha)
		Q[1] = P[1]
		Q[2] = P[0] * math.sin(alpha) + P[2] * math.cos(alpha)

	if (axis == 3):
		Q[0] = P[0] * math.cos(alpha) + P[1] * math.sin(alpha)
		Q[1] = P[1] * math.cos(alpha) - P[0] * math.sin(alpha)
		Q[2] = P[2]
	return Q

#########################################################################################

def rotate(angle, x, y):
	# Rotate the point to correct orientation in XY-plane.

	ha = math.radians(angle)
	px = x * math.cos(ha) - y * math.sin(ha);
	py = x * math.sin(ha) + y * math.cos(ha);
	
	return px, py

#########################################################################################


def dymax_point(tri, lcd, P):
	# converts a Cartesian coordinate (in P) to a Dymaxion point
	# arrays to hold staging variables
	h0 = numpy.zeros([3], 'f')
	h1 = numpy.zeros([3], 'f')

	# In order to rotate the given point into the template spherical
	# triangle, we need the spherical polar coordinates of the center 
	# of the face and one of the face vertices. So set up which vertex 
	# to use.                                                          

	if tri == 1:
		v1 = 1
	elif tri == 2:
		v1 =  1
	elif tri == 3:
		v1 =  1
	elif tri == 4:
		v1 =  1
	elif tri == 5:
		v1 =  1
	elif tri == 6:
		v1 =  2
	elif tri == 7:
		v1 =  3
	elif tri == 8:
		v1 =  3
	elif tri == 9:
		v1 =  4
	elif tri == 10:
		v1 =  4
	elif tri == 11:
		v1 =  5
	elif tri == 12:
		v1 =  5
	elif tri == 13:
		v1 =  6
	elif tri == 14:
		v1 =  2
	elif tri == 15:
		v1 =  2
	elif tri == 16:
		v1 =  8
	elif tri == 17:
		v1 =  9
	elif tri == 18:
		v1 = 10
	elif tri == 19:
		v1 = 11
	elif tri == 20:
		v1 =  8

	h0[:] = P[:] # ensure a copy rather than reference

	h1[:] = V[v1,:]

	hlng, hlat = c_to_s(C[tri]);

	axis = 3
	h0 = r2(axis, hlng, h0)
	h1 = r2(axis, hlng, h1)

	axis = 2
	h0 = r2(axis, hlat, h0)
	h1 = r2(axis, hlat, h1)

	hlng, hlat = c_to_s(h1);
	hlng = hlng - math.radians(90.0)

	axis = 3
	h0 = r2(axis, hlng, h0)
	
	## exact transformation equations

	gz = math.sqrt(1 - h0[0]**2 - h0[1]**2)
	gs = math.sqrt(5 + 2 * math.sqrt(5) ) / ( gz * math.sqrt(15) )

	gxp = h0[0] * gs ;
	gyp = h0[1] * gs ;

	ga1p = 2.0 * gyp / math.sqrt(3.0) + (gel / 3.0)
	ga2p = gxp - (gyp / math.sqrt(3)) +  (gel / 3.0)
	ga3p = (gel / 3.0) - gxp - (gyp / math.sqrt(3))

	ga1 = gt + math.atan( (ga1p - 0.5 * gel) / gdve)
	ga2 = gt + math.atan( (ga2p - 0.5 * gel) / gdve)
	ga3 = gt + math.atan( (ga3p - 0.5 * gel) / gdve)

	gx = 0.5 * (ga2 - ga3)

	gy = (1.0 / (2.0 * math.sqrt(3)) ) * (2 * ga1 - ga2 - ga3)

	# Re-scale so plane triangle edge length is 1.

	x = gx / garc;
	y = gy / garc;

	# rotate and translate to correct position

	if tri == 1:
		x,y = rotate(240.0, x, y)
		px = x + 2.0
		py = y + 7.0 / (2.0 * math.sqrt(3.0))
		
	elif tri == 2:
		x,y = rotate(300.0, x, y)
		px = x + 2.0
		py = y + 5.0 / (2.0 * math.sqrt(3.0))
		
	elif tri == 3:
		x,y = rotate(0.0, x, y)
		px = x + 2.5
		py = y + 2.0 / math.sqrt(3.0)
		
	elif tri == 4: 
		x,y = rotate(60.0, x, y)
		px = x + 3.0;
		py = y + 5.0 / (2.0 * math.sqrt(3.0))
		
	elif tri == 5:
		x,y = rotate(180.0, x, y)
		px = x + 2.5
		py = y + 4.0 * math.sqrt(3.0) / 3.0
		
	elif tri == 6:
		x,y = rotate(300.0, x, y)
		px = x + 1.5
		py = y + 4.0 * math.sqrt(3.0) / 3.0
		
	elif tri == 7:
		x,y = rotate(300.0, x, y);
		px = x + 1.0
		py = y + 5.0 / (2.0 * math.sqrt(3.0))
		
	elif tri == 8:
		x,y = rotate(0.0, x, y)
		px = x + 1.5
		py = y + 2.0 / math.sqrt(3.0)
		
	elif tri == 9:
		if (lcd > 2):
			x,y = rotate(300.0, x, y)
			px = x + 1.5
			py = y + 1.0 / math.sqrt(3.0)
		else:
			x,y = rotate(0.0, x, y)
			px = x + 2.0
			py = y + 1.0 / (2.0 * math.sqrt(3.0))

	elif tri == 10:
		x,y = rotate(60.0, x, y)
		px = x + 2.5
		py = y + 1.0 / math.sqrt(3.0)
		
	elif tri == 11:
		x,y = rotate(60.0, x, y)
		px = x + 3.5
		py = y + 1.0 / math.sqrt(3.0)
		
	elif tri == 12:
		x,y = rotate(120.0, x, y)
		px = x + 3.5
		py = y + 2.0 / math.sqrt(3.0)
		
	elif tri == 13:
		x,y = rotate(60.0, x, y)
		px = x + 4.0
		py = y + 5.0 / (2.0 * math.sqrt(3.0))
		
	elif tri == 14:
		x,y = rotate(0.0, x, y)	
		px = x + 4.0
		py = y + 7.0 / (2.0 * math.sqrt(3.0))
		
	elif tri == 15:
		x,y = rotate(0.0, x, y)
		px = x + 5.0
		py = y + 7.0 / (2.0 * math.sqrt(3.0))
		
	elif tri == 16:
		if (lcd < 4):
			x,y = rotate(60.0, x, y)
			px = x + 0.5
			py = y + 1.0 / math.sqrt(3.0)
		else:
			x,y = rotate(0.0, x, y)
			px = x + 5.5
			py = y + 2.0 / math.sqrt(3.0)
			
	elif tri == 17:
		x,y = rotate(0.0, x, y)
		px = x + 1.0
		py = y + 1.0 / (2.0 * math.sqrt(3.0))
		
	elif tri == 18:
		x,y = rotate(120.0, x, y)
		px = x + 4.0
		py = y + 1.0 / (2.0 * math.sqrt(3.0))
		
	elif tri == 19:
		x,y = rotate(120.0, x, y)
		px = x + 4.5
		py = y + 2.0 / math.sqrt(3.0)
		
	elif tri == 20:
		rotate(300.0, x, y);
		px = x + 5.0
		py = y + 5.0 / (2.0 * math.sqrt(3.0))

	return px, py

#########################################################################################

def ll_to_dx(lon, lat):
	t, p = conv_ll_t_sc(lon, lat)
	v = s_to_c(t, p)
	tri, lcd = s_tri_info(v)
	px, py = dymax_point(tri, lcd, v)
	return px, py, tri, lcd
	
#########################################################################################

def cart_to_dx(v, tri=-1):
	tri, lcd = s_tri_info(v, tri)
	px, py = dymax_point(tri, lcd, v)
	return px, py, tri, lcd
