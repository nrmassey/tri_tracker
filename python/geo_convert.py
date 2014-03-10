#******************************************************************************
#** Program : geo_convert.py
#** Author  : Neil Massey
#** Date    : 01/08/13
#** Purpose : functions to convert between the 3D cartesian representation of
#**           the icosahedron and the lat / lon of the model
#******************************************************************************

import math
from vector_3D import *

#******************************************************************************

def model_to_cart(lon, lat, R=1.0):
	# convert to radians
	lon_r = math.radians(lon)
	lat_r = math.radians((90.0 - lat))

	# convert spherical to cartesian on a unit sphere
	S = vector_3D(R * math.cos(lon_r) * math.sin(lat_r),
				  R * math.sin(lon_r) * math.sin(lat_r),
				  R * math.cos(lat_r))
	return S

#******************************************************************************

def cart_to_model(V):
	# convert a cartesian coordinate to a "model coordinate" - spherical
	# but with pole at 90 N, lons range from -180 to 180
	rad_to_deg = 180.0 / math.pi;
	phi = math.atan2(V[1], V[0]);
	the = math.acos(V[2] / V.mag());

	lon = phi * rad_to_deg;
	lat = 90 - the * rad_to_deg;
	
	return [lon, lat]