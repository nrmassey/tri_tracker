#!/usr/bin/env python

###############################################################################
# Program : plot_all_extrema.py
# Author  : Neil Massey
# Date	  : 08/07/13
# Purpose : Program to plot all extrema in an extrema file
###############################################################################

import sys, os, getopt
from netcdf_file import *

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
from matplotlib.patches import Wedge as Wedge

from extrema_list import *
from vector_3D import *
from geo_convert import *

#pole_latitude=39.25
#pole_longitude=198.0

#pole_latitude=38.0
#pole_longitude=190.0

###############################################################################

def calc_bearing(lon1, lat1, lon2, lat2):

    # METHOD 1
 	nlat1 = math.radians(lat1)
 	nlat2 = math.radians(lat2)
 	dlon  = math.radians(lon2 - lon1)
 	
 	a = math.sin(dlon)*math.cos(nlat2)
 	b = math.cos(nlat1)*math.sin(nlat2) - math.sin(nlat1)*math.cos(nlat2)*math.cos(dlon)
 	ang = math.atan2(a,b)	        	
 	ang = math.degrees(ang)
 	
	return ang;

###############################################################################

def plot_extrema(sp, ex):
	print "Plotting extrema"
	# plot the extrema triangles
	# plot each extrema for each timestep
	n_t = ex.get_n_t_steps()
	earth_r = 6371 * 1000
	t_step = (6.0 * 60 * 60) # seconds per timestep (6 hourly)
	lat_c = (2*math.pi * earth_r)
	for t in range(53, 55):
		color = (1.0*float(t)/n_t,0.0,0.0)
		color2 = (0.0,1.0*float(t)/n_t,0.0)
		# get the extrema for this timestep
		ex_t_step = ex.get_extrema_for_t_step(t)
		for ex_n in range(0, len(ex_t_step)):
			ex_t = ex_t_step[ex_n]
			if len(ex_t.object_list) == 0:
				continue
			lon = ex_t.lon
			lat = ex_t.lat
			if (lon >= 180.0):
				lon -= 360.0
			sp.plot(lon, lat,'ro', ms=2.5, mec='r', zorder=3)
			# plot geostrophic wind if necessary
			if ex_t.steer_x != 0.0 and abs(ex_t.steer_x) < 1e5:
			    # add the geostrophic wind (in sph. coords.) to the position
				circ = 2*math.pi * math.cos(math.radians(lat)) * earth_r
				P_lon = lon + ex_t.steer_x * t_step * 360.0 / circ
				P_lat = lat + ex_t.steer_y * t_step * 180.0 / lat_c
			    # plot arrow between the two points
				sp.arrow(lon, lat, P_lon-lon, P_lat-lat, head_width=0.3, 
				         fc=color2, ec=color2, zorder=1)
				# plot circle of radius equal to the magnitude of the geostrophic wind
				# over six hours
				P3_lon = lon + 500 * 1000 * 360.0/circ
#				sr = plt.Circle((lon,lat),P3_lon-lon,color=color, fill=False)
				# calculate bearing of PR_lon,PR_lat from lon, lat
				ang = calc_bearing(lon, lat, P_lon, P_lat)
				V = math.sqrt((P_lon-lon)**2+(P_lat-lat)**2)
				sr = Wedge((lon,lat),V, -ang, 180-ang, color=color, fill=False,ls='solid')
				plt.gca().add_artist(sr)
			sp.text(lon-0.5, lat, str(t), fontsize=6, zorder=0)
		
###############################################################################

def load_extrema_file(ex_file):
	ex = extrema_list()
	if ex_file != "":		
		if os.path.exists(ex_file):
			fh = open(ex_file, 'rb')
			ex.load(fh)
			fh.close()
			print "Loaded extrema file " + ex_file
		else:
			print "Extrema file: " + ex_file + " not found"
			sys.exit(2)
	return ex

###############################################################################

if __name__ == "__main__":
	ex_file     = ""
	lat_limits  = [-90, 90]
	lon_limits  = [-180, 180]
	
	opts, args = getopt.getopt(sys.argv[1:], 'e:y:x:o:', 
	                           ['extrema_file=', 'lat_limits=', 'lon_limits=', 'out_name='])

	for opt, val in opts:
		if opt in ['--extrema_file','-e']:
			ex_file = val
		if opt in ['--lat_limits',  '-y']:
			y = val.strip("[]").split(",")
			lat_limits = [float(y[0]), float(y[1])]
		if opt in ['--lon_limits',  '-x']:
			x = val.strip("[]").split(",")
			lon_limits = [float(x[0]), float(x[1])]
		if opt in ['--out_name',    '-o']:
			out_name = val

	# load each part in
	if ex_file != "":		
		ex = load_extrema_file(ex_file)

	projection = ccrs.PlateCarree()
	sp0 = plt.subplot(111, projection=projection)
	if ex_file != "":
		plot_extrema(sp0, ex)
	
	sp0.coastlines()
	sp0.gridlines()
	sp0.get_axes().set_extent([lon_limits[0], lon_limits[1], lat_limits[0], lat_limits[1]])
	sp0.set_aspect(1.0)
	plt.savefig(out_name, bbox_inches='tight', dpi=250)
	print "Saved to file: " + out_name
	plt.close()
