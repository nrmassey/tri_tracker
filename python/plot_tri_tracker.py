#!/usr/bin/env python

###############################################################################
# Program : plot_tri_tracker.py
# Author  : Neil Massey
# Date	  : 08/07/13
# Purpose : Program to plot the grid from the regional version of tri_tracker
###############################################################################

import sys, os, getopt
from netcdf_file import *

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from tri_grid import *
from data_store import *
from extrema_list import *
from track_list import *
from draw_continents import *
from geo_convert import *

###############################################################################

def draw_colorbar(min_V, max_V, n_ticks, colmap, format="%i", title="none"):
	# block plot has been drawn so add the colorbar
	fig = plt.gcf()
	cbh = 0.025
	cbb = 0.05
#	self.sp.set_position([self.cbb, self.cbh+2*self.cbb, 
#						  1.0-2*self.cbb, 1.0-self.cbh-3*self.cbb])
	cb = fig.add_subplot("212")
	cb.set_position([cbb, cbb*1.5, 1.0-2*cbb, cbh])
	# set no y ticks or labels and don't draw the x tick lines
	# finally! draw the colorbar
	sc = abs(float(min_V - max_V)) / colmap.N
	w = int(colmap.N/n_ticks)
	xt = []
	for c in range(0, colmap.N, w):
		x = min_V + c * sc
		col = colmap.__call__(c)
		cb.fill([x,x,x+w*sc,x+w*sc],[0,1,1,0],fc=col,ec=col)
		xt.append(x + 0.5*w*sc)
	cb.set_xlim(min_V, max_V)

	# set the x tick values and labels
	cb.yaxis.set_ticks_position('none')
	cb.yaxis.set_ticks([])
	cb.xaxis.set_ticks(xt)
	xtl = [format % float(xt[i]) for i in range(0, len(xt))]   # limit to two d.p.
	cb.xaxis.set_ticklabels(xtl)
	
	if title != "none":
		cb.set_title(title)

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
			P = [cart_to_model(TRI[0]), cart_to_model(TRI[1]), cart_to_model(TRI[2])]
			# don't draw if it goes around the date line
			if abs(P[0][0] - P[1][0]) > 180 or abs(P[1][0] - P[2][0]) > 180 or abs(P[2][0] - P[1][0]) > 180:
				continue
			fc = [1,1,1,0]
			ec = [0,0,0,0.5]
			sp.fill([P[0][0], P[1][0], P[2][0], P[0][0]], \
					[P[0][1], P[1][1], P[2][1], P[0][1]], \
					facecolor=fc, edgecolor=ec, lw=0.25)

###############################################################################

def plot_data(sp, tg, ds, time_step, level, colmap=cm.RdYlBu_r, dzt=1, keep_scale=False, symmetric_cb=True):
	print "Plotting data"
	# plot the mesh and the values on the mesh
	# get the triangles first
	tri_list = tg.get_triangles_at_level(level)
	ds_list = tg.get_ds_indices_at_level(level)
	# see if we should keep the same scale between timesteps
	if not keep_scale:
		max_V = ds.max(time_step, ds_list)
		min_V = ds.min(time_step, ds_list)
	else:
		max_V = ds.max(-1, ds_list)
		min_V = ds.min(-1, ds_list)
	
	max_V = float(int(max_V / 100))
	min_V = float(int(min_V / 100))
	if (symmetric_cb):
		if (abs(min_V) > abs(max_V)):
			max_V = abs(min_V)
		if (abs(max_V) > abs(min_V)):
			min_V = -max_V
	# calculate the denominator in the normalising equation
	trans = max_V - min_V	
	for t in tri_list:
		# get each corner of the triangle and convert to model coordinates (lat/lon)
		TRI = t.get_data()
		P = [cart_to_model(TRI[0]), cart_to_model(TRI[1]), cart_to_model(TRI[2])]
		# don't draw if it goes around the date line
		if abs(P[0][0] - P[1][0]) > 180 or abs(P[1][0] - P[2][0]) > 180 or abs(P[2][0] - P[1][0]) > 180:
			continue
		# get the value from the datastore
		ds_idx = TRI.get_ds_index()
		V = ds[time_step, TRI.get_ds_index()]+0.5
		if numpy.isnan(V):
			continue
		else:
			V = int(V)
		
		if abs(V) >= abs(ds.get_missing_value())*0.99:
			col = "#FFFFFF"
		else:
			V = float(int(V/100))
			# transform between 0 and 1
			vt = (V - min_V) / trans
			# get the colour from the colour map
			if not (numpy.isnan(vt)):
				col = colmap.__call__(int(vt*colmap.N/dzt)*dzt)
			else:
				col = colmap.__call__(0)
			
		sp.fill([P[0][0], P[1][0], P[2][0], P[0][0]], \
				[P[0][1], P[1][1], P[2][1], P[0][1]], \
				facecolor=col, edgecolor=col, lw=0.25, zorder=0)
				
	# draw a colorbar
	draw_colorbar((min_V), (max_V), 12.25, colmap, format="%i", title="MSLP hPa")

###############################################################################

def plot_extrema(sp, tg, ex, time_step, ex_num=-1):
	print "Plotting extrema"
	# plot the extrema triangles
	# get the extrema for this timestep
	ex_t_step = ex.get_extrema_for_t_step(time_step)
	# plot each extrema
	cur = 0;
	for ex_n in range(0, len(ex_t_step)):
		if ex_num != -1 and ex_n != ex_num:
			continue
		ex_t = ex_t_step[ex_n]
		lon = ex_t.lon
		if lon > 180.0: 	# wrap around date line
			lon = lon-360.0
		sp.plot(lon, ex_t.lat, 'b+', ms=5.0)
		# plot the triangles in the object from their labels
		for tl in ex_t.object_list:
			# get the triangle from the grid via its label
			TRI = tg.get_triangle(tl)			
			# get the points
			P = [cart_to_model(TRI[0]), cart_to_model(TRI[1]), cart_to_model(TRI[2])]
			# don't draw if it goes around the date line
			if abs(P[0][0] - P[1][0]) > 180 or abs(P[1][0] - P[2][0]) > 180 or abs(P[2][0] - P[1][0]) > 180:
				continue
			# draw see-through triangle with black border
			fc = [1,1,1,0]
#			ec = [[0,0,0,1], [1,0,0,1]][cur%2]
			ec = [0,0,0,1]
			sp.fill([P[0][0], P[1][0], P[2][0], P[0][0]], \
					[P[0][1], P[1][1], P[2][1], P[0][1]], \
					facecolor=fc, edgecolor=ec, lw=0.25, zorder=2)
		cur+=1
		
###############################################################################

def plot_objects(sp, tg, ds, ex, time_step, level, ex_num=-1):
	print "Plotting objects"
	# plot the extrema triangles
	# get the extrema for this timestep
	ex_t_step = ex.get_extrema_for_t_step(time_step)
	ds_list = tg.get_ds_indices_at_level(level)
	# see if we should keep the same scale between timesteps
	max_V = ds.max(time_step, ds_list)
	min_V = ds.min(time_step, ds_list)
	
	max_V = float(int(max_V / 100))
	min_V = float(int(min_V / 100))
	trans = max_V - min_V
	dzt=1
	# plot each extrema	
	for ex_n in range(0, len(ex_t_step)):
		if ex_num != -1 and ex_n != ex_num:
			continue
		ex_t = ex_t_step[ex_n]
		lon = ex_t.lon
		if lon > 180.0: 	# wrap around date line
			lon = lon-360.0
		sp.plot(lon, ex_t.lat, 'r+', ms=5.0)
		# plot the triangles in the object from their labels
		for tl in ex_t.object_list:
			# get the triangle from the grid via its label
			TRI = tg.get_triangle(tl)
			# get the points
			P = [cart_to_model(TRI[0]), cart_to_model(TRI[1]), cart_to_model(TRI[2])]
			# don't draw if it goes around the date line
			if abs(P[0][0] - P[1][0]) > 180 or abs(P[1][0] - P[2][0]) > 180 or abs(P[2][0] - P[1][0]) > 180:
				continue
			# fill with the color
			# get the value from the datastore
			ds_idx = TRI.get_ds_index()
			V = int(ds[time_step, TRI.get_ds_index()]+0.5)
		
			if abs(V) >= abs(ds.get_missing_value())*0.99:
				col = "#FFFFFF"
			else:
				V = float(int(V/100))
				# transform between 0 and 1
				vt = (V - min_V) / trans
				# get the colour from the colour map
				if not (numpy.isnan(vt)):
					col = colmap.__call__(int(vt*colmap.N/dzt)*dzt)
				else:
					col = colmap.__call__(0)
			ec=[0,0,0,1]
			sp.fill([P[0][0], P[1][0], P[2][0], P[0][0]], \
					[P[0][1], P[1][1], P[2][1], P[0][1]], \
					facecolor=col, edgecolor=ec, lw=0.25, zorder=2)

###############################################################################

def load_mesh_file(mesh_file):
	tg = tri_grid()
	if not os.path.exists(mesh_file):
		print "Mesh file: " + mesh_file + " not found"
		sys.exit(2)
	else:
		fh = open(mesh_file, 'rb')
		tg.load(fh)
		fh.close()
		print "Loaded mesh file " + mesh_file
	return tg
					
###############################################################################

def load_regrid_file(regrid_file):
	ds = data_store()
	if os.path.exists(regrid_file):
		fh = open(regrid_file, 'rb')
		ds.load(fh)
		fh.close()
		print "Loaded regridded data file " + regrid_file
	else:
		print "Regridded file: " + regrid_file + " not found"
		sys.exit(2)
	return ds

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

def load_track_list(trk_file):
	trk = trk_list()
	if os.path.exists(trk_file):
		fh = open(trk_file, 'rb')
		trk.load(fh)
		fh.close()
		print "Loaded track file " + trk_file
	else:
		print "Track file: " + trk_file + " not found"
		sys.exit(2)
	return trk

###############################################################################

def load_original_grid(orig_mesh_file, orig_mesh_var):
	nc_fh = netcdf_file(orig_mesh_file)
	nc_var = nc_fh.variables[orig_mesh_var]
	# get lat / lon dimension - assume 4D variable
	# should do this via looking at the axis info.
	lat_dim = nc_var.dimensions[2]
	lon_dim = nc_var.dimensions[3]
	lon_vals = nc_fh.variables[lon_dim][:]
	lat_vals = nc_fh.variables[lat_dim][:]
	return lat_vals, lon_vals

###############################################################################

def plot_original_grid(sp, lat_vals, lon_vals):
	# plot the original grid lon / lats as dots
	for la in lat_vals:
		for lo in lon_vals:
			# wrap lo val
			if lo > 180.0: 	# wrap around date line
				lo = lo-360.0
			sp.plot(lo, la, 'ko', ms=1)

###############################################################################

def plot_wind_speed(sp, lat, lon, wnd_speed):
	wind_skp = 1
	print "Plotting winds"
	sc = 0.33
	# plot the wind for the lat / lon coordinates
	for th in range(0, lat.shape[0], wind_skp):
		for lm in range(0, lat.shape[1], wind_skp):
			la = lat[th, lm]
			lo = lon[th, lm]
			if lo > 180.0:
				lo -= 360.0
			sp.plot(lo, la, 'k.', markersize = wnd_speed[0, th, lm] * sc, zorder=1, alpha=0.5)
					  
###############################################################################

def load_wind(wind_file, t_step):
	fh = netcdf_file(wind_file)
	lon = fh.variables["global_longitude1"][:,:]
	lat = fh.variables["global_latitude1"][:,:]
	wnd = fh.variables["field50"][t_step]
	print "Max " + str(numpy.max(wnd))
	
	return wnd, lon, lat
	
###############################################################################

if __name__ == "__main__":
	mesh_file   = ""
	regrid_file = ""
	ex_file     = ""
	ex_num		= -1
	trk_file    = ""
	time_step   = 0
	grid_level  = 0
	draw_mesh   = False
	lat_limits  = [-90, 90]
	lon_limits  = [-180, 180]
	orig_mesh_file = ""
	orig_mesh_var = ""	
	wnd_file    = ""
	
	opts, args = getopt.getopt(sys.argv[1:], 'm:o:r:t:e:l:f:x:y:n:g:w:d', 
							   ['mesh_file=', 'out_name=', 'regrid_file=', 't_step=', 
							   	'extrema_file=', 'level=', 'track_file=', 'draw_mesh',
							   	'lat_limits=', 'lon_limits=', 'extrema_num=',
							   	'orig_file=', 'orig_var=',						   	
							   	'wind_file='])

	for opt, val in opts:
		if opt in ['--mesh_file',   '-m']:
			mesh_file = val
		if opt in ['--regrid_file', '-r']:
			regrid_file = val
		if opt in ['--t_step',      '-t']:
			time_step = int(val)
		if opt in ['--out_name',    '-o']:
			out_name = val
		if opt in ['--extrema',     '-e']:
			ex_file = val
		if opt in ['--track_file',  '-f']:
			trk_file = val
		if opt in ['--extrema_num', '-n']:
			ex_num = int(val)
		if opt in ['--level',       '-l']:
			grid_level = int(val)
		if opt in ['--draw_mesh',   '-d']:
			draw_mesh = True
		if opt in ['--lat_limits',  '-y']:
			y = val.strip("[]").split(",")
			lat_limits = [float(y[0]), float(y[1])]
		if opt in ['--lon_limits',  '-x']:
			x = val.strip("[]").split(",")
			lon_limits = [float(x[0]), float(x[1])]
		if opt in ['--orig_file',  '-g']:
			orig_mesh_file = val
		if opt in ['--orig_var',   '-v']:
			orig_mesh_var = val			
		if opt in ['--wind_file',  '-w']:
			wind_file = val

	# load each part in
	tg = load_mesh_file(mesh_file)
	if regrid_file != "":
		ds = load_regrid_file(regrid_file)
	if ex_file != "":		
		ex = load_extrema_file(ex_file)
	if trk_file != "":
		trk = load_track_file(trk_file)
		
	# decide what to do with the data
	sp = plt.subplot("111")
	colmap=cm.RdYlBu_r
	if regrid_file != "":
		plot_data(sp, tg, ds, time_step, grid_level, colmap=colmap, dzt=1, keep_scale=False, symmetric_cb=False)

	if mesh_file != "" and draw_mesh:
		plot_mesh(sp, tg, grid_level)

	if wind_file != "":
		wnd, wnd_lon, wnd_lat = load_wind(wind_file, time_step)
		plot_wind_speed(sp, wnd_lat, wnd_lon, wnd)

	if ex_file != "":
		plot_extrema(sp, tg, ex, time_step, ex_num)
	
	if orig_mesh_file != "" and orig_mesh_var != "":
		# load and plot the original mesh as dots
		lon, lat = load_original_grid(orig_mesh_file, orig_mesh_var)
		plot_original_grid(sp, lon, lat)
	
	draw_continents(sp)
	sp.set_aspect(1.0)
	print "Done"
	sp.set_xlim(lon_limits)
	sp.set_ylim(lat_limits)
	plt.savefig(out_name, bbox_inches='tight', dpi=500)
	print "Saved to file: " + out_name
