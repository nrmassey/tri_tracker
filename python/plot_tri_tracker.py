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
import matplotlib.colors as col
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec

from tri_grid import *
from data_store import *
from extrema_list import *
from track_list import *
from geo_convert import *
from conv_rot_grid import *

#pole_latitude=39.25
#pole_longitude=198.0

#pole_latitude=38.0
#pole_longitude=190.0

###############################################################################

def draw_colorbar(sp1, min_V, max_V, n_ticks, colmap, format="%i", title="none", label_every=1):
	# block plot has been drawn so add the colorbar
	# set no y ticks or labels and don't draw the x tick lines
	xt = []
	n = 1.0/(max_V-min_V)
	w = (max_V - min_V) / n_ticks
	for c in range(int(min_V), int(max_V)+int(w), int(w)):
		col = colmap.__call__((c-min_V)*n)
		sp1.fill([c,c,c+w,c+w],[0,1,1,0],fc=col,ec=col)
		xt.append(c)
	sp1.set_xlim(min_V, max_V)

	# set the x tick values and labels
	sp1.yaxis.set_ticks_position('none')
	sp1.yaxis.set_ticks([])
	sp1.xaxis.set_ticks(xt)
	xtl = []
	for i in range(0, len(xt)):
		if (i % label_every == 0):
			xtl.append(format % float(xt[i]))
		else:
			xtl.append("")
	
	sp1.xaxis.set_ticklabels(xtl)
	
	if title != "none":
		sp1.set_title(title)
	A = (max_V - min_V) * 0.03
	sp1.set_aspect(A)

###############################################################################

def plot_mesh(sp, tg, l, tl=-1, pole_longitude=0.0, pole_latitude=90.0):
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
			if pole_latitude != 90.0:
				P = [glob2rot(P[0][0], P[0][1], pole_longitude, pole_latitude),
					 glob2rot(P[1][0], P[1][1], pole_longitude, pole_latitude),
					 glob2rot(P[2][0], P[2][1], pole_longitude, pole_latitude)]				
			fc = [1,1,1,0]
			ec = [0,0,0,0.5]
			sp.fill([P[0][0], P[1][0], P[2][0], P[0][0]], \
					[P[0][1], P[1][1], P[2][1], P[0][1]], \
					facecolor=fc, edgecolor=ec, lw=0.25)

###############################################################################

def plot_data(sp0, sp1, tg, ds, time_step, level, colmap=cm.RdYlBu_r, dzt=1, keep_scale=False,\
			  symmetric_cb=True, pole_longitude=0.0, pole_latitude=90.0):
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
	print "Plotting data range: ", min_V, max_V

	max_V = 1050
	min_V = 950
#	max_V = 50
#	min_V = -50
	
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
		if pole_latitude != 90.0:
			P = [glob2rot(P[0][0], P[0][1], pole_longitude, pole_latitude),
				 glob2rot(P[1][0], P[1][1], pole_longitude, pole_latitude),
				 glob2rot(P[2][0], P[2][1], pole_longitude, pole_latitude)]
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
			
		sp0.fill([P[0][0], P[1][0], P[2][0], P[0][0]], \
			 	 [P[0][1], P[1][1], P[2][1], P[0][1]], \
			 	 facecolor=col, edgecolor=col, lw=0.25, zorder=0)
				
	# draw a colorbar
	draw_colorbar(sp1, (min_V), (max_V), int(max_V-min_V), colmap, format="%i", label_every=12)#title="MSLP hPa")

###############################################################################

def plot_geostrophic_wind(sp, steer_x, steer_y, lon, lat, pole_longitude=0.0, pole_latitude=90.0):
	# plot geostrophic wind if necessary
	n_hrs_per_tstep = 6
	n_secs_per_tstep = n_hrs_per_tstep * 60 * 60
	earth_r = 6371 * 1000
	lat_c = (2*math.pi * earth_r)
	
	if steer_x != 0.0 and abs(steer_x) < 1e5:
		# add the geostrophic wind (in sph. coords.) to the position
		circ = 2*math.pi * math.cos(math.radians(lat)) * earth_r
		P_lon = lon + steer_x * n_secs_per_tstep * 360.0 / circ
		P_lat = lat + steer_y * n_secs_per_tstep * 180.0 / lat_c
		PG = glob2rot(P_lon, P_lat, pole_longitude, pole_latitude)
		E = glob2rot(lon, lat, pole_longitude, pole_latitude)
		# plot arrow between the two points
		sp.arrow(E[0], E[1], PG[0]-E[0], PG[1]-E[1], head_width=0.3, 
		         fc='r', ec='r', zorder=1)

###############################################################################

def plot_object_triangles(sp, tg, object_list, pole_longitude=0.0, pole_latitude=90.0):
	# plot the triangles in the object from their labels
	for tl in object_list:
		# get the triangle from the grid via its label
		TRI = tg.get_triangle(tl)			
		# get the points
		P = [cart_to_model(TRI[0]), cart_to_model(TRI[1]), cart_to_model(TRI[2])]
		# don't draw if it goes around the date line
		if abs(P[0][0] - P[1][0]) > 180 or abs(P[1][0] - P[2][0]) > 180 or abs(P[2][0] - P[1][0]) > 180:
			continue
		if (pole_latitude != 90.0):
			P = [glob2rot(P[0][0], P[0][1], pole_longitude, pole_latitude),
				 glob2rot(P[1][0], P[1][1], pole_longitude, pole_latitude),
				 glob2rot(P[2][0], P[2][1], pole_longitude, pole_latitude)]				
		# draw see-through triangle with black border
		fc = [1,1,1,0]
		ec = [0,0,0,1]
		sp.fill([P[0][0], P[1][0], P[2][0], P[0][0]], \
				[P[0][1], P[1][1], P[2][1], P[0][1]], \
				facecolor=fc, edgecolor=ec, lw=0.5, zorder=2)

###############################################################################

def plot_extrema(sp, tg, ex, time_step, ex_num=-1, pole_longitude=0.0, pole_latitude=90.0):
	print "Plotting extrema"
	# plot the extrema triangles
	# get the extrema for this timestep
	ex_t_step = ex.get_extrema_for_t_step(time_step)
	# plot each extrema
	cur = 0
	
	for ex_n in range(0, len(ex_t_step)):
		if ex_num != -1 and ex_n != ex_num:
			continue
		ex_t = ex_t_step[ex_n]
		lon = ex_t.lon
		lat = ex_t.lat
		if abs(lon) > 1000 or abs(lat) > 1000:
			# latitude or longitude is the missing value so get the point from the first triangle
			TRI = tg.get_triangle(ex_t.object_list[0])
			C = cart_to_model(TRI.centroid)
			lon = C[0]
			lat = C[1]
		E = glob2rot(lon, lat, pole_longitude, pole_latitude)
		sp.plot(E[0], E[1],'ro', ms=2.5, mec='r', zorder=3)
		plot_object_triangles(sp, tg, ex_t.object_list, pole_longitude, pole_latitude)
#		plot_geostrophic_wind(sp, ex_t.steer_x, ex_t.steer_y, lon, lat, pole_longitude, pole_latitude)

		cur+=1
		
###############################################################################

def plot_objects(sp, tg, ds, ex, time_step, level, ex_num=-1, pole_longitude=0.0, pole_latitude=90.0):
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
	dzt=10
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
			if pole_latitude != 90.0:
				P = [glob2rot(P[0][0], P[0][1], pole_longitude, pole_latitude),
					 glob2rot(P[1][0], P[1][1], pole_longitude, pole_latitude),
					 glob2rot(P[2][0], P[2][1], pole_longitude, pole_latitude)]				
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

def load_track_file(trk_file):
	trk = track_list()
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

def plot_wind_speed(sp, lat, lon, wnd_speed, max_wnd):
	wind_skp = 1
	print "Plotting winds"
	sc = 0.175
	# plot the wind for the lat / lon coordinates
	for th in range(0, lat.shape[0], wind_skp):
		for lm in range(0, lon.shape[0], wind_skp):
			la = lat[th]
			lo = lon[lm]
			if lo > 180.0:
				lo -= 360.0
			a = wnd_speed[0,th,lm] / max_wnd
			sp.plot(lo, la, 'k.', markersize = wnd_speed[0,th,lm] * sc, zorder=1, alpha=a, mec='k')

###############################################################################

def load_wind(wind_file):
	fh = netcdf_file(wind_file)
	lon = fh.variables["longitude1"][:]
	lat = fh.variables["latitude1"][:]
	wnd = fh.variables["field50"][:]
	print "Max " + str(numpy.max(wnd))
	
	return wnd, lon, lat, numpy.max(wnd)
	
###############################################################################

def create_wind_color_map(levels):
    # replicate XWS colours for wind speed
    #
    cmap = ["#ffffff", "#fff4e8", "#ffe1c1", "#ffaa4e", "#ff6d00",
    		"#d33100", "#890000"] #, "#650000", "#390000"]
    ccmap, norm = col.from_levels_and_colors(levels, cmap, 'neither')
    return ccmap, norm

###############################################################################

def plot_wind_field(sp, wnd_lat, wnd_lon, wnd, wnd_max):
	cbar_vals = [x for x in range(0,40,5)]
	cmap, norm = create_wind_color_map(cbar_vals)
	pmesh = sp.pcolormesh(wnd_lon, wnd_lat, wnd, cmap=cmap, vmax=cbar_vals[-1], vmin=cbar_vals[0], norm=norm)
	draw_colorbar(cbar_vals[0], cbar_vals[-1]-5, len(cbar_vals)-2, cmap, format="%i", title="Max wind speed m/s")

###############################################################################

def plot_track(sp, trk, t_step, pole_longitude, pole_latitude):
	print "Plotting tracks"
	# plot the tracks which encompass the current time step
	six_hrs = 6 * 60 * 60
	for ct in range(0, trk.get_n_tracks()):
		c_trk = trk.get_track(ct)
		# check whether the current timestep is in the range of the track
		fn = c_trk.get_track_point(0).frame_number
		np = c_trk.get_n_points()
		if t_step >= fn and t_step < fn + np:
			# plot the track - first point first
			c_tp = c_trk.get_track_point(0)
			P_0 = glob2rot(c_tp.ex.lon, c_tp.ex.lat, pole_longitude, pole_latitude)
			sp.plot(P_0[0], P_0[1], 'ko', ms=2.5)
			earth_r = 6371 * 1000
			circ = 2*math.pi * math.cos(math.radians(c_tp.ex.lat)) * earth_r

			if c_tp.frame_number == t_step:
				V = math.sqrt(c_tp.ex.steer_x**2 + c_tp.ex.steer_y**2)
				if V * six_hrs / 1000.0 < 500:
					R = 500 * 1000
				else:
					R = V * 1.5 * six_hrs
				P3_lon = c_tp.ex.lon + V * 1.5 * six_hrs * 360.0/circ
				P3 = glob2rot(P3_lon, c_tp.ex.lat, pole_longitude, pole_latitude)
				sr = plt.Circle((P_0[0],P_0[1]),P3[0]-P_0[0],color='r', fill=False)
				sp.add_artist(sr)
				plot_geostrophic_wind(sp, c_tp.ex.steer_x, c_tp.ex.steer_y, c_tp.ex.lon, c_tp.ex.lat, pole_longitude, pole_latitude)
				plot_object_triangles(sp, tg, c_tp.ex.object_list, pole_longitude, pole_latitude)
				
			for tp in range(1, np):
				# convert the 2nd point
				c_tp = c_trk.get_track_point(tp)
				P_1 = glob2rot(c_tp.ex.lon, c_tp.ex.lat, pole_longitude, pole_latitude)
				sp.plot(P_1[0], P_1[1], 'ko', ms=2.5)
				
				if c_tp.frame_number == t_step:
					V = math.sqrt(c_tp.ex.steer_x**2 + c_tp.ex.steer_y**2)
					if V * six_hrs / 1000.0 < 500:
						R = 500 * 1000
					else:
						R = V * 1.5 * six_hrs
					P3_lon = c_tp.ex.lon + R * 360.0/circ
					P3 = glob2rot(P3_lon, c_tp.ex.lat, pole_longitude, pole_latitude)
					sr = plt.Circle((P_1[0],P_1[1]),P3[0]-P_1[0],color='r', fill=False)
					sp.add_artist(sr)					
					plot_geostrophic_wind(sp, c_tp.ex.steer_x, c_tp.ex.steer_y, c_tp.ex.lon, c_tp.ex.lat, pole_longitude, pole_latitude)
					plot_object_triangles(sp, tg, c_tp.ex.object_list, pole_longitude, pole_latitude)

				sp.plot([P_0[0], P_1[0]], [P_0[1], P_1[1]], 'k')
				
				P_0 = P_1
				
###############################################################################

if __name__ == "__main__":
	mesh_file   = ""
	regrid_file = ""
	wind_file   = ""
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
	n_steps     = 1
	pole_latitude = 90.0
	pole_longitude = 0.0
	
	opts, args = getopt.getopt(sys.argv[1:], 'm:o:r:t:s:e:l:f:x:y:n:g:w:p:q:d', 
							   ['mesh_file=', 'out_name=', 'regrid_file=', 't_step=', 'n_steps',
							   	'extrema_file=', 'level=', 'track_file=', 'draw_mesh',
							   	'lat_limits=', 'lon_limits=', 'extrema_num=',
							   	'orig_file=', 'orig_var=',						   	
							   	'wind_file=',
							   	'pole_lon=', 'pole_lat='])

	for opt, val in opts:
		if opt in ['--mesh_file',   '-m']:
			mesh_file = val
		if opt in ['--regrid_file', '-r']:
			regrid_file = val
		if opt in ['--t_step',      '-t']:
			time_step = int(val)
		if opt in ['--n_steps',     '-s']:
			n_steps = int(val)
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
		if opt in ['--pole_lon',  '-p']:
			pole_longitude = float(val)
		if opt in ['--pole_lat',  '-q']:
			pole_latitude = float(val)

	# load each part in
	tg = load_mesh_file(mesh_file)
	if regrid_file != "":
		ds = load_regrid_file(regrid_file)
	if ex_file != "":		
		ex = load_extrema_file(ex_file)
	if trk_file != "":
		trk = load_track_file(trk_file)
	if wind_file != "":
		wnd, wnd_lon, wnd_lat, wnd_max = load_wind(wind_file)
		
	if orig_mesh_file != "" and orig_mesh_var != "":
		# load and plot the original mesh as dots
		o_lon, o_lat = load_original_grid(orig_mesh_file, orig_mesh_var)
		
	# decide what to do with the data
	gs = gridspec.GridSpec(7,8)
	
	if pole_longitude != 0.0 and pole_latitude != 90.0:
		projection = ccrs.RotatedPole(pole_latitude=pole_latitude, pole_longitude=pole_longitude)
	else:
		projection = ccrs.PlateCarree()
	colmap=cm.RdYlBu_r
	for t_step in range(time_step, time_step+n_steps):
		sp0 = plt.subplot(gs[0:6,:], projection=projection)
		if regrid_file != "":
			sp1 = plt.subplot(gs[6,1:-1])
			plot_data(sp0, sp1, tg, ds, t_step, grid_level, colmap=colmap, dzt=10, keep_scale=True,\
					  symmetric_cb=False, pole_longitude=pole_longitude, pole_latitude=pole_latitude)

		if mesh_file != "" and draw_mesh:
			plot_mesh(sp0, tg, grid_level, pole_longitude=pole_longitude, pole_latitude=pole_latitude)

		if wind_file != "":
			if regrid_file != "":
				plot_wind_speed(sp0, wnd_lat, wnd_lon, wnd[t_step], wnd_max)
			else:
				plot_wind_field(sp0, wnd_lat, wnd_lon, wnd[t_step,0], wnd_max)

		if ex_file != "":
			plot_extrema(sp0, tg, ex, t_step, ex_num, pole_longitude=pole_longitude, pole_latitude=pole_latitude)
	
		if trk_file != "":
			plot_track(sp0, trk, t_step, pole_longitude=pole_longitude, pole_latitude=pole_latitude)
	    
		if orig_mesh_file != "" and orig_mesh_var != "":
			plot_original_grid(sp0, o_lon, o_lat)

		sp0.coastlines()
		sp0.gridlines()
		sp0.get_axes().set_extent([lon_limits[0], lon_limits[1], lat_limits[0], lat_limits[1]])
		sp0.set_aspect(1.0)
		if n_steps > 0:
			this_out_name = out_name[:-4] + "_%03i" % t_step + ".png"
		else:
			this_out_name = out_name
		plt.savefig(this_out_name, bbox_inches='tight', dpi=125)
		print "Saved to file: " + this_out_name
		plt.close()
