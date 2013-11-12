#! /usr/bin/env python

#******************************************************************************
#** Program : time_mean_regrid_data.py
#** Author  : Neil Massey
#** Date    : 12/08/13
#** Purpose : perform time-averaging on data that has been regridded to the
#**           triangular mesh
#******************************************************************************

import os, sys, getopt
from data_store import *

#******************************************************************************

def average_tri_data(input_file, output_file, avg_period):

	# load the data in
	ds = data_store()
	fh_in = open(input_file, 'rb')
	ds.load(fh_in)
	fh_in.close()
	n_ts = ds.get_n_t_steps()
	n_idxs = ds.get_n_idxs()
	
	# create a new datastore
	mean_ds = data_store()
	mean_ds.mv = ds.mv
	# determine whether there is to be an averaging period
	if avg_period > 1:
		if avg_period == 30 and n_ts != 30:
			avg_period = n_ts
		# set the scaling for each averaging period
		scale = 1.0 / avg_period;
		# create the datastore
		mean_ds.set_size(n_ts/avg_period, n_idxs);
		# loop through the data producing the average
		for t in range(0, n_ts):
			dest_pos = t / avg_period;
			for i in range(0, n_idxs):
				# get the current data value, add the value from the original ds store
				# multiplied by the scaler
				c_val = mean_ds[dest_pos, i]
				t_val = scale * ds[t, i] + c_val
				mean_ds[dest_pos, i] = t_val
	else:
		mean_ds = ds;
	
	# write the data out
	fh_out = open(output_file, 'wb')
	mean_ds.save(fh_out)
	fh_out.close()

#******************************************************************************

if __name__ == "__main__":
	input_file = ""
	output_file = ""
	time_period = 1
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:t:", ["input_file=", "output_file="])
	except:
		print "time_mean_regrid_data -i | --input_file -t | --time_period -o | --output_file"
		sys.exit(2)

	for opt, arg in opts:
		if opt in ['-i', '--input_file']:
			input_file = arg
		if opt in ['-o', '--output_file']:
			output_file = arg
		if opt in ['-t', '--time_period']:
			time_period = int(arg)
			
	average_tri_data(input_file, output_file, time_period)