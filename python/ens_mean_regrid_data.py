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

def ens_avg_tri_data(input_list, output_file):
	# read the input list
	fh_list = open(input_list)
	ens_files = fh_list.readlines()
	fh_list.close()

	scaler = 1.0/len(ens_files)

	# create a new datastore
	mean_ds = data_store()
	
	for e in ens_files:
		# load the data in
		ds = data_store()
		fh_in = open(e.strip('\n'), 'rb')
		ds.load(fh_in)
		fh_in.close()
		n_ts = ds.get_n_t_steps()
		n_idxs = ds.get_n_idxs()
	
		# resize the datastore if this is the first load
		if (mean_ds.get_n_t_steps() == 0):
			mean_ds.set_size(n_ts, n_idxs)
			mean_ds.mv = ds.mv
		# loop through the data producing the average
		for t in range(0, n_ts):
			for i in range(0, n_idxs):
				# get the current data value, add the value from the original ds store
				# multiplied by the scaler
				c_val = mean_ds[t, i]
				t_val = scaler * ds[t, i] + c_val
				mean_ds[t, i] = t_val
	
	# write the data out
	fh_out = open(output_file, 'wb')
	mean_ds.save(fh_out)
	fh_out.close()

#******************************************************************************

if __name__ == "__main__":
	input_list = ""
	output_file = ""
	time_period = 1
	try:
		opts, args = getopt.getopt(sys.argv[1:], "l:o:", ["input_list=", "output_file="])
	except:
		print "time_mean_regrid_data -l | --input_list -o | --output_file"
		sys.exit(2)

	for opt, arg in opts:
		if opt in ['-l', '--input_list']:
			input_list = arg
		if opt in ['-o', '--output_file']:
			output_file = arg
			
	ens_avg_tri_data(input_list, output_file)