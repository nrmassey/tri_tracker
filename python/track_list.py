#******************************************************************************
#** Program : track_list.py
#** Author  : Neil Massey
#** Date    : 31/07/13
#** Purpose : class to represent a track list and to load in from the
#**			  track locating program
#******************************************************************************

from bin_file_utils import *
from extrema_list import extremum
from meta_data import *

class track_point:

	#**************************************************************************

	def __init__(self):
		self.ex = extremum()
		self.frame_number = 0
		self.rules_bf = 0
		
	#**************************************************************************

	def load(self, fh):
		self.frame_number = read_int(fh)
		self.ex.load(fh)
		self.rules_bf = read_int(fh)

#******************************************************************************
		
class track:

	#**************************************************************************

	def __init__(self):
		self.track_points = []

	#**************************************************************************

	def load(self, fh):
		# read number of track points
		trk_pts = read_int(fh)
		assert(trk_pts < 1e7)
		for p in range(0, trk_pts):
			trk_pt = track_point()
			trk_pt.load(fh)
			self.track_points.append(trk_pt)
			
	#**************************************************************************

	def get_n_points(self):
		return len(self.track_points)

	#**************************************************************************
	
	def get_track_point(self, p):
		return self.track_points[p]

#******************************************************************************

class track_list:

	#**************************************************************************

	def __init__(self):
		self.track_list = []
		self.meta_data = {}
		
	#**************************************************************************

	def load(self, fh):
		# read the metadata
		self.meta_data = read_meta_data(fh)
		# load number of tracks
		n_trks = read_int(fh)
		assert(n_trks < 1e7)
		for t in range(0, n_trks):
			trk = track()
			trk.load(fh)
			self.track_list.append(trk)
			
	#**************************************************************************
	
	def get_n_tracks(self):
		return len(self.track_list)

	#**************************************************************************

	def get_meta_data(self):
		return self.meta_data
		
	#**************************************************************************
		
	def get_track(self, t):
		return self.track_list[t]