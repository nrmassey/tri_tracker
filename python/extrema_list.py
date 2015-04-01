#******************************************************************************
#** Program : extrema_list.py
#** Author  : Neil Massey
#** Date    : 31/07/13
#** Purpose : class to represent an extrema list and to load in from the
#**           extrema locating program
#******************************************************************************

from bin_file_utils import *
from meta_data import *

class extremum:

    #**************************************************************************
    
    def __init__(self, ilon=0.0, ilat=0.0, iintensity=0.0, idelta=0.0, \
                 isteer_x = 0.0, isteer_y = 0.0):
        self.lon = ilon
        self.lat = ilat
        self.intensity = iintensity
        self.delta = idelta
        self.steer_x = isteer_x
        self.steer_y = isteer_y
        self.object_list = []
        
    #**************************************************************************

    def load(self, fh):
        self.lon = read_float(fh)
        self.lat = read_float(fh)
        self.intensity = read_float(fh)
        self.delta = read_float(fh)
        self.steer_x = read_float(fh)
        self.steer_y = read_float(fh)
        # read the object list
        n_objs = read_int(fh)
        assert(n_objs < 1e7)
        for o in range(0, n_objs):
            obj_label = read_label(fh)
            self.object_list.append(obj_label)
        
#******************************************************************************

class extrema_list:

    #**************************************************************************

    def __init__(self):
        self.ex_list = []
        self.mv = 2e20
        self.meta_data = {}
        
    #**************************************************************************

    def get_n_t_steps(self):
        return len(self.ex_list)

    #**************************************************************************

    def get_n_extrema(self, t_step):
        return len(self.ex_list[t_step])

    #**************************************************************************

    def get_total_n_extrema(self):
        sum = 0
        for t in range(0, len(self.ex_list)):
            sum += self.get_n_extrema(t)
        return sum

    #**************************************************************************

    def get_extrema_for_t_step(self, t_step):
        return self.ex_list[t_step]

    #**************************************************************************

    def get(self, t_step, ex_n):
        return self.ex_list[t_step][ex_n]

    #**************************************************************************
    
    def get_missing_value(self):
        return self.mv
    
    #**************************************************************************

    def get_meta_data(self):
        return self.meta_data

    #**************************************************************************

    def load(self, fh):
        # read the meta data in
        self.meta_data = read_meta_data(fh)
        # missing value first then number of time steps
        self.mv = read_float(fh)
        n_t_steps = read_int(fh)
        assert(n_t_steps < 1e7)
        for t in range(0, n_t_steps):
            # add a blank list to load into
            self.ex_list.append([])
            # read number of extrema in this timestep
            n_ex = read_int(fh)
            assert(n_ex < 1e7)
            for e in range(0, n_ex):
                ex = extremum()
                ex.load(fh)
                # add to the last list in the extrema list
                self.ex_list[-1].append(ex)
