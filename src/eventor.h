/******************************************************************************
** Program : eventor.h
** Author  : Neil Massey
** Date    : 11/04/14
** Purpose : class to read in regridded data on a triangular mesh and 
**           construct an event set based on that data
******************************************************************************/

#ifndef EVENTOR_H
#define EVENTOR_H

#include <string>
#include <vector>
#include <list>
#include <queue>

#include "extrema_list.h"
#include "event_track_list.h"
#include "ncdata.h"
#include "tri_grid.h"
#include <fstream>

/*****************************************************************************/

class eventor
{
	public:
		eventor(std::vector<std::string> iinput_fname, int imin_per, 
				FP_TYPE imin_len, FP_TYPE imin_dev, FP_TYPE isr,
				int ievent_t_steps, FP_TYPE ioverlap,
				std::vector<std::string> mslp_file_name, std::string mslp_field_name,
				std::vector<std::string> wind_file_name, std::string wind_field_name,
				std::string lsm_file_name, std::string mesh_file);
		~eventor(void);
		void find_events(void);
		void save(std::string output_fname);
		
	private:
		int min_per;
		FP_TYPE min_len;
		FP_TYPE min_dev;
		FP_TYPE sr;
		int event_t_steps;
		FP_TYPE min_overlap;
		int hours_per_t_step;
		
		// extrema related files
		extrema_list ex_list;
		FP_TYPE mv;				// missing value
		
		// functions to build the event set
		void build_first_frame(void);
		void build_event_set(void);
		bool evaluate_candidate_event(steering_extremum* trk_pt, steering_extremum* cand_pt,
									  FP_TYPE& dist, FP_TYPE& overlap);
		void clear_queue(void);
		
		// post processing of tracks
		void trim_tracks(void);
		bool merge_events(void);
		
		// queue to assign extremas
		std::queue<steering_extremum*> ex_queue;
		
		// store data
		event_track_list event_list;
		
		// data for original mslp and wind fields and LSM
		std::vector<ncdata*> mslp_field;
		std::vector<ncdata*> wind_speed_field;
		std::vector<FP_TYPE*> time_field;
		std::vector<int> time_length;
		ncdata* lsm;
		// reference time
		int ref_year, ref_month, ref_day;
		FP_TYPE ref_day_sc, ref_ndays_py;
		
		// functions to get the mslp data / wind speed data / lsm
		FP_TYPE get_mslp_data(int t, int y, int x);
		FP_TYPE get_wind_speed_data(int t, int y, int x);
		FP_TYPE get_lsm(int y, int x);
		FP_TYPE get_time(int t);
		
		// calculate the maximum wind speed of an event
		void calc_max_ws_min_mslp(int event_track_number, int& event_point_index, 
								  FP_TYPE& max_ws, FP_TYPE& min_mslp, int& t_step);
		// calculate the footprint
		FP_TYPE* wind_footprint;
		FP_TYPE* mslp_footprint;
		int footprint_width;
		int footprint_height;
		void calculate_wind_footprint(int event_track_number, int event_point_index,
									  int& n_timesteps, int& n_points, int& central_t_step);
		void save_72hour_footprint(std::string out_name, int event_number);
		void save_time_varying_footprint(std::string out_name, int event_number);
		void calculate_time_step_footprint(int event_track_number, int event_point_index,
										   int& n_points);
		void write_footprint(std::ofstream& fh);
		// triangular mesh
		tri_grid tg;
};

#endif