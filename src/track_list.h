/******************************************************************************
** Program : track_list.h
** Author  : Neil Massey
** Date    : 10/08/09
** Purpose : class that contains track details
******************************************************************************/

#ifndef TRACK_LIST_H
#define TRACK_LIST_H

#include "extrema_list.h"  		// for extrema definition
#include <iostream>
#include <vector>

/*****************************************************************************/

class track_point
{
	friend class track;
	friend class track_list;
	public:
		track_point(void);
		track_point(FP_TYPE iframe);
		void set_point(steering_extremum ipt, FP_TYPE cost, 
					   FP_TYPE c0, FP_TYPE c1, FP_TYPE c2, FP_TYPE c3);
		steering_extremum* get_point(void);
		const FP_TYPE get_frame_number(void) const;

	private:
		steering_extremum pt;	// point
		FP_TYPE frame_number;
		FP_TYPE cost;		// cost of adding to the track
		FP_TYPE c0, c1, c2, c3;	// component costs
};

/*****************************************************************************/

class track
{
	friend class track_list;
	public:
		track(void);
		void add_point(track_point tp);
		void clear(void);	
		const int size(void) const;
		track_point* get_last_track_point(void);
		track_point* get_track_point(int tstep);
		track_point* get_track_point_idx(int idx);
		// candidate point settings
		FP_TYPE get_cand_pt_cost(void);
		steering_extremum* get_cand_pt(void);
		void set_cand_pt(steering_extremum svex_pt, FP_TYPE cost);
		// debug code to enable us to recover the candidate cost contribution
		void set_cand_pt(steering_extremum svex_pt, FP_TYPE cost, FP_TYPE cc0, 
						 FP_TYPE cc1, FP_TYPE cc2, FP_TYPE cc3);
		// consolidate the candidate points
		void consolidate_candidate_point(int t);

		// get length / persistence of the track so that they can be trimmed
		FP_TYPE get_length(void);
		int get_persistence(void);
		// deviation is the distance between the first and last point
		FP_TYPE get_deviation(void);

		// reduce the track to remove static points
		void reduce(void);
		// interpolate the tracks using a spline to obtain a smooth path
		void interpolate(extrema_list* ex_list, int n_spl_pts);	// number of points per spline section

	private:
		std::vector<track_point> tr;
		FP_TYPE length;
		steering_extremum cand_pt;
		FP_TYPE cand_cost;
		// debug code
		FP_TYPE cc0, cc1, cc2, cc3;
};

/*****************************************************************************/

typedef std::vector<track> TR_LIST_TYPE;

/*****************************************************************************/

class track_list
{
	public:
		track_list(void);
		void set_size(int n_tracks);
		const int size(void) const;
		track* get(int n_track);
		const track* get(int n_track) const;
		void add(track trk);
		void consolidate_candidate_points(int t);
		void save(std::string output_fname);
		void save_text(std::string output_fname);
		void load(std::string input_fname);
		META_DATA_TYPE* get_meta_data(void);
		void set_meta_data(META_DATA_TYPE* meta_data);

	private:
		TR_LIST_TYPE tr_list;
		META_DATA_TYPE meta_data;
};

#endif
