/******************************************************************************
** Program : tracker.h
** Author  : Neil Massey
** Date    : 05/08/09
** Purpose : class that connects extrema feature points into a track
******************************************************************************/

#ifndef TRACKER_H
#define TRACKER_H

#include "extrema_list.h"
#include "track_list.h"
#include "ncdata.h"
#include <vector>
#include <string>
#include <stack>

/*****************************************************************************/

class tracker
{
	public:
		tracker(std::vector<std::string> input_fname, int min_per, 
				FP_TYPE min_len, FP_TYPE min_dev, FP_TYPE sr);
		~tracker(void);
		void find_tracks(void);
		void save(std::string output_fname);
		void save_text(std::string output_fname);
		void set_weights(FP_TYPE w0, FP_TYPE w1, FP_TYPE w2, FP_TYPE w3);
		void set_tsteps(int t);
		void set_n_spl_pts(int s);	// number of points in interpolated spline
		void set_ksteps(int k);		// number of timesteps to skip between

	private:
		int min_per;			// minimum persistence of track
		FP_TYPE min_len;		// minimum length
		FP_TYPE min_dev;		// minimum deviation
		FP_TYPE sr;				// search radius between frames
		extrema_list ex_list;
		track_list tr_list;
		FP_TYPE w0, w1, w2, w3;	// weights for the algorithm
		int ts;					// number of permitted timesteps between track points
		int n_spl;				// number of spline pts in interpolant
		int ksteps;				// skip between number of timesteps
		std::vector<std::string> input_fname;
		FP_TYPE mv;				// missing value
		FP_TYPE MAX_COST;		// maximum cost - actually a scaler * the sum of the weights

		// stacks to assign extremas and unassigned extremas
		std::stack<steering_extremum> ex_stack;
		std::stack<steering_extremum> ua_ex_stack;
		void clear_stacks(void);

		// initializes the tracks in the first frame
		void build_first_frame(void);
		int  determine_track_for_point(steering_extremum* svex, int t);
		void assign_point_to_track(steering_extremum* c_svex, int min_tr, int t);
		void add_unassigned_points_as_tracks(int t);

		// the all important cost_function - returns the constituent components
		FP_TYPE cmpt_cost_fn(track* TR, steering_extremum* EX_svex, int t,
                             FP_TYPE& c0, FP_TYPE& c1, FP_TYPE& c2, FP_TYPE& c3);

		// reduce the tracks to remove the timesteps where no movement occurs
		void reduce_tracks(void);
		// trim any tracks that don't meet the minimum persitence / distance
		void trim_tracks(void);
		// reform the tracks so that timesteps t(n+1) = t(n)+1
		void reform_tracks(void);
		// interpolate tracks along a spline
		void interpolate_tracks(void);
};

#endif
