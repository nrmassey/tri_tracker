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
#include <queue>

/*****************************************************************************/

class tracker
{
	public:
		tracker(std::vector<std::string> input_fname, FP_TYPE sr,
				FP_TYPE hrs_per_t_step);
		~tracker(void);
		void find_tracks(void);
		void save(std::string output_fname);
		void save_text(std::string output_fname);

	protected:
		// input variables
		FP_TYPE sr;				// search radius between frames
		FP_TYPE hrs_per_t_step;	// number of hours per timestep
		std::vector<std::string> input_fname;
		FP_TYPE mv;				// missing value

        // output track list
        track_list tr_list;
        // input extrema list
		extrema_list ex_list;

		// queues to assign extremas and unassigned extremas
		std::queue<steering_extremum*> ex_queue;

		// initializes the tracks in the first frame
		void build_first_frame(void);
		// build the queue of extrema for the timestep
		void build_extrema_queue(int timestep);
		// merge tracks using phantom feature points
		void apply_merge_tracks(void);
		
		int  determine_track_for_candidate(steering_extremum* svex, int t);
		bool assign_candidate(steering_extremum* c_svex, int min_tr, int t);
		void add_unassigned_points_as_tracks(int t);

		// the all important search function - returns a bitfield based on whether
		// the rules have been met (in a specific order)
		int apply_rules(track* TR, steering_extremum* EX_svex);
		// check whether a current candidate event should be overwritten with a new one
		bool should_replace_candidate_point(track* TR, track_point* cur_cand,
											 track_point* new_cand);

};

#endif
