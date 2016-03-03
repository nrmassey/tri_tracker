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
                FP_TYPE ov, FP_TYPE hrs_per_t_step, int opt_steps);
        ~tracker(void);
        void find_tracks(void);
        void save(std::string output_fname);
        void save_text(std::string output_fname);

    protected:
        // input variables
        FP_TYPE sr;             // search radius between frames
        FP_TYPE ov;             // overlap percentage now tunable
        FP_TYPE hrs_per_t_step; // number of hours per timestep
        std::vector<std::string> input_fname;
        FP_TYPE mv;             // missing value
        int opt_steps;          // number of optimisation iterations

        // output track list
        track_list tr_list;
        // input extrema list
        extrema_list ex_list;

        // queues to assign extremas and unassigned extremas
        std::queue<steering_extremum> ex_queue;

        // initializes the tracks in the first frame
        void build_first_frame(void);
        // build the queue of extrema for the timestep
        void build_extrema_queue(int timestep);
        
        // optimisation processes
        void add_phantom_points(track* trk_P);
        void remove_phantom_points(track* trk_P);
        void interpolate_phantom_points(track* trk_P);
        void apply_optimise_tracks(void);
        // get a list of overlapping tracks (overlapping in time) for a particular track number
        std::vector<int> get_overlapping_tracks(track* tr_A);
        void apply_merge_tracks(void);

        // functions for deriving the initial tracks
        int  determine_track_for_candidate(track_point& n_cand, int t);
        bool assign_candidate(track_point& n_cand, int min_tr, int t);
        void add_unassigned_points_as_tracks(int t);

        // the all important search function - returns a bitfield based on whether
        // the rules have been met (in a specific order) and the cost of meeting that
        // rule.  All values returned within the candidate point
        void apply_rules(track* TR, track_point& n_cand);
};

#endif
