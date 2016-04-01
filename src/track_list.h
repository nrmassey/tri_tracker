/******************************************************************************
** Program : track_list.h
** Author  : Neil Massey
** Date    : 10/08/09
** Purpose : class that contains track details
******************************************************************************/

#ifndef TRACK_LIST_H
#define TRACK_LIST_H

#include "extrema_list.h"       // for extrema definition
#include <iostream>
#include <vector>

/*****************************************************************************/

class track_point
{
    public:
        track_point(void) : timestep(0), rules_bf(0), cost(2e20) {}
        track_point(steering_extremum ipt, int itimestep, int irules_bf, FP_TYPE icost);
        track_point(const track_point& rhs) : pt(rhs.pt), timestep(rhs.timestep), rules_bf(rhs.rules_bf), cost(rhs.cost) {}
        steering_extremum pt;   // point
        int timestep;
        int rules_bf;
        FP_TYPE cost;
};

/*****************************************************************************/

class track
{
    friend class track_list;
    friend class tracker;
    public:
        track(void);

        // get and set the candidate point      
        void set_candidate_point(track_point cand_pt);
        track_point* get_candidate_point(void);
        
        // consolidate the candidate points - add the candidate points to the
        // end of the track
        void consolidate_candidate_point(void);

        // get various points on the track
        track_point* get_last_track_point(void);
        track_point* get_track_point(int idx);
        
        // get length / persistence / deviation / sum of curvature
        FP_TYPE get_length(void);
        int get_persistence(void);
        FP_TYPE get_deviation(void);
        FP_TYPE get_curvature_sum(void);
        FP_TYPE get_curvature_mean(void);
        FP_TYPE get_curvature_stddev(FP_TYPE mean=2e20);
        
        // get / set deleted
        void set_deleted();
        bool is_deleted();

    private:
        track_point cand_pt;
        std::vector<track_point> tr;
        bool deleted;
};

/*****************************************************************************/

class track_list
{
    public:
        // create
        track_list(void);
        int get_number_of_tracks(void);
        track* get_track(int track_n);
        
        // add and consolidate tracks
        void add_track(track& new_track);
        void consolidate_tracks(void);
        void prune_tracks(void);
        
        // save / load
        void save(std::string output_fname);
        void save_text(std::string output_fname, FP_TYPE sr);
        void load(std::string input_fname);
        
        // meta data for the track
        META_DATA_TYPE* get_meta_data(void);
        void set_meta_data(META_DATA_TYPE* meta_data);

    private:
        std::vector<track> tr_list;
        META_DATA_TYPE meta_data;
};

#endif
