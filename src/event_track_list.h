/******************************************************************************
** Program : event_track_list.h
** Author  : Neil Massey
** Date    : 14/04/14
** Purpose : class that contains track details for a (windstorm) event
******************************************************************************/

#ifndef EVENT_TRACK_LIST_H
#define EVENT_TRACK_LIST_H

#include "extremum.h"
#include <vector>

/*****************************************************************************/

struct event_point
{
	int timestep;
	steering_extremum* svex;
};

/*****************************************************************************/

class event_track
{
	public:
		event_track(void);
		void set_candidate_event(event_point evp, FP_TYPE cand_dist, 
								 FP_TYPE cand_overlap);
		event_point* get_candidate_event(void);
		FP_TYPE get_candidate_distance(void);
		FP_TYPE get_candidate_overlap(void);
		
		// adds candidate point to the end of the track
		void consolidate_candidate_event(void);
		
		// get / set
		event_point* get_last_event(void);
		std::vector<event_point>* get_track(void);
		
		// distance / persistence / deviation
		FP_TYPE get_length(void);
		int get_persistence(void);
		FP_TYPE get_deviation(void);		
		
	private:
		event_point cand_event;
		FP_TYPE cand_dist;
		FP_TYPE cand_overlap;
		
		std::vector<event_point> track;
};

/*****************************************************************************/

class event_track_list
{
	public:
		event_track_list(void);
		event_track* get_event_track(int ev_trk_n);
		int get_number_of_event_tracks(void);
		void add_event_track(event_track& new_event_track);
		void consolidate_event_tracks(void);
		
	private:
		std::vector<event_track> track_list;
};

/*****************************************************************************/

#endif