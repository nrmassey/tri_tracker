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

/*****************************************************************************/

struct event_point
{
	int t_step;
	steering_extremum* svex;
};

typedef std::list<std::list<event_point> > EVENT_LIST_TYPE;

/*****************************************************************************/

class eventor
{
	public:
		eventor(std::vector<std::string> iinput_fname, int imin_per, 
				FP_TYPE imin_len, FP_TYPE imin_dev, FP_TYPE isr,
				int ievent_t_steps, FP_TYPE ioverlap);
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
		
		// extrema related files
		extrema_list ex_list;
		FP_TYPE mv;				// missing value
		
		// functions to build the event set
		void build_first_frame(void);
		void build_event_set(void);
		bool evaluate_candidate_point(steering_extremum* trk_pt, steering_extremum* cand_pt,
									  FP_TYPE& dist, FP_TYPE& overlap);
		void clear_queue(void);
		// queue to assign extremas
		std::queue<steering_extremum*> ex_queue;
		
		// store data
		EVENT_LIST_TYPE event_list;
};

#endif