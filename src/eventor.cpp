/******************************************************************************
** Program : eventor.cpp
** Author  : Neil Massey
** Date    : 11/04/14
** Purpose : class to read in regridded data on a triangular mesh and 
**           construct an event set based on that data
******************************************************************************/

#include "eventor.h"
#include "haversine.h"

/*****************************************************************************/

eventor::eventor(std::vector<std::string> iinput_fname, int imin_per, 
				 FP_TYPE imin_len, FP_TYPE imin_dev, FP_TYPE isr,
				 int ievent_t_steps, FP_TYPE ioverlap) :
				 min_per(imin_per), min_len(imin_len), min_dev(imin_dev),
				 sr(isr), event_t_steps(ievent_t_steps), min_overlap(ioverlap)
{
	for (unsigned int i=0; i<iinput_fname.size(); i++)
	{
	    std::ifstream file_in(iinput_fname[i].c_str(), std::ios::in);
		std::cout << "# Reading data" << std::endl;
		if (!file_in.is_open())
			throw ("File " + iinput_fname[i] + " could not be opened.");
		ex_list.load(iinput_fname[i], mv, i>0);
		file_in.close();
	}
}

/*****************************************************************************/

eventor::~eventor(void)
{
}

/*****************************************************************************/

bool eventor::evaluate_candidate_point(steering_extremum* trk_pt, 
									   steering_extremum* cand_pt,
									   FP_TYPE& dist, FP_TYPE& overlap)
{
	// evaluate the candidate point against the last point in the track
	// calculate distance
	dist = haversine(trk_pt->lon, trk_pt->lat, cand_pt->lon, cand_pt->lat, EARTH_R);
	// calculate overlap
	FP_TYPE n_overlapping_labels = 0.0;
	for (LABEL_STORE::iterator trk_label = trk_pt->object_labels.begin();
		 trk_label != trk_pt->object_labels.end(); trk_label++)
	{
		if (std::find(cand_pt->object_labels.begin(), cand_pt->object_labels.end(), 
				*trk_label) != cand_pt->object_labels.end())
		{
			n_overlapping_labels++;
		}
	}
	overlap = 100 * n_overlapping_labels / trk_pt->object_labels.size();
	// logic check - if percentage overlap is >= min_overlap and dist <= search radius
	if (dist <= sr && overlap >= min_overlap)
		return true;
	else
		return false;
}

/*****************************************************************************/

void eventor::clear_queue(void)
{
	while (!ex_queue.empty())
		ex_queue.pop();
}

/*****************************************************************************/

void eventor::find_events(void)
{
	build_first_frame();
	// loop through the other time steps and the events within them
	for (int t=0; t<ex_list.size(); t++)
	{
		// clear the queue of extremas yet to be added
		clear_queue();
		for (int e=0; e<ex_list.number_of_extrema(t); e++)
		{
			steering_extremum* svex_cand = ex_list.get(t, e);
			ex_queue.push(svex_cand);
		}
		// loop until no assignment has been made
		bool assigned = true;
		while (assigned)
		{
			assigned = false;
			// pop an event off the front of the queue
			steering_extremum* svex_cand = ex_queue.front();
			ex_queue.pop();						
		}
	}
}

/*****************************************************************************/

void eventor::save(std::string output_fname)
{
}

/*****************************************************************************/

void eventor::build_first_frame(void)
{
	const int ff = 95;
	std::cout << ex_list.number_of_extrema(ff) << std::endl;
	for (int e=0; e<ex_list.number_of_extrema(ff); e++)
	{
		// create a new track for each extremum in the first time step
		event_point evp;
		evp.t_step = ff;
		evp.svex = ex_list.get(ff, e);
		// make a list and add to the list
		std::list<event_point> local_event_list;
		local_event_list.push_back(evp);
		event_list.push_back(local_event_list);
	}
}