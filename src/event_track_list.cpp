/******************************************************************************
** Program : event_track_list.cpp
** Author  : Neil Massey
** Date    : 14/04/14
** Purpose : class that contains track details for a (windstorm) event
******************************************************************************/

#include "event_track_list.h"
#include "haversine.h"

/*****************************************************************************/
/** EVENT_TRACK                                                             **/
/*****************************************************************************/

event_track::event_track(void)
{
	// set overlap and distance to dummy values
	cand_dist = 2e20;
	cand_overlap = -1;
}

/*****************************************************************************/

void event_track::set_candidate_event(event_point evp, FP_TYPE icand_dist,
									  FP_TYPE icand_overlap)
{
	cand_event   = evp;
	cand_dist    = icand_dist;
	cand_overlap = icand_overlap;
}

/*****************************************************************************/

event_point* event_track::get_candidate_event(void)
{
	return &cand_event;
}

/*****************************************************************************/

FP_TYPE event_track::get_candidate_distance(void)
{
	return cand_dist;
}

/*****************************************************************************/

FP_TYPE event_track::get_candidate_overlap(void)
{
	return cand_overlap;
}

/*****************************************************************************/
		
void event_track::consolidate_candidate_event(void)
{
	if (cand_dist < 2e20)
		track.push_back(cand_event);
	// reset overlap and distance to dummy values
	cand_dist = 2e20;
	cand_overlap = -1;
}

/*****************************************************************************/
		
event_point* event_track::get_last_event(void)
{
	return &(track.back());
}

/*****************************************************************************/

std::vector<event_point>* event_track::get_track(void)
{
	return &track;
}

/*****************************************************************************/

// distance / persistence / deviation
FP_TYPE event_track::get_length(void)
{
	std::vector<event_point>::iterator it_pt1 = track.begin();
	std::vector<event_point>::iterator it_pt2 = track.begin();
	it_pt2 ++;
	
	FP_TYPE length = 0.0;
	
	while (it_pt2 != track.end())
	{
		steering_extremum* tp1 = it_pt1->svex;
		steering_extremum* tp2 = it_pt2->svex;
		length += haversine(tp1->lon, tp1->lat, tp2->lon, tp2->lat, EARTH_R);
		it_pt1 ++;
		it_pt2 ++;
	}
	return length;
}

/*****************************************************************************/

int event_track::get_persistence(void)
{
	return track.size();
}

/*****************************************************************************/

FP_TYPE event_track::get_deviation(void)
{
	steering_extremum* tp1 = track.front().svex;;
	steering_extremum* tp2 = track.back().svex;
	FP_TYPE d = haversine(tp1->lon, tp1->lat, tp2->lon, tp2->lat, EARTH_R);
	return d;
}

/*****************************************************************************/
/** EVENT_TRACK_LIST                                                        **/
/*****************************************************************************/
		
event_track_list::event_track_list(void)
{
}

/*****************************************************************************/

event_track* event_track_list::get_event_track(int ev_trk_n)
{
	return &(track_list[ev_trk_n]);
}

/*****************************************************************************/

int event_track_list::get_number_of_event_tracks(void)
{
	return track_list.size();
}

/*****************************************************************************/

void event_track_list::add_event_track(event_track& new_event_track)
{
	track_list.push_back(new_event_track);
}

/*****************************************************************************/

void event_track_list::consolidate_event_tracks(void)
{
	for (int e=0; e<track_list.size(); e++)
		track_list[e].consolidate_candidate_event();
}