/******************************************************************************
** Program : eventor.cpp
** Author  : Neil Massey
** Date    : 11/04/14
** Purpose : class to read in regridded data on a triangular mesh and 
**           construct an event set based on that data
******************************************************************************/

#include "eventor.h"
#include "haversine.h"
#include <sstream>
#include <netcdfcpp.h>
#include <ctime>
#include <math.h>
#include <algorithm>
const FP_TYPE w_mv = -2e20;

/*****************************************************************************/

std::string days_since_to_date_string(FP_TYPE days_since, int ref_year, 
									  int ref_month, int ref_day, FP_TYPE day_sc,
									  FP_TYPE ref_ndays_py)
{	
	// get the number of years, months and days in the days_since
	int year, month, day, hours;
	std::string date_string;
	// for a 360 day year we have to do the conversion ourselves
	if (ref_ndays_py == 360.0)
	{
		FP_TYPE remainder;
		// scale first to hours, etc.  Most of the time day_sc should be 1.0
		days_since = days_since * day_sc;
		year = int(days_since)/ref_ndays_py;
		remainder = days_since - year * ref_ndays_py;
		month = int(remainder)/30;
		remainder = remainder - month*30;
		day = int(remainder);
		remainder = remainder-day;
		hours = remainder*24;
	
		// add to the base date
		day += ref_day;
		if (day > 30)
		{
			day = day - 30;
			month += 1;
		}
		month += ref_month;
		if (month > 12)
		{
			month = month - 12;
			year += 1;
		}
		year += ref_year;
	
	}
	// if the calendar has a 365.25 day year then we can use c functions
	// mktime and localtime
	if (ref_ndays_py == 365.25)
	{
		// create the reference time object
		std::tm cur_time;
		memset(&cur_time, 0, sizeof(std::tm));
		cur_time.tm_year = ref_year;
		cur_time.tm_mon  = ref_month-1;
		cur_time.tm_mday = ref_day;
		std::mktime(&cur_time);
		// number of seconds per day
		const int n_secs_pd = 24 * 60 * 60;
		// add the number of days to the time
		cur_time.tm_mday += days_since * day_sc;
		cur_time.tm_hour = (days_since*day_sc - int(days_since*day_sc))*24;
		// correct the fields
		std::mktime(&cur_time);
		year  = cur_time.tm_year;
		month = cur_time.tm_mon+1;
		day   = cur_time.tm_mday;
		hours = cur_time.tm_hour;
	}
	// construct the string
	std::stringstream ss;
	ss.fill('0');
	ss.width(4);
	ss << year << "-";
	ss.width(2);
	ss << month << "-";
	ss.width(2);
	ss << day;
	ss << "T";
	ss.width(2);
	ss << hours << ":";
	ss.width(2);
	ss << 0 << ":";
	ss.width(2);
	ss << 0;
	date_string = ss.str();
	return date_string;
}

/*****************************************************************************/

eventor::eventor(std::vector<std::string> iinput_fname, int imin_per, 
				 FP_TYPE imin_len, FP_TYPE imin_dev, FP_TYPE isr,
				 int ievent_t_steps, FP_TYPE ioverlap, 
				 std::vector<std::string> mslp_file_name, std::string mslp_field_name,
				 std::vector<std::string> wind_file_name, std::string wind_field_name,
				 std::string lsm_file_name, std::string mesh_file) :
				 min_per(imin_per), min_len(imin_len), min_dev(imin_dev),
				 sr(isr), event_t_steps(ievent_t_steps), min_overlap(ioverlap),
				 lsm(NULL), mslp_footprint(NULL), wind_footprint(NULL)
{
	// initialise here but should be command line argument
	hours_per_t_step = 6;
	// check that there are enough netCDF files to go with the extrema files
	if (iinput_fname.size() != mslp_file_name.size() ||
		iinput_fname.size() != wind_file_name.size())
		throw (std::string("Number of input netCDF files do not equal the number of extrema files"));

	// load the extrema
	for (unsigned int i=0; i<iinput_fname.size(); i++)
	{
	    std::ifstream file_in(iinput_fname[i].c_str(), std::ios::in);
		std::cout << "# Reading data" << std::endl;
		if (!file_in.is_open())
			throw (std::string("File " + iinput_fname[i] + " could not be opened."));
		ex_list.load(iinput_fname[i], mv, i>0);
		file_in.close();
	}
	
	// load the MSLP netCDF files
	for (unsigned int i=0; i<mslp_file_name.size(); i++)
	{
		mslp_field.push_back(new ncdata(mslp_file_name[i], mslp_field_name));
	}

	// load the wind speed netCDF files
	for (unsigned int i=0; i<wind_file_name.size(); i++)
	{
		wind_speed_field.push_back(new ncdata(wind_file_name[i], wind_field_name));
	}
	
	// load the time netCDF files
	for (unsigned int i=0; i<mslp_file_name.size(); i++)
	{
		NcFile* nc_file = new NcFile(mslp_file_name[i].c_str());
		// get the time netCDF variable
		NcVar* nc_t_var = nc_file->get_var(mslp_field[0]->get_time_dim_name().c_str());
		int t_len = nc_t_var->get_dim(0)->size();
		FP_TYPE* time_data = new FP_TYPE[t_len];
		nc_t_var->get(time_data, t_len);
		time_field.push_back(time_data);
		time_length.push_back(t_len);
		// if this is the first netCDF file then get the reference time
		if (i==0)
			mslp_field[0]->get_reference_time(ref_year, ref_month, ref_day, ref_day_sc, ref_ndays_py);
	}
	
	// create the wind and mslp footprints
	footprint_width = mslp_field[0]->get_lon_len();
	footprint_height = mslp_field[0]->get_lat_len();
	wind_footprint = new FP_TYPE[footprint_height*footprint_width];
	mslp_footprint = new FP_TYPE[footprint_height*footprint_width];
	
	if (lsm_file_name != "")
		lsm = new ncdata(lsm_file_name, "lsm");
		
	// load the triangular mesh file
	tg.load(mesh_file);
}

/*****************************************************************************/

eventor::~eventor(void)
{
	for (int i=0; i<mslp_field.size(); i++)
		delete mslp_field[i];
	for (int i=0; i<wind_speed_field.size(); i++)
		delete wind_speed_field[i];
	for (int i=0; i<wind_speed_field.size(); i++)
		delete [] time_field[i];
	delete lsm;
	delete [] mslp_footprint;
	delete [] wind_footprint;
}

/*****************************************************************************/

bool eventor::evaluate_candidate_event(steering_extremum* trk_pt, 
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
	std::cout << "# Locating event tracks, timestep: ";
	build_first_frame();
	// loop through the other time steps and the events within them
	for (int t=1; t<ex_list.size(); t++)
	{
		std::cout << t;
        std::cout.flush();
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
			// repeat for as many events currently in the queue
			int qsize = ex_queue.size();
			for (int q=0; q<qsize; q++)
			{
				// pop an event off the front of the queue
				steering_extremum* svex_cand = ex_queue.front();
				ex_queue.pop();	
				// loop over every track so far
				for (int k=0; k<event_list.get_number_of_event_tracks(); k++)
				{
					// get the event track from the list
					event_track* evt = event_list.get_event_track(k);
					// get the last event point from the event_track
					event_point* evp = evt->get_last_event();
					// only assign to tracks where the last timestep was the previous timestep
					if (evp->timestep != t-1)
						continue;
					// otherwise calculate the distance and the overlap
					FP_TYPE svex_cand_dist, svex_cand_overlap;
					evaluate_candidate_event(evp->svex, svex_cand, svex_cand_dist, svex_cand_overlap);
					// check against maximum distance and minimum overlap
					if (svex_cand_dist > sr || svex_cand_overlap < min_overlap)
						continue;
					// get the current candidate event, its distance and overlap
					event_point* cur_cand_event = evt->get_candidate_event();
					FP_TYPE cur_cand_dist = evt->get_candidate_distance();
					
					// if the new candidate distance is less than the current candidate
					// distance then assign and put the old candidate distance back into
					// the queue
					if (svex_cand_dist < cur_cand_dist)
					{
						event_point cand_evp;
						cand_evp.timestep = t;
						cand_evp.svex = svex_cand;
						evt->set_candidate_event(cand_evp, svex_cand_dist, svex_cand_overlap);
						// if it exists then put the current candidate distance back into the queue
						if (cur_cand_dist < 2e20)
							ex_queue.push(cur_cand_event->svex);
						assigned = true;
					}
				}
				// if not assigned then push point back onto queue
				if (not assigned)
					ex_queue.push(svex_cand);				
			}
		}
		// consolidate the tracks that have been added
		event_list.consolidate_event_tracks();
		// unassigned points left?
		for (int q=0; q<ex_queue.size(); q++)
		{
			steering_extremum* svex = ex_queue.front();
			ex_queue.pop();	
			// create a new track with the point
			event_point evp;
			evp.timestep = t;
			evp.svex = svex;
			// make an event_track and add the event
			event_track evt;
			evt.set_candidate_event(evp, 0.0, 0.0);
			evt.consolidate_candidate_event();
			// add the event_track to the event_track_list
			event_list.add_event_track(evt);			
		}
        int e = t;
        if (t == 0)
            std::cout << "\b";
        while (e > 0)
        {
            e = e / 10;
            std::cout << "\b";
        }		
	}
	std::cout << std::endl;
	// merge any events that occur over the same date, at the same location
	while (merge_events()){}
	// remove any short / short lived / slow moving tracks
	trim_tracks();
}

/*****************************************************************************/

void eventor::save(std::string output_fname)
{
	// write the 72 hour footprints
	std::cout << "# Number of events : " << event_list.get_number_of_event_tracks() << std::endl;
	std::cout << "# Calculating 72 hour footprints : ";
	for (int e=0; e<event_list.get_number_of_event_tracks(); e++)
	{
		std::cout << e;
        std::cout.flush();
		// output name construction
		std::stringstream ss;
		ss.width(2);
		ss.fill('0');
		ss << e;
		std::string out_72h_name = output_fname + "_e" + ss.str() + "_72h.csv";
		save_72hour_footprint(out_72h_name, e);
        int p = e;
        if (e == 0)
            std::cout << "\b";
        while (p > 0)
        {
            p = p / 10;
            std::cout << "\b";
        }
	}
	std::cout << std::endl;
	
	// write the time varying footprints
	std::cout << "# Calculating time varying footprints : ";	
	for (int e=0; e<event_list.get_number_of_event_tracks(); e++)
	{
		std::cout << e;
		std::cout.flush();
        int p = e;
        if (e == 0)
            std::cout << "\b";
        while (p > 0)
        {
            p = p / 10;
            std::cout << "\b";
        }		
		std::stringstream ss;
		ss.width(2);
		ss.fill('0');		
		ss << e;
		std::string out_tv_name = output_fname + "_e" + ss.str() + "_tv.csv";
		save_time_varying_footprint(out_tv_name, e);
	}
	std::cout << std::endl;
}

/*****************************************************************************/

void eventor::build_first_frame(void)
{
	const int ff = 0;		// first frame
	for (int e=0; e<ex_list.number_of_extrema(ff); e++)
	{
		// create a new track for each extremum in the first time step
		event_point evp;
		evp.timestep = ff;
		evp.svex = ex_list.get(ff, e);
		// make an event_track and add the event
		event_track evt;
		evt.set_candidate_event(evp, 0.0, 0.0);
		evt.consolidate_candidate_event();
		// add the event_track to the event_track_list
		event_list.add_event_track(evt);
	}
}

/*****************************************************************************/

void eventor::trim_tracks(void)
{
	// reduce the tracks based on minimum persistence, deviation and length
	event_track_list new_event_list;
	for (int e=0; e<event_list.get_number_of_event_tracks(); e++)
	{
		event_track* evt = event_list.get_event_track(e);
		bool add = true;
		add = add && evt->get_persistence() >= min_per;
		add = add && evt->get_length() >= min_len;
		add = add && evt->get_deviation() >= min_dev;

		if (add)
			new_event_list.add_event_track(*evt);
	}
	event_list = new_event_list;
}

/*****************************************************************************/

bool event_times_overlap(std::vector<event_point>* evt_trk_0,
						 std::vector<event_point>* evt_trk_1)
{
	int evt_0s = (*evt_trk_0)[0].timestep;
	int evt_0e = evt_0s + evt_trk_0->size();
	int evt_1s = (*evt_trk_1)[0].timestep;
	int evt_1e = evt_1s + evt_trk_1->size();
	
	// six possible overlap cases:
	// e0:   --   ----  ---    ---     ---  ---
	// e1:  ----   --    ---  ---   ---        ---
	bool overlap = false;
	if (evt_0s >= evt_1s && evt_0e <= evt_1e)
		overlap = true;
	if (evt_0s <= evt_1s && evt_0e >= evt_1e)
		overlap = true;
	if (evt_0s <= evt_1s && evt_0e >= evt_1s)
		overlap = true;
	if (evt_0s >= evt_1s && evt_0e <= evt_1s)
		overlap = true;
//	if (evt_0e+1 == evt_1s)
//		overlap = true;
//	if (evt_0s == evt_1e+1)
//		overlap = true;
	return overlap;
}						 

/*****************************************************************************/

FP_TYPE calculate_overlap_percentage(std::vector<event_point>* evt_pts_0,
									 std::vector<event_point>* evt_pts_1)
{
	// calculate by how much two events overlap during the timesteps which
	// they do overlap
	LABEL_STORE e0_labels;
	LABEL_STORE e1_labels;
	
	// build two lists of the labels in each object
	for (int e0=0; e0 < evt_pts_0->size(); e0++)
	{
		LABEL_STORE svex_labs_e0 = (*evt_pts_0)[e0].svex->object_labels;
		for (LABEL_STORE::iterator it_e0 = svex_labs_e0.begin();
			 it_e0 != svex_labs_e0.end(); it_e0++)
		{
			if (std::find(e0_labels.begin(), e0_labels.end(), *it_e0) ==
						  e0_labels.end())
				e0_labels.push_back(*it_e0);
		}
	}
	
	for (int e1=0; e1 < evt_pts_1->size(); e1++)
	{
		LABEL_STORE svex_labs_e1 = (*evt_pts_1)[e1].svex->object_labels;
		for (LABEL_STORE::iterator it_e1 = svex_labs_e1.begin();
			 it_e1 != svex_labs_e1.end(); it_e1++)
		{
			if (std::find(e1_labels.begin(), e1_labels.end(), *it_e1) ==
						  e1_labels.end())
				e1_labels.push_back(*it_e1);
		}
	}
	
	// now calculate the overlap between the two lists
	FP_TYPE overlap_e1 = 0;
	for (LABEL_STORE::iterator it_e1 = e1_labels.begin(); it_e1 != e1_labels.end();
		 it_e1 ++)
	{
		if (std::find(e0_labels.begin(), e0_labels.end(), *it_e1) != e0_labels.end())
			overlap_e1 += 1;
	}

	FP_TYPE overlap_e0 = 0;
	for (LABEL_STORE::iterator it_e0 = e0_labels.begin(); it_e0 != e0_labels.end();
		 it_e0 ++)
	{
		if (std::find(e1_labels.begin(), e1_labels.end(), *it_e0) != e1_labels.end())
			overlap_e0 += 1;
	}
	
	overlap_e1 = overlap_e1 / e1_labels.size() * 100.0;
	overlap_e0 = overlap_e0 / e0_labels.size() * 100.0;
	
	if (overlap_e1 > overlap_e0)
		return overlap_e1;
	else
		return overlap_e0;
}

/*****************************************************************************/

FP_TYPE calculate_distance_over_track(std::vector<event_point>* evt_pts_0,
									  std::vector<event_point>* evt_pts_1)
{
	FP_TYPE mean_dist = 0.0;
	int n_pts = 0;
	for (int e0=0; e0 < evt_pts_0->size(); e0++)
	{
		FP_TYPE lon_0 = (*evt_pts_0)[e0].svex->lon;
		FP_TYPE lat_0 = (*evt_pts_0)[e0].svex->lat;
		int evt_0 = (*evt_pts_0)[e0].timestep;
		
		for (int e1=0; e1 < evt_pts_1->size(); e1++)
		{
			int evt_1 = (*evt_pts_0)[e0].timestep;
			if (evt_0 == evt_1)
			{
				FP_TYPE lon_1 = (*evt_pts_1)[e1].svex->lon;
				FP_TYPE lat_1 = (*evt_pts_1)[e1].svex->lat;
				FP_TYPE dist = haversine(lon_0, lat_0, lon_1, lat_1, EARTH_R);
				mean_dist += dist / 1000.0;
				n_pts += 1;
			}
		}
	}
	if (n_pts == 0)
		return 2e20;
	else
		return mean_dist / n_pts;
}

/*****************************************************************************/


void merge_event_points(std::vector<event_point>* evt_pts_0,
						std::vector<event_point>* evt_pts_1)
{
	// three ways to merge evt_pts_1 into evt_pts_0
	// 1 - the timestep is less than the first in evt_pts_0 - add the event point to the
	//     beginning of evt_pts_0
	// 2 - the timestep is greater than the last in evt_pts_0 - add the event point to the
	//     end of evt_pts_1
	// 3 - the timestep is contained within the timesteps in evt_pts_0 - add the object
	//     labels from evt_pts_1 to evt_pts_0
	int evt_0s = (*evt_pts_0)[0].timestep;
	int evt_0e = evt_0s + evt_pts_0->size();
	int evt_1s = (*evt_pts_0)[1].timestep;
	
	// this operation will change the size of the vector so do the 3rd case afterwards
	for (int e1=0; e1 < evt_pts_1->size(); e1++)
	{
		int evt_1 = (*evt_pts_1)[e1].timestep;
		if (evt_1 < evt_0s)
			evt_pts_0->insert(evt_pts_0->begin(), (*evt_pts_1)[e1]);
		if (evt_1 > evt_0e)
			evt_pts_0->push_back((*evt_pts_1)[e1]);
	}
	
	// now merge the objects which have timesteps between evt_0s and evt_0e
	for (int e1=0; e1 < evt_pts_1->size(); e1++)
	{
		int evt_1 = (*evt_pts_1)[e1].timestep;
		if (evt_1 >= evt_0s && evt_1 <= evt_0e)
		{
			// get the index into the evt_pts_0 list
			int evt_idx = -1;
			for (int e0=0; e0 < evt_pts_0->size(); e0++)
				if ((*evt_pts_0)[e0].timestep == evt_1)
					evt_idx = e0;
			if (evt_idx == -1)
				continue;
			// get the two steering extremums								
			steering_extremum* svex_e0 = (*evt_pts_0)[evt_idx].svex;
			steering_extremum* svex_e1 = (*evt_pts_1)[e1].svex;
			// now merge the labels
			for (LABEL_STORE::iterator it_e1_labs = svex_e1->object_labels.begin();
				 it_e1_labs != svex_e1->object_labels.end(); it_e1_labs++)
			{
				if (std::find(svex_e0->object_labels.begin(), 
							  svex_e0->object_labels.end(),
							  *it_e1_labs) == svex_e0->object_labels.end())
				{
					svex_e0->object_labels.push_back(*it_e1_labs);
				}
			}
		}
	}
}

/*****************************************************************************/

bool eventor::merge_events(void)
{
	// merge the events together based on three criteria:
	// 1. the date - the events must overlap
	// 2. the distance between the minimum MSLP positions (from the triangle centroid)
	// 3. the value of the MSLP minima
	// 4. the value of the maximum wind speed
	
	// consider every combination of event
	std::cout << "# Merging events : " << event_list.get_number_of_event_tracks() << " events" << std::endl;
	bool merged=false;
	for (int e0=0; e0<event_list.get_number_of_event_tracks(); e0++)
	{
		// get the event track from the event list
		event_track* evt_0 = event_list.get_event_track(e0);
		std::vector<event_point>* evt_pts_0 = evt_0->get_track();
		if (evt_pts_0->size() == 0)
			continue;
		
		for (int e1=e0+1; e1<event_list.get_number_of_event_tracks(); e1++)
		{
			// don't merge the same
			if (e0 == e1)
				continue;
			// get the event track from the event list
			event_track* evt_1 = event_list.get_event_track(e1);
			std::vector<event_point>* evt_pts_1 = evt_1->get_track();
			// already merged objects have had tracks cleared!
			if (evt_pts_1->size() == 0)
				continue;
			bool overlap = event_times_overlap(evt_pts_0, evt_pts_1);
			// get the min mslp and max windspeed for this event
			if (overlap)
			{
				FP_TYPE overlap_pct = calculate_overlap_percentage(evt_pts_0, evt_pts_1);
				FP_TYPE mean_dist = calculate_distance_over_track(evt_pts_0, evt_pts_1);
				if (overlap_pct >= 85.0)
				{
					merge_event_points(evt_pts_0, evt_pts_1);
					merged=true;
					evt_pts_1->clear();
				}
			}
			else
			{
				int evt_0e = (*evt_pts_0)[0].timestep + evt_pts_0->size();
				int evt_1s = (*evt_pts_1)[0].timestep;
				if (evt_0e == evt_1s-1 && 
					(evt_pts_1->size() < event_t_steps || evt_pts_0->size() < event_t_steps))
				{
					FP_TYPE overlap_pct = calculate_overlap_percentage(evt_pts_0, evt_pts_1);
					if (overlap_pct >= 85.0)
					{
						merge_event_points(evt_pts_0, evt_pts_1);
						merged=true;
						evt_pts_1->clear();
					}
				
				}
			}
		}
	}
	// consolidate the events
	event_track_list new_event_list;
	for (int e0=0; e0<event_list.get_number_of_event_tracks(); e0++)
	{
		// get the event track from the event list
		event_track* evt_0 = event_list.get_event_track(e0);
		std::vector<event_point>* evt_pts_0 = evt_0->get_track();
		if (evt_pts_0->size() == 0)
			continue;
		new_event_list.add_event_track(*evt_0);		
	}
	event_list = new_event_list;
	return merged;
}

/*****************************************************************************/

FP_TYPE eventor::get_mslp_data(int t, int y, int x)
{
	// which netCDF data should we use?
	ncdata* this_data = NULL;
	int c_offset = 0;
	for (int nci=0; nci<mslp_field.size(); nci++)
	{
		if (t >= mslp_field[nci]->get_t_len()+c_offset)
			c_offset += mslp_field[nci]->get_t_len();
		else
		{
			this_data = mslp_field[nci];
			break;
		}
	}
	FP_TYPE retval = mv;
	if (this_data != NULL)
	{
		if (t-c_offset < this_data->get_t_len())
			retval = this_data->get_data(x,y,0,t-c_offset);
	}
	return retval;
}

/*****************************************************************************/

FP_TYPE eventor::get_wind_speed_data(int t, int y, int x)
{
	// which netCDF data should we use?
	ncdata* this_data = NULL;
	int c_offset = 0;
	// check the indexes against the size of the arrays
	if (y >= wind_speed_field[0]->get_lat_len() ||
	    x >= wind_speed_field[0]->get_lon_len())
		return mv;
	for (int nci=0; nci<wind_speed_field.size(); nci++)
	{
		if (t >= wind_speed_field[nci]->get_t_len()+c_offset)
			c_offset += wind_speed_field[nci]->get_t_len();
		else
		{
			this_data = wind_speed_field[nci];
			break;
		}
	}
	FP_TYPE retval = mv;
	if (this_data != NULL)
	{
		if (t-c_offset < this_data->get_t_len())
			retval = this_data->get_data(x,y,0,t-c_offset);
	}
	return retval;
}

/*****************************************************************************/

FP_TYPE eventor::get_time(int t)
{
	int c_offset = 0;
	int this_nci = -1;
	for (int nci=0; nci<time_field.size(); nci++)
	{
		if (t >= time_length[nci]+c_offset)
			c_offset += time_length[nci];
		else
		{
			this_nci=nci;
			break;
		}
	}
	FP_TYPE retval = mv;
	if (this_nci != -1)
	{
		if (t-c_offset < time_length[this_nci])
		{
			retval = time_field[this_nci][t-c_offset];
		}
	}
	return retval;
}

/*****************************************************************************/

FP_TYPE eventor::get_lsm(int y, int x)
{
	return lsm->get_data(x,y,0,0);
}

/*****************************************************************************/

void eventor::calc_max_ws_min_mslp(int event_track_number, int& event_point_index, 
								   FP_TYPE& max_ws, FP_TYPE& min_mslp, int& t_step)
{
	// get the event track
	event_track* evt = event_list.get_event_track(event_track_number);
	max_ws = -1.0;	// start with a dummy value
	min_mslp = 2e20;
	int current_evpi = 0;
	// loop over the track
	for (std::vector<event_point>::iterator it_p = evt->get_track()->begin();
		 it_p != evt->get_track()->end(); it_p++)
	{
		// loop over all the labels in the extremum
		for (LABEL_STORE::iterator it_l = it_p->svex->object_labels.begin();
			 it_l != it_p->svex->object_labels.end(); it_l++)
		{
			// get the indices from the triangular grid
			indexed_force_tri_3D* tri = tg.get_triangle(*it_l);
			const std::list<grid_index>* grid_indices = tri->get_grid_indices();
			for (std::list<grid_index>::const_iterator it_gi = grid_indices->begin();
				 it_gi != grid_indices->end(); it_gi++)
			{
				FP_TYPE wind_ws = get_wind_speed_data(it_p->timestep, it_gi->j, it_gi->i);
				// do we need to multiply by the LSM to get values over land only?
				if (lsm != NULL)
				{
					FP_TYPE lsm_val = get_lsm(it_gi->j, it_gi->i);
					wind_ws *= lsm_val;
				}
				if (wind_ws > max_ws)
				{
					max_ws = wind_ws;
					event_point_index = current_evpi;
					t_step = it_p->timestep;
				}
				FP_TYPE mslp_ws = get_mslp_data(it_p->timestep, it_gi->j, it_gi->i);
				if (mslp_ws < min_mslp && mslp_ws > 0.0)
					min_mslp = mslp_ws;
			}
		}
		current_evpi++;
	}
}

/*****************************************************************************/

void eventor::calculate_wind_footprint(int event_track_number, int event_point_index,
									   int& n_timesteps, int& n_points, int& central_t_step)
{
	// first clear the footprints
	for (int j=0; j<footprint_height; j++)
		for (int i=0; i<footprint_width; i++)
		{
			mslp_footprint[j*footprint_width+i] = w_mv;
			wind_footprint[j*footprint_width+i] = w_mv;
		}
	// get the event track
	event_track* evt = event_list.get_event_track(event_track_number);
	int trk_p = evt->get_persistence();
	// get the start and end points of the 72 hour period
	int trk_start_idx = event_point_index - event_t_steps/2;
	int trk_end_idx = event_point_index + event_t_steps/2;
	// check if these are over / under the persistence / 0
	if (trk_end_idx >= trk_p)
	{
		trk_start_idx = trk_p - event_t_steps;
		trk_end_idx = trk_p;
	}
	if (trk_start_idx < 0)
		trk_start_idx = 0;
	if (trk_end_idx - trk_start_idx  < event_t_steps && 
		trk_start_idx + event_t_steps -1 < trk_p)
	{
		trk_end_idx = trk_start_idx + event_t_steps;
	}
	// loop over the track
	std::vector<event_point>* evt_track = evt->get_track();
	for (int t=trk_start_idx; t<trk_end_idx; t++)
	{
		// get the svex from the event point
		steering_extremum* svex = (*evt_track)[t].svex;
		int t_step = (*evt_track)[t].timestep;
		// loop over all the labels in the extremum
		for (LABEL_STORE::iterator it_l = svex->object_labels.begin();
			 it_l != svex->object_labels.end(); it_l++)
		{
			// get the indices from the triangular grid
			indexed_force_tri_3D* tri = tg.get_triangle(*it_l);
			const std::list<grid_index>* grid_indices = tri->get_grid_indices();
			for (std::list<grid_index>::const_iterator it_gi = grid_indices->begin();
				 it_gi != grid_indices->end(); it_gi++)
			{
				FP_TYPE wind_ws = get_wind_speed_data(t_step, it_gi->j, it_gi->i);
				FP_TYPE mslp = get_mslp_data(t_step, it_gi->j, it_gi->i);
				// check if this wind is greater than the current wind in the footprint at
				// these indices
				int idx = it_gi->j*footprint_width+it_gi->i;
				if (wind_ws > wind_footprint[idx] && 
					wind_ws > 0.0 && mslp > 0.0)
				{
					wind_footprint[idx] = wind_ws;
					mslp_footprint[idx] = mslp;
				}
			}
		}
	}
	// calculate the central timestep
	central_t_step = (*evt_track)[trk_start_idx].timestep + (trk_end_idx - trk_start_idx)/2;
	
	// calculate number of timesteps
	n_timesteps = trk_end_idx - trk_start_idx;
	
	// count the number of points
	n_points = 0;
	for (int j=0; j<footprint_height; j++)
		for (int i=0; i<footprint_width; i++)
			if (wind_footprint[j*footprint_width+i] > w_mv)
				n_points++;
}

/*****************************************************************************/

void eventor::write_footprint(std::ofstream& fh)
{
	// write the footprint out - write values out that are not mv
	for (int j=0; j<footprint_height; j++)
		for (int i=0; i<footprint_width; i++)
		{
			int idx = j*footprint_width+i;
			if (wind_footprint[idx] > w_mv)
			{
				// get the true latitude / longitude from the rotated grid
				FP_TYPE tru_lon, tru_lat;
				tru_lon = mslp_field[0]->get_rotated_grid()->get_global_longitude_value(i,j);
				tru_lat = mslp_field[0]->get_rotated_grid()->get_global_latitude_value(i,j);
				// get the grid box corner points
				std::list<FP_TYPE> gb_corners = mslp_field[0]->get_rotated_grid()->get_global_grid_box_values(i,j);
				// write the data
				fh.precision(2);
				fh << std::fixed;
				fh << tru_lon << "," << tru_lat << ",";
				for (std::list<FP_TYPE>::iterator it_gbc = gb_corners.begin(); 
					 it_gbc != gb_corners.end(); it_gbc++)
					fh << *it_gbc << ",";		// order is lon0,lat0, lon1,lat1, lon2,lat2, lon3,lat3
				fh.precision(0);
				fh << mslp_footprint[idx] << ",";
				fh.precision(2);
				fh << wind_footprint[idx] << "," << std::endl;
			}
		}
}

/*****************************************************************************/

void eventor::save_72hour_footprint(std::string out_name, int e)
{
	// Calculate the wind footprints, integrated over a maximum of 72hours
	// and save them
	FP_TYPE max_ws, min_mslp;
	int evp;
	int n_timesteps = 0;
	int n_points = 0;
	int central_t_step;
	int max_tstep = 0;
	calc_max_ws_min_mslp(e, evp, max_ws, min_mslp, max_tstep);
	calculate_wind_footprint(e, evp, n_timesteps, n_points, central_t_step);
	
	// open the file as a text file
	std::ofstream fh;
	fh.open(out_name.c_str());
	// get the value of the central timestep
	FP_TYPE t_val = get_time(max_tstep);
	fh << days_since_to_date_string(t_val, ref_year, ref_month, ref_day, ref_day_sc, ref_ndays_py) << "," << std::endl;
	fh << n_timesteps * hours_per_t_step << "," << std::endl;
	event_track* evt = event_list.get_event_track(e);
//	std::vector<event_point>* evt_track = evt->get_track();
//	fh << (*evt_track)[0].timestep << " " << (*evt_track)[evt_track->size()-1].timestep << "," << std::endl;	
	fh.precision(2);
	fh << std::fixed;
	fh << max_ws << "," << std::endl;
	fh.precision(0);
	fh << min_mslp << "," << std::endl;
	fh << n_points << "," << std::endl;
	write_footprint(fh);
	fh.close();
}

/*****************************************************************************/

void eventor::calculate_time_step_footprint(int event_track_number, int event_point_index,
					 			 			int& n_points)
{
	// first clear the footprints
	for (int j=0; j<footprint_height; j++)
		for (int i=0; i<footprint_width; i++)
		{
			mslp_footprint[j*footprint_width+i] = w_mv;
			wind_footprint[j*footprint_width+i] = w_mv;
		}
	// get the event track
	event_track* evt = event_list.get_event_track(event_track_number);
	// get the svex for the trackpoint
	event_point evp = (*(evt->get_track()))[event_point_index];
	for (LABEL_STORE::iterator it_l = evp.svex->object_labels.begin();
		 it_l != evp.svex->object_labels.end(); it_l++)
	{
		// get the indices from the triangular grid
		indexed_force_tri_3D* tri = tg.get_triangle(*it_l);
		const std::list<grid_index>* grid_indices = tri->get_grid_indices();
		for (std::list<grid_index>::const_iterator it_gi = grid_indices->begin();
			 it_gi != grid_indices->end(); it_gi++)
		{
			FP_TYPE wind_ws = get_wind_speed_data(evp.timestep, it_gi->j, it_gi->i);
			FP_TYPE mslp = get_mslp_data(evp.timestep, it_gi->j, it_gi->i);
			// check if this wind is greater than the current wind in the footprint at
			// these indices			
			int idx = it_gi->j*footprint_width+it_gi->i;
			if (wind_ws > wind_footprint[idx] && 
				wind_ws > 0.0 && mslp > 0.0)
			{
				wind_footprint[idx] = wind_ws;
				mslp_footprint[idx] = mslp;
			}
		}
	}
	// count the number of points
	n_points = 0;
	for (int j=0; j<footprint_height; j++)
		for (int i=0; i<footprint_width; i++)
			if (wind_footprint[j*footprint_width+i] > w_mv)
				n_points++;	
}

/*****************************************************************************/

void eventor::save_time_varying_footprint(std::string out_name, int e)
{
	// write out the timevarying footprints - i.e. the footprint at every timestep
	// over the lifetime of the event
	event_track* evt = event_list.get_event_track(e);
	std::vector<event_point>* evt_track = evt->get_track();

	// get the min mslp, max windspeed and time at which the track starts
	int evp, max_tstep;
	FP_TYPE max_ws, min_mslp;
	calc_max_ws_min_mslp(e, evp, max_ws, min_mslp, max_tstep);
	FP_TYPE start_t_val = get_time((*evt_track)[0].timestep);
	start_t_val = float(int(start_t_val*4))/4;
	
	// open the file as a text file
	std::ofstream fh;
	fh.open(out_name.c_str());
	// write the date out
	fh << days_since_to_date_string(start_t_val, ref_year, ref_month, ref_day, ref_day_sc, ref_ndays_py) << "," << std::endl;
	// write the time period out
	int n_timesteps = evt->get_persistence();
	fh <<  n_timesteps * hours_per_t_step << "," << std::endl;
	// write the max windspeed and minimum mslp out
	fh.precision(2);
	fh << std::fixed;
	fh << max_ws << "," << std::endl;
	fh.precision(0);
	fh << min_mslp << "," << std::endl;
	// write the number of timesteps out
	fh << n_timesteps << "," << std::endl;
	
	// now loop over the track
	for (int p=0; p<evt_track->size(); p++)
	{
		// write the time out
		FP_TYPE t_val = start_t_val + p * 0.25;
		fh << days_since_to_date_string(t_val, ref_year, ref_month, ref_day, ref_day_sc, ref_ndays_py) << "," << std::endl;
		int n_points = 0;
		calculate_time_step_footprint(e, p, n_points);
		fh.precision(2);
		fh << std::fixed;
		// write out the lon, lat, intensity and delta of the pressure minimum
		steering_extremum* svex = (*evt_track)[p].svex;
		fh << svex->lon << "," << svex->lat << ",";
		fh.precision(0);
		fh << svex->intensity << "," << svex->delta << std::endl;
		fh << n_points << "," << std::endl;
		write_footprint(fh);
	}
	fh.close();
}