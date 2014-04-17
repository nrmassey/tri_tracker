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
const FP_TYPE mv = -2e20;

/*****************************************************************************/

std::string days_since_to_date_string(FP_TYPE days_since)
{
	// should be able to pass in base date, and calendar type, but we don't have time!
	int base_year  = 1959;
	int base_month = 12;
	int base_day   = 1;
	
	// get the number of years, months and days in the days_since
	FP_TYPE remainder;
	int year = int(days_since)/360;
	remainder = days_since - year * 360;
	int month = int(remainder)/30;
	remainder = remainder - month*30;
	int day = int(remainder);
	remainder = remainder-day;
	int hours = remainder*24;
	
	// add to the base date
	day += base_day;
	if (day > 30)
	{
		day = day - 30;
		month += 1;
	}
	month += base_month;
	if (month > 12)
	{
		month = month - 12;
		year += 1;
	}
	year += base_year;
	
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
	return ss.str();
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
		NcVar* nc_var = nc_file->get_var("time1");
		int t_len = nc_var->get_dim(0)->size();
		FP_TYPE* time_data = new FP_TYPE[t_len];
		nc_var->get(time_data, t_len);
		time_field.push_back(time_data);
		time_length.push_back(t_len);
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
								   FP_TYPE& max_ws, FP_TYPE& min_mslp)
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
			mslp_footprint[j*footprint_width+i] = mv;
			wind_footprint[j*footprint_width+i] = mv;
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
				if (wind_ws > wind_footprint[it_gi->j*footprint_width+it_gi->i] && 
					wind_ws > 0.0 && mslp > 0.0)
				{
					int idx = it_gi->j*footprint_width+it_gi->i;
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
			if (wind_footprint[j*footprint_width+i] > mv)
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
			if (wind_footprint[idx] > mv)
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
	calc_max_ws_min_mslp(e, evp, max_ws, min_mslp);
	calculate_wind_footprint(e, evp, n_timesteps, n_points, central_t_step);
	
	// open the file as a text file
	std::ofstream fh;
	fh.open(out_name.c_str());
	// get the value of the central timestep
	FP_TYPE t_val = get_time(central_t_step);
	fh << days_since_to_date_string(t_val) << "," << std::endl;
	fh << n_timesteps * hours_per_t_step << "," << std::endl;
	fh.precision(2);
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
			mslp_footprint[j*footprint_width+i] = mv;
			wind_footprint[j*footprint_width+i] = mv;
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
			if (wind_ws > wind_footprint[it_gi->j*footprint_width+it_gi->i] && 
				wind_ws > 0.0 && mslp > 0.0)
			{
				int idx = it_gi->j*footprint_width+it_gi->i;
				wind_footprint[idx] = wind_ws;
				mslp_footprint[idx] = mslp;
			}
		}
	}
	// count the number of points
	n_points = 0;
	for (int j=0; j<footprint_height; j++)
		for (int i=0; i<footprint_width; i++)
			if (wind_footprint[j*footprint_width+i] > mv)
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
	int evp;
	FP_TYPE max_ws, min_mslp;
	calc_max_ws_min_mslp(e, evp, max_ws, min_mslp);
	FP_TYPE start_t_val = get_time((*evt_track)[0].timestep);
	start_t_val = float(int(start_t_val*4))/4;
	
	// open the file as a text file
	std::ofstream fh;
	fh.open(out_name.c_str());
	// write the date out
	fh << days_since_to_date_string(start_t_val) << "," << std::endl;
	// write the time period out
	int n_timesteps = evt->get_persistence();
	fh <<  n_timesteps * hours_per_t_step << ", " << std::endl;
	// write the max windspeed and minimum mslp out
	fh.precision(2);
	fh << max_ws << ", " << std::endl;
	fh.precision(0);
	fh << min_mslp << ", " << std::endl;
	// write the number of timesteps out
	fh << n_timesteps << ", " << std::endl;
	
	// now loop over the track
	for (int p=0; p<evt_track->size(); p++)
	{
		// write the time out
		FP_TYPE t_val = start_t_val + p * 0.25;
		if (!(t_val > mv))	// bogus time value?
			t_val = get_time((*evt_track)[p].timestep+1);
		fh << days_since_to_date_string(t_val) << "," << std::endl;
		int n_points = 0;
		calculate_time_step_footprint(e, p, n_points);
		fh << n_points << ", " << std::endl;
		write_footprint(fh);
	}
	fh.close();
}