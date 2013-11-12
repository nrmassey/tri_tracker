/******************************************************************************
** Program : track_list.cpp
** Author  : Neil Massey
** Date    : 18/07/13
** Purpose : class to hold track data + io operators, regional version
******************************************************************************/

#include "track_list.h"
#include <math.h>
#include <netcdfcpp.h>
#include "haversine.h"
#include "geo_convert.h"
#include "vector_3D.h"
#include "get_bearing.h"
#include "spline.h"
#include "bin_file_utils.h"

/*****************************************************************************/

track_point::track_point(void) : frame_number(0) {}

/*****************************************************************************/

track_point::track_point(FP_TYPE iframe) : frame_number(iframe) {}

/*****************************************************************************/

void track_point::set_point(steering_extremum ipt, FP_TYPE icost, 
							FP_TYPE ic0, FP_TYPE ic1, FP_TYPE ic2, FP_TYPE ic3)
{
    pt = ipt;
	cost = icost;
	c0 = ic0;
	c1 = ic1;
	c2 = ic2;
	c3 = ic3;
}

/*****************************************************************************/

steering_extremum* track_point::get_point(void)
{
	return &pt;
}

/*****************************************************************************/

const FP_TYPE track_point::get_frame_number(void) const
{
	return frame_number;
}

/*****************************************************************************/

track::track(void) : length(-1.0), cand_cost(2e20f)
{
	tr.clear();
}
/*****************************************************************************/

void track::clear(void)
{
	tr.clear();
	length = -1.0;
	cand_cost = 2e20f;
}

/*****************************************************************************/

void track::add_point(track_point tp)
{
	tr.push_back(tp);
}

/*****************************************************************************/

FP_TYPE track::get_length(void)
{
	if (length == -1)	// cache the value
	{
		length = 0.0;
		for (unsigned int i=1; i<tr.size(); i++)
		{
			extremum* tp1 = tr[i].get_point();
			extremum* tp0 = tr[i-1].get_point();
			length += haversine(tp0->lon, tp0->lat, tp1->lon, tp1->lat, EARTH_R);
		}
	}
	return length;
}

/*****************************************************************************/

FP_TYPE track::get_deviation(void)
{
	extremum* tp0 = tr[0].get_point();
	extremum* tp1 = tr[tr.size()-1].get_point();
	FP_TYPE d = haversine(tp0->lon, tp0->lat, tp1->lon, tp1->lat, EARTH_R);
	return d;
}

/*****************************************************************************/

void track::reduce(void)
{
	// interpolate between stationary points
	for (unsigned int tr_i=1; tr_i<tr.size(); tr_i++)
	{
		if(tr[tr_i].get_point()->lon == tr[tr_i-1].get_point()->lon &&
           tr[tr_i].get_point()->lat == tr[tr_i-1].get_point()->lat)
		{
			// copy point
			steering_extremum new_geo = *(tr[tr_i].get_point());
			// interpolate the point to be between the two neighbouring points
			if (tr_i < tr.size()-1)
			{
				new_geo.lon = (tr[tr_i-1].get_point()->lon + tr[tr_i+1].get_point()->lon)/2.0;
				new_geo.lat = (tr[tr_i-1].get_point()->lat + tr[tr_i+1].get_point()->lat)/2.0;
				tr[tr_i].set_point(new_geo, tr[tr_i].cost, tr[tr_i].c0, tr[tr_i].c1, tr[tr_i].c2, tr[tr_i].c3);
			}
			// extrapolate from two end points
			else if (tr_i > 2)
			{
				new_geo.lon = 2*tr[tr_i-1].get_point()->lon - tr[tr_i-2].get_point()->lon;
				new_geo.lat = 2*tr[tr_i-1].get_point()->lat - tr[tr_i-2].get_point()->lat;
				tr[tr_i].set_point(new_geo, tr[tr_i].cost, tr[tr_i].c0, tr[tr_i].c1, tr[tr_i].c2, tr[tr_i].c3);
			}
		}
	}
}

/*****************************************************************************/

int sign(FP_TYPE a)
{
	// return the sign of a number: -1 = negative, 1 = positive (and zero)
	if (a>=0)
		return 1;
	else
		return -1;
}

/*****************************************************************************/

void track::interpolate(extrema_list* ex_list, int n_spl_pts)
{
	// new track after spline interpolation
	std::vector<track_point> new_tr;

	// vectors to store values to pass to the spline
	std::vector<FP_TYPE> t_V;
	std::vector<FP_TYPE> lon_V;
	std::vector<FP_TYPE> lat_V;
	std::vector<FP_TYPE> ints_V;
	std::vector<FP_TYPE> delta_V;

	// build the vectors containing the time, longitude, latitudes, values and deltas
	for (unsigned int i=0; i<tr.size(); i++)
	{
		t_V.push_back(i);
		lon_V.push_back(tr[i].get_point()->lon);
		lat_V.push_back(tr[i].get_point()->lat);
		ints_V.push_back(tr[i].get_point()->intensity);
		delta_V.push_back(tr[i].get_point()->delta);
	}

	// check for the longitude crossing the date line
	for (unsigned int i=1; i<lon_V.size(); i++)
	{
		if (lon_V[i] - lon_V[i-1] < -180)
			lon_V[i] += 360.0;
		if (lon_V[i] - lon_V[i-1] > 180)
			lon_V[i] -= 360.0;
	}

	// construct the splines
	spline lon_SPL(lon_V, t_V, -2e20f);
	spline lat_SPL(lat_V, t_V, -2e20f);
	spline ints_SPL(ints_V, t_V, -2e20f);
	spline delta_SPL(delta_V, t_V, -2e20f);

	// now interpolate the new track
	FP_TYPE t_start = tr[0].get_frame_number();
	for (unsigned int i=0; i<tr.size()-1; i++)
	{
		for (int j=0; j<n_spl_pts; j++)
		{
			// calculate the value to evaluate the splines at
			FP_TYPE v = i + FP_TYPE(j) * 1.0/n_spl_pts;
			// evaluate the splines
			FP_TYPE lon = lon_SPL.evaluate(v);
			// need these for crossing the date line
			if (lon > 360.0)
				lon -= 360.0;
			if (lon < 0.0)
				lon += 360.0;
			FP_TYPE lat = lat_SPL.evaluate(v);
			FP_TYPE ints = ints_SPL.evaluate(v);
			FP_TYPE delta = delta_SPL.evaluate(v);
			// calculate the timestep
			FP_TYPE t = t_start + v;
			// create the point
			track_point tp(t);
			// create the extremum
			steering_extremum svex(lon, lat, ints, delta, 
							 	 tr[i].get_point()->sv_u, 
							 	 tr[i].get_point()->sv_v);
			// set the extremum and add to the track
			tp.set_point(svex, tr[i].cost, tr[i].c0, tr[i].c1, tr[i].c2, tr[i].c3);
			new_tr.push_back(tp);
		}
	}
	// add last point and reassign
	new_tr.push_back(*(get_last_track_point()));
	tr = new_tr;
}

/*****************************************************************************/

int track::get_persistence(void)
{
	return tr.size();
}

/*****************************************************************************/

const int track::size(void) const
{
	return tr.size();
}

/*****************************************************************************/

track_point* track::get_last_track_point(void)
{
	return &(tr[tr.size()-1]);
}

/*****************************************************************************/

track_point* track::get_track_point(int tstep)
{
	int first_time_step = tr[0].get_frame_number();
	return &(tr[tstep-first_time_step]);
}

/*****************************************************************************/

track_point* track::get_track_point_idx(int idx)
{
	return &(tr[idx]);
}

/*****************************************************************************/

FP_TYPE track::get_cand_pt_cost(void)
{
	return cand_cost;
}

/*****************************************************************************/

steering_extremum* track::get_cand_pt(void)
{
	return &cand_pt;
}

/*****************************************************************************/

void track::set_cand_pt(steering_extremum svex_pt, FP_TYPE cost)
{
	cand_pt = svex_pt;
	cand_cost = cost;
}

/*****************************************************************************/

void track::set_cand_pt(steering_extremum svex_pt, FP_TYPE cost,
						FP_TYPE i_cc0, FP_TYPE i_cc1, FP_TYPE i_cc2, FP_TYPE i_cc3)
{
	cand_pt = svex_pt;
	cand_cost = cost;
	cc0 = i_cc0;
	cc1 = i_cc1;
	cc2 = i_cc2;
	cc3 = i_cc3;
}

/*****************************************************************************/

void track::consolidate_candidate_point(int t)
{
	if (cand_cost < 2e20f)	// don't add if no candidate points
	{
		track_point tp(t);
		tp.set_point(cand_pt, cand_cost, cc0, cc1, cc2, cc3);
		add_point(tp);
	}
	cand_cost = 2e20f;	// reset the candidate cost ready for the next frame
}

/*****************************************************************************/

track_list::track_list(void){}

/*****************************************************************************/

void track_list::set_size(int n_tracks)
{
    tr_list.resize(n_tracks);
}

/*****************************************************************************/

const int track_list::size(void) const
{
    return tr_list.size();
}

/*****************************************************************************/

track* track_list::get(int n_track)
{
	return &(tr_list[n_track]);
}

/*****************************************************************************/

const track* track_list::get(int n_track) const
{
	return &(tr_list[n_track]);
}

/*****************************************************************************/

void track_list::add(track trk)
{
	tr_list.push_back(trk);
}

/*****************************************************************************/

void track_list::consolidate_candidate_points(int t)
{
	for (unsigned int tr=0; tr<tr_list.size(); tr++)
		tr_list[tr].consolidate_candidate_point(t);
}

/*****************************************************************************/

META_DATA_TYPE* track_list::get_meta_data(void)
{
	return &meta_data;
}

/*****************************************************************************/

void track_list::set_meta_data(META_DATA_TYPE* in_meta_data)
{
	// copying like this rather than using the copy constructor allows multiple
	// sources of metadata
	for (META_DATA_TYPE::iterator it_md = in_meta_data->begin();
		 it_md != in_meta_data->end(); it_md++)
		meta_data[it_md->first] = it_md->second;
}

/*****************************************************************************/

void track_list::save(std::string output_fname)
{
	if (size() == 0)
		throw(std::string("# No tracks found!"));
	// open a binary file
	std::ofstream out;
	std::cout << "# Saving track data" << std::endl;
	out.open(output_fname.c_str(), std::ios::out | std::ios::binary);	
	if (!out)
		throw(std::string("Saving track data.  File could not be opened or written to: " + output_fname));
	
	// if meta data then write
	if (meta_data.size() != 0)
		write_meta_data(out, meta_data);
	// write number of tracks
	write_int(out, tr_list.size());
	
	// loop through all the tracks
	for (TR_LIST_TYPE::iterator it=tr_list.begin(); it!=tr_list.end(); it++)
	{
		// write the number of track points
		write_int(out, it->size());
		for (int tp=0; tp < it->size(); tp++)
		{
			// get the track point
			track_point* trk_pt = it->get_track_point_idx(tp);
			// get the steering_extremum			
			steering_extremum* svex = trk_pt->get_point();
			// write the frame number first
			write_float(out, trk_pt->frame_number);
			// write the extremum point
			svex->save(out);
			// write the costs
			write_float(out, trk_pt->cost);	// cost of adding to the track
			write_float(out, trk_pt->c0);
			write_float(out, trk_pt->c1);
			write_float(out, trk_pt->c2);
			write_float(out, trk_pt->c3);	// component costs			
		}
	}
	out.close();
}

/*****************************************************************************/

void track_list::save_text(std::string output_fname)
{
	if (size() == 0)
		throw(std::string("# No tracks found!"));
	// open a text file
	std::ofstream out;
	out.open(output_fname.c_str(), std::ios::out);
	if (!out)
		throw(std::string("Saving track data.  File could not be opened or written to: " + output_fname));
	// write number of tracks
	out << tr_list.size() << std::endl;
	
	// loop through all the tracks
	for (int i=0; i<tr_list.size(); i++)
	{
		// write the number of track points
		out << tr_list[i].size() << std::endl;
		for (int tp=0; tp < tr_list[i].size(); tp++)
		{
			// get the track point
			track_point* trk_pt = tr_list[i].get_track_point_idx(tp);
			// get the steering_extremum			
			steering_extremum* svex = trk_pt->get_point();
			// write the frame number first
			out << trk_pt->frame_number << " ";
			// write the extremum point
			svex->save_text(out);
			// write the costs
			out << trk_pt->cost << " " 
				<< trk_pt->c0 << " "
				<< trk_pt->c1 << " "
				<< trk_pt->c2 << " "
				<< trk_pt->c3 << std::endl;	// component costs			
		}
	}
	out.close();		
}

/*****************************************************************************/

void track_list::load(std::string input_fname)
{
	// open a binary file
	std::ifstream in_file;
	std::cout << "# Loading track data" << std::endl;
	in_file.open(input_fname.c_str(), std::ios::in | std::ios::binary);
	if (!in_file)
		throw(std::string("Loading track data.  File could not be opened: " + input_fname));
	
	// load any metadata
	meta_data = read_meta_data(in_file);
	// read number of tracks
	int n_tracks = read_int(in_file);
	
	// loop through all the tracks
	for (int i=0; i<n_tracks; i++)
	{
		// create track
		track trk;
		// read the number of track points
		int n_tp = read_int(in_file);
		for (int tp=0; tp < n_tp; tp++)
		{
			// create the track point
			track_point trk_pt;
			// read the frame number
			trk_pt.frame_number = read_float(in_file);
			// get the steering_extremum			
			steering_extremum* svex = trk_pt.get_point();
			// read the extremum point
			svex->load(in_file);
			// read in the costs
			trk_pt.cost = read_float(in_file);	// cost of adding to the track
			trk_pt.c0 = read_float(in_file);
			trk_pt.c1 = read_float(in_file);
			trk_pt.c2 = read_float(in_file);
			trk_pt.c3 = read_float(in_file);	// component costs
			// add the point to the track
			trk.add_point(trk_pt);
		}
		// add the track to the track list
		add(trk);
	}
	in_file.close();
}