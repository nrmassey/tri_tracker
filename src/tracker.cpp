/******************************************************************************
** Program : tracker.cpp
** Author  : Neil Massey
** Date    : 05/08/09
** Purpose : class that connects extrema feature points into a track
******************************************************************************/

#include "tracker.h"
#include "haversine.h"
#include "set_cout_precision.h"
#include "vector_3D.h"
#include "geo_convert.h"
#include "get_bearing.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "meta_data.h"

/*****************************************************************************/

tracker::tracker(std::vector<std::string> iinput_fname, int imin_per, 
				 FP_TYPE imin_len, FP_TYPE imin_dev, FP_TYPE isr)
		: min_per(imin_per), min_len(imin_len), min_dev(imin_dev), sr(isr), 
		  w0(0.25), w1(0.25), w2(0.625), w3(0.75), ts(1), n_spl(0)
{
	input_fname = iinput_fname;
	for (unsigned int i=0; i<input_fname.size(); i++)
	{
	    std::ifstream file_in(input_fname[i].c_str(), std::ios::in);
		std::cout << "# Reading data" << std::endl;
		if (!file_in.is_open())
			throw ("File " + input_fname[i] + " could not be opened.");
		ex_list.load(input_fname[i], mv, i>0);
		file_in.close();
	}

//	MAX_COST = (w0 + w1 + w2 + w3);
	MAX_COST = 1.0;
//	MAX_COST = 1e3;

	// create and set the meta_data
	std::stringstream ss;
	META_DATA_TYPE meta_data;
	for (unsigned int i=0; i<input_fname.size(); i++)
	{
		ss << "input_file_name_" << i;
		meta_data[ss.str()] = input_fname[i];
	}
	ss << min_per;
	meta_data["minimum_persistence"] = ss.str();
	ss.str(""); ss << min_len;
	meta_data["minimum_length"] = ss.str();
	ss.str(""); ss << min_dev;
	meta_data["minimum_deviation"] = ss.str();
	ss.str(""); ss << sr;
	meta_data["search_radius"] = ss.str();
	ss.str(""); ss << w0;
	meta_data["weight_0_(proximity)"] = ss.str();
	ss.str(""); ss << w1;
	meta_data["weight_1_(intensity)"] = ss.str();
	ss.str(""); ss << w2;
	meta_data["weight_2_(steering_vector)"] = ss.str();
	ss.str(""); ss << w3;
	meta_data["weight_3_(curvature)"] = ss.str();
	ss.str(""); ss << ts;
	meta_data["timesteps_allowed_between_feature_points"] = ss.str();
	ss.str(""); ss << n_spl;
	meta_data["spline_interpolation"] = ss.str();
	ss.str(""); ss << ksteps;				// skip between number of timesteps
	meta_data["timesteps_skipped_between_feature_points"] = ss.str();
	ss.str(""); ss << MAX_COST;
	meta_data["maximum_allowed_in_cost_function"] = ss.str();
	tr_list.set_meta_data(&meta_data);
}

/*****************************************************************************/

tracker::~tracker(void)
{
}

/*****************************************************************************/

void tracker::save(std::string output_fname)
{
	// open the output file
	tr_list.save(output_fname);
}

/*****************************************************************************/

void tracker::save_text(std::string output_fname)
{
	// open the output file
	tr_list.save_text(output_fname);
}

/*****************************************************************************/

void tracker::set_weights(FP_TYPE iw0, FP_TYPE iw1, FP_TYPE iw2, FP_TYPE iw3)
{
	w0 = iw0;
	w1 = iw1;
	w2 = iw2;
	w3 = iw3;
	std::cout << "# Using weights:" << std::endl;
	std::cout << "# " << w0 << " " << w1 << " " << w2 << " " << w3 << std::endl;
}

/*****************************************************************************/

void tracker::set_tsteps(int t)
{
	ts = t;
}

/*****************************************************************************/

void tracker::set_n_spl_pts(int s)
{
	n_spl = s;
}

/*****************************************************************************/

void tracker::set_ksteps(int k)
{
	ksteps = k;
}

/*****************************************************************************/

FP_TYPE dist_cost(track* TR, steering_extremum* EX_svex, FP_TYPE sr)
{
	// calculate cost incurred by the distance from the track to the candidate
	// point
	steering_extremum* TR_svex = TR->get_last_track_point()->get_point();	
	FP_TYPE d = haversine(TR_svex->lon, TR_svex->lat, EX_svex->lon, EX_svex->lat, EARTH_R);
	// check distance - normalise by search radius
	FP_TYPE c = d / sr;				// normalize by search radius
	return c;
}

/*****************************************************************************/

FP_TYPE intensity_cost(track* TR, steering_extremum* EX_svex)
{
	steering_extremum* TR_prev = TR->get_last_track_point()->get_point();
	FP_TYPE c = fabs((TR_prev->intensity - EX_svex->intensity) / TR_prev->intensity);
	return c;
}

/*****************************************************************************/

FP_TYPE geowind_cost(track* TR, steering_extremum* EX_svex)
{
	// Calculate the cost according to the steering term - in this case geostrophic 
	// wind convert coordinates to 3D, add geo wind vector to track point vector, 
	// project back to lat/lon and take the bearing

	// convert the lon / lat of the last point to a 3D vector
	steering_extremum* TR_lst_pt = TR->get_last_track_point()->get_point();
	vector_3D TR_lst_pt_3D = model_to_cart(TR_lst_pt->lon, TR_lst_pt->lat);
	// convert the candidate point
	vector_3D EX_cnd_pt_3D = model_to_cart(EX_svex->lon, EX_svex->lat);
	// get the magnitude of the vector between TR_lst_pt and EX_cnd_pt
	FP_TYPE cnd_mag = (EX_cnd_pt_3D - TR_lst_pt_3D).mag();
	// get the magnitude of the geostrophic wind vector
	vector_3D geo_pt_3D = vector_3D(TR_lst_pt->sv_u, TR_lst_pt->sv_v, 0.0);
	FP_TYPE geo_mag = geo_pt_3D.mag();
	// now create the vector of the last pt plus the geo-strophic wind
	vector_3D proj_pt_3D = TR_lst_pt_3D + geo_pt_3D * (cnd_mag / geo_mag);
	// convert back to latitude / longitude
	FP_TYPE proj_lat, proj_lon;
	proj_pt_3D *= 1.0 / proj_pt_3D.mag();
	cart_to_model(proj_pt_3D, proj_lon, proj_lat);
	// measure the bearing between the projected point and the original pt
	// and subtract from the bearing between the original pt and the candidate pt
	FP_TYPE a = get_bearing(TR_lst_pt->lon, TR_lst_pt->lat, proj_lon, proj_lat);
	FP_TYPE b = get_bearing(TR_lst_pt->lon, TR_lst_pt->lat, EX_svex->lon, EX_svex->lat);
	FP_TYPE c = fabs(b-a);
	return c;
}

/*****************************************************************************/

FP_TYPE curvature_cost(track* TR, steering_extremum* EX_svex, int t)
{
	// Calculate the curvature as the change in bearing between (pt1 -> pt2)
	// and (pt2 -> candidate point)
	int pr = TR->get_persistence();	// last point in the track
	steering_extremum* TR_pt_1 = TR->get_track_point_idx(pr-2)->get_point();
	steering_extremum* TR_pt_2 = TR->get_track_point_idx(pr-1)->get_point();
	FP_TYPE c = 0.0;

	// same point incurs zero cost, rather than undefined
	if (TR_pt_2->lon == EX_svex->lon && TR_pt_2->lat == EX_svex->lat)
		c = 0.0;
	else
	{
		// if the two last points are the same point then the bearing difference 
		// will give a peculiar answer - so use a previous point
		int i=3;
		while (pr-i >= 0 && (TR_pt_1->lon == TR_pt_2->lon && TR_pt_1->lat == TR_pt_2->lat))
		{
			TR_pt_1 = TR->get_track_point_idx(pr-i)->get_point();
			i++;			
		}
		if (pr-i >= 0)
		{
			// check whether they are the same point first - this has zero cost,
			// rather than the undefined cost that the bearing difference will report
			FP_TYPE a = get_bearing(TR_pt_1->lon, TR_pt_1->lat, TR_pt_2->lon, TR_pt_2->lat);
			FP_TYPE b = get_bearing(TR_pt_2->lon, TR_pt_2->lat, EX_svex->lon, EX_svex->lat);
			c = fabs(b - a);
		}
		else
			c = 0.0;
	}
	return c;
}

/*****************************************************************************/

FP_TYPE tracker::cmpt_cost_fn(track* TR, steering_extremum* EX_svex, int t,
	    					  FP_TYPE& c0, FP_TYPE& c1, FP_TYPE& c2, FP_TYPE& c3)
{
	// calculate and return the components of the cost function
	FP_TYPE c = 0;
	c0 = w0 * dist_cost(TR, EX_svex, sr);
	c = c0;
	if (TR->get_persistence() >= 1)
	{
		c1 = w1 * intensity_cost(TR, EX_svex);
		c += c1;
	}
	else
		c1 = 0.0;
	if (TR->get_persistence() >= 1 && EX_svex->sv_u != mv)
	{
		c2 = w2 * geowind_cost(TR, EX_svex);
		c += c2;
	}
	else
		c2 = 0.0;
	if (TR->get_persistence() >= 2)
	{
		c3 = w3 * curvature_cost(TR, EX_svex, t);
		c += c3;
	}
	else
		c3 = 0.0;
	return c;
}

/*****************************************************************************/

int tracker::determine_track_for_point(steering_extremum* svex, int t)
{
	// set up the minimum cost so far
	int min_tr = -1;
	FP_TYPE min_c = 2e20f;
	FP_TYPE c0, c1, c2, c3;
	// loop through each track
	for (int tr=0; tr<tr_list.size(); tr++)
	{
		// check that the track occurs within the permitted number of
		// timesteps
		int lt = tr_list.get(tr)->get_last_track_point()->get_frame_number();
		if (t-lt > ts+ksteps)
			continue; // next track!
		// calculate the cost of adding the extremum to this track
		FP_TYPE c = cmpt_cost_fn(tr_list.get(tr), svex, t, c0, c1, c2, c3);
		// check if this is less than the current candidate cost
		if (c < min_c && c < tr_list.get(tr)->get_cand_pt_cost() &&
			c < MAX_COST)
		{
			min_c = c;
			min_tr = tr;
		}
	}
	return min_tr;
}

/*****************************************************************************/

void tracker::assign_point_to_track(steering_extremum* c_svex, int min_tr, int t)
{
	// now the minimum track location has been found - add to the track
	// check first whether a track was found
	if (min_tr == -1)
	{
		// not found so add to unassigned stack
		ua_ex_stack.push(*c_svex);
	}
	else
	{
		// otherwise - check whether an extremum has already been 
		// assigned to the track
		// if it has then put the old extremum back into the stack
		if (tr_list.get(min_tr)->get_cand_pt_cost() < 2e20f)
		{
			steering_extremum* geo_old_ex = tr_list.get(min_tr)->get_cand_pt();					
			ex_stack.push(*geo_old_ex);
		}
		// get the cost - have to get each component again
		FP_TYPE c0, c1, c2, c3;
		FP_TYPE c = cmpt_cost_fn(tr_list.get(min_tr), c_svex, t, c0, c1, c2, c3);				
		// set as candidate point
		tr_list.get(min_tr)->set_cand_pt(*c_svex, c, c0, c1, c2, c3);
	}
}

/*****************************************************************************/

void tracker::clear_stacks(void)
{
	// clear the two stacks used
	while (!ex_stack.empty())
		ex_stack.pop();
	while (!ua_ex_stack.empty())
		ua_ex_stack.pop();
}

/*****************************************************************************/

void tracker::add_unassigned_points_as_tracks(int t)
{
	while (!ua_ex_stack.empty())
	{
		// get the unassigned extremum and pop it from the stack
		steering_extremum ua_svex = ua_ex_stack.top();
		ua_ex_stack.pop();
		// create a new track point
		track trk;
		track_point tp(t);
		tp.set_point(ua_svex, 0.0, 0.0, 0.0, 0.0, 0.0);
		trk.add_point(tp);
		tr_list.add(trk);
	}
}

/*****************************************************************************/

void tracker::find_tracks(void)
{
	std::cout << "# Locating tracks, timestep: ";
	// create the tracks from the first frame's worth of points
	build_first_frame();

	// loop through all the timesteps
	for (int t=ksteps; t<ex_list.size(); t+=ksteps)
	{
		std::cout << t;
        std::cout.flush();
		// build a stack of the extrema
		clear_stacks();
		for (int e=0; e<ex_list.number_of_extrema(t); e++)
		{
			steering_extremum* svex = ex_list.get(t,e);
			ex_stack.push(*svex);
		}
		// loop until no extrema left
		while (!ex_stack.empty())
		{
			// pop off an extremum to compare against 
			steering_extremum c_svex = ex_stack.top();
			ex_stack.pop();
			// determine which track it should be assigned to and assign it
			int min_tr = determine_track_for_point(&c_svex, t);
			assign_point_to_track(&c_svex, min_tr, t);
		}
		// consolidate the candidate points - i.e. add the cand pts to the end
		// of the track
		tr_list.consolidate_candidate_points(t);
		
		// create new tracks for the remaining (unassigned) points
		add_unassigned_points_as_tracks(t);
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
	// reduce the tracks to timesteps that have had movement within them
	reduce_tracks();
	// cull any tracks that are too short in length and time
	trim_tracks();
	// interpolate the tracks along a spline for the number of points passed in
	interpolate_tracks();
}

/*****************************************************************************/

void tracker::build_first_frame(void)
{
	tr_list.set_size(ex_list.number_of_extrema(0));
	// create a track for each extrema that exists in timestep 0
	for (int i=0; i<ex_list.number_of_extrema(0); i++)
	{
		track_point tp(0);
		// get the extremum
		steering_extremum svex(*(ex_list.get(0, i)));
		// add to the list
		tp.set_point(svex, 0.0, 0.0, 0.0, 0.0, 0.0);
		tr_list.get(i)->add_point(tp);
	}
}

/*****************************************************************************/

void tracker::reduce_tracks(void)
{
	// remove duplicate lon / lat values from the tracks - in prep for
	// interpolation later
	std::cout << "# Reducing tracks" << std::endl;
	for (int tr_i=0; tr_i<tr_list.size(); tr_i++)
		tr_list.get(tr_i)->reduce();

	// split tracks that have a long gap between timesteps
	track_list new_tr_list;
	FP_TYPE const gts = 4;	// should supply this as a parameter as it's time dependent
	for (int tr_list_i=0; tr_list_i<tr_list.size(); tr_list_i++)
	{
		track new_tr;
		track* tr = tr_list.get(tr_list_i);
		if (tr->size() == 0)
			continue;
		for (int tr_i=0; tr_i<tr->size()-1; tr_i++)
		{
			track_point* tp0 = tr->get_track_point_idx(tr_i);
			track_point* tp1 = tr->get_track_point_idx(tr_i+1);
			if (tp1->get_frame_number() - tp0->get_frame_number() >= gts*ksteps)
			{
				new_tr_list.add(new_tr);
				new_tr.clear();
			}
			else
				new_tr.add_point(*tp0);
		}
		new_tr.add_point(*(tr->get_last_track_point()));
        new_tr_list.add(new_tr);
    }

	// remove any orphaned tracks
	track_list new_tr_list2;
	new_tr_list2.set_meta_data(tr_list.get_meta_data());
	for (int tr_i=0; tr_i<new_tr_list.size(); tr_i++)
	{
		if (new_tr_list.get(tr_i)->size() > 0)
			new_tr_list2.add(*(new_tr_list.get(tr_i)));
	}
	tr_list = new_tr_list2;
}

/*****************************************************************************/

void tracker::trim_tracks(void)
{
	std::cout << "# Trimming tracks" << std::endl;
	// remove any tracks that don't meet the minimum persistence or length
	// criteria
	track_list new_tr_list;
	new_tr_list.set_meta_data(tr_list.get_meta_data());

	for (int tr_i=0; tr_i<tr_list.size(); tr_i++)
	{
		track* tr = tr_list.get(tr_i);
		// calculate the latitude scaling for the curvature
		bool add = true;
		add = add && tr->get_persistence() >= min_per;
		add = add && tr->get_length() >= min_len;
		add = add && tr->get_deviation() >= min_dev;

		if (add)
			new_tr_list.add(*tr);
	}
	tr_list = new_tr_list;
}

/*****************************************************************************/

void tracker::interpolate_tracks(void)
{
	if (n_spl != 0)
	{
		std::cout << "# Interpolating tracks" << std::endl;
		for (int tr_i=0; tr_i<tr_list.size(); tr_i++)
			tr_list.get(tr_i)->interpolate(&ex_list, n_spl);
	}
}
