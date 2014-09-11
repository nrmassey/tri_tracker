/******************************************************************************
** Program : maxima_locator.cpp
** Author  : Neil Massey
** Date    : 10/07/13
** Purpose : class inherited from extrema_locator that searches for maxima 
**           in data regridded onto a regional hierarchical triangular mesh
******************************************************************************/

#include "maxima_locator.h"
#include <sstream>
#include <math.h>
#include "geo_convert.h"

/*****************************************************************************/

maxima_locator::maxima_locator(void)
			   :extrema_locator()
{
}

/*****************************************************************************/
				   
maxima_locator::~maxima_locator(void)
{
}

/*****************************************************************************/

void maxima_locator::parse_arg_string(std::string method_string)
{
	// args for minima locator are:
	// 1 = threshold
	// 2 = min delta
	// 3 = min number of surrounding triangles that are greater than current triangle
	
	// get the first bracket
	int c_pos = method_string.find_first_of("(")+1;
	int e_pos = method_string.find(")", c_pos);
	char dummy;
	// check whether the text is "help" and print the method parameters if it is
	if (method_string.substr(c_pos, e_pos-c_pos) == "help")
		throw(std::string("minima parameters = (threshold, minimum delta, minimum number of surrounding triangles)"));
	std::stringstream stream(method_string.substr(c_pos, e_pos-c_pos));
	stream >> thresh >> dummy >> max_delta >> dummy >> min_sur_tri;	// dummy reads the comma
	// add to the meta data
	std::stringstream ss;
	meta_data["method"] = "minima";
	ss << thresh;
	meta_data["threshold"] = ss.str();
	ss.str(""); ss << max_delta;
	meta_data["minimum_delta"] = ss.str();
	ss.str(""); ss << min_sur_tri;
	meta_data["minimum_number_of_surrounding_triangles"] = ss.str();
}

/*****************************************************************************/

bool maxima_locator::is_extrema(indexed_force_tri_3D* tri, int t_step)
{
	// args in base class are encoded as follows:
	// thresh = threshold
	// max_delta = max delta
	// min_sur_tri = min number of surrounding triangles that are greater than current triangle
	// get the value of the triangle
	
	// build a list of surrounding triangles first
	const LABEL_STORE* adj_tri_labels = tri->get_adjacent_labels(adj_type);
	// create the tri list and add to it
	std::list<indexed_force_tri_3D*> adj_tris;
	for (LABEL_STORE::const_iterator adj_it = adj_tri_labels->begin();
		 adj_it != adj_tri_labels->end(); adj_it++)
	{
		indexed_force_tri_3D* a_tri = tg.get_triangle(*adj_it);
		adj_tris.push_back(a_tri);
	}

	// now find the maxima	
	FP_TYPE tri_val = ds.get_data(t_step, tri->get_ds_index());
	// if it's the missing value then return false
	if (tri_val == ds.get_missing_value())
		return false;
	// if it's lower than the threshold then return false
	if (tri_val < thresh)
		return false;
	// number of surrounding triangles that are less than current triangle
	int n_st = 0;
	// loop through all the adjacent triangles
	for (std::list<indexed_force_tri_3D*>::iterator it = adj_tris.begin();
		 it != adj_tris.end(); it++)
	{
		// get the value of the adjacent triangle
		FP_TYPE adj_val = ds.get_data(t_step, (*it)->get_ds_index());
		// if it's the missing value then continue onto next one
		if (adj_val == ds.get_missing_value())
			continue;
		// difference between middle triangle and this surrounding triangle
		FP_TYPE diff = tri_val - adj_val;
		// if the diff is +ve then the tri_val is greater than the adj_val
		n_st += (diff >= 0.0) ? 1 : 0;
	}
	return (n_st >= min_sur_tri);
}

/*****************************************************************************/

bool maxima_locator::is_in_object(indexed_force_tri_3D* O_TRI, 
							  	  indexed_force_tri_3D* C_TRI, int t_step)
{
	// O_TRI - original label
	// T_TRI - current triangle label
	// C_TRI - candidate label - triangle being tested for inclusion
	bool is_in = true;
	// get the candidate triangle value
	FP_TYPE cl_v = ds.get_data(t_step, C_TRI->get_ds_index());
	FP_TYPE ol_v = ds.get_data(t_step, O_TRI->get_ds_index());
	// quick checks
	is_in = is_in && cl_v != ds.get_missing_value();
	is_in = is_in && cl_v >= ol_v;
	is_in = is_in && cl_v > thresh;
	return is_in;	
}

/*****************************************************************************/

void maxima_locator::calculate_object_position(int o, int t)
{
	steering_extremum* svex = ex_list.get(t, o);
	LABEL_STORE* object_labels = &(svex->object_labels);
	// find min / max of values in the object
	FP_TYPE min_v = 2e20f;	// minimum value in object
	FP_TYPE max_v = 0;		// maximum value in object
	get_min_max_values(min_v, max_v, o, t);
	
	// position vector in Cartesian coordinates
	vector_3D P;
	// get the position and weight for each triangle in the object
	for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
		 it_ll != object_labels->end(); it_ll++)	// tri labels
	{
		// get the triangle and its centroid
		indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
		vector_3D C = c_tri->centroid();
		// get the data value from the datastore
		FP_TYPE V = ds.get_data(t, c_tri->get_ds_index());
		// if this value equals the min value then this is the object position
		if (V == max_v)
			P = C;
	}
	// put the values back in the geo extremum
	FP_TYPE lon, lat;
	cart_to_model(P, lon, lat);
	svex->lon = lon;
	svex->lat = lat;
}

/*****************************************************************************/

void maxima_locator::calculate_object_intensity(int o, int t)
{
	// get the extremum - the lat and lon will have been set already
	steering_extremum* svex = ex_list.get(t, o);
	LABEL_STORE* object_labels = &(svex->object_labels);
	
	// find min / max of values in the object
	FP_TYPE min_v = 2e20f;	// minimum value in object
	FP_TYPE max_v = 0;		// maximum value in object
	get_min_max_values(min_v, max_v, o, t);		
	svex->intensity = max_v;
}

/*****************************************************************************/

void maxima_locator::calculate_object_delta(int o, int t)
{
	// calculate the delta as the absolute difference between the minimum and
	// the maximum value of the object
	// find min / max of values in the object
	FP_TYPE min_v = 2e20f;	// minimum value in object
	FP_TYPE max_v = -2e20f;		// maximum value in object
	get_min_max_values(min_v, max_v, o, t);	
	steering_extremum* svex = ex_list.get(t, o);
	svex->delta = fabs(max_v - min_v);
}