/******************************************************************************
** Program : minima_locator.cpp
** Author  : Neil Massey
** Date    : 10/07/13
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded onto a regional hierarchical triangular mesh
******************************************************************************/

#include "minima_locator.h"
#include <sstream>
#include <math.h>

/*****************************************************************************/

minima_locator::minima_locator(void)
			   :extrema_locator()
{
}

/*****************************************************************************/
				   
minima_locator::~minima_locator(void)
{
}

/*****************************************************************************/

void minima_locator::parse_arg_string(std::string method_string)
{
	// args for minima locator are:
	// 1 = threshold
	// 2 = min delta
	// 3 = min number of surrounding triangles that are greater than current triangle
	// 4 = tolerance for region growing routine
	
	// get the first bracket
	int c_pos = method_string.find_first_of("(")+1;
	int e_pos = method_string.find(")", c_pos);
	char dummy;
	// check whether the text is "help" and print the method parameters if it is
	if (method_string.substr(c_pos, e_pos-c_pos) == "help")
		throw(std::string("minima parameters = (threshold, minimum delta, minimum number of surrounding triangles, tolerance in object finding)"));
	std::stringstream stream(method_string.substr(c_pos, e_pos-c_pos));
	stream >> thresh >> dummy >> min_delta >> dummy >> min_sur_tri >> dummy >> TOL;	// dummy reads the comma
	// add to the meta data
	std::stringstream ss;
	meta_data["method"] = "minima";
	ss << thresh;
	meta_data["threshold"] = ss.str();
	ss.str(""); ss << min_delta;
	meta_data["minimum_delta"] = ss.str();
	ss.str(""); ss << min_sur_tri;
	meta_data["minimum_number_of_surrounding_triangles"] = ss.str();
	ss.str(""); ss << TOL;
	meta_data["tolerance_in_gradient_ascent"] = ss.str();
}

/*****************************************************************************/

bool minima_locator::is_extrema(indexed_force_tri_3D* tri, int t_step)
{
	// args in base class are encoded as follows:
	// thresh = threshold
	// min_delta = min delta
	// min_sur_tri = min number of surrounding triangles that are greater than current triangle
	// TOL = tolerance for region growing routine
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

	// now find the minima
	
	FP_TYPE tri_val = ds.get_data(t_step, tri->get_ds_index());
	// if it's the missing value then return false
	if (tri_val == ds.get_missing_value())
		return false;
	// if it's higher than the threshold then return false
	if (tri_val > thresh)
		return false;
	// number of surrounding triangles that are greater than current triangle
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
		// if the diff is -ve then the tri_val is less than the adj_val
		n_st += (diff <= 0.0) ? 1 : 0;
	}
	return (n_st >= min_sur_tri);
}

/*****************************************************************************/

bool minima_locator::is_in_object(indexed_force_tri_3D* O_TRI, 
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
	is_in = is_in && cl_v <= ol_v;
	is_in = is_in && cl_v < thresh;
	return is_in;	
}

/*****************************************************************************/

bool minima_locator::process_data(void)
{
	return true;
}

/*****************************************************************************/

FP_TYPE minima_locator::calculate_point_weight(FP_TYPE V, FP_TYPE min_v, FP_TYPE max_v)
{
	// for minima the weight is 1.0 - (V-min_v) / (max_v-min_v)
	return 1.0 - (V-min_v) / (max_v-min_v);
}

/*****************************************************************************/

indexed_force_tri_3D* minima_locator::get_original_triangle(int o, int t)
{
	// get the original triangle - i.e. the one at the centre of the object
	// for the minima this is defined as the triangle with the lowest value
	FP_TYPE c_val = 2e20f;
	indexed_force_tri_3D* o_tri = NULL;
	LABEL_STORE* object_labels = &(ex_list.get(t, o)->object_labels);
	// if there are no labels do not try to find the min/max
	if (object_labels->size() == 0)
		return NULL;
	// get the missing value
	FP_TYPE mv = ds.get_missing_value();
	// loop over all the triangles in the object
	for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
		 it_ll != object_labels->end(); it_ll++)	// tri indices
	{
		// get the triangle
		indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
		FP_TYPE val = ds.get_data(t, c_tri->get_ds_index());
		// find the min and max values
		if (fabs(val) < 0.99*fabs(mv) && val < c_val)
		{
			c_val = val;
			o_tri = c_tri;
		}
	}
	return o_tri;
}