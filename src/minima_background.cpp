/******************************************************************************
** Program : minima_background.cpp
** Author  : Neil Massey
** Date    : 07/08/13
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded, and removing the background field
******************************************************************************/

#include "minima_background.h"
#include "haversine.h"
#include "geo_convert.h"
#include <sstream>
#include <math.h>

/******************************************************************************/

inline bool is_mv(FP_TYPE V, FP_TYPE mv)
{
	return !(fabs(V) < fabs(mv*0.99));
}

/******************************************************************************/

FP_TYPE contour_data(FP_TYPE V, FP_TYPE C)
{
	return FP_TYPE(int(V/C) * C);
}

/******************************************************************************/

minima_background::minima_background(void) : extrema_locator(), bck_field_ds(NULL)
{
}

/******************************************************************************/

minima_background::~minima_background(void)
{
	delete data_minus_bck;
	delete bck_field_ds;
}

/******************************************************************************/

void minima_background::locate(void)
{
	process_data();
	find_extrema();
	refine_objects();
	find_objects();
	trim_objects();
	merge_objects();
	ex_points_from_objects();
}

/******************************************************************************/

void minima_background::parse_arg_string(std::string method_string)
{
	// arguments are:
	// arg[0] = file to read background
	// arg[1] = averaging period to take background field over
	// arg[2] = mesh level to detect large scale minima at
	// arg[3] = contour level
	// arg[4] = minimum delta
	// parameters for minima location with background removal
	// get the first bracket
	int c_pos = method_string.find_first_of("(")+1;
	int e_pos = method_string.find(")", c_pos);
	char dummy;
	if (method_string.substr(c_pos, e_pos-c_pos) == "help")
		throw(std::string("minima_back parameters = (file name to take background field from, averaging period of background field, contour value, min delta)"));
	
	int b_pos = method_string.find(",", c_pos);
	bck_field_file = method_string.substr(c_pos, b_pos-c_pos);
	std::stringstream stream(method_string.substr(b_pos+1, e_pos-b_pos));
	stream >> bck_avg_period >> dummy
		   >> min_mesh_lvl >> dummy
		   >> contour_value >> dummy
		   >> min_delta;
		   
	// add to the metadata
	std::stringstream ss;
	meta_data["method"] = "minima_back";
	meta_data["background_file"] = bck_field_file;
	ss.str(""); ss << bck_avg_period;
	meta_data["background_averaging_period"] = ss.str();
	ss.str(""); ss << min_mesh_lvl;
	meta_data["minima_mesh_level"] = ss.str();
	ss.str(""); ss << contour_value;
	meta_data["contour_value"] = ss.str();
	ss.str(""); ss << min_delta;
	meta_data["minimum_delta"] = ss.str();
}

/******************************************************************************/

void add_child_labels_to_object(steering_extremum* svex, QT_TRI_NODE* c_qt_tri,
								int max_lvl)
{
	if (c_qt_tri->get_level() == max_lvl)
		svex->object_labels.push_back(c_qt_tri->get_data()->get_label());
	else if (!c_qt_tri->is_leaf())
	{
		for (int i=0; i<4; i++)
			add_child_labels_to_object(svex, c_qt_tri->get_child(i), max_lvl);
	}
}

/******************************************************************************/

void minima_background::find_extrema(void)
{
	// overloaded this function as it does something quite different to the
	// standard function:
	// 1. minima are located at a level lower than that of the maximum grid level
	// 2. this extremum's children are also marked as minima
	// resize the extrema list to be the correct size for the number of timesteps
	std::cout << "# Locating extrema, timestep: ";
	ex_list.set_size(data_minus_bck->get_number_of_time_steps());
	// get a list of all the triangles at the required level
	std::list<QT_TRI_NODE*> tris = tg.get_triangles_at_level(min_mesh_lvl);
	// get the missing value
	FP_TYPE mv = data_minus_bck->get_missing_value();
	// repeat over all timesteps
	FP_TYPE sum_ex = 0.0;
	for (int t=0; t<data_minus_bck->get_number_of_time_steps(); t++)
	{
		tstep_out_begin(t);	
		// repeat over all the nodes in this level
		for (std::list<QT_TRI_NODE*>::iterator it = tris.begin();
			 it != tris.end(); it++)
		{
			indexed_force_tri_3D* c_tri = (*it)->get_data();
			if (is_extrema(c_tri, t))
			{
				// add to the extrema list, via the object list - the location and
				// value of the svex is currently just filled with the missing value
				steering_extremum svex(mv, mv, mv, mv, mv, mv);
				add_child_labels_to_object(&svex, (*it), grid_level);
				// now add to the extrema list
				ex_list.add(t, svex);
				sum_ex += 1.0;
			}
		}
		tstep_out_end(t);
	}
	std::cout << " Number of extrema: " << sum_ex << " ";
	std::cout << std::endl;
}

/******************************************************************************/

bool minima_background::is_extrema(indexed_force_tri_3D* tri, int t_step)
{
	// get the data for the triangle for this time step and contour it
	FP_TYPE tri_val = contour_data(data_minus_bck->get_data(t_step, tri->get_ds_index()),
								   contour_value);
	// if it's the missing value then return false
	FP_TYPE mv = contour_data(data_minus_bck->get_missing_value(), contour_value);
	if (fabs(tri_val) >= fabs(0.99*mv))
		return false;

	int n_st = 0;	
	// loop through all the adjacent triangles
	const LABEL_STORE* tri_adj_labels = tri->get_adjacent_labels(adj_type);
	// nearest neighbour search of all adjacent triangles for a minimum
	for (LABEL_STORE::const_iterator tri_adj_it = tri_adj_labels->begin();
		 tri_adj_it != tri_adj_labels->end(); tri_adj_it++)
	{
		// get the triangle from the label
		indexed_force_tri_3D* tri_adj = tg.get_triangle(*tri_adj_it);
		// get the value of the adjacent triangle		
		FP_TYPE tri_adj_val = contour_data(data_minus_bck->get_data(t_step, tri_adj->get_ds_index()),
		                                   contour_value);
		// if it's the missing value then continue onto next one - but count it as greater than current
		if (fabs(tri_adj_val) >= fabs(0.99*mv))
		{
			n_st += 1;
			continue;
		}
		// if the middle triangle is less than or equal to this surrounding triangle
		if (tri_val <= min_delta && tri_val <= tri_adj_val)
			n_st += 1;
	}
	int min_sur = tri_adj_labels->size();
	return (n_st >= min_sur);
}

/******************************************************************************/

bool minima_background::is_in_object(indexed_force_tri_3D* O_TRI, 
				  					 indexed_force_tri_3D* C_TRI, int t_step)
{
	// O_TRI - original triangle
	// C_TRI - candidate triangle - triangle being tested for inclusion
	bool is_in = false;

	// get the candidate triangle value
//	FP_TYPE cl_v = contour_data(data_minus_bck->get_data(t_step, C_TRI->get_ds_index()), contour_value);
	FP_TYPE cl_v = data_minus_bck->get_data(t_step, C_TRI->get_ds_index());
	FP_TYPE ol_v = contour_data(data_minus_bck->get_data(t_step, O_TRI->get_ds_index()), contour_value);

	// quick check	
	is_in = cl_v < (ol_v + contour_value);  // within 1 contour
	// not the mv
	is_in = is_in && fabs(cl_v) <= fabs(0.99 * data_minus_bck->get_missing_value());
	// less than the minimum delta
	is_in = is_in && (ol_v <= min_delta);
	// less than 1000km radius
	is_in = is_in && tg.distance_between_triangles(O_TRI->get_label(), C_TRI->get_label())/1000.0 < 1000.0;
	return is_in;
}

/******************************************************************************/

void minima_background::calculate_background_field(void)
{
	std::cout << "# Calculating background field" << std::endl;
	
	// load in the field first
	bck_field_ds = new data_store();
	bck_field_ds->load(bck_field_file);
	// get the size of the current datastore
	int n_ts = bck_field_ds->get_number_of_time_steps();
	int n_idx = bck_field_ds->get_number_of_indices();
	
	// do we need to take a mean?
	if (bck_avg_period > 1 && n_ts >= bck_avg_period)
	{
		// create a new background field
		data_store* new_bck_field_ds = new data_store();
		// set the scaling for each averaging period
		FP_TYPE scale = 1.0 / bck_avg_period;
		// create the datastore
		new_bck_field_ds->set_size(n_ts/bck_avg_period, n_idx);
		new_bck_field_ds->set_missing_value(bck_field_ds->get_missing_value());
		// loop through the data producing the average
		for (int t=0; t<n_ts; t++)
		{
			int dest_pos = t / bck_avg_period;
			for (int i=0; i<n_idx;  i++)
			{
				// get the current data value, add the value from the original ds store
				// multiplied by the scaler
				FP_TYPE c_val = new_bck_field_ds->get_data(dest_pos, i);
				FP_TYPE t_val = scale * bck_field_ds->get_data(t, i) + c_val;
				new_bck_field_ds->set_data(dest_pos, i, t_val);
			}
		}
		// assign the bck_field to be the new_field and delete the old one
		data_store* old_bck_field_ds = bck_field_ds;
		bck_field_ds = new_bck_field_ds;
		delete old_bck_field_ds;
	}
}

/******************************************************************************/

bool minima_background::process_data(void)
{
	// process the data
	// this removes the background field from the data
	// the background field may be averaged over a time period - i.e. remove
	// the monthly average / 5 day average or just daily average
	// the difference is taken between the triangle in the data and the corresponding
	// triangle in the background field.  This is repeated for all levels so that data
	// is subtracted from the background field at a number of resolutions
	
	// first calculate the background field
	calculate_background_field();

	// get details about the input data
	int n_ts = ds.get_number_of_time_steps();
	FP_TYPE mv = ds.get_missing_value();	
	int bck_field_nts = bck_field_ds->get_number_of_time_steps();

	// create the storage for the result
	data_minus_bck = new data_store();	
	data_minus_bck->set_size(n_ts, ds.get_number_of_indices());
	data_minus_bck->set_missing_value(mv);
	
	std::cout << "# Processing data" << std::endl;

	// we want to process every level, not just the grid level
	for (int l=0; l<tg.get_max_level(); l++)
	{
		// get the triangles for this level
		std::list<QT_TRI_NODE*> tris_qn = tg.get_triangles_at_level(l);
		for (std::list<QT_TRI_NODE*>::iterator it_qt = tris_qn.begin();
			 it_qt != tris_qn.end(); it_qt++)
		{
			// get the index into the tri-grid
			int ds_idx = (*it_qt)->get_data()->get_ds_index();
			// use this index over every timestep to subtract the background field
			for (int t=0; t<n_ts; t++)
			{
				FP_TYPE dV = ds.get_data(t, ds_idx);
				// check for stepping outside of the array
				int bck_t = t/bck_avg_period;
				if (bck_t >= bck_field_nts)
					bck_t = bck_field_nts-1;
				FP_TYPE bV = bck_field_ds->get_data(bck_t, ds_idx);
				FP_TYPE V = mv;
				if (!(is_mv(dV, mv) || is_mv(bV, mv)))
					V =  dV - bV;
				data_minus_bck->set_data(t, ds_idx, V);
			}
		}
	}
	
	// option to save the output - build the filename first
	std::string out_fname = ds_fname.substr(0, ds_fname.size()-4)+"_bck.rgd";	
	data_minus_bck->save(out_fname);
	return true;
}

/******************************************************************************/

void minima_background::trim_objects(void)
{
	std::cout << "# Trimming objects, timestep: ";
	FP_TYPE sum_o = 0.0;
	for (int t=0; t<ex_list.size(); t++)
	{
		tstep_out_begin(t);
		int o_s = ex_list.number_of_extrema(t);
		sum_o += o_s;
		for (int o1=0; o1<o_s; o1++)
		{
			// remove those greater than minimum delta
			FP_TYPE min_vd, max_vd;
			get_min_max_values_delta(min_vd, max_vd, o1, t);
			if (min_vd > min_delta)
			{
				ex_list.get(t, o1)->object_labels.clear();// delete!
				sum_o -= 1;
			}
		}
		tstep_out_end(t);	
	}
	std::cout << " Number of objects: " << sum_o << " ";
	std::cout << std::endl;
}

/******************************************************************************/

void minima_background::refine_objects(void)
{
	// refine the object (after the initial merge) so that any triangle which is not
	// the minimum value is removed from this initial object
	std::cout << "# Refining extrema, timestep: ";
	FP_TYPE mv = contour_data(data_minus_bck->get_missing_value(), contour_value);
	for (int t=0; t<ex_list.size(); t++)
	{
		tstep_out_begin(t);
		int o_s = ex_list.number_of_extrema(t);
		for (int o1=0; o1<o_s; o1++)
		{
			FP_TYPE min_v, max_v;
			get_min_max_values_delta(min_v, max_v, o1, t);
			min_v = contour_data(min_v, contour_value);
			LABEL_STORE* object_labels = &(ex_list.get(t, o1)->object_labels);
			LABEL_STORE new_labels;
			for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
				 it_ll != object_labels->end(); it_ll++)	// tri indices
			{
				indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
				FP_TYPE val = data_minus_bck->get_data(t, c_tri->get_ds_index());
				val = contour_data(val, contour_value);
				if (val < min_v + contour_value && val <= min_delta && !is_mv(val, mv))
					new_labels.push_back(*it_ll);
			}
			ex_list.get(t, o1)->object_labels = new_labels;
		}
		tstep_out_end(t);		
	}
	std::cout << std::endl;
}

/******************************************************************************/

void minima_background::get_min_max_values_delta(FP_TYPE& min_v, FP_TYPE& max_v, 
										 	     int o, int t)
{
	// get the minimum and maximum values for an object containing a number
	// of labels
	min_v = 2e20f;
	max_v = -2e20f;
	LABEL_STORE* object_labels = &(ex_list.get(t, o)->object_labels);
	// if there are no labels do not try to find the min/max
	if (object_labels->size() == 0)
		return;
	// get the missing value
	FP_TYPE mv = data_minus_bck->get_missing_value();
	// loop over all the triangles in the object
	for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
		 it_ll != object_labels->end(); it_ll++)	// tri indices
	{
		// get the triangle
		indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
		FP_TYPE val = data_minus_bck->get_data(t, c_tri->get_ds_index());
		// find the min and max values
		if (fabs(val) < 0.99*fabs(mv) && val > max_v && fabs(val) < 1e10)
			max_v = val;
		if (fabs(val) < 0.99*fabs(mv) && val < min_v && fabs(val) < 1e10)
			min_v = val;
	}	
}

/******************************************************************************/

void minima_background::calculate_object_position(int o, int t)
{
	steering_extremum* svex = ex_list.get(t, o);
	LABEL_STORE* object_labels = &(svex->object_labels);
	// find min / max of values in the object
	FP_TYPE min_v = 2e20f;	// minimum value in object
	FP_TYPE max_v = 0;		// maximum value in object
	get_min_max_values_delta(min_v, max_v, o, t);
	
	// get all the points within the contour
	min_v = contour_data(min_v, contour_value);
	
	// position vector in Cartesian coordinates
	vector_3D P;
	FP_TYPE sum_V = 0.0;
	// get the position and weight for each triangle in the object
	for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
		 it_ll != object_labels->end(); it_ll++)	// tri labels
	{
		// get the triangle and its centroid
		indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
		vector_3D C = c_tri->centroid();
		// get the data value from the datastore
		FP_TYPE V = data_minus_bck->get_data(t, c_tri->get_ds_index());
		// the value is the weight
		P += C * fabs(V);
		sum_V += fabs(V);
	}
	// divide by the weights and project to the sphere again
	if (sum_V > 0.0)
		P *= 1.0 / (P.mag() * sum_V);
	// put the values back in the geo extremum
	FP_TYPE lon, lat;
	cart_to_model(P, lon, lat);
	svex->lon = lon;
	svex->lat = lat;
}

/*****************************************************************************/

void minima_background::calculate_object_intensity(int o, int t)
{
	// find min / max of values in the object
	FP_TYPE min_v = 2e20f;		// minimum value in object
	FP_TYPE max_v = -2e20f;		// maximum value in object
	get_min_max_values(min_v, max_v, o, t);		
	steering_extremum* svex = ex_list.get(t, o);
	svex->intensity = min_v;
}

/*****************************************************************************/

void minima_background::calculate_object_delta(int o, int t)
{
	// calculate the delta as the absolute difference between the minimum
	// and the background field
	// find min / max of values in the object
	FP_TYPE min_v = 2e20f;	// minimum value in object
	FP_TYPE max_v = -2e20f;		// maximum value in object
	get_min_max_values_delta(min_v, max_v, o, t);	
	steering_extremum* svex = ex_list.get(t, o);
	svex->delta = min_v;
}