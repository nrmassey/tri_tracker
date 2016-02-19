/******************************************************************************
** Program : minima_back_wind.cpp
** Author  : Neil Massey
** Date    : 16/09/13
** Purpose : class inherited from minima_background that searches for minima 
**           in data regridded, after the removal of the background field
**           and then uses the wind gradient to expand objects to contain a
**           wind storm object
******************************************************************************/

#include "minima_back_wind.h"
#include <sstream>
#include <algorithm>
#include <math.h>
#include "concentric_shell.h"
#include "haversine.h"
#include "geo_convert.h"
#include "read_from_string.h"

/******************************************************************************/

minima_back_wind::minima_back_wind(void) : minima_background(),
										   wind_spd_field(NULL)
{
	perim=4;
}

/******************************************************************************/

minima_back_wind::~minima_back_wind(void)
{
	delete wind_spd_field;
}

/******************************************************************************/

void minima_back_wind::locate(void)
{
	process_data();
	find_extrema();
	if (grid_level > 5)
	{
		merge_objects();
		refine_objects();
	}
	find_objects();
	merge_objects();
	expand_objects();
	merge_objects();
	ex_points_from_objects();
}

/******************************************************************************/

bool minima_back_wind::process_data(void)
{
	wind_spd_field = new ncdata(wind_file_name, wind_spd_field_name);
	return minima_background::process_data();
}

/******************************************************************************/

void minima_back_wind::parse_arg_string(std::string method_string)
{
	// arguments are:
	// arg[0] = file to read background
	// arg[1] = level to use as the background field
	// arg[2] = averaging period to take background field over
	// arg[3] = mesh level to detect large scale minima at
	// arg[4] = contour value
	// arg[5] = minimum delta
	// arg[6] = file name containing wind field
	// arg[7] = field name of wind in file
	// arg[8] = percentile of value to use for threshold of winds
	// parameters for minima location with background removal
	// get the first bracket
	int c_pos = method_string.find_first_of("(")+1;
	int e_pos = method_string.find(")", c_pos);
	char dummy;
	if (method_string.substr(c_pos, e_pos-c_pos) == "help")
		throw(std::string("minima_back_wind parameters = (file name to take background field from, level in mesh to use as background field, averaging period of background field, contour value)"));
	
	// get the background field file
	c_pos = read_from_string(method_string, c_pos, ",", bck_field_file);
	// form a stream with the rest of the string in it
	std::stringstream stream(method_string.substr(c_pos, e_pos-c_pos));
	// read the other variables
	stream >> bck_mesh_lvl >> dummy 
		   >> bck_avg_period >> dummy
		   >> min_mesh_lvl >> dummy
		   >> contour_value >> dummy
		   >> min_delta >> dummy;
	// go to the end of the stream
	c_pos = c_pos + stream.tellg();

	// get the substring of the rest of the string
	method_string = method_string.substr(c_pos, e_pos-c_pos);
	c_pos = 0;
	c_pos = read_from_string(method_string, c_pos, ",", wind_file_name);
	c_pos = read_from_string(method_string, c_pos, ",", wind_spd_field_name);
	
	std::stringstream stream2(method_string.substr(c_pos, e_pos));	
	stream2 >> wind_thresh_value >> dummy
			>> wind_high_value;
	
	// add to the metadata
	std::stringstream ss;
	meta_data["method"] = "minima_back_wind";
	meta_data["background_file"] = bck_field_file;
	ss << bck_mesh_lvl;
	meta_data["background_mesh_level"] = ss.str();
	ss.str(""); ss << bck_avg_period;
	meta_data["background_averaging_period"] = ss.str();
	ss.str(""); ss << min_mesh_lvl;
	meta_data["minima_mesh_level"] = ss.str();
	ss.str(""); ss << contour_value;
	meta_data["contour_value"] = ss.str();
	ss.str(""); ss << min_delta;
	meta_data["minimum_delta"] = ss.str();
	meta_data["wind_file_name"] = wind_file_name;
	meta_data["wind_var_name"] = wind_spd_field_name;	
	ss.str(""); ss << wind_thresh_value;
	meta_data["wind_thresh_value"] = ss.str();
	ss.str(""); ss << wind_high_value;
}
				  
/******************************************************************************/

FP_TYPE calc_percentile(std::vector<FP_TYPE>* wnd_val, FP_TYPE ptile)
{
	// calculate the percentile of the wind values
	int v_size = wnd_val->size();
	int pos_0 = v_size * (ptile/100);
	std::vector<FP_TYPE>::iterator it_val = wnd_val->begin();
	std::advance(it_val, pos_0);
	FP_TYPE val_0 = *it_val;
	if (it_val != wnd_val->end())
		it_val ++;
	else
		it_val --;
	FP_TYPE val_1 = *it_val;
	FP_TYPE ptile_val = 0.5*(val_0 + val_1);
	it_val = wnd_val->end();
	return ptile_val;
}

/******************************************************************************/

bool minima_back_wind::is_in_object(indexed_force_tri_3D* O_TRI, 
				  					indexed_force_tri_3D* C_TRI, int t_step)
{
	// O_TRI - original triangle
	// C_TRI - candidate triangle - triangle being tested for inclusion
	// restrict to 1000km
	bool is_in = true;
	is_in = is_in && calculate_triangle_distance(O_TRI, C_TRI) < 1000.0;
	is_in = is_in && minima_background::is_in_object(O_TRI, C_TRI, t_step);
	return is_in;
}

/******************************************************************************/

FP_TYPE minima_back_wind::calculate_average_wind(indexed_force_tri_3D* TRI, int t_step)
{
	const std::list<grid_index>* c_tri_idxs = TRI->get_grid_indices();
	FP_TYPE c_wnd_spd = 0.0;
	// candidate triangle wind speed
	for (std::list<grid_index>::const_iterator c_t_idx = c_tri_idxs->begin();
		 c_t_idx != c_tri_idxs->end(); c_t_idx++)
		c_wnd_spd += wind_spd_field->get_data(c_t_idx->i, c_t_idx->j, 0, t_step);
		
	if (c_tri_idxs->size() == 0)
		return -1;
	else
		c_wnd_spd /= c_tri_idxs->size();
	return c_wnd_spd;
}

/******************************************************************************/

bool minima_back_wind::wind_test(indexed_force_tri_3D* O_TRI, 
								 indexed_force_tri_3D* C_TRI, int t_step)
{
	// distance check - distance from original triangle to current triangle
	// should be less than 1000km
	const FP_TYPE dist_thresh = 1000;
	FP_TYPE c_wnd_spd = calculate_average_wind(C_TRI, t_step);
	bool wnd_test = c_wnd_spd >= wind_thresh_value;
	FP_TYPE dist = calculate_triangle_distance(O_TRI, C_TRI);
	wnd_test = wnd_test && dist < dist_thresh;
	// how much over the wind_thresh? - special case for very strong winds on
	// edge of 1000km radius
	if (c_wnd_spd >= wind_high_value && dist < 1500)
		wnd_test = true;
	return wnd_test;
}

/******************************************************************************/

void minima_back_wind::expand_objects(void)
{
	std::cout << "# Expanding objects, timestep: ";
	concentric_shell c_shell;
	LABEL_STORE shell_in_object;
	for (int t=0; t<ds.get_number_of_time_steps(); t++)
	{
		tstep_out_begin(t);
		for (int e=0; e<ex_list.number_of_extrema(t); e++)
		{
			// reset shell
			c_shell.clear();
			// get the extremum for this object
			steering_extremum* svex = ex_list.get(t, e);
			// check to see whether a label actually exists
			if (svex->object_labels.size() == 0)
				continue;
			// get the "original triangle" for this object 
			// - i.e. the one at the centre of the object
			indexed_force_tri_3D* O_TRI = get_original_triangle(e, t);
			// continue until no more triangles are added to the object
			c_shell.calculate_inner_ring(&tg, svex);	// calculate the initial shape
			// calculate the 1st concentric shell
			c_shell.calculate(&tg, svex);
			bool keep_adding_to_obj = true;
			while (keep_adding_to_obj)
			{
				keep_adding_to_obj = false;
				LABEL_STORE* c_shell_labs = c_shell.get_labels();
				shell_in_object.clear();
				// loop over all the labels
				for (LABEL_STORE::iterator it_c_shell_labs = c_shell_labs->begin();
					 it_c_shell_labs != c_shell_labs->end(); it_c_shell_labs++)
				{
					// get the candidate triangle
					indexed_force_tri_3D* C_TRI = tg.get_triangle(*it_c_shell_labs);
					FP_TYPE V = ds.get_data(t, C_TRI->get_ds_index());
					if (wind_test(O_TRI, C_TRI, t) && fabs(V) < fabs(ds.get_missing_value()*0.99) &&
					    V > ds.get_data(t, O_TRI->get_ds_index()) )
					{
						keep_adding_to_obj = true;
						// add to object but only if not already added
						if (std::find(svex->object_labels.begin(), svex->object_labels.end(), 
							*it_c_shell_labs) == svex->object_labels.end())
						{
							svex->object_labels.push_back(*it_c_shell_labs);
						}
						if (std::find(shell_in_object.begin(), shell_in_object.end(), 
							*it_c_shell_labs) == shell_in_object.end())
						{
							shell_in_object.push_back(*it_c_shell_labs);							
						}						
					}
				}
				c_shell.recalculate(&tg, &shell_in_object);
			}			
		}
		tstep_out_end(t);
	}
	std::cout << std::endl;	
}

/*****************************************************************************/

void minima_back_wind::calculate_object_position(int o, int t)
{
	steering_extremum* svex = ex_list.get(t, o);
	LABEL_STORE* object_labels = &(svex->object_labels);
	// find min / max of values in the object
	FP_TYPE min_v = 2e20;

	// position vector in Cartesian coordinates
	vector_3D P;
	FP_TYPE sum_w = 0.0;
	// get the position and weight for each triangle in the object
	for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
		 it_ll != object_labels->end(); it_ll++)	// tri labels
	{
		// get the triangle and its centroid
		indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
		// get the data value from the datastore and the centroid
		FP_TYPE V = ds.get_data(t, c_tri->get_ds_index());
		vector_3D C = c_tri->centroid();
		if (V < min_v && fabs(V) < 20000)
		{
			P = C;
			min_v = V;
		}
	}

	// put the values back in the geo extremum
	FP_TYPE lon, lat;
	cart_to_model(P, lon, lat);
	svex->lon = lon;
	svex->lat = lat;
}

/*****************************************************************************/

void minima_back_wind::calculate_object_intensity(int o, int t)
{
	// object intensity is calculated as a weighted sum of the intensities of
	// the triangles in the object.  The weight is defined as 1.0 - dist/max_dist
	// where dist is the distance between the lat and lon of the feature point
	// and max_dist is the maximum distance of all the objects
	
	// get the extremum - the lat and lon will have been set already
	steering_extremum* svex = ex_list.get(t, o);
	LABEL_STORE* object_labels = &(svex->object_labels);
	
	// find min / max of values in the object
	FP_TYPE min_v = 2e20f;	// minimum value in object
	
		// get the position and weight for each triangle in the object
	for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
		 it_ll != object_labels->end(); it_ll++)	// tri labels
	{
		// get the data value from the datastore
		indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
		FP_TYPE V = ds.get_data(t, c_tri->get_ds_index());
		if (V < min_v && fabs(V) < 20000)
		{
			min_v = V;
		}
	}
	
	svex->intensity = min_v;
}
