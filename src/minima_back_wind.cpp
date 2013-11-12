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

/******************************************************************************/

minima_back_wind::minima_back_wind(void) : wind_spd_field(NULL),
										   minima_background()
{
	previous_object = -1;
}

/******************************************************************************/

minima_back_wind::~minima_back_wind(void)
{
	delete wind_spd_field;
}

/******************************************************************************/

bool minima_back_wind::process_data(void)
{
	// load the netcdf data and then call the base class process_data
	wind_spd_field = new ncdata(wind_file_name, wind_spd_field_name);
	// calculate the distribution of wind speeds per timestep	
	for (int t=0; t<wind_spd_field->get_t_len(); t++)
	{
		// calculate the wind speed at each grid point for this time step
		std::vector<FP_TYPE> wind_for_t_step;
		for (int i=0; i<wind_spd_field->get_lon_len(); i++)
		{
			for (int j=0; j<wind_spd_field->get_lat_len(); j++)
			{
				FP_TYPE wnd_speed = wind_spd_field->get_data(i,j,0,t);
				wind_for_t_step.push_back(wnd_speed);
			}
		}
		std::sort(wind_for_t_step.begin(), wind_for_t_step.end());
		// create the storage
		std::vector<FP_TYPE> percentile_winds;
		// now calculate the distribution in 1% percentiles
		int s = wind_for_t_step.size();
		for (int p=0; p<100; p++)
		{
			float ptile = float(p)/100.0;
			int idx = int(float(s) * ptile);
			float frac_part = float(s) * ptile - int(float(s) * ptile);
			float value = wind_for_t_step[idx] * (1.0-frac_part) + wind_for_t_step[idx+1] * frac_part;
			percentile_winds.push_back(value);
		}
		all_wind_distr.push_back(percentile_winds);
	}
	return minima_background::process_data();
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
	if (grid_level > 5)
		trim_objects();
	merge_objects();
	ex_points_from_objects();
	expand_objects();
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
}

/******************************************************************************/

FP_TYPE calc_percentile(std::vector<FP_TYPE>& V, FP_TYPE p)
{
	int s = V.size();
	int idx_5 = float(s) * p;
	float frac_part = float(s) * p - int(float(s) * p);
	float value = V[idx_5] * (1.0-frac_part) + V[idx_5+1] * frac_part;
	return value;
}

/******************************************************************************/

bool minima_back_wind::wind_test(indexed_force_tri_3D* O_TRI, 
								 indexed_force_tri_3D* C_TRI, int t_step)
{
	const std::list<grid_index>* c_tri_idxs = C_TRI->get_grid_indices();
	FP_TYPE c_wnd_spd = 0.0;
	// candidate triangle wind speed
	for (std::list<grid_index>::const_iterator c_t_idx = c_tri_idxs->begin();
		 c_t_idx != c_tri_idxs->end(); c_t_idx++)
	{
		c_wnd_spd += wind_spd_field->get_data(c_t_idx->i, c_t_idx->j, 0, t_step);
	}
	if (c_tri_idxs->size() == 0.0)
		c_wnd_spd = -1.0;
	else
		c_wnd_spd /= c_tri_idxs->size();
	
	bool wnd_test = false;
	wnd_test = c_wnd_spd >= 0.0 && c_wnd_spd >= all_wind_distr[t_step][80];
	return wnd_test;
}

/******************************************************************************/

void minima_back_wind::expand_objects(void)
{
	//
	std::cout << "# Expanding objects, timestep: ";
	concentric_shell c_shell;
	LABEL_STORE shell_in_object;
	for (int t=0; t<ds.get_number_of_time_steps(); t++)
	{
		tstep_out_begin(t);
		for (int e=0; e<ex_list.number_of_extrema(t); e++)
		{
			c_shell.clear();
			// get the extremum for this object
			steering_extremum* svex = ex_list.get(t, e);
			// check to see whether a label actually exists
			if (svex->object_labels.size() == 0)
				continue;
			// get the "original triangle" for this object 
			// - i.e. the one at the centre of the object
			indexed_force_tri_3D* O_TRI = get_original_triangle(e, t);
			bool keep_adding_to_obj = true;
			// continue until no more triangles are added to the object
			c_shell.calculate_inner_ring(&tg, svex);	// calculate the initial shape
			// calculate the 1st concentric shell
			c_shell.calculate(&tg, svex);
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
					if (wind_test(O_TRI, C_TRI, t))
					{
						keep_adding_to_obj = true;
						// add to object but only if not already added
						if (std::find(svex->object_labels.begin(), svex->object_labels.end(), 
							*it_c_shell_labs) == svex->object_labels.end())
						{
							svex->object_labels.push_back(*it_c_shell_labs);
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

/******************************************************************************/