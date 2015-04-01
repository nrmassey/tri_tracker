/******************************************************************************
** Program : geo_wind_vector.cpp
** Author  : Neil Massey
** Date    : 05/08/13
** Purpose : Calculate the geostrophic wind from a file passed in.  Inherits
**           from steering vector class
** Modified: 18/08/14 - updated to use spherical coordinate version of 
**           geostrophic wind equation and to directly calculated the 
**           derivatives from the 5 point splines
******************************************************************************/

#include "geo_wind_vector.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include "haversine.h"
#include "spline.h"

/*****************************************************************************/

const FP_TYPE f_deg_to_rad = M_PI/180.0;

/*****************************************************************************/

geo_wind_vector::geo_wind_vector(void) : geopot_ht(NULL)
{
}

/*****************************************************************************/

geo_wind_vector::~geo_wind_vector(void)
{
	delete geopot_ht;
}

/*****************************************************************************/

void geo_wind_vector::parse_arg_string(std::string arg_string)
{
	// parse the argument string.  Format is:
	// geostropic(file_name, var_name, z_level)
	std::string file_name, var_name;
	int z_level;
	
	// get the filename
	int c_pos = arg_string.find_first_of("(")+1;
	int e_pos = arg_string.find(",", c_pos);
	file_name = arg_string.substr(c_pos, e_pos-c_pos);
	
	// get the varname
	c_pos = e_pos+1;
	e_pos = arg_string.find(",", c_pos);
	var_name = arg_string.substr(c_pos, e_pos-c_pos);

	// get the level
	c_pos = e_pos+1;
	e_pos = arg_string.find(")", c_pos);
	std::stringstream stream(arg_string.substr(c_pos, e_pos-c_pos));
	stream >> z_level;
		
	// have the file name, var name and level so create the nc data
	geopot_ht = new ncdata(file_name, var_name);
	
	// add the meta_data
	meta_data["steering_wind_method"] = "geostrophic wind";
	meta_data["steering_wind_file_name"] = file_name;
	meta_data["steering_wind_var_name"] = var_name;
	std::stringstream ss;
	ss << z_level;
	meta_data["steering_wind_z_level"] = ss.str();
}

/*****************************************************************************/

void get_lon_lat_idx(ncdata* nc_data, int lon_idx, int lat_idx,
                     int lon_off, int lat_off, int& n_lon_idx, int& n_lat_idx)
{
    // calculate new longitude and latitude indices based on the supplied
    // lon lat and the offsets
    n_lon_idx = lon_idx + lon_off;
    n_lat_idx = lat_idx + lat_off;

    // wrap around the North Pole if necessary and not a rotated grid!
    if (n_lat_idx < 0)
    {
    	if (nc_data->has_rotated_grid())
    	{
    		n_lat_idx = 0;
    	}
    	else
    	{
	        n_lat_idx = -n_lat_idx - 1;
    	    n_lon_idx += n_lon_idx / 2;
    	}
    }

    // wrap around the South Pole if necessary
    int lat_len = nc_data->get_lat_len();
    if (n_lat_idx >= lat_len)
    {
    	if (nc_data->has_rotated_grid())
    	{
    		n_lat_idx = lat_len-1;
    	}
    	else
    	{
	        n_lat_idx = 2 * lat_len - n_lat_idx - 1;
    	    n_lon_idx += n_lon_idx / 2;
    	}
    }

    // wrap around the date line
    int lon_len = nc_data->get_lon_len();
    if (n_lon_idx < 0)
    {
    	if (nc_data->has_rotated_grid())
    		n_lon_idx = 0;
    	else
	        n_lon_idx = lon_len + n_lon_idx;
    }
    if (n_lon_idx >= lon_len)
    {
    	if (nc_data->has_rotated_grid())
    		n_lon_idx = lon_len - 1;
	    else
	        n_lon_idx = n_lon_idx - lon_len;
	}
}

/*****************************************************************************/

spline form_lon_spline(ncdata* nc_data, int lon_idx, int lat_idx, int z, int t)
{
	const int spl_size = 10;
	int spl2 = spl_size / 2;
    // construct a 5 pt spline from the data in the lon direction
    std::vector<FP_TYPE> lon_vals(spl_size, 0.0);
	std::vector<FP_TYPE> x_vals(spl_size);
    for (int i=-spl2; i<spl2+1; i++)
    {
        int lat_i, lon_i;
        get_lon_lat_idx(nc_data, lon_idx, lat_idx, i, 0, lon_i, lat_i);
        lon_vals[i+spl2] = nc_data->get_data(lon_i, lat_i, z, t);
        if (nc_data->has_rotated_grid())
	        x_vals[i+spl2] = nc_data->get_rotated_grid()->get_global_longitude_value(lon_i, lat_i) * f_deg_to_rad;
        else
	        x_vals[i+spl2] = nc_data->get_lon_from_idx(lon_i)* f_deg_to_rad;
    }
    spline lon_spline(lon_vals, x_vals, nc_data->get_missing_value());
    return lon_spline;
}

/*****************************************************************************/

spline form_lat_spline(ncdata* nc_data, int lon_idx, int lat_idx, int z, int t)
{
	const int spl_size = 5;
	int spl2 = spl_size / 2;
    // construct a 5 pt spline from the data in the lat direction
    std::vector<FP_TYPE> lat_vals(spl_size, 0.0);
	std::vector<FP_TYPE> y_vals(spl_size);
	
    for (int j=-spl2; j<spl2+1; j++)
    {
        int lat_i, lon_i;
        get_lon_lat_idx(nc_data, lon_idx, lat_idx, 0, j, lon_i, lat_i);
        // put the data values in backwards
        lat_vals[spl2-j] = nc_data->get_data(lon_i, lat_i, z, t);
        if (nc_data->has_rotated_grid())
	        y_vals[spl2-j] = nc_data->get_rotated_grid()->get_global_latitude_value(lon_i, lat_i) * f_deg_to_rad;
        else
	        y_vals[spl2-j] = nc_data->get_lat_from_idx(lat_i) * f_deg_to_rad;
    }
    spline lat_spline(lat_vals, y_vals, nc_data->get_missing_value());
    return lat_spline;
}

/*****************************************************************************/

void calc_geo_wind(ncdata* gph_data, int t, int gph_z, int lon_idx, int lat_idx,
                   FP_TYPE mv, FP_TYPE& u, FP_TYPE& v)
{
    // calculate the geostrophic wind - spherical coordinates
	// determine the scale factor - do we need to convert geopotential height
	// to the geopotential?
	FP_TYPE g_sc = 1.0;
	if (gph_data->get_units() == "m")
		g_sc = 9.81;
	
    spline lon_spline = form_lon_spline(gph_data, lon_idx, lat_idx, gph_z, t);
    spline lat_spline = form_lat_spline(gph_data, lon_idx, lat_idx, gph_z, t);

    FP_TYPE lat_r;
    FP_TYPE lon;
    
    if (gph_data->has_rotated_grid())
    {
    	lat_r = gph_data->get_rotated_grid()->get_global_latitude_value (lon_idx, lat_idx) * f_deg_to_rad;
    	lon   = gph_data->get_rotated_grid()->get_global_longitude_value(lon_idx, lat_idx) * f_deg_to_rad;
    }
    else
    {
	    lat_r = gph_data->get_lat_from_idx(lat_idx) * f_deg_to_rad;
	    lon   = gph_data->get_lon_from_idx(lon_idx) * f_deg_to_rad;
	}

    // calculate the Coriolis parameter
    const FP_TYPE O = 7.292e-5;
    FP_TYPE f = 2 * O * sin(lat_r);
    // radius of the Earth
    const FP_TYPE R = 6371 * 1000;
    // calculate the u and v of the wind using the derivatives from the spline
    // u and v are now in spherical coordinates
    u = -g_sc / (f*R) * lat_spline.evaluate_dx(lat_r); 
    v = g_sc / (f*R*cos(lat_r)) * lon_spline.evaluate_dx(lon);
}

/*****************************************************************************/

void geo_wind_vector::calculate_steering_vector(tri_grid* tg, 
		 										steering_extremum* svex, int t, 
												FP_TYPE mv)
{
	LABEL_STORE* object_labels = &(svex->object_labels);	
	// geostrophic wind is calculated as an average of the geostrophic wind
	// for each triangle in the object
	FP_TYPE sv_u = 0.0, sv_v = 0.0;
	FP_TYPE sum_sv_u = 0.0, sum_sv_v = 0.0;
	int n_sv_pts = 0;
	// loop through the triangle objects
	for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
		 it_ll != object_labels->end(); it_ll++)
	{
		// get the triangle
		indexed_force_tri_3D* c_tri = tg->get_triangle(*it_ll);
		// get the indices
		const std::list<grid_index>* g_idx = c_tri->get_grid_indices();
		for (std::list<grid_index>::const_iterator it_g_idx = g_idx->begin();
			 it_g_idx != g_idx->end(); it_g_idx++)
		{
			// use the x and y index into the original grid from the triangle to
			// get the data from the original grid and calculate the geostrophic wind
			calc_geo_wind(geopot_ht, t, z_level, it_g_idx->i, it_g_idx->j, mv, sv_u, sv_v);
		    if (sv_u != mv && sv_v != mv && isfinite(sv_u) && isfinite(sv_v))
		    {
		    	sum_sv_u += sv_u;
		    	sum_sv_v += sv_v;
		    	n_sv_pts++;
		    }
        }
	}
	if (n_sv_pts != 0)
	{
		svex->sv_u = sum_sv_u / n_sv_pts;
		svex->sv_v = sum_sv_v / n_sv_pts;
	}
}
