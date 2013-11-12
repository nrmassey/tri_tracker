/******************************************************************************
** Program : geo_wind_vector.cpp
** Author  : Neil Massey
** Date    : 05/08/13
** Purpose : Calculate the geostrophic wind from a file passed in.  Inherits
**           from steering vector class
******************************************************************************/

#include "geo_wind_vector.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include "haversine.h"
#include "spline.h"

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

    // wrap around the North Pole if necessary
    if (n_lat_idx < 0)
    {
        n_lat_idx = -n_lat_idx - 1;
        n_lon_idx += n_lon_idx / 2;
    }

    // wrap around the South Pole if necessary
    int lat_len = nc_data->get_lat_len();
    if (n_lat_idx >= lat_len)
    {
        n_lat_idx = 2 * lat_len - n_lat_idx - 1;
        n_lon_idx += n_lon_idx / 2;
    }

    // wrap around the date line
    int lon_len = nc_data->get_lon_len();
    if (n_lon_idx < 0)
        n_lon_idx = lon_len + n_lon_idx;
    if (n_lon_idx >= lon_len)
        n_lon_idx = n_lon_idx - lon_len;
}

/*****************************************************************************/

spline form_lon_spline(ncdata* nc_data, int lon_idx, int lat_idx, int z, int t)
{
    // construct a 5 pt spline from the data in the lon direction
    std::vector<FP_TYPE> lon_vals(5, 0.0);
    for (int i=-2; i<3; i++)
    {
        int lat_i, lon_i;
        get_lon_lat_idx(nc_data, lon_idx, lat_idx, i, 0, lon_i, lat_i);
        lon_vals[i+2] = nc_data->get_data(lon_i, lat_i, z, t);
    }
	// check to see if any of the values are the missing value
	const int l[5] = {4,3,2,3,4};
	for (int i=0; i<5; i++)
	{
		if (lon_vals[i] == nc_data->get_missing_value())
		{
			int dir = 1;
			if (i > 2)
				dir = -1;
			for (int j=1; j<l[i]; j++)
			{
				int p = i + j*dir;
				if (lon_vals[p] != nc_data->get_missing_value())
				{
					lon_vals[i] = lon_vals[p];
					break;
				}
			}
		}
	}
	std::vector<FP_TYPE> x_vals(5);
	for (int i=0; i<5; i++)
		x_vals[i] = i;
    spline lon_spline(lon_vals, x_vals, nc_data->get_missing_value());
    return lon_spline;
}

/*****************************************************************************/

spline form_lat_spline(ncdata* nc_data, int lon_idx, int lat_idx, int z, int t)
{
    // construct a 5 pt spline from the data in the lat direction
    std::vector<FP_TYPE> lat_vals(5, 0.0);
    for (int j=-2; j<3; j++)
    {
        int lat_i, lon_i;
        get_lon_lat_idx(nc_data, lon_idx, lat_idx, 0, j, lon_i, lat_i);
        lat_vals[j+2] = nc_data->get_data(lon_i, lat_i, z, t);
    }
    // check to see if any of the values are the missing value
    const int l[5] = {4,3,2,3,4};
    for (int j=0; j<5; j++)
    {
        if (lat_vals[j] == nc_data->get_missing_value())
        {
            int dir = 1;
            if (j > 2)
                dir = -1;
            for (int k=1; k<l[j]; k++)
            {
                int p = j + k*dir;
                if (lat_vals[p] != nc_data->get_missing_value())
                {
                    lat_vals[j] = lat_vals[p];
                    break;
                }
            }
        }
    }
	std::vector<FP_TYPE> y_vals(5);
	for (int i=0; i<5; i++)
		y_vals[i] = i;

    spline lat_spline(lat_vals, y_vals, nc_data->get_missing_value());
    return lat_spline;
}

/*****************************************************************************/

FP_TYPE calc_five_pt_difference(spline& spl, FP_TYPE mv)
{
    // calculate centred difference from the 5 point spline
    // set h to 1 to use the actual values
    const FP_TYPE h = 0.1;
    const FP_TYPE y = 2.0;

    // calculate the terms in the difference
    FP_TYPE t1 = -spl.evaluate(y+2*h);   // first term -S(y+2h)
    FP_TYPE t2 = 8*spl.evaluate(y+h);    // 8S(y+h)
    FP_TYPE t3 = -8*spl.evaluate(y-h);   // -8S(y-h)
    FP_TYPE t4 = spl.evaluate(y-2*h);    // S(y-2h)

	FP_TYPE t;
	if (t1 == mv || t2 == mv || t3 == mv || t4 == mv)
		t = mv;
	else
		t = (t1+t2+t3+t4) / (12*h);
    return t;
}

/*****************************************************************************/

void calc_geo_wind(ncdata* gph_data, int t, int gph_z,
				   int lon_idx, int lat_idx,
                   FP_TYPE& u, FP_TYPE& v)
{
    // calculate the geostrophic wind
	// determine the scale factor - do we need to convert geopotential height
	// to the geopotential?
	FP_TYPE g_sc = 1.0;
	if (gph_data->get_units() == "m")
		g_sc = 9.81;
	
    spline lon_spline = form_lon_spline(gph_data, lon_idx, lat_idx, gph_z, t);
    spline lat_spline = form_lat_spline(gph_data, lon_idx, lat_idx, gph_z, t);

    FP_TYPE lat = gph_data->get_lat_from_idx(lat_idx);
    FP_TYPE lon = gph_data->get_lon_from_idx(lon_idx);

    // calculate the derivatives
	FP_TYPE mv = gph_data->get_missing_value();
    FP_TYPE dphi_dx = calc_five_pt_difference(lon_spline, mv);
    FP_TYPE dphi_dy = calc_five_pt_difference(lat_spline, mv);
    FP_TYPE f = 2 * 7.292e-5 * sin(lat * M_PI/180.0);

    // calculate the width and height of the grid boxes for this latitude
    FP_TYPE dy = haversine(0.0, lat-gph_data->get_lat_d()/2.0,
                          0.0, lat+gph_data->get_lat_d()/2.0, EARTH_R);
    FP_TYPE dx = haversine(lon, lat, lon+gph_data->get_lon_d(), lat, EARTH_R);

	if (dphi_dx == mv || dphi_dy == mv)
		u = v = mv;
	else if (f != 0.0)
	{
	    u = g_sc/f * dphi_dy * 1.0/dy;		// latitude running opposite way to equations
    	v = g_sc/f * dphi_dx * 1.0/dx;
	}
	else
		u = v = 0.0;
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
			calc_geo_wind(geopot_ht, t, z_level, it_g_idx->i, it_g_idx->j, sv_u, sv_v);
		    if (sv_u != mv && sv_v != mv)
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
