/******************************************************************************
** Program : mslp_wind_vector.cpp
** Author  : Neil Massey
** Date    : 05/08/13
** Purpose : Calculate the geostrophic wind from the MSLP of a file passed in.
**           Inherits from steering vector class
** Modified: 18/08/14 - updated to use spherical coordinate version of 
**           geostrophic wind equation and to directly calculated the 
**           derivatives from the 5 point splines
******************************************************************************/

#include "mslp_wind_vector.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include "haversine.h"
#include "geo_convert.h"
#include "spline.h"

/*****************************************************************************/

const FP_TYPE f_deg_to_rad = M_PI/180.0;

/*****************************************************************************/

mslp_wind_vector::mslp_wind_vector(void) : mslp(NULL)
{
    max_tri_level = -1;
}

/*****************************************************************************/

mslp_wind_vector::~mslp_wind_vector(void)
{
    delete mslp;
}

/*****************************************************************************/

void mslp_wind_vector::parse_arg_string(std::string arg_string)
{
    // parse the argument string.  Format is:
    // mslp_wind(file_name, var_name)
    std::string file_name, var_name;
    
    // get the filename
    int c_pos = arg_string.find_first_of("(")+1;
    int e_pos = arg_string.find(",", c_pos);
    file_name = arg_string.substr(c_pos, e_pos-c_pos);
    
    // get the varname
    c_pos = e_pos+1;
    e_pos = arg_string.find(")", c_pos);
    var_name = arg_string.substr(c_pos, e_pos-c_pos);

    // have the file name, var name and level so create the nc data
    mslp = new ncdata(file_name, var_name);
    
    // add the meta_data
    meta_data["steering_wind_method"] = "mslp geostrophic wind";
    meta_data["steering_wind_file_name"] = file_name;
    meta_data["steering_wind_var_name"] = var_name;
}

/*****************************************************************************/

extern void get_lon_lat_idx(ncdata* nc_data, int lon_idx, int lat_idx,
                            int lon_off, int lat_off, int& n_lon_idx, int& n_lat_idx);
extern spline form_lon_spline(ncdata* nc_data, int lon_idx, int lat_idx, int z, int t);
extern spline form_lat_spline(ncdata* nc_data, int lon_idx, int lat_idx, int z, int t);

/*****************************************************************************/

void calc_mslp_wind(ncdata* mslp_data, int t, int lon_idx, int lat_idx,
                    FP_TYPE mv, FP_TYPE& u, FP_TYPE& v)
{
    // calculate the geostrophic wind - spherical coordinates, from MSLP
    
    spline lon_spline = form_lon_spline(mslp_data, lon_idx, lat_idx, 0, t);
    spline lat_spline = form_lat_spline(mslp_data, lon_idx, lat_idx, 0, t);

    FP_TYPE lat_r;
    FP_TYPE lon;
    
    if (mslp_data->has_rotated_grid())
    {
        lat_r = mslp_data->get_rotated_grid()->get_global_latitude_value (lon_idx, lat_idx) * f_deg_to_rad;
        lon   = mslp_data->get_rotated_grid()->get_global_longitude_value(lon_idx, lat_idx) * f_deg_to_rad;
    }
    else
    {
        lat_r = mslp_data->get_lat_from_idx(lat_idx) * f_deg_to_rad;
        lon   = mslp_data->get_lon_from_idx(lon_idx) * f_deg_to_rad;
    }

    // calculate the Coriolis parameter
    const FP_TYPE O = 7.292e-5;
    FP_TYPE f = 2 * O * sin(lat_r);
    // radius of the Earth
    const FP_TYPE R = 6371 * 1000;
    // reciprocal of density of dry air km/m^3
    const FP_TYPE roh = 1.29;
    // calculate the u and v of the wind using the derivatives from the spline
    // u and v are now in spherical coordinates
    u =  -1.0 / (roh*f*R) * lat_spline.evaluate_dx(lat_r); 
    v =  1.0 / (roh*f*R*cos(lat_r)) * lon_spline.evaluate_dx(lon);
}

/*****************************************************************************/

void mslp_wind_vector::calculate_steering_vector(tri_grid* tg, 
                                                steering_extremum* svex, int t, 
                                                FP_TYPE mv)
{
    // find where the svex occurs
    vector_3D C = model_to_cart(svex->lon, svex->lat);
    if (max_tri_level == -1)
        max_tri_level = tg->get_max_level();
    LABEL mid_L = tg->get_triangle_for_point(&C, max_tri_level);
    
    // now get the triangle and the adjacent labels
    indexed_force_tri_3D* c_tri = tg->get_triangle(mid_L);
    const LABEL_STORE* object_labels = c_tri->get_adjacent_labels(POINT);
    // geostrophic wind is calculated as an average of the geostrophic wind
    // for each triangle in the adjacent triangle list
    FP_TYPE sv_u = 0.0, sv_v = 0.0;
    FP_TYPE sum_sv_u = 0.0, sum_sv_v = 0.0;
    int n_sv_pts = 0;
    // calculate the geo wind for the first triangle
    // get the triangle
    // get the indices
    const std::list<grid_index>* g_idx = c_tri->get_grid_indices();
    for (std::list<grid_index>::const_iterator it_g_idx = g_idx->begin();
         it_g_idx != g_idx->end(); it_g_idx++)
    {
        // use the x and y index into the original grid from the triangle to
        // get the data from the original grid and calculate the geostrophic wind
        calc_mslp_wind(mslp, t, it_g_idx->i, it_g_idx->j, mv, sv_u, sv_v);
        if (sv_u != mv && sv_v != mv && isfinite(sv_u) && isfinite(sv_v))
        {
            sum_sv_u += sv_u;
            sum_sv_v += sv_v;
            n_sv_pts++;
        }
    }
    
    // loop through the adjacent triangle
    for (LABEL_STORE::const_iterator it_ll = object_labels->begin(); 
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
            calc_mslp_wind(mslp, t, it_g_idx->i, it_g_idx->j, mv, sv_u, sv_v);
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
