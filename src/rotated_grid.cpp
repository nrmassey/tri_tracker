/******************************************************************************
** Program : rotated_grid.cpp
** Author  : Neil Massey
** Date    : 05/06/13
** Purpose : class to represent a rotated grid
******************************************************************************/

#include "rotated_grid.h"
#include "Rot2Global.h"
#include "ncdata.h"

/*****************************************************************************/

rotated_grid::rotated_grid(FP_TYPE i_rotated_pole_lat, FP_TYPE i_rotated_pole_lon,
	 			  			class ncdata* p_parent_data) :
	 			  			rotated_pole_longitude(i_rotated_pole_lon),
	 			  			rotated_pole_latitude(i_rotated_pole_lat),
	 			  			global_longitude_values(NULL),
	 			  			global_latitude_values(NULL),
	 			  			parent_data(p_parent_data)
{
	// after the variables are initialised, calculate the global coordinates
	// from the rotated pole coordinates
	calculate_global_coordinates();
}

/*****************************************************************************/

rotated_grid::~rotated_grid(void)
{
	delete [] global_latitude_values;
	delete [] global_longitude_values;	
}

/*****************************************************************************/

FP_TYPE rotated_grid::get_global_latitude_value(int i, int j)
{
	int off = j*parent_data->get_lon_len() + i;
	return global_latitude_values[off];
}

/*****************************************************************************/

FP_TYPE rotated_grid::get_global_longitude_value(int i, int j)
{
	int off = j*parent_data->get_lon_len() + i;
	return global_longitude_values[off];
}

/*****************************************************************************/

std::list<FP_TYPE> rotated_grid::get_global_grid_box_values(int i, int j)
{
	// return an 8-long vector containing the four grid box corners
	// in lon0,lat0, lon1,lat1, lon2,lat2, lon3,lat3 order
	// box is defined in a clockwise order
	
	// get the midpoint and calculate
	FP_TYPE f_mid_lat = parent_data->get_lat_s() + j*parent_data->get_lat_d();
	FP_TYPE f_mid_lon = parent_data->get_lon_s() + i*parent_data->get_lon_d();
	FP_TYPE f_half_lat_d = parent_data->get_lat_d() * 0.5;
	FP_TYPE f_half_lon_d = parent_data->get_lon_d() * 0.5;
	std::list<FP_TYPE> output_list;
	
	// offsets from the mid point in terms of multiples of half_lon_d, half_lat_d
	const FP_TYPE x_map[4] = {-1.0,1.0,1.0,-1.0};
	const FP_TYPE y_map[4] = {-1.0,-1.0,1.0,1.0};
	
	for (int x=0; x<4; x++)
	{
		// calculate the grid box corner in the rotated grid
		FP_TYPE rot_lon = f_mid_lon + x_map[x] * f_half_lon_d;
		FP_TYPE rot_lat = f_mid_lat + y_map[x] * f_half_lat_d;
		// convert to global lon / lat
		FP_TYPE global_lon, global_lat;
		Rot2Global(rot_lat, rot_lon,
				   rotated_pole_latitude, rotated_pole_longitude,
				   global_lat, global_lon);
		// append to the vector
		output_list.push_back(global_lon);
		output_list.push_back(global_lat);
	}
	return output_list;
}

/*****************************************************************************/

FP_TYPE rotated_grid::get_rotated_pole_longitude(void)
{
	return rotated_pole_longitude;
}

/*****************************************************************************/

FP_TYPE rotated_grid::get_rotated_pole_latitude(void)
{
	return rotated_pole_latitude;
}

/*****************************************************************************/

void rotated_grid::calculate_global_coordinates(void)
{
	// create the storage for the global coordinates
	int size = parent_data->get_lon_len() * parent_data->get_lat_len();
	global_latitude_values = new FP_TYPE[size];
	global_longitude_values = new FP_TYPE[size];
	// loop through the latitude and longitude values - lats first
	for (int j=0; j<parent_data->get_lat_len(); j++)
	{
		for (int i=0; i<parent_data->get_lon_len(); i++)
		{
			FP_TYPE f_global_lat, f_global_lon;
			FP_TYPE f_rot_lat = parent_data->get_lat_s() + j*parent_data->get_lat_d();
			FP_TYPE f_rot_lon = parent_data->get_lon_s() + i*parent_data->get_lon_d();
			Rot2Global(f_rot_lat, f_rot_lon,
					   rotated_pole_latitude, rotated_pole_longitude,
					   f_global_lat, f_global_lon);
			// offset is j*ncols + i
			int off = j*parent_data->get_lon_len() + i;
			global_latitude_values[off] = f_global_lat;
			global_longitude_values[off] = f_global_lon;
		}
	}
}