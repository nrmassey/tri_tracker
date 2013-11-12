/******************************************************************************
** Program : rotated_grid.h
** Author  : Neil Massey
** Date    : 05/06/13
** Purpose : class to represent a rotated grid
******************************************************************************/

#ifndef ROTATED_GRID_H
#define ROTATED_GRID_H

class rotated_grid
{
	public:
		rotated_grid(FP_TYPE rotated_pole_lat, FP_TYPE rotated_pole_lon,
					 class ncdata* parent_data);
		~rotated_grid(void);
		FP_TYPE get_global_latitude_value(int i, int j);
		FP_TYPE get_global_longitude_value(int i, int j);
	private:
		// this function fills in the global_latitude_values and global_longitude_values
		void calculate_global_coordinates(void);
		FP_TYPE rotated_pole_longitude;
		FP_TYPE rotated_pole_latitude;
		FP_TYPE* global_longitude_values;
		FP_TYPE* global_latitude_values;
		class ncdata* parent_data;
};

#endif