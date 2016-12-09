/******************************************************************************
** Program : ncdata.h
** Author  : Neil Massey
** Date    : 02/09/09
** Purpose : class to represent netCDF data and carry out common functions
******************************************************************************/

#ifndef NC_DATA_H
#define NC_DATA_H

#include <string>
#include "rotated_grid.h"
#include "field_data.h"
#include <netcdf>

class ncdata
{
	public:
		ncdata(std::string file_name, std::string var_name);
		~ncdata(void);
		FP_TYPE get_data(int lon_idx, int lat_idx, int z_idx, int t_idx);
		FP_TYPE get_data(FP_TYPE lon, FP_TYPE lat, int z_idx, int t_idx);
		field_data get_field(void);
		field_data get_field(int z_idx, int t_idx);
		
		FP_TYPE get_lon_s(void) { return lon_s; }
		FP_TYPE get_lat_s(void) { return lat_s; }
		FP_TYPE get_lon_d(void) { return lon_d; }
		FP_TYPE get_lat_d(void) { return lat_d; }
		int get_lon_len(void) { return lon_l; }
		int get_lat_len(void) { return lat_l; }
		int get_t_len(void) { return t_len; }
		FP_TYPE get_t_s(void) { return t_s; }
		FP_TYPE get_t_d(void) { return t_d; }
		int get_lon_idx(FP_TYPE lon);
		int get_lat_idx(FP_TYPE lat);
		FP_TYPE get_lon_from_idx(int lon_idx);
		FP_TYPE get_lat_from_idx(int lat_idx);
		void get_grid_details(FP_TYPE& olon_s, FP_TYPE& olat_s, FP_TYPE& olon_d,
							  FP_TYPE& olat_d, int& olon_l, int& olat_l);
		FP_TYPE get_missing_value(void) { return mv; }
		bool has_rotated_grid(void) {return p_rotated_grid != NULL;}
		rotated_grid* get_rotated_grid(void) {return p_rotated_grid;}
		std::string get_units(void);
		std::string get_file_name(void);
		std::string get_var_name(void);
		std::string get_time_dim_name(void);
		void get_reference_time(int& year, int& month, int& day, FP_TYPE& day_scale, FP_TYPE& n_days_py);

	private:
		netCDF::NcFile* nc_file;
		netCDF::NcVar*  nc_var;
		FP_TYPE lon_s, lat_s;       // start of longitude / latitude
		FP_TYPE lon_d, lat_d;       // delta of lon / lat
		int     lon_l, lat_l;        // length of longitude / latitude
		FP_TYPE t_s, t_d;           // time start / delta
		// dimensions for lat / lon / time / level
		int lon_dim, lat_dim;
		int t_dim, t_len, z_dim;
		FP_TYPE mv;                 // missing value
		std::string fname, vname;   // useful for metadata
		FP_TYPE* current_data;      // current field of data
		int c_z, c_t;               // valid z_level for current_data
		rotated_grid* p_rotated_grid;   // rotated grid of variable - if it exists
};

#endif
