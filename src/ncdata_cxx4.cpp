/******************************************************************************
** Program : ncdata.cpp
** Author  : Neil Massey
** Date    : 02/09/09
** Purpose : class to represent netCDF data and carry out common functions
** Modified: 24/10/14 - to use new netcdf cpp libraries 4.2
******************************************************************************/

#include "ncdata.h"
#include <exception>
#include <sstream>
#include <vector>

/*****************************************************************************/

int get_dim_pos(NcVar nc_var, std::string search_name)
{
    int dim_id = -1;
    for (int i=0; i<nc_var.getDimCount(); i++)
    {
        NcDim dim = nc_var.getDim(i);
        std::string dim_name = dim.getName();
        if (dim_name.find(search_name) != dim_name.npos)
        {
            dim_id = i;
            break;
        }
    }
    return dim_id;
}

/*****************************************************************************/

void get_dimension_pos(NcVar nc_var, int& t_dim, int& z_dim,
                       int& lon_dim, int& lat_dim)
{
	// time dimension - try "time" and "t"
    t_dim = get_dim_pos(nc_var, "time");
	if (t_dim == -1) t_dim = get_dim_pos(nc_var, "t");
    lon_dim = get_dim_pos(nc_var, "lon");
    lat_dim = get_dim_pos(nc_var, "lat");
    if (lon_dim == -1 || lat_dim == -1)
		throw std::string("netCDF variable is not associated with a latitude or longitude dimension");
    z_dim = get_dim_pos(nc_var, "z");
    // first try "height"
    if (z_dim == -1) z_dim = get_dim_pos(nc_var, "height");
	// next try "p"
    if (z_dim == -1) z_dim = get_dim_pos(nc_var, "p");
	// next try "surface"
	if (z_dim == -1) z_dim = get_dim_pos(nc_var, "surface");
	// next try "lev" (ERA_Interim)
	if (z_dim == -1) z_dim = get_dim_pos(nc_var, "lev");
}

/*****************************************************************************/

void get_grid_spacing(NcFile nc_file, NcVar nc_var,
                      FP_TYPE& lon_s, FP_TYPE& lat_s,
                      FP_TYPE& lon_d, FP_TYPE& lat_d,
                      int&    lon_l, int&    lat_l)
{
    // get the start, spacing and length of the longitude dimension
    int lon_id = get_dim_pos(nc_var, "lon");
    int lat_id = get_dim_pos(nc_var, "lat");
    // get the variables of the lon / lat
    NcVar lon_var = nc_file.getVar(nc_var.getDim(lon_id).getName());
    NcVar lat_var = nc_file.getVar(nc_var.getDim(lat_id).getName());

    // get the start and deltas
    FP_TYPE* lon_vals;
    FP_TYPE* lat_vals;
    lon_var.getVar(lon_vals);
    lat_var.getVar(lat_vals);
    lon_s = lon_vals[0];
    lat_s = lat_vals[0];
    lon_d = lon_vals[1] - lon_s;
    lat_d = lat_vals[1] - lat_s;

    // get the size from the dimension
    lon_l = nc_var.getDim(lon_id).getSize();
    lat_l = nc_var.getDim(lat_id).getSize();
}

/*****************************************************************************/

ncdata::ncdata(std::string file_name, std::string var_name) :
			   current_data(NULL), c_z(-1), c_t(-1), p_rotated_grid(NULL)
{
	fname = file_name;
	vname = var_name;
	nc_file = NcFile(file_name.c_str(), NcFile::read);
	std::cout << "# Loading netCDF file " << file_name << std::endl;
	// get the netCDF variable
	nc_var = nc_file.getVar(var_name.c_str());
	// get the grid spacing
	get_grid_spacing(nc_file, nc_var, lon_s, lat_s, lon_d, lat_d, lon_l, lat_l);
	// get the indexes
	get_dimension_pos(nc_var, t_dim, z_dim, lon_dim, lat_dim);
	// if no time dimension then set length to 1
	if (t_dim == -1)
		t_len = 1;
	else
		// otherwise get the size
		t_len = nc_var.getDim(t_dim).getSize();

	current_data = new FP_TYPE [lon_l * lat_l];
	mv = 2e20;
	// get the missing value
	try
	{
		FP_TYPE* att_vars;
		NcVarAtt fv_att = nc_var.getAtt("_FillValue");
		fv_att.getValues(att_vars);
		mv = att_vars[0];
	}
	catch (...) 
	{
	}
	
	// see if this variable has a grid mapping attribute
	try
	{
		NcVarAtt gm_att = nc_var.getAtt("grid_mapping");
		// get the name of the variable with the grid mapping and then get the variable
		std::string gm_name;
		gm_att.getValues(gm_name);
		NcVar gm_var = nc_file.getVar(gm_name);
		// in this variable there are the rotated pole lat and lon
		NcVarAtt gmv_att_lat = gm_var.getAtt("grid_north_pole_latitude");
		NcVarAtt gmv_att_lon = gm_var.getAtt("grid_north_pole_longitude");
		float* rotated_pole_lat_vals;
		gmv_att_lat.getValues(rotated_pole_lat_vals);
		float* rotated_pole_lon_vals;
		gmv_att_lon.getValues(rotated_pole_lon_vals);
		float rotated_pole_lat = rotated_pole_lat_vals[0];
		float rotated_pole_lon = rotated_pole_lon_vals[0];
		// create the rotated grid
		p_rotated_grid = new rotated_grid(rotated_pole_lat, rotated_pole_lon,
										  this);
	}
	catch (...)
	{
	}
}

/*****************************************************************************/

ncdata::~ncdata(void)
{
	delete [] current_data;
	delete p_rotated_grid;
}

/*****************************************************************************/

float ncdata::get_data(int lon_idx, int lat_idx, int z_idx, int t_idx)
{
	float value;
	if (t_dim != -1 && z_dim != -1)	// 4D data
	{
		// check whether the data has been cached or not
		if (c_z != z_idx || c_t != t_idx)
		{
			// get the full field of data
			std::vector<size_t> s(4);
			std::vector<size_t> c(4);
			s[z_dim] = z_idx;
			s[t_dim] = t_idx;
			c[lon_dim] = lon_l;
			c[lat_dim] = lat_l;
			nc_var.getVar(s,c,current_data);
			c_z = z_idx;
			c_t = t_idx;
		}
		// now calculate the index into the field
		value = current_data[lat_idx*lon_l + lon_idx];
	}
	else if (t_dim != -1)			// 3D data with time component
	{
		if (c_t != t_idx)
		{
			std::vector<size_t> s(3);
			std::vector<size_t> c(3);
			s[t_dim] = t_idx;
			c[lon_dim] = lon_l;
			c[lat_dim] = lat_l;
			nc_var.getVar(s,c,current_data);
			c_t = t_idx;
		}
		value = current_data[lat_idx*lon_l + lon_idx];
	}
	else if (z_dim != -1)			// 3D data with z component
	{
		if (c_z != z_idx)
		{
			std::vector<size_t> s(3);
			std::vector<size_t> c(3);
			s[z_dim] = z_idx;
			c[lon_dim] = lon_l;
			c[lat_dim] = lat_l;
			nc_var.getVar(s,c,current_data);
			c_z = z_idx;
		}
		value = current_data[lat_idx*lon_l + lon_idx];
	}
	else							// 2D data
	{
		if (c_t == -1)
		{
			std::vector<size_t> s(2);
			std::vector<size_t> c(2);
			s[lon_dim] = lon_idx;
			s[lat_dim] = lat_idx;
			nc_var.getVar(s,c,current_data);
			c_t = 0;
		}
		value = current_data[lat_idx*lon_l + lon_idx];
	}

	return value;
}

/*****************************************************************************/

float ncdata::get_data(FP_TYPE lon, FP_TYPE lat, int z_idx, int t_idx)
{
	// calculate the lon and lat indices from the start and spacing
	int lon_idx = int((lon - lon_s) / lon_d + 0.5);
	int lat_idx = int((lat - lat_s) / lat_d + 0.5);

	return get_data(lon_idx, lat_idx, z_idx, t_idx);
}

/*****************************************************************************/

int ncdata::get_lon_idx(FP_TYPE lon)
{
	int idx = int((lon - lon_s) / lon_d + 0.5);
	if (idx < 0)
		idx = int((lon+360 - lon_s) / lon_d + 0.5);
	return idx;
}

/*****************************************************************************/

int ncdata::get_lat_idx(FP_TYPE lat)
{
	return int((lat - lat_s) / lat_d + 0.5);
}

/*****************************************************************************/

FP_TYPE ncdata::get_lon_from_idx(int lon_idx)
{
	
	return lon_s + lon_idx * lon_d;
}

/*****************************************************************************/

FP_TYPE ncdata::get_lat_from_idx(int lat_idx)
{
	return lat_s + lat_idx * lat_d;
}

/*****************************************************************************/

void ncdata::get_grid_details(FP_TYPE& olon_s, FP_TYPE& olat_s, FP_TYPE& olon_d, 
					  		  FP_TYPE& olat_d, int& olon_l, int& olat_l)
{
	olon_s = lon_s;
	olat_s = lat_s;
	olon_d = lon_d;
	olat_d = lat_d;
	olon_l = lon_l;
	olat_l = lat_l;
}

/*****************************************************************************/

std::string ncdata::get_units(void)
{
	// get a string representing the units
	std::string units = "none";
	try
	{
		NcVarAtt u_att = nc_var.getAtt("units");
		u_att.getValues(units);
	}
	catch (...)
	{
	}
	return units;
}

/*****************************************************************************/

std::string ncdata::get_file_name(void)
{
	return fname;
}

/*****************************************************************************/

std::string ncdata::get_var_name(void)
{
	return vname;
}

/*****************************************************************************/

std::string ncdata::get_time_dim_name(void)
{
	return nc_var.getDim(t_dim).getName();
}

/*****************************************************************************/

void ncdata::get_reference_time(int& year, int& month, int& day, FP_TYPE& day_sc, FP_TYPE& n_days_py)
{
    // get the reference time from the netcdf file and variable
    try
    {
		// get the time variable
		NcVar time_var = nc_file.getVar(nc_var.getDim(t_dim).getName());

		NcVarAtt tu_att = time_var.getAtt("units");
		char dummy;
		std::string scale;
		std::string dummy_string;
		std::string time_string;
		std::string time_units;
		tu_att.getValues(time_units);
		// parse the string - format is <scale> since <year>-<month>-<day> 00:00:00
		std::stringstream stream(time_units);
		stream >> scale >> dummy_string >> time_string;
		std::stringstream (time_string.substr(0,4)) >> year;
		std::stringstream (time_string.substr(5,2)) >> month;
		std::stringstream (time_string.substr(8,2)) >> day;
		// if the scale is "days since" - day_sc is 1, if it's "hours since" - day_sc is 1.0/24
		if (scale == "days")
			day_sc = 1.0;
		if (scale == "hours")
			day_sc = 1.0/24.0;		
		// get the number of days per year
		NcVarAtt dpy_att = time_var.getAtt("calendar");
		std::string cal_type;
		dpy_att.getValues(cal_type);
		if (cal_type == "standard")
			n_days_py = 365.25;
		if (cal_type == "360_day")
			n_days_py = 360.0;
	}
	catch (...)
	{
	}
}
