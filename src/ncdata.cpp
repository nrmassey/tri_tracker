/******************************************************************************
** Program : ncdata.cpp
** Author  : Neil Massey
** Date    : 02/09/09
** Purpose : class to represent netCDF data and carry out common functions
******************************************************************************/

#include "ncdata.h"
#include <netcdf>
#include <exception>
#include <sstream>
#include <map>

using namespace netCDF;

/*****************************************************************************/

int get_dim_pos(NcVar* nc_var, std::string search_name)
{
    int dim_id = -1;
    for (int i=0; i<nc_var->getDimCount(); i++)
    {
        NcDim dim = nc_var->getDim(i);
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

void get_dimension_pos(NcVar* nc_var, int& t_dim, int& z_dim,
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

void get_grid_spacing(NcFile* nc_file, NcVar* nc_var,
                      FP_TYPE& lon_s, FP_TYPE& lat_s,
                      FP_TYPE& lon_d, FP_TYPE& lat_d,
                      int&     lon_l, int&     lat_l,
                      FP_TYPE& t_s,   FP_TYPE& t_d)
{
    // get the start, spacing and length of the longitude dimension
    int lon_id, lat_id, z_id, t_id;
    get_dimension_pos(nc_var, t_id, z_id, lon_id, lat_id);
    // get the variables of the lon / lat
    NcVar lon_var = nc_file->getVar(nc_var->getDim(lon_id).getName());
    NcVar lat_var = nc_file->getVar(nc_var->getDim(lat_id).getName());

    // get the start and deltas - this very different in new netCDF cpp4
    // create a single dimension vector with an index to the very first element in the array
    std::vector<size_t> idx(1);
    idx[0] = 0;
    lon_var.getVar(idx, &lon_s);
    lat_var.getVar(idx, &lat_s);
    // iterate the index to point to the next latitude value
    idx[0] += 1;
    lon_var.getVar(idx, &lon_d);
    lon_d -= lon_s;					// get the difference between the 1st and 0th index
    lat_var.getVar(idx, &lat_d);
    lat_d -= lat_s;

    // get the size from the dimension
    lon_l = nc_var->getDim(lon_id).getSize();
    lat_l = nc_var->getDim(lat_id).getSize();
    
    // get the time start and delta - first get the variable
    NcVar t_var = nc_file->getVar(nc_var->getDim(t_id).getName());
    idx[0] = 0;
    t_var.getVar(idx, &t_s);
    // get the spacing if time dimension length is > 1
    if (t_var.getDim(0).getSize() > 1)
    {
	    idx[0] += 1;
    	t_var.getVar(idx, &t_d);
	    t_d -= t_s;
	}
	else
	{
		t_d = 0.0;
	}
}

/*****************************************************************************/

ncdata::ncdata(std::string file_name, std::string var_name) :
			   nc_file(NULL), nc_var(NULL), current_data(NULL), 
			   c_z(-1), c_t(-1), p_rotated_grid(NULL)
{
	fname = file_name;
	vname = var_name;
	try
	{
		nc_file = new NcFile(file_name, NcFile::read);
	}
	catch(exceptions::NcException& e)
	{
	    e.what();
		throw std::string("file " + file_name + " not found or not a netCDF file.");
	}
		
	std::cout << "# Loading netCDF file " << file_name << std::endl;
	// get the netCDF variable
	try
	{
		nc_var = new NcVar(nc_file->getVar(var_name));
	}
	catch(exceptions::NcException& e)
	{
		e.what();
		throw std::string("var " + var_name + " not found in file " + file_name);
	}
	// get the grid spacing
	get_grid_spacing(nc_file, nc_var, lon_s, lat_s, lon_d, lat_d, lon_l, lat_l, t_s, t_d);
	// get the indexes
	get_dimension_pos(nc_var, t_dim, z_dim, lon_dim, lat_dim);
	// if no time dimension then set length to 1
	if (t_dim == -1)
		t_len = 1;
	else
		// otherwise get the size
		t_len = nc_var->getDim(t_dim).getSize();

	current_data = new FP_TYPE [lon_l * lat_l];
	// get the missing value
	try
	{
		NcVarAtt fv_att = nc_var->getAtt("_FillValue");
		fv_att.getValues(&mv);
	}
	catch(exceptions::NcException& e)
	{
		mv = 2e20;
	}
	
	// see if this variable has a grid mapping attribute
	try
	{
		NcVarAtt gm_att = nc_var->getAtt("grid_mapping");
		// get the name of the variable with the grid mapping and then get the variable
		std::string gm_name;
		gm_att.getValues(gm_name);
		NcVar gm_var = nc_file->getVar(gm_name);
		
		// in this variable there are the rotated pole lat and lon
		NcVarAtt gmv_att_lat = gm_var.getAtt("grid_north_pole_latitude");
		NcVarAtt gmv_att_lon = gm_var.getAtt("grid_north_pole_longitude");
		FP_TYPE rotated_pole_lat, rotated_pole_lon;
		gmv_att_lat.getValues(&rotated_pole_lat);
		gmv_att_lon.getValues(&rotated_pole_lon);
		// create the rotated grid
		p_rotated_grid = new rotated_grid(rotated_pole_lat, rotated_pole_lon,
										  this);
	}
	catch(exceptions::NcException)
	{
	    p_rotated_grid = NULL;
	}
}

/*****************************************************************************/

ncdata::~ncdata(void)
{
	delete [] current_data;
	delete p_rotated_grid;
	delete nc_file;
	delete nc_var;
}

/*****************************************************************************/

FP_TYPE ncdata::get_data(int lon_idx, int lat_idx, int z_idx, int t_idx)
{
	FP_TYPE value;
	if (t_dim != -1 && z_dim != -1)	// 4D data
	{
		// check whether the data has been cached or not
		if (c_z != z_idx || c_t != t_idx)
		{
			// get the full field of data - set the position and counts
			std::vector<size_t> s(4);
			std::vector<size_t> c(4);
			s[lon_dim] = 0;
			s[lat_dim] = 0;
			s[z_dim] = z_idx;
			s[t_dim] = t_idx;
			
			c[lon_dim] = lon_l;
			c[lat_dim] = lat_l;
			c[z_dim] = 1;
			c[t_dim] = 1;
			
			// get the data
			nc_var->getVar(s, c, current_data);
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

			// get the full field of data - set the position and counts
			std::vector<size_t> s(3);
			std::vector<size_t> c(3);
			s[lon_dim] = 0;
			s[lat_dim] = 0;
			s[t_dim] = t_idx;
			
			c[lon_dim] = lon_l;
			c[lat_dim] = lat_l;
			c[t_dim] = 1;
			
			// get the data
			nc_var->getVar(s, c, current_data);
			c_t = t_idx;
		}
		value = current_data[lat_idx*lon_l + lon_idx];
	}
	else if (z_dim != -1)			// 3D data with z component
	{
		if (c_z != z_idx)
		{
			// get the full field of data - set the position and counts
			std::vector<size_t> s(3);
			std::vector<size_t> c(3);
			s[lon_dim] = 0;
			s[lat_dim] = 0;
			s[z_dim] = z_idx;
			
			c[lon_dim] = lon_l;
			c[lat_dim] = lat_l;
			c[t_dim] = 1;
			
			// get the data
			nc_var->getVar(s, c, current_data);
			c_z = z_idx;
		}
		value = current_data[lat_idx*lon_l + lon_idx];
	}
	else							// 2D data
	{
		if (c_t == -1)
		{
			// get the full field of data - set the position and counts
			std::vector<size_t> s(2);
			std::vector<size_t> c(2);
			s[lon_dim] = 0;
			s[lat_dim] = 0;
			
			c[lon_dim] = lon_l;
			c[lat_dim] = lat_l;
			
			// get the data
			nc_var->getVar(s, c, current_data);
			c_t = 0;
		}
		value = current_data[lat_idx*lon_l + lon_idx];
	}

	return value;
}

/*****************************************************************************/

FP_TYPE ncdata::get_data(FP_TYPE lon, FP_TYPE lat, int z_idx, int t_idx)
{
	// calculate the lon and lat indices from the start and spacing
	int lon_idx = int((lon - lon_s) / lon_d + 0.5);
	int lat_idx = int((lat - lat_s) / lat_d + 0.5);

	return get_data(lon_idx, lat_idx, z_idx, t_idx);
}

/*****************************************************************************/

field_data ncdata::get_field(int z_idx, int t_idx)
{
    field_data out_data(lat_l, lon_l, 0.0);
    FP_TYPE* px_out_data = out_data.get();
    
	std::vector<size_t> s(4);
	std::vector<size_t> c(4);
	s[lon_dim] = 0;
	s[lat_dim] = 0;
	s[z_dim] = z_idx;
	s[t_dim] = t_idx;
	
	c[lon_dim] = lon_l;
	c[lat_dim] = lat_l;
	c[z_dim] = 1;
	c[t_dim] = 1;
    
    nc_var->getVar(s, c, px_out_data);
    return out_data;
}

/*****************************************************************************/

field_data ncdata::get_field(void)
{
    field_data out_data(lat_l, lon_l, 0.0);
    FP_TYPE* px_out_data = out_data.get();
	// get all the data
	nc_var->getVar(px_out_data);
    return out_data;
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

std::string get_var_units(NcVar* nc_var)
{
	// get a string representing the units
	std::string units = "none";
	std::map<std::string, NcVarAtt> all_atts = nc_var->getAtts();
	for (std::map<std::string, NcVarAtt>::iterator i = all_atts.begin();
	     i != all_atts.end(); i++)
	{
		if ((*i).first == "units")
		{
 			(*i).second.getValues(units);
 			break;
		}
	}
	return units;
}

/*****************************************************************************/

std::string ncdata::get_units(void)
{
	return get_var_units(nc_var);
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
	return nc_var->getDim(t_dim).getName();
}

/*****************************************************************************/

void ncdata::get_reference_time(int& year, int& month, int& day, FP_TYPE& day_sc, FP_TYPE& n_days_py)
{
    // get the reference time from the netcdf file and variable
    try
    {
		NcVar time_var = nc_file->getVar(nc_var->getDim(t_dim).getName());
		std::string time_units = get_var_units(&time_var);
		
		char dummy;
		std::string scale;
		std::string dummy_string;
		std::string time_string;
		if (time_units != "none")
		{
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
		}
		// get the number of days per year
		NcVarAtt dpy_att = time_var.getAtt("calendar");
		std::string cal_type;
		dpy_att.getValues(cal_type);
		if (cal_type == "standard")
			n_days_py = 365.25;
		if (cal_type == "360_day")
			n_days_py = 360.0;
	}
	catch(exceptions::NcException& e)
	{
	    e.what();
		throw std::string("Could not find time variable, or time variable is ill-formatted");
	}
}
