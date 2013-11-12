/******************************************************************************
** Program : ncdata.cpp
** Author  : Neil Massey
** Date    : 02/09/09
** Purpose : class to represent netCDF data and carry out common functions
******************************************************************************/

#include "ncdata.h"
#include <netcdfcpp.h>
#include <exception>

/*****************************************************************************/

int get_dim_pos(NcVar* nc_var, std::string search_name)
{
    int dim_id = -1;
    for (int i=0; i<nc_var->num_dims(); i++)
    {
        NcDim* dim = nc_var->get_dim(i);
        std::string dim_name = dim->name();
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
}

/*****************************************************************************/

void get_grid_spacing(NcFile* nc_file, NcVar* nc_var,
                      FP_TYPE& lon_s, FP_TYPE& lat_s,
                      FP_TYPE& lon_d, FP_TYPE& lat_d,
                      int&    lon_l, int&    lat_l)
{
    // get the start, spacing and length of the longitude dimension
    int lon_id = get_dim_pos(nc_var, "lon");
    int lat_id = get_dim_pos(nc_var, "lat");
    // get the variables of the lon / lat
    NcVar* lon_var = nc_file->get_var(nc_var->get_dim(lon_id)->name());
    NcVar* lat_var = nc_file->get_var(nc_var->get_dim(lat_id)->name());

    // get the start and deltas
    lon_s = lon_var->as_float(0);
    lat_s = lat_var->as_float(0);
    lon_d = lon_var->as_float(1) - lon_s;
    lat_d = lat_var->as_float(1) - lat_s;

    // get the size from the dimension
    lon_l = nc_var->get_dim(lon_id)->size();
    lat_l = nc_var->get_dim(lat_id)->size();
}

/*****************************************************************************/

ncdata::ncdata(std::string file_name, std::string var_name) :
			   nc_file(NULL), nc_var(NULL), current_data(NULL), 
			   c_z(-1), c_t(-1), p_rotated_grid(NULL)
{
	fname = file_name;
	vname = var_name;
	nc_file = new NcFile(file_name.c_str());
	if (!nc_file->is_valid())
	{
		throw std::string("file " + file_name + " not found or not a netCDF file.");
	}
	std::cout << "# Loading netCDF file" << std::endl;
	// get the netCDF variable
	nc_var = nc_file->get_var(var_name.c_str());
	// get the grid spacing
	get_grid_spacing(nc_file, nc_var, lon_s, lat_s, lon_d, lat_d, lon_l, lat_l);
	// get the indexes
	get_dimension_pos(nc_var, t_dim, z_dim, lon_dim, lat_dim);
	// if no time dimension then set length to 1
	if (t_dim == -1)
		t_len = 1;
	else
		// otherwise get the size
		t_len = nc_var->get_dim(t_dim)->size();

	current_data = new FP_TYPE [lon_l * lat_l];
	mv = 2e20;
	// get the missing value
	{
		NcError error(NcError::silent_nonfatal);
		NcAtt* fv_att = nc_var->get_att("_FillValue");
		if (fv_att != 0)
			mv = fv_att->as_float(0);
		delete fv_att;
	}
	
	// see if this variable has a grid mapping attribute
	{
		NcError error(NcError::silent_nonfatal);	// don't die if the attribute not found
		NcAtt* gm_att = nc_var->get_att("grid_mapping");
		if (gm_att != 0)
		{
			// get the name of the variable with the grid mapping and then get the variable
			char* gm_name = gm_att->as_string(0);
			NcVar* gm_var = nc_file->get_var(gm_name);
			// in this variable there are the rotated pole lat and lon
			NcAtt* gmv_att_lat = gm_var->get_att("grid_north_pole_latitude");
			NcAtt* gmv_att_lon = gm_var->get_att("grid_north_pole_longitude");
			float rotated_pole_lat = gmv_att_lat->as_float(0);
			float rotated_pole_lon = gmv_att_lon->as_float(0);
			// create the rotated grid
			p_rotated_grid = new rotated_grid(rotated_pole_lat, rotated_pole_lon,
											  this);
			// delete attribute objects
			delete gm_att;
			delete gmv_att_lat;
			delete gmv_att_lon;
		}
	}
}

/*****************************************************************************/

ncdata::~ncdata(void)
{
	delete [] current_data;
	delete nc_file;
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
			long s[4] = {0,0,0,0};
			long c[4] = {1,1,1,1};
			s[z_dim] = z_idx;
			s[t_dim] = t_idx;
			c[lon_dim] = lon_l;
			c[lat_dim] = lat_l;
			nc_var->set_cur(s);
			nc_var->get(current_data, c);
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
			long s[3] = {0,0,0};
			long c[3] = {1,1,1};
			s[t_dim] = t_idx;
			c[lon_dim] = lon_l;
			c[lat_dim] = lat_l;
			nc_var->set_cur(s);
			nc_var->get(current_data, c);
			c_t = t_idx;
		}
		value = current_data[lat_idx*lon_l + lon_idx];
	}
	else if (z_dim != -1)			// 3D data with z component
	{
		if (c_z != z_idx)
		{
			long s[3] = {0,0,0};
			long c[3] = {1,1,1};
			s[z_dim] = z_idx;
			c[lon_dim] = lon_l;
			c[lat_dim] = lat_l;
			nc_var->get(current_data, c);
			nc_var->set_cur(s);
			c_z = z_idx;
		}
		value = current_data[lat_idx*lon_l + lon_idx];
	}
	else							// 2D data
	{
		if (c_t == -1)
		{
			long s[2] = {0,0};
			long c[2] = {1,1};
			s[lon_dim] = lon_idx;
			s[lat_dim] = lat_idx;
			nc_var->set_cur(s);
			nc_var->get(current_data, c);
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
	int i_n_atts = nc_var->num_atts();
	for (int i=0; i<i_n_atts; i++)
	{
		NcAtt* nc_att = nc_var->get_att(i);
		if (std::string(nc_att->name()) == "units")
		{
			units = std::string(nc_att->as_string(0));
			break;
		}
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
