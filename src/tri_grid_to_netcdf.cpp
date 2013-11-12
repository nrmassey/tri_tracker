/******************************************************************************
** Program : tri_grid_to_netcdf.cpp
** Author  : Neil Massey
** Date    : 03/07/13
** Purpose : convert the output of the triangular regridding algorithm back 
**           into a netCDF file - mostly for debugging purposes
**           netCDF grid will be created from 90 to -90 latitude, -180 to 180 
**           longitude with spacing defined by the user
******************************************************************************/

#include <tclap/CmdLine.h>
#include <string>
#include <fstream>
#include <netcdfcpp.h>

#include "tri_grid.h"
#include "data_store.h"
#include "geo_convert.h"
#include <math.h>

/*****************************************************************************/

int get_lat_idx(FP_TYPE lat, FP_TYPE lat_s, FP_TYPE lat_d)
{
	return int((lat_s - lat) / lat_d + 0.5);
}

/*****************************************************************************/

int get_lon_idx(FP_TYPE lon, FP_TYPE lon_s, FP_TYPE lon_d)
{
	int idx = int((lon - lon_s) / lon_d + 0.5);
	if (idx < 0)
		idx = int((lon+360 - lon_s) / lon_d + 0.5);
	return idx;
}

/*****************************************************************************/

void save(std::string nc_fname, std::string nc_var_name, FP_TYPE lon_d, 
		  FP_TYPE lat_d, int t_s, FP_TYPE* data, FP_TYPE mv)
{
	// create file
	NcFile out_file(nc_fname.c_str(), NcFile::Replace);
	// check creation
	if (!out_file.is_valid())
		throw(std::string("# Could not save to file: ") + nc_fname);
	// calculate latitude / longitude dimensions
	int lon_s = int(360.0 / lon_d + 0.5);
	int lat_s = int(180 / lat_d + 0.5);
	// add the dimensions
	NcDim* t_dim = out_file.add_dim("t", t_s);
	NcDim* lon_dim = out_file.add_dim("longitude", lon_s);
	NcDim* lat_dim = out_file.add_dim("latitude", lat_s);
	// create the vectors
	FP_TYPE* lon = new FP_TYPE[lon_s];
	FP_TYPE* lat = new FP_TYPE[lat_s];
	int* td = new int[t_s];
	// fill lon
	lon[0] = 0.0;
	for (int l=1; l<lon_s; l++)
		lon[l] = lon[l-1]+lon_d;
	// fill lat
	lat[0] = 90.0;
	for (int l=1; l<lat_s; l++)
		lat[l] = lat[l-1]-lat_d;
	// fill t = just numbers
	for (int t=0; t<t_s; t++)
		td[t] = t;
	// add the dimension variables and write their data
	NcVar* t_var = out_file.add_var("t", ncFloat, t_dim);
	t_var->put(td, t_s);
	NcVar* lon_var = out_file.add_var("longitude", ncFloat, lon_dim);
	lon_var->put(lon, lon_s);
	NcVar* lat_var = out_file.add_var("latitude", ncFloat, lat_dim);
	lat_var->put(lat, lat_s);
	// create the actual variable
	NcVar* data_var = out_file.add_var(nc_var_name.c_str(), ncFloat, t_dim, lat_dim, lon_dim);
	data_var->add_att("missing_value", mv);
	data_var->set_cur(0,0,0);
	data_var->put(data, t_s, lat_s, lon_s);
	// delete the created memory
	delete [] lon;
	delete [] lat;
	delete [] td;
		
	out_file.close();
}
/*****************************************************************************/

void tri_grid_to_netcdf(tri_grid& tg, data_store& ds, int mesh_level, 
						FP_TYPE lon_d, FP_TYPE lat_d, std::string nc_fname, 
						std::string nc_var_name)
{
	// calculate the size of the data and create it
	int lon_s = int(360.0 / lon_d + 0.5);
	int lat_s = int(180 / lat_d + 0.5);
	int t_s = ds.get_number_of_time_steps();
	FP_TYPE* data = new FP_TYPE[t_s*lon_s*lat_s];
	// fill with missing
	for (int i=0; i<t_s*lon_s*lat_s; i++)
		data[i] = ds.get_missing_value();
	// get all the triangles at the level of the mesh we want to output
	std::list<QT_TRI_NODE*> tri_node_list = tg.get_triangles_at_level(mesh_level);
	
	for (std::list<QT_TRI_NODE*>::iterator it = tri_node_list.begin();
		 it != tri_node_list.end(); it++)
	{
		// get the triangle centroid
		vector_3D c = (*it)->get_data()->centroid();
		// convert to lat / lon
		FP_TYPE lon, lat;
		cart_to_model(c, lon, lat);
		// get the triangle datastore index
		int tgt_idx = (*it)->get_data()->get_ds_index();
		int lon_idx = get_lon_idx(lon, 0,  lon_d);
		int lat_idx = get_lat_idx(lat, 90, lat_d);
		// loop through each timestep
		FP_TYPE v;
		for (int t=0; t<t_s; t++)
		{
			// get the value
			v = ds.get_data(t, tgt_idx);
			// calculate output index position
			//            < field       > + < line        > + <cell >
			int out_idx	= t*(lon_s*lat_s) + (lon_s*lat_idx) + lon_idx;
			data[out_idx] = v;
		}
	}
	save(nc_fname, nc_var_name, lon_d, lat_d, t_s, data, ds.get_missing_value());
	delete [] data;
}

/*****************************************************************************/

int main(int argc, char** argv)
{
	// command line arguments
	std::string mesh_fname;		// mesh from gen_grid filename
	std::string regrid_fname;	// output file from regrid filename
	FP_TYPE		lon_d;			// longitude spacing
	FP_TYPE		lat_d;			// latitude spacing
	std::string nc_fname;		// name of netCDF file to write to
	std::string nc_var_name;	// name of variable in netCDF file
	int mesh_level;				// level of the mesh for which to output the data

	std::cout << "#### tri_grid_to_netcdf" << std::endl;

	try
	{
		TCLAP::CmdLine cmd("Output data from triangular mesh regridding to a regularly gridded netCDF file.");
		TCLAP::ValueArg<FP_TYPE> lon_d_arg("s", "lon_d", "Longitude spacing in netCDF output file", false, 1.0, "float");
		cmd.add(lon_d_arg);
		TCLAP::ValueArg<FP_TYPE> lat_d_arg("t", "lat_d", "Latitude spacing in netCDF output file", false, 1.0, "float");
		cmd.add(lat_d_arg);
		TCLAP::ValueArg<std::string> nc_fname_arg("o", "nc_file", "Name of netCDF output file", true, "", "string");
		cmd.add(nc_fname_arg);
		TCLAP::ValueArg<std::string> nc_var_name_arg("v", "nc_var", "Name of variable in netCDF output file", true, "", "string");
		cmd.add(nc_var_name_arg);
		TCLAP::ValueArg<int> mesh_level_arg("l", "mesh_level", "Level of mesh for which to output data", false, 1, "integer");
		cmd.add(mesh_level_arg);
		TCLAP::ValueArg<std::string> mesh_fname_arg("m", "mesh_file", "Name of file containing mesh generated by gen_grid", true, "", "string");
		cmd.add(mesh_fname_arg);
		TCLAP::ValueArg<std::string> regrid_fname_arg("i", "regrid_file", "Name of file containing regridded data generated by regrid", true, "", "string");
		cmd.add(regrid_fname_arg);
		cmd.parse(argc, argv);
		
		// get the values
		lon_d = lon_d_arg.getValue();
		lat_d = lat_d_arg.getValue();
		nc_fname = nc_fname_arg.getValue();
		nc_var_name = nc_var_name_arg.getValue();
		mesh_level = mesh_level_arg.getValue();
		mesh_fname = mesh_fname_arg.getValue();
		regrid_fname = regrid_fname_arg.getValue();
	}
	catch (TCLAP::ArgException &e)  // catch exceptions
	{
		std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	catch (char const* m)           // general error message exception
	{
		std::cerr << "Error: " << m << std::endl;
		return 1;
	}
	catch (std::string &s)          // string error message
	{
		std::cerr << "Error: " << s << std::endl;
		return 1;
	}
	
	// now load the mesh and data into the tri grid and data store respectively
	try
	{
		tri_grid tg;
		tg.load(mesh_fname);
		data_store ds;
		ds.load(regrid_fname);
		tri_grid_to_netcdf(tg, ds, mesh_level, lon_d, lat_d, nc_fname, nc_var_name);
		std::cout << "Saved data to file: " << nc_fname.c_str() << std::endl;
	}
	catch(std::string &s)
	{
		std::cerr << s << std::endl;
		return 1;
	}
	return 0;
}