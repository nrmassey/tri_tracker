/******************************************************************************
** Program : gen_grid.cpp
** Author  : Neil Massey
** Date    : 15/03/13
** Purpose : program to generate and output triangular, equal area, equally
**           spaced and hierarchical mesh, using a regularly spaced grid as
**           an input.  This includes regional grids and may result in a non-
**           uniform mesh being output, including triangles which have not been
**           split from the original icosahedron.
******************************************************************************/

#include <iostream>
#include <fstream>
#include <math.h>
#include <tclap/CmdLine.h>
#include <string>

#include "ncdata.h"
#include "tri_grid.h"

int main(int argc, char** argv)
{
	// set up the command line arguments
	int max_I, init_S, max_pts, max_levs;
	std::string out_name = "";
	std::string in_name = "";
	std::string var_name = "";
	bool stats = false;
	bool text_out = false;
	bool save_text = false;
	std::cout << "#### gen_grid" << std::endl;

	try
	{
		TCLAP::CmdLine cmd("Generate an equally spaced, equal area, potentially irregular, triangular grid mesh", ' ', "0.2");
		// define the arguments - backwards
		TCLAP::SwitchArg text_out_arg("T", "text", "Output grid in text format, as well as the binary format", cmd, false);
		TCLAP::SwitchArg lat_lon_out_arg("L", "lat_lon", "Output grid in text latitude / longitude coordinates, as well as the binary format", cmd, false);
		TCLAP::SwitchArg stats_arg("S", "stats", "Show statistics for the mesh", cmd, false);
		TCLAP::ValueArg<int> init_shp_arg("i", "init_shape", "Initialisation shape: 0 = ICOSAHEDRON, 1 = OCTAHEDRON, 2 = DYMAXION(tm)", false, 0, "int", cmd);
		TCLAP::ValueArg<int> max_its_arg("I", "max_its", "Maximum number of iterations in equalisation phase", false, 1e3, "int", cmd);
		TCLAP::ValueArg<int> max_pts_arg("n", "max_pts", "Maximum number of points from the original grid allowed per triangle", false, 0, "int", cmd);
		TCLAP::ValueArg<int> max_levs_arg("l", "max_level", "Maximum number of levels allowed in the grid", false, -1, "int", cmd);
		TCLAP::ValueArg<std::string> out_name_arg("o", "output", "Grid output file name", true, "", "string", cmd);
		TCLAP::ValueArg<std::string> var_name_arg("v", "var_name", "Variable name in source grid file", true, "", "string", cmd);
		TCLAP::ValueArg<std::string> in_name_arg("f", "grid_input", "Input of source grid file name", true, "", "string", cmd);

		// Parse the command line arguments
		cmd.parse(argc, argv);

		// get the values
		text_out = text_out_arg.getValue();
		save_text = lat_lon_out_arg.getValue();
		stats = stats_arg.getValue();
		out_name = out_name_arg.getValue();
		var_name = var_name_arg.getValue();
		in_name = in_name_arg.getValue();
		max_I = max_its_arg.getValue();
		init_S = init_shp_arg.getValue();
		max_pts = max_pts_arg.getValue();
		max_levs = max_levs_arg.getValue();		
	}
	catch (TCLAP::ArgException &e)  // catch exceptions
	{
		std::cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
		return 1;
	}
	catch (const char* m)           // general error message exception
	{
		std::cerr << "Error: " << m << std::endl;
		return 1;
	}
	catch (...)
	{
		return 1;
	}
	
	try
	{
		ncdata nc_data(in_name, var_name);
		tri_grid tg;
		tg.initialize(SHAPE(init_S), &nc_data, max_pts, max_levs, max_I);
		tg.save(out_name);
		std::cout << "# Saved mesh to file: " << out_name.c_str() << std::endl;
		if (save_text)
			tg.save_text(out_name+"lat_lon.txt");
	}
	catch (std::string s)
	{
		std::cerr << "# Error: " << s.c_str() << std::endl;
		return 1;
	}
	return 0;
}
