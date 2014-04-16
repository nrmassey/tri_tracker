/******************************************************************************
** Program : event_set.cpp
** Author  : Neil Massey
** Date    : 11/04/14
** Purpose : program to read in regridded data on a triangular mesh and 
**           construct an event set based on that data
******************************************************************************/

#include <tclap/CmdLine.h>
#include <string>
#include <iostream>
#include "eventor.h"
#include "read_from_string.h"

int main(int argc, char** argv)
{
	std::cout << "#### event set." << std::endl;
	std::vector<std::string> input_fname;
	std::string output_fname;
	std::vector<std::string> mslp_fname;
	std::vector<std::string> wind_fname;
	std::string mslp_field;
	std::string wind_field;
	std::string lsm_file;
	std::string mesh_file;
	FP_TYPE min_per;		// minimum persistence of a track (frames)
	FP_TYPE min_tl;			// minimum track length (km)
	FP_TYPE search_rad;		// search radius (km)
	FP_TYPE min_dev;		// minimum deviation
	FP_TYPE overlap;		// minimum percentage of overlapping objects
	int event_t_steps;		// number of timesteps in integrated event 
							// (e.g. 6 hourly data, 72 hour event - event_t_steps = 12)
							
	try
	{
		TCLAP::CmdLine cmd("Build an event set using extrema located from the extrema program");
		TCLAP::MultiArg<std::string> mslp_fname_arg("m", "mslp_file", "Original MSLP files used in tracking and field name. MSLP values for the event set will be taken from this file, rather than from the extrema file", true, "string", cmd);		
		TCLAP::ValueArg<std::string> mslp_field_arg("n", "mslp_field", "Field name of MSLP in original MSLP file.", false, "field8", "string", cmd);
		TCLAP::MultiArg<std::string> wind_fname_arg("w", "wind_file", "Original wind files used in tracking and field name. 10m wind values for the event set will be taken from this file, rather than from the extrema file", true, "string", cmd);
		TCLAP::ValueArg<std::string> wind_field_arg("x", "wind_field", "Field name of 10m wind variable in original wind file.", false, "field50", "string", cmd);
		TCLAP::ValueArg<std::string> lsm_file_arg("s", "lsm_file", "LSM mask for the model.  So that maximum winds can be computed over land only", false, "", "string", cmd);
		TCLAP::ValueArg<FP_TYPE> overlap_arg("v", "overlap", "Minimum percentage of overlapping triangles between objects in adjacent frames", false, 10.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<int> event_t_steps_arg("t", "event_timesteps", "Number of timesteps to use in integrated event set", false, 12, "integer", cmd);
		TCLAP::ValueArg<FP_TYPE> sr_arg("r", "search_radius", "Radius, in km, to search when locating tracks", false, 500.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<FP_TYPE> min_tl_arg("l", "min_length", "Minimum track length in km", false, 0.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<FP_TYPE> min_dev_arg("d", "min_dev", "Minimum track deviation from first to last point in km", false, 0.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<int> min_per_arg("p", "min_per", "Minimum persistence of the track in frames / timesteps", false, 1, "integer", cmd);
		TCLAP::ValueArg<std::string> mesh_file_arg("g", "mesh_file", "Mesh used in generation of extrema.", true, "", "string", cmd);
		TCLAP::ValueArg<std::string> output_fname_arg("o", "output", "Output file name", true, "", "string", cmd);
		TCLAP::MultiArg<std::string> input_fname_arg("i", "input", "Input file names - extrema files created by ./extrema", true, "string", cmd);
		cmd.parse(argc, argv);

		overlap = overlap_arg.getValue();
		event_t_steps = event_t_steps_arg.getValue();
		mslp_fname = mslp_fname_arg.getValue();
		mslp_field = mslp_field_arg.getValue();
		wind_fname = wind_fname_arg.getValue();
		wind_field = wind_field_arg.getValue();
		lsm_file = lsm_file_arg.getValue();
		mesh_file = mesh_file_arg.getValue();
		search_rad = sr_arg.getValue() * 1000.0;	// convert to metres
		min_tl = min_tl_arg.getValue() * 1000.0;	// convert to metres
		min_dev = min_dev_arg.getValue() * 1000.0;
		min_per = min_per_arg.getValue();
		output_fname = output_fname_arg.getValue();
		input_fname = input_fname_arg.getValue();
		
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
	
	try
	{
		eventor EV(input_fname, min_per, min_tl, min_dev, search_rad, event_t_steps, overlap,
				   mslp_fname, mslp_field, wind_fname, wind_field, lsm_file, mesh_file);
		EV.find_events();
		EV.save(output_fname);
	}
	catch(std::string &s)
	{
		std::cerr << s << std::endl;
	}

}