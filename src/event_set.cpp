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

int main(int argc, char** argv)
{
	std::cout << "#### event set." << std::endl;
	std::vector<std::string> input_fname;
	std::string output_fname;
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
		TCLAP::ValueArg<FP_TYPE> overlap_arg("v", "overlap", "Minimum percentage of overlapping triangles between objects in adjacent frames", false, 10.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<int> event_t_steps_arg("t", "event_timesteps", "Number of timesteps to use in integrated event set", false, 12, "integer", cmd);
		TCLAP::ValueArg<FP_TYPE> sr_arg("r", "search_radius", "Radius, in km, to search when locating tracks", false, 500.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<FP_TYPE> min_tl_arg("l", "min_length", "Minimum track length in km", false, 0.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<FP_TYPE> min_dev_arg("d", "min_dev", "Minimum track deviation from first to last point in km", false, 0.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<int> min_per_arg("p", "min_per", "Minimum persistence of the track in frames / timesteps", false, 1, "integer", cmd);
		TCLAP::ValueArg<std::string> output_fname_arg("o", "output", "Output file name", true, "", "string", cmd);
		TCLAP::MultiArg<std::string> input_fname_arg("i", "input", "Input file name", true, "string", cmd);
		cmd.parse(argc, argv);

		overlap = overlap_arg.getValue();
		event_t_steps = event_t_steps_arg.getValue();
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
		eventor EV(input_fname, min_per, min_tl, min_dev, search_rad, event_t_steps, overlap);
		EV.find_events();
	}
	catch(std::string &s)
	{
		std::cerr << s << std::endl;
	}

}