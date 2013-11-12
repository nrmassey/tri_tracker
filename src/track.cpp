/******************************************************************************
** Program : track.cpp
** Author  : Neil Massey
** Date    : 05/08/09
** Purpose : program to read in regridded data on a triangular mesh and detect
**           extrema within that data
******************************************************************************/

#include <tclap/CmdLine.h>
#include <string>
#include <iostream>
#include "tracker.h"

int main(int argc, char** argv)
{
	std::cout << "#### track" << std::endl;
	std::vector<std::string> input_fname;
	std::string output_fname;
	FP_TYPE min_per;		// minimum persistence of a track (frames)
	FP_TYPE min_tl;		// minimum track length (km)
	FP_TYPE search_rad;	// search radius (km)
	FP_TYPE min_dev;		// minimum deviation
	int t,s,k;
	FP_TYPE w0, w1, w2, w3;
	bool text_out = false;

	try
	{
		TCLAP::CmdLine cmd("Track extrema located using extrema program");

		TCLAP::ValueArg<int> k_arg("k", "k", "Number of timesteps in file to skip each tracking iteration", false, 1, "integer", cmd);
		TCLAP::ValueArg<int> s_arg("s", "s", "Number of points to interpolate between identified points in final track", false, 0, "integer", cmd);
		TCLAP::ValueArg<int> t_arg("t", "t", "Number of permissable timesteps between track points", false, 1, "integer", cmd);
		TCLAP::ValueArg<FP_TYPE> w0_arg("0", "w0", "Weight for distance in cost function", false, 1e6, "FP_TYPE", cmd);
		TCLAP::ValueArg<FP_TYPE> w1_arg("1", "w1", "Weight for intensity comparison in cost function", false, 0.25, "FP_TYPE", cmd);
		TCLAP::ValueArg<FP_TYPE> w2_arg("2", "w2", "Weight for geostrophic wind in cost function", false, 0.625, "FP_TYPE", cmd);
		TCLAP::ValueArg<FP_TYPE> w3_arg("3", "w3", "Weight for curvature in cost function", false, 0.75, "FP_TYPE", cmd);
		TCLAP::SwitchArg text_out_arg("T", "text", "Output grid in text format, as well as the binary format", cmd, false);	
		TCLAP::ValueArg<FP_TYPE> sr_arg("r", "search_radius", "Radius, in km, to search when locating tracks", false, 500.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<FP_TYPE> min_tl_arg("l", "min_length", "Minimum track length in km", false, 0.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<FP_TYPE> min_dev_arg("d", "min_dev", "Minimum track deviation from first to last point in km", false, 0.0, "FP_TYPE", cmd);
		TCLAP::ValueArg<int> min_per_arg("p", "min_per", "Minimum persistence of the track in frames / timesteps", false, 1, "integer", cmd);
		TCLAP::ValueArg<std::string> output_fname_arg("o", "output", "Output file name", true, "", "string", cmd);
		TCLAP::MultiArg<std::string> input_fname_arg("i", "input", "Input file name", true, "string", cmd);
		cmd.parse(argc, argv);

		w0 = w0_arg.getValue();
		w1 = w1_arg.getValue();
		w2 = w2_arg.getValue();
		w3 = w3_arg.getValue();

		t = t_arg.getValue();
		s = s_arg.getValue();
		k = k_arg.getValue();
		if (k<=0) k=1;

		search_rad = sr_arg.getValue() * 1000.0;	// convert to metres
		min_tl = min_tl_arg.getValue() * 1000.0;	// convert to metres
		min_dev = min_dev_arg.getValue() * 1000.0;
		min_per = min_per_arg.getValue();
		output_fname = output_fname_arg.getValue();
		input_fname = input_fname_arg.getValue();
		text_out = text_out_arg.getValue();
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
		tracker TR(input_fname, min_per, min_tl, min_dev, search_rad);
		if (w0 != 1e6)
			TR.set_weights(w0, w1, w2, w3);
		TR.set_tsteps(t);
		TR.set_n_spl_pts(s);
		TR.set_ksteps(k);
		TR.find_tracks();
		TR.save(output_fname);
		std::cout << "# Saved to file: " << output_fname << std::endl;		
		if (text_out)
			TR.save_text(output_fname + ".txt");
	}
	catch(std::string &s)
	{
		std::cerr << s << std::endl;
	}
}
