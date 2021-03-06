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
    FP_TYPE search_rad;         // search radius (km)
    FP_TYPE hrs_per_ts;         // hours per timestep
    FP_TYPE overlap;            // percentage overlap
    bool text_out = false;
    int  opt_steps=1;           // number of optimisation iterations

    try
    {
        TCLAP::CmdLine cmd("Track extrema located using extrema program");
        TCLAP::ValueArg<int> opt_steps_arg("s", "opt_step", "Number of iterations to run the optimiser for", false, 1, "int", cmd);
        
        TCLAP::ValueArg<FP_TYPE> hrs_arg("f", "hours", "Hours per timestep in original data", false, 24, "FP_TYPE", cmd);
        TCLAP::ValueArg<FP_TYPE> ov_arg("v", "overlap", "Overlap between objects between timesteps", false, 50, "FP_TYPE", cmd);
        TCLAP::SwitchArg text_out_arg("T", "text", "Output tracks in text format, as well as the binary format", cmd, false);   
        TCLAP::ValueArg<FP_TYPE> sr_arg("r", "search_radius", "Radius, in km, to search when locating tracks", false, 500.0, "FP_TYPE", cmd);
        TCLAP::ValueArg<std::string> output_fname_arg("o", "output", "Output file name", true, "", "string", cmd);
        TCLAP::MultiArg<std::string> input_fname_arg("i", "input", "Input file name", true, "string", cmd);
        cmd.parse(argc, argv);

        opt_steps = opt_steps_arg.getValue();
        hrs_per_ts = hrs_arg.getValue();
        search_rad = sr_arg.getValue() * 1000.0;    // convert to metres
        overlap = ov_arg.getValue();
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
        tracker TR(input_fname, search_rad, overlap, hrs_per_ts, opt_steps);
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
