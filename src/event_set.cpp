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
#include "read_from_string.h"
#include "event_creator.h"

int main(int argc, char** argv)
{
    std::cout << "#### event set." << std::endl;
    std::string input_track;        // input track filename
    std::string output_fname;       // output event filename
    
    std::string mslp_fname;         // original MSLP filename
    std::string wind_fname;         // original 10m windspeed filename
    std::string precip_fname;       // original precipitation filename
    
    std::string mslp_field;         // variable name of MSLP in netCDF file
    std::string wind_field;         // variable name of wind in netCDF file
    std::string precip_field;       // variable name of precipitation in netCDF file
    
    std::string pop_fname;          // file name of population map (to calculate loss function)
    std::string pop_field;          // variable name of population in population map

    std::string remap_fname;        // file name for remap file to calculate 3s gust from max windspeed

    std::string lsm_fname;          // file name of land sea mask (to calculate max wind speed / power over land)
    std::string lsm_field;          // variable name of lsm in land sea mask file
    
    FP_TYPE search_rad;             // search radius (km)
    int event_t_steps;              // minimum number of timesteps in integrated event 
                                    // (e.g. 6 hourly data, 24 hour event - event_t_steps = 5)
    
    try
    {
        TCLAP::CmdLine cmd("Build an event set using tracks located from the tracker program");
        TCLAP::ValueArg<std::string> mslp_fname_arg("m", "mslp_file", "Original MSLP file used in tracking and field name. MSLP values for the event set will be taken from this file, rather than from the extrema file.", true, "", "string", cmd);       
        TCLAP::ValueArg<std::string> mslp_field_arg("M", "mslp_field", "Field name of MSLP in original MSLP file.", false, "field8", "string", cmd);
        TCLAP::ValueArg<std::string> wind_fname_arg("w", "wind_file", "Original wind file used in tracking and field name. 10m wind values for the event set will be taken from this file, rather than from the extrema file.", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> wind_field_arg("W", "wind_field", "Field name of 10m wind variable in original wind file.", false, "field50", "string", cmd);
        TCLAP::ValueArg<std::string> precip_fname_arg("p", "precip_file", "Original precipitation file used in tracking and field name. Precipitation values for the event set will be taken from this file, rather than from the extrema file.", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> precip_field_arg("P", "precip_field", "Field name of precipitation variable in original wind file.", false, "field90", "string", cmd);

        TCLAP::ValueArg<std::string> pop_fname_arg("q", "pop_file", "Population density file to be used in calculating loss function of storm.", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> pop_field_arg("Q", "pop_field", "Field name of population density variable in population density file.", false, "pop", "string", cmd);
        
        TCLAP::ValueArg<std::string> remap_fname_arg("e", "remap_file", "File used to remap max windspeed to 3s gust", true, "", "string", cmd);

        TCLAP::ValueArg<std::string> lsm_fname_arg("l", "lsm_file", "Land sea mask file to be used in calculating max wind power over land.", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> lsm_field_arg("L", "lsm_field", "Field name of land sea mask variable in land sea mask file.", false, "lsm", "string", cmd);

        TCLAP::ValueArg<int> event_t_steps_arg("t", "event_timesteps", "Minimum number of timesteps to use in integrated event set.", false, 12, "integer", cmd);
        TCLAP::ValueArg<FP_TYPE> sr_arg("r", "search_radius", "Radius, in km, to search when locating tracks.", false, 500.0, "FP_TYPE", cmd);

        TCLAP::ValueArg<std::string> output_fname_arg("o", "output", "Output prefix for event file names", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> input_fname_arg("i", "input", "Input file name - track file created by ./track", true, "", "string", cmd);
        cmd.parse(argc, argv);

        mslp_fname = mslp_fname_arg.getValue();
        mslp_field = mslp_field_arg.getValue();
        
        wind_fname = wind_fname_arg.getValue();
        wind_field = wind_field_arg.getValue();
        
        precip_fname = precip_fname_arg.getValue();
        precip_field = precip_field_arg.getValue();

        pop_fname = pop_fname_arg.getValue();
        pop_field = pop_field_arg.getValue();

        remap_fname = remap_fname_arg.getValue();

        lsm_fname = lsm_fname_arg.getValue();
        lsm_field = lsm_field_arg.getValue();
        
        event_t_steps = event_t_steps_arg.getValue();
        search_rad = sr_arg.getValue() * 1000.0;    // convert to metres

        output_fname = output_fname_arg.getValue();
        input_track = input_fname_arg.getValue();
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
    
        event_creator EV(input_track,
                         mslp_fname,   mslp_field,
                         wind_fname,   wind_field,
                         precip_fname, precip_field,
                         pop_fname,    pop_field,
                         remap_fname,
                         lsm_fname,    lsm_field,
                         search_rad,   event_t_steps);
        EV.find_events();
        EV.save_events(output_fname);
    }
    catch(std::string &s)
    {
        std::cerr << s << std::endl;
    }

}