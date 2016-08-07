/******************************************************************************
** Program : event_creator.h
** Author  : Neil Massey
** Date    : 05/08/16
** Purpose : class that creates the windstorm events from the original netcdf
**           files and the track file produced from the regrid->extrema->track
**           tool chain
******************************************************************************/

#ifndef <EVENT_CREATOR_H>
#define <EVENT_CREATOR_H>

#include <string>

/*****************************************************************************/

class event_creator
{
    public:
        event_creator(std::string input_track,  std::string output_fname,
                      std::string mslp_fname,   std::string mslp_field,
                      std::string wind_fname,   std::string wind_field,
                      std::string precip_fname, std::string precip_field,
                      std::string pop_fname,    std::string pop_field,
                      std::string remap_fname,
                      std::string lsm_fname,    std::string lsm_field,
                      FP_TYPE search_rad,       int event_t_steps);
        
    private:
};

#endif