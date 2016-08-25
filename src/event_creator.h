/******************************************************************************
** Program : event_creator.h
** Author  : Neil Massey
** Date    : 05/08/16
** Purpose : class that creates the windstorm events from the original netcdf
**           files and the track file produced from the regrid->extrema->track
**           tool chain
******************************************************************************/

#ifndef EVENT_CREATOR_H
#define EVENT_CREATOR_H

#include <string>
#include <list>
#include "ncdata.h"
#include "track_list.h"
#include "event.h"

/*****************************************************************************/

class event_creator
{
    public:
        event_creator(std::string input_track,
                      std::string mslp_fname,   std::string mslp_field,
                      std::string wind_fname,   std::string wind_field,
                      std::string precip_fname, std::string precip_field,
                      std::string pop_fname,    std::string pop_field,
                      std::string remap_fname,
                      std::string lsm_fname,    std::string lsm_field,
                      FP_TYPE search_rad,       int event_t_steps);
        ~event_creator(void);
        void find_events(void);     // function to start search procedure
        void save_events(std::string output_prefix); // write out the events
        
    private:
        // function to write event to
        void write_event(std::string out_fname, event* evt);
    
        // input track
        track_list tr_list;
        
        // ncdata objects for all of the netcdf input files
        ncdata mslp_data;
        ncdata wind_data;
        ncdata precip_data;
        ncdata pop_data;
        ncdata remap_slope, remap_icept;
        ncdata lsm_data;
        
        // time data
        int ref_year, ref_month, ref_day;
        FP_TYPE ref_day_scale, ref_n_days_py;
        
        // search radius
        float sr;
        // event time steps per day (i.e. 4 for 6hr data)
        int   evt;
        
        // output filename
        std::string out_fname;
        
        // list to store the events in
        std::list<event*> event_list;
};

#endif