/******************************************************************************
** Program : event.h
** Author  : Neil Massey
** Date    : 23/08/16
** Purpose : class that stores the windstorm event and its associated field
**           data
******************************************************************************/

#include "field_data.h"
#include <list>
#include "ncdata.h"
#include "netcdfcpp.h"

class event
{
    friend class event_creator;
    public:
        event(int x_len, int y_len);
        void set_ref_data(ncdata* ref_data);  // copy the pointer to the mslp data
        
    private:
        field_data mslp;                // minimum mslp for the event
        field_data wind_max;            // 10 metre max 6 hourly wind speed
        field_data wind_gust;           // 10 metre gust - calculated from regression
        field_data precip;              // precipitation
        field_data loss;                // output from the loss function
        
        std::list<FP_TYPE> track_lon;   // position of track points
        std::list<FP_TYPE> track_lat;   // position of track points
        std::list<int> timestep;        // timestep of track point
        
        ncdata* ref_data;               // pointer to mslp data to get the grid definition from
        const FP_TYPE mv;               // missing value
        const ncbyte scale_mv;          // scaled missing value
        
        // scaling values for byte packing the data
        FP_TYPE mslp_offset, mslp_scale;
        FP_TYPE wind_offset, wind_scale;
        FP_TYPE precip_offset, precip_scale;
        FP_TYPE loss_offset, loss_scale;
};