/******************************************************************************
** Program : event_creator.cpp
** Author  : Neil Massey
** Date    : 05/08/16
** Purpose : class that creates the windstorm events from the original netcdf
**           files and the track file produced from the regrid->extrema->track
**           tool chain
******************************************************************************/

#include "event_creator.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include <sys/stat.h>
#include "netcdfcpp.h"

/*****************************************************************************/

event_creator::event_creator(std::string input_track,
                             std::string mslp_fname,   std::string mslp_field,
                             std::string wind_fname,   std::string wind_field,
                             std::string precip_fname, std::string precip_field,
                             std::string pop_fname,    std::string pop_field,
                             std::string remap_fname,
                             std::string lsm_fname,    std::string lsm_field,
                             FP_TYPE search_rad,       int event_t_steps)
              : mslp_data(mslp_fname, mslp_field),
                wind_data(wind_fname, wind_field),
                precip_data(precip_fname, precip_field),
                pop_data(pop_fname, pop_field),
                remap_slope(remap_fname, "slope"),
                remap_icept(remap_fname, "intercept"),
                lsm_data(lsm_fname, lsm_field),
                sr(search_rad), evt(event_t_steps)
{
    // load the track list
    tr_list.load(input_track);
}

/*****************************************************************************/

event_creator::~event_creator(void)
{
    for (std::list<event*>::iterator it_evt_lst = event_list.begin();
         it_evt_lst != event_list.end(); it_evt_lst++)
    {
        delete *it_evt_lst;
    }
}

/*****************************************************************************/

void event_creator::find_events(void)
{
    // find the events by looping over the tracks
    std::cout << "# Creating event set, track number: ";
    // loop through all the tracks
    for (int tn=0; tn<tr_list.get_number_of_tracks(); tn++)
    {
        std::cout << tn;
        std::cout.flush();
        // get the track
        track* c_trk = tr_list.get_track(tn);
        // check the persistent is equal to or greater than evt
        if (c_trk->get_persistence() >= evt)
        {
            // create the event
            event* new_event = new event(mslp_data.get_lon_len(), mslp_data.get_lat_len());
            // set the mslp field as the reference data
            new_event->set_ref_data(&mslp_data);
            // loop over the track points
            for (int tp=0; tp < c_trk->get_persistence(); tp++)
            {
                // get the track point
                track_point* trk_pt = c_trk->get_track_point(tp);
                // get the indices from the track point
                int lon_idx = mslp_data.get_lon_idx(trk_pt->pt.lon);
                int lat_idx = mslp_data.get_lat_idx(trk_pt->pt.lat);
                
                // append the lon, lat and timestep to the new_event
                new_event->track_lon.push_back(trk_pt->pt.lon);
                new_event->track_lat.push_back(trk_pt->pt.lat);
                new_event->timestep.push_back(trk_pt->timestep);
            }
            // push the new event onto the event list
            event_list.push_back(new_event);
        }
        // track number counter
        int e = tn;
        if (tn == 0)
            std::cout << "\b";
        while (e > 0)
        {
            e = e / 10;
            std::cout << "\b";
        }
    }
    std::cout << std::endl;
    std::cout << "# Number of events: " << event_list.size() << std::endl;
}

/*****************************************************************************/

std::string get_event_start_date(ncdata* ref_data, int timestep)
{
    using namespace boost::gregorian;
    using namespace boost::posix_time;
    // build the filename containing the start date of the event
    int year, month, day;
    FP_TYPE day_scale, n_days_py;
    // first get the reference time and create a posix time from it
    ref_data->get_reference_time(year, month, day, day_scale, n_days_py);
    // create a boost posix time of the date and time start
    ptime ref_date_time(date(year,month,day),time_duration(0,0,0));
    // do the arithmetic on the time
    // number of hours that event occurs at since the reference
    time_duration n_hrs = hours(int((timestep * ref_data->get_t_d() + ref_data->get_t_s())*24+0.5));
    // create the time at which the event starts
    ptime trk_date_time = ref_date_time + n_hrs;
    
    // write out into output string
    std::ostringstream output_string;
    
    output_string << to_iso_extended_string(trk_date_time.date())
                  << "T" << std::setw(2) << std::setfill('0')
                  << trk_date_time.time_of_day().hours() << ":"
                  << std::setw(2) << std::setfill('0')
                  << trk_date_time.time_of_day().minutes() << ":"
                  << std::setw(2) << std::setfill('0')
                  << trk_date_time.time_of_day().seconds();
    return output_string.str();
}

/*****************************************************************************/

bool file_exists(const std::string& fpath)
{
    struct stat buf;
    int result = stat(fpath.c_str(), &buf);
    return result == 0;
}

/*****************************************************************************/

void event_creator::save_events(std::string output_prefix)
{
    // write out the events
    // these variables are needed to get the year, month and day of the beginning
    // of the event so that a sensible filename can be derived
    
    for (std::list<event*>::iterator it_evt = event_list.begin();
         it_evt != event_list.end(); it_evt++)
    {
        // get the event
        event* cur_evt = *it_evt;
        // get the reference time for this event
        std::string start_date = get_event_start_date(cur_evt->ref_data, cur_evt->timestep.front());
        // create the output file name
        std::string out_fname = output_prefix + "_" + start_date;
        // check whether the file exists or not and iterate until a unique filename is
        // created.  This is needed for a number of events which start on the same day
        char c_it = '0';
        std::string out_path = out_fname + "_" + c_it + ".nc";
        while (file_exists(out_path))
        {
            c_it++;
            out_path = out_fname + "_" + c_it + ".nc";
        }
        // write out the event
        write_event(out_path, cur_evt);
    }
}

/*****************************************************************************/

ncbyte* pack_data(FP_TYPE* data, int y_len, int x_len, 
                  FP_TYPE offset, FP_TYPE scale, 
                  FP_TYPE mv, ncbyte scale_mv)
{
    // pack the floating point data into a byte - this will reduce the size
    // of the ensemble by 75%
    int data_size = y_len * x_len;
    // create the output data - need to delete this later
    ncbyte* packed_data = new ncbyte[data_size];
    // loop over the data
    for (int i=0; i<data_size; i++)
    {
        // replace mv with scaled mv
        if (data[i] == mv)
        {
            data[i] = scale_mv;
        }
        else
        {
            FP_TYPE v = (data[i]- offset) * 1.0 / scale;
            v = (float)((int)(v));
            packed_data[i] = (ncbyte)(v);
            std::cout << packed_data[i] << " ";
        }
    }
    std::cout << std::endl;
    return packed_data;
}

/*****************************************************************************/

void event_creator::write_event(std::string out_fname, event* evt)
{
    // write a single event out.  Each event consists of a number of elements:
    // 1. Minimum MSLP field (2D)
    // 2. Maximum wind @ 10m field (2D)
    // 3. Wind gust @ 10m field (2D)
    // 4. Precipitation field (2D)
    // 5. Loss value field (2D)
    // 6. Max loss (single data point)
    // 7. Max wind power (single data point)
    // 8. Track longitudes (1D)
    // 9. Track latitudes (1D)
    // 10.Track timesteps (1D)
    
    // create the file for writing
    NcFile x_nc_file(out_fname.c_str(), NcFile::New);
    if (!x_nc_file.is_valid())
    {
        std::cerr << "Error: Output file: " << out_fname <<
                     " already exists or could not be opened" << std::endl;
        return;
    }
    // get the reference data from the event
    ncdata* ref_data = evt->ref_data;
    // create the longitude and latitude dimensions
    NcDim* px_rot_lat_dim = x_nc_file.add_dim("latitude",
                                              ref_data->get_lat_len());
    NcDim* px_rot_lon_dim = x_nc_file.add_dim("longitude",
                                              ref_data->get_lon_len());
    // create the longitude and latitude variables
    NcVar* px_rot_lat_var = x_nc_file.add_var("latitude", ncFloat, px_rot_lat_dim);
    NcVar* px_rot_lon_var = x_nc_file.add_var("longitude", ncFloat, px_rot_lon_dim);
    // add the attributes
    px_rot_lat_var->add_att("standard_name", "grid_latitude");
    px_rot_lon_var->add_att("standard_name", "grid_longitude");
    px_rot_lat_var->add_att("units", "degrees");
    px_rot_lon_var->add_att("units", "degrees");
    px_rot_lat_var->add_att("axis", "Y");
    px_rot_lon_var->add_att("axis", "X");

    // if this is a rotated grid then add the grid mapping variable
    NcVar* px_grid_map_var = NULL;
    if (ref_data->has_rotated_grid())
    {
        px_grid_map_var = x_nc_file.add_var("rotated_pole", ncChar, 0L);
        px_grid_map_var->add_att("grid_mapping_name", "rotated_latitude_longitude");
        px_grid_map_var->add_att("grid_north_pole_latitude",
                                 ref_data->get_rotated_grid()->get_rotated_pole_latitude());
        px_grid_map_var->add_att("grid_north_pole_longitude",
                                 ref_data->get_rotated_grid()->get_rotated_pole_longitude());
    }

    // add the values for the longitude and latitude (on the rotated grid)
    for (int y=0; y < ref_data->get_lat_len(); y++)
    {
        FP_TYPE c_lat = ref_data->get_lat_s() + ref_data->get_lat_d()*y;
        px_rot_lat_var->set_cur(y);
        px_rot_lat_var->put(&c_lat, 1L);
    }
    for (int x=0; x < ref_data->get_lon_len(); x++)
    {
        FP_TYPE c_lon = ref_data->get_lon_s() + ref_data->get_lon_d()*x;
        px_rot_lon_var->set_cur(x);
        px_rot_lon_var->put(&c_lon, 1L);
    }
    

    // write the tracks out
    // first number of track points dimensions
    int n_trk_pts = evt->timestep.size();
    
    NcDim* px_trk_pts_dim = x_nc_file.add_dim("track_point", n_trk_pts);
    // create the variables for the track
    NcVar* px_trk_time_var = x_nc_file.add_var("track_time", ncFloat, px_trk_pts_dim);
    NcVar* px_trk_lat_var = x_nc_file.add_var("track_latitude", ncFloat, px_trk_pts_dim);
    NcVar* px_trk_lon_var = x_nc_file.add_var("track_longitude", ncFloat, px_trk_pts_dim);
    // add the attributes for the track lon and lat
    px_trk_lat_var->add_att("standard_name", "latitude");
    px_trk_lon_var->add_att("standard_name", "longitude");
    px_trk_lat_var->add_att("units", "degrees_north");
    px_trk_lon_var->add_att("units", "degrees_east");

    // add the attributes for the time variable
    int year, month, day;
    FP_TYPE day_scale, n_days_py;
    // first get the reference time
    ref_data->get_reference_time(year, month, day, day_scale, n_days_py);
    // output to string stream
    std::ostringstream time_units;
    time_units << "days since " << year << "-" << month << "-" << day << " 00:00:00";
    // add the variable
    px_trk_time_var->add_att("standard_name", "time");
    px_trk_time_var->add_att("units", time_units.str().c_str());
    
    // write the track lats and lons out
    std::list<FP_TYPE>::iterator track_lon_it = evt->track_lon.begin();   // position of track points
    std::list<FP_TYPE>::iterator track_lat_it = evt->track_lat.begin();   // position of track points
    std::list<int>::iterator     timestep_it  = evt->timestep.begin();    // timestep of track point

    for (int t=0; t < n_trk_pts; t++)
    {
        FP_TYPE time = (*timestep_it * ref_data->get_t_d() + ref_data->get_t_s());
        px_trk_time_var->set_cur(t);
        px_trk_time_var->put(&time, 1L);
        
        px_trk_lon_var->set_cur(t);
        px_trk_lon_var->put(&(*track_lon_it), 1L);

        px_trk_lat_var->set_cur(t);
        px_trk_lat_var->put(&(*track_lat_it), 1L);
    
        track_lon_it++;
        track_lat_it++;
        timestep_it++;
    }
    
    // write the fields out need to pack them first
    // create each var first
    NcVar* px_mslp_var = x_nc_file.add_var("mslp", ncByte, px_rot_lat_dim, px_rot_lon_dim);
    NcVar* px_wind_var = x_nc_file.add_var("wind_max", ncByte, px_rot_lat_dim, px_rot_lon_dim);
    NcVar* px_gust_var = x_nc_file.add_var("wind_gust", ncByte, px_rot_lat_dim, px_rot_lon_dim);
    NcVar* px_loss_var = x_nc_file.add_var("loss", ncByte, px_rot_lat_dim, px_rot_lon_dim);
    NcVar* px_precip_var = x_nc_file.add_var("precip", ncByte, px_rot_lat_dim, px_rot_lon_dim);

    // add the attributes
    // rotated grid mapping
    px_mslp_var->add_att("grid_mapping", "rotated_pole");
    px_wind_var->add_att("grid_mapping", "rotated_pole");
    px_gust_var->add_att("grid_mapping", "rotated_pole");
    px_loss_var->add_att("grid_mapping", "rotated_pole");
    px_precip_var->add_att("grid_mapping", "rotated_pole");

    // missing value
    px_mslp_var->add_att("_FillValue", ncbyte(evt->scale_mv));
    px_wind_var->add_att("_FillValue", ncbyte(evt->scale_mv));
    px_gust_var->add_att("_FillValue", ncbyte(evt->scale_mv));
    px_loss_var->add_att("_FillValue", ncbyte(evt->scale_mv));
    px_precip_var->add_att("_FillValue", ncbyte(evt->scale_mv));
    
    // minimum / maximum values (depending on variables)
    px_mslp_var->add_att("minimum", evt->mslp.get_min());
    FP_TYPE wind_max = evt->wind_max.get_max();
    px_wind_var->add_att("maximum", wind_max);
    px_wind_var->add_att("maximum_power", wind_max*wind_max*wind_max);
    FP_TYPE gust_max = evt->wind_gust.get_max();
    px_wind_var->add_att("maximum", gust_max);
    px_wind_var->add_att("maximum_power", gust_max*gust_max*gust_max);
    px_precip_var->add_att("maximum", evt->precip.get_max());
    px_loss_var->add_att("maximum", evt->loss.get_max());
    
    // standard names and units
    px_mslp_var->add_att("standard_name", "air_pressure_at_sea_level");
    px_mslp_var->add_att("long_name", "mean sea level pressure");
    px_mslp_var->add_att("units", "Pa");
    px_wind_var->add_att("standard_name", "wind_speed");
    px_wind_var->add_att("long_name", "wind speed at 10m");
    px_wind_var->add_att("units", "m s-1");
    px_gust_var->add_att("standard_name", "wind_speed_of_gust");
    px_gust_var->add_att("long_name", "3s peak gust wind speed at 10m");
    px_gust_var->add_att("units", "m s-1");
    // loss has no standard name or units but it does have a long name
    px_loss_var->add_att("long_name", "estimated losses");
    px_loss_var->add_att("units", "1.0");
    px_precip_var->add_att("standard_name", "precipitation_amount");
    px_precip_var->add_att("long_name", "daily accumulated precipitation amount");
    px_precip_var->add_att("units", "kg m-2");

    // byte packing offsets and scales
    px_mslp_var->add_att("add_offset", evt->mslp_offset);
    px_mslp_var->add_att("scaling_factor", evt->mslp_scale);
    px_mslp_var->add_att("_Unsigned", "true");
    px_wind_var->add_att("add_offset", evt->wind_offset);
    px_wind_var->add_att("scaling_factor", evt->wind_scale);
    px_wind_var->add_att("_Unsigned", "true");
    px_gust_var->add_att("add_offset", evt->wind_offset);
    px_gust_var->add_att("scaling_factor", evt->wind_scale);
    px_gust_var->add_att("_Unsigned", "true");
    px_loss_var->add_att("add_offset", evt->loss_offset);
    px_loss_var->add_att("scaling_factor", evt->loss_scale);
    px_loss_var->add_att("_Unsigned", "true");
    px_precip_var->add_att("add_offset", evt->precip_offset);
    px_precip_var->add_att("scaling_factor", evt->precip_scale);
    px_precip_var->add_att("_Unsigned", "true");

    // byte pack the data
    ncbyte* px_byte_mslp = pack_data(evt->mslp.get(), 
                                     ref_data->get_lat_len(), ref_data->get_lon_len(), 
                                     evt->mslp_offset, evt->mslp_scale,
                                     evt->mv, evt->scale_mv);

    ncbyte* px_byte_wind = pack_data(evt->wind_max.get(), 
                                     ref_data->get_lat_len(), ref_data->get_lon_len(), 
                                     evt->wind_offset, evt->wind_scale,
                                     evt->mv, evt->scale_mv);

    ncbyte* px_byte_gust = pack_data(evt->wind_gust.get(), 
                                     ref_data->get_lat_len(), ref_data->get_lon_len(), 
                                     evt->wind_offset, evt->wind_scale,
                                     evt->mv, evt->scale_mv);

    ncbyte* px_byte_loss = pack_data(evt->loss.get(), 
                                     ref_data->get_lat_len(), ref_data->get_lon_len(), 
                                     evt->loss_offset, evt->loss_scale,
                                     evt->mv, evt->scale_mv);

    ncbyte* px_byte_precip = pack_data(evt->precip.get(), 
                                     ref_data->get_lat_len(), ref_data->get_lon_len(), 
                                     evt->precip_offset, evt->precip_scale,
                                     evt->mv, evt->scale_mv);
    
    // write the data
    px_mslp_var->put(px_byte_mslp, ref_data->get_lat_len(), ref_data->get_lon_len());
    px_wind_var->put(px_byte_wind, ref_data->get_lat_len(), ref_data->get_lon_len());
    px_gust_var->put(px_byte_gust, ref_data->get_lat_len(), ref_data->get_lon_len());
    px_loss_var->put(px_byte_loss, ref_data->get_lat_len(), ref_data->get_lon_len());
    px_precip_var->put(px_byte_precip, ref_data->get_lat_len(), ref_data->get_lon_len());

    x_nc_file.add_att("Conventions", "CF-1.6");
    x_nc_file.close();
    
    // delete the allocated memory
    delete [] px_byte_mslp;
    delete [] px_byte_wind;
    delete [] px_byte_gust;
    delete [] px_byte_loss;
    delete [] px_byte_precip;
}