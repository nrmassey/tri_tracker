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
#include <list>
#include <vector>
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include <sys/stat.h>
#include <netcdf>
#include <netcdf.h>
#include "haversine.h"
#include "Rot2Global.h"

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

std::list<std::vector<int> > get_indices_in_circle(int ox, int oy, int max_x, 
                                                   int max_y, int rad)
{
    // use circle equation to determine which grid points are with rad grid
    // box radius of the origin (ox, oy)
    // error checking handled by passing in max_x and maxY
    std::list<std::vector<int> > out_list_of_indices;
    int R = rad*rad;
    for (int y = -rad; y <= rad; y++)
        for (int x = -rad; x <= rad; x++)
            if ((x*x + y*y <= R) and 
                (ox+x > 0) and (ox+x < max_x) and
                (oy+y > 0) and (oy+y < max_y))
            {
                std::vector<int> P(2);
                P[0] = ox+x;
                P[1] = oy+y;
                out_list_of_indices.push_back(P);
            }
    return out_list_of_indices;
}

/*****************************************************************************/

int calc_grid_box_radius(FP_TYPE sr, ncdata* input_data)
{
    // calculate the search radius in number of grid boxes
    // this is for a rotated grid where each grid box is (roughly) equal length
    FP_TYPE x_len = haversine(0.0, 0.0, input_data->get_lon_d(), 0.0, EARTH_R);
    // search radius is supplied in kms and EARTH_R is defined in kms
    int gsr = int(sr / x_len + 0.5);        // round up
    return gsr;
}

/*****************************************************************************/

FP_TYPE calc_grid_box_area(ncdata* input_data)
{
    FP_TYPE x_len = haversine(0.0, 0.0, input_data->get_lon_d(), 0.0, EARTH_R);
    FP_TYPE y_len = haversine(0.0, 0.0, 0.0, input_data->get_lat_d(), EARTH_R);
    return x_len * y_len;
}

/*****************************************************************************/

void event_creator::find_events(void)
{
    // find the events by looping over the tracks
    std::cout << "# Creating event set, track number: ";
    // calculate the search radius in grid boxes
    int gsr = calc_grid_box_radius(sr, &mslp_data);

    // get the slope and intercept for the wind gust
    field_data remap_slope_field = remap_slope.get_field();
    field_data remap_icept_field = remap_icept.get_field();
    field_data pop_field = pop_data.get_field();

    // get the pole latitude and longitude
    FP_TYPE plon, plat;
    if (mslp_data.has_rotated_grid())
    {
        plon = mslp_data.get_rotated_grid()->get_rotated_pole_longitude();
        plat = mslp_data.get_rotated_grid()->get_rotated_pole_latitude();
    }
    else
    {
        plon = 0.0;
        plat = 90.0;
    }

    // get the scaling factor for the timesteps for the wind data and the precip data
    // this allows precip and wind data at less than 6 hourly timesteps to be used
    // for example daily precip data
    FP_TYPE precip_t_step_scale = float(precip_data.get_t_len()) / mslp_data.get_t_len();
    FP_TYPE wind_t_step_scale = float(wind_data.get_t_len()) / mslp_data.get_t_len();

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
            FP_TYPE mv = new_event->mv;
            // loop over the track points
            for (int tp=0; tp < c_trk->get_persistence(); tp++)
            {
                // get the track point
                track_point* trk_pt = c_trk->get_track_point(tp);
                
                // get the indices from the track point
                // first convert to rotated grid if necessary
                FP_TYPE pt_lon, pt_lat;
                if (mslp_data.has_rotated_grid())
                {
                    Global2Rot(trk_pt->pt.lat, trk_pt->pt.lon, plat, plon,
                               pt_lat, pt_lon);
                }
                else
                {
                    pt_lon = trk_pt->pt.lon;
                    pt_lat = trk_pt->pt.lat;
                }
                int lon_idx = mslp_data.get_lon_idx(pt_lon);
                int lat_idx = mslp_data.get_lat_idx(pt_lat);
                int tstep = trk_pt->timestep;
                
                // append the lon, lat and timestep to the new_event
                new_event->track_lon.push_back(trk_pt->pt.lon);
                new_event->track_lat.push_back(trk_pt->pt.lat);
                new_event->timestep.push_back(trk_pt->timestep);
                
                // build the circle of indices
                std::list<std::vector<int> > circ_idxs = get_indices_in_circle(lon_idx, lat_idx, 
                                                                   mslp_data.get_lon_len(), 
                                                                   mslp_data.get_lat_len(), gsr);
                
                for (std::list<std::vector<int> >::iterator it_idx = circ_idxs.begin();
                     it_idx != circ_idxs.end(); it_idx++)
                {
                    // unpack the indices
                    int c_x, c_y;
                    c_x = (*it_idx)[0];
                    c_y = (*it_idx)[1];
                    // get the mslp value from the source data
                    FP_TYPE src_mslp_V = mslp_data.get_data(c_x, c_y, 0, tstep);
                    // get the mslp value from the target data
                    FP_TYPE tgt_mslp_V = new_event->mslp.get(c_x, c_y);
                    // if the src mslp is < tgt mslp or tgt mslp is mv
                    if ((tgt_mslp_V == mv) or (src_mslp_V < tgt_mslp_V))
                        new_event->mslp.set(c_x, c_y, src_mslp_V);
                        
                    // do the same for the wind data - need a consistency check with the
                    // mslp data first
                    if (tstep < int(float(wind_data.get_t_len()) * wind_t_step_scale))
                    {
                        FP_TYPE src_wind_V = wind_data.get_data(c_x, c_y, 0, int(float(tstep)*wind_t_step_scale));
                        FP_TYPE tgt_wind_V = new_event->wind_max.get(c_x, c_y);
                        if ((tgt_wind_V == mv) or (src_wind_V > tgt_wind_V))
                            new_event->wind_max.set(c_x, c_y, src_wind_V);
                    
                        // do the same for the precip data with the timestep check as well
                        FP_TYPE src_precip_V = precip_data.get_data(c_x, c_y, 0, int(float(tstep)*precip_t_step_scale));
                        FP_TYPE tgt_precip_V = new_event->precip.get(c_x, c_y);
                        if ((tgt_precip_V == mv) or (src_wind_V > tgt_precip_V))
                            new_event->precip.set(c_x, c_y, src_precip_V);
                    }
                }
                // calculate the wind gust - first copy the wind_max to wind_gust
                new_event->wind_gust.copy_ip(new_event->wind_max);
                // multiply by slope and add intercept (in place)
                new_event->wind_gust.mult_ip(remap_slope_field, mv);
                new_event->wind_gust.add_ip(remap_icept_field, mv);
                
                // calculate the loss using a simple loss model of population*wind power
                new_event->loss.copy_ip(new_event->wind_gust);      // copy wind gust
                // calculate wind power (w^3) by multiplying twice
                new_event->loss.mult_ip(new_event->wind_gust, mv);
                new_event->loss.mult_ip(new_event->wind_gust, mv);
                // multiply by population per km^2
                new_event->loss.mult_ip(pop_field, mv);
                // multiply by a scalar
                new_event->loss.mult_ip(1e-6,mv);
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
    
    // output string
    std::ostringstream output_string;

    // need the capability to do two separate strings - one for 365.25 day
    // data (e.g. PRECIS, ERA-I) and one for 360 day data (e.g. wah data, UM data)
    if (n_days_py == 365.25)
    {
        // create a boost posix time of the date and time start
        ptime ref_date_time(date(year, month, day), time_duration(0,0,0));
        // do the arithmetic on the time
        // number of hours that event occurs at since the reference
        time_duration n_hrs = hours(int((timestep * ref_data->get_t_d() + ref_data->get_t_s())*24+0.5));
        // create the time at which the event starts
        ptime trk_date_time = ref_date_time + n_hrs;
        // write out into output string    
        output_string << to_iso_extended_string(trk_date_time.date())
                      << "T" << std::setw(2) << std::setfill('0')
                      << trk_date_time.time_of_day().hours() << "-"
                      << std::setw(2) << std::setfill('0')
                      << trk_date_time.time_of_day().minutes() << "-"
                      << std::setw(2) << std::setfill('0')
                      << trk_date_time.time_of_day().seconds();
    }
    else if (n_days_py == 360.0)
    {
        // calculate number of days since the reference time
        FP_TYPE n_evt_days = timestep * ref_data->get_t_d() + ref_data->get_t_s();
        // create the number of days of the reference time and add on the event days
        FP_TYPE n_total_days = year * 360 + (month-1) * 30 + (day-1) + n_evt_days;
        // calculate the year
        int pyear = int(n_total_days / 360.0);
        n_total_days -= float(pyear) * 360;
        // calculate the month
        int pmonth = int(n_total_days/30);
        n_total_days -= float(pmonth) * 30;
        // calculate the day
        int pday = int(n_total_days) + 1;
        n_total_days -= float(pday-1);
        // calculate the hour
        int hrs = int(n_total_days*24);
        // write out the string
        output_string << std::setw(4) << std::setfill('0') << pyear << "-"
                      << std::setw(2) << std::setfill('0') << pmonth + 1 << "-"
                      << std::setw(2) << std::setfill('0') << pday << "T"
                      << std::setw(2) << std::setfill('0') << hrs << "-"
                      << "00-00";
    }
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

    // require the lsm to determine largest gust over land
    field_data lsm_field = lsm_data.get_field();

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
        write_event(out_path, cur_evt, &lsm_field);
    }
}

/*****************************************************************************/

BYTE* pack_data(FP_TYPE* data, int y_len, int x_len, 
                FP_TYPE offset, FP_TYPE scale, 
                FP_TYPE mv, BYTE scale_mv)
{
    // pack the floating point data into a byte - this will reduce the size
    // of the ensemble by 75%
    int data_size = y_len * x_len;
    // create the output data - need to delete this later
    BYTE* packed_data = new BYTE[data_size];
    // loop over the data
    for (int i=0; i<data_size; i++)
    {
        // replace mv with scaled mv
        if (data[i] == mv)
        {
            packed_data[i] = scale_mv;
        }
        else
        {
            FP_TYPE v = (data[i]- offset) * 1.0 / scale;
            int v1 = (int)(v);
            packed_data[i] = (BYTE)(v1);
        }
    }
    return packed_data;
}

/*****************************************************************************/

void event_creator::write_event(std::string out_fname, event* evt, field_data* lsm_field)
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
    netCDF::NcFile* px_nc_file = NULL;
    try
    {
        px_nc_file = new netCDF::NcFile(out_fname.c_str(), netCDF::NcFile::replace, 
                                        netCDF::NcFile::nc4);
    }   
    catch(netCDF::exceptions::NcException& e)
    {
        e.what();
        throw std::string("file " + out_fname + " could not be written to or already exists.");
    }
    // CF compliant file
    px_nc_file->putAtt("Conventions", "CF-1.6");

    // get the reference data from the event
    ncdata* ref_data = evt->ref_data;
    // create the longitude and latitude dimensions
    netCDF::NcDim x_rot_lat_dim = px_nc_file->addDim("latitude", ref_data->get_lat_len());
    netCDF::NcDim x_rot_lon_dim = px_nc_file->addDim("longitude", ref_data->get_lon_len());
    // create the longitude and latitude variables
    netCDF::NcVar x_rot_lat_var = px_nc_file->addVar("latitude", netCDF::ncFloat, x_rot_lat_dim);
    netCDF::NcVar x_rot_lon_var = px_nc_file->addVar("longitude", netCDF::ncFloat, x_rot_lon_dim);
    // add the attributes
    x_rot_lat_var.putAtt("standard_name", "grid_latitude");
    x_rot_lon_var.putAtt("standard_name", "grid_longitude");
    x_rot_lat_var.putAtt("units", "degrees");
    x_rot_lon_var.putAtt("units", "degrees");
    x_rot_lat_var.putAtt("axis", "Y");
    x_rot_lon_var.putAtt("axis", "X");

    // if this is a rotated grid then add the grid mapping variable
    if (ref_data->has_rotated_grid())
    {
        // NRM - note that the empty vector is a bit of a hack borrowed from the next
        // version of C++ netCDF
	    netCDF::NcVar x_grid_map_var = px_nc_file->addVar("rotated_pole", netCDF::ncChar,
	                                                      std::vector<netCDF::NcDim>());
        x_grid_map_var.putAtt("grid_mapping_name", "rotated_latitude_longitude");
        x_grid_map_var.putAtt("grid_north_pole_latitude", netCDF::ncFloat,
                              ref_data->get_rotated_grid()->get_rotated_pole_latitude());
        x_grid_map_var.putAtt("grid_north_pole_longitude", netCDF::ncFloat,
                              ref_data->get_rotated_grid()->get_rotated_pole_longitude());
    }

    // write the tracks out
    // first number of track points dimensions
    int n_trk_pts = evt->timestep.size();
    
    netCDF::NcDim x_trk_pts_dim = px_nc_file->addDim("track_point", n_trk_pts);
    // create the variables for the track
    netCDF::NcVar x_trk_time_var = px_nc_file->addVar("track_time", netCDF::ncFloat, x_trk_pts_dim);
    netCDF::NcVar x_trk_lat_var = px_nc_file->addVar("track_latitude", netCDF::ncFloat, x_trk_pts_dim);
    netCDF::NcVar x_trk_lon_var = px_nc_file->addVar("track_longitude", netCDF::ncFloat, x_trk_pts_dim);
    // add the attributes for the track lon and lat
    x_trk_lat_var.putAtt("standard_name", "latitude");
    x_trk_lon_var.putAtt("standard_name", "longitude");
    x_trk_lat_var.putAtt("units", "degrees_north");
    x_trk_lon_var.putAtt("units", "degrees_east");

    // add the attributes for the time variable
    int year, month, day;
    FP_TYPE day_scale, n_days_py;
    // first get the reference time
    ref_data->get_reference_time(year, month, day, day_scale, n_days_py);
    // output to string stream
    std::ostringstream time_units;
    time_units << "days since " << year << "-" << month << "-" << day << " 00:00:00";
    // add the attributes
    x_trk_time_var.putAtt("standard_name", "time");
    x_trk_time_var.putAtt("units", time_units.str());
    // add the calendar type
    if (n_days_py == 365.25)
        x_trk_time_var.putAtt("calendar", "standard");
    else if (n_days_py == 360.0)
        x_trk_time_var.putAtt("calendar", "360_day");
        
    // write the fields out need to pack them first
    // create each var first
    // need a vector of the dimensions
    std::vector<netCDF::NcDim> rot_dims(2);
    rot_dims[0] = x_rot_lat_dim;
    rot_dims[1] = x_rot_lon_dim;
    NC_SBYTE ncSbyte;
    
    netCDF::NcVar x_mslp_var = px_nc_file->addVar("mslp", ncSbyte, rot_dims);
    netCDF::NcVar x_wind_var = px_nc_file->addVar("wind_max", ncSbyte, rot_dims);
    netCDF::NcVar x_gust_var = px_nc_file->addVar("wind_gust", ncSbyte, rot_dims);
    netCDF::NcVar x_loss_var = px_nc_file->addVar("loss", ncSbyte, rot_dims);

    netCDF::NcVar x_precip_var = px_nc_file->addVar("precip", ncSbyte, rot_dims);

    // add the attributes
    // rotated grid mapping
    x_mslp_var.putAtt("grid_mapping", "rotated_pole");
    x_wind_var.putAtt("grid_mapping", "rotated_pole");
    x_gust_var.putAtt("grid_mapping", "rotated_pole");
    x_loss_var.putAtt("grid_mapping", "rotated_pole");
    x_precip_var.putAtt("grid_mapping", "rotated_pole");

    // missing value
    BYTE b_scale_mv = evt->scale_mv;
    x_mslp_var.putAtt("_FillValue", ncSbyte, b_scale_mv);
    x_wind_var.putAtt("_FillValue", ncSbyte, b_scale_mv);
    x_gust_var.putAtt("_FillValue", ncSbyte, b_scale_mv);
    x_loss_var.putAtt("_FillValue", ncSbyte, b_scale_mv);
    x_precip_var.putAtt("_FillValue", ncSbyte, b_scale_mv);
    
    // minimum / maximum values (depending on variables)
    FP_TYPE mslp_min = evt->mslp.get_min(evt->mv);
	x_mslp_var.putAtt("minimum", netCDF::ncFloat,  mslp_min);
    
    // wind maximum over land
    field_data wind_max_lsm(ref_data->get_lat_len(), ref_data->get_lon_len());
    wind_max_lsm.copy_ip(evt->wind_max);
    wind_max_lsm.mult_ip(*lsm_field, evt->mv);
    FP_TYPE wind_max = wind_max_lsm.get_max(evt->mv);
    x_wind_var.putAtt("maximum", netCDF::ncFloat, wind_max);
    FP_TYPE wind_power = wind_max*wind_max*wind_max;
    x_wind_var.putAtt("maximum_power", netCDF::ncFloat, wind_power);
    
    // calculate the gust maximum - which is multiplied by the lsm to get maximum over land
    field_data gust_max_lsm(ref_data->get_lat_len(), ref_data->get_lon_len());
    gust_max_lsm.copy_ip(evt->wind_gust);
    gust_max_lsm.mult_ip(*lsm_field, evt->mv);
    FP_TYPE gust_max = gust_max_lsm.get_max(evt->mv);
    x_gust_var.putAtt("maximum", netCDF::ncFloat, gust_max);
    FP_TYPE gust_power = gust_max*gust_max*gust_max;
    x_gust_var.putAtt("maximum_power", netCDF::ncFloat, gust_power);
    
    // precip and loss
    FP_TYPE max_precip = evt->precip.get_max(evt->mv);
    x_precip_var.putAtt("maximum", netCDF::ncFloat, max_precip);
    FP_TYPE max_loss = evt->loss.get_max(evt->mv);
    x_loss_var.putAtt("maximum", netCDF::ncFloat, max_loss);
    FP_TYPE sum_loss = evt->loss.get_sum(evt->mv);
    x_loss_var.putAtt("sum", netCDF::ncFloat, sum_loss);
    
    // standard names and units
    x_mslp_var.putAtt("standard_name", "air_pressure_at_sea_level");
    x_mslp_var.putAtt("long_name", "mean sea level pressure");
    x_mslp_var.putAtt("units", "Pa");
    x_wind_var.putAtt("standard_name", "wind_speed");
    x_wind_var.putAtt("long_name", "wind speed at 10m");
    x_wind_var.putAtt("units", "m s-1");
    x_gust_var.putAtt("standard_name", "wind_speed_of_gust");
    x_gust_var.putAtt("long_name", "3s peak gust wind speed at 10m");
    x_gust_var.putAtt("units", "m s-1");
    // loss has no standard name but it does have a long name and units
    x_loss_var.putAtt("long_name", "estimated losses");
    x_loss_var.putAtt("units", "1e6 m3 s-3 persons km^-2");
    
    x_precip_var.putAtt("standard_name", "precipitation_flux");
    x_precip_var.putAtt("long_name", "total precipitation flux");
    x_precip_var.putAtt("units", "kg m-2 s-1");

    // byte packing offsets and scales
    x_mslp_var.putAtt("add_offset", netCDF::ncFloat, evt->mslp_offset);
    x_mslp_var.putAtt("scale_factor", netCDF::ncFloat, evt->mslp_scale);

    x_wind_var.putAtt("add_offset", netCDF::ncFloat, evt->wind_offset);
    x_wind_var.putAtt("scale_factor", netCDF::ncFloat, evt->wind_scale);

    x_gust_var.putAtt("add_offset", netCDF::ncFloat, evt->wind_offset);
    x_gust_var.putAtt("scale_factor", netCDF::ncFloat, evt->wind_scale);

    x_loss_var.putAtt("add_offset", netCDF::ncFloat, evt->loss_offset);
    x_loss_var.putAtt("scale_factor", netCDF::ncFloat, evt->loss_scale);

    x_precip_var.putAtt("add_offset", netCDF::ncFloat, evt->precip_offset);
    x_precip_var.putAtt("scale_factor", netCDF::ncFloat, evt->precip_scale);

    // byte pack the data
    BYTE* px_byte_mslp = pack_data(evt->mslp.get(), 
                                   ref_data->get_lat_len(), ref_data->get_lon_len(), 
                                   evt->mslp_offset, evt->mslp_scale,
                                   evt->mv, evt->scale_mv);

    BYTE* px_byte_wind = pack_data(evt->wind_max.get(), 
                                   ref_data->get_lat_len(), ref_data->get_lon_len(), 
                                   evt->wind_offset, evt->wind_scale,
                                   evt->mv, evt->scale_mv);

    BYTE* px_byte_gust = pack_data(evt->wind_gust.get(), 
                                   ref_data->get_lat_len(), ref_data->get_lon_len(), 
                                   evt->wind_offset, evt->wind_scale,
                                   evt->mv, evt->scale_mv);

    BYTE* px_byte_loss = pack_data(evt->loss.get(), 
                                   ref_data->get_lat_len(), ref_data->get_lon_len(), 
                                   evt->loss_offset, evt->loss_scale,
                                   evt->mv, evt->scale_mv);

    BYTE* px_byte_precip = pack_data(evt->precip.get(), 
                                     ref_data->get_lat_len(), ref_data->get_lon_len(), 
                                     evt->precip_offset, evt->precip_scale,
                                     evt->mv, evt->scale_mv);

    // add the values for the longitude and latitude (on the rotated grid)
    std::vector<size_t> p(1);		// require vectors for addressing in new netCDF
    std::vector<size_t> c(1);
    p[0] = 0;
    c[0] = 1;
    for (int y=0; y < ref_data->get_lat_len(); y++)
    {
        FP_TYPE c_lat = ref_data->get_lat_s() + ref_data->get_lat_d()*y;
        x_rot_lat_var.putVar(p, c, &c_lat);
        p[0] += 1;
    }
    p[0] = 0;
    for (int x=0; x < ref_data->get_lon_len(); x++)
    {
        FP_TYPE c_lon = ref_data->get_lon_s() + ref_data->get_lon_d()*x;
        x_rot_lon_var.putVar(p, c, &c_lon);
        p[0] += 1;
    }

	// end file definition
	nc_enddef(px_nc_file->getId());

    // write the track lats and lons out
    std::list<FP_TYPE>::iterator track_lon_it = evt->track_lon.begin();   // position of track points
    std::list<FP_TYPE>::iterator track_lat_it = evt->track_lat.begin();   // position of track points
    std::list<int>::iterator     timestep_it  = evt->timestep.begin();    // timestep of track point

	// p & c are still available but need to be reset
	p[0] = 0;
	c[0] = 1;
    for (int t=0; t < n_trk_pts; t++)
    {
        FP_TYPE time = (*timestep_it * ref_data->get_t_d() + ref_data->get_t_s());
        x_trk_time_var.putVar(p,c,&time);
        
        x_trk_lon_var.putVar(p,c,&(*track_lon_it));

        x_trk_lat_var.putVar(p,c, &(*track_lat_it));
    
        track_lon_it++;
        track_lat_it++;
        timestep_it++;
        p[0]++;
    }

    // write the data
    // dimensions count
    std::vector<size_t> field_start(2);
    std::vector<size_t> field_count(2);
    field_start[0] = field_start[1] = 0;
    field_count[0] = ref_data->get_lat_len();
    field_count[1] = ref_data->get_lon_len();
    
    x_mslp_var.putVar(field_start, field_count, px_byte_mslp);
    x_wind_var.putVar(field_start, field_count, px_byte_wind);
    x_gust_var.putVar(field_start, field_count, px_byte_gust);
    x_loss_var.putVar(field_start, field_count, px_byte_loss);
    
    x_precip_var.putVar(field_start, field_count, px_byte_precip);
    
    delete px_nc_file;
}