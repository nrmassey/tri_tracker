/******************************************************************************
** Program : minima_background.cpp
** Author  : Neil Massey
** Date    : 07/08/13
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded, and removing the background field
******************************************************************************/

#include "minima_background.h"
#include "haversine.h"
#include "geo_convert.h"
#include <sstream>
#include <math.h>

/******************************************************************************/

minima_background::minima_background(void) : minima_processed(), bck_field_ds(NULL)
{
}

/******************************************************************************/

minima_background::~minima_background(void)
{
    delete bck_field_ds;
}

/******************************************************************************/

void minima_background::parse_arg_string(std::string method_string)
{
    // arguments are:
    // arg[0] = file to read background
    // arg[1] = averaging period to take background field over
    // arg[2] = contour level
    // arg[3] = minimum delta
    // parameters for minima location with background removal
    // get the first bracket
    int c_pos = method_string.find_first_of("(")+1;
    int e_pos = method_string.find(")", c_pos);
    char dummy;
    if (method_string.substr(c_pos, e_pos-c_pos) == "help")
        throw(std::string("minima_back parameters = (file name to take background field from, averaging period of background field, contour value, min delta)"));
    
    int b_pos = method_string.find(",", c_pos);
    bck_field_file = method_string.substr(c_pos, b_pos-c_pos);
    std::stringstream stream(method_string.substr(b_pos+1, e_pos-b_pos));
    stream >> bck_avg_period >> dummy
           >> contour_value >> dummy
           >> min_delta;
    // add to the metadata
    std::stringstream ss;
    meta_data["method"] = "minima_back";
    meta_data["background_file"] = bck_field_file;
    ss.str(""); ss << bck_avg_period;
    meta_data["background_averaging_period"] = ss.str();
    ss.str(""); ss << contour_value;
    meta_data["contour_value"] = ss.str();
    ss.str(""); ss << min_delta;
    meta_data["minimum_delta"] = ss.str();
}

/******************************************************************************/

void minima_background::calculate_background_field(void)
{
    std::cout << "# Calculating background field" << std::endl;
    
    // load in the field first
    bck_field_ds = new data_store();
    bck_field_ds->load(bck_field_file);
    // get the size of the current datastore
    int n_ts = bck_field_ds->get_number_of_time_steps();
    int n_idx = bck_field_ds->get_number_of_indices();
    // do we need to take a mean?
    if (bck_avg_period > 1 && n_ts >= bck_avg_period)
    {
        // create a new background field
        data_store* new_bck_field_ds = new data_store();
        // set the scaling for each averaging period
        FP_TYPE scale = 1.0 / bck_avg_period;
        // create the datastore
        new_bck_field_ds->set_size(n_ts/bck_avg_period, n_idx);
        new_bck_field_ds->set_missing_value(bck_field_ds->get_missing_value());
        // loop through the data producing the average
        for (int t=0; t<n_ts; t++)
        {
            int dest_pos = t / bck_avg_period;
            for (int i=0; i<n_idx;  i++)
            {
                // get the current data value, add the value from the original ds store
                // multiplied by the scaler
                FP_TYPE c_val = new_bck_field_ds->get_data(dest_pos, i);
                FP_TYPE t_val = scale * bck_field_ds->get_data(t, i) + c_val;
                new_bck_field_ds->set_data(dest_pos, i, t_val);
            }
        }
        // assign the bck_field to be the new_field and delete the old one
        data_store* old_bck_field_ds = bck_field_ds;
        bck_field_ds = new_bck_field_ds;
        delete old_bck_field_ds;
    }
}

/******************************************************************************/

bool minima_background::process_data(void)
{
    // process the data
    // this removes the background field from the data
    // the background field may be averaged over a time period - i.e. remove
    // the monthly average / 5 day average or just daily average
    // the difference is taken between the triangle in the data and the corresponding
    // triangle in the background field.  This is repeated for all levels so that data
    // is subtracted from the background field at a number of resolutions
    
    // first calculate the background field
    calculate_background_field();

    // get details about the input data
    int n_ts = ds.get_number_of_time_steps();
    FP_TYPE mv = ds.get_missing_value();    
    int bck_field_nts = bck_field_ds->get_number_of_time_steps();

    // create the storage for the result
    data_processed->set_size(n_ts, ds.get_number_of_indices());
    data_processed->set_missing_value(mv);
    
    std::cout << "# Processing data" << std::endl;

    // we want to process every level, not just the grid level
    for (int l=0; l<tg.get_max_level(); l++)
    {
        // get the triangles for this level
        std::list<QT_TRI_NODE*> tris_qn = tg.get_triangles_at_level(l);
        for (std::list<QT_TRI_NODE*>::iterator it_qt = tris_qn.begin();
             it_qt != tris_qn.end(); it_qt++)
        {
            // get the index into the tri-grid
            int ds_idx = (*it_qt)->get_data()->get_ds_index();
            // use this index over every timestep to subtract the background field
            for (int t=0; t<n_ts; t++)
            {
                FP_TYPE dV = ds.get_data(t, ds_idx);
                // check for stepping outside of the array
                int bck_t = t/bck_avg_period;
                if (bck_t >= bck_field_nts)
                    bck_t = bck_field_nts-1;
                FP_TYPE bV = bck_field_ds->get_data(bck_t, ds_idx);
                FP_TYPE V = mv;
                if (!(is_mv(dV, mv) || is_mv(bV, mv)))
                    V =  dV - bV;
                data_processed->set_data(t, ds_idx, V);
            }
        }
    }
    // option to save the output - build the filename first
    std::string out_fname = ds_fname.substr(0, ds_fname.size()-4)+"_bck.rgd";
    data_processed->save(out_fname);

    return true;
}
