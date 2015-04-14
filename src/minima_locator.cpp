/******************************************************************************
** Program : minima_locator.cpp
** Author  : Neil Massey
** Date    : 10/07/13
** Modified: 09/09/14
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded onto a regional hierarchical triangular mesh
******************************************************************************/

#include "minima_locator.h"
#include <sstream>
#include <math.h>
#include "geo_convert.h"

/******************************************************************************/

minima_locator::minima_locator(void)
               :minima_processed()
{
}

/*****************************************************************************/

minima_locator::~minima_locator(void)
{
}

/*****************************************************************************/

void minima_locator::parse_arg_string(std::string method_string)
{
    // args for minima locator are:
    // arg[0] = pole background value
    // arg[1] = equator background value
    // arg[2] = contour level
    // arg[3] = min delta
    
    // get the first bracket
    int c_pos = method_string.find_first_of("(")+1;
    int e_pos = method_string.find(")", c_pos);
    char dummy;
    // check whether the text is "help" and print the method parameters if it is
    if (method_string.substr(c_pos, e_pos-c_pos) == "help")
        throw(std::string("minima parameters = (background value, contour value, minimum delta)"));
    std::stringstream stream(method_string.substr(c_pos, e_pos-c_pos));
    stream >> pole_bck >> dummy >> equ_bck >> dummy >> contour_value >> dummy >> min_delta ;    // dummy reads the comma
    // add to the meta data
    std::stringstream ss;
    meta_data["method"] = "minima";
    ss.str(""); ss << pole_bck;
    meta_data["pole_background_value"] = ss.str();
    ss.str(""); ss << equ_bck;
    meta_data["equator_background_value"] = ss.str();
    ss.str(""); ss << contour_value;
    meta_data["contour_value"] = ss.str();
    ss.str(""); ss << min_delta;
    meta_data["minimum_delta"] = ss.str();
}

/*****************************************************************************/

bool minima_locator::process_data(void)
{
    const FP_TYPE deg_to_rad = M_PI / 180.0;

    // get the information for the current datastore
    int n_ts = ds.get_number_of_time_steps();
    FP_TYPE mv = ds.get_missing_value();
    // create a new datastore
    data_processed->set_size(n_ts, ds.get_number_of_indices());
    data_processed->set_missing_value(mv);
    
    std::cout << "# Processing data" << std::endl;
    
     // get the triangles at the extrema detection level
    // we only have to process this level
    std::list<QT_TRI_NODE*> tris_qn = tg.get_triangles_at_level(grid_level);
    for (std::list<QT_TRI_NODE*>::iterator it_qt = tris_qn.begin();
         it_qt != tris_qn.end(); it_qt++)
    {
        // get the index into the tri-grid
        indexed_force_tri_3D* TRI = (*it_qt)->get_data();
        int ds_idx = TRI->get_ds_index();

        // get the latitude / longitude of this triangle's centroid
        FP_TYPE tri_lat, tri_lon;
        cart_to_model(TRI->centroid(), tri_lon, tri_lat);

        // calculate the background value for this triangle
        FP_TYPE bV = ((equ_bck - pole_bck) * cos(tri_lat*deg_to_rad) + pole_bck);

        // use this index over every timestep to subtract the background field
        for (int t=0; t<n_ts; t++)
        {
            // get the data value
            FP_TYPE dV = ds.get_data(t, ds_idx);
            // default to missing
            FP_TYPE V = mv;
            // if not missing the value is the data value minus the background
            if (!(is_mv(dV, mv)))
                V =  dV - bV;
            data_processed->set_data(t, ds_idx, V);
        }
    }
    // option to save the output - build the filename first
    std::string out_fname = ds_fname.substr(0, ds_fname.size()-4)+"_min.rgd";   
    data_processed->save(out_fname);

    return true;
}
