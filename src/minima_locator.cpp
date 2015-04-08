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

bool is_mv(FP_TYPE V, FP_TYPE mv)
{
    return !(fabs(V) < fabs(mv*0.99));
}

/******************************************************************************/

minima_locator::minima_locator(void)
               :extrema_locator()
{
}

/*****************************************************************************/

minima_locator::~minima_locator(void)
{
}

/*****************************************************************************/

FP_TYPE minima_locator::trans_val(FP_TYPE val, vector_3D centroid)
{
    // transform the data to a contoured version which has the pressure gradient
    // removed
    FP_TYPE tri_lat, tri_lon;
    cart_to_model(centroid, tri_lon, tri_lat);

    const FP_TYPE deg_to_rad = M_PI / 180.0;
    FP_TYPE bck_val = ((equ_bck - pole_bck) * cos(tri_lat*deg_to_rad) + pole_bck);
    FP_TYPE val_2 = val - bck_val;
    FP_TYPE C = FP_TYPE(int(val_2/contour_value) * contour_value);
    return C;
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

bool minima_locator::is_extrema(indexed_force_tri_3D* tri, int t_step)
{
    // args in base class are encoded as follows:
    // thresh = threshold
    // min_delta = min delta
    // min_sur_tri = min number of surrounding triangles that are greater than current triangle
    // get the value of the triangle
    
    // build a list of surrounding triangles first
    const LABEL_STORE* tri_adj_labels = tri->get_adjacent_labels(adj_type);
    FP_TYPE mv = ds.get_missing_value();
    // get the value and transform the data to the contoured version minus the background field
    FP_TYPE tri_val = ds.get_data(t_step, tri->get_ds_index());
    FP_TYPE tri_val_C = trans_val(tri_val, tri->centroid());
    // if it's the missing value then return false
    if (is_mv(tri_val, mv))
        return false;
    // number of surrounding triangles that are greater than current triangle
    int n_st = 0;
    // loop through all the adjacent triangles
    for (LABEL_STORE::const_iterator tri_adj_it = tri_adj_labels->begin();
         tri_adj_it != tri_adj_labels->end(); tri_adj_it++)
    {
        // get the triangle from the label
        indexed_force_tri_3D* tri_adj = tg.get_triangle(*tri_adj_it);
        // get the value of the adjacent triangle
        FP_TYPE adj_val = ds.get_data(t_step, tri_adj->get_ds_index());
        // if it's the missing value then continue onto next one
        if (is_mv(adj_val, mv))
        {
            n_st += 1;
            continue;
        }
        // transform the data
        FP_TYPE adj_val_C = trans_val(adj_val, tri->centroid());
        
        // tri val less than surrounding triangles and less than min_delta
        if (tri_val_C <= min_delta && tri_val_C <= adj_val_C)
            n_st += 1;
    }
    return (n_st >= tri_adj_labels->size());
}

/*****************************************************************************/

bool minima_locator::is_in_object(indexed_force_tri_3D* O_TRI, 
                                  indexed_force_tri_3D* C_TRI, int t_step)
{
    // O_TRI - original triangle
    // C_TRI - candidate triangle - triangle being tested for inclusion
    bool is_in = false;

    // get the candidate triangle value
    FP_TYPE cl_v = ds.get_data(t_step, C_TRI->get_ds_index());
    FP_TYPE ol_v = ds.get_data(t_step, O_TRI->get_ds_index());

    // transform the data
    FP_TYPE cl_v_C = trans_val(cl_v, C_TRI->centroid());
    FP_TYPE ol_v_C = trans_val(ol_v, O_TRI->centroid());

    // quick check  
    is_in = cl_v_C <= ol_v_C;  // within 1 contour
    // not the mv
    is_in = is_in && fabs(cl_v) <= fabs(0.99 * ds.get_missing_value());
    // less than the minimum delta
    is_in = is_in && (cl_v_C <= min_delta);
    // less than 1000km radius
    is_in = is_in && tg.distance_between_triangles(O_TRI->get_label(), C_TRI->get_label())/1000.0 < 1000.0;
    return is_in;
}

/******************************************************************************/

void minima_locator::calculate_object_position(int o, int t)
{
    steering_extremum* svex = ex_list.get(t, o);
    LABEL_STORE* object_labels = &(svex->object_labels);
    
    // get all the points within the contour
    FP_TYPE mv = ds.get_missing_value();

    // position vector in Cartesian coordinates
    vector_3D P;
    FP_TYPE sum_V = 0.0;
    // get the position and weight for each triangle in the object
    for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
         it_ll != object_labels->end(); it_ll++)    // tri labels
    {
        // get the triangle and its centroid
        indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
        vector_3D C = c_tri->centroid();
        // get the data value from the datastore
        FP_TYPE V = ds.get_data(t, c_tri->get_ds_index());
        // the value is the weight
        if (!is_mv(V, mv))
        {
            P += C * fabs(V);
            sum_V += fabs(V);
        }
    }
    // divide by the weights and project to the sphere again
    if (sum_V > 0.0)
    {
        P *= 1.0 / (P.mag() * sum_V);
        // put the values back in the geo extremum
        FP_TYPE lon, lat;
        cart_to_model(P, lon, lat);
        svex->lon = lon;
        svex->lat = lat;
    }
    else
    {
        svex->lon = mv;
        svex->lat = mv;
    }
}

/*****************************************************************************/

void minima_locator::calculate_object_intensity(int o, int t)
{
    // get the extremum - the lat and lon will have been set already
    steering_extremum* svex = ex_list.get(t, o);
    
    // find min / max of values in the object
    FP_TYPE min_v = 2e20f;  // minimum value in object
    FP_TYPE max_v = 0;      // maximum value in object
    get_min_max_values(min_v, max_v, o, t);     
    svex->intensity = min_v;
}

/*****************************************************************************/

void minima_locator::calculate_object_delta(int o, int t)
{
    // calculate the delta as the absolute difference between the minimum and
    // the maximum value of the object
    // find min / max of values in the object
    FP_TYPE min_v = 2e20f;  // minimum value in object
    FP_TYPE max_v = -2e20f;     // maximum value in object
    get_min_max_values(min_v, max_v, o, t); 
    steering_extremum* svex = ex_list.get(t, o);
    svex->delta = max_v - min_v;
}