/******************************************************************************
** Program : minima_processed.h
** Author  : Neil Massey
** Date    : 10/04/15
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded and then processed in some way.  This is a
**           virtual class and needs to be inherited from.  It was written
**           as many feature identification routines have the same functions
**           after the data processing.  So all you have to do is write the
**           data processing function (and overload anything as necessary)
******************************************************************************/

#include "minima_processed.h"
#include "haversine.h"
#include "geo_convert.h"
#include "spline.h"
#include <sstream>
#include <math.h>

/******************************************************************************/

minima_processed::minima_processed(void) : extrema_locator(), data_processed(NULL)
{
    data_processed = new data_store();
}

/******************************************************************************/

minima_processed::~minima_processed(void)
{
    delete data_processed;
}

/******************************************************************************/

void minima_processed::locate(void)
{
    process_data();
    find_extrema();
    find_objects();
    merge_objects();
    ex_points_from_objects();
}

/******************************************************************************/

bool minima_processed::is_mv(FP_TYPE V, FP_TYPE mv)
{
    return !(fabs(V) < fabs(mv*0.99));
}

/******************************************************************************/

FP_TYPE minima_processed::contour_data(FP_TYPE V, FP_TYPE C)
{
    FP_TYPE c_val = FP_TYPE(int(V/C) * C);
    return c_val;
}

/******************************************************************************/

bool minima_processed::is_extrema(indexed_force_tri_3D* tri, int t_step)
{
    // get the data for the triangle for this time step and contour it
    FP_TYPE tri_val = data_processed->get_data(t_step, tri->get_ds_index());

    // if it's the missing value then return false
    FP_TYPE mv = data_processed->get_missing_value();
    if (is_mv(tri_val, mv))
        return false;

    // contour the data minus the background
    FP_TYPE tri_val_C = contour_data(tri_val, contour_value);

    // quickest check
    if (tri_val_C > min_delta)
        return false;

    int n_st = 0;
    // loop through all the adjacent triangles
    const LABEL_STORE* tri_adj_labels = tri->get_adjacent_labels(adj_type);
    // nearest neighbour search of all adjacent triangles for a minimum
    for (LABEL_STORE::const_iterator tri_adj_it = tri_adj_labels->begin();
         tri_adj_it != tri_adj_labels->end(); tri_adj_it++)
    {
        // get the triangle from the label
        indexed_force_tri_3D* tri_adj = tg.get_triangle(*tri_adj_it);
        // get the value of the adjacent triangle       
        FP_TYPE tri_adj_val = data_processed->get_data(t_step, tri_adj->get_ds_index());
        // if it's the missing value then continue onto next one - but count it as greater than current
        if (is_mv(tri_adj_val, mv))
        {
            n_st += 1;
            continue;
        }
        FP_TYPE tri_adj_val_C = contour_data(tri_adj_val, contour_value);
        // if the middle triangle is less than or equal to this surrounding triangle
        if (tri_val_C <= tri_adj_val_C)
            n_st += 1;
    }
    int min_sur = tri_adj_labels->size();
    return (n_st >= min_sur);
}

/******************************************************************************/

bool minima_processed::is_in_object(indexed_force_tri_3D* O_TRI, 
                                    indexed_force_tri_3D* C_TRI, int t_step)
{
    // O_TRI - original triangle
    // C_TRI - candidate triangle - triangle being tested for inclusion
    bool is_in = true;

    // get the candidate triangle value
    FP_TYPE cl_v = data_processed->get_data(t_step, C_TRI->get_ds_index());
    FP_TYPE ol_v = data_processed->get_data(t_step, O_TRI->get_ds_index());
    
    FP_TYPE cl_v_C = contour_data(cl_v, contour_value);
    FP_TYPE ol_v_C = contour_data(ol_v, contour_value);

    // not the mv
    FP_TYPE mv = data_processed->get_missing_value();
    is_in = is_in && !is_mv(cl_v, mv);
    // less than the minimum delta
    is_in = is_in && (cl_v_C <= min_delta);
    // less than 1000km radius
    is_in = is_in && tg.distance_between_triangles(O_TRI->get_label(), C_TRI->get_label())/1000.0 < 1000.0;

    if (is_in)
    {
        // saddle point check
        LABEL_STORE path = tg.get_path(O_TRI->get_label(), C_TRI->get_label(), grid_level);
        // create storage for spline
        std::vector<FP_TYPE> y_vals(path.size(), 0.0);
        std::vector<FP_TYPE> x_vals(path.size(), 0.0);
        // loop through and recover the values, adding to the arrays
        int i = 0;
        for (LABEL_STORE::iterator it_ll = path.begin(); 
             it_ll != path.end(); it_ll++) 
        {
            // get the triangle from the label
            indexed_force_tri_3D* I_TRI = tg.get_triangle(*it_ll);
            // get the value
            FP_TYPE il_v = data_processed->get_data(t_step, I_TRI->get_ds_index());
            x_vals[i] = (float)(i);
            y_vals[i] = il_v;
            i += 1;
        }
        // form the spline
        spline path_spline = spline(y_vals, x_vals, mv);
        // evaluate the 2nd differential
        FP_TYPE path_d2x = path_spline.evaluate_d2x(path.size()-1);
        is_in &= path_d2x > 0;
    }
    
    return is_in;
}

/******************************************************************************/

void minima_processed::get_min_max_values_processed(FP_TYPE& min_v, FP_TYPE& max_v, 
                                                    int o, int t)
{
    // get the minimum and maximum values for an object containing a number
    // of labels
    min_v = 2e20f;
    max_v = -2e20f;
    LABEL_STORE* object_labels = &(ex_list.get(t, o)->object_labels);
    // if there are no labels do not try to find the min/max
    if (object_labels->size() == 0)
        return;
    // get the missing value
    FP_TYPE mv = data_processed->get_missing_value();
    // loop over all the triangles in the object
    for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
         it_ll != object_labels->end(); it_ll++)    // tri indices
    {
        // get the triangle
        indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
        FP_TYPE val = data_processed->get_data(t, c_tri->get_ds_index());
        // find the min and max values
        if (!is_mv(val, mv) && val > max_v && fabs(val) < 1e10)
            max_v = val;
        if (!is_mv(val, mv) && val < min_v && fabs(val) < 1e10)
            min_v = val;
    }   
}

/******************************************************************************/

void minima_processed::calculate_object_position(int o, int t)
{
    steering_extremum* svex = ex_list.get(t, o);
    LABEL_STORE* object_labels = &(svex->object_labels);
    // find min / max of values in the object
    FP_TYPE min_v = 2e20f;  // minimum value in object
    FP_TYPE max_v = 0;      // maximum value in object
    get_min_max_values_processed(min_v, max_v, o, t);
    
    // get all the points within the contour
    min_v = contour_data(min_v, contour_value);
    FP_TYPE mv = data_processed->get_missing_value();

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
        FP_TYPE V = data_processed->get_data(t, c_tri->get_ds_index());
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

void minima_processed::calculate_object_intensity(int o, int t)
{
    // find min / max of values in the object
    FP_TYPE min_v = 2e20f;      // minimum value in object
    FP_TYPE max_v = -2e20f;     // maximum value in object
    get_min_max_values(min_v, max_v, o, t);
    steering_extremum* svex = ex_list.get(t, o);
    svex->intensity = min_v;
}

/*****************************************************************************/

void minima_processed::calculate_object_delta(int o, int t)
{
    // calculate the delta as the absolute difference between the minimum
    // and the background field
    // find min / max of values in the object
    FP_TYPE min_v = 2e20f;  // minimum value in object
    FP_TYPE max_v = -2e20f;     // maximum value in object
    get_min_max_values_processed(min_v, max_v, o, t);   
    steering_extremum* svex = ex_list.get(t, o);
    svex->delta = min_v;
}