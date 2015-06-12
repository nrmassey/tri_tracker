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
    refine_extrema();
    find_objects();
    merge_objects();
    trim_objects();
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

    int n_less = 0;
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
            n_less += 1;
            continue;
        }
        FP_TYPE tri_adj_val_C = contour_data(tri_adj_val, contour_value);
        // if the middle triangle is less than or equal to this surrounding triangle
        if (tri_val_C <= tri_adj_val_C)
            n_less += 1;
    }
    int min_sur = tri_adj_labels->size();
    bool is_min = (n_less >= min_sur);
    return is_min;
}

/******************************************************************************/

void minima_processed::refine_extrema(void)
{
    // refine the extrema to the lowest grid level first
    extrema_locator::refine_extrema();
    // now construct a new set of labels where the values are within one
    // contour of the minimum value
    for (int t=0; t<ds.get_number_of_time_steps(); t++)
    {
        tstep_out_begin(t);
        for (int e=0; e<ex_list.number_of_extrema(t); e++)
        {
            steering_extremum* svex = ex_list.get(t, e);
            // get the min max value
            FP_TYPE min_v, max_v;
            get_min_max_values_processed(min_v, max_v, e, t);
            min_v = contour_data(min_v, contour_value);
            // get the labels
            LABEL_STORE* object_labels = &(ex_list.get(t, e)->object_labels);
            // construct a new label list
            LABEL_STORE new_label_list;
            // loop over all labels
            for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
                 it_ll != object_labels->end(); it_ll++)
            {
                // get the value of the triangle
                indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
                FP_TYPE val = contour_data(data_processed->get_data(t, c_tri->get_ds_index()), contour_value);
                // add to new labels only if val is <= min_v
                if (val <= min_v)
                    new_label_list.push_back(*it_ll);
            }
            svex->object_labels.clear();
            svex->object_labels = new_label_list;
        }
        tstep_out_end(t);
    }
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
    is_in &= !is_mv(cl_v, mv);
    // less than the minimum delta
    is_in &= (cl_v_C <= min_delta);
    // less than 1000km radius
    is_in &= tg.distance_between_triangles(O_TRI->get_label(), C_TRI->get_label())/1000.0 < 1000.0;
    // within one contour
    is_in &= cl_v < ol_v_C + contour_value;
    
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
    FP_TYPE min_v = 2e20f;      // minimum value in object
    FP_TYPE max_v = -2e20f;     // maximum value in object
    get_min_max_values_processed(min_v, max_v, o, t);   
    steering_extremum* svex = ex_list.get(t, o);
    svex->delta = min_v;
}

/*****************************************************************************/

void minima_processed::trim_objects(void)
{
    std::cout << "# Trimming objects, timestep: ";
    FP_TYPE sum_o = 0.0;
    // get the surface area of a triangle
    steering_extremum* svex = ex_list.get(0, 0);
    if (svex == NULL)
    {
        std::cout << std::endl;
        return;
    }
    indexed_force_tri_3D* c_tri = tg.get_triangle(svex->object_labels[0]);
    FP_TYPE surf_A = c_tri->surface_area();
    
    for (int t=0; t<ex_list.size(); t++)
    {
        tstep_out_begin(t);
        int o_s = ex_list.number_of_extrema(t);
        sum_o += o_s;
        for (int o1=0; o1<o_s; o1++)
        {
            // remove those greater than minimum delta
            FP_TYPE min_vd, max_vd;
            get_min_max_values_processed(min_vd, max_vd, o1, t);
            if (min_vd > min_delta)
            {
                ex_list.get(t, o1)->object_labels.clear();// delete!
                sum_o -= 1;
            }
            // remove those with a surface area > 5,000km^2
            if (ex_list.get(t, o1)->object_labels.size()*surf_A > (5000)*(5000))
                ex_list.get(t, o1)->object_labels.clear();// delete!

        }
        tstep_out_end(t);   
    }
    std::cout << " Number of objects: " << sum_o << " ";
    std::cout << std::endl;
}

/*****************************************************************************/

void minima_processed::smooth_processed_data(int start_level)
{
    std::cout << "# Smoothing processed data." << std::endl;
    // create new smoothed data output
    data_store* data_smooth = new data_store();
    int n_ts = ds.get_number_of_time_steps();
    FP_TYPE mv = ds.get_missing_value();
    // create a new datastore
    data_smooth->set_size(n_ts, ds.get_number_of_indices());
    data_smooth->set_missing_value(mv);
    // smooth all the data below the extrema detection level
    for (int l=start_level; l<tg.get_max_level(); l++)
    {
        // get all the triangles at the current level
        std::list<QT_TRI_NODE*> tris_qn = tg.get_triangles_at_level(l);
        for (std::list<QT_TRI_NODE*>::iterator it_qt = tris_qn.begin();
             it_qt != tris_qn.end(); it_qt++)
        {
            // get the index into the datastore for this triangle
            indexed_force_tri_3D* TRI = (*it_qt)->get_data();
            int ds_idx = TRI->get_ds_index();
            // get the adjacent triangles
            const LABEL_STORE* tri_adj_labels = TRI->get_adjacent_labels(adj_type);
            // loop over every timestep
            for (int t=0; t<n_ts; t++)
            {
                // start the sum and sum of weights
                FP_TYPE sum_V = data_processed->get_data(t, ds_idx);
                FP_TYPE sum_W = 1.0;
                // check for missing value
                if (sum_V == mv)
                {
                    sum_V = 0.0;
                    sum_W = 0.0;
                }
                // scaling for surrounding triangles
                FP_TYPE S = 1.0/12;
                // loop over every adjacent triangle
                for (LABEL_STORE::const_iterator tri_adj_it = tri_adj_labels->begin();
                     tri_adj_it != tri_adj_labels->end(); tri_adj_it++)
                {
                    // get the adjacent triangle, the index and the value from the processed data
                    indexed_force_tri_3D* A_TRI = tg.get_triangle(*tri_adj_it);
                    int ds_a_idx = A_TRI->get_ds_index();
                    FP_TYPE V = data_processed->get_data(t, ds_a_idx);
                    if (V != mv)
                    {
                        sum_V += S*V;
                        sum_W += S;
                    }
                }
                // assign the value
                if (sum_W == 0.0)
                    data_smooth->set_data(t, ds_idx, mv);
                else
                    data_smooth->set_data(t, ds_idx, sum_V/sum_W);
            }
        }
    }
    // assign the smoothed data to the processed data
    delete data_processed;
    data_processed = data_smooth;
}