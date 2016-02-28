/******************************************************************************
** Program : extrema_locator.cpp
** Author  : Neil Massey
** Date    : 10/07/13
** Purpose : class that searches for extrema in data regridded onto a 
**           regional hierarchical triangular mesh - provided by regrid
**           this class is abstract and should be inherited from.  See
**           minima_locator.h and maxima_locator.h for an example
******************************************************************************/

#include <iostream>
#include <math.h>
#include <sstream>
#include "extrema_locator.h"
#include "geo_convert.h"
#include "haversine.h"
#include "concentric_shell.h"
#include "read_from_string.h"

/*****************************************************************************/

void extrema_locator::tstep_out_begin(int t)
{
    std::cout << t;
    std::cout.flush();
}

/*****************************************************************************/

void extrema_locator::tstep_out_end(int t)
{
    int e = t;
    if (t == 0)
        std::cout << "\b";
    while (e > 0)
    {
        e = e / 10;
        std::cout << "\b";
    }
}

/*****************************************************************************/

extrema_locator::extrema_locator(void) : sv(NULL)
{
}

/*****************************************************************************/

extrema_locator::~extrema_locator(void)
{
}

/*****************************************************************************/

void extrema_locator::locate(void)
{
    max_merge_dist = 500.0 * 1000.0;        // 500 km
    find_extrema();
    refine_extrema();
    find_objects();
    split_objects();
    merge_objects();
    ex_points_from_objects();
}

/*****************************************************************************/

void extrema_locator::save(std::string output_fname, bool save_text)
{
    int n_ex = 0;
    for (int t=0; t<ex_list.size(); t++)
        for (int o=0; o<ex_list.number_of_extrema(t); o++)
            if (!ex_list.get(t,0)->deleted)
                n_ex += 1;
    std::cout << "# Number of extrema: " << n_ex << std::endl;
    ex_list.set_meta_data(&meta_data);

    if (sv != NULL)
        ex_list.set_meta_data(sv->get_meta_data());
    ex_list.save(output_fname, ds.get_missing_value());
    if (save_text)
    {
        std::cout << output_fname + ".txt" << std::endl;
        ex_list.save_text(output_fname+".txt", &tg);
    }
}

/*****************************************************************************/

void extrema_locator::set_steering_vector(steering_vector* isv)
{
    sv = isv;
}

/*****************************************************************************/

void extrema_locator::set_inputs(std::string input_fname, std::string mesh_fname,
                                 int i_extrema_level, ADJACENCY i_adj_type)
{
    // load the data
    ds_fname = input_fname;
    ds.load(input_fname);
    tg.load(mesh_fname);
    // set the other inputs
    extrema_level = i_extrema_level;
    adj_type = i_adj_type;
    
    // create the meta data
    std::stringstream ss;
    meta_data["mesh_file_name"] = mesh_fname;
    meta_data["input_file_name"] = input_fname;
    ss << i_extrema_level;
    meta_data["extrema_grid_level"] = ss.str();
    meta_data["adjacency_type"] = i_adj_type == POINT ? "point" : "edge";
}

/*****************************************************************************/

void extrema_locator::find_extrema(void)
{
    // resize the extrema list to be the correct size for the number of timesteps
    ex_list.set_size(ds.get_number_of_time_steps());
    std::cout << "# Locating extrema, timestep: ";
    // get a list of all the triangles at the required level
    std::list<QT_TRI_NODE*> tris = tg.get_triangles_at_level(extrema_level);
    // get the missing value
    FP_TYPE mv = ds.get_missing_value();
    // repeat over all timesteps
    for (int t=0; t<ds.get_number_of_time_steps(); t++)
    {
        tstep_out_begin(t); 
        // repeat over all the nodes in this level
        for (std::list<QT_TRI_NODE*>::iterator it = tris.begin();
             it != tris.end(); it++)
        {
            indexed_force_tri_3D* c_tri = (*it)->get_data();
            if (is_extrema(c_tri, t))
            {
                // add to the extrema list, via the object list - the location and
                // value of the svex is currently just filled with the missing value
                steering_extremum svex(mv, mv, mv, mv, false, mv, mv);
                svex.object_labels.push_back(c_tri->get_label());
                // now add to the extrema list
                ex_list.add(t, svex);
            }
        }
        tstep_out_end(t);
    }
    std::cout << std::endl;
}

/*****************************************************************************/

void extrema_locator::get_leaf_node_labels(QT_TRI_NODE* c_tri_node, LABEL_STORE& label_list, int max_level)
{
    // we know how to form the labels at the child node - add 1 to 4 to the
    // end of the label and repeat this n times, where n is the number of
    // levels between the extrema detection level and the max level
    // doing this recursively makes for a very elegant program
    
    LABEL c_label = c_tri_node->get_data()->get_label();
    if (c_label.max_level < max_level)   // not at max level yet
    {
        // go to the child level
        for (int i=0; i<4; i++)
            if (c_tri_node->get_child(i) != NULL)
                get_leaf_node_labels(c_tri_node->get_child(i), label_list, max_level);
    }
    else
        label_list.push_back(c_label);
}

/*****************************************************************************/


void extrema_locator::refine_extrema(void)
{
    // refine the extrema to the maximum grid level by replacing the labels
    // of the extrema with the leaf node labels
    // we just have to check that the triangle of the leaf node actually
    // exists by checking for a NULL pointer returned from tg.get_triangle()
    
    // first check that the extrema level is not the maximum grid level
    int max_lev = tg.get_max_level() - 1;
    if (extrema_level == max_lev)
        return;
    
    std::cout << "# Refining extrema, timestep: ";
    // repeat for all detected extrema
    for (int t=0; t<ds.get_number_of_time_steps(); t++)
    {
        tstep_out_begin(t);
        for (int e=0; e<ex_list.number_of_extrema(t); e++)
        {
            steering_extremum* svex = ex_list.get(t, e);
            // get the labels at the leaf node descended from this label
            // should only be one object label at this stage
            LABEL_STORE label_list;
            QT_TRI_NODE* tri_node = tg.get_triangle_node(svex->object_labels[0]);
            get_leaf_node_labels(tri_node, label_list, max_lev);
            // clear the current object labels
            svex->object_labels.clear();
            svex->object_labels = label_list;
        }
        tstep_out_end(t);
    }
    std::cout << std::endl;
}

/*****************************************************************************/

void extrema_locator::find_objects(void)
{
    // at this point the extrema list contains the labels of the triangles
    // which have extrema points.  We want to grow these points into objects
    // which encompass (for example) the entire low pressure system
    std::cout << "# Locating objects from extrema, timestep: ";
    concentric_shell c_shell;
    LABEL_STORE shell_in_object;
    for (int t=0; t<ds.get_number_of_time_steps(); t++)
    {
        tstep_out_begin(t);
        for (int e=0; e<ex_list.number_of_extrema(t); e++)
        {
            c_shell.clear();
            // get the extremum for this object
            steering_extremum* svex = ex_list.get(t, e);
            // check to see whether a label actually exists
            if (svex->deleted)
                continue;
            // get the "original triangle" for this object 
            // - i.e. the one at the centre of the object
            indexed_force_tri_3D* O_TRI = tg.get_triangle(svex->object_labels[0]);
            bool keep_adding_to_obj = true;
            // continue until no more triangles are added to the object
            c_shell.calculate_inner_ring(&tg, svex);    // calculate the initial shape
            // calculate the 1st concentric shell
            c_shell.calculate(&tg, svex);
            while (keep_adding_to_obj)
            {
                keep_adding_to_obj = false;
                LABEL_STORE* c_shell_labs = c_shell.get_labels();
                shell_in_object.clear();
                // loop over all the labels
                for (LABEL_STORE::iterator it_c_shell_labs = c_shell_labs->begin();
                     it_c_shell_labs != c_shell_labs->end(); it_c_shell_labs++)
                {
                    // get the candidate triangle
                    indexed_force_tri_3D* C_TRI = tg.get_triangle(*it_c_shell_labs);
                    if (is_in_object(O_TRI, C_TRI, t))
                    {
                        keep_adding_to_obj = true;
                        // add to object but only if not already added
                        if (std::find(svex->object_labels.begin(), svex->object_labels.end(), 
                            *it_c_shell_labs) == svex->object_labels.end())
                        {
                            svex->object_labels.push_back(*it_c_shell_labs);
                            shell_in_object.push_back(*it_c_shell_labs);
                        }
                    }
                }
                c_shell.recalculate(&tg, &shell_in_object);
            }
            // set deleted if size is 0 - i.e. no triangles found
            if (svex->object_labels.size() == 0)
                svex->deleted = true;
        }
        tstep_out_end(t);
    }
    std::cout << std::endl;
}

/*****************************************************************************/

bool extrema_locator::objects_share_nodes(const LABEL_STORE* o1, 
                                          const LABEL_STORE* o2)
{
    // are two labels equal in the object lists
    for (LABEL_STORE::const_iterator it_o1 = o1->begin();
         it_o1 != o1->end(); it_o1++)
    {
        if (std::find(o2->begin(), o2->end(), *it_o1) != o2->end())
            return true;
    }
    return false;
}

/*****************************************************************************/
void extrema_locator::split_objects(void)
{
    // split the objects into constituent objects - i.e. one triangle per object
    extrema_list new_ex_list;
    // set the new extrema list to have the same number of timesteps as the current one
    new_ex_list.set_size(ex_list.size());
    std::cout << "# Splitting objects, timestep: ";
    for (int t=0; t<ex_list.size(); t++)
    {
        tstep_out_begin(t);
        // loop over every object for this timestep
        int o_s = ex_list.number_of_extrema(t);
        for (int o1=0; o1<o_s; o1++)
        {
           // loop over every label in this object
            steering_extremum* svex = ex_list.get(t, o1);
            // check to see whether a label actually exists
            if (svex == NULL || svex->deleted)
                continue;
            // create one object for every triangle in the current object
            for (int o2=0; o2<svex->object_labels.size(); o2++)
            {
                // create a new svex with just one label
                steering_extremum new_svex(svex->lon, svex->lat, 
                                           svex->intensity, svex->delta, 
                                           false,
                                           svex->sv_u, svex->sv_v);
                // add the label to the object
                new_svex.object_labels.push_back(svex->object_labels[o2]);
                // add the object (with one label / triangle) to the new
                // extrema list
                new_ex_list.add(t, new_svex);
            }
        }
        tstep_out_end(t);
    }
    // assign the extrema list to the new extrema list
    ex_list = new_ex_list;
    std::cout << std::endl;
}

/*****************************************************************************/

void extrema_locator::merge_objects(void)
{
    // merge objects together - two or more objects may have emerged from
    // extrema located close together
    std::cout << "# Merging objects, timestep: ";
    for (int t=0; t<ex_list.size(); t++)
    {
        tstep_out_begin(t);
        int o_s = ex_list.number_of_extrema(t);
        // calculate the concentric shells for each object
        std::vector<concentric_shell> obj_c_shells;
        for (int o1=0; o1<o_s; o1++)
        {
            concentric_shell o_c_shell;
            o_c_shell.calculate(&tg, ex_list.get(t, o1));
            obj_c_shells.push_back(o_c_shell);
        }
        for (int o1=0; o1<o_s; o1++)
        {           
            LABEL_STORE* o1_shell_labs = obj_c_shells[o1].get_labels();
            steering_extremum* obj_1 = ex_list.get(t, o1);
            if (obj_1->deleted) // deleted object as above
                continue;
            // get the labels in the object as well as the shell
            LABEL_STORE* o1_labs = &(obj_1->object_labels);
            for (int o2=o1+1; o2<o_s; o2++)
            {
                steering_extremum* obj_2 = ex_list.get(t, o2);
                if (obj_2->deleted)
                    continue;
                LABEL_STORE* o2_shell_labs = obj_c_shells[o2].get_labels();
                // get the labels for both objects
                LABEL_STORE* o2_labs = &(obj_2->object_labels);
                
                // is a label in the shell found in the object?
                // test at different levels - 1st shells overlap?
                bool test = false;
                // 1st test object shells overlap
                if (!test) test = objects_share_nodes(o1_shell_labs, o2_shell_labs);
                // 2nd - 1st object shell overlaps with 2nd object
                if (!test) test = objects_share_nodes(o1_shell_labs, o2_labs);
                // 3rd - 2nd object shell overlaps with 1st object
                if (!test) test = objects_share_nodes(o1_labs, o2_shell_labs);
                // 4th - labels in 1st object overlaps 2nd object
                if (!test) test = objects_share_nodes(o1_labs, o2_labs);

                if (test)
                {
                    // create two sets of new labels - one for each object
                    LABEL_STORE new_obj_1_L;
                    LABEL_STORE new_obj_2_L;
                    LABEL_STORE new_shell;
                    
                    // get the first object label in object 2
                    LABEL obj_2_centre = (*o2_labs)[0];
                    
                    // add object 2 labels to new object 2 labels
                    for (LABEL_STORE::iterator it_o2 = o2_labs->begin(); 
                         it_o2 != o2_labs->end(); it_o2++)
                        new_obj_2_L.push_back(*it_o2);

                    // loop through object 1's labels
                    for (LABEL_STORE::iterator it_o1 = o1_labs->begin(); 
                         it_o1 != o1_labs->end(); it_o1++)
                    {
                        // measure the distance between object 2's central triangle
                        // and the new triangle
                        FP_TYPE dist = tg.distance_between_triangles(obj_2_centre, *it_o1);
                        // if the distance is less than 1000km then add to object 2
                        if (dist < max_merge_dist)
                        {
                            // if it's not already in the object
                            if (find(new_obj_2_L.begin(), new_obj_2_L.end(), *it_o1) == new_obj_2_L.end())
                            {
                                new_obj_2_L.push_back(*it_o1);
                                // get the new labels to add to the shell
                                new_shell.push_back(*it_o1);
                            }
                        }
                        else
                        // otherwise add to new object 1's labels
                            new_obj_1_L.push_back(*it_o1);
                    }
                    // set labels
                    obj_2->object_labels = new_obj_2_L;
                    // if new obj 1 labels are not empty then add to object 1
                    // otherwise delete object 1
                    if (new_obj_1_L.size() > 1)
                        obj_1->object_labels = new_obj_1_L;
                    else
                        obj_1->deleted = true;
                    // recalculate the shell
                    obj_c_shells[o2].recalculate(&tg, &new_shell);
                }
            }
            // sort mostly for debugging purposes
//          sort(o1_labs->begin(), o1_labs->end());
        }
        tstep_out_end(t);
    }
    std::cout << std::endl;
}

/*****************************************************************************/

void extrema_locator::get_min_max_values(FP_TYPE& min_v, FP_TYPE& max_v, 
                                         int o, int t)
{
    // get the minimum and maximum values for an object containing a number
    // of labels
    min_v = 2e20f;
    max_v = -2e20f;
    LABEL_STORE* object_labels = &(ex_list.get(t, o)->object_labels);
    // if there are no labels do not try to find the min/max
    if (ex_list.get(t, o)->deleted)
        return;
    // get the missing value
    FP_TYPE mv = ds.get_missing_value();
    // loop over all the triangles in the object
    for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
         it_ll != object_labels->end(); it_ll++)    // tri indices
    {
        // get the triangle
        indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
        FP_TYPE val = ds.get_data(t, c_tri->get_ds_index());
        // find the min and max values
        if (fabs(val) < 0.99*fabs(mv) && val > max_v && fabs(val) < 1e10)
            max_v = val;
        if (fabs(val) < 0.99*fabs(mv) && val < min_v && fabs(val) < 1e10)
            min_v = val;
    }
}

/*****************************************************************************/

void extrema_locator::calculate_steering_vector(int o, int t)
{
    // get the extremum first and the list of labels in the object
    steering_extremum* svex = ex_list.get(t, o);
    FP_TYPE mv = ds.get_missing_value();
    if (svex != NULL)
        sv->calculate_steering_vector(&tg, svex, t, mv);
}

/*****************************************************************************/

void extrema_locator::ex_points_from_objects(void)
{
    // get the triangle from the tri_grid, via the label
    std::cout << "# Generating points from objects, timestep ";
    for (int t=0; t<ex_list.size(); t++)    // time
    {
        tstep_out_begin(t);
        for (int o=0; o<ex_list.number_of_extrema(t); o++)  // objects
        {
            if (!ex_list.get(t,o)->deleted)
            {
                calculate_object_position(o, t);
                if (sv != NULL)
                    calculate_steering_vector(o, t);
                calculate_object_delta(o, t);
                calculate_object_intensity(o, t);
            }
        }
        tstep_out_end(t);
    }
    std::cout << std::endl;
}