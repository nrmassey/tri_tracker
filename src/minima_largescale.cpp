/******************************************************************************
** Program : minima_largescale.cpp
** Author  : Neil Massey
** Date    : 07/04/15
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded, after the removal of the large scale flow.
**           The large scale flow is defined as the regridded data at a lower
**           mesh resolution
******************************************************************************/

#include "minima_largescale.h"
#include <sstream>
#include <math.h>
#include <vector>

/******************************************************************************/

minima_largescale::minima_largescale(void) : minima_processed()
{
}

/******************************************************************************/

minima_largescale::~minima_largescale(void)
{
}

/******************************************************************************/

void minima_largescale::parse_arg_string(std::string method_string)
{
    // arguments are:
    // arg[0] = large scale features mesh level (e.g. 2)
    // arg[1] = contour level
    // arg[2] = minimum delta
    // parameters for minima location with large scale removal
    // get the first bracket
    int c_pos = method_string.find_first_of("(")+1;
    int e_pos = method_string.find(")", c_pos);
    char dummy;
    if (method_string.substr(c_pos, e_pos-c_pos) == "help")
        throw(std::string("minima_largescale parameters = (mesh level of large scale features, contour value, min delta)"));
    
    std::stringstream stream(method_string.substr(c_pos, e_pos-c_pos));
    stream >> ls_msh_lvl >> dummy >> contour_value >> dummy >> min_delta ;    // dummy reads the comma
    
    // add to the metadata
    std::stringstream ss;
    meta_data["method"] = "minima_largescale";
    ss.str(""); ss << ls_msh_lvl;
    meta_data["large_scale_mesh_level"] = ss.str();
    ss.str(""); ss << contour_value;
    meta_data["contour_value"] = ss.str();
    ss.str(""); ss << min_delta;
    meta_data["minimum_delta"] = ss.str();
}

/******************************************************************************/

bool minima_largescale::process_data(void)
{
    // get the information for the current datastore
    int n_ts = ds.get_number_of_time_steps();
    FP_TYPE mv = ds.get_missing_value();
    // create a new datastore
    data_processed->set_size(n_ts, ds.get_number_of_indices());
    data_processed->set_missing_value(mv);
    
    std::cout << "# Processing data." << std::endl;
    
    // get the triangles at the extrema detection level
    // we only have to process this level and below
/*    for (int l=grid_level; l<tg.get_max_level(); l++)
    {
        // calculate number of levels to go up the mesh to find the large scale flow
        int n_up = l - ls_msh_lvl;
        std::list<QT_TRI_NODE*> tris_qn = tg.get_triangles_at_level(l);
        for (std::list<QT_TRI_NODE*>::iterator it_qt = tris_qn.begin();
             it_qt != tris_qn.end(); it_qt++)
        {
            // get the index into the tri-grid
            indexed_force_tri_3D* TRI = (*it_qt)->get_data();
            int ds_idx = TRI->get_ds_index();
            
            // get the parent triangle
            QT_TRI_NODE* qt_p_node = (*it_qt);
            for (int f=0; f<=n_up; f++)
                qt_p_node = qt_p_node->get_parent();
            indexed_force_tri_3D* P_TRI = qt_p_node->get_data();
            
            // we want to do a distance weighted mean for the background value so get
            // the list of adjacent triangles for the parent triangle
            const LABEL_STORE* p_tri_adj_labels = P_TRI->get_adjacent_labels(adj_type);

            // two lists, one of indices and one of weights
            std::vector<int> p_ds_indices;
            std::vector<FP_TYPE> p_ds_weights;

            // add the first - weight is 1000/distance
            FP_TYPE distance = tg.distance_between_triangles(TRI->get_label(), P_TRI->get_label());
            FP_TYPE D = 1000.0;
            if (distance < 1000.0)    // prevent a div by zero or over weighting
                distance = 1000.0;
            FP_TYPE W0 = D/(distance);
            int I0 = P_TRI->get_ds_index();
            p_ds_indices.push_back(I0);
            p_ds_weights.push_back(W0);

            // now do this for each triangle in the tri list
            for (LABEL_STORE::const_iterator p_tri_adj_it = p_tri_adj_labels->begin();
                 p_tri_adj_it != p_tri_adj_labels->end(); p_tri_adj_it++)
            {
                // get the distance
                distance = tg.distance_between_triangles(TRI->get_label(), *p_tri_adj_it);
                // get the triangle's ds index
                int p_adj_ds_index = tg.get_triangle(*p_tri_adj_it)->get_ds_index();
                FP_TYPE W1 = D/(distance);
                p_ds_indices.push_back(p_adj_ds_index);
                p_ds_weights.push_back(W1);
            }

             
            // use this index over every timestep to subtract the background field
            for (int t=0; t<n_ts; t++)
            {
                // get the extrema level value
                FP_TYPE dV = ds.get_data(t, ds_idx);
                // default to missing
                FP_TYPE V = mv;
                // if not missing
                if (is_mv(dV, mv))
                    continue;
                // otherwise do the weighted sum of the parent triangle values
                // get the parent value
                FP_TYPE sum_pV = 0.0;
                FP_TYPE sum_w = 0.0;
                for (int c=0; c<p_ds_indices.size(); c++)
                {
                    FP_TYPE pV = ds.get_data(t, p_ds_indices[c]);
                    FP_TYPE W = p_ds_weights[c];
                    if (!is_mv(pV, mv))
                    {
                        sum_pV += pV * W;
                        sum_w += W;
                    }
                }
                if (sum_w > 0)
                    V =  dV - sum_pV / sum_w;
                else
                    V = mv;
                data_processed->set_data(t, ds_idx, V);
            }
        }
    }
    // smooth the data
    smooth_processed_data(grid_level);*/
    create_smoothed_largescale();
    // option to save the output - build the filename first
    std::string out_fname = ds_fname.substr(0, ds_fname.size()-4);
    std::stringstream ss;
    ss << "_L" << tg.get_max_level()-1 << "_E" << grid_level << "_S" << ls_msh_lvl << "_ls.rgd";
    out_fname += ss.str();
    data_processed->save(out_fname);
    return true;
}

/******************************************************************************/

void minima_largescale::create_smoothed_largescale(void)
{
    // create the largescale background field by smoothing the field at 
    // level e-nup to levels l, where e is the level to detect extrema at,
    // nup is the number of levels to go up from e to get to the large scale
    // field at level s (nup = e-s) and l are levels below e (l > e)
    
    // create new smoothed data output
    data_store* data_smooth = new data_store();
    int n_ts = ds.get_number_of_time_steps();
    FP_TYPE mv = ds.get_missing_value();
    // create a new datastore
    data_smooth->set_size(n_ts, ds.get_number_of_indices());
    data_smooth->set_missing_value(mv);
    
    // get and loop through the triangles at the large scale flow level
    std::list<QT_TRI_NODE*> tris_qn = tg.get_triangles_at_level(ls_msh_lvl);
    for (std::list<QT_TRI_NODE*>::iterator it_qt = tris_qn.begin();
         it_qt != tris_qn.end(); it_qt++)
    {
        LABEL SL = (*it_qt)->get_data()->get_label();
        // get the corner triangles at one level below
        LABEL_STORE corner_tris = tg.get_corner_child_triangles(SL, 1);
        // for each triangle, get the adjacent triangles to the corner tris
        for (int cti=0; cti<corner_tris.size(); cti++)
        {
            // get the starting label
            LABEL CL = corner_tris[cti];
            const LABEL_STORE* adj_corner_tris = tg.get_triangle(CL)->get_adjacent_labels(POINT);
            // now get the parents of these adjacent corner triangles.  Only add each parent once
            // this builds a list of triangles at the large scale flow level which surround the 
            // corner of the triangle.  Build these as ds indices and weights
            std::vector<int> ds_indices;
            std::vector<FP_TYPE> ds_weights;
            // add the first index with weight of 2 (twice as much as the surrounding ones
            indexed_force_tri_3D* C_TRI = tg.get_triangle(CL);
            ds_indices.push_back(C_TRI->get_ds_index());
            ds_weights.push_back(2.0);
            
            for (LABEL_STORE::const_iterator it_adj_corner_tris  = adj_corner_tris->begin();
                                             it_adj_corner_tris != adj_corner_tris->end();
                                             it_adj_corner_tris++)
            {
                // get the parent of the corner triangles
                indexed_force_tri_3D* p_adj_corner_tri = tg.get_triangle_node(CL)->get_parent()->get_data();
                // get the ds index label
                int ds_idx_adj_CT = p_adj_corner_tri->get_ds_index();
                // add if not already added
                if (std::find(ds_indices.begin(), ds_indices.end(), ds_idx_adj_CT) == ds_indices.end())
                {
                    ds_indices.push_back(ds_idx_adj_CT);
                    ds_weights.push_back(1.0);
                }
            }
            // we now have a list of indices and weights that we can use, through time and the mesh
            // levels to assign the large scale flow value to the corner triangles at higher resolutions
            // loop over the time
            // propagate this value down the mesh
            for (int l=1; l<=tg.get_max_level()-ls_msh_lvl; l++) // number of levels to descend
            {
                // get the corner tris
                LABEL_STORE child_corner_tris = tg.get_corner_child_triangles(SL, l);
                // get the corner tri that corresponds to the current corner we are in
                LABEL child_corner_tri = child_corner_tris[cti];
                // get the ds index
                int ds_idx = tg.get_triangle(child_corner_tri)->get_ds_index();
                for (int t=0; t<ds.get_number_of_time_steps(); t++)
                {
                    // calculate the triangle value
                    FP_TYPE sum_V = 0.0;
                    FP_TYPE sum_W = 0.0;
                    // weighted sum of surrounding triangles
                    for (int v=0; v<ds_indices.size(); v++)
                    {
                        FP_TYPE V1 = ds.get_data(t, ds_indices[v]);
                        if (V1 != mv)
                        {
                            sum_V += V1 * ds_weights[v];
                            sum_W += ds_weights[v];
                        }
                    }
                    FP_TYPE V = sum_V / sum_W;
                    // set the value
                    if (sum_W != 0.0)
                        data_smooth->set_data(t, ds_idx, V);
                }
            }
        }
        // now add the value for the centroid - this is just the value of the large scale triangle
        // need the ds index
        int SL_ds_idx = tg.get_triangle(SL)->get_ds_index();
        // propagate the value down the mesh
        for (int l=1; l<tg.get_max_level()-ls_msh_lvl; l++) // number of levels to descend
        {
            LABEL child_centroid_tri = tg.get_centroid_child_triangle(SL, l);
            if (child_centroid_tri.label != -1)
            {
                // get the index
                int ds_idx = tg.get_triangle(child_centroid_tri)->get_ds_index();
                for (int t=0; t<ds.get_number_of_time_steps(); t++)
                {
                    // get the value
                    FP_TYPE V = ds.get_data(t, SL_ds_idx);
                    // set the value
                    data_smooth->set_data(t, ds_idx, V);
                }
            }
        }
    }
    delete data_processed;
    data_processed = data_smooth;
}
