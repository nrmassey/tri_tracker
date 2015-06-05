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
    
    std::cout << "# Processing data" << std::endl;
    
    // get the triangles at the extrema detection level
    // we only have to process this level and below
    for (int l=grid_level; l<tg.get_max_level(); l++)
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
            if (distance == 0.0)    // prevent a div by zero
                distance = 1.0;

            FP_TYPE W0 = 1000.0/distance;
//            FP_TYPE W0 = 1.0;
            int I0 = P_TRI->get_ds_index();
            FP_TYPE max_w = -1.0; 
            FP_TYPE min_w = 2e20; 
            // now do this for each triangle in the tri list
            for (LABEL_STORE::const_iterator p_tri_adj_it = p_tri_adj_labels->begin();
                 p_tri_adj_it != p_tri_adj_labels->end(); p_tri_adj_it++)
            {
                // get the distance
                distance = tg.distance_between_triangles(TRI->get_label(), *p_tri_adj_it);
                // get the triangle's ds index
                int p_adj_ds_index = tg.get_triangle(*p_tri_adj_it)->get_ds_index();
                FP_TYPE W1 = 1000.0/distance;
//                FP_TYPE W1 = 1.0;
                p_ds_indices.push_back(p_adj_ds_index);
                p_ds_weights.push_back(W1);
                if (W1 > max_w)
                    max_w = W1;
                if (W1 < min_w)
                    min_w = W1;
            }
            // restrict W0 to the maximum weight to prevent parent triangles that
            // are very close to the centroid of the triangle from dominating in the
            // removal of the field
            if (W0 > max_w)
                W0 = max_w;

            p_ds_indices.push_back(I0);
            p_ds_weights.push_back(W0);
            
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
    smooth_processed_data();
    // option to save the output - build the filename first
    std::string out_fname = ds_fname.substr(0, ds_fname.size()-4);
    std::stringstream ss;
    ss << "_L" << tg.get_max_level()-1 << "_E" << grid_level << "_S" << ls_msh_lvl << "_ls.rgd";
    out_fname += ss.str();
    data_processed->save(out_fname);
    return true;
}

/*****************************************************************************/

void minima_largescale::smooth_processed_data(void)
{
    std::cout << "# Smoothing processed data.";
    // create new smoothed data output
    data_store* data_smooth = new data_store();
    int n_ts = ds.get_number_of_time_steps();
    FP_TYPE mv = ds.get_missing_value();
    // create a new datastore
    data_smooth->set_size(n_ts, ds.get_number_of_indices());
    data_smooth->set_missing_value(mv);
    // smooth all the data below the extrema detection level
    for (int l=grid_level; l<tg.get_max_level(); l++)
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

    std::cout << std::endl;
}
