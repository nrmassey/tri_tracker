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
    
    // calculate number of levels to go up the mesh to find the large scale flow
    int n_up = grid_level - ls_msh_lvl;
     // get the triangles at the extrema detection level
    // we only have to process this level
    std::list<QT_TRI_NODE*> tris_qn = tg.get_triangles_at_level(grid_level);
    for (std::list<QT_TRI_NODE*>::iterator it_qt = tris_qn.begin();
         it_qt != tris_qn.end(); it_qt++)
    {
        // get the index into the tri-grid
        indexed_force_tri_3D* TRI = (*it_qt)->get_data();
        int ds_idx = TRI->get_ds_index();
        // get the n parent triangle by manipulating the label of the current triangle
        long int ls_label_int = (TRI->get_label().label) % (long int)(pow(10, n_up+2));
        // calculate division needed to go up to the level
        LABEL ls_label = LABEL(ls_label_int, ls_msh_lvl);
        // get the parent triangle
        indexed_force_tri_3D* P_TRI = tg.get_triangle(ls_label);

        // we want to do a distance weighted mean for the background value so get
        // the list of adjacent triangles for the parent triangle
        const LABEL_STORE* p_tri_adj_labels = P_TRI->get_adjacent_labels(adj_type);
        
        // two lists, one of indices and one of weights
        std::vector<int> p_ds_indices;
        std::vector<FP_TYPE> p_ds_weights;
        
        // add the first - weight is 1000/distance
        FP_TYPE distance = tg.distance_between_triangles(TRI->get_label(), ls_label);
        p_ds_indices.push_back(P_TRI->get_ds_index());
        p_ds_weights.push_back(1000.0/distance);

        // now do this for each triangle in the tri list
        for (LABEL_STORE::const_iterator p_tri_adj_it = p_tri_adj_labels->begin();
             p_tri_adj_it != p_tri_adj_labels->end(); p_tri_adj_it++)
        {
            // get the distance
            distance = tg.distance_between_triangles(TRI->get_label(), *p_tri_adj_it);
            // get the triangle's ds index
            int p_adj_ds_index = tg.get_triangle(*p_tri_adj_it)->get_ds_index();
            p_ds_indices.push_back(p_adj_ds_index);
            p_ds_weights.push_back(1000.0/distance);
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
    // option to save the output - build the filename first
    std::string out_fname = ds_fname.substr(0, ds_fname.size()-4)+"_ls.rgd";   
    data_processed->save(out_fname);

    return true;
}
