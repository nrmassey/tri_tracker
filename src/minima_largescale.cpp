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
#include "geo_convert.h"

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
    stream >> ls_msh_lvl >> dummy >> contour_value >> dummy >> min_delta;    // dummy reads the comma
    
    // add to the metadata
    std::stringstream ss;
    meta_data["method"] = "minima_largescale";
    ss.str(""); ss << ls_msh_lvl;
    meta_data["large_scale_mesh_level"] = ss.str();
    ss.str(""); ss << contour_value;
    meta_data["contour_value"] = ss.str();
    ss.str(""); ss << min_delta;
    meta_data["min_delta"] = ss.str();
}

/******************************************************************************/

bool minima_largescale::process_data(void)
{
    // get the information for the current datastore
    int n_ts = ds.get_number_of_time_steps();
    int n_idxs = ds.get_number_of_indices();
    FP_TYPE mv = ds.get_missing_value();
    
    // create a new datastore
    data_processed->set_size(n_ts, n_idxs);
    data_processed->set_missing_value(mv);
    
    std::cout << "# Processing data" << std::endl;
    // create the smoothed largescale data
    data_store* data_smooth = create_smoothed_largescale();
    // now subtract the smoothed largescale data from the input data
    for (int t=0; t<n_ts; t++)
    {
        for (int d=0; d<n_idxs; d++)
        {
            // subtract the foreground (actual) data from the background (smoothed largescale)
            FP_TYPE fV = ds.get_data(t,d);
            FP_TYPE bV = data_smooth->get_data(t,d);
            FP_TYPE V = fV - bV;
            // set the data in the data processed data_store
            if (bV != mv && fV != mv)
                data_processed->set_data(t,d,V);
            else
                data_processed->set_data(t,d,mv);
        }
    }
//     // option to save the output - build the filename first
//     std::string out_fname = ds_fname.substr(0, ds_fname.size()-4);
//     std::stringstream ss;
//     // large scale field subtracted from data at extrema location level
//     ss << "_L" << tg.get_max_level()-1 << "_E" << extrema_level << "_S" << ls_msh_lvl << "_ls_del.rgd";
//     data_processed->save(out_fname+ss.str());
//     // just the smoothed large scale / background field
//     ss.str("");
//     ss << "_L" << tg.get_max_level()-1 << "_E" << extrema_level << "_S" << ls_msh_lvl << "_ls.rgd";
//     data_smooth->save(out_fname+ss.str());
    
    // clear the smoothed data
    delete data_smooth;
    
    // smooth the processed data below the large scale mesh level
    smooth_processed_data(ls_msh_lvl+1);
    return true;
}

/******************************************************************************/

int triangles_share_points(force_tri_3D* tri_a, force_tri_3D* tri_b)
{
    // record how many points the two triangles share.  NB - these triangles can
    // be at different levels in the mesh.
    int n_shared_pts = 0;
    // loop through each point in triangle A
    for (int i=0; i<3; i++)
        // loop through each point in triangle B
        for (int j=0; j<3; j++)
        {
            // note - it is sufficient to check that they share an index into
            // the point cloud, due to the method of subdividing the grid
            if ((*tri_a)[i] == (*tri_b)[j])
                n_shared_pts++;
        }
        
    // matching points will be detected twice
    return n_shared_pts;
}

/******************************************************************************/

void smooth_single_triangle(tri_grid* tg, QT_TRI_NODE* tri_qn, data_store* data_smooth)
{
    // debug option - whether we should smooth or just copy the triangle value
    bool smooth_on = true;
    // get all the adjacent triangles
    const LABEL_STORE* adj_labs = tri_qn->get_data()->get_adjacent_labels(POINT);
    FP_TYPE mv = data_smooth->get_missing_value();
    // for each child triangle, except the centroid
    for (int c=0; c<3; c++)
    {
        // get the child triangle node
        QT_TRI_NODE* child_qn = tri_qn->get_child(c);
        // check that it is not NULL
        if (child_qn == NULL)
            continue;
        indexed_force_tri_3D* child_tri = child_qn->get_data();
        // now build a list of indices into the data store where the adjacent
        // triangles at the tri_qn level share a point with the triangle at
        // the child_qn level
        std::vector<int> ds_indices;        // storage for ds indices
        
        if (smooth_on)
        {
            for (LABEL_STORE::const_iterator it_adj_labs = adj_labs->begin();
                 it_adj_labs != adj_labs->end(); it_adj_labs++)
            {
                // get the triangle from the tri grid
                indexed_force_tri_3D* adj_tri = tg->get_triangle(*it_adj_labs);
                // check for NULL triangle
                if (adj_tri == NULL)
                    continue;
                // check if they share a vertex
                if (triangles_share_points(child_tri, adj_tri) > 0)
                    // if they do then add this adjacent triangle at the level above 
                    // this child triangle to the list of ds indices
                    ds_indices.push_back(adj_tri->get_ds_index());
            }
        }
        // also add the current triangles ds index
        ds_indices.push_back(tri_qn->get_data()->get_ds_index());
        
        // now loop over each timestep and do the averaging into the child
        // triangle ds index
        int child_ds_index = child_tri->get_ds_index();
        int n_ds = ds_indices.size();
        for (int t=0; t<data_smooth->get_number_of_time_steps(); t++)
        {
            // do the averaging by looping over each ds index
            FP_TYPE sum = 0.0;
            int N = 0;
            // check first that there are ds_indices to loop over
            if (n_ds == 0)
                data_smooth->set_data(t, child_ds_index, mv);
            else
            {
                for (int d=0; d<n_ds; d++)
                {
                    FP_TYPE V = data_smooth->get_data(t, ds_indices[d]);
                    if (V != mv)
                    {
                        sum += V;
                        N += 1;
                    }
                }
                if (N > 0)
                    data_smooth->set_data(t, child_ds_index, sum / N);
                else
                    data_smooth->set_data(t, child_ds_index, mv);
            }
        }
    }
    // now assign the centroid to have just the value of the mean of the
    // 3 triangles above and the value of the parent triangle
    QT_TRI_NODE* child_qn = tri_qn->get_child(3);
    if (child_qn != NULL)
    {
        // get a list of the ds indices
        std::vector<int> ds_indices;
        if (smooth_on)
        {
            for (int c=0; c<3; c++)
            {
                // get the child triangle node
                QT_TRI_NODE* child_adj_qn = tri_qn->get_child(c);
                if (child_adj_qn != NULL)
                {
                    ds_indices.push_back(child_adj_qn->get_data()->get_ds_index());
                }
            }
        }
        // add the parent triangle
        ds_indices.push_back(tri_qn->get_data()->get_ds_index());

        for (int t=0; t<data_smooth->get_number_of_time_steps(); t++)
        {
            FP_TYPE sum_tris = 0.0;
            int N = 0;
            for (int d=0; d<ds_indices.size(); d++)
            {
                FP_TYPE V = data_smooth->get_data(t, ds_indices[d]);
                if (V != mv)
                {
                    sum_tris += V;
                    N += 1;
                }
            }
            if (N > 0)
                data_smooth->set_data(t, child_qn->get_data()->get_ds_index(), sum_tris / N);
            else
                data_smooth->set_data(t, child_qn->get_data()->get_ds_index(), mv);
        }
    }
}

/******************************************************************************/

data_store* minima_largescale::create_smoothed_largescale(void)
{
    // create the largescale background field by smoothing the field at 
    // level e-nup to levels l, where e is the level to detect extrema at,
    // nup is the number of levels to go up from e to get to the large scale
    // field at level s (nup = e-s) and l are levels below e (l > e)
    
    // create new smoothed data output
    data_store* data_smooth = new data_store();
    // copy the data from the non-processed datastore to the smoothed datastore
    data_smooth->copy(ds);

    // loop through each level from the large scale level to the max_level
    for (int l = ls_msh_lvl; l < tg.get_max_level(); l++)
    {
        // get all the triangles at this level
        std::list<QT_TRI_NODE*> tl = tg.get_triangles_at_level(l);
        // loop over all of these triangles
        for (std::list<QT_TRI_NODE*>::iterator it_tl = tl.begin(); it_tl != tl.end(); it_tl++)
        {
            // smooth this triangle
            smooth_single_triangle(&tg, *it_tl, data_smooth);
        }
    }
    return data_smooth;
}
