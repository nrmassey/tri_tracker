/******************************************************************************
** Program : regridder.cpp
** Author  : Neil Massey
** Date    : 28/06/13
** Purpose : class that takes a regional triangular grid and a netCDF file   
**           the netCDF file is regridded onto the mesh
******************************************************************************/

#include "regridder.h"
#include <iostream>
#include <stdlib.h>
#include <sstream>

/*****************************************************************************/

regridder::regridder(const std::string mesh_file_name, const std::string nc_fname, 
                     const std::string nc_vname, const int iz_level,
                     const int ip_method) 
          : nc_data(nc_fname, nc_vname), z_level(iz_level),
            p_method(ip_method)
{
    // open the grid file
    tg.load(mesh_file_name);
    // create the storage needed for the regridded values
    int n_tris = tg.get_max_ds_index();
    int n_tsteps = nc_data.get_t_len();
    
    // create the storage
    ds = new data_store();
    ds->set_size(n_tsteps, n_tris);
    ds->set_missing_value(nc_data.get_missing_value());
    // build the meta data
    META_DATA_TYPE meta_data;
    std::stringstream ss;

    meta_data["mesh_file_name"] = mesh_file_name;
    meta_data["input_file_name"] = nc_fname;
    meta_data["input_var_name"] = nc_vname;
    ss << iz_level;
    meta_data["z_level"] = ss.str();
    switch(ip_method)
    {
        case 0: meta_data["parent_method"] = "mean"; break;
        case 1: meta_data["parent_method"] = "min"; break;
        case 2: meta_data["parent_method"] = "max"; break;
    }
    ds->set_meta_data(&meta_data);
}

/*****************************************************************************/

regridder::~regridder(void)
{
    delete ds;
}

/*****************************************************************************/

void regrid_node(QT_TRI_NODE* current_node, ncdata* nc_data, int z_level, int t, 
                 data_store* ds, int p_method)
{
    // regridding whilst walking the tree
    // get the list of indices from the current triangle
    const std::list<grid_index>* grid_indices = current_node->get_data()->get_grid_indices();
    FP_TYPE sum = 0.0;
    FP_TYPE min = 2e20;
    FP_TYPE max = -2e20;
    int n = 0;
    FP_TYPE sum_W = 0.0;
    
    // loop through and recover the values from the nc_data
    for (std::list<grid_index>::const_iterator it = grid_indices->begin();
         it != grid_indices->end(); it++)
    {
#ifdef DEBUG
        if (it->i > nc_data.get_lon_len() || it->j > nc_data.get_lat_len())
            throw "Index out of range - incorrect grid used?";
#endif
        FP_TYPE d = nc_data->get_data(it->i, it->j, z_level, t);
        sum += d * it->W;
        sum_W += it->W;
        if (d < min)
            min = d;
        if (d > max)
            max = d;
        n ++;
    }
        
    // complete the average - check for no data and assign mv if non
    FP_TYPE val;
    if (n == 0)
        val = nc_data->get_missing_value();
    else
    {
        sum = sum / sum_W;
        // place the data - assign what the user wishes, either the mean, min or max
        switch(p_method) // user can choose: 0=mean, 1=min, 2=max
        {
            case 0: val = sum; break;
            case 1: val = min; break;
            case 2: val = max; break;
        }
    }
    
    ds->set_data(t, current_node->get_data()->get_ds_index(), val);
    
    // regrid any child nodes
    for (int i=0; i<4; i++)
        if (current_node->get_child(i) != NULL)
            regrid_node(current_node->get_child(i), nc_data, z_level, t, ds, p_method);
}

/*****************************************************************************/

void regridder::regrid(void)
{
    std::cout << "# Regridding, time step: ";
    // loop through all time steps
    for (int t=0; t<nc_data.get_t_len(); t++)
    {
        std::cout << t;
        std::cout.flush();
        // loop through all the base triangles
        for (int i=0; i<tg.get_number_of_base_tris(); i++)
        {
            // walk the tree doing the regridding as we go
            QT_TRI_NODE* root_node = tg.get_base_tri(i);
            regrid_node(root_node, &nc_data, z_level, t, ds, p_method);
        }
        // timestep counter
        int e = t;
        if (t == 0)
            std::cout << "\b";
        while (e > 0)
        {
            e = e / 10;
            std::cout << "\b";
        }
    }
    std::cout << std::endl;
}

/*****************************************************************************/

void regridder::save(std::string out_fname)
{
    ds->save(out_fname);
}

/*****************************************************************************/

void regridder::save_text(std::string out_fname)
{
    ds->save_text(out_fname);
}
