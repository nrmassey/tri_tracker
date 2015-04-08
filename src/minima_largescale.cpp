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

/******************************************************************************/

minima_largescale::minima_largescale(void) : extrema_locator()
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
    
    // calculate number of levels to go up the mesh to find the large scale flow
    n_up = grid_level - ls_msh_lvl;
    
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

FP_TYPE minima_largescale::get_val(indexed_force_tri_3D* TRI, int t)
{
     // get the n parent triangle by manipulating the label of the current
     // triangle
     long int label_int = (TRI->get_label().label) % (long int)(pow(10, n_up+2));
     // calculate division needed to go up to the level
     LABEL parent_label = LABEL(label_int, ls_msh_lvl);
     // get the parent triangle
     indexed_force_tri_3D* P_TRI = tg.get_triangle(parent_label);
     // subtract the parent value from the triangle value and contour
     FP_TYPE val = ds.get_data(TRI->get_ds_index(), t) - ds.get_data(P_TRI->get_ds_index(), t);
     FP_TYPE C_val = FP_TYPE(int(val/contour_value) * contour_value);
     return C_val;
}

/******************************************************************************/

void minima_largescale::calculate_object_position(int o, int t)
{
}

/******************************************************************************/

void minima_largescale::calculate_object_intensity(int o, int t)
{
}

/******************************************************************************/

void minima_largescale::calculate_object_delta(int o, int t)
{
}

/******************************************************************************/

bool minima_largescale::is_extrema(indexed_force_tri_3D* tri, int t_step)
{
    return false;
}

/******************************************************************************/

bool minima_largescale::is_in_object(indexed_force_tri_3D* O_TRI, 
                                     indexed_force_tri_3D* C_TRI, int t_step)
{
    return false;
}

/*********************************************************************/

void get_min_max_values_delta(FP_TYPE& min, FP_TYPE& max, int o, int t)
{
}