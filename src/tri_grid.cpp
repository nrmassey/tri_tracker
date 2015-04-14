/******************************************************************************
** Program : tri_grid.cpp
** Author  : Neil Massey
** Date    : 10/06/13
** Purpose : class to represent a triangular grid for the region
******************************************************************************/

#include "tri_grid.h"
#include "ncdata.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <sstream>
#include "geo_convert.h"
#include "bin_file_utils.h"
#include "haversine.h"
#include "Rot2Global.h"

/******************************************************************************/
// alias for adding points
#define PC point_cloud_instance
#define PC_ADD_ROT(X,Y,Z,T) point_cloud_instance.add_point(force_point(ROT_PT(vector_3D(X,Y,Z), T)))
#define PC_ADD(X,Y,Z)       point_cloud_instance.add_point(force_point(X,Y,Z))
#define TRI_ADD(P,X,Y,Z,S)  triangles.push_back(new QT_TRI(indexed_force_tri_3D((P), X, Y, Z, S)))
#define PC_ADD_SPLIT(X)     point_cloud_instance.add_point(X)

/******************************************************************************/

tri_grid::tri_grid(void)
{
}

/******************************************************************************/

tri_grid::~tri_grid(void)
{
    for (unsigned int i=0; i<triangles.size(); i++)
        delete triangles[i];
}

/******************************************************************************/

void tri_grid::initialize(SHAPE initial_shape, ncdata* nc_input_data,
                          int max_levs, int max_its, int perim)
{
    max_lev_input = max_levs;
    // initialise the eight or twelve quad trees with the icosahedron or octahedron
    // vertices
    assert(initial_shape == ICOSAHEDRON || initial_shape == OCTAHEDRON || initial_shape == DYMAXION);
    create_shape(initial_shape, 1.0);
    assign_points_from_grid(nc_input_data, perim);
    split_triangles(max_levs, nc_input_data);
    if (max_its != 0)
        PC.equalise(max_its, get_max_level());
    build_adjacency_maps();
    build_ds_indices();
    // build the metadata
    std::stringstream ss;

    meta_data["input_grid_file_name"] = nc_input_data->get_file_name();
    meta_data["input_grid_var_name"] = nc_input_data->get_file_name();
    switch(initial_shape)
    {
        case ICOSAHEDRON:
            meta_data["initial_shape"] = "icosahedron";
        break;
        case OCTAHEDRON:
            meta_data["initial_shape"] = "octahedron";
        break;
        case DYMAXION:
            meta_data["initial_shape"] = "dymaxion";
        break;
    }
    ss << max_levs;
    meta_data["maximum_levels"] = ss.str();
    ss.str(""); ss << max_its;
    meta_data["maximum_iterations"] = ss.str();
}

/*****************************************************************************/

META_DATA_TYPE* tri_grid::get_meta_data(void)
{
    return &meta_data;
}

/*****************************************************************************/

indexed_force_tri_3D* tri_grid::get_triangle(LABEL label)
{
    // get the triangle from the label.
    return get_triangle_node(label)->get_data();
}

/*****************************************************************************/

QT_TRI_NODE* tri_grid::get_triangle_node(LABEL label)
{
    // get the triangle node from the label.
    // the label is now completely numeric and has the following format
    
    // two least significant figures : parent triangle (0->19)
    // next sig fig : level 1 triangle (100,200,300,400)
    // next sig fig : level 2 triangle (1000,2000,3000,4000) etc.
    // Add these together to get the label: 2318, 1100, etc.
    
    // first two powers of ten are the first triangle
    long int tri_number = label.label % 100;
    // loop through the rest of the label, taking the requisite child
    QT_TRI_NODE* current = triangles[tri_number]->get_root();
    // these two values keep track of which significant figure we are at
    long int mp1 = 1000;
    long int mp2 = 100;
    for (unsigned int i=0; i<label.max_level; i++)
    {
        long int child_number = (label.label % mp1)/mp2;
        current = current->get_child(child_number-1);
        assert(current != NULL);
        mp1 *= 10;
        mp2 *= 10;
    }
    return current;
}

/*****************************************************************************/

void get_max_ds_index_node(QT_TRI_NODE* current, int& max_tgt_idx)
{
    if (current->get_data()->get_ds_index() > max_tgt_idx)
        max_tgt_idx = current->get_data()->get_ds_index();
    for (int i=0; i<4; i++)
        if (current->get_child(i) != NULL)
            get_max_ds_index_node(current->get_child(i), max_tgt_idx);
}

/*****************************************************************************/

int tri_grid::get_max_ds_index(void)
{
    int max_tgt_idx = -1;
    for(unsigned int i=0; i<triangles.size(); i++)
    {
        get_max_ds_index_node(triangles[i]->get_root(), max_tgt_idx);
    }
    return max_tgt_idx;
}

/******************************************************************************/

int tri_grid::get_max_level(void)
{
    std::vector<int> depths(triangles.size());
    for (unsigned int tri=0; tri<triangles.size(); tri++)
        depths[tri] = triangles[tri]->get_max_level();
    std::sort(depths.begin(), depths.end());
    return depths[triangles.size()-1];
}

/*****************************************************************************/

bool tri_grid::point_in_tri(const vector_3D* P, const force_tri_3D* T)
{
    // project the point to the plane
    const FP_TYPE TOL = 0.0;
    bool t1a = (PC[(*T)[0]].xp(PC[(*T)[1]])).dp(*P) > (0-TOL);
    bool t2a = (PC[(*T)[1]].xp(PC[(*T)[2]])).dp(*P) > (0-TOL);
    bool t3a = (PC[(*T)[2]].xp(PC[(*T)[0]])).dp(*P) > (0-TOL);
    return t1a && t2a && t3a;
}

/****************************************************************************/

vector_3D ROT_PT(const vector_3D& P, const vector_3D& T)
{
    // rotate around the x-axis
    vector_3D X1;
    FP_TYPE CX = cos(T[0]);
    FP_TYPE SX = sin(T[0]);
    X1[0] = P[0];
    X1[1] = P[1] * CX - P[2] * SX;
    X1[2] = P[1] * SX + P[2] * CX;

    // rotate around the y-axis
    vector_3D X2;
    FP_TYPE CY = cos(T[1]);
    FP_TYPE SY = sin(T[1]);
    X2[0] = X1[0] * CY + X1[2] * SY;
    X2[1] = X1[1];
    X2[2] = X1[2] * CY - X1[0] * SY;
    
    // rotate around the z-axis
    vector_3D X3;
    FP_TYPE CZ = cos(T[2]);
    FP_TYPE SZ = sin(T[2]);
    X3[0] = X2[0] * CZ - X2[1] * SY;
    X3[1] = X2[0] * SZ + X2[1] * CZ;
    X3[2] = X2[2];
    
    return X3;
}

/*****************************************************************************/

void tri_grid::create_shape(SHAPE initial_shape, FP_TYPE R)
{
    // create the initial shape
    switch (initial_shape)
    {
        case ICOSAHEDRON:
        {
            // icosahedron defined by 12 points and 20 faces
            // calculate the common values for the corners
            FP_TYPE b   = cos(M_PI/5.0);
            FP_TYPE c   = sqrt(1+4*b*b);
            FP_TYPE phi = (2 * cos(M_PI/5.0)) / c;
            FP_TYPE the = 1.0 / c;

            // create the vertices
            // rotate the original points so that one triangle has its
            // centroid at the pole
            vector_3D T(M_PI/8.0, 0.0, 0.0);

            int idx_0  = PC_ADD_ROT( 0.0,  phi,  the, T);
            int idx_1  = PC_ADD_ROT( 0.0, -phi,  the, T);
            int idx_2  = PC_ADD_ROT( 0.0,  phi, -the, T);
            int idx_3  = PC_ADD_ROT( 0.0, -phi, -the, T);
            int idx_4  = PC_ADD_ROT( the,  0.0,  phi, T);
            int idx_5  = PC_ADD_ROT(-the,  0.0,  phi, T);
            int idx_6  = PC_ADD_ROT( the,  0.0, -phi, T);
            int idx_7  = PC_ADD_ROT(-the,  0.0, -phi, T);
            int idx_8  = PC_ADD_ROT( phi,  the,  0.0, T);
            int idx_9  = PC_ADD_ROT(-phi,  the,  0.0, T);
            int idx_10 = PC_ADD_ROT( phi, -the,  0.0, T);
            int idx_11 = PC_ADD_ROT(-phi, -the,  0.0, T);

            // create the faces as root nodes in the quadtree
            TRI_ADD(&PC, idx_1,  idx_3,  idx_10, LABEL(0L, 0));
            TRI_ADD(&PC, idx_4,  idx_1,  idx_10, LABEL(1L, 0));
            TRI_ADD(&PC, idx_8,  idx_4,  idx_10, LABEL(2L, 0));
            TRI_ADD(&PC, idx_6,  idx_8,  idx_10, LABEL(3L, 0));
            TRI_ADD(&PC, idx_6,  idx_10, idx_3,  LABEL(4L, 0));
            TRI_ADD(&PC, idx_3,  idx_1,  idx_11, LABEL(5L, 0));
            TRI_ADD(&PC, idx_5,  idx_11, idx_1,  LABEL(6L, 0));
            TRI_ADD(&PC, idx_1,  idx_4,  idx_5,  LABEL(7L, 0));
            TRI_ADD(&PC, idx_0,  idx_5,  idx_4,  LABEL(8L, 0));
            TRI_ADD(&PC, idx_4,  idx_8,  idx_0,  LABEL(9L, 0));

            TRI_ADD(&PC, idx_2,  idx_0,  idx_8,  LABEL(10L, 0));
            TRI_ADD(&PC, idx_8,  idx_6,  idx_2,  LABEL(11L, 0));
            TRI_ADD(&PC, idx_7,  idx_2,  idx_6,  LABEL(12L, 0));
            TRI_ADD(&PC, idx_6,  idx_3,  idx_7,  LABEL(13L, 0));
            TRI_ADD(&PC, idx_11, idx_7,  idx_3,  LABEL(14L, 0));
            TRI_ADD(&PC, idx_11, idx_5,  idx_9,  LABEL(15L, 0));
            TRI_ADD(&PC, idx_5,  idx_0,  idx_9,  LABEL(16L, 0));
            TRI_ADD(&PC, idx_0,  idx_2,  idx_9,  LABEL(17L, 0));
            TRI_ADD(&PC, idx_9,  idx_2,  idx_7,  LABEL(18L, 0));
            TRI_ADD(&PC, idx_11, idx_9,  idx_7,  LABEL(19L, 0));

            break;
        }
        case OCTAHEDRON:
        {
            // octahedron defined by 6 points and 8 faces
            // four points around the centre
            vector_3D T(M_PI/8.0, 0.0, 0.0);

            int idx_0 = PC_ADD_ROT(-1.0,  1.0,  0.0, T);
            int idx_1 = PC_ADD_ROT( 1.0,  1.0,  0.0, T);
            int idx_2 = PC_ADD_ROT( 1.0, -1.0,  0.0, T);
            int idx_3 = PC_ADD_ROT(-1.0, -1.0,  0.0, T);
            // two poles
            int idx_4 = PC_ADD_ROT( 0.0,  0.0,  1.0, T);
            int idx_5 = PC_ADD_ROT( 0.0,  0.0, -1.0, T);

            // create the faces
            TRI_ADD(&PC, idx_4,  idx_0,  idx_1, LABEL(0L, 0));
            TRI_ADD(&PC, idx_5,  idx_0,  idx_1, LABEL(1L, 0));
            TRI_ADD(&PC, idx_4,  idx_1,  idx_2, LABEL(2L, 0));
            TRI_ADD(&PC, idx_5,  idx_1,  idx_2, LABEL(3L, 0));
            TRI_ADD(&PC, idx_4,  idx_2,  idx_3, LABEL(4L, 0));
            TRI_ADD(&PC, idx_5,  idx_2,  idx_3, LABEL(5L, 0));
            TRI_ADD(&PC, idx_4,  idx_3,  idx_0, LABEL(6L, 0));
            TRI_ADD(&PC, idx_5,  idx_3,  idx_0, LABEL(7L, 0));

            break;
        }
        case DYMAXION:
        {
            // Buckminster Fuller's Dymaxion (tm) map.  Icosahedron but with specific
            // coordinates.  Taken from:
            // Gray, Robert W., Exact Transformation Equations For Fuller's World Map, Cartographica, 32(3): 17-25, 1995.

            int idx_1  = PC_ADD( 0.420152426708710003,  0.078145249402782959,  0.904082550615019298);
            int idx_2  = PC_ADD( 0.995009439436241649, -0.091347795276427931,  0.040147175877166645);
            int idx_3  = PC_ADD( 0.518836730327364437,  0.835420380378235850,  0.181331837557262454);
            int idx_4  = PC_ADD(-0.414682225320335218,  0.655962405434800777,  0.630675807891475371);
            int idx_5  = PC_ADD(-0.515455959944041808, -0.381716898287133011,  0.767200992517747538);
            int idx_6  = PC_ADD( 0.355781402532944713, -0.843580002466178147,  0.402234226602925571);
            int idx_7  = PC_ADD( 0.414682225320335218, -0.655962405434800777, -0.630675807891475371);
            int idx_8  = PC_ADD( 0.515455959944041808,  0.381716898287133011, -0.767200992517747538);
            int idx_9  = PC_ADD(-0.355781402532944713,  0.843580002466178147, -0.402234226602925571);
            int idx_10 = PC_ADD(-0.995009439436241649,  0.091347795276427931, -0.040147175877166645);
            int idx_11 = PC_ADD(-0.518836730327364437, -0.835420380378235850, -0.181331837557262454);
            int idx_12 = PC_ADD(-0.420152426708710003, -0.078145249402782959, -0.904082550615019298);
            
            // create the faces
            TRI_ADD(&PC, idx_1,  idx_2,  idx_3,  LABEL(0L, 0));
            TRI_ADD(&PC, idx_1,  idx_3,  idx_4,  LABEL(1L, 0));
            TRI_ADD(&PC, idx_1,  idx_4,  idx_5,  LABEL(2L, 0));
            TRI_ADD(&PC, idx_1,  idx_5,  idx_6,  LABEL(3L, 0));
            TRI_ADD(&PC, idx_1,  idx_6,  idx_2,  LABEL(4L, 0)); // reversed
            TRI_ADD(&PC, idx_2,  idx_8,  idx_3,  LABEL(5L, 0)); // reversed
            TRI_ADD(&PC, idx_8,  idx_9,  idx_3,  LABEL(6L, 0)); // reversed
            TRI_ADD(&PC, idx_9,  idx_4,  idx_3,  LABEL(7L, 0)); // reversed
            TRI_ADD(&PC, idx_10, idx_4,  idx_9,  LABEL(8L, 0)); // reversed
            TRI_ADD(&PC, idx_5,  idx_4,  idx_10, LABEL(9L, 0)); // reversed

            TRI_ADD(&PC, idx_5,  idx_10, idx_11, LABEL(10L, 0)); // reversed
            TRI_ADD(&PC, idx_5,  idx_11, idx_6,  LABEL(11L, 0)); // reversed
            TRI_ADD(&PC, idx_11, idx_7,  idx_6,  LABEL(12L, 0)); // reversed
            TRI_ADD(&PC, idx_7,  idx_2,  idx_6,  LABEL(13L, 0)); // reversed
            TRI_ADD(&PC, idx_8,  idx_2,  idx_7,  LABEL(14L, 0)); // reversed
            TRI_ADD(&PC, idx_12, idx_9,  idx_8,  LABEL(15L, 0));
            TRI_ADD(&PC, idx_12, idx_10, idx_9,  LABEL(16L, 0)); // reversed
            TRI_ADD(&PC, idx_12, idx_11, idx_10, LABEL(17L, 0));
            TRI_ADD(&PC, idx_12, idx_7,  idx_11, LABEL(18L, 0)); // reversed
            TRI_ADD(&PC, idx_12, idx_8,  idx_7,  LABEL(19L, 0));
            break;
        }
    }
}

/*****************************************************************************/

void tri_grid::assign_points_from_grid(ncdata* nc_input_data, int perim)
{
    // relative coordinates for the 5 points of a grid box
    const int n_gp = 5;
    const FP_TYPE x_pts[5] = {0,-1,1,1,-1};
    const FP_TYPE y_pts[5] = {0,-1,-1,1,1};
    FP_TYPE lon_d = nc_input_data->get_lon_d();
    FP_TYPE lat_d = nc_input_data->get_lat_d();
    FP_TYPE lon_s = nc_input_data->get_lon_s();
    FP_TYPE lat_s = nc_input_data->get_lat_s();
    
    FP_TYPE rot_pole_lon, rot_pole_lat;
    if (nc_input_data->has_rotated_grid())
    {
         rot_pole_lon = nc_input_data->get_rotated_grid()->get_rotated_pole_longitude();
         rot_pole_lat = nc_input_data->get_rotated_grid()->get_rotated_pole_latitude();
    }
    
    // get the coordinates from the ncdata file and add to the quad trees
    for (int j=perim; j<nc_input_data->get_lat_len()-perim; j++)
    {
        for (int i=perim; i<nc_input_data->get_lon_len()-perim; i++)
        {
            FP_TYPE lon_b = lon_s + i*lon_d;
            FP_TYPE lat_b = lat_s + j*lat_d;
            for (int k=0; k<n_gp; k++)
            {
                // get the lon as the corners of a grid box
                FP_TYPE lon = lon_b + x_pts[k] * lon_d * 0.49999;
                FP_TYPE lat = lat_b + y_pts[k] * lat_d * 0.49999;
                // extra conversion for rotated grid
                if (nc_input_data->has_rotated_grid())
                {
                    FP_TYPE glob_lon, glob_lat;
                    Rot2Global(lat, lon, rot_pole_lat, rot_pole_lon, 
                               glob_lat, glob_lon);
                    lat = glob_lat;
                    lon = glob_lon;
                }
                // calculate area of triangle
                FP_TYPE A_gb = grid_box_area(lat_b, lon_d, lat_d);
                // convert to Cartesian coordinates
                vector_3D cart_coords = model_to_cart(lon, lat);
                // determine which triangle (head node) this point is in
                for (unsigned int tri=0; tri<triangles.size(); tri++)
                {
                    if (point_in_tri(&cart_coords, triangles[tri]->get_root()->get_data()))
                    {
                        // calculate triangle surface area
                        FP_TYPE A_T = triangles[tri]->get_root()->get_data()->surface_area();
                        FP_TYPE W = A_gb / A_T;
                        triangles[tri]->get_root()->get_data()->add_index(i,j, cart_coords, W);
                    }
                }
            }
        }
    }
}

/*****************************************************************************/

void tri_grid::distribute_grid_indices_to_children(QT_TRI_NODE* triangle, ncdata* nc_input_data)
{
    // check whether the triangles should contain one of the indexed points
    const std::list<grid_index>* grid_index_list = triangle->get_data()->get_grid_indices();
    
    // corners and mid-point of grid box
    const int n_gp = 5;
    const FP_TYPE x_pts[5] = {0,-1,1,1,-1};
    const FP_TYPE y_pts[5] = {0,-1,-1,1,1};

    // grid spacing info
    FP_TYPE lon_d = nc_input_data->get_lon_d();
    FP_TYPE lat_d = nc_input_data->get_lat_d();
    FP_TYPE lon_s = nc_input_data->get_lon_s();
    FP_TYPE lat_s = nc_input_data->get_lat_s();
    
    // check for rotated grid
    FP_TYPE rot_pole_lon, rot_pole_lat;
    if (nc_input_data->has_rotated_grid())
    {
         rot_pole_lon = nc_input_data->get_rotated_grid()->get_rotated_pole_longitude();
         rot_pole_lat = nc_input_data->get_rotated_grid()->get_rotated_pole_latitude();
    }

    // loop over all grid points
    for (std::list<grid_index>::const_iterator gii = grid_index_list->begin();
         gii != grid_index_list->end(); gii++)
    {
        // calculate the centre point of the grid box
        FP_TYPE lon_b = lon_s + gii->i*lon_d;
        FP_TYPE lat_b = lat_s + gii->j*lat_d;
        // area of the grid box
        FP_TYPE A_gb = grid_box_area(lat_b, lon_d, lat_d);

        // get the actual triangle from the quad tree
        // check for point inclusion of current grid index
        bool tri_found = false;
        for (int i=0; i<4; i++)
        {
            // need to check each child triangle
            indexed_force_tri_3D* child_tri = triangle->get_child(i)->get_data();
            // get the surface area of the triangle
            FP_TYPE A_T = child_tri->surface_area();
            // pit = point in triangle, false at the moment
            int n_pit = 0;
            // loop over the centre and 4 corner points of the grid box
            for (int k=0; k<n_gp; k++)
            {
                // get the lon and lat as the corners of a grid box
                FP_TYPE lon = lon_b + x_pts[k] * lon_d * 0.49999;
                FP_TYPE lat = lat_b + y_pts[k] * lat_d * 0.49999;
                // extra conversion for rotated grid
                if (nc_input_data->has_rotated_grid())
                {
                    FP_TYPE glob_lon, glob_lat;
                    Rot2Global(lat, lon, rot_pole_lat, rot_pole_lon, 
                               glob_lat, glob_lon);
                    lat = glob_lat;
                    lon = glob_lon;
                }
                // convert to Cartesian coordinates
                vector_3D cart_coords = model_to_cart(lon, lat);

                // check to see whether either the grid box centre, or each grid box corner
                // is in the triangle
                if (point_in_tri(&(cart_coords), child_tri))
                    n_pit+=1;
            }
            if (n_pit > 0)
            {
                // calculate the weight
                FP_TYPE W = n_pit*0.2*A_gb/A_T;
                child_tri->add_index(gii->i, gii->j, gii->cart_coord, W);
                tri_found = true;
            }
        }
        // if not added then add to nearest
        FP_TYPE min_dist = 1e10;
        int min_tri_idx = -1;
        if (not tri_found)
        {
            for (int i=0; i<4; i++)
            {
                // calculate distance between point and centroid
                indexed_force_tri_3D* child_tri = triangle->get_child(i)->get_data();
                FP_TYPE dist = (child_tri->centroid() - gii->cart_coord).mag();
                if (dist < min_dist)
                    min_tri_idx = i;
            }
            // add to the child tri the point is nearest to
            triangle->get_child(min_tri_idx)->get_data()->add_index(gii->i, gii->j, gii->cart_coord, 1.0);
        }
    }
}

/*****************************************************************************/

void tri_grid::split_triangle(QT_TRI_NODE* triangle, int max_levs, ncdata* nc_input_data)
{
    // split the triangle indicated by PTI (parent triangle index) into four
    // new triangles and insert the 3 new force points into the point cloud 
    // and the triangle definition into the HTMA containing the references to 
    // the force points
    
    // if the number of indices is zero then do not split
    if (triangle->get_data()->get_number_of_indices() == 0)
        return;
    
    int iM[3];      // indices to new points in the point cloud
    for (int i=0; i<3; i++)
    {
        int j = i<2 ? i+1 : 0;
        // get the two points of the current triangle edge and create a mid point
        force_point MP = (PC[(*(triangle->get_data()))[i]] + PC[(*(triangle->get_data()))[j]]) * 0.5;
        // normalise the mid point
        MP *= 1.0 / MP.mag();
        // add to the point cloud
        iM[i] = PC_ADD_SPLIT(MP);
    }
    // new points have been created - now add the four new triangles
    indexed_force_tri_3D* ift = triangle->get_data();
    // get the current label
    LABEL c_label = ift->get_label();
    // calculate the multiplier
    long int mp = (long int)(pow(10, c_label.max_level+2));
    triangle->add_child(indexed_force_tri_3D(&PC, (*ift)[0],     iM[0], iM[2], LABEL(c_label.label+1*mp, c_label.max_level+1)));
    triangle->add_child(indexed_force_tri_3D(&PC,     iM[0], (*ift)[1], iM[1], LABEL(c_label.label+2*mp, c_label.max_level+1)));
    triangle->add_child(indexed_force_tri_3D(&PC,     iM[1], (*ift)[2], iM[2], LABEL(c_label.label+3*mp, c_label.max_level+1)));
    triangle->add_child(indexed_force_tri_3D(&PC,     iM[0],     iM[1], iM[2], LABEL(c_label.label+4*mp, c_label.max_level+1)));
    
    distribute_grid_indices_to_children(triangle, nc_input_data);
    
    // check whether each child should be split again
    for (int i=0; i<4; i++)
    {
        if (triangle->get_level() < (max_levs-1))
            // recursive split
            split_triangle(triangle->get_child(i), max_levs, nc_input_data);
    }
}

/*****************************************************************************/

void tri_grid::split_triangles(int max_levs, ncdata* nc_input_data)
{
    std::cout << "# Splitting triangles" << std::endl;
    // loop through each base triangle and split until the level of the triangle
    // is less than the max level
    for (unsigned int tri=0; tri<triangles.size(); tri++)
    {
        // check first whether it should be split or not
        if (triangles[tri]->get_root()->get_data()->get_number_of_indices() > 0)
            split_triangle(triangles[tri]->get_root(), max_levs, nc_input_data);  // NB - this function is recursive
    }
}

/*****************************************************************************/

void tri_grid::build_adjacency_maps(void)
{
    std::cout << "# Building adjacency map level: ";
    for (int lev=0; lev<get_max_level(); lev++)
    {
        std::list<QT_TRI_NODE* > global_tri_list = get_triangles_at_level(lev);
        // now loop over each triangle and consider each other triangle as well
        std::list<QT_TRI_NODE* >::iterator it_1;
        std::list<QT_TRI_NODE* >::iterator it_2;
        std::cout << lev;
        for (it_1 = global_tri_list.begin(); it_1 != global_tri_list.end(); it_1++)
        {
            for (it_2 = global_tri_list.begin(); it_2 != global_tri_list.end(); it_2++)
            {
                int n_shared_pts = 0;           // number of shared points
                if (it_1 == it_2)   // don't check against itself
                    continue;           
                // check each point in the triangle reference by it_index_1 against
                // each point in the triangle it_index_2
                for (int i=0; i<3; i++)
                    for (int j=0; j<3; j++)
                        if ((*(*it_1)->get_data())[i] == (*(*it_2)->get_data())[j])
                            n_shared_pts ++;    // sufficient just to check indices match
                if (n_shared_pts >= 1)  // point adjacency              
                    (*it_1)->get_data()->add_adjacency((*it_2)->get_data()->get_label(), POINT);
                if (n_shared_pts >= 2)  // edge adjacency
                    (*it_1)->get_data()->add_adjacency((*it_2)->get_data()->get_label(), EDGE);
            }
        }
        if (lev < 10)
            std::cout << "\b";
        if (lev >= 10 and lev < 100)
            std::cout << "\b\b";
        if (lev >= 100 and lev < 1000)
            std::cout <<"\b\b\b";
    }
    std::cout << std::endl;
}

/*****************************************************************************/

void build_ds_index_node(QT_TRI_NODE* current_node, int& current_index)
{
    // set the datastore indices into a simple array for fast and efficient data storage
    current_node->get_data()->set_ds_index(current_index);
    current_index ++;
    for (int i=0; i<4; i++)
    {
        if (current_node->get_child(i) != NULL)
            build_ds_index_node(current_node->get_child(i), current_index);
    }
}

/*****************************************************************************/

void tri_grid::build_ds_indices(void)
{
    std::cout << "# Building datastore indices" << std::endl;
    int current_index = 0;
    for (std::vector<QT_TRI* >::iterator it = triangles.begin(); 
         it != triangles.end(); it++)
        build_ds_index_node((*it)->get_root(), current_index);
}

/*****************************************************************************/

int tri_grid::get_number_of_base_tris(void)
{
    return triangles.size();
}

/*****************************************************************************/

QT_TRI_NODE* tri_grid::get_base_tri(int n)
{
    return triangles[n]->get_root();
}
/*****************************************************************************/

std::list<QT_TRI_NODE*> tri_grid::get_triangles_at_level(int level)
{
    // get the first triangle to create a list
    std::list<QT_TRI_NODE*> tri_node_list;
    // now extend the list to contain the extra triangles
    for (unsigned int i=0; i<triangles.size(); i++)
    {
        std::list<QT_TRI_NODE*> this_tri_list = triangles[i]->get_all_nodes_at_level(level);
        tri_node_list.splice(tri_node_list.end(), this_tri_list);
    }
    return  tri_node_list;
}

/*****************************************************************************/
        
void save_node(std::ofstream& out, QT_TRI_NODE* current_node)
{
    current_node->get_data()->save(out);
    // save the children recursively
    for (int i=0; i<4; i++)
        if (current_node->get_child(i) != NULL)
            save_node(out, current_node->get_child(i));
}

/*****************************************************************************/

void tri_grid::save(std::string filename)
{
    std::cout << "# Saving mesh" << std::endl;
    std::ofstream out;
    out.open(filename.c_str(), std::ios::out | std::ios::binary);
    if (!out)
        throw(std::string("Saving mesh.  File could not be opened or written to: " + filename));
    // write out the meta data
    write_meta_data(out, meta_data);
    // write out the point cloud
    point_cloud_instance.save(out);
    // write out number of triangle roots
    write_int(out, triangles.size());
    // write out all the triangles
    for (unsigned int tri=0; tri<triangles.size(); tri++)
        save_node(out, triangles[tri]->get_root());
    out.close();
}

/*****************************************************************************/

indexed_force_tri_3D* load_node(std::ifstream& in, QT_TRI_NODE* current_node, 
                                indexed_force_tri_3D* current_tri,
                                point_cloud* pPC)
{
    // assign the current data
    (*current_node->get_data()) = *current_tri;
    // get the next triangle in the file
    indexed_force_tri_3D* next_tri = new indexed_force_tri_3D;
    next_tri->load(in, pPC);
    // check the current triangle label against the previous triangle label
    if (next_tri->get_label().size() == current_tri->get_label().size() || in.eof())
        // occur at same level so do not create any children
        return next_tri;
    if (next_tri->get_label().size() > current_tri->get_label().size())
    {
        // do not occur at same level so add children and start the recursion
        for (int i=0; i<4; i++)
        {
            // add child and assign the next triangle to it
            current_node = current_node->add_child(indexed_force_tri_3D());
            next_tri = load_node(in, current_node, next_tri, pPC);
            // need to go back up to the parent node so as to create the siblings correctly
            // otherwise we would be reproducing with our own siblings
            if (current_node->get_parent() != NULL)
                current_node = current_node->get_parent();
        }
    }
    // control will never reach here but will return next_tri anyway to avoid warnings
    return next_tri;
}

/*****************************************************************************/

void tri_grid::load(std::string filename)
{
    std::cout << "# Loading mesh" << std::endl;
    // do this as a text file first
    std::ifstream in;
    in.open(filename.c_str(), std::ios::in | std::ios::binary);
    if (!in)
        throw(std::string("Loading mesh.  File could not be opened ") + filename);
    meta_data = read_meta_data(in);
    // read in the point cloud
    point_cloud_instance.load(in);
    // read in the number of triangle roots
    int n_tris = read_int(in);
    // create the quad trees for the triangle roots
    // start with the first triangle
    indexed_force_tri_3D* current_tri = new indexed_force_tri_3D;
    current_tri->load(in, &point_cloud_instance);
    // now loop through the rest of the triangles
    for (int i=0; i<(int)(n_tris); i++)
    {
        triangles.push_back(new quadtree<indexed_force_tri_3D>);
        // get the current node as the first root node
        QT_TRI_NODE* current_node = triangles[i]->get_root();
        current_tri = load_node(in, current_node, current_tri, &point_cloud_instance);
    }
    in.close();
}

/*****************************************************************************/

void tri_grid::save_text(std::string filename)
{
    // write out the latitude and longitude coordinates for the max_level
    // open the file
    std::ofstream out;
    out.open(filename.c_str(), std::ios::out);
    if (!out)
        throw(std::string("Saving lat / lon mesh.  File could not be opened or written to: " + filename));
    // loop through the triangles
    std::list<QT_TRI_NODE* > tri_list = get_triangles_at_level(max_lev_input);
    std::list<QT_TRI_NODE* >::iterator it_tri_list;
    for (it_tri_list = tri_list.begin(); it_tri_list != tri_list.end(); it_tri_list++)
    {
        // write the label, followed by the three lon / lat pairs
        indexed_force_tri_3D* c_tri = (*it_tri_list)->get_data();
        // label
        out << c_tri->get_label().label << " ";
        // lon / lat pairs for each vertex point
        for (int v=0; v<3; v++)
        {
            // get the current vertex
            force_point& c_vert = point_cloud_instance[(*c_tri)[v]];
            // transform to lon / lat
            FP_TYPE lon, lat;
            cart_to_model(c_vert, lon, lat);
            // write
            out << lon << " " << lat << " ";
        }
        // write the data store grid index
        out << c_tri->get_ds_index() << " ";
        // write source grid indices
        const std::list<grid_index>* g_idxs = c_tri->get_grid_indices();
        std::list<grid_index>::const_iterator g_it;
        for (g_it = g_idxs->begin(); g_it != g_idxs->end(); g_it++)
        {
            out << g_it->i << " " << g_it->j << " ";
        }
        out << std::endl;
    }
    out.close();
    std::cout << "# Saved lon / lat grid for level " << max_lev_input << " to file: " << filename << std::endl;
}

/*****************************************************************************/

LABEL tri_grid::get_triangle_for_point(vector_3D* P, int level)
{
    // determine which of the 20 base triangles the point is in
    int base_tri = -1;
    for (unsigned int tri=0; tri<triangles.size(); tri++)
        if (point_in_tri(P, triangles[tri]->get_root()->get_data()))
        {
            base_tri = tri;
            break;
        }
    // not found so find by distance
    if (base_tri == -1)
    {
        FP_TYPE min_dist = 2e20;
        int min_base = -1;
        FP_TYPE P_lon, P_lat;
        cart_to_model(*P, P_lon, P_lat);
        for (unsigned int tri=0; tri<triangles.size(); tri++)
        {
            FP_TYPE T_lon, T_lat;
            cart_to_model(triangles[tri]->get_root()->get_data()->centroid(), T_lon, T_lat);
            FP_TYPE dist = haversine(P_lon, P_lat, T_lon, T_lat, 1.0);
            if (dist < min_dist)
            {
                dist = min_dist;
                min_base = tri;
            }
        }   
        base_tri = min_base;
    }
    // get the base node
    QT_TRI_NODE* current_node = get_base_tri(base_tri);
    bool found = false;
    int current_level = 0;
    while (not found)
    {
        if (current_node->is_leaf() or current_level == level)
        {
            found = true;
            break;
        }
        else
        {
            bool found_in_child = false;
            current_level += 1;
            for (int i=0; i<4; i++)
            {
                QT_TRI_NODE* current_child = current_node->get_child(i);
                if (point_in_tri(P, current_child->get_data()))
                {
                    current_node = current_child;
                    found_in_child = true;
                    break;
                }
            }
            // not found by point inclusion - find by minimum distance
            if (not found_in_child)
            {
                FP_TYPE min_dist = 2e20;
                int min_child = -1;
                FP_TYPE P_lon, P_lat;
                cart_to_model(*P, P_lon, P_lat);
                for (int i=0; i<4; i++)
                {
                    FP_TYPE T_lon, T_lat;
                    cart_to_model(current_node->get_child(i)->get_data()->centroid(), T_lon, T_lat);
                    FP_TYPE dist = haversine(P_lon, P_lat, T_lon, T_lat, 1.0);
                    if (dist < min_dist)
                    {
                        dist = min_dist;
                        min_child = i;
                    }
                }
                current_node = current_node->get_child(min_child);
            }
        }
    }
    return current_node->get_data()->get_label();
}

/*****************************************************************************/

FP_TYPE tri_grid::distance_between_triangles(LABEL SL, LABEL EL)
{
    // get each triangle
    indexed_force_tri_3D* tri_SL = get_triangle(SL);
    indexed_force_tri_3D* tri_EL = get_triangle(EL);
    // get the centroid
    vector_3D vec_SL = tri_SL->centroid();
    vector_3D vec_EL = tri_EL->centroid();
    // convert the points to model coordinates (lat-lon)
    FP_TYPE lon_SL, lat_SL;
    FP_TYPE lon_EL, lat_EL;
    cart_to_model(vec_SL, lon_SL, lat_SL);
    cart_to_model(vec_EL, lon_EL, lat_EL);
    // get the distance
    FP_TYPE dist = haversine(lon_SL, lat_SL, lon_EL, lat_EL, EARTH_R);
    return dist;
}

/*****************************************************************************/

LABEL_STORE tri_grid::get_path(LABEL SL, LABEL EL, int level, int resolution)
{
    // create a path of labels between the starting label (SL) and the 
    // ending label (EL)
    // In the feature recognition routines, SL is the triangle being tested for
    // inclusion (candidate triangle - CL ) and EL is the original triangle at 
    // the centre of the feature (OL)

    // the path is derived by taking a line between the SL and EL and interpolating
    // along the line at points equal to the minimum of the distance between SL and
    // its neighbours.  The point inclusion test is then used to determine which
    // triangle each point is in
    
    // This algorithm has the benefit of terminating always - other path finding
    // methods I have tried, based on the adjacency list, have not!

    // if the level is -1 then set it to be the maximum level
    if (level == -1) level = get_max_level();

    // what is the minimum distance between SL and it's adjacent labels?
    indexed_force_tri_3D* s_tri = get_triangle(SL);
    LABEL_STORE path;
    FP_TYPE min_dist = 2e20;
    const LABEL_STORE* s_adj_labs = s_tri->get_adjacent_labels(POINT);
    for (LABEL_STORE::const_iterator it_adj_lab = s_adj_labs->begin();
         it_adj_lab != s_adj_labs->end(); it_adj_lab++)
    {
        indexed_force_tri_3D* a_tri = get_triangle(*it_adj_lab);
        vector_3D V = s_tri->centroid() - a_tri->centroid();
        FP_TYPE dist = V.mag();
        if (dist < min_dist)
            min_dist = dist;
    }
    
    // form the SL to EL vector
    indexed_force_tri_3D* e_tri = get_triangle(EL);
    vector_3D se_vec = s_tri->centroid() - e_tri->centroid();
    // loop until we reach the EL triangle  
    FP_TYPE n_pts = se_vec.mag()/(min_dist*resolution); // how many points ?    
    vector_3D add_vec = se_vec / n_pts;     // what to add each iteration (vector)
    for (int i=0; i<int(n_pts+0.5)+1; i++)
    {
        vector_3D c_vec = s_tri->centroid() - add_vec * i;
        // normalise
        c_vec = c_vec / c_vec.mag();
        // point inclusion
        LABEL tri_label = get_triangle_for_point(&c_vec, level);
        // only add if not already added
        if (std::find(path.begin(), path.end(), tri_label) == path.end())
            path.push_back(tri_label);
    }
    return path;
}
