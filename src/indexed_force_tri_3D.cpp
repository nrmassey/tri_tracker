/******************************************************************************
** Program : indexed_force_tri_3D.cpp
** Author  : Neil Massey
** Date    : 11/06/13
** Purpose : class to hold 3D triangle information with reference to points in
**           the point cloud and reference to indices in the original grid
******************************************************************************/

#include "indexed_force_tri_3D.h"
#include "bin_file_utils.h"

/*****************************************************************************/

indexed_force_tri_3D::indexed_force_tri_3D(void) : force_tri_3D()
{
}

/*****************************************************************************/

indexed_force_tri_3D::indexed_force_tri_3D(point_cloud* pPC, int ip0, int ip1, int ip2, LABEL i_label)
					 :force_tri_3D(pPC, ip0, ip1, ip2), label(i_label)
{
}

/*****************************************************************************/

indexed_force_tri_3D::indexed_force_tri_3D(const force_tri_3D& rhs)
                     :force_tri_3D(rhs)
{
}

/*****************************************************************************/

void indexed_force_tri_3D::add_index(int ii, int ij, vector_3D cart_coord)
{
	// add a grid index to this triangle
	grid_index gdx;
	gdx.i = ii;
	gdx.j = ij;
	gdx.cart_coord = cart_coord;
	grid_indices.push_back(gdx);
}

/*****************************************************************************/

void indexed_force_tri_3D::set_label(LABEL i_label)
{
	label = i_label;
}

/*****************************************************************************/

LABEL indexed_force_tri_3D::get_label(void)
{
	return label;
}

/*****************************************************************************/

int indexed_force_tri_3D::get_number_of_indices(void)
{
	return grid_indices.size();
}

/*****************************************************************************/

const std::list<grid_index>* indexed_force_tri_3D::get_grid_indices(void) const
{
	return &grid_indices;
}

/*****************************************************************************/

void indexed_force_tri_3D::add_adjacency(LABEL label, ADJACENCY type)
{
	adjacency[type].push_back(label);
}

/*****************************************************************************/

const LABEL_STORE* indexed_force_tri_3D::get_adjacent_labels(ADJACENCY type) const
{
	return &(adjacency[type]);
}

/*****************************************************************************/

int indexed_force_tri_3D::get_ds_index(void)
{
	return ds_index;
}

/*****************************************************************************/

void indexed_force_tri_3D::set_ds_index(int ti)
{
	ds_index = ti;
}
	
/*****************************************************************************/

void indexed_force_tri_3D::save(std::ofstream& out)
{
	// write the label string
	write_label(out, label);
	// write the point cloud indices for each vertex
	for (int t=0; t<3; t++)
		write_int(out, p[t]);
	// write the target index
	write_int(out, ds_index);
	// write the index list - first the length
	write_int(out, grid_indices.size());
	// now write all the indices
	for (std::list<grid_index>::iterator it=grid_indices.begin();
		 it!=grid_indices.end(); it++)
	{
		write_int(out, it->i);
		write_int(out, it->j);
		write_vector(out, it->cart_coord);
	}
	// write the point and edge adjacency list
	for (int a=0; a<2; a++)
	{
		write_int(out, adjacency[a].size());
		for (LABEL_STORE::iterator it=adjacency[a].begin();
			 it!= adjacency[a].end(); it++)
			write_label(out, *it);
	}
}

/*****************************************************************************/

void indexed_force_tri_3D::load(std::ifstream& in, point_cloud* pPC)
{
	// read the label in and set it
	LABEL label = read_label(in);
	set_label(label);
	if (in.eof())
		return;

	// set the point cloud
	ppoint_cloud_instance = pPC;
	// read the point cloud indices for each vertex
	for (int t=0; t<3; t++)
		p[t] = read_int(in);
	// calculate the centroid
	calculate_centroid();
	// read the target index
	ds_index = read_int(in);
	// read the length of index list in
	int n_idx = read_int(in);
	for (int i=0; i<n_idx; i++)
	{
		grid_index gr_idx;
		gr_idx.i = read_int(in);
		gr_idx.j = read_int(in);
		gr_idx.cart_coord = read_vector(in);
		grid_indices.push_back(gr_idx);
	}
	// read the point and adjacency list
	for (int a=0; a<2; a++)
	{
		// get the size first
		int s = read_int(in);
		for (int i=0; i<s; i++)
			adjacency[a].push_back(read_label(in));
	}
}