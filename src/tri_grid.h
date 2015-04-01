/******************************************************************************
** Program : tri_grid.h
** Author  : Neil Massey
** Date    : 10/06/13
** Purpose : class to represent a triangular grid for the region
******************************************************************************/

#ifndef TRI_GRID_H
#define TRI_GRID_H

#include "indexed_force_tri_3D.h"
#include "point_cloud.h"
#include "quadtree.h"
#include "meta_data.h"

enum SHAPE{ICOSAHEDRON=0, OCTAHEDRON, DYMAXION};
#define QT_TRI quadtree<indexed_force_tri_3D>
#define QT_TRI_NODE qt_node<indexed_force_tri_3D>

class tri_grid
{
	public:
		tri_grid(void);
		~tri_grid(void);
		
		// set up initial triangles
		void initialize(SHAPE initial_shape, class ncdata* nc_input_data,
						int max_levs, int max_its, int perim);
		indexed_force_tri_3D* get_triangle(LABEL label);
		QT_TRI_NODE* get_triangle_node(LABEL label);
		int get_max_level(void);
		int get_max_ds_index(void);	// for the size of the regridded data array
		// get the number of base triangles
		int get_number_of_base_tris(void);
		// get the root node of a base triangle
		QT_TRI_NODE* get_base_tri(int n);
		// get all triangles at a particular level
		std::list<QT_TRI_NODE*> get_triangles_at_level(int level);
		// get a path between two triangles as a vector of labels
		LABEL_STORE get_path(LABEL SL, LABEL EL, int resolution=1);
		// get the label of the triangle which the point lies in
		LABEL get_triangle_for_point(vector_3D* P);
		// get the distance between two triangles in the tri-grid
		FP_TYPE distance_between_triangles(LABEL SL, LABEL EL);
		// get the metadata
		META_DATA_TYPE* get_meta_data(void);
		// input / output functions
		void save(std::string filename);
		// fast_load controls whether duplicates are allowed in the point cloud		
		void load(std::string filename);
		// save as lat/lon text file
		void save_text(std::string filename);
		
	private:
		/***********************************************************************/
		void create_shape(SHAPE initial_shape, FP_TYPE R);
		void assign_points_from_grid(ncdata* nc_input_data, int perim);
		void split_triangle(QT_TRI_NODE* triangle, int max_levs, ncdata* nc_input_data);
		void split_triangles(int max_levs, ncdata* nc_input_data);
		void build_adjacency_maps(void);
		void build_ds_indices(void);
		bool point_in_tri(const vector_3D* P, const force_tri_3D* T);
		void distribute_grid_indices_to_children(QT_TRI_NODE* triangle, ncdata* nc_input_data);
		
		/***********************************************************************/
		std::vector<QT_TRI* > triangles;
		point_cloud point_cloud_instance;
		int max_lev_input;
		META_DATA_TYPE meta_data;
		/***********************************************************************/
};

#endif