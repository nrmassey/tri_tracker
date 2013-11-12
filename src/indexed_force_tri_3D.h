/******************************************************************************
** Program : indexed_force_tri_3D.h
** Author  : Neil Massey
** Date    : 11/06/13
** Purpose : class to hold 3D triangle information with reference to points in
**           the point cloud and reference to indices in the original grid
******************************************************************************/

#ifndef INDEXED_FORCE_TRI_3D_H
#define INDEXED_FORCE_TRI_3D_H

#include <vector>
#include <list>
#include <fstream>
#include "force_tri_3D.h"

enum ADJACENCY{POINT=0, EDGE};

struct grid_index
{
	// a struct to store the grid index so that it can be stored in a list
	int i;
	int j;
	vector_3D cart_coord;
};

class LABEL
{
	public:
		LABEL(void) : label(-1), max_level(-1) {}
		LABEL(long int ilabel, int imax_level) : label(ilabel), max_level(imax_level) {}
		bool operator==(LABEL& rhs) {return label == rhs.label;}
		bool operator==(const LABEL& rhs) const {return label == rhs.label;}
		bool operator<(LABEL& rhs) {return label < rhs.label;}
		bool operator>(LABEL& rhs) {return label > rhs.label;}
		bool operator<(const LABEL& rhs) const {return label < rhs.label;}
		bool operator>(const LABEL& rhs) const {return label > rhs.label;}
		int size(void){ return max_level; }
		long int label;
		int max_level;
};

typedef std::vector<LABEL> LABEL_STORE;

class indexed_force_tri_3D : public force_tri_3D
{
	public:
		indexed_force_tri_3D(void);
		indexed_force_tri_3D(point_cloud* pPC, int ip0, int ip1, int ip2, LABEL label);
		indexed_force_tri_3D(const force_tri_3D& rhs);
		void add_index(int ii, int ij, vector_3D cart_coord);
		void set_label(LABEL label);
		LABEL get_label(void);
		const std::list<grid_index>* get_grid_indices(void) const;
		int  get_number_of_indices(void);
		void set_ds_index(int ti);
		int  get_ds_index(void);
		void add_adjacency(LABEL label, ADJACENCY type);
		const LABEL_STORE* get_adjacent_labels(ADJACENCY type) const;
		void save(std::ofstream& out);
		void load(std::ifstream& in, point_cloud* pPC);
		
	private:
		std::list<grid_index> grid_indices;	// indices into the grid - each triangle may
											// have more than one
		LABEL label;						// the labelling system is a quick way to navigate to a triangle 
											// and provide adjacency relationships.  it is a long int for
											// speed of comparisons
		LABEL_STORE adjacency[2];			// list of labels for adjacency maps.  Two maps - POINT and EDGE
		int ds_index;						// index into an array to store the regridded data
};

#endif