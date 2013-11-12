/******************************************************************************
** Program : point_cloud.h
** Author  : Neil Massey
** Date    : 24/08/11
** Purpose : class to hold all the force points used in the triangular mesh
******************************************************************************/

#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include "force_point.h"
#include <vector>

class point_cloud
{
	public:
		point_cloud(void);
		int add_point(force_point FP, bool add_dups=false);
		void equalise(int max_its, int n_levels=0);
        force_point& operator[](int i);
        int size(void) {return point_cloud_store.size();}
		void save(std::ofstream& out);
		void load(std::ifstream& in);

	private:
		std::vector<force_point> point_cloud_store;
};

#endif
