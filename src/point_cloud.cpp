/******************************************************************************
** Program : point_cloud.h
** Author  : Neil Massey
** Date    : 24/08/11
** Purpose : class to hold all the force points used in the triangular mesh
******************************************************************************/

#include "point_cloud.h"
#include <iostream>
#include "bin_file_utils.h"
#define PCT std::vector<force_point>

/*****************************************************************************/

point_cloud::point_cloud(void)
{
}

/*****************************************************************************/

force_point& point_cloud::operator[](int i)
{
	return point_cloud_store[i];
}

/*****************************************************************************/

int point_cloud::add_point(force_point FP, bool add_dups)
{
	// by not checking whether duplicates are added, the loading of the
	// mesh can be massively quickened

	if (!add_dups)
	{
		// see if the point already exists in the point cloud
		for (unsigned int i=0; i<point_cloud_store.size(); i++)
			if (FP == point_cloud_store[i])
				return i;
	}

	// not found so add to the back of the vector and return the last index
	point_cloud_store.push_back(FP);
	return point_cloud_store.size() - 1;
}

/*****************************************************************************/

void point_cloud::equalise(int max_its, int n_levels)
{
	PCT::iterator it_p;
	PCT::iterator it_q;
	// distance, force, radius vector
    FP_TYPE d, f;
	vector_3D r;

	std::cout << "# Equalising triangles, number of triangle points: " 
			  << point_cloud_store.size() << std::endl;
	std::cout << "# Iteration number: ";

	FP_TYPE F;
	// prevent too much adjustment per iteration
	if (n_levels != 0)
		F = 1.0 / n_levels;
	else
		F = 1.0;

	int c = 0;
	// equalise the distance between points in the point cloud
	while (c < max_its)
	{
		c++;
		std::cout << c;
		std::cout.flush();
		// apply the forces from one point to the next
		for (it_p = point_cloud_store.begin(); it_p != point_cloud_store.end(); it_p++)
		{
			for (it_q = point_cloud_store.begin(); it_q != point_cloud_store.end(); it_q++)
			{
				r = *it_p - *it_q;
				d = r.mag();
				// check they're not the same point
				if (d > 0.0)
				{
					// repel the points apart by adding a force vector
					// to each point proportional to 1/2 the distance from
					// each other point
					f = F*0.5/d;
					it_p->add_force(r*f);
					it_q->add_force(r*-f);
				}
			}
			it_p->apply_forces();
		}
		int e = c;
		while (e > 0)
		{
			e = e / 10;
			std::cout << "\b";
		}
	}
	std::cout << std::endl;
}

/*****************************************************************************/

void point_cloud::save(std::ofstream& out)
{
	// save the size of the point cloud
	write_int(out, point_cloud_store.size());
	// save the point cloud to the (binary) stream
	for (PCT::iterator it_p = point_cloud_store.begin();
		 it_p != point_cloud_store.end(); it_p++)
		write_vector(out, *it_p);
}

/*****************************************************************************/

void point_cloud::load(std::ifstream& in)
{
	// read the size of the point cloud
	int size = read_int(in);
	// read <size> number of vectors and add to the point cloud
	point_cloud_store.resize(size);
	for (int i=0; i<size; i++)
	{
		vector_3D pt = read_vector(in);
		point_cloud_store[i] = pt;
	}
}
