/******************************************************************************
** Program : force_tri_3D.h
** Author  : Neil Massey
** Date    : 23/08/11
** Purpose : class to hold 3D triangle information where the points are merely
**           indices to force points which occur in the point cloud
******************************************************************************/

#ifndef FORCE_TRI_3D_H
#define FORCE_TRI_3D_H

#include "point_cloud.h"
#include "vector_3D.h"

class force_tri_3D
{
	public:
		force_tri_3D(void);
		force_tri_3D(point_cloud* pPC, int ip0, int ip1, int ip2);
		force_tri_3D(const force_tri_3D& rhs);
		int& operator[](int i);
		const int& operator[](int i) const;
		inline point_cloud* get_point_cloud(void) const {return ppoint_cloud_instance;}

		bool operator==(const force_tri_3D& rhs) const;
		bool operator!=(const force_tri_3D& rhs) const;

		void set(point_cloud* pPC, int ip0, int ip1, int ip2);

		FP_TYPE surface_area(FP_TYPE radius=1.0) const;       // surface area of the triangle
		inline vector_3D centroid(void) const {return c0;}
		vector_3D plane_normal(void) const; // calculate the normal to the plane defined by the tri$

	protected:
		int p[3];							// 3 indices to the vertices
		vector_3D c0;						// centroid - pre calc when triangle is created
		point_cloud* ppoint_cloud_instance;	// pointer reference to point cloud
		void calculate_centroid(void);		//
};

#endif
