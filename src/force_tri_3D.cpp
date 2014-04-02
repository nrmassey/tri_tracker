/******************************************************************************
** Program : force_tri_3D.cpp
** Author  : Neil Massey
** Date    : 23/08/11
** Purpose : class to hold 3D triangle information with force point vectors
** Reference : Renka, R.J. 1984.  Interpolation of data on the surface of a sphere.  
**             ACM Trans. Math Softw. 10, 4 (Dec), 437-439.
******************************************************************************/

#include <assert.h>
#include <math.h>
#include "force_tri_3D.h"

#define PPC ppoint_cloud_instance

/*****************************************************************************/

force_tri_3D::force_tri_3D(void)
			 :PPC(NULL)
{
	p[0] = p[1] = p[2] = -1;
}

/*****************************************************************************/

force_tri_3D::force_tri_3D(point_cloud* pPC, int ip0, int ip1, int ip2)
{
	assert(pPC != NULL);
	set(pPC, ip0, ip1, ip2);
}

/*****************************************************************************/

force_tri_3D::force_tri_3D(const force_tri_3D& rhs)
{
	if (rhs.PPC != NULL)
		set(rhs.PPC, rhs.p[0], rhs.p[1], rhs.p[2]);
	else
		PPC = NULL;
}

/*****************************************************************************/

int& force_tri_3D::operator[](int i)
{
	assert(i >=0 && i < 3);
	return p[i];
}

/*****************************************************************************/

const int& force_tri_3D::operator[](int i) const
{
	assert(i >=0 && i < 3);
	return p[i];
}

/*****************************************************************************/

void force_tri_3D::calculate_centroid(void)
{
	// calculate the centroid of a 3D triangle
	const FP_TYPE THIRD = 1.0/3.0;
	if (p[0] != -1 && p[1] != -1 && p[2] != -1)
	{
		c0[0] = THIRD*((*PPC)[p[0]][0] + (*PPC)[p[1]][0] + (*PPC)[p[2]][0]);
		c0[1] = THIRD*((*PPC)[p[0]][1] + (*PPC)[p[1]][1] + (*PPC)[p[2]][1]);
		c0[2] = THIRD*((*PPC)[p[0]][2] + (*PPC)[p[1]][2] + (*PPC)[p[2]][2]);
	}
}

/*****************************************************************************/

void force_tri_3D::set(point_cloud* pPC, int ip0, int ip1, int ip2)
{
	assert(pPC != NULL);
	p[0] = ip0;
	p[1] = ip1;
	p[2] = ip2;
	PPC = pPC;
	calculate_centroid();
}

/*****************************************************************************/

bool force_tri_3D::operator==(const force_tri_3D& rhs) const
{
	// if one point and the plane normal match then assume triangles are
	// equal.  NB - this only applies when triangles are distributed over
	// the sphere.
	bool share_pt = false;
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			share_pt = share_pt || (operator[](i) == rhs[j]);
	bool pn_equal = (plane_normal() == rhs.plane_normal()) ||
                    (plane_normal() == rhs.plane_normal() * -1.0);
	return share_pt && pn_equal;
}

/*****************************************************************************/

bool force_tri_3D::operator!=(const force_tri_3D& rhs) const
{
	return !(*this == rhs);
}

/*****************************************************************************/

FP_TYPE force_tri_3D::surface_area(FP_TYPE radius) const
{
	// get the points
	vector_3D p0 = (*PPC)[p[0]];
	vector_3D p1 = (*PPC)[p[1]];
	vector_3D p2 = (*PPC)[p[2]];

	// multiply by the radius	
	p0 = p0 * radius;
	p1 = p1 * radius;
	p2 = p2 * radius;

	// calculate the area
	vector_3D P = p1 - p0;
	vector_3D Q = p2 - p0;
	vector_3D C = P.xp(Q);
	return 0.5 * C.mag();
}

/*****************************************************************************/

vector_3D force_tri_3D::plane_normal(void) const
{
	// calculate the normal to the plane defined by the 3 pts of the triangle
	vector_3D n = ((*PPC)[p[1]] - (*PPC)[p[0]]).xp((*PPC)[p[2]] - (*PPC)[p[0]]);
	return n / n.mag();
}
