/******************************************************************************
** Program : force_point.cpp
** Author  : Neil Massey
** Date    : 25/06/09
** Purpose : class to hold points used in force repelling sphere equalisation
**           routine
******************************************************************************/

#include "force_point.h"
#include <stdlib.h>

/*****************************************************************************/

force_point::force_point(void) : vector_3D(0.0, 0.0, 0.0), nf(0)
{
	fv.zero();
}


/*****************************************************************************/

force_point::force_point(const FP_TYPE a, const FP_TYPE b, const FP_TYPE c)
			: vector_3D(a, b, c), nf(0)
{
	fv.zero();
}

/*****************************************************************************/

force_point::force_point(const vector_3D& v)
			: vector_3D(v), nf(0)
{
	fv.zero();
}

/*****************************************************************************/

force_point::force_point(const force_point& rhs) 
			: vector_3D(rhs), nf(rhs.nf), fv(rhs.fv)
{
}

/*****************************************************************************/

void force_point::reset_forces(void)
{
	fv.zero();
	nf = 0;
}

/*****************************************************************************/

void force_point::normalise(void)
{
	FP_TYPE l = 1.0/mag();
	x *= l;
	y *= l;
	z *= l;
}

/*****************************************************************************/

void force_point::add_force(const vector_3D& F)
{
	fv += F;
	nf += 1;
}

/*****************************************************************************/

void force_point::apply_forces(void)
{
	x += fv[0] / nf;
	y += fv[1] / nf;
	z += fv[2] / nf;
	normalise();
	reset_forces();
}

/*****************************************************************************/
