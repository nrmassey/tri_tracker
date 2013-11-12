/******************************************************************************
** Program : force_point.h
** Author  : Neil Massey
** Date    : 25/06/09
** Purpose : class to hold points used in force repelling sphere equalisation
**			 routine
******************************************************************************/

#ifndef FORCE_POINT_H
#define FORCE_POINT_H

#include "vector_3D.h"

class force_point : public vector_3D
{
	public:
		force_point(void);
		force_point(const vector_3D& v);
		force_point(const FP_TYPE a, const FP_TYPE b, const FP_TYPE c);
		force_point(const force_point& rhs);

		void reset_forces(void);
		void normalise(void);
		void add_force(const vector_3D& F);
		void apply_forces(void);

	private:
		int		  nf;	// number of forces applied this iteration
		vector_3D fv;	// current force vector
};

#endif
