/******************************************************************************
** Program : spline.h
** Author  : Neil Massey
** Date    : 24/09/09
** Purpose : Spline class with control points that vary in t
******************************************************************************/

#ifndef SPLINE_H
#define SPLINE_H

#include <vector>

class spline
{
	public:
		// constructor, list of values and points
		spline(std::vector<FP_TYPE> v, std::vector<FP_TYPE> x, FP_TYPE mv);
		// evaluate the spline for pt t
		FP_TYPE evaluate(FP_TYPE t);

	private:
		void calc_coefficients(std::vector<FP_TYPE> v);
		std::vector<FP_TYPE> a;
		std::vector<FP_TYPE> b;
		std::vector<FP_TYPE> c;
		std::vector<FP_TYPE> d;
		std::vector<FP_TYPE> x;
		int find_piece(FP_TYPE t);
		FP_TYPE mv;
};

#endif
