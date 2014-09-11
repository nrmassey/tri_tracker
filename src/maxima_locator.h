/******************************************************************************
** Program : maxima_locator.h
** Author  : Neil Massey
** Date    : 10/07/13
** Purpose : class inherited from extrema_locator that searches for maxima 
**           in data regridded onto a regional hierarchical triangular mesh
******************************************************************************/

#ifndef MAXIMA_LOCATOR_H
#define MAXIMA_LOCATOR_H

#include "extrema_locator.h"

/*****************************************************************************/

class maxima_locator : public extrema_locator
{
	public:
		maxima_locator(void);
		~maxima_locator(void);
		virtual void parse_arg_string(std::string method_string);
		
	protected:
		// Virtual functions that require overloading
		virtual void calculate_object_position(int o, int t);
		virtual void calculate_object_intensity(int o, int t);
		virtual void calculate_object_delta(int o, int t);
		virtual bool is_extrema(indexed_force_tri_3D* tri, int t_step);
		virtual bool is_in_object(indexed_force_tri_3D* O_TRI, 
								  indexed_force_tri_3D* C_TRI, int t_step);
		FP_TYPE thresh, max_delta;
		int min_sur_tri;
};

#endif