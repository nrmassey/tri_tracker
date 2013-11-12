/******************************************************************************
** Program : maxima_locator.cpp
** Author  : Neil Massey
** Date    : 10/07/13
** Purpose : class inherited from extrema_locator that searches for maxima 
**           in data regridded onto a regional hierarchical triangular mesh
******************************************************************************/

#include "maxima_locator.h"
#include <math.h>

maxima_locator::maxima_locator(void)
			   :extrema_locator()
{
}

/*****************************************************************************/
				   
maxima_locator::~maxima_locator(void)
{
}

/*****************************************************************************/

void maxima_locator::parse_arg_string(std::string method)
{
}

/*****************************************************************************/

bool maxima_locator::is_extrema(indexed_force_tri_3D* tri, int t_step)
{
	return true;
}

/*****************************************************************************/

bool maxima_locator::is_in_object(indexed_force_tri_3D* O_TRI, 
							  	  indexed_force_tri_3D* C_TRI, int t_step)
{
	
	return true;
}

/*****************************************************************************/

bool maxima_locator::process_data(void)
{
	return true;
}

/*****************************************************************************/

FP_TYPE maxima_locator::calculate_point_weight(FP_TYPE V, FP_TYPE min_v, FP_TYPE max_v)
{
}

/*****************************************************************************/

indexed_force_tri_3D* maxima_locator::get_original_triangle(int o, int t)
{
	// get the original triangle - i.e. the one at the centre of the object
	// for the maxima this is defined as the triangle with the highest value
	FP_TYPE c_val = -2e20f;
	indexed_force_tri_3D* o_tri = NULL;
	LABEL_STORE* object_labels = &(ex_list.get(t, o)->object_labels);
	// if there are no labels do not try to find the min/max
	if (object_labels->size() == 0)
		return NULL;
	// get the missing value
	FP_TYPE mv = ds.get_missing_value();
	// loop over all the triangles in the object
	for (LABEL_STORE::iterator it_ll = object_labels->begin(); 
		 it_ll != object_labels->end(); it_ll++)	// tri indices
	{
		// get the triangle
		indexed_force_tri_3D* c_tri = tg.get_triangle(*it_ll);
		FP_TYPE val = ds.get_data(t, c_tri->get_ds_index());
		// find the min and max values
		if (fabs(val) < 0.99*fabs(mv) && val > c_val)
		{
			c_val = val;
			o_tri = c_tri;
		}
	}
	return o_tri;
}