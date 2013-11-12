/******************************************************************************
** Program : steering_vector.h
** Author  : Neil Massey
** Date    : 05/08/13
** Purpose : Pure virtual class to be overloaded to calculate the steering
**           vectors to be used in the 3rd term of the cost function
******************************************************************************/

#ifndef STEERING_VECTOR_H
#define STEERING_VECTOR_H

#include <string>
#include "tri_grid.h"
#include "extremum.h"
#include "meta_data.h"

class steering_vector
{
	public:
		// sets the arguments to the steering vector calculator - these are
		// passed in from the command line
		virtual void parse_arg_string(std::string arg_string)=0;
		// calculate the steering vector from gridded data
		// tg = tri_grid, svex = steering extremum, t = timestep,
		// sv_u = steering vector u component, sv_v = steering vector v cmpt.
		virtual void calculate_steering_vector(tri_grid* tg, 
											   steering_extremum* svex, int t,
											   FP_TYPE mv=2e20)=0;
		// meta data operators
		META_DATA_TYPE* get_meta_data(void) {return &meta_data;}
		
	protected:
		META_DATA_TYPE meta_data;
};

#endif