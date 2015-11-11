/******************************************************************************
** Program : mslp_wind_vector.h
** Author  : Neil Massey
** Date    : 11/11/15
** Purpose : Calculate the geostrophic wind from the MSLP field of a file passed
**           in.  Inherits from steering vector class
******************************************************************************/

#ifndef MSLP_WIND_VECTOR_H
#define MSLP_WIND_VECTOR_H

#include <string>
#include "steering_vector.h"
#include "ncdata.h"

class mslp_wind_vector : public steering_vector
{
    public:
        mslp_wind_vector(void);
        ~mslp_wind_vector(void);
        // sets the arguments to the steering vector calculator - these are
        // passed in from the command line
        void parse_arg_string(std::string arg_string);
        // calculate the steering vector from gridded data
        // t = timestep, i = x index in grid, j = y index in grid,
        // sv_u = steering vector u component, sv_v = steering vector v cmpt.
        void calculate_steering_vector(tri_grid* tg, 
                                       steering_extremum* svex, int t, 
                                       FP_TYPE mv=2e20);
    private:
        ncdata* mslp;
        int max_tri_level;
};

#endif