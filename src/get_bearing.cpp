/******************************************************************************
** Program : get_bearing.cpp
** Author  : Neil Massey
** Date    : 23/09/09
** Purpose : function to return bearing between two points
******************************************************************************/

#include "get_bearing.h"
#include <math.h>
#include <iostream>

FP_TYPE get_bearing(FP_TYPE lon1, FP_TYPE lat1, FP_TYPE lon2, FP_TYPE lat2)
{
    // get the bearing between two lon / lat pairs
    const FP_TYPE deg_to_rad = M_PI / 180.0;
    const FP_TYPE rad_to_deg = 180.0 / M_PI;
   
    FP_TYPE nlat1 = lat1 * deg_to_rad;
    FP_TYPE nlat2 = lat2 * deg_to_rad;
    FP_TYPE dlon  = (lon2 - lon1) * deg_to_rad;
    
    FP_TYPE a = sin(nlat1) * sin(nlat2) + cos(nlat1) * cos(nlat2) * cos(dlon);
    // get the arccos and return
    FP_TYPE ang = acos(a);
    ang = ang * rad_to_deg;
    return ang;  // returned in degrees;
}
