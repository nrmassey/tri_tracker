/******************************************************************************
** Program : get_bearing.cpp
** Author  : Neil Massey
** Date    : 23/09/09
** Purpose : function to return bearing between two points
******************************************************************************/

#include "get_bearing.h"
#include <math.h>
#include <iostream>

/*****************************************************************************/

FP_TYPE get_bearing(FP_TYPE lon1, FP_TYPE lat1, FP_TYPE lon2, FP_TYPE lat2)
{
    // get the bearing between two lon / lat pairs
    const FP_TYPE deg_to_rad = M_PI / 180.0;
    const FP_TYPE rad_to_deg = 180.0 / M_PI;
   
    FP_TYPE nlat1 = lat1 * deg_to_rad;
    FP_TYPE nlat2 = lat2 * deg_to_rad;
    FP_TYPE dlon  = (lon2 - lon1) * deg_to_rad;
    
    FP_TYPE a = sin(dlon) * cos(nlat2);
    FP_TYPE b = cos(nlat1) * sin(nlat2) - sin(nlat1) * cos(nlat2) * cos(dlon);
    // get the arccos and return
    FP_TYPE ang = atan2(a, b);
    ang = ang * rad_to_deg;
    return ang;  // returned in degrees;
}

/******************************************************************************/

int get_sector(FP_TYPE b)
{
    // get the sector of the bearing.  bearing and sector are related as:
    //       -180/180
    //           |       
    //       2   |   1   
    //           |       
    // -90 ------+------ 90
    //           |       
    //       3   |   0   
    //           |       
    //           0
    
    int S = 0;
    
    if ((b >= 0.0) && (b < 90.0))
        S = 0;
    else if ((b >= 90) && (b <= 180))
        S = 1;
    else if ((b >= -180) && (b < -90))
        S = 2;
    else if ((b >= -90) && (b < 0))
        S = 3;
    return S;
}

/******************************************************************************/

FP_TYPE get_curvature(FP_TYPE lon0, FP_TYPE lat0, FP_TYPE lon1, FP_TYPE lat1, FP_TYPE lon2, FP_TYPE lat2)
{
    // Calculate the curvature as the change in bearing between (pt1 -> pt2)
    // and (pt2 -> pt3)
    FP_TYPE c = 0.0;

    // same point incurs zero cost, rather than undefined
    if (lon1 == lon2 && lat1 == lat2)
        c = 0.0;
    else
    {
        // calculate the difference in the bearings
        FP_TYPE a = get_bearing(lon0, lat0, lon1, lat1);
        FP_TYPE b = get_bearing(lon1, lat1, lon2, lat2);
        //  get the sector for each bearing
        int Sa = get_sector(a);
        int Sb = get_sector(b);
        if (Sa == 1 && Sb == 2)
            c = b - a + 360;
        else if (Sa == 2 && Sb == 1)
            c = b - a - 360;
        else
            c = b - a;
    }
    return fabs(c);
}