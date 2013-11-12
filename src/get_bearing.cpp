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
    FP_TYPE nlat1 = lat1 * deg_to_rad;
    FP_TYPE nlat2 = lat2 * deg_to_rad;
    FP_TYPE dlon  = (lon2 - lon1) * deg_to_rad;
    FP_TYPE a = sin(nlat1)*sin(nlat2) + cos(nlat1)*cos(nlat2)*cos(dlon);
    FP_TYPE ang = 0.0;
//    if (a >= -1.0 and a <= 1.0)
        ang = acos(a);
    return ang;
}

/*****************************************************************************/

FP_TYPE ang_between(const vector_3D& V1, const vector_3D& V2)
{
	FP_TYPE mag1 = V1.mag();
	FP_TYPE mag2 = V2.mag();
	FP_TYPE a = 0.0;
	if (!(mag1 == 0.0 || mag2 == 0.0))
//		a = acos(V2.dp(V1) / (mag1 * mag2));
	// switch to using atan2 as it has a sign associated with it
		a = atan2(V2[1],V2[0]) - atan2(V1[1],V1[0]);
	return a;
}

/*****************************************************************************/

FP_TYPE ang_between(const vector_3D& V1, const vector_3D& V2, const vector_3D& V3)
{
	// measure angle between two vectors V1->V2 & V2->V3
	vector_3D vd1 = V2 - V1;
	vector_3D vd2 = V3 - V2;
	FP_TYPE a = 0.0;
	FP_TYPE mag1 = vd1.mag();
	FP_TYPE mag2 = vd2.mag();
	if (!(mag1 == 0.0 || mag2 == 0.0))
//		a = acos(vd2.dp(vd1) / (mag1 * mag2));
	// switch to using atan2 as it has a sign associated with it
		a = atan2(vd2[1], vd2[0]) - atan2(vd1[1], vd1[0]);
	return a;
}
