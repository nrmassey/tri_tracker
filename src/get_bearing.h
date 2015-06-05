/******************************************************************************
** Program : get_bearing.h
** Author  : Neil Massey
** Date    : 23/09/09
** Purpose : function to return bearing between two or three points
******************************************************************************/

#ifndef GET_BEARING_H
#define GET_BEARING_H

#include "vector_3D.h"

FP_TYPE get_bearing(FP_TYPE lon1, FP_TYPE lat1, FP_TYPE lon2, FP_TYPE lat2);
int get_sector(FP_TYPE b);
FP_TYPE get_curvature(FP_TYPE lon0, FP_TYPE lat0, FP_TYPE lon1, FP_TYPE lat1, FP_TYPE lon2, FP_TYPE lat2);

#endif
