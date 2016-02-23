#ifndef HAVERSINE_H
#define HAVERSINE_H

#include "vector_3D.h"

const FP_TYPE EARTH_R = 6371 * 1000;

FP_TYPE haversine(FP_TYPE lon1, FP_TYPE lat1, FP_TYPE lon2, FP_TYPE lat2, FP_TYPE R);
FP_TYPE haversine(vector_3D& point_1, vector_3D& point_2, FP_TYPE R);

#endif
