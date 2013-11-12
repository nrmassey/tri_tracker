/******************************************************************************
** Program : geo_convert.h
** Author  : Neil Massey
** Date    : 08/07/09
** Purpose : convert model coordinates (lon, lat) to xyz cartesian coords
******************************************************************************/

#ifndef GEO_CONVERT_H
#define GEO_CONVERT_H

#include "vector_3D.h"

vector_3D model_to_cart(FP_TYPE lon, FP_TYPE lat);
FP_TYPE grid_box_area(FP_TYPE lat, FP_TYPE lon_d, FP_TYPE lat_d);
void cart_to_model(const vector_3D& V, FP_TYPE& lon, FP_TYPE& lat);
vector_3D model_to_norm(FP_TYPE lon, FP_TYPE lat);
void norm_to_model(const vector_3D& V, FP_TYPE& lon, FP_TYPE& lat);

#endif
