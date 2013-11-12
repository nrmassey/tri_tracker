/******************************************************************************
** Program : geo_convert.cpp
** Author  : Neil Massey
** Date    : 08/05/09
** Purpose : convert model coords to cartesian
******************************************************************************/

#include <math.h>
#include "geo_convert.h"

vector_3D model_to_cart(FP_TYPE lon, FP_TYPE lat)
{
	// convert to radians
	const FP_TYPE R = 1.0;
	const FP_TYPE deg_to_rad = M_PI / 180.0;
	FP_TYPE lon_r = lon * deg_to_rad;
	// translate latitude first
	FP_TYPE lat_r = (90.0 - lat) * deg_to_rad;

	// convert spherical to cartesian on a unit sphere
	vector_3D S;
	S[0] = R * cos(lon_r) * sin(lat_r);
	S[1] = R * sin(lon_r) * sin(lat_r);
	S[2] = R * cos(lat_r);

	return S;
}

/*****************************************************************************/

void cart_to_model(const vector_3D& V, FP_TYPE& lon, FP_TYPE& lat)
{
	// convert a cartesian coordinate to a "model coordinate" - spherical
	// but with pole at 90 N
	const FP_TYPE rad_to_deg = 180.0 / M_PI;
	FP_TYPE phi = atan2(V[1], V[0]);
	FP_TYPE the = acos(V[2] / V.mag());

	lon = phi * rad_to_deg;
	if (lon < 0.0)	// adjust to 0..360 from -180..180
		lon += 360.0;
	lat = 90 - the * rad_to_deg;
}

/*****************************************************************************/

FP_TYPE grid_box_area(FP_TYPE lat, FP_TYPE lon_d, FP_TYPE lat_d)
{
	// calculate the area of a grid box of the above dimensions, on a unit
	// sphere - NB only longitude spacing is important
	const FP_TYPE deg_to_rad = M_PI / 180.0;
	FP_TYPE lat_r = lat * deg_to_rad;
	FP_TYPE lon_dr = lon_d * deg_to_rad;
	FP_TYPE lat_dr = 0.5 * lat_d * deg_to_rad;

	return lon_dr * fabs(sin(lat_r-lat_dr) - sin(lat_r+lat_dr));
}

/*****************************************************************************/

vector_3D model_to_norm(FP_TYPE lon, FP_TYPE lat)
{
	// return normalized co-ordinates (0..1) depicting lon / lat on a 2D
	// plane - origin is the South Pole
	vector_3D C;
	C[0] = lon / 360.0;
	C[1] = (lat + 90.0) / 180.0;
	C[2] = 0.0;
	return C;
}

/*****************************************************************************/

void norm_to_model(const vector_3D& V, FP_TYPE& lon, FP_TYPE& lat)
{
	lon = V[0] * 360.0;
	lat = V[1] * 180.0 - 90.0;
}

