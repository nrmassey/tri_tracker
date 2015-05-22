#include "haversine.h"
#include <math.h>
#include <iostream>

FP_TYPE haversine(FP_TYPE lon1, FP_TYPE lat1, FP_TYPE lon2, FP_TYPE lat2, FP_TYPE R)
{
	// quick check to prevent very small distances at the poles
	if (lat1 == 90.0)
		lat1 = 90.5;
	if (lat1 == -90.0)
		lat1 = -90.5;
	if (lat2 == 90.0)
		lat2 = 90.5;
	if (lat2 == -90.0)
		lat2 = -90.5;
	// convert lon / lat into radians
	const FP_TYPE deg_to_rad = M_PI / 180.0;
	FP_TYPE lon1_r = lon1 * deg_to_rad;
	FP_TYPE lat1_r = lat1 * deg_to_rad;
	FP_TYPE lon2_r = lon2 * deg_to_rad;
	FP_TYPE lat2_r = lat2 * deg_to_rad;

	// differences
	FP_TYPE d_lon = lon2_r - lon1_r;
	FP_TYPE d_lat = lat2_r - lat1_r;

	// component parts of equation
	FP_TYPE a = sin(d_lat/2.0);
	FP_TYPE b = sin(d_lon/2.0);
	FP_TYPE c = a*a + cos(lat1_r) * cos(lat2_r) * b*b;

	// don't scupper the arcsin
	if (c > 1.0) c = 1.0;
	FP_TYPE d = R * 2 * asin(sqrt(c));
	return d;
}
