// File    : Rot2Global.cpp
// Author  : Neil Massey, Tolu Aina, Simon Wilson
// Date    : 11/12/08
// Purpose : Functions to convert rotated grids to global grids

#include <math.h>

void Rot2Global(FP_TYPE  f_rot_lat,  FP_TYPE f_rot_lon,
                FP_TYPE  f_pole_lat, FP_TYPE f_pole_lon,
                FP_TYPE& f_global_lat, FP_TYPE& f_global_lon)
{
	// convert a rotated grid 
	FP_TYPE f_pi = M_PI;
    FP_TYPE f_deg_to_rad = f_pi/180.0;

    // Make sure rotlon is between 0 and 360
    while (f_rot_lon >= 360.0) f_rot_lon -= 360.0;
    while (f_rot_lon <    0.0) f_rot_lon += 360.0;
    // Make sure pole_lon is between 0 and 360
    while (f_pole_lon >= 360.0) f_pole_lon -= 360.0;
    while (f_pole_lon <    0.0) f_pole_lon += 360.0;

    // Convert inputs to radians
    f_rot_lon *= f_deg_to_rad;
    f_rot_lat *= f_deg_to_rad;
    f_pole_lon *= f_deg_to_rad;
    f_pole_lat *= f_deg_to_rad;

    // Amount rotated about 180E meridian
    FP_TYPE f_sock;
    if (f_pole_lon == 0.0)
        f_sock = 0.0;
    else
        f_sock = f_pole_lon - f_pi;

    FP_TYPE f_cpart = cos(f_rot_lon) * cos(f_rot_lat);
    FP_TYPE f_x = cos(f_pole_lat) * f_cpart + sin(f_pole_lat) * sin(f_rot_lat);

    if (f_x >=  1.0) f_x =  1.0;
    if (f_x <= -1.0) f_x = -1.0;
    f_global_lat = asin(f_x);

    FP_TYPE f_t1 = -1.0 * cos(f_pole_lat) * sin(f_rot_lat);
    FP_TYPE f_t2 = sin(f_pole_lat) * f_cpart;

    f_x = (f_t1 + f_t2) / cos(f_global_lat);
    if (f_x >=  1.0) f_x =  1.0;
    if (f_x <= -1.0) f_x = -1.0;
    f_global_lon = -1.0 * acos(f_x);

    // Make sure rotlon is between -PI and PI    
    while (f_rot_lon < -1*f_pi) f_rot_lon += 2.0*f_pi;
    while (f_rot_lon >    f_pi) f_rot_lon -= 2.0*f_pi;

    if (f_rot_lon >= 0.0 && f_rot_lon <= f_pi) f_global_lon *= -1.0 ;
    f_global_lon += f_sock;

    // Convert back to degrees
    f_global_lon /= f_deg_to_rad;
    f_global_lat /= f_deg_to_rad;
}

void Global2Rot(FP_TYPE  f_global_lat, FP_TYPE f_global_lon,
                FP_TYPE  f_pole_lat, FP_TYPE f_pole_lon,
                FP_TYPE& f_rot_lat, FP_TYPE& f_rot_lon)
{
	// convert a rotated grid 
	FP_TYPE f_pi = M_PI;
    FP_TYPE f_deg_to_rad = f_pi/180.0;

	// Make sure global lon is between 0 and 360
	while (f_global_lon >= 360.0) f_global_lon -= 360.0;
	while (f_global_lon < 0.0) f_global_lon += 360.0;

    // Make sure pole_lon is between 0 and 360
	while (f_pole_lon >= 360.0) f_pole_lon -= 360.0;
	while (f_pole_lon < 0.0) f_pole_lon += 360.0;

	// Convert inputs to radians
    f_global_lon *= f_deg_to_rad;
    f_global_lat *= f_deg_to_rad;
    f_pole_lon *= f_deg_to_rad;
    f_pole_lat *= f_deg_to_rad;

	// Amount rotated about 180E meridian
	FP_TYPE f_sock;
	if (f_pole_lon == 0.0)
		f_sock = 0.0;
	else
		f_sock = f_pole_lon - f_pi;

	// Need to get the screw in range -pi to pi
	FP_TYPE f_screw = f_global_lon - f_sock;
	while (f_screw < -1.0 * f_pi) f_screw += 2.0 * f_pi;
	while (f_screw > f_pi) f_screw -= 2.0 * f_pi;
	FP_TYPE f_bpart = cos(f_screw) * cos(f_global_lat);

	FP_TYPE f_x = (-1.0 * cos(f_pole_lat) * f_bpart) + (sin(f_pole_lat) * sin(f_global_lat));
	if (f_x >=  1.0) f_x =  1.0;
	if (f_x <= -1.0) f_x = -1.0;
	FP_TYPE f_lat2 = asin(f_x);

	FP_TYPE t1 = cos(f_pole_lat) * sin(f_global_lat);
	FP_TYPE t2 = sin(f_pole_lat) * f_bpart;

	f_x = (t1 + t2) / cos(f_lat2);

	if (f_x >=  1.0) f_x =  1.0;
	if (f_x <= -1.0) f_x = -1.0;
	FP_TYPE f_lon2 = -1.0 * acos(f_x);

	if (f_screw >= 0.0 and f_screw <= f_pi)
		f_lon2 *= -1.0;

	// Convert back to degrees
	f_rot_lat = f_lat2 / f_deg_to_rad;
	f_rot_lon = f_lon2 / f_deg_to_rad;
}
