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
