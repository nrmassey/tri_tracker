//*****************************************************************************
//* File    : Rot2Global.h
//* Author  : Neil Massey, Tolu Aina, Simon Wilson
//* Date    : 11/12/08
//* Purpose : Functions to convert rotated grids to global grids
//*****************************************************************************

#ifndef ROT_2_GLOBAL
#define ROT_2_GLOBAL

void Rot2Global(FP_TYPE  f_rot_lat,  FP_TYPE f_rot_lon,
                FP_TYPE  f_pole_lat, FP_TYPE f_pole_lon,
                FP_TYPE& f_global_lat, FP_TYPE& f_global_lon);

#endif
