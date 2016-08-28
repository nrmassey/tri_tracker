/******************************************************************************
** Program : event.h
** Author  : Neil Massey
** Date    : 23/08/16
** Purpose : class that stores the windstorm event and its associated field
**           data
******************************************************************************/

#include "event.h"

/*****************************************************************************/

event::event(int x_len, int y_len) : mv(2e20), scale_mv(-128)
{
    // set the sizes of all the fields
    mslp.set_size(x_len, y_len, mv);
    wind_max.set_size(x_len, y_len, mv);
    wind_gust.set_size(x_len, y_len, mv);
    precip.set_size(x_len, y_len, mv);
    loss.set_size(x_len, y_len, mv);
    
    // calculate the scales and offsets
    // formulas are (reserving -128 as the missing value)
    // scale_factor = (max-min) / (2^n-2)   n=8 so 2^n-2 = 254
    // add_offset   = min-scalefactor
    
    const int N = 254;
    
    // mslp values between 910 to 1080 hPa
    const FP_TYPE mslp_min = 910.0e2;
    const FP_TYPE mslp_max = 1080.0e2;
    mslp_scale  = (mslp_max - mslp_min) / N;
    mslp_offset = (mslp_max + mslp_min) / 2;
    
    // wind (and gust) have values between 0 and 70 ms-1
    const FP_TYPE wind_max = 100.0;
    wind_scale  = (wind_max - 0.0) / N;
    wind_offset = 0.0;
    
    // precip flux has values between 0 and 0.005 kg m-2 s-1
    const FP_TYPE precip_max = 0.005;
    precip_scale  = (precip_max - 0.0) / N;
    precip_offset = 0;
    
    // not sure what loss is at the moment!
    const FP_TYPE loss_max = 1e12;
    loss_scale  = (loss_max - 0.0) / N;
    loss_offset = 0;
}

/*****************************************************************************/

void event::set_ref_data(ncdata* pref_data)
{
    ref_data = pref_data;
}