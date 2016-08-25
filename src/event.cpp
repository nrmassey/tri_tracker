/******************************************************************************
** Program : event.h
** Author  : Neil Massey
** Date    : 23/08/16
** Purpose : class that stores the windstorm event and its associated field
**           data
******************************************************************************/

#include "event.h"

/*****************************************************************************/

event::event(int x_len, int y_len) : mv(2e20), scale_mv(0),
       mslp_offset(85000.0), mslp_scale(100.0),
       wind_offset(0.0), wind_scale(0.5),
       precip_offset(0.0), precip_scale(0.5),
       loss_offset(0.0), loss_scale(100.0)
{
    // set the sizes of all the fields
    mslp.set_size(x_len, y_len, mv);
    wind_max.set_size(x_len, y_len, mv);
    wind_gust.set_size(x_len, y_len, mv);
    precip.set_size(x_len, y_len, mv);
    loss.set_size(x_len, y_len, mv);
}

/*****************************************************************************/

void event::set_ref_data(ncdata* pref_data)
{
    ref_data = pref_data;
}