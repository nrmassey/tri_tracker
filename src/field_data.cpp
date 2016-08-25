/******************************************************************************
** Program : field_data.cpp
** Author  : Neil Massey
** Date    : 19/06/16
** Purpose : class to represent a 2D field of data and carry out common 
**           functions, such as multiplication, min and max
******************************************************************************/

#include "field_data.h"
#include <stddef.h>
#include <cstring>
#include <algorithm>

/*****************************************************************************/

field_data::field_data(void) : data(NULL), x_len(0), y_len(0)
{
}

/*****************************************************************************/

field_data::field_data(int ix_len, int iy_len, const FP_TYPE v) 
           : x_len(ix_len), y_len(iy_len)
{
    data = new FP_TYPE[x_len*y_len];
    // set to value
    std::fill(data, data+x_len*y_len, v);
}

/*****************************************************************************/

field_data::~field_data(void)
{
    delete [] data;
}

/*****************************************************************************/

void field_data::set_size(int ix_len, int iy_len, const FP_TYPE v)
{
    x_len = ix_len;
    y_len = iy_len;
    if (data != NULL)
        delete [] data;
    data = new FP_TYPE[x_len*y_len];
    std::fill(data, data+x_len*y_len, v);
}

/*****************************************************************************/

void field_data::max_ip(field_data& rhs)
{
    // set the field values to the maximum of the two fields
    for (int i=0; i<x_len*y_len; i++)
        if (rhs.data[i] > data[i])
            data[i] = rhs.data[i];
}

/*****************************************************************************/

void field_data::min_ip(field_data& rhs)
{
    // set the field values to the minimum of the two fields
    for (int i=0; i<x_len*y_len; i++)
        if (rhs.data[i] < data[i])
            data[i] = rhs.data[i];
}

/*****************************************************************************/
        
void field_data::mult_ip(field_data& rhs)
{
    // set the field values to the multiplication of the two fields
    for (int i=0; i<x_len*y_len; i++)
        data[i] = data[i] * rhs.data[i];
}

/*****************************************************************************/

void field_data::div_ip(field_data& rhs)
{
    // set the field value to the field divided by the right hand field
    for (int i=0; i<x_len*y_len; i++)
        data[i] = data[i] / rhs.data[i];
}

/*****************************************************************************/

void field_data::add_ip(field_data& rhs)
{
    // in place addition of the two data fields
    for (int i=0; i<x_len*y_len; i++)
        data[i] = data[i] + rhs.data[i];
}

/*****************************************************************************/

void field_data::sub_ip(field_data& rhs)
{
    for (int i=0; i<x_len*y_len; i++)
        data[i] = data[i] - rhs.data[i];
}

/*****************************************************************************/

FP_TYPE field_data::get_min(void)
{
    FP_TYPE c_min = 2e20;
    for (int i=0; i<x_len*y_len; i++)
        if (data[i] < c_min)
            c_min = data[i];
    return c_min;
}

/*****************************************************************************/

FP_TYPE field_data::get_max(void)
{
    FP_TYPE c_max = -2e20;
    for (int i=0; i<x_len*y_len; i++)
        if (data[i] > c_max)
            c_max = data[i];
    return c_max;
}

/*****************************************************************************/

FP_TYPE field_data::get(int x, int y)
{
    return data[y*x_len + x];
}

/*****************************************************************************/

FP_TYPE* field_data::get(void)
{
    return data;
}

/*****************************************************************************/

void field_data::set(FP_TYPE* rhs_data)
{
    // copy the data
    memcpy(data, rhs_data, x_len*y_len*sizeof(FP_TYPE));
}
