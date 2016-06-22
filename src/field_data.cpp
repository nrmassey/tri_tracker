/******************************************************************************
** Program : field_data.cpp
** Author  : Neil Massey
** Date    : 19/06/16
** Purpose : class to represent a 2D field of data and carry out common 
**           functions, such as multiplication, min and max
******************************************************************************/

#include "field_data.h"

field_data::field_data(int x_len, int y_len)
{
    data = new FP_TYPE[x_len*y_len];
}

/*****************************************************************************/

field_data::~field_data(void)
{
    delete [] data;
}

/*****************************************************************************/
