/******************************************************************************
** Program : field_data.h
** Author  : Neil Massey
** Date    : 19/06/16
** Purpose : class to represent a 2D field of data and carry out common 
**           functions, such as multiplication, min and max
******************************************************************************/

#ifndef FIELD_DATA_H
#define FIELD_DATA_H

class field_data
{
    public:
        field_data(int x_len, int y_len);
        ~field_data(void);
        
        // in place functions only supported so far
        void max_ip(field_data& rhs);
        void min_ip(field_data& rhs);
        
        void mult_ip(field_data& rhs);  // result is (this) = (this)*rhs
        void div_ip(field_data& rhs);   // result is (this) = (this)/rhs
        void add_ip(field_data& rhs);   // result is (this) = (this)+rhs
        void sub_ip(field_data& rhs);   // result is (this) = (this)-rhs
    
        FP_TYPE get(int x, int y);      // return the value at x,y
        void set(FP_TYPE* rhs_data);    // set the values for the field
    
    private:
        FP_TYPE* data;
};

#endif