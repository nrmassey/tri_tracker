/******************************************************************************
** Program : minima_processed.h
** Author  : Neil Massey
** Date    : 10/04/15
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded and then processed in some way.  This is a
**           virtual class and needs to be inherited from.  It was written
**           as many feature identification routines have the same functions
**           after the data processing.  So all you have to do is write the
**           data processing function (and overload anything as necessary)
******************************************************************************/

#ifndef MINIMA_PROCESSED_H
#define MINIMA_PROCESSED_H

#include "extrema_locator.h"

/*****************************************************************************/

class minima_processed : public extrema_locator
{
    public:
        minima_processed(void);
        ~minima_processed(void);   
        // Virtual functions that require overloading
        virtual void parse_arg_string(std::string method_string)=0; 
        virtual void locate(void);

    protected:
    
        // Virtual functions that require overloading
        virtual void calculate_object_position(int o, int t);
        virtual void calculate_object_intensity(int o, int t);
        virtual void calculate_object_delta(int o, int t);
        virtual bool is_extrema(indexed_force_tri_3D* tri, int t_step);
        virtual bool is_in_object(indexed_force_tri_3D* O_TRI, 
                                  indexed_force_tri_3D* C_TRI, int t_step);
                                  
        /*********************************************************************/
        
        void get_min_max_values_processed(FP_TYPE& min, FP_TYPE& max, int o, int t);
        virtual bool process_data(void)=0; // pure virtual so have to inherit from the class.

        /*********************************************************************/
        
        bool is_mv(FP_TYPE V, FP_TYPE mv);
        FP_TYPE contour_data(FP_TYPE V, FP_TYPE C);
        FP_TYPE contour_value;
        FP_TYPE min_delta;

        /*********************************************************************/
        
        data_store* data_processed;
};

#endif