/******************************************************************************
** Program : minima_largescale.h
** Author  : Neil Massey
** Date    : 07/04/15
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded, after the removal of the large scale flow.
**           The large scale flow is defined as the regridded data at a lower
**           mesh resolution
******************************************************************************/

#ifndef MINIMA_LARGESCALE_H
#define MINIMA_LARGESCALE_H

#include "extrema_locator.h"

/*****************************************************************************/

class minima_largescale : public extrema_locator
{
    public:
        minima_largescale(void);
        ~minima_largescale(void);   
        // Virtual functions that require overloading
        virtual void parse_arg_string(std::string method_string);       
        
    protected:
    
        // Virtual functions that require overloading
        virtual void calculate_object_position(int o, int t);
        virtual void calculate_object_intensity(int o, int t);
        virtual void calculate_object_delta(int o, int t);
        virtual bool is_extrema(indexed_force_tri_3D* tri, int t_step);
        virtual bool is_in_object(indexed_force_tri_3D* O_TRI, 
                                  indexed_force_tri_3D* C_TRI, int t_step);
                                  
        /*********************************************************************/
        
        void get_min_max_values_delta(FP_TYPE& min, FP_TYPE& max, int o, int t);
        FP_TYPE get_val(indexed_force_tri_3D* TRI, int t);

        /*********************************************************************/
        
        int ls_msh_lvl; // mesh level at which large scale flow is determined to be
        int n_up;       // number of levels to traverse up the mesh to the large scale
        FP_TYPE contour_value, min_delta;
};

#endif