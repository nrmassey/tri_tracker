/******************************************************************************
** Program : minima_locator.h
** Author  : Neil Massey
** Date    : 10/07/13
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded onto a regional hierarchical triangular mesh
******************************************************************************/

#ifndef MINIMA_LOCATOR_H
#define MINIMA_LOCATOR_H

#include "minima_processed.h"

/*****************************************************************************/

class minima_locator : public minima_processed
{
    public:
        minima_locator(void);
        ~minima_locator(void);
        virtual void parse_arg_string(std::string method_string);
        
    protected:
    
        /*********************************************************************/

        virtual bool process_data(void);

        /*********************************************************************/

        // Virtual functions that require overloading
        FP_TYPE trans_val(FP_TYPE val, vector_3D centroid);
        FP_TYPE pole_bck, equ_bck;
};

#endif