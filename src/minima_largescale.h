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

#include "minima_processed.h"

/*****************************************************************************/

class minima_largescale : public minima_processed
{
    public:
        minima_largescale(void);
        ~minima_largescale(void);   
        // Virtual functions that require overloading
        virtual void parse_arg_string(std::string method_string);       

    protected:
            
        /*********************************************************************/
        
        bool process_data(void);
        void smooth_processed_data(void);

        /*********************************************************************/
        
        int ls_msh_lvl; // mesh level at which large scale flow is determined to be
};

#endif