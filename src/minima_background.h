/******************************************************************************
** Program : minima_background.h
** Author  : Neil Massey
** Date    : 07/08/13
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded, after the removal of the background field
******************************************************************************/

#ifndef MINIMA_BACKGROUND_H
#define MINIMA_BACKGROUND_H

#include "minima_processed.h"

/*****************************************************************************/

class minima_background : public minima_processed
{
    public:
        minima_background(void);
        ~minima_background(void);   
        // Virtual functions that require overloading
        virtual void parse_arg_string(std::string method_string);       
        
    protected:
                                      
        /*********************************************************************/

        void calculate_background_field(void);
        virtual bool process_data(void);

        /*********************************************************************/
        
        std::string bck_field_file;
        int bck_avg_period;
        data_store* bck_field_ds;
};

#endif