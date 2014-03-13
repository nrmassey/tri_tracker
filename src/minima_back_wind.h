/******************************************************************************
** Program : minima_back_wind.h
** Author  : Neil Massey
** Date    : 16/09/13
** Purpose : class inherited from minima_background that searches for minima 
**           in data regridded, after the removal of the background field
**           and then uses the wind gradient to expand objects to contain a
**           wind storm object
******************************************************************************/

#ifndef MINIMA_BACK_WIND_H
#define MINIMA_BACK_WIND_H

#include "minima_background.h"
#include "ncdata.h"
#include <list>

/*****************************************************************************/

class minima_back_wind : public minima_background
{
	public:
		minima_back_wind(void);
		virtual ~minima_back_wind(void);	
		// Virtual functions that require overloading
		virtual void parse_arg_string(std::string method_string);		
		virtual void locate(void);
		virtual bool process_data(void);
		virtual bool is_in_object(indexed_force_tri_3D* O_TRI, 
								  indexed_force_tri_3D* C_TRI, int t_step);
		
	protected:

		/*********************************************************************/

		std::string wind_file_name;
		std::string wind_spd_field_name;
		ncdata* wind_spd_field;
		int ptile_thresh;
		FP_TYPE wind_thresh_value;

		/*********************************************************************/
		
		bool wind_test(indexed_force_tri_3D* O_TRI, indexed_force_tri_3D* C_TRI, int t_step);
		void expand_objects(void);
		virtual void calculate_object_position(int o, int t);
		virtual void calculate_object_intensity(int o, int t);
				
		/*********************************************************************/		
};

#endif