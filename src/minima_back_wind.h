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
#include <vector>

/*****************************************************************************/

class minima_back_wind : public minima_background
{
	public:
		minima_back_wind(void);
		~minima_back_wind(void);	
		// Virtual functions that require overloading
		virtual void parse_arg_string(std::string method_string);		
		virtual void locate(void);
		virtual bool process_data(void);		
		
	protected:

		/*********************************************************************/

		std::string wind_file_name;
		std::string wind_spd_field_name;
		ncdata* wind_spd_field;
		int previous_object;
		std::vector<std::vector<FP_TYPE> > all_wind_distr;

		/*********************************************************************/
		
		bool wind_test(indexed_force_tri_3D* O_TRI, indexed_force_tri_3D* C_TRI, int t_step);
		void expand_objects(void);
		
		/*********************************************************************/		
};

#endif