/******************************************************************************
** Program : minima_background.h
** Author  : Neil Massey
** Date    : 07/08/13
** Purpose : class inherited from extrema_locator that searches for minima 
**           in data regridded, after the removal of the background field
******************************************************************************/

#ifndef MINIMA_BACKGROUND_H
#define MINIMA_BACKGROUND_H

#include "extrema_locator.h"

/*****************************************************************************/

class minima_background : public extrema_locator
{
	public:
		minima_background(void);
		~minima_background(void);	
		// Virtual functions that require overloading
		virtual void parse_arg_string(std::string method_string);		
		virtual void locate(void);
		
	protected:
	
		// Virtual functions that require overloading
		virtual void calculate_object_position(int o, int t);
		virtual void calculate_object_intensity(int o, int t);
		virtual void calculate_object_delta(int o, int t);
		virtual bool is_extrema(indexed_force_tri_3D* tri, int t_step);
		virtual bool is_in_object(indexed_force_tri_3D* O_TRI, 
								  indexed_force_tri_3D* C_TRI, int t_step);
		virtual void find_extrema(void);
								  
		/*********************************************************************/

		void calculate_background_field(void);
		void trim_objects(void);
		void refine_objects(void);
		bool process_data(void);
		void get_min_max_values_delta(FP_TYPE& min, FP_TYPE& max, int o, int t);		

		/*********************************************************************/
		
		std::string bck_field_file;
		int bck_avg_period;
		int min_mesh_lvl;
		FP_TYPE contour_value, min_delta;
		data_store* bck_field_ds;
		data_store* data_minus_bck;
};

#endif