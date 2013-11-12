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
		virtual bool is_extrema(indexed_force_tri_3D* tri, int t_step);
		virtual bool is_in_object(indexed_force_tri_3D* O_TRI, 
								  indexed_force_tri_3D* C_TRI, int t_step);
		virtual bool process_data(void);
		virtual FP_TYPE calculate_point_weight(FP_TYPE V, FP_TYPE min_v, FP_TYPE max_v);
		virtual void locate(void);
		virtual indexed_force_tri_3D* get_original_triangle(int o, int t);
		
	protected:

		/*********************************************************************/

		void calculate_background_field(void);
		void contour_data(void);
		void smooth_data(void);
		void get_child_values(QT_TRI_NODE* c_node, int max_depth, int t_step,
							  std::vector<FP_TYPE>& child_values);
		void trim_objects(void);
		void refine_objects(void);

		/*********************************************************************/
		
		std::string bck_field_file;
		int bck_mesh_lvl, bck_avg_period;
		int min_mesh_lvl;
		FP_TYPE contour_value, min_delta;
		data_store* bck_field_ds;
};

#endif