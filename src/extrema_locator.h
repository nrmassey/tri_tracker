/******************************************************************************
** Program : extrema_locator.h
** Author  : Neil Massey
** Date    : 10/07/13
** Purpose : class that searches for extrema in data regridded onto a 
**           regional hierarchical triangular mesh - provided by regrid
**           this class is abstract and should be inherited from.  See
**           minima_locator.h and maxima_locator.h for an example
******************************************************************************/

#ifndef EXTREMA_LOCATOR_H
#define EXTREMA_LOCATOR_H

#include <list>
#include "indexed_force_tri_3D.h"
#include "tri_grid.h"
#include "data_store.h"
#include "extrema_list.h"
#include "ncdata.h"
#include "steering_vector.h"
#include "meta_data.h"

/*****************************************************************************/

class extrema_locator
{
	public:
		extrema_locator(void);
		virtual ~extrema_locator(void);
		virtual void locate(void);
		void save(std::string output_fname, bool save_text=false);
		void set_steering_vector(steering_vector* sv);
		void set_inputs(std::string input_fname, std::string mesh_fname,
					    int grid_level, ADJACENCY adj_type);
		void calculate_steering_vector(int o, int t);
		
		// Virtual functions that require overloading
		virtual void parse_arg_string(std::string method_string) = 0;
		virtual bool is_extrema(indexed_force_tri_3D* tri, int t_step) = 0;
		virtual bool is_in_object(indexed_force_tri_3D* O_TRI, 
						  		  indexed_force_tri_3D* C_TRI, int t) = 0;
		virtual bool process_data(void) = 0;
		virtual FP_TYPE calculate_point_weight(FP_TYPE V, FP_TYPE min_v, FP_TYPE max_v) = 0;
		virtual indexed_force_tri_3D* get_original_triangle(int o, int t) = 0;
		
	protected:	
		/*********************************************************************/
		// control variables
		int grid_level;
		ADJACENCY adj_type;
		
		/*********************************************************************/
		// mesh and data storage
		tri_grid tg;
		data_store ds;
		extrema_list ex_list;
		steering_vector* sv;		// steering vector class
		META_DATA_TYPE meta_data;	// add the meta data as we parse the methods and steering vector
		
		/*********************************************************************/

		void find_extrema(void);
		void find_objects(void);
		void merge_objects(void);
		void ex_points_from_objects(void);
		void get_min_max_values(FP_TYPE& min, FP_TYPE& max, int o, int t);
		virtual void calculate_object_position(int o, int t);
		virtual void calculate_object_intensity(int o, int t);
		virtual void calculate_object_delta(int o, int t);
		bool objects_share_nodes(const LABEL_STORE* o1,
								 const LABEL_STORE* o2);
		void tstep_out_begin(int t);
		void tstep_out_end(int t);
		
		int read_from_string(std::string str, int c_pos, const char* delim,
							 std::string& ret_str);
};

#endif