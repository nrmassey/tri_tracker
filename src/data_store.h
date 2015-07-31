/******************************************************************************
** Program : data_store.h
** Author  : Neil Massey
** Date    : 03/07/13
** Purpose : class to store the data regridded onto the regional triangular
**           grid.  NB - each triangle in the tri_grid contains a direct
**           method of accessing this store via the get_target_index method
******************************************************************************/

#ifndef DATA_STORE_H
#define DATA_STORE_H

#include <string>
#include "meta_data.h"

class data_store
{
	public:
		data_store(void);
		data_store(const data_store& rhs);
		~data_store(void);
		void operator=(const data_store& rhs);
		void set_size(int n_t_steps, int n_idxs);
		void set_data(int t_step, int idx, FP_TYPE value);
		void zero(void);		// set all to zero
		void copy(const data_store& rhs);
		int get_number_of_time_steps(void) const;
		int get_number_of_indices(void) const;
		FP_TYPE get_data(int t_step, int idx);
		void set_missing_value(FP_TYPE mv);
		FP_TYPE get_missing_value(void) const;
		void load(std::string f_name);
		void save(std::string f_name);
		void save_text(std::string f_name);
		void set_meta_data(META_DATA_TYPE* meta_data);
		META_DATA_TYPE* get_meta_data(void);
		
	private:
		FP_TYPE* data;
		FP_TYPE missing_value;
		int n_t_steps;
		int n_idxs;
		META_DATA_TYPE meta_data;
};

#endif