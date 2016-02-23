/******************************************************************************
** Program : extrema_list.h
** Author  : Neil Massey
** Date    : 06/08/09
** Purpose : class to hold extrema data + io operators
******************************************************************************/

#ifndef EXTREMA_LIST_H
#define EXTREMA_LIST_H

#include <vector>
#include <iostream>
#include <string>

#include "extremum.h"
#include "tri_grid.h"
#include "meta_data.h"

typedef std::vector<std::vector<steering_extremum> > EX_LIST_TYPE;

/*****************************************************************************/

class extrema_list
{
	public:
		extrema_list(void);
		extrema_list(int n_tsteps);
		void set_size(int n_tsteps);
		const int size(void) const;
		const int number_of_extrema(int t) const;
		steering_extremum* get(int t_step, int ex_n);
		const steering_extremum* get(int t_step, int ex_n) const;
		void add(int t_step, steering_extremum ex);
		void remove(int t_step, int ex_n);
		void save(std::string output_fname, FP_TYPE mv);
		void save_text(std::string output_fname, tri_grid* tg);
		void load(std::string input_fname, FP_TYPE& mv, bool append=false);
		META_DATA_TYPE* get_meta_data(void);
		void set_meta_data(META_DATA_TYPE* meta_data);

	private:
		EX_LIST_TYPE ex_list;
		META_DATA_TYPE meta_data;
};

#endif
