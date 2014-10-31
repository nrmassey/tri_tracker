/******************************************************************************
** Program : extrema_list.cpp
** Author  : Neil Massey
** Date    : 06/08/09
** Purpose : class to hold extrema data + io operators
******************************************************************************/

#include "extrema_list.h"
#include <fstream>
#include <iostream>
#include "vector_3D.h"
#include "geo_convert.h"
#include "bin_file_utils.h"

/*****************************************************************************/

extrema_list::extrema_list(void)
{
}

/*****************************************************************************/

extrema_list::extrema_list(int n_tsteps)
{
	ex_list.resize(n_tsteps);
}

/*****************************************************************************/

void extrema_list::set_size(int n_tsteps)
{
	ex_list.resize(n_tsteps);
}

/*****************************************************************************/

const int extrema_list::size(void) const
{
	return ex_list.size();
}

/*****************************************************************************/

const int extrema_list::number_of_extrema(int t) const
{
	return ex_list[t].size();
}

/*****************************************************************************/

steering_extremum* extrema_list::get(int t_step, int ex_n)
{
	return &(ex_list[t_step][ex_n]);
}

/*****************************************************************************/

const steering_extremum* extrema_list::get(int t_step, int ex_n) const
{
	return &(ex_list[t_step][ex_n]);
}

/*****************************************************************************/

void extrema_list::add(int t_step, steering_extremum ex)
{
	ex_list[t_step].push_back(ex);
}

/*****************************************************************************/

void extrema_list::remove(int t_step, int ex_n)
{
	ex_list[t_step][ex_n].object_labels.clear();		// set as removed
}

/*****************************************************************************/

void extrema_list::consolidate(FP_TYPE mv)
{
	// remove objects from the list that have been marked as having no
	// size - either the object label list is size 0 or the longitude
	// is the missing value
	bool stop = false;
	while (!stop)
	{
		stop = true;
		for (int t=0; t<size(); t++)
			for (int e=0; e<number_of_extrema(t); e++)
				if (ex_list[t][e].lon == mv || 
					ex_list[t][e].object_labels.size() == 0)
				{					
					ex_list[t].erase(ex_list[t].begin()+e);
					stop = false;
					break;
				}
	}
}

/*****************************************************************************/

META_DATA_TYPE* extrema_list::get_meta_data(void)
{
	return &meta_data;
}

/*****************************************************************************/

void extrema_list::set_meta_data(META_DATA_TYPE* in_meta_data)
{
	// copying like this rather than using the copy constructor allows multiple
	// sources of metadata
	for (META_DATA_TYPE::iterator it_md = in_meta_data->begin();
		 it_md != in_meta_data->end(); it_md++)
		meta_data[it_md->first] = it_md->second;
}

/*****************************************************************************/

void extrema_list::save(std::string output_fname, FP_TYPE mv)
{
	// save extrema to a binary format file
	std::ofstream out;
	std::cout << "# Saving extrema data" << std::endl;
	out.open(output_fname.c_str(), std::ios::out | std::ios::binary);	
	if (!out)
		throw(std::string("Saving extrema data.  File could not be opened or written to: " + output_fname));
	// write metadata
	if (meta_data.size() != 0)
		write_meta_data(out, meta_data);
	// write missing value
	write_float(out, mv);
	// write number of timesteps
	write_int(out, size());	
	for (int t=0; t<size(); t++)
	{
		// write number of extrema for this timestep
		write_int(out, number_of_extrema(t));
		for (int e=0; e<number_of_extrema(t); e++)
		{
			steering_extremum* svex = get(t, e);
			svex->save(out);
		}
	}
	out.close();
}

/*****************************************************************************/

void extrema_list::save_text(std::string output_fname, tri_grid* tg) 
{
	// save the extrema to a text file
	std::ofstream out;
	out.open(output_fname.c_str(), std::ios::out);
	if (!out)
		throw(std::string("Saving extrema data as text file.  File could not be opened or written to: " + output_fname));	
	for (int t=0; t<size(); t++)
	{
		for (int e=0; e<number_of_extrema(t); e++)
		{
			steering_extremum* svex = get(t, e);
			out << t << " ";
			svex->save_text(out);
			out << std::endl;
		}
	}
	out.close();
}

/*****************************************************************************/

void extrema_list::load(std::string input_fname, FP_TYPE& mv, bool append)
{
	// read extrema in a binary format file
	std::ifstream in_file;
	std::cout << "# Loading extrema data" << std::endl;
	in_file.open(input_fname.c_str(), std::ios::in | std::ios::binary);
	if (!in_file)
        throw(std::string("Loading extrema data.  File could not be opened: ") + input_fname);
    
    // read meta data
    meta_data = read_meta_data(in_file);
    // read the missing value in
    mv = read_float(in_file);
    // read the number of timesteps in
    int n_t_steps = read_int(in_file);
    // this offset allows for appending extrema to the end of the list
	int offset = 0;

	// set the size of the list
	if (append)
	{
		offset = ex_list.size();
		ex_list.resize(ex_list.size() + n_t_steps);
	}
	else
		ex_list.resize(n_t_steps);
	
	// now loop through the number of timesteps
	for (int t=0; t<n_t_steps; t++)
	{
		// read in the number of extrema for this time step
		int n_ex = read_int(in_file);
		// loop over the number of extrema
		for (int e=0; e<n_ex; e++)
		{
			steering_extremum svex;
			svex.load(in_file);
		 	add(t+offset, svex);
		}
	}
}
