/******************************************************************************
** Program : data_store.cpp
** Author  : Neil Massey
** Date    : 03/07/13
** Purpose : class to store the data regridded onto the regional triangular
**           grid.  NB - each triangle in the tri_grid contains a direct
**           method of accessing this store via the get_target_index method
******************************************************************************/

#include "data_store.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include "bin_file_utils.h"

/*****************************************************************************/

data_store::data_store() : data(NULL), missing_value(2e20)
{
}

/*****************************************************************************/

data_store::data_store(const data_store& rhs)
{
    *this = rhs;
}

/*****************************************************************************/

data_store::~data_store()
{
    delete [] data;
}

/*****************************************************************************/

void data_store::operator=(const data_store& rhs)
{
    missing_value = rhs.missing_value;
    set_size(rhs.get_number_of_time_steps(), rhs.get_number_of_indices());
    int size = rhs.get_number_of_time_steps() * rhs.get_number_of_indices();
    memcpy(data, rhs.data, size * sizeof(FP_TYPE));
    meta_data = rhs.meta_data;
}

/*****************************************************************************/

void data_store::set_size(int in_t_steps, int in_idxs)
{
    n_t_steps = in_t_steps;
    n_idxs = in_idxs;
    if (data != NULL)
        delete [] data;
    data = new FP_TYPE[n_t_steps * n_idxs];
    memset(data, static_cast<int>(0.0f), n_t_steps * n_idxs * sizeof(FP_TYPE));
}

/*****************************************************************************/

void data_store::zero(void)
{
    memset(data, static_cast<int>(0.0f), n_t_steps * n_idxs * sizeof(FP_TYPE)); 
}

/*****************************************************************************/

void data_store::copy(const data_store& rhs)
{
    // copy the data and info from one datastore to the next
    // get the info
    int n_ts = rhs.get_number_of_time_steps();
    int n_idxs = rhs.get_number_of_indices();
    FP_TYPE mv = rhs.get_missing_value();
    // set the sizes / info
    set_size(n_ts, n_idxs);
    set_missing_value(mv);
    // copy the data
    for (int i=0; i<n_ts * n_idxs; i++)
        data[i] = rhs.data[i];
}

/*****************************************************************************/

void data_store::set_data(int t_step, int idx, FP_TYPE value)
{
    int pos = t_step * n_idxs + idx;
#ifdef DEBUG
    if (pos > n_t_steps * n_idxs)
        throw "Index out of range - incorrect grid used?";
#endif
    data[pos] = value;
}

/*****************************************************************************/

int data_store::get_number_of_time_steps(void) const
{
    return n_t_steps;
}

/*****************************************************************************/

int data_store::get_number_of_indices(void) const
{
    return n_idxs;
}

/*****************************************************************************/

FP_TYPE data_store::get_data(int t_step, int idx)
{
    int pos = t_step * n_idxs + idx;
    return data[pos];
}

/*****************************************************************************/

void data_store::set_missing_value(FP_TYPE mv)
{
    missing_value = mv;
    // set all the data store to mv
    memset(data, static_cast<int>(mv), n_t_steps * n_idxs * sizeof(FP_TYPE));   
}

/*****************************************************************************/

FP_TYPE data_store::get_missing_value(void) const
{
    return missing_value;
}

/*****************************************************************************/

void data_store::set_meta_data(META_DATA_TYPE* in_meta_data)
{
    for (META_DATA_TYPE::iterator it_md = in_meta_data->begin();
         it_md != in_meta_data->end(); it_md++)
        meta_data[it_md->first] = it_md->second;
}

/*****************************************************************************/

META_DATA_TYPE* data_store::get_meta_data(void)
{
    return &meta_data;
}

/*****************************************************************************/


void data_store::load(std::string f_name)
{
    std::cout << "# Loading regridded data" << std::endl;
    // delete data if it already exists
    if (data)
        delete [] data;
    std::ifstream in;
    // open the file in binary mode
    in.open(f_name.c_str(), std::ios::in | std::ios::binary);
    if (!in)
        throw(std::string("Loading regridded data.  File could not be opened ") + f_name);
    meta_data = read_meta_data(in);
    // read the number of time steps and indices
    n_t_steps = read_int(in);
    n_idxs = read_int(in);
    missing_value = read_float(in);
    // create the data storage
    data = new FP_TYPE[n_t_steps * n_idxs];
    // read in the data
    in.read(reinterpret_cast<char*>(data), sizeof(FP_TYPE)*n_t_steps*n_idxs);
    in.close();
}

/*****************************************************************************/

void data_store::save(std::string f_name)
{
    std::cout << "# Saving regridded data to file: " << f_name << std::endl;
    std::ofstream out;
    out.open(f_name.c_str(), std::ios::out | std::ios::binary);
    if (!out)
        throw(std::string("Saving regridded data.  File could not be opened or written to: " + f_name));
    if (meta_data.size() != 0)
        write_meta_data(out, meta_data);
    // write out the number of time steps and indices
    write_int(out, n_t_steps);
    write_int(out, n_idxs);
    write_float(out, missing_value);
    // write the big array!
    out.write(reinterpret_cast<char*>(data), sizeof(FP_TYPE)*n_t_steps*n_idxs);
    out.close();
}

/*****************************************************************************/

void data_store::save_text(std::string f_name)
{
    std::ofstream out;
    out.open(f_name.c_str(), std::ios::out | std::ios::binary);
    if (!out)
        throw(std::string("Saving regridded data as text file.  File could not be opened or written to: " + f_name));
    // write out the number of time steps
    out << n_t_steps << std::endl;
    // write out the number of indices
    out << n_idxs << std::endl;
    // now write out the array
    for (int t=0; t<n_t_steps; t++)
    {
        for (int d=0; d<n_idxs; d++)
            out << data[t * n_idxs + d] << " ";
        out << std::endl;
    }
}