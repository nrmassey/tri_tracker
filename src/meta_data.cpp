/******************************************************************************
** Program : meta_data.cpp
** Author  : Neil Massey
** Date    : 11/11/13
** Purpose : functions to add metadata to the pure binary file output
******************************************************************************/

#include "meta_data.h"
#include "bin_file_utils.h"

/*****************************************************************************/

void write_meta_data(std::ofstream& out, META_DATA_TYPE& meta_data)
{
	// write the metadata header
	out << 'M'; out << 'E'; out << 'T'; out <<'A';	// takes 4 bytes
	// write out the number of key, value pairs
	write_int(out, meta_data.size());
	// write out the pairs of strings as the metadata
	for (META_DATA_TYPE::iterator it_md = meta_data.begin(); 
		 it_md != meta_data.end(); it_md++)
	{
		std::string key = it_md->first;
		std::string value = it_md->second;
		write_string(out, key);
		write_string(out, value);
	}
}

/*****************************************************************************/

META_DATA_TYPE read_meta_data(std::ifstream& in)
{
	META_DATA_TYPE md_in;
	// read the first four bytes
	char m,e,t,a;
	in >> m >> e >> t >> a;
	if (m == 'M' && e == 'E' && t == 'T' && a == 'A')
	{
		// load the meta data
		int n_meta_items = read_int(in);
		for (int i=0; i<n_meta_items; i++)
		{
			std::string value;
			std::string key;
			key = read_string(in);
			value = read_string(in);
			md_in[key] = value;
		}
	}
	else
	{
		// otherwise rewind four bytes
		int c_pos = in.tellg();
		in.seekg(c_pos-4);
	}
	return md_in;
}