/******************************************************************************
** Program : meta_data.h
** Author  : Neil Massey
** Date    : 11/11/13
** Purpose : functions to add metadata to the pure binary file output
******************************************************************************/

#ifndef META_DATA_H
#define META_DATA_H

#include <map>
#include <string>
#include <fstream>

typedef std::map<std::string, std::string> META_DATA_TYPE;

void write_meta_data(std::ofstream& out, META_DATA_TYPE& meta_data);
META_DATA_TYPE read_meta_data(std::ifstream& in);

#endif