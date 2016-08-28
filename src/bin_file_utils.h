/******************************************************************************
** Program : bin_file_utils.h
** Author  : Neil Massey
** Date    : 26/06/13
** Purpose : functions to read and write values to a pure binary file
******************************************************************************/

#ifndef BIN_FILE_UTILS_H
#define BIN_FILE_UTILS_H

#include <string>
#include <fstream>
#include "vector_3D.h"
#include "indexed_force_tri_3D.h"

void write_string(std::ofstream& out, std::string string);
std::string read_string(std::ifstream& in);

#if FP_TYPE==float
	void write_float(std::ofstream& out, float val);
	float read_float(std::ifstream& in);
#else
	void write_float(std::ofstream& out, double val);
	double read_float(std::ifstream& in);
#endif

void write_int(std::ofstream& out, int val);
int read_int(std::ifstream& in);

void write_label(std::ofstream& out, LABEL val);
LABEL read_label(std::ifstream& in);

void write_vector(std::ofstream& out, vector_3D& vec);
vector_3D read_vector(std::ifstream& in);

#endif