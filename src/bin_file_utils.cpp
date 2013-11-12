/******************************************************************************
** Program : bin_file_utils.cpp
** Author  : Neil Massey
** Date    : 26/06/13
** Purpose : functions to read and write values to a pure binary file
******************************************************************************/

#include "bin_file_utils.h"

/*****************************************************************************/


void write_string(std::ofstream& out, std::string string)
{
	// write out length of string and each character in the string
	write_int_as_byte(out, string.size());
	for (unsigned int c=0; c<string.size();c++)
		out << string.c_str()[c];
}

/*****************************************************************************/

std::string read_string(std::ifstream& in)
{
	std::string in_string = "";
	// read in length of the string
	int sl = read_int_as_byte(in);
	for (int sc=0; sc<sl; sc++)
	{
		char c;
		in >> c;
		in_string += c;
	}
	return in_string;
}

/*****************************************************************************/

#if FP_TYPE==float
void write_float(std::ofstream& out, float val)
{
	out.write(reinterpret_cast<char*>(&val), sizeof(val));
}

float read_float(std::ifstream& in)
{
	float val;
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	return val;
}

/*****************************************************************************/

#else
void write_float(std::ofstream& out, double val)
{
	out.write(reinterpret_cast<char*>(&val), sizeof(val));
}

double read_float(std::ifstream& in)
{
	double val;
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	return val;
}

#endif

/*****************************************************************************/

void write_int(std::ofstream& out, int val)
{
	out.write(reinterpret_cast<char*>(&val), sizeof(val));
}

/*****************************************************************************/

void write_int_as_byte(std::ofstream& out, int val)
{
	out.write(reinterpret_cast<char*>(&val), 1);
}

/*****************************************************************************/

int read_int(std::ifstream& in)
{
	int val;
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	return val;
}

/*****************************************************************************/

int read_int_as_byte(std::ifstream& in)
{
	char val;
	in.read(&val, 1);
	return static_cast<int>(val);
}

/*****************************************************************************/

LABEL read_label(std::ifstream& in)
{
	LABEL val;
	in.read(reinterpret_cast<char*>(&(val.label)), sizeof(val.label));
	val.max_level = read_int_as_byte(in);
	return val;
}

/*****************************************************************************/

void write_label(std::ofstream& out, LABEL val)
{
	// write out the label's index and max_level
	out.write(reinterpret_cast<char*>(&(val.label)), sizeof(val.label));
	write_int_as_byte(out, val.max_level);
}

/*****************************************************************************/

void write_vector(std::ofstream& out, vector_3D& vec)
{
	for (int v=0; v<3; v++)
		write_float(out, vec[v]);
}

/*****************************************************************************/

vector_3D read_vector(std::ifstream& in)
{
	double vi[3];
	for (int v=0; v<3; v++)
		vi[v] = read_float(in);
	vector_3D vec(vi[0], vi[1], vi[2]);
	return vec;
}

/*****************************************************************************/