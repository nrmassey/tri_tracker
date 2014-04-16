/******************************************************************************
** Program : read_from_string.cpp
** Author  : Neil Massey
** Date    : 15/04/14
** Purpose : read from a string upto the delimiter and return the position
**           in the string of the next character
******************************************************************************/

#include <string>

int read_from_string(std::string str, int c_pos, const char* delim,
					 std::string& ret_str)
{
	// useful function for reading a filename, field_name, etc. from the string
	// takes current string position, returns next string position
	// delim is delimiting character to use
	// returns in ret_str
	int b_pos = str.find(delim, c_pos);
	ret_str = str.substr(c_pos, b_pos-c_pos);
	return b_pos+1;	
}