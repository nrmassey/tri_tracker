#ifndef READ_FROM_STRING_H
#define READ_FROM_STRING_H

#include <string>

int read_from_string(std::string str, int c_pos, const char* delim, std::string& ret_str);

#endif