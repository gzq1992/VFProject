#include "stdafx.h"
#include <vector>

std::string SearchBody(int argc, char* argv[],const char* fname);

std::vector<std::string> ParameterConvert(char* argv[], std::string Gender);

float FindCloset(float x, float table[], int LEN);