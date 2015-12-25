// SearchBody.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <cstring>
#include <string>
#include <vector>
using namespace std;
void SearchBody(int argc, char* argv[]);
//int _tmain(int argc, _TCHAR* argv[])
int main(int argc, char* argv[])
{
	if (argc != 6)
	{
		cerr << "Input parameters wrong" << endl;
		return 0;
	}
	string body_filepath = "H:\\VirtualFittingDataset\\Body\\";
	string s;
	s = string(argv[1]);
	if (s == "Male")
	{
		body_filepath = body_filepath + "Male\\Male_Body_Mat\\";
	}
	else
	{
		if (s == "Female")
		{
			body_filepath = body_filepath + "Female\\Female_Body_Mat\\";
		}
		else
		{
			cerr << "Gender is not clear" << endl;
			return 0;
		}
	}
	vector<string> Vstr[4];
	for (int i = 2; i < argc;++i)
	{
	}
	


	return 0;
}

