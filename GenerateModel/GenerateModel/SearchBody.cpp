#include "stdafx.h"
#include "SearchBody.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<fstream>

using std::string;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::ofstream;

string SearchBody(int argc, char* argv[], const char* fname)
{
	string filename;
	if (argc != 7)
	{
		cerr << "Input parameters wrong" << endl;
		return filename;
	}
	string body_filepath = "E:\\VirtualFitting\\generateModel\\VirtualFittingDataset\\SCAPE_Body_Database\\";
	string s(argv[1]);
	if (s == "Male")
	{
		body_filepath = body_filepath + "Male\\Male_Body_Mat\\Male_Body_";
	}
	else
	{
		if (s == "Female")
		{
			body_filepath = body_filepath + "Female\\Female_Body_Mat\\Female_Body_";
		}
		else
		{
			cerr << "Gender is not clear" << endl;
			return filename;
		}
	}

	vector<string> Vstr = ParameterConvert(argv, s);
	filename = body_filepath + "Height" + Vstr[0] + "_ChestLine" + Vstr[1] + "_WaistLine" + Vstr[2] + "_HipLine" + Vstr[3] + ".mat";

	ofstream out(fname);
	out << Vstr[0] << ' ' << Vstr[1] << ' ' << Vstr[2] << ' ' << Vstr[3] << ' ' << endl;

	return filename;

}

vector<string> ParameterConvert(char* argv[], string Gender)
{
	float Height[5], Chest[3], Waist[3], Hip[3];
	if (Gender == "Male")
	{
		float a[5] = { 160, 165, 170, 175, 180 };
		memcpy(Height, a, 5 * sizeof(float));
		float b[3] = { 95, 105, 115 };
		memcpy(Chest, b, 3 * sizeof(float));
		float c[3] = { 85, 95, 105 };
		memcpy(Waist, c, 3 * sizeof(float));
		float d[3] = { 90, 97, 104 };
		memcpy(Hip, d, 3 * sizeof(float));
	}
	if (Gender == "Female")
	{
		float a[5] = { 150, 155, 160, 165, 170 };
		memcpy(Height, a, 5 * sizeof(float));
		float b[3] = { 84, 92, 100 };
		memcpy(Chest, b, 3 * sizeof(float));
		float c[3] = { 70, 80, 90 };
		memcpy(Waist, c, 3 * sizeof(float));
		float d[3] = { 92, 101, 110 };
		memcpy(Hip, d, 3 * sizeof(float));
	}


	vector<string> result;
	char buffer[8];

	float val = atof(argv[2]);
	memset(buffer, 0, 8 * sizeof(char));
	_itoa_s(static_cast<int>(FindCloset(val, Height, 5)), buffer, 8, 10);
	result.push_back(string(buffer));

	val = atof(argv[3]);
	memset(buffer, 0, 8 * sizeof(char));
	_itoa_s(static_cast<int>(FindCloset(val, Chest, 3)), buffer, 8, 10);
	result.push_back(string(buffer));

	val = atof(argv[4]);
	memset(buffer, 0, 8 * sizeof(char));
	_itoa_s(static_cast<int>(FindCloset(val, Waist, 3)), buffer, 8, 10);
	result.push_back(string(buffer));

	val = atof(argv[5]);
	memset(buffer, 0, 8 * sizeof(char));
	_itoa_s(static_cast<int>(FindCloset(val, Hip, 3)), buffer, 8, 10);
	result.push_back(string(buffer));

	return result;
};

float FindCloset(float x, float table[], int LEN)
{

	if (x <= table[0])
	{
		return table[0];
	}
	if (x >= table[LEN - 1])
	{
		return table[LEN - 1];
	}
	for (int i = 1; i < LEN; ++i)
	{
		if (x <= table[i])
		{
			if (x - table[i - 1] < table[i] - x)
				return table[i - 1];
			else
				return table[i];
		}
	}

}
