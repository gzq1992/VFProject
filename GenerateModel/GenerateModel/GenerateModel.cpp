// GenerateModel.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <string>
#include <iostream>
#include "SearchBody.h"
#include "Spring.h"
#include <Eigen/Core>
#include "HumanModel.h"
#include "Boundary.h"
using namespace std;
using namespace Eigen;
//int _tmain(int argc, _TCHAR* argv[])

int main(int argc,char* argv[])
{
#if 1
	

	string UserID(argv[6]);
	string HeadFile_Path = "E:\\VirtualFitting\\generateModel\\VirtualFittingDataset\\Head_Database\\";
	
	
	string HeadFile = HeadFile_Path + UserID + ".mat";
	string OutputFile_Path = "E:\\VirtualFitting\\generateModel\\VirtualFittingDataset\\User_Database\\";
	string logfile = OutputFile_Path + UserID + ".txt";
	string BodyFile = SearchBody(argc, argv,logfile.c_str());
	string OutputFile = OutputFile_Path + UserID + ".obj";

	string BoundaryFile("E:\\VirtualFitting\\generateModel\\VirtualFittingDataset\\MatlabFile\\Boundary.mat");
	//HumanModel.cpp line458
	Boundary boundary(BoundaryFile.c_str());

	HumanModel User;
	HumanModel::HumanBody body;
	User.ReadMatlabFile(BodyFile.c_str(), body);

	HumanModel::HumanHead head;
	User.ReadMatlabFile(HeadFile.c_str(), head);
	User.ConvergeHeadAndBody(head, body, boundary);
	User.ComputeNormals();
	//HumanModel.cpp line522
	User.SaveModel2Obj(OutputFile.c_str());
#endif

	return 0;
}

