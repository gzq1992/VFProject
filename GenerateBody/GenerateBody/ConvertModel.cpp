// ConvertModel.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <mat.h>
#include "BaseInfo.h"
#include "AdjcentVertices.h"
#include "Boundary.h"
#include <Eigen/SVD>
#include <fstream>
#include<Eigen/Dense>
#include <string>
#include <stdio.h>
#include <string.h>
#include <windows.h>
#include "HumanModel.h"
#include <vector>
#include <io.h>
using namespace Eigen;
using namespace std;

#if 0
void getFiles(string path, vector<string>& files)
{
	//文件句柄  
	long   hFile = 0;
	//文件信息  
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			//如果是目录,迭代之  
			//如果不是,加入列表  
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}
#endif


int _tmain(int argc, _TCHAR* argv[])
{
#if 1
	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	string HeadFilePath = "G:\\TestDataset\\Head1\\MatFile\\MergedHead";
	string ResultFilePath("G:\\TestDataset\\ConvergeResult1\\ConvergeModelZhy");
	char Index[5];
	string s;
	string headfilename;
	string resultfilename;

	const int CountN = 1;
	for (int i = 1; i < CountN + 1; ++i)
	{
		HumanModel::HumanBody body;
		HumanModel a;
		a.ReadMatlabFile("G:\\TestDataset\\Body1\\MatFile\\Body.mat", body);

		Boundary NowBounday("MatlabFile\\Boundary.mat");

		HumanModel::HumanHead head;
		sprintf_s(Index, "%d", i);
		s = Index;
		headfilename = HeadFilePath + s + ".mat";
		resultfilename = ResultFilePath + s + ".obj";
		cout << headfilename << endl;
		const char* fhead = headfilename.c_str();

		a.ReadMatlabFile(fhead, head);


		a.ConvergeHeadAndBody(head, body,NowBounday);

		const char* fresult = resultfilename.c_str();

		a.SaveModel2Obj(fresult,const_cast<MatrixXd&>(a.GetVertices()),const_cast<MatrixXi&>(a.GetFaces()));
		
	}
#endif
	QueryPerformanceCounter(&t2);
	cout<< (t2.QuadPart - t1.QuadPart)*1.0 / tc.QuadPart<<endl;
	return 0;
}