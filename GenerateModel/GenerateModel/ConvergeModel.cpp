// ConvergeModel.cpp : 定义控制台应用程序的入口点。
//

// ConvertModel.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <mat.h>
#include "BaseInfo.h"
#include "AdjcentVertices.h"
#include "Boundary.h"
#include <Eigen/SVD>
#include <fstream>

//#include<Eigen/Dense>
#include <Eigen/SparseCore>
#include <string>
#include <stdio.h>
#include <string.h>
#include <windows.h>
#include "HumanModel.h"
#include <vector>
#include <io.h>
#include "Spring.h"
using namespace Eigen;
using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{
	DWORD timeStart = GetTickCount();



#if 0
	Spring Man("MatlabFile\\Male_SCAPE_Parameters.mat");
	double male_parameters[9] = { 0.55, 0.98, 0.58, 1.8, 0.97, 0.98, 0.82, 0.96, 1.03 };
	MatrixXd Tman(4, 4);
	Tman << 0.7978, 0.6025, -0.0245, -0.0910,
		0.0110, 0.0261, 0.9996, 0.4003,
		0.6029, -0.7977, 0.0142, -0.0092,
		0, 0, 0, 1;
	//MatrixXd mantheta(15,3);
	//mantheta.setZero();
	Man.GenerateModel(male_parameters, 9, Tman);
	Man.SaveModel2Obj("MatlabFile\\ManTest.obj");
	MatrixXd ManBodyVertex = Man.GetBodyVertices();
	MatrixXi ManBodyFaces = Man.GetBodyFaces();
	MatrixXd ManBodyColor = Man.GetBodyColor();

	HumanModel::HumanBody manbody;
	manbody.SetVertices(ManBodyVertex);
	manbody.SetFaces(ManBodyFaces);
	manbody.SetColors(ManBodyColor);
#endif



	Spring Woman("MatlabFile\\Female_SCAPE_Parameters.mat");
	double female_parameters[9] = { 0.52,0.81,0.48,1.55,0.81,0.67,0.64,0.67,0.86 };
	MatrixXd Twoman(4, 4);
	Twoman << 1, -0.09, -0.01, -0.04,
		0.01, 0.01, 1, 0.55,
		-0.09, -1, 0.01, 0.03,
		0, 0, 0, 1;
	//MatrixXd womantheta(15, 3);
	//womantheta.setZero();
	//womantheta.row(12) << 30, 30, 0;
	//womantheta.row(13) << -30, -30, 0;
	//cout << womantheta << endl;
	//Twoman.setIdentity(Twoman.rows(), Twoman.cols());
	//cout << Twoman << endl;
	Woman.GenerateModel(female_parameters, 9, Twoman);
	MatrixXd v = Woman.GetBodyVertices();
	const MatrixXi& Color = Woman.GetBodyFaces();
	Woman.SaveModel2Obj("MatlabFile\\WomanTest.obj",v,Color);

	MatrixXd WomanBodyVertex = Woman.GetBodyVertices();
	MatrixXi WomanBodyFaces = Woman.GetBodyFaces();
	MatrixXd WomanBodyColor = Woman.GetBodyColor();

	HumanModel::HumanBody womanbody;
	womanbody.SetVertices(WomanBodyVertex);
	womanbody.SetFaces(WomanBodyFaces);
	womanbody.SetColors(WomanBodyColor);

#if 1
	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	string HeadFilePath = "G:\\TestDataset\\Head1\\MatFile\\HeadWoman";
	string ResultFilePath("G:\\TestDataset\\ConvergeResult1\\ConvergeModelZhy");
	char Index[5];
	string s;
	string headfilename;
	string resultfilename;

	const int CountN = 1;
	for (int i = 1; i < CountN + 1; ++i)
	{
		//HumanModel::HumanBody body;
		HumanModel a;
		//a.ReadMatlabFile("G:\\TestDataset\\Body1\\MatFile\\Body.mat", body);
		Boundary NowBounday("MatlabFile\\Boundary.mat");

		HumanModel::HumanHead head;
		sprintf_s(Index, "%d", i);
		s = Index;
		headfilename = HeadFilePath + s + ".mat";
		resultfilename = ResultFilePath + s + ".obj";
		cout << headfilename << endl;
		const char* fhead = headfilename.c_str();

		a.ReadMatlabFile(fhead, head);


		a.ConvergeHeadAndBody(head, womanbody, NowBounday);

		const char* fresult = resultfilename.c_str();

		a.SaveModel2Obj(fresult, const_cast<MatrixXd&>(a.GetVertices()), const_cast<MatrixXi&>(a.GetFaces()));

	}
#endif
	QueryPerformanceCounter(&t2);
	cout << (t2.QuadPart - t1.QuadPart)*1.0 / tc.QuadPart << endl;
	DWORD timeEnd = GetTickCount();
	cout << "t : " << timeEnd - timeStart << endl;
	return 0;
}
