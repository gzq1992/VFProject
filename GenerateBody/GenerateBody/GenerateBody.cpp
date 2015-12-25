// GenerateBody.cpp : 定义控制台应用程序的入口点。

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




	Spring Man("MatlabFile\\Male_SCAPE_Parameters.mat");
	MatrixXd Tman(4, 4);
	Tman << 0.7978, 0.6025, -0.0245, -0.0910,
		0.0110, 0.0261, 0.9996, 0.4003,
		0.6029, -0.7977, 0.0142, -0.0092,
		0, 0, 0, 1;
	double male_parameters[9] = { 0 };
	string matfilepath = "MatlabFile\\Male\\Male_Body_Mat\\";
	string objfilepath = "MatlabFile\\Male\\Male_Body_Obj\\";
	string fullfilepath = "MatlabFile\\Male\\Male_Full\\";
	string matfilename = matfilepath;
	string objfilename = objfilepath;
	string fullfilename;
	string s_filename;
	char Index[4];
	string s;
#if 0
	for (int height_index = 0; height_index < 5; ++height_index)
	{
		double height = 160 + 5 * height_index;
		//cout << height << endl;
		male_parameters[3] = height / 100;
		male_parameters[0] = 0.3246 * height / 100;
		male_parameters[1] = 0.5156 * height / 100;
		male_parameters[2] = 0.3082 * height / 100;
		//sprintf_s(Index, "%d", static_cast<int>(height * 100));
		//s = Index;
		//s_filename = "Male_Body_Height" + s + "_";
		for (int chestsize_index = 0; chestsize_index < 3;++chestsize_index)
		{
			double chestsize = 95 + chestsize_index * 10;
			male_parameters[4] = chestsize / 100;
			male_parameters[5] = 0.9741 * chestsize / 100;
			//sprintf_s(Index, "%d", static_cast<int>(chestsize * 100));
			//s = Index;
			//s_filename = s_filename + "ChestLine" + s + "_";
			for (int waistsize_index = 0; waistsize_index < 3;++waistsize_index)
			{
				double waistsize = 85 + waistsize_index * 10;

				male_parameters[6] = waistsize / 100;
				male_parameters[7] = 1.0297 * waistsize / 100;
				//sprintf_s(Index, "%d", static_cast<int>(waistsize * 100));
				//s = Index;
				//s_filename = s_filename + "WaistLine" + s + "_";
				for (int hipsize_index = 0; hipsize_index < 3;++hipsize_index)
				{
					double hipsize = 90 + hipsize_index * 7;

					male_parameters[8] = hipsize / 100;

					sprintf_s(Index, "%d", static_cast<int>(height));
					s = Index;
					s_filename = "Male_Body_Height" + s + "_";

					sprintf_s(Index, "%d", static_cast<int>(chestsize));
					s = Index;
					s_filename = s_filename + "ChestLine" + s + "_";

					sprintf_s(Index, "%d", static_cast<int>(waistsize));
					s = Index;
					s_filename = s_filename + "WaistLine" + s + "_";
					
					sprintf_s(Index, "%d", static_cast<int>(hipsize));
					s = Index;
					matfilename = matfilepath + s_filename + "HipLine" + s + ".mat";
					objfilename = objfilepath + s_filename + "HipLine" + s + ".obj";
					fullfilename = fullfilepath + s_filename + "HipLine" + s + ".obj";

					Man.GenerateModel(male_parameters, 9, Tman);
					MatrixXd ManBodyVertex = Man.GetBodyVertices();
					MatrixXi ManBodyFaces = Man.GetBodyFaces();

					Man.SaveModel2Obj(fullfilename.c_str());
					Man.SaveModel2Obj(objfilename.c_str(), ManBodyVertex, ManBodyFaces);
					Man.WriteMatFile(matfilename.c_str(), ManBodyVertex, ManBodyFaces);
				}
			}
		}
	}



	//double male_parameters[9] = { 0.55, 0.98, 0.58, 1.8, 0.97, 0.98, 0.82, 0.96, 1.03 };

	//MatrixXd mantheta(15,3);
	//mantheta.setZero();
#if 0
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

#endif

#if 1
	Spring Woman("MatlabFile\\Female_SCAPE_Parameters.mat");
	MatrixXd Twoman(4, 4);
	Twoman << 1, -0.09, -0.01, -0.04,
		0.01, 0.01, 1, 0.55,
		-0.09, -1, 0.01, 0.03,
		0, 0, 0, 1;
	double female_parameters[9] = {0};
	matfilepath = "MatlabFile\\Female\\Female_Body_Mat\\";
	objfilepath = "MatlabFile\\Female\\Female_Body_Obj\\";
	fullfilepath = "MatlabFile\\Female\\Female_Full\\";
	matfilename = matfilepath;
	objfilename = objfilepath;
	
	for (int height_index = 0; height_index < 5; ++height_index)
	{
		double height = 150 + 5 * height_index;
		female_parameters[3] = height / 100;
		female_parameters[0] = 0.3307 * height / 100;
		female_parameters[1] = 0.5251 * height / 100;
		female_parameters[2] = 0.2968 * height / 100;
		//sprintf_s(Index, "%d", static_cast<int>(height * 100));
		//s = Index;
		//s_filename = "Female_Body_Height" + s + "_";
		for (int chestsize_index = 0; chestsize_index < 3; ++chestsize_index)
		{
			//double chestsize = 0.84 + chestsize_index * 0.08;
			double chestsize = 84 + chestsize_index * 8;

			female_parameters[4] = chestsize / 100;
			female_parameters[5] = 0.8901 * chestsize / 100;
			//sprintf_s(Index, "%d", static_cast<int>(chestsize * 100));
			//s = Index;
			//s_filename = s_filename + "ChestLine" + s + "_";
			for (int waistsize_index = 0; waistsize_index < 3; ++waistsize_index)
			{
				//double waistsize = 0.7 + waistsize_index * 0.1;
				double waistsize = 70 + waistsize_index * 10;
				female_parameters[6] = waistsize / 100;
				female_parameters[7] = 1.2115 * waistsize / 100;
				//sprintf_s(Index, "%d", static_cast<int>(waistsize * 100));
				//s = Index;
				//s_filename = s_filename + "WaistLine" + s + "_";
				for (int hipsize_index = 0; hipsize_index < 3; ++hipsize_index)
				{
					double hipsize = 92 + hipsize_index * 9;
					female_parameters[8] = hipsize / 100;

					sprintf_s(Index, "%d", static_cast<int>(height));
					s = Index;
					s_filename = "Female_Body_Height" + s + "_";

					sprintf_s(Index, "%d", static_cast<int>(chestsize));
					s = Index;
					s_filename = s_filename + "ChestLine" + s + "_";

					sprintf_s(Index, "%d", static_cast<int>(waistsize));
					s = Index;
					s_filename = s_filename + "WaistLine" + s + "_";

					sprintf_s(Index, "%d", static_cast<int>(hipsize));
					s = Index;
					matfilename = matfilepath + s_filename + "HipLine" + s + ".mat";
					objfilename = objfilepath + s_filename + "HipLine" + s + ".obj";
					fullfilename = fullfilepath + s_filename + "HipLine" + s + ".obj";

					Woman.GenerateModel(female_parameters, 9, Twoman);
					MatrixXd WomanBodyVertex = Woman.GetBodyVertices();
					MatrixXi WomanBodyFaces = Woman.GetBodyFaces();

					Woman.SaveModel2Obj(objfilename.c_str(), WomanBodyVertex, WomanBodyFaces);
					Woman.WriteMatFile(matfilename.c_str(), WomanBodyVertex, WomanBodyFaces);
					Woman.SaveModel2Obj(fullfilename.c_str());
				}
			}
		}
	}
#endif

#if 0
	Spring Woman("MatlabFile\\Female_SCAPE_Parameters.mat");
	double female_parameters[9] = { 0.52, 0.81, 0.48, 1.55, 0.81, 0.67, 0.64, 0.67, 0.86 };
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
	Woman.SaveModel2Obj("MatlabFile\\WomanTest.obj", v, Color);

	MatrixXd WomanBodyVertex = Woman.GetBodyVertices();
	MatrixXi WomanBodyFaces = Woman.GetBodyFaces();
	MatrixXd WomanBodyColor = Woman.GetBodyColor();

	//Woman.WriteMatFile("MatlabFile\\datatest.mat", WomanBodyVertex, WomanBodyFaces);



	HumanModel::HumanBody womanbody;
	womanbody.SetVertices(WomanBodyVertex);
	womanbody.SetFaces(WomanBodyFaces);
	womanbody.SetColors(WomanBodyColor);
#endif


#if 0
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
	QueryPerformanceCounter(&t2);
	cout << (t2.QuadPart - t1.QuadPart)*1.0 / tc.QuadPart << endl;
#endif
	DWORD timeEnd = GetTickCount();
	cout << "t : " << timeEnd - timeStart << endl;
	return 0;
}
