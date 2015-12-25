#include "stdafx.h"
#include "HumanModel.h"
#include "Boundary.h"
#include <Eigen/SVD>
#include <fstream>
#include <iostream>
#include <engine.h>
#include <windows.h>
using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::RowVector3d;


HumanModel::HumanModel()
{}

HumanModel::HumanModel(HumanHead& HeadPart, HumanHead& BodyPart, Boundary& NowBoundary)
{
	ConvergeHeadAndBody(HeadPart, BodyPart,NowBoundary);
}

HumanModel::HumanModel(MatrixXd& Vertices, MatrixXi& Faces)
{
	SetHumanModel(Vertices, Faces);
}

HumanModel::~HumanModel()
{}

void HumanModel::SetHumanModel(MatrixXd& Vertices, MatrixXi& Faces)
{
	HumanFullVertices = Vertices;
	HumanFullFaces = Faces;
	NumberOfVertices = Vertices.rows();
	NumberOfFaces = Vertices.cols();
}


void HumanModel::ReadMatlabFile(const char* filename,HumanPart& part)
{
	typedef double MatDataFormat;
	MATFile* pmatFile = matOpen(filename,"r");
	if (pmatFile == NULL)
	{
		std::cerr << "Not find matlab file" << std::endl;
		//ReadFlag = false;
		return;
	}
	//else
	//	ReadFlag = true;
	mxArray* pVertices = matGetVariable(pmatFile,"Vertices");
	MatDataFormat* pdVertices = static_cast<MatDataFormat*> (mxGetData(pVertices));
	size_t M = mxGetM(pVertices);
	size_t N = mxGetN(pVertices);
	part.Vertices.resize(M, N);
	//for (size_t i = 0; i < M;++i)
	//{
	//	for (size_t j = 0; j < N;++j)
	//	{
	//		part.Vertices(i, j) = (pdVertices[j*M + i]);
	//	}
	//}

	memcpy(part.Vertices.data(), pdVertices, sizeof(MatDataFormat)* M * N); // zhy
	mxFree(pdVertices);
	part.NumberOfVertices = M;

	mxArray* pFace = matGetVariable(pmatFile, "Faces");
	MatDataFormat* pdFace = static_cast<MatDataFormat*> (mxGetData(pFace));
	M = mxGetM(pFace);
	N = mxGetN(pFace);
	part.Faces.resize(M, N);
	//for (size_t i = 0; i < M; ++i)
	//{
	//	for (size_t j = 0; j < N; ++j)
	//	{
	//		part.Faces(i, j) = (pdFace[j*M + i]);
	//	}
	//}
	Eigen::MatrixXd temp;
	temp.resize(M, N);
	memcpy(temp.data(), pdFace, sizeof(double)* M * N); // zhy
	part.Faces = temp.cast<int>();
	mxFree(pdFace);
	part.NumberOfFaces = M;
#if 0
	mxArray* pColor = matGetVariable(pmatFile, "Color");
	MatDataFormat* pdColor = static_cast<MatDataFormat*> (mxGetData(pColor));
	M = mxGetM(pColor);
	N = mxGetN(pColor);
	part.Colors.resize(M, N);
	//for (size_t i = 0; i < M; ++i)
	//{
	//	for (size_t j = 0; j < N; ++j)
	//	{
	//		part.Colors(i, j) = (pdColor[j*M + i]);
	//	}
	//}
	memcpy(part.Colors.data(), pdColor, sizeof(MatDataFormat)* M * N); // zhy
	mxFree(pdColor);
#endif
	matClose(pmatFile);
}

#if 0
void HumanModel::ReadObjFile(const char* filename,HumanPart& part)
{
	const size_t NumberOfVertices = 1443;
	const size_t NumberOfFaces = 2852;
	using std::ifstream;
	ifstream file;
	file.open(filename);

	char ch;
	double x(0), y(0), z(0);
	int a(0), b(0), c(0);
	part.Vertices.resize(NumberOfVertices, 3);
	part.Faces.resize(NumberOfFaces, 3);
	part.Colors.resize(NumberOfVertices, 3);

	size_t vi = 0,fi=0;
	while (!file.eof())
	{
		file >> ch;
		if (ch == 'v')
		{
			file >> part.Vertices(vi, 0) >> part.Vertices(vi, 1) >> part.Vertices(vi, 2) >> part.Colors(vi, 0) >> part.Colors(vi, 1) >> part.Colors(vi, 2);
			++vi;
			ch = 0;
		}
		else
		{
			if (ch == 'f')
			{
				file >> part.Faces(fi, 0) >> part.Faces(fi, 1) >> part.Faces(fi, 2);
				++fi;
				ch = 0;
			}
		}
	}
	file.close();
}
#endif


//const MatrixXd& HumanModel::GetVertices() const
//{
//	return HumanFullVertices;
//}

//const MatrixXi& HumanModel::GetFaces() const
//{
//	return HumanFullFaces;
//}

//MatrixXi HumanModel::GetBoundaryVertices(PartOfHuman Indication)
//{}

//size_t HumanModel::GetNumberOfVertices() const
//{
//	return NumberOfVertices;
//}

//size_t HumanModel::GetNumberOfFaces() const
//{
//	return NumberOfFaces;
//}

void HumanModel::ConvergeHeadAndBody(HumanHead& HeadPart, HumanHead& BodyPart, Boundary& NowBoundary)
{
	ConvergeHeadAndBodyStuff(HeadPart, BodyPart,NowBoundary);
}

void HumanModel::SaveModel2Obj(const char* filename)
{
	if (VNormal.cols() == 0)
		ComputeNormals();
	SaveModel2Obj(filename, HumanFullVertices, HumanFullFaces, VNormal);
}

void HumanModel::SaveModel2Obj(const char* filename, Eigen::MatrixXd& Vertices, Eigen::MatrixXi& Faces, MatrixXd& Normals)
{
	std::ofstream out;
	out.open(filename);
	size_t NowNumberOfVertices = Vertices.rows();
	size_t NowNumberOfFaces = Faces.rows();

	std::string buff;
	for (int i = 0; i < NowNumberOfVertices; ++i)
	{

		buff.append("v");
		for (int j = 0; j < Vertices.cols(); ++j)
		{
			buff.append(" " + std::to_string(Vertices(i, j)));
		}
		buff.append("\n");

		buff.append("vn");
		for (int j = 0; j < Normals.cols(); ++j)
		{
			buff.append(" " + std::to_string(Normals(i, j)));
		}
		buff.append("\n");
		//for (int j = 0; j < Colors.cols(); ++j)
		//{
		//	buff.append(" " + std::to_string(Colors(i, j)));
		//}
	}
	for (int i = 0; i < NowNumberOfFaces; ++i)
	{
		buff.append("f");
		for (int j = 0; j < Faces.cols(); ++j)
		{
			buff.append(" " + std::to_string(Faces(i, j)));
		}
		buff.append("\n");
	}
	/*size_t NowNumberOfNormal = VNormal.rows();
	for (int i = 0; i < NowNumberOfNormal; ++i)
	{
		buff.append("vn");
		for (int j = 0; j < VNormal.cols(); ++j)
		{
			buff.append(" " + std::to_string(VNormal(i, j)));
		}
		buff.append("\n");
	}*/
	out.write(buff.c_str(), buff.size());
	out.close();
}

double HumanModel::AdjustScale(Eigen::MatrixXd& HeadBoundaryVertices, Eigen::MatrixXd& BodyBoundaryVertices)
{
	size_t N = HeadBoundaryVertices.rows();
	double ScaleSum = 0;
	MatrixXd Line(1, 3);
	double LengthOfLineHead = 0;
	double LengthOfLineBody = 0;
	size_t ccc = 0;
#if 0
	std::cout << HeadBoundaryVertices << std::endl;
	std::cout << BodyBoundaryVertices << std::endl;
#endif
	for (size_t idx = 1; idx < N;++idx)
	{
		LengthOfLineHead = (HeadBoundaryVertices.row(idx) - HeadBoundaryVertices.row(idx - 1)).norm();
		//std::cout << HeadBoundaryVertices.row(idx) << std::endl << HeadBoundaryVertices.row(idx - 1) << std::endl;

		LengthOfLineBody = (BodyBoundaryVertices.row(idx) - BodyBoundaryVertices.row(idx - 1)).norm();
		//std::cout << Line << std::endl;
		if (LengthOfLineHead != 0)
		{
			ScaleSum += (LengthOfLineBody / LengthOfLineHead);
			++ccc;
		}
		//std::cout << LengthOfLineHead << std::endl;
	}
	return ScaleSum/ccc;
}

MatrixXd HumanModel::RigidTransformation(Eigen::MatrixXd& target, Eigen::MatrixXd& source)
{
	MatrixXd tar_center(1, 3);
	MatrixXd src_center(1, 3);
	//std::cout << "target:" << std::endl << target << std::endl;
	//std::cout << "source:" << std::endl << source << std::endl;
	for (size_t i = 0; i < 3; ++i)
	{
		tar_center(i) = (target.col(i)).sum() / target.rows();
		src_center(i) = (source.col(i)).sum() / source.rows();
	}
	MatrixXd src_new(source.rows(), source.cols());
	MatrixXd tar_new(target.rows(), target.cols());
	src_new = source - MatrixXd::Ones(source.rows(), 1) * tar_center;
	tar_new = target - MatrixXd::Ones(target.rows(), 1) * tar_center;
	//MatrixXd variance = tar_new.adjoint() * src_new;
	Eigen::JacobiSVD<MatrixXd> svd(tar_new.adjoint() * src_new, Eigen::ComputeFullU | Eigen::ComputeFullV);
	MatrixXd R(3, 3);
	R = svd.matrixU() * (svd.matrixV().adjoint());
	MatrixXd t(3, 1);
	t = tar_center.adjoint() - R * src_center.adjoint();
	MatrixXd M(3, 4);
	M.col(0) = R.col(0);
	M.col(1) = R.col(1);
	M.col(2) = R.col(2);
	M.col(3) = t.col(0);
	//std::cout << M << std::endl;
	return M;
}


void HumanModel::TransformVertices(MatrixXd& Vertices, MatrixXd R, MatrixXd t)
{
	Vertices = (R*Vertices.adjoint()).adjoint() + Eigen::MatrixXd::Ones(Vertices.rows(), 1) * t.adjoint();
}

MatrixXd HumanModel::InterpolationForConverge(Eigen::MatrixXd CenterVertex_H, MatrixXd& CenterVertex_B, MatrixXd& Neighbors)
{
	size_t NumberOfNeighbors = Neighbors.rows();
	MatrixXd Res(1, 3);
	if (NumberOfNeighbors < 4)
	{
		for (int i = 0; i < 3; ++i)
			Res(i) = (Neighbors.col(i)).sum() / NumberOfNeighbors;
		return Res;
	}
	MatrixXd u = Neighbors.col(0);
	MatrixXd v = Neighbors.col(1);
	MatrixXd s = Neighbors.col(2);
	MatrixXd A(NumberOfNeighbors, 6);
	A.col(0).setOnes();
	A.col(1) = u.col(0);	//any problem??
	A.col(2) = v.col(0);
	A.col(3) = (u.cwiseProduct(u)).col(0);
	A.col(4) = (u.cwiseProduct(v)).col(0);
	A.col(5) = (v.cwiseProduct(v)).col(0);
	MatrixXd a = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(s.col(0));
	//std::cout << a << std::endl;
	//a = (A.adjoint() * A).inverse() * A.adjoint() * s.col(0);
	MatrixXd VertexC(1, 3);
	VertexC = (CenterVertex_H + CenterVertex_B) / 2;
	MatrixXd VerCoef(1, 6);
	VerCoef(0) = 1;
	VerCoef(1) = VertexC(0);
	VerCoef(2) = VertexC(1);
	VerCoef(3) = VertexC(0) * VertexC(0);
	VerCoef(4) = VertexC(0) * VertexC(1);
	VerCoef(5) = VertexC(1) * VertexC(1);
	MatrixXd res = VerCoef * a;
	MatrixXd VertexV = VertexC;
	VertexV(2) = res(0);
	if ((VertexC - VertexV).norm() > 0.01)
		return VertexC;
	else
		return VertexV;
}

void HumanModel::ConvergeHeadAndBodyStuff(HumanHead& HeadPart, HumanHead& BodyPart, Boundary& NowBoundary)
{
	size_t NumberOfBoundaryVertex = NowBoundary.GetNumberOfBoundaryVertex();
	MatrixXd BoundaryVerticesInHead(NumberOfBoundaryVertex, 3);
	MatrixXd BoundaryVerticesInBody(NumberOfBoundaryVertex, 3);
	int BoundaryIndexInHead = 0;
	int BoundaryIndexInBody = 0;
	for (int idx = 0; idx < NumberOfBoundaryVertex;++idx)
	{
		BoundaryIndexInHead = NowBoundary.GetBoundaryVertexIndex(idx, Boundary::Head) - 1;
		BoundaryIndexInBody = NowBoundary.GetBoundaryVertexIndex(idx, Boundary::Body) - 1;
		for (size_t idy = 0; idy < 3;++idy)
		{
			BoundaryVerticesInHead(idx, idy) = HeadPart.Vertices(BoundaryIndexInHead, idy);
			BoundaryVerticesInBody(idx, idy) = BodyPart.Vertices(BoundaryIndexInBody, idy);
		}
		/*std::cout <<"Head: " << BoundaryVerticesInHead.row(idx) << std::endl;
		std::cout << "Body: " << BoundaryVerticesInBody.row(idx) << std::endl;*/
	}
	double scale = AdjustScale(BoundaryVerticesInHead,BoundaryVerticesInBody);
	//std::cout << '1' <<std::endl<< BoundaryVerticesInHead - BoundaryVerticesInBody << std::endl;
	//std::cout << "BoundaryVerticesInHead:" << std::endl << BoundaryVerticesInHead << std::endl;
	BoundaryVerticesInHead *= scale;
	
	MatrixXd T = RigidTransformation(BoundaryVerticesInBody, BoundaryVerticesInHead);

	//std::cout << T << std::endl << scale << std::endl;
	MatrixXd t(3, 1);
	t = T.col(3);

	MatrixXd R(3, 3);

	R.col(0) = T.col(0);
	R.col(1) = T.col(1);
	R.col(2) = T.col(2);
	//std::cout << R << std::endl;
	//std::cout << BoundaryVerticesInHead.rows() << ' ' << BoundaryVerticesInHead.cols() << std::endl;
	BoundaryVerticesInHead = (R * BoundaryVerticesInHead.adjoint()).adjoint() + Eigen::MatrixXd::Ones(BoundaryVerticesInHead.rows(), 1) * t.adjoint();	//已经乘过scale,/(ㄒoㄒ)/~~这个bug调太久！！！！！
	R *= scale;

	//std::cout << R << std::endl;
	//std::cout << '2' << std::endl << BoundaryVerticesInHead - BoundaryVerticesInBody << std::endl;
	//std::cout << t << std::endl;

	TransformVertices(HeadPart.Vertices,R,t);
	
	//SaveModel2Obj("MatlabFile\\1.obj", HeadPart.Vertices, HeadPart.Faces);

	MatrixXd NeighborVertices;
	VertexAndNeighbors head;
	VertexAndNeighbors body;
	int IndexInHead = 0;
	int IndexInBody = 0;
	MatrixXd CenterVertexInHead(1, 3);
	MatrixXd CenterVertexInBody(1, 3);
	MatrixXd AllInterpolationVertex(NumberOfBoundaryVertex,3);
	for (int idx = 0; idx < NumberOfBoundaryVertex;++idx)
	{
		head = NowBoundary.GetVertexAndNeighbors(idx, Boundary::Head);
		body = NowBoundary.GetVertexAndNeighbors(idx, Boundary::Body);
		NeighborVertices.resize(head.NumberOfNeighborVertex + body.NumberOfNeighborVertex, 3);
		for (int i = 0; i < head.NumberOfNeighborVertex;++i)
		{
			IndexInHead = (head.IndexOfNeighborVertex)[i] - 1;
			NeighborVertices.row(i) = (HeadPart.Vertices).row(IndexInHead);
		}
		CenterVertexInHead = BoundaryVerticesInHead.row(idx);		//应该没取错，应该是已经变换了的坐标
		for (int j = 0; j < body.NumberOfNeighborVertex;++j)
		{
			IndexInBody = (body.IndexOfNeighborVertex)[j] - 1;
			NeighborVertices.row(head.NumberOfNeighborVertex + j) = (BodyPart.Vertices).row(IndexInBody);
		}
		CenterVertexInBody = BoundaryVerticesInBody.row(idx);
		//std::cout << CenterVertexInHead << std::endl << CenterVertexInBody << std::endl;
		AllInterpolationVertex.row(idx) = InterpolationForConverge(CenterVertexInHead, CenterVertexInBody, NeighborVertices);
		IndexInHead = NowBoundary.GetBoundaryVertexIndex(idx, Boundary::Head) - 1;
		IndexInBody = NowBoundary.GetBoundaryVertexIndex(idx, Boundary::Body) - 1;
		HeadPart.Vertices.row(IndexInHead) = AllInterpolationVertex.row(idx);
		BodyPart.Vertices.row(IndexInBody) = AllInterpolationVertex.row(idx);
	}
	AdjustHeadEdgeVertex(HeadPart, NowBoundary);
	PlaceHeadAndBodyToModel(HeadPart, BodyPart,NowBoundary);
}

void HumanModel::PlaceHeadAndBodyToModel(HumanHead& HeadPart, HumanBody& BodyPart, Boundary& NowBoundary)
{
	NumberOfVertices = HeadPart.NumberOfVertices + BodyPart.NumberOfVertices - NowBoundary.GetNumberOfBoundaryVertex();	//Head和Body中有重合点（边界点），如果不合并，边界处会有问题
	HumanFullVertices.resize(NumberOfVertices, 3);
	NumberOfFaces = HeadPart.NumberOfFaces + BodyPart.NumberOfFaces;
	HumanFullFaces.resize(NumberOfFaces, 3);
	//Colors.resize(NumberOfVertices, 3);
	for (size_t i = 0; i < HeadPart.NumberOfVertices; ++i)
	{
		HumanFullVertices.row(i) = HeadPart.Vertices.row(i);
		//Colors.row(i) = HeadPart.Colors.row(i);
	}
	for (size_t i = 0; i < HeadPart.NumberOfFaces; ++i)
	{
		HumanFullFaces.row(i) = HeadPart.Faces.row(i);
	}

	MatrixXi& BodyMaptable = NowBoundary.GetBodyMapTable(HeadPart.NumberOfVertices, BodyPart.NumberOfVertices);
	for (size_t i = 0; i < BodyPart.NumberOfVertices; ++i)
	{
		HumanFullVertices.row(BodyMaptable(i) - 1) = BodyPart.Vertices.row(i);
		//Colors.row(BodyMaptable(i) - 1) = BodyPart.Colors.row(i);
	}
	MatrixXi a;
	for (size_t i = 0; i < BodyPart.NumberOfFaces; ++i)
	{
		a = BodyPart.Faces.row(i);
		a(0) = BodyMaptable(a(0) - 1);
		a(1) = BodyMaptable(a(1) - 1);
		a(2) = BodyMaptable(a(2) - 1);
		HumanFullFaces.row(HeadPart.NumberOfFaces + i) = a.row(0);
	}
}

void HumanModel::AdjustHeadEdgeVertex(HumanHead& HeadPart, Boundary& NowBoundary)
{
	NowBoundary.ReadBoundaryWeightFile("E:\\VirtualFitting\\generateModel\\VirtualFittingDataset\\MatlabFile\\BoundaryWeight.mat");
	const MatrixXd& Weight = NowBoundary.GetBoundaryWeight();
#if 1
	for (int i = 0; i < Weight.rows(); ++i)
	{
		HeadPart.Vertices.row(static_cast<int>(Weight(i, 0) - 1)) = Eigen::RowVector3d(0, 0, 0);
		for (int j = 0; j < 3;++j)
		{
			if (static_cast<int>(Weight(i, j + 1)) == 0)
				break;
			HeadPart.Vertices.row(static_cast<int>(Weight(i, 0) - 1)) += (Weight(i, j + 4) * HeadPart.Vertices.row(static_cast<int>(Weight(i, j + 1) - 1)));
		}		
	}
#endif
}



//void HumanModel::SetHumanPart(const MatrixXd& Vertices, HumanPart& part)
//{
//	HumanPart.Vertices
//}
//
//void HumanModel::SetHumanPart(const MatrixXi& Faces, HumanPart& part)
//{}
//
//void HumanModel::SetHumanPart(const MatrixXd& Vertices, const MatrixXi& Faces, const MatrixXd& Colors, HumanPart& part)
//{}


void HumanModel::ComputeNormals()
{
	_ComputeFaceNormal();
	_ComputeVertexNormal();
}

MatrixXd HumanModel::MatrixRowCross(MatrixXd& left, MatrixXd& right)
{
	MatrixXd res(1, 3);
	res(0) = left(1)*right(2) - left(2)*right(1);
	res(1) = left(2)*right(0) - left(0)*right(2);
	res(2) = left(0)*right(1) - left(1)*right(0);
	return res;
}

void HumanModel::_ComputeFaceNormal()
{
	FNormal.resize(NumberOfFaces, 3);
	MatrixXd u;
	MatrixXd x,y,z,t1,t2;

	for (size_t i_f_row = 0; i_f_row < NumberOfFaces; ++i_f_row)
	{
		x = HumanFullVertices.row(HumanFullFaces(i_f_row, 0) - 1);
		y = HumanFullVertices.row(HumanFullFaces(i_f_row, 1) - 1) - x;
		z = HumanFullVertices.row(HumanFullFaces(i_f_row, 2) - 1) - x;
		u = MatrixRowCross(y,z);
		u = u / u.norm();
		FNormal.row(i_f_row) = u;
	}
}

void HumanModel::_ComputeVertexNormal()
{
	const char filename[] = "E:\\VirtualFitting\\generateModel\\VirtualFittingDataset\\MatlabFile\\Faces_Neighbors.mat";
	MATFile* pmatFile = matOpen(filename, "r");
	if (pmatFile == NULL)
	{
		std::cerr << "Not find Vertex_Neighbors.mat file" << std::endl;
		return;
	}
	mxArray* pNeighbors = matGetVariable(pmatFile, "Neighbors");
	double* pdNeighbors = static_cast<double*>(mxGetData(pNeighbors));
	size_t M = mxGetM(pNeighbors);
	size_t N = mxGetN(pNeighbors);
	MatrixXd Neighbors;
	Neighbors.resize(M, N);
	memcpy(Neighbors.data(), pdNeighbors, sizeof(double)* M * N);
	mxFree(pdNeighbors);

	VNormal.resize(NumberOfVertices, 3);
	MatrixXd u(1, 3);
	for (size_t i_v_row = 0; i_v_row < NumberOfVertices; ++i_v_row)
	{
		u.setZero();
		for (size_t index_iter = 0; index_iter < Neighbors(i_v_row, 19); ++index_iter)
			u = u + FNormal.row(static_cast<int>(Neighbors(i_v_row, index_iter) - 1));
		u = u / u.norm();
		VNormal.row(i_v_row) = u;
	}
}