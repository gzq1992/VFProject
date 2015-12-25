#include "stdafx.h"
#include "Spring.h"
//#include <vector>
#include <iostream>
#include <fstream>
#include <Eigen/SparseCholesky>

//#include <vector>
#include <windows.h>
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::MatrixXi;
using Eigen::SimplicialLDLT;
using std::cerr;
using std::cout;
using std::endl;


Spring::Spring(char* SCAPE_Parameters_Filename)
{
	Init(SCAPE_Parameters_Filename);
}


Spring::~Spring()
{
	delete[] SCAPE_Template.Paramters_2_Shape;
	delete[] SCAPE_Template.SCAPE_P;
	SCAPE_Template.Paramters_2_Shape = NULL;
	SCAPE_Template.SCAPE_P = NULL;
}

void Spring::SaveModel2Obj(const char* filename)
{
	SaveModel2Obj(filename, vout, SCAPE_Template.Faces);
}

void Spring::SaveModel2Obj(const char* filename, const MatrixXd& Vertices, const MatrixXi& Faces)
{
#if 0
	std::ofstream out;
	out.open(filename);
	size_t NumberOfVertices = SCAPE::Number_Of_Vertices;
	size_t NumberOfFaces = SCAPE::Number_Of_Faces;
	std::string buff;

	for (size_t i = 0; i < NumberOfVertices; ++i)
	{
		buff.append("v");
		for (int j = 0; j < vout.cols(); ++j)
		{
			buff.append(" " + std::to_string(vout(i, j)));
		}
		buff.append("\n");
	}
	for (size_t i = 0; i < NumberOfFaces; ++i)
	{
		buff.append("f");
		for (int j = 0; j < SCAPE_Template.Faces.cols(); ++j)
		{
			buff.append(" " + std::to_string(SCAPE_Template.Faces(i, j)));
		}
		buff.append("\n");
	}
	out.write(buff.c_str(), buff.size());
	out.close();
#endif
	std::ofstream out;
	out.open(filename);
	size_t NumberOfVertices = Vertices.rows();
	size_t NumberOfFaces = Faces.rows();
	std::string buff;

	for (size_t i = 0; i < NumberOfVertices; ++i)
	{
		buff.append("v");
		for (int j = 0; j < Vertices.cols(); ++j)
		{
			buff.append(" " + std::to_string(Vertices(i, j)));
		}
		buff.append("\n");
	}
	for (size_t i = 0; i < NumberOfFaces; ++i)
	{
		buff.append("f");
		for (int j = 0; j < Faces.cols(); ++j)
		{
			buff.append(" " + std::to_string(Faces(i, j)));
		}
		buff.append("\n");
	}
	out.write(buff.c_str(), buff.size());
	out.close();
}

#if 0
void Spring::ReadSegmentRef(char* RefFilename)
{
	MATFile* pmatFile = matOpen(RefFilename, "r");
	if (pmatFile == NULL)
	{
		cerr << "Not find " << RefFilename << endl;
		return;
	}

	mxArray* ptr_Index_BodyVertex_InSCAPE = matGetVariable(pmatFile, "Index_BodyVertex_InSCAPE");
	Read2DMatData(ptr_Index_BodyVertex_InSCAPE, SegmentIndex.Index_BodyVertex_InSCAPE);
	
	mxArray* ptr_Index_HeadVertex_InSCAPE = matGetVariable(pmatFile, "Index_HeadVertex_InSCAPE");
	Read2DMatData(ptr_Index_HeadVertex_InSCAPE, SegmentIndex.Index_HeadVertex_InSCAPE);

	mxArray* ptr_Index_BoundaryVertex_InBody = matGetVariable(pmatFile, "Index_BoundaryVertex_InBody");
	Read2DMatData(ptr_Index_BoundaryVertex_InBody, SegmentIndex.Index_BoundaryVertex_InBody);

	mxArray* ptr_Index_BoundaryVertex_InHead = matGetVariable(pmatFile, "Index_BoundaryVertex_InHead");
	Read2DMatData(ptr_Index_BoundaryVertex_InHead, SegmentIndex.Index_BoundaryVertex_InHead);
	
	matClose(pmatFile);
}
#endif

void Spring::Init(char* SCAPE_Parameters_Filename)
{
	typedef double MatDataFormat;
	MATFile* pmatFile = matOpen(SCAPE_Parameters_Filename,"r");
	if (pmatFile == NULL)
	{
		cerr << "Not find " << SCAPE_Parameters_Filename << endl;
		return;
	}
	mxArray* ptr_Paramters_2_Shape = matGetVariable(pmatFile, "Parameter_2_Shape");
	Read3DMatData(ptr_Paramters_2_Shape,SCAPE_Template.Paramters_2_Shape,SCAPE::Paramters_2_Shape_sizex,SCAPE::Paramters_2_Shape_sizey,SCAPE::Number_Of_Faces);

	mxArray* ptr_Face_In_Bone = matGetVariable(pmatFile, "Face_In_Bone");
	Read2DMatData(ptr_Face_In_Bone, SCAPE_Template.Face_In_Bone);
	//cout << SCAPE_Template.Face_In_Bone.row(3) << endl << SCAPE_Template.Face_In_Bone.row(5) << endl;

	mxArray* ptr_Template = matGetVariable(pmatFile, "Template");
	Read2DMatData(ptr_Template, SCAPE_Template.Template);

	mxArray* ptr_SCAPE_P = matGetVariable(pmatFile, "SCAPE_P");
	Read3DMatData(ptr_SCAPE_P, SCAPE_Template.SCAPE_P,SCAPE::SCAPE_P_sizex,SCAPE::SCAPE_P_sizey,SCAPE::Number_Of_Faces);

	mxArray* ptr_SCAPE_M = matGetVariable(pmatFile, "SCAPE_M_Index");
	ReadSparseData(ptr_SCAPE_M, SCAPE_Template.SCAPE_M, SCAPE::SCAPE_M_rows, SCAPE::SCAPE_M_cols);

	mxArray* ptr_SCAPE_N = matGetVariable(pmatFile, "SCAPE_N_Index");
	SCAPE_Template.SCAPE_N.resize(SCAPE::SCAPE_N_rows, SCAPE::SCAPE_N_cols);
	ReadSparseData(ptr_SCAPE_N, SCAPE_Template.SCAPE_N, SCAPE::SCAPE_N_rows, SCAPE::SCAPE_N_cols);

	mxArray* ptr_Faces = matGetVariable(pmatFile, "Faces");
	Read2DMatData(ptr_Faces, SCAPE_Template.Faces);

	mxArray* ptr_Bone_Of_Face = matGetVariable(pmatFile, "Bone_Of_Face");
	Read2DMatData(ptr_Bone_Of_Face, SCAPE_Template.Bone_Of_Face);

	mxArray* ptr_Bone_Of_Vertex = matGetVariable(pmatFile, "Bone_Of_Vertex");
	Read2DMatData(ptr_Bone_Of_Vertex, SCAPE_Template.Bone_Of_Vertex);

	mxArray* ptr_Bone_Conjunctions = matGetVariable(pmatFile, "Bone_Conjunctions");
	Read2DMatData(ptr_Bone_Conjunctions, SCAPE_Template.Bone_Conjunctions);
	//cout << endl << SCAPE_Template.Bone_Conjunctions << endl << endl;

	mxArray* ptr_T = matGetVariable(pmatFile, "T");
	Read3DMatData(ptr_T, SCAPE_Template.T, 4, 4, SCAPE::Number_Of_Bones);

	mxArray* ptr_Joints = matGetVariable(pmatFile, "Joints");
	Read2DMatData(ptr_Joints, SCAPE_Template.Joints);

	mxArray* ptr_Body = matGetVariable(pmatFile, "Body");
	Read2DMatData(ptr_Body, SCAPE_Template.Body);

	mxArray* ptr_Index_BodyVertex_InSCAPE = matGetVariable(pmatFile, "Index_BodyVertex_InSCAPE");
	Read2DMatData(ptr_Index_BodyVertex_InSCAPE, SCAPE_Template.Index_BodyVertex_InSCAPE);

	mxArray* ptr_BodyFaces = matGetVariable(pmatFile, "BodyFaces");
	Read2DMatData(ptr_BodyFaces, SCAPE_Template.BodyFaces);

	mxArray* ptr_BodyColor = matGetVariable(pmatFile, "BodyColor");
	Read2DMatData(ptr_BodyColor, SCAPE_Template.BodyColor);

	matClose(pmatFile);
}




void Spring::Read2DMatData(mxArray* MatData, MatrixXd& C2DMatrixData)
{
	MatDataFormat* ptr_data = static_cast<MatDataFormat*> (mxGetData(MatData));
	size_t M = mxGetM(MatData);
	size_t N = mxGetN(MatData);
	C2DMatrixData.resize(M, N);
	//for (size_t i = 0; i < N;++i)
	//{
	//	size_t x = i * M;
	//	for (size_t j = 0; j < M; ++j)
	//	{
	//		C2DMatrixData(j, i) = (ptr_data[x + j]);
	//	}
	//}
	memcpy(C2DMatrixData.data(), ptr_data, sizeof(double)* M * N);
	mxFree(ptr_data);
	//cout << "M:" << M << ' ' << "N:" << N << endl;
}

void Spring::Read2DMatData(mxArray* MatData, MatrixXi& C2DMatrixData)
{
	MatDataFormat* ptr_data = static_cast<MatDataFormat*> (mxGetData(MatData));
	size_t M = mxGetM(MatData);
	size_t N = mxGetN(MatData);
	C2DMatrixData.resize(M, N);
	//for (size_t i = 0; i < N; ++i)
	//{
	//	size_t x = i * M;
	//	for (size_t j = 0; j < M; ++j)
	//	{
	//		C2DMatrixData(j, i) = static_cast<int>(ptr_data[x + j]);
	//	}
	//}
	//cout << C2DMatrixData.row(1) << endl;
	Eigen::MatrixXd temp;
	temp.resize(M, N);
	memcpy(temp.data(), ptr_data, sizeof(double)* M * N);
	C2DMatrixData = temp.cast<int>();
	mxFree(ptr_data);

}


void Spring::Read3DMatData(mxArray* MatData, MatrixXd*& C3DMatrixData, size_t size_x, size_t size_y, size_t size_z)
{
	MatDataFormat* dataptr = static_cast<MatDataFormat*> (mxGetData(MatData));
	C3DMatrixData = new MatrixXd[size_z];
	MatrixXd ReadData(size_x, size_y);
	size_t offset = 0;
	//for (size_t k = 0; k < size_z;++k)
	//{
	//	for (size_t j = 0; j < size_y;++j)
	//	{
	//		for (size_t i = 0; i < size_x;++i)
	//		{
	//			ReadData(i, j) = dataptr[offset];
	//			++offset;
	//		}
	//	}
	//	C3DMatrixData[k] = ReadData;
	//}
	for (size_t k = 0; k < size_z; ++k)
	{
		C3DMatrixData[k].resize(size_x, size_y);
		memcpy(C3DMatrixData[k].data(), dataptr + size_x * size_y * k, size_x * size_y * sizeof(double));
	}
	mxFree(dataptr);
}

void Spring::ReadSparseData(mxArray* MatData, SparseMatrix<double>& CSparseData,size_t rows,size_t cols)
{
	MatDataFormat* dataptr = static_cast<MatDataFormat*>(mxGetData(MatData));
	size_t M = mxGetM(MatData);
	size_t N = mxGetN(MatData);
	int* r = new int[N];
	int* c = new int[N];
	double* value = new double[N];
	size_t offset = 0;
	for (size_t i = 0; i < N;++i)
	{
		r[i] = static_cast<int> (dataptr[offset]);
		++offset;
		c[i] = static_cast<int> (dataptr[offset]);
		++offset;
		value[i] = dataptr[offset];
		++offset;
	}
	mxFree(dataptr);
	CSparseData.resize(static_cast<int>(rows),static_cast<int>(cols));
	Triplet<double>* ptr_triplet = new Triplet<double>[N];
	for (size_t i = 0; i < N;++i)
	{
		ptr_triplet[i] = Triplet<double>(r[i], c[i], value[i]);
	}
	CSparseData.setFromTriplets(ptr_triplet, ptr_triplet + N);

	delete[] ptr_triplet;
	delete[] r;
	delete[] c;
	delete[] value;
	ptr_triplet = NULL;
	r = NULL;
	c = NULL;
	value = NULL;
}

MatrixXd Spring::SCAPE_vector_from_motion(MatrixXd& Vertices, MatrixXi& Faces, MatrixXi& Bone_Of_Face, MatrixXi& Bone_Conjunctions, MatrixXd& Joints, size_t face, MatrixXd& Theta)
{
	MatrixXd T(1, 7);
	T(6) = 1;
	int bone = Bone_Of_Face(face) - 1;
	if (Bone_Conjunctions(bone,4) == 1)
	{
		for (int i = 0; i < 3;++i)
		{
			T(i) = Theta(Bone_Conjunctions(bone, 0) - 1, i);
			T(i + 3) = Theta(Bone_Conjunctions(bone, 0) - 1, i);
		}
		return T;
	}

	if (Bone_Conjunctions(bone, 4) == 2)
	{
		for (int i = 0; i < 3; ++i)
		{
			T(i) = Theta(Bone_Conjunctions(bone, 0) - 1, i);
			T(i + 3) = Theta(Bone_Conjunctions(bone, 1) - 1, i);
		}
		return T;
	}
	MatrixXd distance(Bone_Conjunctions(bone, 4), 1);
	MatrixXi v = Faces.row(face);
	MatrixXd center = Vertices.row(v(0) - 1) + Vertices.row(v(1) - 1) + Vertices.row(v(2) - 1);
	center = center / 3;
	for (int i = 0; i < Bone_Conjunctions(bone, 4);++i)
	{
		distance(i) = (Joints.row(Bone_Conjunctions(bone, i) - 1) - center).norm();
	}

	int min_location = 0;
	int min_location2 = 1;
	double min_val = distance(0), min_val2 = distance(1);
	if (distance(0) > distance(1))
	{
		min_location = 1;
		min_location2 = 0;
		min_val = distance(1);
		min_val2 = distance(0);
	}
	for (int i = 2; i < distance.rows();++i)
	{
		if (distance(i) < distance(min_location))
		{
			min_location2 = min_location;
			min_location = i;
			min_val2 = min_val;
			min_val = distance(i);
		}
		else
		{
			if (distance(i) < distance(min_location2))
			{
				min_location2 = i;
				min_val2 = distance(i);
			}
		}
	}
	for (int i = 0; i < 3; ++i)
	{
		T(i) = Theta(Bone_Conjunctions(bone, min_location2) - 1, i);
		T(i + 3) = Theta(Bone_Conjunctions(bone, min_location) - 1, i);
	}
	return T;
}


Spring& Spring::GenerateModel(double parameters[], size_t NumberOfParameters,MatrixXd T, MatrixXd& theta)
{
	if (parameters != NULL)
	{
		GenerateSpringModel(parameters, NumberOfParameters, T, theta);
	}
	return *this;
}

Spring& Spring::GenerateModel(double parameters[], size_t NumberOfParameters, MatrixXd T)
{
	MatrixXd theta(15, 3);
	theta.setZero();
	GenerateModel(parameters, NumberOfParameters, T, theta);
	return *this;
}


void Spring::GenerateSpringModel(double p[], size_t NumberOfParameters, MatrixXd T,MatrixXd& theta)
{

	MatrixXd Parameters(NumberOfParameters + 1, 1);
	for (size_t i = 0; i < NumberOfParameters; ++i)
	{
		Parameters(i) = p[i];
	}
	Parameters(NumberOfParameters) = 1;
	
	MatrixXd* S = new MatrixXd[SCAPE::Number_Of_Faces];

	for (size_t face = 0; face < SCAPE::Number_Of_Faces;++face)
	{
		MatrixXd sv = SCAPE_Template.Paramters_2_Shape[face] * Parameters;
		MatrixXd temp(3, 3);
		for (size_t i = 0; i < 9;++i)
		{
			temp(i) = sv(i);
		}
		S[face] = temp.adjoint();
	}

	MatrixXd bx(SCAPE::Number_Of_Faces * 2 + 1, 1);
	MatrixXd by(SCAPE::Number_Of_Faces * 2 + 1, 1);
	MatrixXd bz(SCAPE::Number_Of_Faces * 2 + 1, 1);
	
	//LARGE_INTEGER t1, t2, tc;
	//QueryPerformanceFrequency(&tc);
	//QueryPerformanceCounter(&t1);
	MatrixXd v1,v2,v3;
	MatrixXd Q_(3, 3);
	//MatrixXd theta(15, 3);
	//theta.setZero();
	MatrixXd pv, q, Q, v2_, v3_;
	for (size_t bone = 0; bone < SCAPE::Number_Of_Bones; ++bone)
	{
		MatrixXd R(3, 3);
		for (size_t i = 0; i < 3; ++i)
		{
			for (size_t j = 0; j < 3; ++j)
			{
				R(i, j) = SCAPE_Template.T[bone](i, j);
			}
		}
		for (size_t face = 0; face < SCAPE::Number_Of_Faces; ++face)
		{
			if (SCAPE_Template.Face_In_Bone(face, bone) == 0)
				continue;
			size_t tri = face;
#if 1

			v1 = SCAPE_Template.Template.row(SCAPE_Template.Faces(tri, 0) - 1);
			v2 = SCAPE_Template.Template.row(SCAPE_Template.Faces(tri, 1) - 1) - v1;
			v3 = SCAPE_Template.Template.row(SCAPE_Template.Faces(tri, 2) - 1) - v1;
#endif
			pv = SCAPE_vector_from_motion(SCAPE_Template.Template, SCAPE_Template.Faces, SCAPE_Template.Bone_Of_Face, SCAPE_Template.Bone_Conjunctions, SCAPE_Template.Joints, tri, theta);
			q = SCAPE_Template.SCAPE_P[tri] * pv.adjoint();

			for (int i = 0; i < 9; ++i)
			{
				Q_(i) = q(i);
			}
#if 1
			Q = Q_.adjoint();
			v2_ = R * S[tri] * Q * v2.adjoint();
			v3_ = R * S[tri] * Q * v3.adjoint();
#endif

			bx(tri * 2) = v2_(0);
			bx(tri * 2 + 1) = v3_(0);

			by(tri * 2) = v2_(1);
			by(tri * 2 + 1) = v3_(1);

			bz(tri * 2) = v2_(2);
			bz(tri * 2 + 1) = v3_(2);
		}
	}

	//QueryPerformanceCounter(&t2);
	//cout << (t2.QuadPart - t1.QuadPart)*1.0 / tc.QuadPart << endl;
	
	v1.resize(1, 4);
	for (int i = 0; i < 3;++i)
	{
		v1(i) = SCAPE_Template.Template(6316, i);
	}
	v1(3) = 1;

	MatrixXd v1_ = SCAPE_Template.T[SCAPE_Template.Bone_Of_Vertex(6316)] * v1.adjoint();
	MatrixXd v1_3(3,1);
	for (int i = 0; i < 3;++i)
	{
		v1_3(i) = v1_(i) / v1_(3);
	}

	bx(SCAPE::Number_Of_Faces * 2) = v1_3(0);
	by(SCAPE::Number_Of_Faces * 2) = v1_3(1);
	bz(SCAPE::Number_Of_Faces * 2) = v1_3(2);

	//QueryPerformanceFrequency(&tc);
	//QueryPerformanceCounter(&t1);
	SimplicialLDLT<SparseMatrix<double>> solver;
	solver.compute(SCAPE_Template.SCAPE_M);
	if (solver.info()!=Eigen::Success)
	{
		cerr << "Decomposition failed" << endl;
		return;
	}


	vout.resize(SCAPE::Number_Of_Vertices, 3);
	vout.col(0) = solver.solve(SCAPE_Template.SCAPE_N * bx);
	vout.col(1) = solver.solve(SCAPE_Template.SCAPE_N * by);
	vout.col(2) = solver.solve(SCAPE_Template.SCAPE_N * bz);
	//QueryPerformanceCounter(&t2);
	//cout << (t2.QuadPart - t1.QuadPart)*1.0 / tc.QuadPart << endl;
	MatrixXd _vout(SCAPE::Number_Of_Vertices, 4);
	for (int i = 0; i < 3;++i)
	{
		_vout.col(i) = vout.col(i);
	}
	_vout.col(3).setOnes();
	MatrixXd _vout_ = (T * _vout.adjoint()).adjoint();
	for (int i = 0; i < 3;++i)
	{
		vout.col(i) = _vout_.col(i);
	}
}

#if 0
void Spring::ComputeNormal(MatrixXd& vertex, MatrixXi& faces)
{
	size_t nface = faces.rows();
	size_t nvert = faces.rows();
	normal.resize(3,nvert);
	normal.setZero();
	MatrixXd normalf(SCAPE::Number_Of_Faces,3);
	MatrixXd x(3, 1);
	MatrixXd y(3, 1);
	for (size_t idf = 0; idf < SCAPE::Number_Of_Faces; ++idf)
	{
		x = vertex.row(faces(idf, 1)) - vertex.row(faces(idf, 0));
		y = vertex.row(faces(idf, 2)) - vertex.row(faces(idf, 0));
		normalf(idf, 0) = x(1) * y(2) - x(2) * y(1);
		normalf(idf, 1) = x(2) * y(0) - x(0) * y(2);
		normalf(idf, 2) = x(0) * y(1) - x(1) * y(0);
	}
}
#endif

void Spring::WriteMatFile(const char* filename, MatrixXd& Data, char* varName)
{
	size_t M = Data.rows();
	size_t N = Data.cols();

	double* outData = new double[M * N];
	for (size_t j = 0; j < N; ++j)
	{
		for (size_t i = 0; i < M; ++i)
		{
			outData[j*M + i] = Data(i, j);
		}
	}
	MATFile* pmatfile = matOpen(filename, "w");
	mxArray *pMxArray = mxCreateDoubleMatrix(M, N, mxREAL);
	memcpy(mxGetPr(pMxArray), outData, M*N*sizeof(double));
	matPutVariable(pmatfile, varName, pMxArray);
	matClose(pmatfile);
	delete[] outData;
	outData = NULL;
}

void Spring::WriteMatFile(const char* filename, MatrixXd& Vertices, MatrixXi& Faces)
{
	MATFile* pmatfile = matOpen(filename, "w");
	mxArray *pMxArray = NULL;
	size_t M = Vertices.rows();
	size_t N = Vertices.cols();
	double* outData = new double[M * N];
	for (size_t j = 0; j < N; ++j)
	{
		for (size_t i = 0; i < M; ++i)
		{
			outData[j*M + i] = Vertices(i, j);
		}
	}
	pMxArray = mxCreateDoubleMatrix(M, N, mxREAL);
	memcpy(mxGetPr(pMxArray), outData, M*N*sizeof(double));
	matPutVariable(pmatfile, "Vertices", pMxArray);
	delete[] outData;

	M = Faces.rows();
	N = Faces.cols();
	outData = new double[M * N];
	for (size_t j = 0; j < N; ++j)
	{
		for (size_t i = 0; i < M; ++i)
		{
			outData[j*M + i] = static_cast<double>(Faces(i, j));
		}
	}
	pMxArray = mxCreateDoubleMatrix(M, N, mxREAL);
	memcpy(mxGetPr(pMxArray), outData, M*N*sizeof(double));
	matPutVariable(pmatfile, "Faces", pMxArray);
	delete[] outData;
	outData = NULL;

	matClose(pmatfile);
}


#if 0
void Spring::SegmentHeadAndBody()
{
	
}
#endif

const MatrixXd& Spring::GetModelVertice()
{
	return vout;
}

MatrixXd Spring::GetBodyVertices()
{
	MatrixXi& Index = SCAPE_Template.Index_BodyVertex_InSCAPE;
	MatrixXd BodyVertex(Index.rows(), 3);	
	for (int i = 0; i < BodyVertex.rows();++i)
	{
		BodyVertex.row(i) = vout.row(Index(i) - 1);
	}
	return BodyVertex;
}

MatrixXd& Spring::GetBodyColor()
{
	return SCAPE_Template.BodyColor;
}

const MatrixXi& Spring::GetBodyFaces()
{
	return SCAPE_Template.BodyFaces;
}

void Spring::SetBodyColor(MatrixXd& Colors)
{
	SCAPE_Template.BodyColor = Colors;
}