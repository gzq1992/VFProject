#include"stdafx.h"
#include"BaseInfo.h"

using Eigen::MatrixX3i;
using Eigen::MatrixXi;
using std::pair;
using std::make_pair;


MatrixXi& BaseInfo::GetHeadFace()
{
	if ((HeadFace.rows() == 0) || (HeadFace.cols() == 0))
		ReadFaceAndBoundary();
	return HeadFace;
}

MatrixXi& BaseInfo::GetHeadBoundaryVertexIndex()
{
	if ((HeadBoundaryVertexIndex.rows() == 0) || (HeadBoundaryVertexIndex.cols() == 0))
		ReadFaceAndBoundary();
	return HeadBoundaryVertexIndex;
}

MatrixXi& BaseInfo::GetBodyFace()
{
	if ((BodyFace.rows() == 0) || (BodyFace.cols() == 0))
		ReadFaceAndBoundary();
	return BodyFace;
}

MatrixXi& BaseInfo::GetBodyBoundaryVertexIndex()
{
	if ((BodyBoundaryVertexIndex.cols() == 0) || (BodyBoundaryVertexIndex.rows() == 0))
		ReadFaceAndBoundary();
	return BodyBoundaryVertexIndex;
}


BaseInfo& BaseInfo::ReadAllVariable()
{
	ReadFaceAndBoundary();
	return *this;
}

BaseInfo::MatrixSize BaseInfo::GetHeadFaceSize()
{
	return make_pair(HeadFace.rows(), HeadFace.cols());
}

BaseInfo::MatrixSize BaseInfo::GetHeadBoundaryVertexIndexSize()
{
	return make_pair(HeadBoundaryVertexIndex.rows(), HeadBoundaryVertexIndex.cols());
}

BaseInfo::MatrixSize BaseInfo::GetBodyFaceSize()
{
	return make_pair(BodyFace.rows(), BodyFace.cols());
}

BaseInfo::MatrixSize BaseInfo::GetBodyBoundaryVertexIndexSize()
{
	return make_pair(BodyBoundaryVertexIndex.rows(), BodyBoundaryVertexIndex.cols());
}

void BaseInfo::DataTest()
{
	std::cout << HeadFace.cols() << '\t' << HeadFace.rows() << std::endl
		<< BodyFace.cols() << '\t' << BodyFace.rows() << std::endl
		<< HeadBoundaryVertexIndex.cols() << '\t' << HeadBoundaryVertexIndex.rows() << std::endl
		<< BodyBoundaryVertexIndex.cols() << '\t' << BodyBoundaryVertexIndex.rows() << std::endl;
	std::cout << "HeadFace:" << std::endl << std::endl << HeadFace << std::endl;
	std::cout << "HeadBoundaryVertexIndex:" << std::endl << std::endl << HeadBoundaryVertexIndex << std::endl;
	std::cout << "BodyFace:" << std::endl << std::endl << BodyFace << std::endl;
	std::cout << "BodyBoundaryVertexIndex:" << std::endl << std::endl << BodyBoundaryVertexIndex << std::endl;
}


void BaseInfo::ReadData(mxArray* pData, MatrixXi& Data)
{
	double* initA = NULL;
	initA = static_cast<double *> (mxGetData(pData));
	size_t M = mxGetM(pData);
	size_t N = mxGetN(pData);
	Data = MatrixXi(M,N);
	for (size_t i = 0; i < M; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			Data(i, j) = static_cast<int>(initA[j*M + i]);
		}
	}
}

void BaseInfo::ReadFaceAndBoundary()
{
	MATFile* pmatFile = matOpen(FileName, "r");
	if (pmatFile == NULL)
	{
		std::cerr << "FaceAndBoundary file not exist" << std::endl;
		return;
	}
	mxArray* pHeadBoundaryIndex = matGetVariable(pmatFile, "HeadBoundary");
	//ReadHeadBoundaryVertexIndex(pHeadBoundaryIndex);
	if (pHeadBoundaryIndex == NULL)
		std::cerr << "HeadBoundary doesn't esixt in the Matlab data" << std::endl;
	else
		ReadData(pHeadBoundaryIndex, HeadBoundaryVertexIndex);

	//const char* HeadFaceVarName = "HeadFace";
	mxArray* pHeadFace = matGetVariable(pmatFile, "HeadFace");
	if (pHeadFace == NULL)
		std::cerr << "HeadFace doesn't esixt in the Matlab data" << std::endl;
	else
		ReadData(pHeadFace, HeadFace);

	//const char* BodyBoundaryIndexVarName = "BodyBoundary";
	mxArray* pBodyBoundaryIndex = matGetVariable(pmatFile, "BodyBoundary");
	if (pBodyBoundaryIndex == NULL)
		std::cerr << "BodyBoundary doesn't esixt in the Matlab data" << std::endl;
	else
		ReadData(pBodyBoundaryIndex, BodyBoundaryVertexIndex);

	//const char* BodyFaceVarName = "BodyFace";
	mxArray* pBodyFace = matGetVariable(pmatFile, "BodyFace");
	//ReadBodyFace(pBodyFace);
	if (pBodyFace == NULL)
		std::cerr << "BodyFace doesn't esixt in the Matlab data" << std::endl;
	else
		ReadData(pBodyFace, BodyFace);

	matClose(pmatFile);

}


