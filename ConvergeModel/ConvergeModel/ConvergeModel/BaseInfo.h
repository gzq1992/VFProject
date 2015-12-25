#ifndef _BASEINFO_H
#define _BASEINFO_H
#include "stdafx.h"
#include <utility>
//#include <Eigen/Core>

class BaseInfo{
public:
	BaseInfo(char* FaceAndBoundaryFileName) :FileName(FaceAndBoundaryFileName)
	{}
	typedef size_t Rows;
	typedef size_t Cols;
	typedef std::pair<Rows, Cols> MatrixSize;
	Eigen::MatrixXi& GetHeadFace();
	Eigen::MatrixXi& GetBodyFace();
	Eigen::MatrixXi& GetHeadBoundaryVertexIndex();
	Eigen::MatrixXi& GetBodyBoundaryVertexIndex();
	BaseInfo& ReadAllVariable();
	MatrixSize GetHeadFaceSize();
	MatrixSize GetHeadBoundaryVertexIndexSize();
	MatrixSize GetBodyFaceSize();
	MatrixSize GetBodyBoundaryVertexIndexSize();
	void DataTest();
private:
	char* FileName;
	Eigen::MatrixXi HeadFace;
	Eigen::MatrixXi BodyFace;
	Eigen::MatrixXi HeadBoundaryVertexIndex;
	Eigen::MatrixXi BodyBoundaryVertexIndex;
	void ReadData(mxArray* pData, Eigen::MatrixXi& Data);
	void ReadFaceAndBoundary();
#if 0
	void ReadHeadFace(mxArray* pHeadFace);
	void ReadHeadBoundaryVertexIndex(mxArray* pHeadBoundary);
	void ReadBodyFace(mxArray* pBodyFace);
	void ReadBodyBoundaryVertexIndex(mxArray* pBodyBoundary);
#endif
};


#endif