#ifndef _BOUNDARY_H
#define _BOUNDARY_H
#include "stdafx.h"

struct  VertexAndNeighbors{
	int IndexOfCenterVertex;
	size_t NumberOfNeighborVertex;
	int* IndexOfNeighborVertex;
	VertexAndNeighbors() :IndexOfCenterVertex(0), NumberOfNeighborVertex(0), IndexOfNeighborVertex(NULL){}
};

class Boundary
{
public:
	enum ModelPart{Head,Body};
	typedef double MatDataFormat;
	//bool ReadFlag;
	Boundary();
	Boundary(const char* filename);
	Boundary(Eigen::MatrixXi& BoundaryVertexIndex, Eigen::MatrixXi& Face,ModelPart HeadOrBody);
	Boundary(Eigen::MatrixXi& BoundaryVertexIndexInHead, Eigen::MatrixXi& FaceOfHead, Eigen::MatrixXi& BoundaryVertexIndexInBody, Eigen::MatrixXi& FaceOfBody);
	Boundary(const Boundary& rhs);
	Boundary& operator=(const Boundary& rhs);
	~Boundary();
	void SetBoundary(Eigen::MatrixXi& BoundaryVertexIndex, Eigen::MatrixXi& Face, ModelPart HeadOrBody, size_t NumsOfVertex = 0);
	void SetBoundary(Eigen::MatrixXi& BoundaryVertexIndexInHead, Eigen::MatrixXi& FaceOfHead, Eigen::MatrixXi& BoundaryVertexIndexInBody, Eigen::MatrixXi& FaceOfBody, size_t NumsOfVertex = 0);
	const int& GetBoundaryVertexIndex(int CenterVertexNo, ModelPart HeadOrBody) const;
	int* const GetNeighborVertexIndex(int CenterVertexNo, ModelPart HeadOrBody) const;
	const VertexAndNeighbors& GetVertexAndNeighbors(int CenterVertexNo, ModelPart HeadOrBody) const;
	const size_t GetNumberOfNeighbors(int CenterVertexNo,ModelPart HeadOfBody);
	const size_t GetNumberOfBoundaryVertex() const;
	Eigen::MatrixXi& GetBodyMapTable(size_t NumberOfHeadVertex,size_t NumberOfBodyVertices);
	void ReadBoundaryWeightFile(const char* filename);
	const Eigen::MatrixXd& GetBoundaryWeight() const;
	const Eigen::MatrixXi& GetBodyMapTable() const;
	const VertexAndNeighbors* GetVertexAndNeighbors(ModelPart HeadOrBody) const;
private:
	size_t NumberOfBoundaryVertex;
	VertexAndNeighbors* VertexAndNeighborsInHead;
	VertexAndNeighbors* VertexAndNeighborsInBody;
	Eigen::MatrixXi BodyMapTable;
	Eigen::MatrixXd BoundaryWeight;
	void ReadMatlabFile(const char* filename);
	void ReadmxArrayData(mxArray* pCenterVertex, mxArray* pNeighborVertex, VertexAndNeighbors* VertexAndNeighborObject);
	void SearchNeighbors(Eigen::MatrixXi& BoundaryVertexIndex, Eigen::MatrixXi& Faces, ModelPart HeadOrBody);
	void DestroyAllNeigborVertx(ModelPart HeadOrBody);
	void ConstructIndexMap(size_t NumberOfHeadVertex, size_t NumberOfBodyVertex);
	void CopyAssignmentAndCopyConstructBehivar(const Boundary& rhs);
};
#endif