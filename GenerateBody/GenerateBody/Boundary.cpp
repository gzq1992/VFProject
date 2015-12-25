#include "stdafx.h"
#include "Boundary.h"
#include <set>
using Eigen::MatrixXi;
using Eigen::MatrixXd;
using std::set;

Boundary::Boundary() : NumberOfBoundaryVertex(0), VertexAndNeighborsInHead(NULL), VertexAndNeighborsInBody(NULL)
{}

Boundary::Boundary(char* filename) : NumberOfBoundaryVertex(0), VertexAndNeighborsInHead(NULL), VertexAndNeighborsInBody(NULL)
{
	ReadMatlabFile(filename);
}

Boundary::Boundary(Eigen::MatrixXi& BoundaryVertexIndex, Eigen::MatrixXi& Face,ModelPart HeadOrBody)
{
	SearchNeighbors(BoundaryVertexIndex, Face, HeadOrBody);
}

Boundary::Boundary(Eigen::MatrixXi& BoundaryVertexIndexInHead, Eigen::MatrixXi& FaceOfHead, Eigen::MatrixXi& BoundaryVertexIndexInBody, Eigen::MatrixXi& FaceOfBody)
{
	SetBoundary(BoundaryVertexIndexInHead, FaceOfHead, BoundaryVertexIndexInBody, FaceOfBody);
}

Boundary::Boundary(const Boundary& rhs)
{
	CopyAssignmentAndCopyConstructBehivar(rhs);
}

Boundary& Boundary::operator=(const Boundary& rhs)
{
	DestroyAllNeigborVertx(Head);
	DestroyAllNeigborVertx(Body);
	CopyAssignmentAndCopyConstructBehivar(rhs);
	return *this;
}



void Boundary::CopyAssignmentAndCopyConstructBehivar(const Boundary& rhs)
{
	if (this == &rhs)
		return;

	NumberOfBoundaryVertex = rhs.GetNumberOfBoundaryVertex();
	BodyMapTable = rhs.GetBodyMapTable();
	BoundaryWeight = rhs.GetBoundaryWeight();
	VertexAndNeighborsInHead = new VertexAndNeighbors[NumberOfBoundaryVertex];
	VertexAndNeighborsInBody = new VertexAndNeighbors[NumberOfBoundaryVertex];
	for (int k = 0; k < 2;++k)
	{
		const VertexAndNeighbors* ptr_rhs_VertexAndNeighbors = NULL;
		VertexAndNeighbors* ptr_this_VertexAndNeighbors = NULL;
		if (k == 0)
		{
			ptr_rhs_VertexAndNeighbors = rhs.GetVertexAndNeighbors(Head);
			ptr_this_VertexAndNeighbors = VertexAndNeighborsInHead;
		}
		else
		{
			ptr_rhs_VertexAndNeighbors = rhs.GetVertexAndNeighbors(Body);
			ptr_this_VertexAndNeighbors = VertexAndNeighborsInBody;
		}

		memcpy(ptr_this_VertexAndNeighbors, ptr_rhs_VertexAndNeighbors, NumberOfBoundaryVertex*sizeof(VertexAndNeighbors));
		for (size_t i = 0; i < NumberOfBoundaryVertex; ++i)
		{
			ptr_this_VertexAndNeighbors[i].IndexOfNeighborVertex = new int[ptr_this_VertexAndNeighbors[i].NumberOfNeighborVertex];
			memcpy(ptr_this_VertexAndNeighbors[i].IndexOfNeighborVertex, ptr_rhs_VertexAndNeighbors[i].IndexOfNeighborVertex, ptr_this_VertexAndNeighbors[i].NumberOfNeighborVertex*sizeof(int));
		}
#if 0
		for (size_t i = 0; i < NumberOfBoundaryVertex; ++i)
		 {
			ptr_this_VertexAndNeighbors[i] = ptr_rhs_VertexAndNeighbors[i];
			ptr_this_VertexAndNeighbors[i].IndexOfNeighborVertex = new int[ptr_this_VertexAndNeighbors[i].NumberOfNeighborVertex];
			for (size_t j = 0; j < ptr_this_VertexAndNeighbors[i].NumberOfNeighborVertex; ++j)
			{
				ptr_this_VertexAndNeighbors[i].IndexOfNeighborVertex[j] = ptr_rhs_VertexAndNeighbors[i].IndexOfNeighborVertex[j];
			}
		 }
#endif
	}	
}

Boundary::~Boundary()
{
#if 0
	for (int i = 0; i < NumberOfBoundaryVertex; ++i)
	{
		delete[](VertexAndNeighborsInHead[i].IndexOfNeighborVertex);
		delete[](VertexAndNeighborsInBody[i].IndexOfNeighborVertex);
	}
#endif
	DestroyAllNeigborVertx(Head);
	DestroyAllNeigborVertx(Body);

	delete[]VertexAndNeighborsInHead;
	delete[]VertexAndNeighborsInBody;
}

void Boundary::SetBoundary(Eigen::MatrixXi& BoundaryVertexIndex, Eigen::MatrixXi& Face, ModelPart HeadOrBody, size_t NumsOfVertex)
{
	if (NumsOfVertex == 0)
		NumsOfVertex = BoundaryVertexIndex.rows();
	DestroyAllNeigborVertx(HeadOrBody);
	if (Head == HeadOrBody)
	{
		delete[]VertexAndNeighborsInHead;
		VertexAndNeighborsInHead = new VertexAndNeighbors[NumsOfVertex];
	}
	if (Body == HeadOrBody)
	{
		delete[]VertexAndNeighborsInBody;
		VertexAndNeighborsInBody = new VertexAndNeighbors[NumsOfVertex];
	}
	NumberOfBoundaryVertex = NumsOfVertex;
	SearchNeighbors(BoundaryVertexIndex, Face, HeadOrBody);
}

void Boundary::SetBoundary(Eigen::MatrixXi& BoundaryVertexIndexInHead, Eigen::MatrixXi& FaceOfHead, Eigen::MatrixXi& BoundaryVertexIndexInBody, Eigen::MatrixXi& FaceOfBody, size_t NumsOfVertex)
{
	SetBoundary(BoundaryVertexIndexInHead, FaceOfHead, Head, NumsOfVertex);
	SetBoundary(BoundaryVertexIndexInBody, FaceOfBody, Body, NumsOfVertex);
}

const int& Boundary::GetBoundaryVertexIndex(int CenterVertexNo, ModelPart HeadOrBody) const
{
	if (Head == HeadOrBody)
		return VertexAndNeighborsInHead[CenterVertexNo].IndexOfCenterVertex;
	else
		return VertexAndNeighborsInBody[CenterVertexNo].IndexOfCenterVertex;
}

int* const Boundary::GetNeighborVertexIndex(int CenterVertexNo, ModelPart HeadOrBody) const
{
	if (Head == HeadOrBody)
		return VertexAndNeighborsInHead[CenterVertexNo].IndexOfNeighborVertex;
	else
		return VertexAndNeighborsInBody[CenterVertexNo].IndexOfNeighborVertex;
}

const VertexAndNeighbors& Boundary::GetVertexAndNeighbors(int CenterVertexNo, ModelPart HeadOrBody) const
{
	if (Head == HeadOrBody)
		return VertexAndNeighborsInHead[CenterVertexNo];
	else
		return VertexAndNeighborsInBody[CenterVertexNo];
}

const size_t Boundary::GetNumberOfNeighbors(int CenterVertexNo, ModelPart HeadOrBody)
{
	if (Head == HeadOrBody)
		return VertexAndNeighborsInHead[CenterVertexNo].NumberOfNeighborVertex;
	else
		return VertexAndNeighborsInBody[CenterVertexNo].NumberOfNeighborVertex;
}

const size_t Boundary::GetNumberOfBoundaryVertex() const
{
	return NumberOfBoundaryVertex;
}

const MatrixXi& Boundary::GetBodyMapTable() const
{
	return BodyMapTable;
}

const VertexAndNeighbors* Boundary::GetVertexAndNeighbors(ModelPart HeadOrBody) const
{
	if (HeadOrBody == Head)
		return VertexAndNeighborsInHead;
	else
		return VertexAndNeighborsInBody;
}




void Boundary::ReadMatlabFile(char* filename)
{
	MATFile* pmatFile = matOpen(filename, "r");
	if (pmatFile == NULL)
	{
		std::cerr << "Not find Boundary file" << std::endl;
		return;
	}

	mxArray* pNumberOfBoundaryVertex = matGetVariable(pmatFile, "NumberOfBoundaryVertex");
	MatDataFormat* pdNumberOfBoundaryVertex = static_cast<MatDataFormat*>(mxGetData(pNumberOfBoundaryVertex));
	NumberOfBoundaryVertex = static_cast<int> (*pdNumberOfBoundaryVertex);
	VertexAndNeighborsInHead = new VertexAndNeighbors[NumberOfBoundaryVertex];
	VertexAndNeighborsInBody = new VertexAndNeighbors[NumberOfBoundaryVertex];

	mxFree(pdNumberOfBoundaryVertex);

	mxArray* pCenterVertexInHead = matGetVariable(pmatFile, "CenterVertexInHead");
#if 0
	MatDataFormat* pdCenterVertexInHead = static_cast<MatDataFormat*> (mxGetData(pCenterVertexInHead));
	int M = mxGetM(pCenterVertexInHead);
	int N = mxGetN(pCenterVertexInHead);
	for (int idx = 0; idx < NumberOfBoundaryVertex; ++idx)
		VertexAndNeighborsInHead[idx].IndexOfCenterVertex = pdCenterVertexInHead[idx];
#endif
	mxArray* pCenterVertexInBody = matGetVariable(pmatFile, "CenterVertexInBody");
	mxArray* pNeighborVertexInHead = matGetVariable(pmatFile, "NeighborVertexInHead");
	mxArray* pNeighborVertexInBody = matGetVariable(pmatFile, "NeighborVertexInBody");
	ReadmxArrayData(pCenterVertexInHead, pNeighborVertexInHead, VertexAndNeighborsInHead);
	ReadmxArrayData(pCenterVertexInBody, pNeighborVertexInBody, VertexAndNeighborsInBody);

	matClose(pmatFile);
}

void Boundary::ReadmxArrayData(mxArray* pCenterVertex,mxArray* pNeighborVertex,VertexAndNeighbors* VertexAndNeighborObject)//for read int data,not suitable for double data
{
	MatDataFormat* pdCenterVertex = static_cast<MatDataFormat*> (mxGetData(pCenterVertex));
	MatDataFormat* pdNeighborVertex = static_cast<MatDataFormat*> (mxGetData(pNeighborVertex));
	size_t MCenter = mxGetM(pCenterVertex);
	size_t NCenter = mxGetN(pCenterVertex);
	for (size_t i = 0; i < MCenter;++i)
	{
		VertexAndNeighborObject[i].IndexOfCenterVertex = static_cast<int> (pdCenterVertex[i]);
	}
	size_t MNeighbor = mxGetM(pNeighborVertex);
	size_t NNeighbor = mxGetN(pNeighborVertex);
	size_t Offset = (NNeighbor - 1) * MNeighbor;
	for (size_t i = 0; i < MCenter;++i)
	{
		VertexAndNeighborObject[i].NumberOfNeighborVertex = static_cast<int> (pdNeighborVertex[Offset + i]);
		VertexAndNeighborObject[i].IndexOfNeighborVertex = new int[VertexAndNeighborObject[i].NumberOfNeighborVertex];
		for (int j = 0; j < VertexAndNeighborObject[i].NumberOfNeighborVertex;++j)
		{
			VertexAndNeighborObject[i].IndexOfNeighborVertex[j] = static_cast<int> (pdNeighborVertex[j*MNeighbor+i]);
		}
	}
	mxFree(pdCenterVertex);
	mxFree(pdNeighborVertex);
}

void Boundary::SearchNeighbors(Eigen::MatrixXi& BoundaryVertexIndex, Eigen::MatrixXi& Faces, ModelPart HeadOrBody)
{
	size_t RowsOfFaces = Faces.rows();
	size_t ColsOfFaces = Faces.cols();
	VertexAndNeighbors* ptr = NULL;
	if (Head == HeadOrBody)
	{
		ptr = VertexAndNeighborsInHead;
	}
	if (Body == HeadOrBody)
	{
		ptr = VertexAndNeighborsInBody;
	}
	for (size_t idx = 0; idx < NumberOfBoundaryVertex; ++idx)
	{
		set<int> Neighbors;
		for (size_t idy = 0; idy < RowsOfFaces; ++idy)
		{
			if (BoundaryVertexIndex(idx) == Faces(idy,0))
			{
				Neighbors.insert(Faces(idy, 1));
				Neighbors.insert(Faces(idy, 2));
			}
			else
			{
				if (BoundaryVertexIndex(idx) == Faces(idy,1))
				{
					Neighbors.insert(Faces(idy, 0));
					Neighbors.insert(Faces(idy, 2));
				}
				else
				{
					if (BoundaryVertexIndex(idx) == Faces(idy,2))
					{
						Neighbors.insert(Faces(idy, 0));
						Neighbors.insert(Faces(idy, 1));
					}
				}
			}
		}
		ptr[idx].IndexOfCenterVertex = BoundaryVertexIndex(idx);
		ptr[idx].NumberOfNeighborVertex = Neighbors.size();
		ptr[idx].IndexOfNeighborVertex = new int[Neighbors.size()];
		int* p = ptr[idx].IndexOfNeighborVertex;
		set<int>::iterator iter = Neighbors.cbegin();
		for (size_t i = 0; i < Neighbors.size(); ++i)
		{
			p[i] = *iter;
			++iter;
		}
	}
}

void Boundary::DestroyAllNeigborVertx(ModelPart HeadOrBody)
{
	VertexAndNeighbors* VertexAndNeighborObject = NULL;
	if (Head == HeadOrBody)
		VertexAndNeighborObject = VertexAndNeighborsInHead;
	if (Body == HeadOrBody)
		VertexAndNeighborObject = VertexAndNeighborsInBody;
	for (int i = 0; i < NumberOfBoundaryVertex; ++i)
	{
		VertexAndNeighborObject[i].NumberOfNeighborVertex = 0;
		delete[](VertexAndNeighborObject[i].IndexOfNeighborVertex);
		VertexAndNeighborObject[i].IndexOfNeighborVertex= NULL;
	}
}



void Boundary::ConstructIndexMap(size_t NumberOfHeadVertex, size_t NumberOfBodyVertex)
{
	MatrixXi BoundaryMap(NumberOfBoundaryVertex, 2);
	for (int i = 0; i < NumberOfBoundaryVertex; ++i)
	{
		BoundaryMap(i, 0) = VertexAndNeighborsInHead[i].IndexOfCenterVertex;
		BoundaryMap(i, 1) = VertexAndNeighborsInBody[i].IndexOfCenterVertex;
	}

	BodyMapTable.resize(NumberOfBodyVertex, 1);
	int loc = 1;
	bool flag = false;
	int count = 0;

	for (int i = 1; i <= NumberOfBodyVertex;++i)
	{		
		flag = false;
		if ((i >= 10799) && (i <= 11072))
		{
			count = 0;
			while (count < NumberOfBoundaryVertex)
			{
				if (i == BoundaryMap(count, 1))
				{
					flag = true;
					break;
				}
				++count;
			}
		}
		if (flag == false)
		{
			BodyMapTable(i - 1) = static_cast<int>(loc + NumberOfHeadVertex);
			++loc;
		}
		else
		{
			BodyMapTable(i - 1) = BoundaryMap(count, 0);
		}		
	}

}





MatrixXi& Boundary::GetBodyMapTable(size_t NumberOfHeadVertex, size_t NumberOfBodyVertices)
{
	ConstructIndexMap(NumberOfHeadVertex,NumberOfBodyVertices);
	return BodyMapTable;
}

void Boundary::ReadBoundaryWeightFile(char* filename)
{
	MATFile* pmatFile = matOpen(filename, "r");
	if (pmatFile == NULL)
	{
		std::cerr << "Not find BoundaryWeight file" << std::endl;
		return;
	}
	mxArray* pBoundaryWeight = matGetVariable(pmatFile, "BoundaryWeight");
	MatDataFormat* pdBoundaryWeight = static_cast<MatDataFormat*>(mxGetData(pBoundaryWeight));
	size_t M = mxGetM(pBoundaryWeight);
	size_t N = mxGetN(pBoundaryWeight);
	BoundaryWeight.resize(M, N);
	//for (size_t i = 0; i < M;++i)
	//{
	//	for (size_t j = 0; j < N;++j)
	//	{
	//		BoundaryWeight(i, j) = pdBoundaryWeight[j*M + i];
	//	}
	//}
	memcpy(BoundaryWeight.data(), pdBoundaryWeight, sizeof(double)* M * N);
	mxFree(pdBoundaryWeight);
}

const MatrixXd& Boundary::GetBoundaryWeight() const
{
	return BoundaryWeight;
}
