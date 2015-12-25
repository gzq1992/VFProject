#include "stdafx.h"
#include "AdjcentVertices.h"
#include <string.h>
using Eigen::MatrixXi;

AdjcentVertices::AdjcentVertices(char* name)
{
	FileName = new char[strlen(name) + 1];
	strcpy_s(FileName, strlen(name) + 1, name);
}

AdjcentVertices::~AdjcentVertices()
{
	DestructorBehivar();
}

void AdjcentVertices::DestructorBehivar()
{
	delete[] FileName;
	for (int idx = 0; idx < NumberOfBoundaryVertices; ++idx)
	{
		delete[](pAdjCellData[idx].Vertex_Neighbors_Head);
		pAdjCellData[idx].Vertex_Neighbors_Head = NULL;
		delete[](pAdjCellData[idx].Vertex_Neighbors_Body);
		pAdjCellData[idx].Vertex_Neighbors_Body = NULL;
	}
	delete[]pAdjCellData;
}

AdjcentVertices::AdjcentVertices(const AdjcentVertices& rhs)
{
	CopyAssignAndConstructorBehivar(rhs);
}

AdjcentVertices& AdjcentVertices::operator=(const AdjcentVertices&rhs)
{
	if (this == &rhs)
		return * this;
	DestructorBehivar();
	CopyAssignAndConstructorBehivar(rhs);
	return *this;
}

void AdjcentVertices::CopyAssignAndConstructorBehivar(const AdjcentVertices& rhs)
{
	if (this == &rhs)
		return;
	FileName = new char[strlen(rhs.FileName) + 1];
	strcpy_s(FileName, strlen(rhs.FileName) + 1, rhs.FileName);
	NumberOfBoundaryVertices = rhs.NumberOfBoundaryVertices;
	pAdjCellData = new AdjCell[NumberOfBoundaryVertices];
	for (int i = 0; i < NumberOfBoundaryVertices; ++i)
	{
		pAdjCellData[i] = rhs[i];
		pAdjCellData[i].Vertex_Neighbors_Head = new int[pAdjCellData[i].Number_Neighbors_Head];
		pAdjCellData[i].Vertex_Neighbors_Body = new int[pAdjCellData[i].Number_Neighbors_Body];

		memcpy(pAdjCellData[i].Vertex_Neighbors_Head, rhs[i].Vertex_Neighbors_Head, pAdjCellData[i].Number_Neighbors_Head * sizeof(int));
		memcpy(pAdjCellData[i].Vertex_Neighbors_Body, rhs[i].Vertex_Neighbors_Body, pAdjCellData[i].Number_Neighbors_Body * sizeof(int));

#if 0
		for (int j = 0; j < pAdjCellData[i].Number_Neighbors_Head; ++j)
		{
			pAdjCellData[i].Vertex_Neighbors_Head[j] = rhs[i].Vertex_Neighbors_Head[j];
		}
		for (int j = 0; j < pAdjCellData[i].Number_Neighbors_Body; ++j)
		{
			pAdjCellData[i].Vertex_Neighbors_Body[j] = rhs[i].Vertex_Neighbors_Body[j];
		}
#endif
	}
}




AdjcentVertices::AdjCell* AdjcentVertices::GetAdjcentData()
{
	if (FileName == NULL)
	{
		std::cerr << "Filename may be wrong" << std::endl;
		return NULL;
	}
	if (pAdjCellData == NULL)
		ReadAdjcentVertices();
	if (pAdjCellData == NULL)
	{
		std::cerr << "Some errors occured in data read and convert" << std::endl;
		return NULL;
	}
	return pAdjCellData;
}

AdjcentVertices::AdjCell& AdjcentVertices::operator [](int x) const
{
	//if (pAdjCellData == NULL)
	//	ReadAdjcentVertices();
	return *(pAdjCellData+x);
}

#if 0
void AdjcentVertices::AdjTest()
{
	if (FileName == NULL)
	{
		std::cerr << "Filename may be wrong" << std::endl;
		return;
	}
	if (pAdjCellData == NULL)
		ReadAdjcentVertices();
	if (pAdjCellData == NULL)
	{
		std::cerr << "Some errors occured in data read and convert" << std::endl;
		return;
	}
	std::cout << "NumberOfBoundaryVertices:" << NumberOfBoundaryVertices << std::endl;
	std::cout << "CenterVertexInHead and CenterVertexInBody" << std::endl;
	AdjCell* p = pAdjCellData;
	for (int i = 0; i < NumberOfBoundaryVertices;++i)
	{
		p = pAdjCellData + i;
		std::cout << i << ' ' << p->CenterVertexInHead << ' ' << p->CenterVertexInBody << ' ' << p->Number_Neighbors_Head << ' ' << p->CenterVertexInBody << std::endl;
	}
}
#endif

void AdjcentVertices::ReadData(mxArray* pData, Eigen::MatrixXi& Data)
{
	double* initA = NULL;
	initA = static_cast<double *> (mxGetData(pData));
	size_t M = mxGetM(pData);
	size_t N = mxGetN(pData);
	Data = MatrixXi(M, N);
	for (size_t i = 0; i < M; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			Data(i, j) = static_cast<int>(initA[j*M + i]);
		}
	}
}

void AdjcentVertices::DataConvert(int* DataPos, int& Count, Eigen::MatrixXi& Data, int RowNo)
{
	Count = Data(RowNo,Data.cols() -1);
	DataPos = new int[Count];	
	for (int idx = 0; idx < Count; ++idx)
		DataPos[idx] = Data(RowNo, idx);
}

void AdjcentVertices::ReadAdjcentVertices()
{
	MATFile* pmatFile = matOpen(FileName, "r");
	if (pmatFile == NULL)
	{
		std::cerr << "FaceAndBoundary file not exist" << std::endl;
		return;
	}

	MatrixXi CenterVertexInHead;
	mxArray* pCenterVertexInHead = matGetVariable(pmatFile, "CenterVertexInHead");
	if (pCenterVertexInHead == NULL)
		std::cerr << "CenterVertexInHead doesn't esixt in the Matlab data" << std::endl;
	else
		ReadData(pCenterVertexInHead, CenterVertexInHead);
	mxFree(pCenterVertexInHead);

	MatrixXi Neighbor_Vertex_Head;
	mxArray* pNeighbor_Vertex_Head = matGetVariable(pmatFile, "Neighbor_Vertex_Head");
	if (pNeighbor_Vertex_Head == NULL)
		std::cerr << "Neighbor_Vertex_Head doesn't esixt in the Matlab data" << std::endl;
	else
		ReadData(pNeighbor_Vertex_Head, Neighbor_Vertex_Head);
	mxFree(pNeighbor_Vertex_Head);

	MatrixXi CenterVertexInBody;
	mxArray* pCenterVertexInBody = matGetVariable(pmatFile,"CenterVertexInBody");
	if (pCenterVertexInBody == NULL)
		std::cerr << "CenterVertexInBody doesn't esixt in the Matlab data" << std::endl;
	else
		ReadData(pCenterVertexInBody, CenterVertexInBody);
	mxFree(pCenterVertexInBody);

	MatrixXi Neighbor_Vertex_Body;
	mxArray* pNeighbor_Vertex_Body = matGetVariable(pmatFile, "Neighbor_Vertex_Body");
	if (pNeighbor_Vertex_Body == NULL)
		std::cerr << "Neighbor_Vertex_Body doesn't esixt in the Matlab data" << std::endl;
	else
		ReadData(pNeighbor_Vertex_Body, Neighbor_Vertex_Body);
	mxFree(pNeighbor_Vertex_Body);

	matClose(pmatFile);
	NumberOfBoundaryVertices = static_cast<int>(CenterVertexInHead.rows());
	pAdjCellData = new AdjInfo[NumberOfBoundaryVertices];
	for (int i = 0; i < NumberOfBoundaryVertices;++i)
	{
		pAdjCellData[i].CenterVertexInHead = CenterVertexInHead(i);
		pAdjCellData[i].CenterVertexInBody = CenterVertexInBody(i);
		DataConvert(pAdjCellData[i].Vertex_Neighbors_Head, pAdjCellData[i].Number_Neighbors_Head, Neighbor_Vertex_Head, i);
		DataConvert(pAdjCellData[i].Vertex_Neighbors_Body, pAdjCellData[i].Number_Neighbors_Body, Neighbor_Vertex_Body, i);
	}
}