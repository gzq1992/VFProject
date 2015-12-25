#ifndef _ADJCENT_VERTICES_H
#define _ADJCENT_VERTICES_H
#include "stdafx.h"

struct AdjInfo
{
	int CenterVertexInHead;
	int CenterVertexInBody;
	int Number_Neighbors_Head;
	int* Vertex_Neighbors_Head;
	int Number_Neighbors_Body;
	int* Vertex_Neighbors_Body;
	AdjInfo() :CenterVertexInHead(0), CenterVertexInBody(0), Number_Neighbors_Body(0), Number_Neighbors_Head(0), Vertex_Neighbors_Body(NULL), Vertex_Neighbors_Head(NULL){}
};
class AdjcentVertices{
public:
	typedef AdjInfo AdjCell;
	AdjcentVertices(char* name = "MatlabFile\\AdjcentVertexM.mat");
	AdjcentVertices(const AdjcentVertices& rhs);
	~AdjcentVertices();
	AdjcentVertices& operator=(const AdjcentVertices& rhs);
	AdjCell* GetAdjcentData();
	AdjCell& operator[](int) const;
	//void AdjTest();
private:
	char* FileName;
	AdjCell* pAdjCellData;
	int NumberOfBoundaryVertices;
	void ReadData(mxArray* pData, Eigen::MatrixXi& Data);
	void DataConvert(int* DataPos,int& Count, Eigen::MatrixXi& Data, int RowNo);
	void ReadAdjcentVertices();
	void CopyAssignAndConstructorBehivar(const AdjcentVertices& rhs);
	void DestructorBehivar();
};

#endif