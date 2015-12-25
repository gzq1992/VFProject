#ifndef _HUMAN_MODEL_H
#define _HUMAN_MODEL_H
#include "stdafx.h"
#include "Boundary.h"
#if 0
struct HumanHead{
	Eigen::MatrixXi HeadVertices;
	Eigen::MatrixXi HeadFaces;
	size_t NumberOfHeadVertices;
	size_t NumberOfHeadFaces;
};
struct HumanBody{
	Eigen::MatrixXi BodyVertices;
	Eigen::MatrixXi BodyFaces;
	size_t NumberOfBodyVertices;
	size_t NumberOfBodyFaces;
};
#endif

struct HumanPart{
	Eigen::MatrixXd Vertices;
	Eigen::MatrixXi Faces;
	Eigen::MatrixXd Colors;
	size_t NumberOfVertices;
	size_t NumberOfFaces;
	HumanPart() :NumberOfFaces(0), NumberOfVertices(0){}
	void SetVertices(const Eigen::MatrixXd& _Vertices)
	{
		Vertices = _Vertices;
		NumberOfVertices = _Vertices.rows();
	}

	void SetFaces(const Eigen::MatrixXi& _Faces)
	{
		Faces = _Faces;
		NumberOfFaces = _Faces.rows();
	}

	void SetColors(const Eigen::MatrixXd& _Colors)
	{
		Colors = _Colors;
	}
};


class HumanModel
{
public:
	typedef HumanPart HumanHead;
	typedef HumanPart HumanBody;
	enum PartOfHuman{Head,Body,Full};
	HumanModel();
	HumanModel(Eigen::MatrixXd& Vertices, Eigen::MatrixXi& Faces);
	HumanModel(HumanHead& HeadPart, HumanHead& BodyPart,Boundary& NowBoundary);
	~HumanModel();
	void SetHumanModel(Eigen::MatrixXd& Vertices, Eigen::MatrixXi& Faces);		
	//设置完整对象的点和面
	//void ReadSegmentBodyInfo(char* filename);
	//void SetHumanPart(const Eigen::MatrixXd& Vertices, HumanPart& part);
	//void SetHumanPart(const Eigen::MatrixXi& Faces, HumanPart& part);
	//void SetHumanPart(const Eigen::MatrixXd& Vertices, const Eigen::MatrixXi& Faces, const Eigen::MatrixXd& Colors, HumanPart& part);
	void ReadMatlabFile(const char* filename,HumanPart& part);
	//用于读取mat格式的head和body的点面信息
	void ReadObjFile(const char* filename,HumanPart& part,size_t NumberOfVertices,size_t NumberOfFaces);
	//读取obj格式的head和body
	const Eigen::MatrixXd& GetVertices() const;
	//读取当前model的点信息	HumanFullVertices
	const Eigen::MatrixXi& GetFaces() const;
	//读取当前model的面信息	HumanFullFaces
	size_t GetNumberOfVertices() const;
	//读取当前model的点数	NumberOfVertices
	size_t GetNumberOfFaces() const;
	//读取当前model的面数	NumberOfFaces
	void ConvergeHeadAndBody(HumanHead& HeadPart, HumanHead& BodyPart, Boundary& NowBoundary);
	//void TextureOnHead(Eigen::MatrixXd& Colors);
	//void TextureOnBody(Eigen::MatrixXd& Colors);
	//融合头部和身体，NowBoundary主要包含边缘点的信息
	void SaveModel2Obj(const char* filename,Eigen::MatrixXd& Vertices,Eigen::MatrixXi& Faces);
	//将模型点和面保存成OBJ格式
	//static const int NModel;
	//static Eigen::MatrixXd Time;
	//static int count;
private:
	Eigen::MatrixXd HumanFullVertices;
	//模型所有的点（融合后）
	Eigen::MatrixXi HumanFullFaces;
	Eigen::MatrixXd Colors;
	//模型所有的面（融合后）
	size_t NumberOfVertices;
	size_t NumberOfFaces;

	double AdjustScale(Eigen::MatrixXd& HeadBoundaryVertices, Eigen::MatrixXd& BodyBoundaryVertices);
	//调整头部的尺度
	Eigen::MatrixXd RigidTransformation(Eigen::MatrixXd& HeadBoundaryVertices, Eigen::MatrixXd& BodyBoundaryVertices);
	//对齐头部和身体
	void TransformVertices(Eigen::MatrixXd& Vertices,Eigen::MatrixXd R, Eigen::MatrixXd t);
	//根据R,t变换点的坐标
	Eigen::MatrixXd InterpolationForConverge(Eigen::MatrixXd CenterVertex_H,Eigen::MatrixXd& CenterVertex_B,Eigen::MatrixXd& Neighbors );
	//对边界点插值平滑
	void ConvergeHeadAndBodyStuff(HumanHead& HeadPart, HumanHead& BodyPart, Boundary& NowBoundary);
	//融合
	void PlaceHeadAndBodyToModel(HumanHead& HeadPart, HumanBody& BodyPart, Boundary& NowBoundary);
	//重排所有点和面（主要是用于合并边缘点）
	void AdjustHeadEdgeVertex(HumanHead& HeadPart,Boundary& NowBoundary);
};

#endif