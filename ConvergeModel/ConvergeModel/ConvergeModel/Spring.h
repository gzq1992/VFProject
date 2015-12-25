#ifndef _SPRING_H
#define _SPRING_H
#include "stdafx.h"
#include <vector>
#include <Eigen/SparseCore>

struct SCAPE 
{
	Eigen::MatrixXd* Paramters_2_Shape;
	Eigen::MatrixXi Face_In_Bone;
	Eigen::MatrixXd Template;
	Eigen::MatrixXd* SCAPE_P;
	Eigen::SparseMatrix<double> SCAPE_M;
	Eigen::SparseMatrix<double> SCAPE_N;
	Eigen::MatrixXi Faces;
	Eigen::MatrixXi Bone_Of_Face;
	Eigen::MatrixXi Bone_Of_Vertex;
	Eigen::MatrixXi Bone_Conjunctions;
	Eigen::MatrixXd* T;
	Eigen::MatrixXd Joints;
	Eigen::MatrixXd Body;
	Eigen::MatrixXi Index_BodyVertex_InSCAPE;
	Eigen::MatrixXi BodyFaces;
	Eigen::MatrixXd BodyColor;
	static const size_t Number_Of_Faces = 25000;
	static const size_t Number_Of_Bones = 16;
	static const size_t Number_Of_Vertices = 12500;
	static const size_t Paramters_2_Shape_sizex = 9;
	static const size_t Paramters_2_Shape_sizey = 10;
	static const size_t SCAPE_P_sizex = 9;
	static const size_t SCAPE_P_sizey = 7;
	static const size_t SCAPE_M_rows = 12500;
	static const size_t SCAPE_M_cols = 12500;
	static const size_t SCAPE_N_rows = 12500;
	static const size_t SCAPE_N_cols = 50001;
};

//struct SegmentRef
//{
//	Eigen::MatrixXi Index_BodyVertex_InSCAPE;
//	Eigen::MatrixXi Index_HeadVertex_InSCAPE;
//	Eigen::MatrixXi Index_BoundaryVertex_InBody;
//	Eigen::MatrixXi Index_BoundaryVertex_InHead;
//};

//struct SegmentBodyData
//{
//	Eigen::MatrixXd BodyVertex;
//	Eigen::MatrixXd HeadVertex;
//	Eigen::MatrixXd BodyBoundaryVertex;
//};

class Spring{
public:
	typedef double MatDataFormat;
	Spring(char* SCAPE_Parameters_Filename);
	//Spring(const SCAPE& rhs);
	~Spring();
	Spring& GenerateModel(double parameters[],size_t NumberOfParameters,Eigen::MatrixXd T, Eigen::MatrixXd& theta);
	Spring& GenerateModel(double parameters[], size_t NumberOfParameters, Eigen::MatrixXd T);
	const Eigen::MatrixXd& GetModelVertice();
	Eigen::MatrixXd GetBodyVertices();
	Eigen::MatrixXd& GetBodyColor();
	const Eigen::MatrixXi& GetBodyFaces();
	void SetBodyColor(Eigen::MatrixXd& Colors);
	//void ReadSegmentRef(char* RefFilename);
	void SaveModel2Obj(const char* filename);
	void SaveModel2Obj(const char* filename,const Eigen::MatrixXd& Vertices, const Eigen::MatrixXi& Faces);
	void WriteMatFile(const char* filename, Eigen::MatrixXd& Data,char* varName);
	void WriteMatFile(const char* filename, Eigen::MatrixXd& Vertices, Eigen::MatrixXi& Faces);
private:
	SCAPE SCAPE_Template;
	//SegmentRef SegmentIndex;
	Eigen::MatrixXd vout;
	//Eigen::MatrixXd normal;
	Spring(const SCAPE& rhs);
	Spring& operator=(const SCAPE& rhs);
	void Init(char* SCAPE_Parameters_Filename);
	void Read2DMatData(mxArray*, Eigen::MatrixXd&);
	void Read2DMatData(mxArray*, Eigen::MatrixXi&);
	void Read3DMatData(mxArray*, Eigen::MatrixXd* &,size_t,size_t,size_t);
	void ReadSparseData(mxArray*, Eigen::SparseMatrix<double>&,size_t,size_t);
	void GenerateSpringModel(double parameters[], size_t NumberOfParameters, Eigen::MatrixXd T, Eigen::MatrixXd& theta);
	Eigen::MatrixXd SCAPE_vector_from_motion(Eigen::MatrixXd& Vertices, Eigen::MatrixXi& Faces, Eigen::MatrixXi& Bone_Of_Face, Eigen::MatrixXi& Bone_Conjunctions, Eigen::MatrixXd& Joints, size_t face, Eigen::MatrixXd& Theta);
	//void SegmentHeadAndBody();

	//void ComputeNormal(Eigen::MatrixXd&, Eigen::MatrixXi&);
	
};


#endif