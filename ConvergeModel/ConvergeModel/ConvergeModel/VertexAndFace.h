#ifndef _VERTEX_AND_FACE_H
#define _VERTEX_AND_FACE_H
typedef unsigned int UINT;
struct Vertex{
	double Vertex_X;
	double Vertex_Y;
	double Vertex_Z;
	Vertex(int x, int y, int z) :Vertex_X(x), Vertex_Y(y), Vertex_Z(z){}
	//Vertex(int* p) :Vertex_X(*p), Vertex_Y(*(p + 1)), Vertex_Z(*(p + 2)){}
};

struct Face{
	UINT Vertex_Index_Start;
	UINT Vertex_Index_Middle;
	UINT Vertex_Index_End;
	Face(UINT id1, UINT id2, UINT id3) :Vertex_Index_Start(id1), Vertex_Index_Middle(id2), Vertex_Index_End(id3){}
};
#endif