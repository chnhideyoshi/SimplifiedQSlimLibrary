#ifndef BASE_H
#define BASE_H
#define  _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <algorithm>
struct Point3d
{
public:
	float X;
	float Y;
	float Z;
	Point3d() 
	{
		X = 0;
		Y = 0;
		Z = 0;
	}
	~Point3d()
	{
	}
	Point3d(float x, float y, float z) 
	{
		this->X = x;
		this->Y = y;
		this->Z = z;
	}
};
struct Triangle
{
public :
	int P0Index;
	int P1Index;
	int P2Index;
	Triangle(int p0index, int p1index, int p2index)
	{
		this->P0Index=p0index;
		this->P1Index=p1index;
		this->P2Index=p2index;
	}
	Triangle()
	{
		P0Index=-1;
		P1Index=-1;
		P2Index=-1;
	}
};
struct Vector
{
public:
	float X;
	float Y;
	float Z;
	Vector() 
	{
		X = 0;
		Y = 0;
		Z = 0;
	}
	~Vector()
	{
	}
	Vector(float x, float y, float z) 
	{
		this->X = x;
		this->Y = y;
		this->Z = z;
	}
};
Vector CaculateTriangleNormal(Point3d& p0, Point3d& p1, Point3d& p2)
{
	Vector Normal;
	float v1x = p1.X - p0.X;
	float v1y = p1.Y - p0.Y;
	float v1z = p1.Z - p0.Z;
	float v2x = p2.X - p1.X;
	float v2y = p2.Y - p1.Y;
	float v2z = p2.Z - p1.Z;
	Normal.X= v1y * v2z - v1z * v2y;
	Normal.Y = v1z * v2x - v1x * v2z;
	Normal.Z = v1x * v2y - v1y * v2x;
	float len = (float)sqrt(Normal.X * Normal.X + Normal.Y * Normal.Y + Normal.Z * Normal.Z);
	if (len == 0)
	{
	}
	else
	{
		Normal.X /= len;
		Normal.Y /= len;
		Normal.Z /= len;
	}
	return Normal;
}
class Mesh
{
public:
	std::vector<Point3d> Vertices;
	std::vector<Triangle> Faces;
	std::vector<Vector> FaceNormals;
	std::vector<Vector> VertexNormals;
	std::vector<std::vector<int>*> AdjacentFacesPerVertex;
	std::vector<std::vector<int>*> AdjacentVerticesPerVertex;
	Mesh()
	{
	}
	~Mesh()
	{
	}
	int AddVertex(Point3d& toAdd)
	{
		int index = Vertices.size();
		Vertices.push_back(toAdd);
		return index;
	}
	int AddFace(Triangle& tri)
	{
		int index = Faces.size();
		Faces.push_back(tri);
		return index;
	}
	void CaculateFaceNormals()
	{
		FaceNormals.reserve(Faces.size());
		for(size_t i=0;i<Faces.size();i++)
		{
			Point3d& p0=Vertices[Faces[i].P0Index];
			Point3d& p1=Vertices[Faces[i].P1Index];
			Point3d& p2=Vertices[Faces[i].P2Index];
			Vector v=CaculateTriangleNormal(p0,p1,p2);
			FaceNormals.push_back(v);
		}
	}
	void CaculateVertexNormals()
	{
		if(FaceNormals.size()==0)
			CaculateFaceNormals();
		if(AdjacentFacesPerVertex.size()==0)
			CaculateAdjacentFacesPerVertex();
		VertexNormals.reserve(Vertices.size());
		for (size_t i = 0; i < Vertices.size(); i++)
		{
			float sumx = 0;
			float sumy = 0;
			float sumz = 0;
			std::vector<int>& tlist = *(AdjacentFacesPerVertex[i]);
			if(tlist.size()!=0)
			{
				for (size_t j = 0; i < tlist.size(); j++)
				{
					sumx += FaceNormals[tlist[j]].X;
					sumy += FaceNormals[tlist[j]].Y;
					sumz += FaceNormals[tlist[j]].Z;
				}
				VertexNormals.push_back(Vector(sumx / tlist.size(), sumy / tlist.size(), sumz /tlist.size()));
			}
			else
			{
				VertexNormals.push_back(Vector(0,0,0));
			}
		}
	}
	void CaculateAdjacentFacesPerVertex()
	{
		AdjacentFacesPerVertex.reserve(Vertices.size());
		for (size_t i = 0; i < Vertices.size(); i++)
		{
			std::vector<int>* list=new std::vector<int>();
			list->reserve(4);
			AdjacentFacesPerVertex.push_back(list);
		}
		for (size_t i = 0; i < Faces.size(); i++)
		{
			Triangle& t = Faces[i];
			std::vector<int> *t0list= AdjacentFacesPerVertex[t.P0Index];
			std::vector<int> *t1list= AdjacentFacesPerVertex[t.P1Index];
			std::vector<int> *t2list= AdjacentFacesPerVertex[t.P2Index];
			t0list->push_back(i);
			t1list->push_back(i);
			t2list->push_back(i);
		}
	}
	void CaculateAdjacentVerticesPerVertex()
	{
		AdjacentVerticesPerVertex.reserve(Vertices.size());
		for (size_t i = 0; i < Vertices.size(); i++)
		{
			std::vector<int>* list=new std::vector<int>();
			list->reserve(4);
			AdjacentVerticesPerVertex.push_back(list);
		}
		for (size_t i = 0; i < Faces.size(); i++)
		{
			Triangle &t = Faces[i];
			std::vector<int> *p0list= AdjacentVerticesPerVertex[t.P0Index];
			std::vector<int> *p1list= AdjacentVerticesPerVertex[t.P1Index];
			std::vector<int> *p2list= AdjacentVerticesPerVertex[t.P2Index];
			if (std::find(p0list->begin(), p0list->end(), t.P1Index)==p0list->end())
				p0list->push_back(t.P1Index);
			if (std::find(p0list->begin(), p0list->end(), t.P2Index)==p0list->end())
				p0list->push_back(t.P2Index);
			if (std::find(p1list->begin(), p1list->end(), t.P0Index)==p1list->end())
				p1list->push_back(t.P0Index);
			if (std::find(p1list->begin(), p1list->end(), t.P2Index)==p1list->end())
				p1list->push_back(t.P2Index);
			if (std::find(p2list->begin(), p2list->end(), t.P0Index)==p2list->end())
				p2list->push_back(t.P0Index);
			if (std::find(p2list->begin(), p2list->end(), t.P1Index)==p2list->end())
				p2list->push_back(t.P1Index);
		}
	}
	void LaplacianSmooth(int time)
	{
		if(AdjacentVerticesPerVertex.size()==0)
			CaculateAdjacentVerticesPerVertex();
		Point3d* tempPos=new Point3d[Vertices.size()];
		for(int k=0;k<time;k++)
		{
			for (size_t i = 0; i < Vertices.size(); i++)
			{
				float xav = 0;
				float yav = 0;
				float zav = 0;
				std::vector<int>& adjlist =*(AdjacentVerticesPerVertex[i]);
				size_t adjcount=adjlist.size();
				if(adjcount==0)
					continue;
				for (size_t j = 0; j < adjcount; j++)
				{
					xav += Vertices[adjlist[j]].X;
					yav += Vertices[adjlist[j]].Y;
					zav += Vertices[adjlist[j]].Z;
				}
				xav /= adjcount;
				yav /= adjcount;
				zav /= adjcount;
				tempPos[i].X = xav;
				tempPos[i].Y = yav;
				tempPos[i].Z = zav;
			}
			for (size_t i = 0; i < Vertices.size(); i++)
			{
				Vertices[i].X = tempPos[i].X;
				Vertices[i].Y = tempPos[i].Y;
				Vertices[i].Z = tempPos[i].Z;
			}
		}
		delete[] tempPos;
	}
};
class PlyManager
{
public:
	static void Output(Mesh& mesh,const char* filename);
	static void ReadFile(Mesh& mesh,const char* fileName);
private:
	static void AWriteV(FILE* sw, double v1, double v2, double v3,unsigned char r,unsigned char g,unsigned char b);
	static void AWriteF(FILE* sw, long i1, long i2, long i3);
};
void PlyManager::AWriteV(FILE* file, double v1, double v2, double v3,unsigned char r,unsigned char g,unsigned char b)
{
	fprintf(file,"%.2f %.2f %.2f %d %d %d\n",v1,v2,v3,r,g,b);
}
void PlyManager::AWriteF(FILE* file, long i1, long i2, long i3)
{
	fprintf(file,"%d %d %d %d\n",3,i1,i2,i3);
}
void PlyManager::ReadFile(Mesh& mesh,const char* fileName)
{
	long vcount=0;
	long fcount=0;
	FILE * nfile = fopen(fileName,"r");
	fscanf(nfile,"ply\nformat ascii 1.0\ncomment VCGLIB generated\nelement vertex %d\n",&vcount);
	fscanf(nfile,"property float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement face %d\n",&fcount);
	fscanf(nfile,"property list int int vertex_indices\nend_header\n");
	float v1=0,v2=0,v3=0;
	int r=0,g=0,b=0;
	long i1=0,i2=0,i3=0;
	for(long i=0;i<vcount;i++)
	{
		fscanf(nfile,"%f %f %f %d %d %d\n",&v1,&v2,&v3,&r,&g,&b);
		//qym 
		//mesh.AddVertex(Point3d(v1,v2,v3));
		Point3d p3d(v1,v2,v3);
		mesh.AddVertex(p3d);
	}
	for(long j=0;j<fcount;j++)
	{
		fscanf(nfile,"3 %d %d %d\n",&i1,&i2,&i3);
		Triangle t(i1,i2,i3);
		mesh.AddFace(t);
	}
	fclose(nfile);
}
void PlyManager::Output(Mesh& mesh,const char* filename)
{
	FILE * nfile = fopen(filename,"wb");
	fprintf(nfile,"ply\n");
	fprintf(nfile,"format ascii 1.0\n");
	fprintf(nfile,"comment VCGLIB generated\n");
	fprintf(nfile,"element vertex %d\n",mesh.Vertices.size());
	fprintf(nfile,"property float x\n");
	fprintf(nfile,"property float y\n");
	fprintf(nfile,"property float z\n");
	fprintf(nfile,"property uchar red\n");
	fprintf(nfile,"property uchar green\n");
	fprintf(nfile,"property uchar blue\n");
	fprintf(nfile,"element face %d\n",mesh.Faces.size());
	fprintf(nfile,"property list int int vertex_indices\n");
	fprintf(nfile,"end_header\n");
	for(size_t i=0;i<mesh.Vertices.size();i++)
	{
		AWriteV(nfile,mesh.Vertices[i].X,mesh.Vertices[i].Y,mesh.Vertices[i].Z,255,255,255);
	}
	for(size_t i=0;i<mesh.Faces.size();i++)
	{
		AWriteF(nfile,mesh.Faces[i].P0Index,mesh.Faces[i].P1Index,mesh.Faces[i].P2Index);
	}
	fclose(nfile);
}
#endif