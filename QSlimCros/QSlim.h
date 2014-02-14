#ifndef QSLIM_H
#define QSLIM_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <limits.h>
using namespace std;
#define M_PI 3.14159265358979323846
#define FEQ_EPS 1e-6
#define FEQ_EPS2 1e-12
inline bool FEQ(double a,double b,double eps=FEQ_EPS) { return fabs(a-b)<eps; }
inline bool FEQ(float a,float b,float eps=FEQ_EPS) { return fabsf(a-b)<eps; }
#define MIN(a,b) (((a)>(b))?(b):(a))
#define MAX(a,b) (((a)>(b))?(a):(b))
typedef double real;
enum Axis {X=0, Y=1, Z=2, W=3};
template<class T>
inline T min(T a, T b) { return (a < b)?a:b; }
template<class T>
inline T max(T a, T b) { return (a > b)?a:b; }
class Vec2 {
private:
	real elt[2];
protected:
	inline void copy(const Vec2& v);
public:
	Vec2(real x=0, real y=0) { elt[0]=x; elt[1]=y; }
	Vec2(const Vec2& v) { copy(v); }
	Vec2(const real *v) { elt[0]=v[0]; elt[1]=v[1]; }
	real& operator()(int i)       { return elt[i]; }
	real  operator()(int i) const { return elt[i]; }
	real& operator[](int i)       { return elt[i]; }
	real  operator[](int i) const { return elt[i]; }
	real *raw()             { return elt; }
	const real *raw() const { return elt; }
	inline bool operator==(const Vec2& v) const;
	inline bool operator!=(const Vec2& v) const;
	inline void set(real x, real y) { elt[0]=x; elt[1]=y; }
	inline Vec2& operator=(const Vec2& v);
	inline Vec2& operator+=(const Vec2& v);
	inline Vec2& operator-=(const Vec2& v);
	inline Vec2& operator*=(real s);
	inline Vec2& operator/=(real s);
	inline Vec2 operator+(const Vec2& v) const;
	inline Vec2 operator-(const Vec2& v) const;
	inline Vec2 operator-() const;
	inline Vec2 operator*(real s) const;
	inline Vec2 operator/(real s) const;
	inline real operator*(const Vec2& v) const;
};
inline void Vec2::copy(const Vec2& v)
{
	elt[0]=v.elt[0]; elt[1]=v.elt[1];
}
inline bool Vec2::operator==(const Vec2& v) const
{
	real dx=elt[X]-v[X],  dy=elt[Y]-v[Y];
	return (dx*dx + dy*dy) < FEQ_EPS2;
}
inline bool Vec2::operator!=(const Vec2& v) const
{
	real dx=elt[X]-v[X],  dy=elt[Y]-v[Y];
	return (dx*dx + dy*dy) > FEQ_EPS2;
}
inline Vec2& Vec2::operator=(const Vec2& v)
{
	copy(v);
	return *this;
}
inline Vec2& Vec2::operator+=(const Vec2& v)
{
	elt[0] += v[0];   elt[1] += v[1];
	return *this;
}
inline Vec2& Vec2::operator-=(const Vec2& v)
{
	elt[0] -= v[0];   elt[1] -= v[1];
	return *this;
}
inline Vec2& Vec2::operator*=(real s)
{
	elt[0] *= s;   elt[1] *= s;
	return *this;
}
inline Vec2& Vec2::operator/=(real s)
{
	elt[0] /= s;   elt[1] /= s;
	return *this;
}
inline Vec2 Vec2::operator+(const Vec2& v) const
{
	return Vec2(elt[0]+v[0], elt[1]+v[1]);
}
inline Vec2 Vec2::operator-(const Vec2& v) const
{
	return Vec2(elt[0]-v[0], elt[1]-v[1]);
}
inline Vec2 Vec2::operator-() const
{
	return Vec2(-elt[0], -elt[1]);
}
inline Vec2 Vec2::operator*(real s) const
{
	return Vec2(elt[0]*s, elt[1]*s);
}
inline Vec2 Vec2::operator/(real s) const
{
	return Vec2(elt[0]/s, elt[1]/s);
}
inline real Vec2::operator*(const Vec2& v) const
{
	return elt[0]*v[0] + elt[1]*v[1];
}
inline Vec2 operator*(real s, const Vec2& v) { return v*s; }
inline real norm(const Vec2& v) { return sqrt(v[0]*v[0] + v[1]*v[1]); }
inline real norm2(const Vec2& v) { return v[0]*v[0] + v[1]*v[1]; }
inline real length(const Vec2& v) { return norm(v); }
inline real unitize(Vec2& v)
{
	real l=norm2(v);
	if( l!=1.0 && l!=0.0 )
	{
		l = sqrt(l);
		v /= l;
	}
	return l;
}
class Vec3 {
private:
	real elt[3];
protected:
	inline void copy(const Vec3& v);
public:
	Vec3(real x=0, real y=0, real z=0) { elt[0]=x; elt[1]=y; elt[2]=z; }
	Vec3(const Vec3& v) { copy(v); }
	Vec3(const real *v) { elt[0]=v[0]; elt[1]=v[1]; elt[2]=v[2]; }
	real& operator()(int i)       { return elt[i]; }
	real  operator()(int i) const { return elt[i]; }
	real& operator[](int i)       { return elt[i]; }
	real  operator[](int i) const { return elt[i]; }
	real *raw()             { return elt; }
	const real *raw() const { return elt; }
	inline bool operator==(const Vec3& v) const;
	inline bool operator!=(const Vec3& v) const;
	inline void set(real x, real y, real z) { elt[0]=x; elt[1]=y; elt[2]=z; }
	inline Vec3& operator=(const Vec3& v);
	inline Vec3& operator+=(const Vec3& v);
	inline Vec3& operator-=(const Vec3& v);
	inline Vec3& operator*=(real s);
	inline Vec3& operator/=(real s);
	inline Vec3 operator+(const Vec3& v) const;
	inline Vec3 operator-(const Vec3& v) const;
	inline Vec3 operator-() const;
	inline Vec3 operator*(real s) const;
	inline Vec3 operator/(real s) const;
	inline real operator*(const Vec3& v) const;
	inline Vec3 operator^(const Vec3& v) const;
};
inline void Vec3::copy(const Vec3& v)
{
	elt[0]=v.elt[0]; elt[1]=v.elt[1]; elt[2]=v.elt[2];
}
inline bool Vec3::operator==(const Vec3& v) const
{
	real dx=elt[X]-v[X],  dy=elt[Y]-v[Y],  dz=elt[Z]-v[Z];
	return (dx*dx + dy*dy + dz*dz) < FEQ_EPS2;
}
inline bool Vec3::operator!=(const Vec3& v) const
{
	real dx=elt[X]-v[X],  dy=elt[Y]-v[Y],  dz=elt[Z]-v[Z];
	return (dx*dx + dy*dy + dz*dz) > FEQ_EPS2;
}
inline Vec3& Vec3::operator=(const Vec3& v)
{
	copy(v);
	return *this;
}
inline Vec3& Vec3::operator+=(const Vec3& v)
{
	elt[0] += v[0];   elt[1] += v[1];   elt[2] += v[2];
	return *this;
}
inline Vec3& Vec3::operator-=(const Vec3& v)
{
	elt[0] -= v[0];   elt[1] -= v[1];   elt[2] -= v[2];
	return *this;
}
inline Vec3& Vec3::operator*=(real s)
{
	elt[0] *= s;   elt[1] *= s;   elt[2] *= s;
	return *this;
}
inline Vec3& Vec3::operator/=(real s)
{
	elt[0] /= s;   elt[1] /= s;   elt[2] /= s;
	return *this;
}
inline Vec3 Vec3::operator+(const Vec3& v) const
{
	return Vec3(elt[0]+v[0], elt[1]+v[1], elt[2]+v[2]);
}
inline Vec3 Vec3::operator-(const Vec3& v) const
{
	return Vec3(elt[0]-v[0], elt[1]-v[1], elt[2]-v[2]);
}
inline Vec3 Vec3::operator-() const
{
	return Vec3(-elt[0], -elt[1], -elt[2]);
}
inline Vec3 Vec3::operator*(real s) const
{
	return Vec3(elt[0]*s, elt[1]*s, elt[2]*s);
}
inline Vec3 Vec3::operator/(real s) const
{
	return Vec3(elt[0]/s, elt[1]/s, elt[2]/s);
}
inline real Vec3::operator*(const Vec3& v) const
{
	return elt[0]*v[0] + elt[1]*v[1] + elt[2]*v[2];
}
inline Vec3 Vec3::operator^(const Vec3& v) const
{
	Vec3 w( elt[1]*v[2] - v[1]*elt[2],
		-elt[0]*v[2] + v[0]*elt[2],
		elt[0]*v[1] - v[0]*elt[1] );
	return w;
}
inline Vec3 operator*(real s, const Vec3& v) { return v*s; }
inline real norm(const Vec3& v)
{
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
inline real norm2(const Vec3& v)
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}
inline real length(const Vec3& v) { return norm(v); }
inline real unitize(Vec3& v)
{
	real l=norm2(v);
	if( l!=1.0 && l!=0.0 )
	{
		l = sqrt(l);
		v /= l;
	}
	return l;
}
class Vec4 {
private:
	real elt[4];
protected:
	inline void copy(const Vec4& v);
public:
	Vec4(real x=0, real y=0, real z=0, real w=0) {
		elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=w;
	}
	Vec4(const Vec4& v) { copy(v); }
	Vec4(const Vec3& v,real w) {elt[0]=v[0];elt[1]=v[1];elt[2]=v[2];elt[3]=w;}
	Vec4(const real *v) { elt[0]=v[0]; elt[1]=v[1]; elt[2]=v[2]; elt[3]=v[3]; }
	real& operator()(int i)       { return elt[i]; }
	real  operator()(int i) const { return elt[i]; }
	real& operator[](int i)             { return elt[i]; }
	const real& operator[](int i) const { return elt[i]; }
	real *raw()             { return elt; }
	const real *raw() const { return elt; }
	inline bool operator==(const Vec4&) const;
	inline bool operator!=(const Vec4&) const;
	inline void set(real x, real y, real z, real w){
		elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=w;
	}
	inline Vec4& operator=(const Vec4& v);
	inline Vec4& operator+=(const Vec4& v);
	inline Vec4& operator-=(const Vec4& v);
	inline Vec4& operator*=(real s);
	inline Vec4& operator/=(real s);
	inline Vec4 operator+(const Vec4& v) const;
	inline Vec4 operator-(const Vec4& v) const;
	inline Vec4 operator-() const;
	inline Vec4 operator*(real s) const;
	inline Vec4 operator/(real s) const;
	inline real operator*(const Vec4& v) const;
};
inline void Vec4::copy(const Vec4& v)
{
	elt[0]=v.elt[0]; elt[1]=v.elt[1]; elt[2]=v.elt[2]; elt[3]=v.elt[3];
}
inline bool Vec4::operator==(const Vec4& v) const
{
	real dx=elt[X]-v[X],  dy=elt[Y]-v[Y],  dz=elt[Z]-v[Z],  dw=elt[W]-v[W];
	return (dx*dx + dy*dy + dz*dz + dw*dw) < FEQ_EPS2;
}
inline bool Vec4::operator!=(const Vec4& v) const
{
	real dx=elt[X]-v[X],  dy=elt[Y]-v[Y],  dz=elt[Z]-v[Z],  dw=elt[W]-v[W];
	return (dx*dx + dy*dy + dz*dz + dw*dw) > FEQ_EPS2;
}
inline Vec4& Vec4::operator=(const Vec4& v)
{
	copy(v);
	return *this;
}
inline Vec4& Vec4::operator+=(const Vec4& v)
{
	elt[0] += v[0];   elt[1] += v[1];   elt[2] += v[2];   elt[3] += v[3];
	return *this;
}
inline Vec4& Vec4::operator-=(const Vec4& v)
{
	elt[0] -= v[0];   elt[1] -= v[1];   elt[2] -= v[2];   elt[3] -= v[3];
	return *this;
}
inline Vec4& Vec4::operator*=(real s)
{
	elt[0] *= s;   elt[1] *= s;   elt[2] *= s;  elt[3] *= s;
	return *this;
}
inline Vec4& Vec4::operator/=(real s)
{
	elt[0] /= s;   elt[1] /= s;   elt[2] /= s;   elt[3] /= s;
	return *this;
}
inline Vec4 Vec4::operator+(const Vec4& v) const
{
	return Vec4(elt[0]+v[0], elt[1]+v[1], elt[2]+v[2], elt[3]+v[3]);
}
inline Vec4 Vec4::operator-(const Vec4& v) const
{
	return Vec4(elt[0]-v[0], elt[1]-v[1], elt[2]-v[2], elt[3]-v[3]);
}
inline Vec4 Vec4::operator-() const
{
	return Vec4(-elt[0], -elt[1], -elt[2], -elt[3]);
}
inline Vec4 Vec4::operator*(real s) const
{
	return Vec4(elt[0]*s, elt[1]*s, elt[2]*s, elt[3]*s);
}
inline Vec4 Vec4::operator/(real s) const
{
	return Vec4(elt[0]/s, elt[1]/s, elt[2]/s, elt[3]/s);
}
inline real Vec4::operator*(const Vec4& v) const
{
	return elt[0]*v[0] + elt[1]*v[1] + elt[2]*v[2] + elt[3]*v[3];
}
inline Vec4 operator*(real s, const Vec4& v) { return v*s; }
inline Vec4 cross(const Vec4& a, const Vec4& b, const Vec4& c)
{
	Vec4 result;

	real d1 = (b[Z] * c[W]) - (b[W] * c[Z]);
	real d2 = (b[Y] * c[W]) - (b[W] * c[Y]);
	real d3 = (b[Y] * c[Z]) - (b[Z] * c[Y]);
	real d4 = (b[X] * c[W]) - (b[W] * c[X]);
	real d5 = (b[X] * c[Z]) - (b[Z] * c[X]);
	real d6 = (b[X] * c[Y]) - (b[Y] * c[X]);

	result[X] = - a[Y] * d1 + a[Z] * d2 - a[W] * d3;
	result[Y] =   a[X] * d1 - a[Z] * d4 + a[W] * d5;
	result[Z] = - a[X] * d2 + a[Y] * d4 - a[W] * d6;
	result[W] =   a[X] * d3 - a[Y] * d5 + a[Z] * d6;

	return result;
}
inline real norm(const Vec4& v)
{
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
}
inline real norm2(const Vec4& v)
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
}
inline real length(const Vec4& v) { return norm(v); }
inline real unitize(Vec4& v)
{
	real l=norm2(v);
	if( l!=1.0 && l!=0.0 )
	{
		l = sqrt(l);
		v /= l;
	}
	return l;
}
template<class T>
class array {
protected:
	T *data;
	int len;
public:
	array() { data=NULL; len=0; }
	array(int l) { init(l); }
	~array() { free(); }
	inline void init(int l);
	inline void free();
	inline void resize(int l);
	inline T& ref(int i);
	inline T& operator[](int i) { return data[i]; }
	inline T& operator()(int i) { return ref(i); }
	inline const T& ref(int i) const;
	inline const T& operator[](int i) const { return data[i]; }
	inline const T& operator()(int i) const { return ref(i); }
	inline int length() const { return len; }
	inline int maxLength() const { return len; }
};
template<class T>
inline void array<T>::init(int l)
{
	data = new T[l];
	len = l;
}
template<class T>
inline void array<T>::free()
{
	if( data )
	{
		delete[] data;
		data = NULL;
	}
}
template<class T>
inline T& array<T>::ref(int i)
{
	return data[i];
}
template<class T>
inline const T& array<T>::ref(int i) const
{
	return data[i];
}
template<class T>
inline void array<T>::resize(int l)
{
	T *old = data;
	data = new T[l];
	data = (T *)memcpy(data,old,MIN(len,l)*sizeof(T));
	len = l;
	delete[] old;
}
extern Vec3 randomPoint(const Vec3&, const Vec3&);  // on segment
extern Vec3 randomPoint(const Vec3&, const Vec3&, const Vec3&); // in triangle
extern real triangleArea(const Vec3&, const Vec3&, const Vec3&);
extern real triangleCompactness(const Vec3&, const Vec3&, const Vec3&);
class Bounds
{
public:
	Vec3 min, max;
	Vec3 center;
	real radius;
	unsigned int points;
	Bounds() { reset(); }
	void reset();
	void addPoint(const Vec3&);
	void complete();
};
class Plane
{
	Vec3 n;
	real d;
public:
	Plane() : n(0,0,1) { d=0; } // -- this will define the XY plane
	Plane(const Vec3& p, const Vec3& q, const Vec3& r) { calcFrom(p,q,r); }
	Plane(const array<Vec3>& verts) { calcFrom(verts); }
	Plane(const Plane& p) { n=p.n; d=p.d; }
	void calcFrom(const Vec3& p, const Vec3& q, const Vec3& r);
	void calcFrom(const array<Vec3>&);
	bool isValid() const { return n[X]!=0.0 || n[Y]!=0.0 || n[Z]!= 0.0; }
	void markInvalid() { n[X] = n[Y] = n[Z] = 0.0; }
	real distTo(const Vec3& p) const { return n*p + d; }
	const Vec3& normal() const { return n; }
	void coeffs(real *a, real *b, real *c, real *dd) const {
		*a=n[X]; *b=n[Y]; *c=n[Z]; *dd=d;
	}
	Vec4 coeffs() const { return Vec4(n,d); }
};
class Face3
{
protected:
	Plane P;
private:
	void recalcPlane() { P.calcFrom(vertexPos(0),vertexPos(1),vertexPos(2)); }
	void recalcPlane(const Vec3& a,const Vec3& b,const Vec3& c)
	{ P.calcFrom(a,b,c); }
public:
	Face3(const Vec3& a,const Vec3& b,const Vec3& c)
		: P(a,b,c)
	{ }
	virtual const Vec3& vertexPos(int i) const = 0;
	virtual void vertexPos(int i, const Vec3&) = 0; 
	const Plane& plane() { if(!P.isValid()) recalcPlane();   return P;}
	void invalidatePlane() { P.markInvalid(); }
	real distTo(const Vec3& p) const;
	real area();
};
inline real random1()
{
	//#if defined(WIN32) || defined(GFX_USE_RAND)
	return (real)rand() / (real)LONG_MAX;
	//#else
	//	return (real)random() / (real)LONG_MAX;
	//#endif
}
extern Vec3 randomPoint(const Vec3& v1, const Vec3& v2);
extern Vec3 randomPoint(const Vec3& v1, const Vec3& v2, const Vec3& v3);
extern  real triangleArea(const Vec3& v1, const Vec3& v2, const Vec3& v3);
#define ROOT3 1.732050807568877
#define FOUR_ROOT3 6.928203230275509
extern  real triangleCompactness(const Vec3& v1, const Vec3& v2, const Vec3& v3);
extern real __gfx_hoppe_dist(const Face3& f, const Vec3& v);
#define DOTP(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
template<class T>
class buffer : public array<T> {
protected:
	int fill;
public:
	buffer() { init(8); }
	buffer(int l) { init(l); }
	inline void init(int l) { array<T>::init(l); fill=0; }
	inline int add(const T& t);
	inline void reset();
	inline int find(const T&);
	inline T remove(int i);
	inline int addAll(const buffer<T>& buf);
	inline void removeDuplicates();
	inline int length() const { return fill; }
	inline int maxLength() const { return len; }
};
template<class T>
inline int buffer<T>::add(const T& t)
{
	if( fill == len )
		resize( len*2 );

	data[fill] = t;

	return fill++;
}
template<class T>
inline void buffer<T>::reset()
{
	fill = 0;
}
template<class T>
inline int buffer<T>::find(const T& t)
{
	for(int i=0;i<fill;i++)
		if( data[i] == t )
			return i;

	return -1;
}
template<class T>
inline T buffer<T>::remove(int i)
{
	fill--;
	T temp = data[i];
	data[i] = data[fill];
	return temp;
}
template<class T>
inline int buffer<T>::addAll(const buffer<T>& buf)
{
	for(int i=0; i<buf.fill; i++)
		add(buf(i));

	return fill;
}
template<class T>
inline void buffer<T>::removeDuplicates()
{
	for(int i=0; i<fill; i++)
	{
		for(int j=i+1; j<fill; )
		{
			if( data[j] == data[i] )
				remove(j);
			else
				j++;
		}
	}
}
#define NOT_IN_HEAP -47
class Heapable
{
private:
	int token;
public:
	Heapable() { notInHeap(); }
	inline int isInHeap() { return token!=NOT_IN_HEAP; }
	inline void notInHeap() { token = NOT_IN_HEAP; }
	inline int getHeapPos() { return token; }
	inline void setHeapPos(int t) { token=t; }
};
class heap_node {
public:
	float import;
	Heapable *obj;

	heap_node() { obj=NULL; import=0.0; }
	heap_node(Heapable *t, float i=0.0) { obj=t; import=i; }
	heap_node(const heap_node& h) { import=h.import; obj=h.obj; }
};
class Heap : public array<heap_node> {
	int size;
	void swap(int i, int j);
	int parent(int i) { return (i-1)/2; }
	int left(int i) { return 2*i+1; }
	int right(int i) { return 2*i+2; }
	void upheap(int i);
	void downheap(int i);
public:
	Heap() { size=0; }
	Heap(int s) : array<heap_node>(s) { size=0; }
	void insert(Heapable *, float);
	void update(Heapable *, float);
	heap_node *extract();
	heap_node *top() { return size<1 ? (heap_node *)NULL : &ref(0); }
	heap_node *kill(int i);
};

class NPrim
{
public:
	int uniqID;
	inline bool isValid()     { return uniqID >= 0; }
	inline void markInvalid() { if( uniqID>=0 ) uniqID = -uniqID-1; }
	inline void markValid()   { if( uniqID<0  ) uniqID = -uniqID-1; }
	inline int  validID()     { return (uniqID<0)?(-uniqID-1):uniqID; }
};
class NTaggedPrim : public NPrim
{
public:
	int tempID;
	inline void untag() { tempID = 0; }
	inline void tag(int t=1) { tempID = t; }
	inline bool isTagged() { return tempID!=0; }
};
class Vertex;
class Edge;
class Face;
typedef buffer<Vertex *> vert_buffer;
typedef buffer<Edge *> edge_buffer;
typedef buffer<Face *> face_buffer;
extern void untagFaceLoop(Vertex *v);
extern void collectFaceLoop(Vertex *v, face_buffer& faces);
#define EDGE_BOGUS 0
#define EDGE_BORDER 1
#define EDGE_MANIFOLD 2
#define EDGE_NONMANIFOLD 3
extern int classifyEdge(Edge *);
#define VERTEX_INTERIOR 0
#define VERTEX_BORDER 1
#define VERTEX_BORDER_ONLY 2
extern int classifyVertex(Vertex *);
class Vertex : public Vec3, public NTaggedPrim
{
	edge_buffer edge_uses;
public:
	Vertex(real x, real y, real z) : Vec3(x, y, z), edge_uses(6) {
	}
	void kill();
	edge_buffer& edgeUses() { return edge_uses; }
	void linkEdge(Edge *);
	void unlinkEdge(Edge *);
	void remapTo(Vertex *v);
};
class Edge : public NPrim
{
private:
	Vertex *v1;
	face_buffer *face_uses;
	Edge *twin;
	Edge(Edge *twin, Vertex *org); // the twin constructor
public:
	Edge(Vertex *, Vertex *);
	~Edge();
	Vertex *org()  { return v1;       }
	Vertex *dest() { return twin->v1; }
	Edge *sym()    { return twin;     }
	void kill();
	face_buffer& faceUses() { return *face_uses; }
	void linkFace(Face *);
	void unlinkFace(Face *);
	void remapEndpoint(Vertex *from, Vertex *to);
	void remapTo(Edge *);
};
class Face : public Face3, public NTaggedPrim
{
	Edge *edges[3];
public:
	Face(Edge *, Edge *, Edge *);
	const Vec3& vertexPos(int i) const { return *edges[i]->org(); }
	void vertexPos(int, const Vec3&) {
		//fatal_error("Face: can't directly set vertex position.");
	}
	Vertex *vertex(int i) { return edges[i]->org(); }
	const Vertex *vertex(int i) const { return edges[i]->org(); }
	Edge *edge(int i)               { return edges[i];        }
	void kill();
	void remapEdge(Edge *from, Edge *to);
};
class Mat4
{
private:
	Vec4 row[4];
protected:
	inline void copy(const Mat4& m);
	inline Vec4 col(int i) const
	{ return Vec4(row[0][i],row[1][i],row[2][i],row[3][i]); }
public:
	static Mat4 identity;
	static Mat4 zero;
	static Mat4 unit;
	static Mat4 trans(real,real,real);
	static Mat4 scale(real,real,real);
	static Mat4 xrot(real); //
	static Mat4 yrot(real); // Arguments are in radians
	static Mat4 zrot(real); //
	Mat4() { copy(zero); }
	Mat4(const Vec4& r0,const Vec4& r1,const Vec4& r2,const Vec4& r3)
	{ row[0]=r0; row[1]=r1; row[2]=r2; row[3]=r3; }
	Mat4(const Mat4& m) { copy(m); }
	real& operator()(int i, int j)       { return row[i][j]; }
	real  operator()(int i, int j) const { return row[i][j]; }
	const Vec4& operator[](int i) const { return row[i]; }
	inline int operator==(const Mat4&);
	inline Mat4& operator=(const Mat4& m) { copy(m); return *this; }
	inline Mat4& operator+=(const Mat4& m);
	inline Mat4& operator-=(const Mat4& m);
	inline Mat4& operator*=(real s);
	inline Mat4& operator/=(real s);
	inline Mat4 operator+(const Mat4& m) const;
	inline Mat4 operator-(const Mat4& m) const;
	inline Mat4 operator-() const;
	inline Mat4 operator*(real s) const;
	inline Mat4 operator/(real s) const;
	Mat4 operator*(const Mat4& m) const;
	inline Vec4 operator*(const Vec4& v) const; // [x y z w]
	inline Vec3 operator*(const Vec3& v) const; // [x y z w]
	real det() const;
	Mat4 transpose() const;
	Mat4 adjoint() const;
	real inverse(Mat4&) const;
	real cramerInverse(Mat4&) const;
	friend ostream& operator<<(ostream&, const Mat4&);
	friend istream& operator>>(istream&, Mat4&);
};
inline void Mat4::copy(const Mat4& m)
{
	row[0] = m.row[0]; row[1] = m.row[1];
	row[2] = m.row[2]; row[3] = m.row[3];
}
inline int Mat4::operator==(const Mat4& m)
{
	return row[0]==m.row[0] &&
		row[1]==m.row[1] &&
		row[2]==m.row[2] &&
		row[3]==m.row[3] ;
}
inline Mat4& Mat4::operator+=(const Mat4& m)
{
	row[0] += m.row[0]; row[1] += m.row[1];
	row[2] += m.row[2]; row[3] += m.row[3];
	return *this;
}
inline Mat4& Mat4::operator-=(const Mat4& m)
{
	row[0] -= m.row[0]; row[1] -= m.row[1];
	row[2] -= m.row[2]; row[3] -= m.row[3];
	return *this;
}
inline Mat4& Mat4::operator*=(real s)
{
	row[0] *= s; row[1] *= s; row[2] *= s; row[3] *= s;
	return *this;
}
inline Mat4& Mat4::operator/=(real s)
{
	row[0] /= s; row[1] /= s; row[2] /= s; row[3] /= s;
	return *this;
}
inline Mat4 Mat4::operator+(const Mat4& m) const
{
	return Mat4(row[0]+m.row[0],
		row[1]+m.row[1],
		row[2]+m.row[2],
		row[3]+m.row[3]);
}
inline Mat4 Mat4::operator-(const Mat4& m) const
{
	return Mat4(row[0]-m.row[0],
		row[1]-m.row[1],
		row[2]-m.row[2],
		row[3]-m.row[3]);
}
inline Mat4 Mat4::operator-() const
{
	return Mat4(-row[0], -row[1], -row[2], -row[3]);
}
inline Mat4 Mat4::operator*(real s) const
{
	return Mat4(row[0]*s, row[1]*s, row[2]*s, row[3]*s);
}
inline Mat4 Mat4::operator/(real s) const
{
	return Mat4(row[0]/s, row[1]/s, row[2]/s, row[3]/s);
}
inline Vec4 Mat4::operator*(const Vec4& v) const
{
	return Vec4(row[0]*v, row[1]*v, row[2]*v, row[3]*v);
}
inline Vec3 Mat4::operator*(const Vec3& v) const
{
	Vec4 u=Vec4(v,1);
	real w=row[3]*u;

	if(w==0.0)
		return Vec3(row[0]*u, row[1]*u, row[2]*u);
	else
		return Vec3(row[0]*u/w, row[1]*u/w, row[2]*u/w);
}
template<class T>
class array3
{
protected:
	T *data;
	int w, h, d;
public:
	array3() { data=NULL; w=h=d=0; }
	array3(int w, int h, int d) { init(w,h,d); }
	~array3() { free(); }
	inline void init(int w, int h, int d);
	inline void free();
	inline T& ref(int i, int j, int k);
	inline T& operator()(int i,int j, int k) { return ref(i,j, k); }
	inline int width() const { return w; }
	inline int height() const { return h; }
	inline int depth() const { return d; }
};
template<class T>
inline void array3<T>::init(int width, int height, int depth)
{
	w = width;
	h = height;
	d = depth;
	data = new T[w*h*d];
}
template<class T>
inline void array3<T>::free()
{
	if( data )
	{
		delete[] data;
		data = NULL;
		w = h = d = 0;
	}

}
template<class T>
inline T& array3<T>::ref(int i, int j, int k)
{
	return data[k*w*h + j*w + i];
}
class ProxGrid_Cell : public buffer<Vec3 *>
{
public:
	ProxGrid_Cell(): buffer<Vec3 *>(8)
	{ }
};
class ProxGrid
{
	array3<ProxGrid_Cell> cells;
	int xdiv, ydiv, zdiv;
	real cellsize;
	real cellsize2;
	Vec3 min, max;
	void cell_for_point(const Vec3&, int *i, int *j, int *k);
	void maybe_collect_points(Vec3 *v, buffer<Vec3 *>& close,
		ProxGrid_Cell& cell);
public:
	ProxGrid(const Vec3& min, const Vec3& max, real dist);
	~ProxGrid() { cells.free(); }
	void addPoint(Vec3 *);
	void removePoint(Vec3 *);
	void proximalPoints(Vec3 *, buffer<Vec3 *>&);
};

extern int face_target;
extern real error_tolerance;
extern bool will_use_plane_constraint;
extern bool will_use_vertex_constraint;
extern bool will_preserve_boundaries;
extern bool will_preserve_mesh_quality;
extern bool will_constrain_boundaries;
extern real boundary_constraint_weight;
extern bool will_weight_by_area;
#define PLACE_ENDPOINTS 0
#define PLACE_ENDORMID  1
#define PLACE_LINE      2
#define PLACE_OPTIMAL   3
extern int placement_policy;
extern real pair_selection_tolerance;

class Model //: public SMF_Model
{
protected:
	vert_buffer vertices;
	edge_buffer edges;
	face_buffer faces;
private:
	void maybeFixFace(Face *);
public:
	Model() { }
	Bounds bounds;
	int validVertCount;
	int validEdgeCount;
	int validFaceCount;
	Vertex *vertex(int i) { return vertices(i); }
	Edge *edge(int i) { return edges(i); }
	Face *face(int i) { return faces(i); }
	int vertCount() { return vertices.length(); }
	int edgeCount() { return edges.length();    }
	int faceCount() { return faces.length();    }
	vert_buffer& allVertices() { return vertices; }
	edge_buffer& allEdges()    { return edges;    }
	face_buffer& allFaces()    { return faces;    }
	Vertex   *newVertex(real x=0.0, real y=0.0, real z=0.0);
	Edge     *newEdge(Vertex *,Vertex *);
	Face *newFace(Vertex *, Vertex *, Vertex *);
	void killVertex(Vertex *);
	void killEdge(Edge *);
	void killFace(Face *);
	void reshapeVertex(Vertex *, real, real, real);
	void remapVertex(Vertex *from, Vertex *to);
	void contract(Vertex *v1, Vertex *v2, const Vec3& to,face_buffer& changed);
	void removeDegeneracy(face_buffer& changed);
	void contractionRegion(Vertex *v1, Vertex *v2, face_buffer& changed);
	void contractionRegion(Vertex *v1,const vert_buffer& vertices,face_buffer& changed);
	int in_Vertex(const Vec3&);
	int in_Face(int v1, int v2, int v3);
	Vec3 synthesizeNormal(Vertex *);
};
extern Model M0;
extern Mat4 quadrix_vertex_constraint(const Vec3&);
extern Mat4 quadrix_plane_constraint(real a, real b, real c, real d);
extern Mat4 quadrix_plane_constraint(Face& T);
extern Mat4 quadrix_plane_constraint(const Vec3& n, real);
extern Mat4 quadrix_plane_constraint(const Vec3&, const Vec3&, const Vec3&);
extern real quadrix_evaluate_vertex(const Vec3& v, const Mat4& K);
extern bool check_for_discontinuity(Edge *);
extern Mat4 quadrix_discontinuity_constraint(Edge *, const Vec3&);
extern Mat4 quadrix_discontinuity_constraint(Edge *);
extern bool quadrix_find_local_fit(const Mat4& Q,const Vec3& v1, const Vec3& v2,Vec3& candidate);
extern bool quadrix_find_line_fit(const Mat4& Q,const Vec3& v1, const Vec3& v2,Vec3& candidate);
extern bool quadrix_find_best_fit(const Mat4& Q, Vec3& candidate);
extern real quadrix_pair_target(const Mat4& Q,Vertex *v1,Vertex *v2,Vec3& candidate);

extern Vertex *decimate_last_v0;
extern Vertex *decimate_last_v1;
extern bool decimate_quadric(Vertex *v, Mat4& Q);
extern real decimate_min_error();
extern real decimate_max_error(Model& m);
extern real decimate_error(Vertex *);
extern void decimate_contract(Model& m);
extern void decimate_init(Model& m, real limit);



#endif