#include "QSlim.h"

Vec3 randomPoint(const Vec3& v1, const Vec3& v2)
{
	real a = random1();

	return a*v1 + (1-a)*v2;
}
Vec3 randomPoint(const Vec3& v1, const Vec3& v2, const Vec3& v3)
{
	real b1 = 1 - sqrt( 1-random1() );
	real b2 = (1-b1) * random1();
	real b3 = 1 - b1 - b2;

	return b1*v1 + b2*v2 + b3*v3;
}
real triangleArea(const Vec3& v1, const Vec3& v2, const Vec3& v3)
{
	Vec3 a = v2 - v1;
	Vec3 b = v3 - v1;

	return 0.5 * length(a ^ b);
}
real triangleCompactness(const Vec3& v1, const Vec3& v2, const Vec3& v3)
{
	real L1 = norm2(v2-v1);
	real L2 = norm2(v3-v2);
	real L3 = norm2(v1-v3);

	return FOUR_ROOT3 * triangleArea(v1,v2,v3) / (L1+L2+L3);
}

void Heap::swap(int i,int j)
{
	heap_node tmp = ref(i);
	ref(i) = ref(j);
	ref(j) = tmp;
	ref(i).obj->setHeapPos(i);
	ref(j).obj->setHeapPos(j);
}
void Heap::upheap(int i)
{
	if( i==0 ) return;

	if( ref(i).import > ref(parent(i)).import ) {
		swap(i,parent(i));
		upheap(parent(i));
	}
}
void Heap::downheap(int i)
{
	if (i>=size) return;        // perhaps just extracted the last

	int largest = i,
		l = left(i),
		r = right(i);

	if( l<size && ref(l).import > ref(largest).import ) largest = l;
	if( r<size && ref(r).import > ref(largest).import ) largest = r;

	if( largest != i ) {
		swap(i,largest);
		downheap(largest);
	}
}
void Heap::insert(Heapable *t,float v)
{
	if( size == maxLength() )
	{
		cerr << "NOTE: Growing heap from " << size << " to " << 2*size << endl;
		resize(2*size);
	}

	int i = size++;

	ref(i).obj = t;
	ref(i).import = v;

	ref(i).obj->setHeapPos(i);

	upheap(i);
}
void Heap::update(Heapable *t,float v)
{
	int i = t->getHeapPos();

	if( i >= size )
	{
		cerr << "WARNING: Attempting to update past end of heap!" << endl;
		return;
	}
	else if( i == NOT_IN_HEAP )
	{
		cerr << "WARNING: Attempting to update object not in heap!" << endl;
		return;
	}

	float old=ref(i).import;
	ref(i).import = v;

	if( v<old )
		downheap(i);
	else
		upheap(i);
}
heap_node *Heap::extract()
{
	if( size<1 ) return 0;

	swap(0,size-1);
	size--;

	downheap(0);

	ref(size).obj->notInHeap();

	return &ref(size);
}
heap_node *Heap::kill(int i)
{
	if( i>=size )
		cerr << "WARNING: Attempt to delete invalid heap node." << endl;

	swap(i, size-1);
	size--;
	ref(size).obj->notInHeap();

	if( ref(i).import < ref(size).import )
		downheap(i);
	else
		upheap(i);


	return &ref(size);
}

Mat4 Mat4::identity(Vec4(1,0,0,0),Vec4(0,1,0,0),Vec4(0,0,1,0),Vec4(0,0,0,1));
Mat4 Mat4::zero(Vec4(0,0,0,0),Vec4(0,0,0,0),Vec4(0,0,0,0),Vec4(0,0,0,0));
Mat4 Mat4::unit(Vec4(1,1,1,1),Vec4(1,1,1,1),Vec4(1,1,1,1),Vec4(1,1,1,1));
Mat4 Mat4::trans(real x, real y, real z)
{
	return Mat4(Vec4(1,0,0,x),
		Vec4(0,1,0,y),
		Vec4(0,0,1,z),
		Vec4(0,0,0,1));
}
Mat4 Mat4::scale(real x, real y, real z)
{
	return Mat4(Vec4(x,0,0,0),
		Vec4(0,y,0,0),
		Vec4(0,0,z,0),
		Vec4(0,0,0,1));
}
Mat4 Mat4::xrot(real a)
{
	real c = cos(a);
	real s = sin(a);

	return Mat4(Vec4(1, 0, 0, 0),
		Vec4(0, c,-s, 0),
		Vec4(0, s, c, 0),
		Vec4(0, 0, 0, 1));
}
Mat4 Mat4::yrot(real a)
{
	real c = cos(a);
	real s = sin(a);
	return Mat4(Vec4(c, 0, s, 0),Vec4(0, 1, 0, 0),Vec4(-s,0, c, 0),Vec4(0, 0, 0, 1));
}
Mat4 Mat4::zrot(real a)
{
	real c = cos(a);
	real s = sin(a);
	return Mat4(Vec4(c,-s, 0, 0),Vec4(s, c, 0, 0),Vec4(0, 0, 1, 0),Vec4(0, 0, 0, 1));
}
Mat4 Mat4::operator*(const Mat4& m) const
{
	Mat4 A;
	int i,j;
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			A(i,j) = row[i]*m.col(j);

	return A;
}
real Mat4::det() const
{
	return row[0] * cross(row[1], row[2], row[3]);
}
Mat4 Mat4::transpose() const
{
	return Mat4(col(0), col(1), col(2), col(3));
}
Mat4 Mat4::adjoint() const
{
	Mat4 A;
	A.row[0] = cross( row[1], row[2], row[3]);
	A.row[1] = cross(-row[0], row[2], row[3]);
	A.row[2] = cross( row[0], row[1], row[3]);
	A.row[3] = cross(-row[0], row[1], row[2]);
	return A;
}
real Mat4::cramerInverse(Mat4& inv) const
{
	Mat4 A = adjoint();
	real d = A.row[0] * row[0];

	if( d==0.0 )
		return 0.0;

	inv = A.transpose() / d;
	return d;
}
#define SWAP(a, b, t)   {t = a; a = b; b = t;}
real Mat4::inverse(Mat4& B) const
{
	Mat4 A(*this);
	int i, j, k;
	real max, t, det, pivot;
	for (i=0; i<4; i++)                 /* put identity matrix in B */
		for (j=0; j<4; j++)
			B(i, j) = (real)(i==j);
	det = 1.0;
	for (i=0; i<4; i++) {               /* eliminate in column i, below diag */
		max = -1.;
		for (k=i; k<4; k++)             /* find pivot for column i */
			if (fabs(A(k, i)) > max) {
				max = fabs(A(k, i));
				j = k;
			}
			if (max<=0.) return 0.;         /* if no nonzero pivot, PUNT */
			if (j!=i) {                     /* swap rows i and j */
				for (k=i; k<4; k++)
					SWAP(A(i, k), A(j, k), t);
				for (k=0; k<4; k++)
					SWAP(B(i, k), B(j, k), t);
				det = -det;
			}
			pivot = A(i, i);
			det *= pivot;
			for (k=i+1; k<4; k++)           /* only do elems to right of pivot */
				A(i, k) /= pivot;
			for (k=0; k<4; k++)
				B(i, k) /= pivot;
			for (j=i+1; j<4; j++) {         /* eliminate in rows below i */
				t = A(j, i);                /* we're gonna zero this guy */
				for (k=i+1; k<4; k++)       /* subtract scaled row i from row j */
					A(j, k) -= A(i, k)*t;   /* (ignore k<=i, we know they're 0) */
				for (k=0; k<4; k++)
					B(j, k) -= B(i, k)*t;
			}
	}
	for (i=4-1; i>0; i--) {             /* eliminate in column i, above diag */
		for (j=0; j<i; j++) {           /* eliminate in rows above i */
			t = A(j, i);                /* we're gonna zero this guy */
			for (k=0; k<4; k++)         /* subtract scaled row i from row j */
				B(j, k) -= B(i, k)*t;
		}
	}
	return det;
}

ProxGrid::ProxGrid(const Vec3& lo, const Vec3& hi, real dist)
{
	cellsize = dist;
	cellsize2 = dist*dist;
	min = lo;
	max = hi;
	xdiv = (int)ceil((max[X] - min[X])/dist);
	ydiv = (int)ceil((max[Y] - min[Y])/dist);
	zdiv = (int)ceil((max[Z] - min[Z])/dist);
	cells.init(xdiv, ydiv, zdiv);
}
void ProxGrid::cell_for_point(const Vec3& v,int *i_out,int *j_out,int *k_out)
{
	int i = (int)floor((v[X] - min[X]) / cellsize);
	int j = (int)floor((v[Y] - min[Y]) / cellsize);
	int k = (int)floor((v[Z] - min[Z]) / cellsize);
	if( i==xdiv ) i--;
	if( j==ydiv ) j--;
	if( k==zdiv ) k--;
	*i_out = i;
	*j_out = j;
	*k_out = k;
}
void ProxGrid::addPoint(Vec3 *v)
{
	int i, j, k;
	cell_for_point(*v, &i, &j, &k);
	ProxGrid_Cell& cell = cells(i,j,k);
	cell.add(v);
}
void ProxGrid::removePoint(Vec3 *v)
{
	int i, j, k;
	cell_for_point(*v, &i, &j, &k);
	ProxGrid_Cell& cell = cells(i,j,k);
	int index = cell.find(v);
	if( index >= 0 )
		cell.remove(index);
	else
		cerr << "WARNING: ProxGrid -- removing non-member point." << endl;
}
void ProxGrid::maybe_collect_points(Vec3 *v, buffer<Vec3 *>& close,ProxGrid_Cell& cell)
{
	for(int i=0; i<cell.length(); i++)
	{
		Vec3 *u = cell(i);

		if( u!=v && norm2(*u - *v) < cellsize2 )
			close.add(u);
	}
}
void ProxGrid::proximalPoints(Vec3 *v, buffer<Vec3 *>& close)
{
	int i, j, k;
	cell_for_point(*v, &i, &j, &k);
	for(int dk=-1; dk<2; dk++)
		for(int dj=-1; dj<2; dj++)
			for(int di=-1; di<2; di++)
			{
				if( i+di>=0 && j+dj>=0 && k+dk>=0
					&& i+di<cells.width()
					&& j+dj<cells.height()
					&& k+dk<cells.depth() )
				{
					maybe_collect_points(v, close, cells(i+di, j+dj, k+dk));
				}
			}
}



static real Distance2(real x[3], real *y)
{
	real a, b, c;
	a = x[0] - y[0];  b = x[1] - y[1];  c = x[2] - y[2];
	return (a * a + b * b + c * c);
}
static void interp(real *proj, const real *p1, const real *p2,const real *p3, const real *bary)
{
	proj[0] = p1[0] * bary[0] + p2[0] * bary[1] + p3[0] * bary[2];
	proj[1] = p1[1] * bary[0] + p2[1] * bary[1] + p3[1] * bary[2];
	proj[2] = p1[2] * bary[0] + p2[2] * bary[1] + p3[2] * bary[2];
}
static void Projecth(const real *v1, const real *v2, const real *v3,real *bary)
{
	int     i;
	real   vvi[3],vppi[3];
	real   d12sq, don12,d2, mind2, a;
	real   proj[3];
	real   pf[3][3];
	real   ba[3];
	real   cba[3];

	mind2 = 1e30;
	interp(proj,v1,v2,v3,bary);
	pf[0][0] = v1[0]; pf[0][1] = v1[1]; pf[0][2] = v1[2];
	pf[1][0] = v2[0]; pf[1][1] = v2[1]; pf[1][2] = v2[2];
	pf[2][0] = v3[0]; pf[2][1] = v3[1]; pf[2][2] = v3[2];

	ba[0] = bary[0]; ba[1] = bary[1]; ba[2] = bary[2];

	for (i = 0; i < 3; i++){
		if (ba[(i+2)%3] >= 0) continue;
		/* project proj onto segment pf[(i+0)%3]--pf[(i+1)%3]  */

		vvi[0] = pf[(i+1) % 3][0] - pf[i][0];
		vvi[1] = pf[(i+1) % 3][1] - pf[i][1];
		vvi[2] = pf[(i+1) % 3][2] - pf[i][2];

		vppi[0] = proj[0] - pf[i][0];
		vppi[1] = proj[1] - pf[i][1];
		vppi[2] = proj[2] - pf[i][2];

		d12sq = DOTP(vvi, vvi);
		don12 = DOTP(vvi, vppi);

		if (don12<=0) {
			d2 = Distance2(pf[i], proj);
			if (d2 >= mind2) continue;
			mind2=d2; cba[i]=1; cba[(i+1)%3]=0; cba[(i+2)%3]=0;
		}
		else {
			if (don12 >= d12sq) {
				d2 = Distance2(pf[(i+1)%3], proj);
				if (d2>=mind2) continue;
				mind2=d2; cba[i]=0; cba[(i+1)%3]=1; cba[(i+2)%3]=0;
			}
			else {
				a = don12/d12sq;
				cba[i]=1-a; cba[(i+1)%3]=a; cba[(i+2)%3]=0;
				break;
			}
		}
	}

	bary[0] = cba[0]; bary[1] = cba[1]; bary[2] = cba[2];
}
static void ProjectPtri(const real *point,  const real *v1, const real *v2,  const real *v3, real *bary)
{
	int    i;
	real  localv2[3], localv3[3], vpp1[3];
	real  v22,v33,v23,v2pp1,v3pp1;
	real  a1,a2,a3,denom;

	for (i = 0; i < 3; i++){
		localv2[i] = v2[i] - v1[i];
		localv3[i] = v3[i] - v1[i];
		vpp1[i] = point[i] - v1[i];
	}

	v22   = DOTP(localv2, localv2);
	v33   = DOTP(localv3, localv3);
	v23   = DOTP(localv2, localv3);
	v2pp1 = DOTP(localv2, vpp1);
	v3pp1 = DOTP(localv3, vpp1);

	if (!v22) v22=1;        /* recover if v2==0 */
	if (!v33) v33=1;        /* recover if v3==0 */

	denom = ( v33 - v23 * v23 / v22);
	if (!denom) {
		a2 = a3 = 1.0/3.0;    /* recover if v23*v23==v22*v33 */
	}
	else {
		a3=(v3pp1-v23/v22*v2pp1)/denom;
		a2=(v2pp1-a3*v23)/v22;
	}
	a1 = 1 - a2 - a3;

	bary[0] = a1; bary[1] = a2; bary[2] = a3;

	if ((a1 < 0) || (a2 < 0) || (a3 < 0)){
		Projecth(v1,v2,v3,bary);
		return;
	}
}

real __gfx_hoppe_dist(const Face3& f, const Vec3& v)
{
	Vec3 bary;

	const Vec3& v0 = f.vertexPos(0);
	const Vec3& v1 = f.vertexPos(1);
	const Vec3& v2 = f.vertexPos(2);

	ProjectPtri(v.raw(), v0.raw(), v1.raw(), v2.raw(), bary.raw());

	Vec3 p = bary[X]*v0 + bary[Y]*v1 + bary[Z]*v2;

	Vec3 diff = v - p;

	return diff*diff;
}

void Bounds::reset()
{
	min[X] = min[Y] = min[Z] = HUGE;
	max[X] = max[Y] = max[Z] = -HUGE;

	center[X] = center[Y] = center[Z] = 0.0;
	radius = 0.0;

	points = 0;
}
void Bounds::addPoint(const Vec3& v)
{
	if( v[X] < min[X] ) min[X] = v[X];
	if( v[Y] < min[Y] ) min[Y] = v[Y];
	if( v[Z] < min[Z] ) min[Z] = v[Z];

	if( v[X] > max[X] ) max[X] = v[X];
	if( v[Y] > max[Y] ) max[Y] = v[Y];
	if( v[Z] > max[Z] ) max[Z] = v[Z];


	center += v;

	points++;
}
void Bounds::complete()
{
	center /= (real)points;
	Vec3 R1 = max-center;
	Vec3 R2 = min-center;
	radius = MAX(length(R1), length(R2));
}
void Plane::calcFrom(const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
	Vec3 v1 = p2-p1;
	Vec3 v2 = p3-p1;

	n = v1 ^ v2;
	unitize(n);

	d = -n*p1;
}
void Plane::calcFrom(const array<Vec3>& verts)
{
	n[X] = n[Y] = n[Z] = 0.0;

	int i;
	for(i=0; i<verts.length()-1; i++)
	{
		const Vec3& cur = verts[i];
		const Vec3& next = verts[i+1];

		n[X] += (cur[Y] - next[Y]) * (cur[Z] + next[Z]);
		n[Y] += (cur[Z] - next[Z]) * (cur[X] + next[X]);
		n[Z] += (cur[X] - next[X]) * (cur[Y] + next[Y]);
	}

	const Vec3& cur = verts[verts.length()-1];
	const Vec3& next = verts[0];
	n[X] += (cur[Y] - next[Y]) * (cur[Z] + next[Z]);
	n[Y] += (cur[Z] - next[Z]) * (cur[X] + next[X]);
	n[Z] += (cur[X] - next[X]) * (cur[Y] + next[Y]);

	unitize(n);

	d = -n*verts[0];
}
real Face3::area()
{
	return triangleArea(vertexPos(0),vertexPos(1),vertexPos(2));
}
real Face3::distTo(const Vec3& v) const
{
	return __gfx_hoppe_dist(*this, v);
}

void Vertex::kill()
{
	markInvalid();
	edge_uses.reset();
}
void Vertex::linkEdge(Edge *e)
{
	edge_uses.add(e);
}
void Vertex::unlinkEdge(Edge *e)
{
	int index = edge_uses.find(e);
	edge_uses.remove(index);
	if( edge_uses.length() <= 0 )
		kill();
}
void Vertex::remapTo(Vertex *v)
{
	if( v != this )
	{
		for(int i=0; i<edge_uses.length(); i++)
		{
			edge_uses(i)->remapEndpoint(this, v);
		}
		kill();
	}
}
Edge::Edge(Vertex *a, Vertex *b)
{
	v1 = a;
	v1->linkEdge(this);
	face_uses = new buffer<Face *>(2);
	twin = new Edge(this, b);
}
Edge::Edge(Edge *sibling, Vertex *org)
{
	v1 = org;
	v1->linkEdge(this);
	face_uses = sibling->face_uses;
	twin = sibling;
}
Edge::~Edge()
{
	if( twin )
	{
		face_uses->free();
		delete face_uses;
		twin->twin = NULL;
		delete twin;
	}
}
void Edge::kill()
{
	if( isValid() )
	{
		org()->unlinkEdge(this);
		dest()->unlinkEdge(sym());
		markInvalid();
		twin->markInvalid();
		face_uses->reset();
	}
}
void Edge::linkFace(Face *face)
{
	face_uses->add(face);
}
void Edge::unlinkFace(Face *face)
{
	int index = face_uses->find(face);
	face_uses->remove(index);
	if( face_uses->length() == 0 )
		kill();
}
void Edge::remapTo(Edge *e)
{
	if( e != this )
	{
		for(int i=0; i<face_uses->length(); i++)
		{
			(*face_uses)(i)->remapEdge(this, e);
		}
		kill();
	}
}
void Edge::remapEndpoint(Vertex *from, Vertex *to)
{
	if( org()==from )
	{
		v1 = to;
		to->linkEdge(this);
	}
	else if( dest()==from )
	{
		this->twin->v1=NULL;
		twin->v1 = to;
		to->linkEdge(twin);
	}
	else
	{
		cerr << "WARNING remapEndpoint: Illegal endpoint." << endl;
	}
	for(int i=0; i<face_uses->length(); i++)
	{
		face_uses->ref(i)->invalidatePlane();
	}
}
Face::Face(Edge *e0, Edge *e1, Edge *e2): Face3(*e0->org(), *e1->org(), *e2->org())
{
	edges[0] = e0;
	edges[1] = e1;
	edges[2] = e2;
	edges[0]->linkFace(this);
	edges[1]->linkFace(this);
	edges[2]->linkFace(this);
}
void Face::kill()
{
	if( isValid() )
	{
		if( edge(0)->isValid() )
			edge(0)->unlinkFace(this);

		if( edge(1)->isValid() )
			edge(1)->unlinkFace(this);

		if( edge(2)->isValid() )
			edge(2)->unlinkFace(this);

		markInvalid();
	}
}
void Face::remapEdge(Edge *from, Edge *to)
{
	for(int i=0; i<3; i++)
	{
		if( edges[i] == from )
		{
			edges[i] = to;
			to->linkFace(this);
		}
		else if( edges[i] == from->sym() )
		{
			edges[i] = to->sym();
			to->sym()->linkFace(this);
		}
	}

	invalidatePlane();
}
void untagFaceLoop(Vertex *v)
{
	edge_buffer& edges = v->edgeUses();
	for(int j=0; j<edges.length(); j++)
	{	    
		face_buffer& faces = edges(j)->faceUses();
		for(int k=0; k<faces.length(); k++)
			faces(k)->untag();
	}
}
void collectFaceLoop(Vertex *v, face_buffer& loop)
{
	edge_buffer& edges = v->edgeUses();

	for(int j=0; j<edges.length(); j++)
	{	    
		face_buffer& faces = edges(j)->faceUses();
		for(int k=0; k<faces.length(); k++)
			if( !faces(k)->isTagged() )
			{
				loop.add(faces(k));
				faces(k)->tag();
			}
	}
}
int classifyEdge(Edge *e)
{
	int cls = e->faceUses().length();
	if( cls>3 ) cls=3;
	return cls;
}
int classifyVertex(Vertex *v)
{
	int border_count = 0;
	const edge_buffer& edges = v->edgeUses();

	for(int i=0; i<edges.length(); i++)
		if( classifyEdge(edges(i)) == 1 )
			border_count++;

	if( border_count == edges.length() )
		return VERTEX_BORDER_ONLY;
	else 
		return (border_count > 0);
}


double setup_time, init_time, slim_time, write_time;
int face_target = 0;
real error_tolerance = HUGE;
bool will_use_plane_constraint = true;
bool will_use_vertex_constraint = false;
bool will_preserve_boundaries = false;
bool will_preserve_mesh_quality = false;
bool will_constrain_boundaries = false;
real boundary_constraint_weight = 1.0;
bool will_weight_by_area = false;
int placement_policy = PLACE_OPTIMAL;
real pair_selection_tolerance = 0.0;
Model M0;


int Model::in_Vertex(const Vec3& p)
{
	Vertex *v = newVertex(p[X], p[Y], p[Z]);
	bounds.addPoint(p);
	return vertCount() - 1;
}
int Model::in_Face(int a, int b, int c)
{
	Vertex *v1 = vertices(a);
	Vertex *v2 = vertices(b);
	Vertex *v3 = vertices(c);

	Face *t = newFace(v1, v2, v3);

	return faceCount() - 1;
}
Vec3 Model::synthesizeNormal(Vertex *v)
{
	Vec3 n(0,0,0);
	int n_count = 0;

	edge_buffer& edges = v->edgeUses();
	for(int i=0; i<edges.length(); i++)
	{
		face_buffer& faces = edges(i)->faceUses();

		for(int j=0; j<faces.length(); j++)
		{
			n += faces(j)->plane().normal();
			n_count++;
		}
	}

	if( n_count )
		n /= (real)n_count;
	else
	{
		cerr << "Vertex with no normals!!: " << v->uniqID;
		cerr << " / " << v->tempID << endl;
	}
	return n;
}
Vertex *Model::newVertex(real x, real y, real z)
{
	Vertex *v = new Vertex(x, y, z);
	v->uniqID = vertices.add(v);
	validVertCount++;

	return v;
}
Edge *Model::newEdge(Vertex *a, Vertex *b)
{
	Edge *e = new Edge(a, b);

	e->uniqID = edges.add(e);
	e->sym()->uniqID = e->uniqID;

	validEdgeCount++;

	return e;
}
static Edge *get_edge(Model *m, Vertex *org, Vertex *v)
{
	edge_buffer& edge_uses = org->edgeUses();

	for(int i=0; i<edge_uses.length(); i++)
		if( edge_uses(i)->dest() == v )
			return edge_uses(i);

	Edge *e = m->newEdge(org, v);

	return e;
}
Face *Model::newFace(Vertex *v1, Vertex *v2, Vertex *v3)
{
	Edge *e0 = get_edge(this, v1, v2);   // v1->edgeTo(m, v2);
	Edge *e1 = get_edge(this, v2, v3);   // v2->edgeTo(m, v3);
	Edge *e2 = get_edge(this, v3, v1);   // v3->edgeTo(m, v1);
	Face *t = new Face(e0, e1, e2);
	t->uniqID = faces.add(t);
	validFaceCount++;

	return t;
}
void Model::killVertex(Vertex *v)
{
	if( v->isValid() )
	{
		v->kill();
		validVertCount--;
	}
}
void Model::killEdge(Edge *e)
{
	if( e->isValid() )
	{
		e->kill();
		validEdgeCount--;
	}
}
void Model::killFace(Face *f)
{
	if( f->isValid() )
	{
		f->kill();
		validFaceCount--;
	}
}
void Model::reshapeVertex(Vertex *v, real x, real y, real z)
{
	v->set(x, y, z);
}
void Model::remapVertex(Vertex *from, Vertex *to)
{
	from->remapTo(to);
}
void Model::maybeFixFace(Face *F)
{
	Vertex *v0=F->vertex(0); Vertex *v1=F->vertex(1); Vertex *v2=F->vertex(2);
	Edge *e0 = F->edge(0); Edge *e1 = F->edge(1); Edge *e2 = F->edge(2);

	bool a=(v0 == v1),  b=(v0 == v2),  c=(v1 == v2);

	if( a && c )
	{
		// This triangle has been reduced to a point
		killEdge(e0);
		killEdge(e1);
		killEdge(e2);

		killFace(F);
	}
	//
	// In the following 3 cases, the triangle has become an edge
	else if( a )
	{
		killEdge(e0);
		e1->remapTo(e2->sym());
		killFace(F);
	}
	else if( b )
	{
		killEdge(e2);
		e0->remapTo(e1->sym());
		killFace(F);
	}
	else if( c )
	{
		killEdge(e1);
		e0->remapTo(e2->sym());
		killFace(F);
	}
	else
	{
		// This triangle remains non-degenerate
	}
}
void Model::removeDegeneracy(face_buffer& changed)
{
	for(int i=0; i<changed.length(); i++)
		maybeFixFace(changed(i));
}
void Model::contractionRegion(Vertex *v1,const vert_buffer& vertices,face_buffer& changed)
{
	changed.reset();
	int i;
	untagFaceLoop(v1);
	for(i=0; i<vertices.length(); i++)
		untagFaceLoop(vertices(i));
	collectFaceLoop(v1, changed);
	for(i=0; i<vertices.length(); i++)
		collectFaceLoop(vertices(i), changed);
}
void Model::contractionRegion(Vertex *v1, Vertex *v2, face_buffer& changed)
{
	changed.reset();
	untagFaceLoop(v1);
	untagFaceLoop(v2);
	collectFaceLoop(v1, changed);
	collectFaceLoop(v2, changed);
}
void Model::contract(Vertex *v1, Vertex *v2, const Vec3& to,face_buffer& changed)
{
	contractionRegion(v1, v2, changed);
	reshapeVertex(v1, to[X], to[Y], to[Z]);
	v2->remapTo(v1);
	removeDegeneracy(changed);
}

Mat4 quadrix_vertex_constraint(const Vec3& v)
{
	Mat4 L(Mat4::identity);

	L(0,3) = -v[0];
	L(1,3) = -v[1];
	L(2,3) = -v[2];
	L(3,3) = v*v;

	L(3,0) = L(0,3);
	L(3,1) = L(1,3);
	L(3,2) = L(2,3);

	return L;
}
Mat4 quadrix_plane_constraint(real a, real b, real c, real d)
{
	Mat4 K(Mat4::zero);

	K(0,0) = a*a;   K(0,1) = a*b;   K(0,2) = a*c;  K(0,3) = a*d;
	K(1,0) =K(0,1); K(1,1) = b*b;   K(1,2) = b*c;  K(1,3) = b*d;
	K(2,0) =K(0,2); K(2,1) =K(1,2); K(2,2) = c*c;  K(2,3) = c*d;
	K(3,0) =K(0,3); K(3,1) =K(1,3); K(3,2) =K(2,3);K(3,3) = d*d;

	return K;
}
Mat4 quadrix_plane_constraint(const Vec3& n, real d)
{
	return quadrix_plane_constraint(n[X], n[Y], n[Z], d);
}
Mat4 quadrix_plane_constraint(Face& T)
{
	const Plane& p = T.plane();
	real a,b,c,d;
	p.coeffs(&a, &b, &c, &d);

	return quadrix_plane_constraint(a, b, c, d);
}
Mat4 quadrix_plane_constraint(const Vec3& v1, const Vec3& v2, const Vec3& v3)
{
	Plane P(v1,v2,v3);
	real a,b,c,d;
	P.coeffs(&a, &b, &c, &d);
	return quadrix_plane_constraint(a, b, c, d);
}
real quadrix_evaluate_vertex(const Vec3& v, const Mat4& K)
{
	real x=v[X], y=v[Y], z=v[Z];
	return x*x*K(0,0) + 2*x*y*K(0,1) + 2*x*z*K(0,2) + 2*x*K(0,3)
		+ y*y*K(1,1)   + 2*y*z*K(1,2) + 2*y*K(1,3)
		+ z*z*K(2,2)   + 2*z*K(2,3)
		+ K(3,3);
}
static bool is_border(Edge *e )
{
	return classifyEdge(e) == EDGE_BORDER;
}
bool check_for_discontinuity(Edge *e)
{
	return is_border(e);
}
Mat4 quadrix_discontinuity_constraint(Edge *edge, const Vec3& n)
{
	Vec3& org = *edge->org();
	Vec3& dest = *edge->dest();
	Vec3 e = dest - org;

	Vec3 n2 = e ^ n;
	unitize(n2);

	real d = -n2 * org;
	return quadrix_plane_constraint(n2, d);
}
Mat4 quadrix_discontinuity_constraint(Edge *edge)
{
	Mat4 D(Mat4::zero);

	face_buffer& faces = edge->faceUses();

	for(int i=0; i<faces.length(); i++)
		D += quadrix_discontinuity_constraint(edge,faces(i)->plane().normal());

	return D;
}
bool quadrix_find_local_fit(const Mat4& K,const Vec3& v1, const Vec3& v2,	Vec3& candidate)
{

	Vec3 v3 = (v1 + v2) / 2;

	bool try_midpoint = placement_policy > PLACE_ENDPOINTS;

	real c1 = quadrix_evaluate_vertex(v1, K);
	real c2 = quadrix_evaluate_vertex(v2, K);
	real c3;
	if( try_midpoint ) c3 = quadrix_evaluate_vertex(v3, K);

	if( c1<c2 )
	{
		if( try_midpoint && c3<c1 )
			candidate=v3;
		else
			candidate=v1;
	}
	else
	{
		if( try_midpoint && c3<c2 )
			candidate=v3;
		else
			candidate=v2;
	}

	return true;
}
bool quadrix_find_line_fit(const Mat4& Q,const Vec3& v1, const Vec3& v2,Vec3& candidate)
{
	Vec3 d = v1-v2;

	Vec3 Qv2 = Q*v2;
	Vec3 Qd  = Q*d;

	real denom = 2*d*Qd;

	if( denom == 0.0 )
		return false;

	real a = (d*Qv2 + v2*Qd) / denom;

	if( a<0.0 ) a=0.0;
	if( a>1.0 ) a=1.0;


	candidate = a*d + v2;
	return true;
}
bool quadrix_find_best_fit(const Mat4& Q, Vec3& candidate)
{
	Mat4 K = Q;
	K(3,0) = K(3,1) = K(3,2) = 0.0;  K(3,3) = 1;
	Mat4 M;
	real det = K.inverse(M);
	if( FEQ(det, 0.0, 1e-12) )
		return false;
	candidate[X] = M(0,3);
	candidate[Y] = M(1,3);
	candidate[Z] = M(2,3);
	return true;
}
real quadrix_pair_target(const Mat4& Q,Vertex *v1,Vertex *v2,Vec3& candidate)
{
	int policy = placement_policy;

	//
	// This analytic boundary preservation isn't really necessary.  The
	// boundary constraint quadrics are quite effective.  But, I've left it
	// in anyway.
	//
	if( will_preserve_boundaries )
	{
		int c1 = classifyVertex(v1);
		int c2 = classifyVertex(v2);

		if( c1 > c2 )
		{
			candidate = *v1;
			return quadrix_evaluate_vertex(candidate, Q);
		}
		else if( c2 > c1 )
		{
			candidate = *v2;
			return quadrix_evaluate_vertex(candidate, Q);
		}
		else if( c1>0 && policy>PLACE_LINE )
			policy = PLACE_LINE;

		//if( policy == PLACE_OPTIMAL ) assert(c1==0 && c2==0);
	}

	switch( policy )
	{
	case PLACE_OPTIMAL:
		if( quadrix_find_best_fit(Q, candidate) )
			break;

	case PLACE_LINE:
		if( quadrix_find_line_fit(Q, *v1, *v2, candidate) )
			break;

	default:
		quadrix_find_local_fit(Q, *v1, *v2, candidate);
		break;
	}

	return quadrix_evaluate_vertex(candidate, Q);
}


class pair_info : public Heapable
{
public:
	Vertex *v0, *v1;
	Vec3 candidate;
	real cost;
	pair_info(Vertex *a,Vertex *b) { v0=a; v1=b; cost=HUGE; }
	bool isValid() { return v0->isValid() && v1->isValid(); }
};
typedef buffer<pair_info *> pair_buffer;
class vert_info
{
public:
	pair_buffer pairs;
	Mat4 Q;
	real norm;
	vert_info() : Q(Mat4::zero) { pairs.init(2); norm=0.0; }
};
int will_draw_pairs = 0;
static Heap *heap;
static array<vert_info> vinfo;
static real proximity_limit;
static inline vert_info& vertex_info(Vertex *v)
{
	return vinfo(v->validID());
}
static bool check_for_pair(Vertex *v0, Vertex *v1)
{
	const pair_buffer& pairs = vertex_info(v0).pairs;

	for(int i=0; i<pairs.length(); i++)
	{
		if( pairs(i)->v0==v1 || pairs(i)->v1==v1 )
			return true;
	}

	return false;
}
static pair_info *new_pair(Vertex *v0, Vertex *v1)
{
	vert_info& v0_info = vertex_info(v0);
	vert_info& v1_info = vertex_info(v1);

	pair_info *pair = new pair_info(v0,v1);
	v0_info.pairs.add(pair);
	v1_info.pairs.add(pair);

	return pair;
}
static void delete_pair(pair_info *pair)
{
	vert_info& v0_info = vertex_info(pair->v0);
	vert_info& v1_info = vertex_info(pair->v1);
	v0_info.pairs.remove(v0_info.pairs.find(pair));
	v1_info.pairs.remove(v1_info.pairs.find(pair));
	if( pair->isInHeap() )
		heap->kill(pair->getHeapPos());
	delete pair;
}
static bool pair_is_valid(Vertex *u, Vertex *v)
{
	return norm2(*u - *v) < proximity_limit;
}
static int predict_face(Face& F, Vertex *v1, Vertex *v2, Vec3& vnew,Vec3& f1, Vec3& f2, Vec3& f3)
{
	int nmapped = 0;
	if( F.vertex(0) == v1 || F.vertex(0) == v2 )
	{ f1 = vnew;  nmapped++; }
	else f1 = *F.vertex(0);
	if( F.vertex(1) == v1 || F.vertex(1) == v2 )
	{ f2 = vnew;  nmapped++; }
	else f2 = *F.vertex(1);
	if( F.vertex(2) == v1 || F.vertex(2) == v2 )
	{ f3 = vnew;  nmapped++; }
	else f3 = *F.vertex(2);
	return nmapped;
}
#define MESH_INVERSION_PENALTY 1e9
static real pair_mesh_penalty(Model& M, Vertex *v1, Vertex *v2, Vec3& vnew)
{
	static face_buffer changed;
	changed.reset();
	M.contractionRegion(v1, v2, changed);
	real Nmin = 0;
	for(int i=0; i<changed.length(); i++)
	{
		Face& F = *changed(i);
		Vec3 f1, f2, f3;
		int nmapped = predict_face(F, v1, v2, vnew, f1, f2, f3);
		if( nmapped < 2 )
		{
			Plane Pnew(f1, f2, f3);
			real delta =  Pnew.normal() * F.plane().normal();

			if( Nmin > delta ) Nmin = delta;
		}
	}
	if( Nmin < 0.0 )
		return MESH_INVERSION_PENALTY;
	else
		return 0.0;
}
static void compute_pair_info(pair_info *pair)
{
	Vertex *v0 = pair->v0;
	Vertex *v1 = pair->v1;
	vert_info& v0_info = vertex_info(v0);
	vert_info& v1_info = vertex_info(v1);
	Mat4 Q = v0_info.Q + v1_info.Q;
	real norm = v0_info.norm + v1_info.norm;
	pair->cost = quadrix_pair_target(Q, v0, v1, pair->candidate);
	if( will_weight_by_area )
		pair->cost /= norm;
	if( will_preserve_mesh_quality )
		pair->cost += pair_mesh_penalty(M0, v0, v1, pair->candidate);
	if( pair->isInHeap() )
	{
		heap->update(pair, (float)-pair->cost);
	}
	else
	{
		heap->insert(pair, (float)-pair->cost);
	}
}
static void do_contract(Model& m, pair_info *pair)
{
	Vertex *v0 = pair->v0;  Vertex *v1 = pair->v1;
	vert_info& v0_info = vertex_info(v0);
	vert_info& v1_info = vertex_info(v1);
	Vec3 vnew = pair->candidate;
	int i;
	v0_info.Q += v1_info.Q;
	v0_info.norm += v1_info.norm;
	static face_buffer changed;
	changed.reset();
	m.contract(v0, v1, vnew, changed);
	delete_pair(pair);
	for(i=0; i<v0_info.pairs.length(); i++)
	{
		pair_info *p = v0_info.pairs(i);
		compute_pair_info(p);
	}
	static pair_buffer condemned(6); // collect condemned pairs for execution
	condemned.reset();
	for(i=0; i<v1_info.pairs.length(); i++)
	{
		pair_info *p = v1_info.pairs(i);

		Vertex *u;
		if( p->v0 == v1 )      u = p->v1;
		else if( p->v1 == v1)  u = p->v0;
		else cerr << "YOW!  This is a bogus pair." << endl;
		if( !check_for_pair(u, v0) )
		{
			p->v0 = v0;
			p->v1 = u;
			v0_info.pairs.add(p);
			compute_pair_info(p);
		}
		else
			condemned.add(p);
	}
	for(i=0; i<condemned.length(); i++)
		delete_pair(condemned(i));
	v1_info.pairs.reset(); // safety precaution
}
bool decimate_quadric(Vertex *v, Mat4& Q)
{
	if( vinfo.length() > 0 )
	{
		Q = vinfo(v->uniqID).Q;
		return true;
	}
	else
		return false;
}
void decimate_contract(Model& m)
{
	heap_node *top;
	pair_info *pair;

	for(;;)
	{
		top = heap->extract();
		if( !top ) return;
		pair = (pair_info *)top->obj;

		//
		// This may or may not be necessary.  I'm just not quite
		// willing to assume that all the junk has been removed from the
		// heap.
		if( pair->isValid() )
			break;

		delete_pair(pair);
	}

	do_contract(m, pair);

	//if( logfile && (selected_output&OUTPUT_COST) )
	//	*logfile << "#$cost " << m.validFaceCount << " "
	//	<< pair->cost << endl;

	M0.validVertCount--;  // Attempt to maintain valid vertex information
}
real decimate_error(Vertex *v)
{
	vert_info& info = vertex_info(v);

	real err = quadrix_evaluate_vertex(*v, info.Q);

	if( will_weight_by_area )
		err /= info.norm;

	return err;
}
real decimate_min_error()
{
	heap_node *top;
	pair_info *pair;

	for(;;)
	{
		top = heap->top();
		if( !top ) return -1.0;
		pair = (pair_info *)top->obj;

		if( pair->isValid() )
			break;

		top = heap->extract();
		delete_pair(pair);
	}

	return pair->cost;
}
real decimate_max_error(Model& m)
{
	real max_err = 0;
	for(int i=0; i<m.vertCount(); i++)
		if( m.vertex(i)->isValid() )
		{
			max_err = MAX(max_err, decimate_error(m.vertex(i)));
		}
		return max_err;
}
void decimate_init(Model& m, real limit)
{
	int i,j;
	vinfo.init(m.vertCount());
	cout << "  Decimate:  Distributing shape constraints." << endl;
	if( will_use_vertex_constraint )
		for(i=0; i<m.vertCount(); i++)
		{
			Vertex *v = m.vertex(i);
			if( v->isValid() )
				vertex_info(v).Q = quadrix_vertex_constraint(*v);
		}

		for(i=0; i<m.faceCount(); i++)
			if( m.face(i)->isValid() )
			{
				if( will_use_plane_constraint )
				{
					Mat4 Q = quadrix_plane_constraint(*m.face(i));
					real norm = 0.0;

					if( will_weight_by_area )
					{
						norm = m.face(i)->area();
						Q *= norm;
					}

					for(j=0; j<3; j++)
					{
						vertex_info(m.face(i)->vertex(j)).Q += Q;
						vertex_info(m.face(i)->vertex(j)).norm += norm;

					}
				}
			}

			if( will_constrain_boundaries )
			{
				cout << "  Decimate:  Accumulating discontinuity constraints." << endl;
				for(i=0; i<m.edgeCount(); i++)
					if( m.edge(i)->isValid() && check_for_discontinuity(m.edge(i)) )
					{
						Mat4 B = quadrix_discontinuity_constraint(m.edge(i));
						real norm = 0.0;

						if( will_weight_by_area )
						{
							norm = norm2(*m.edge(i)->org() - *m.edge(i)->dest());
							B *= norm;
						}

						B *= boundary_constraint_weight;
						vertex_info(m.edge(i)->org()).Q += B;
						vertex_info(m.edge(i)->org()).norm += norm;
						vertex_info(m.edge(i)->dest()).Q += B;
						vertex_info(m.edge(i)->dest()).norm += norm;
					}
			}

			cout << "  Decimate:  Allocating heap." << endl;
			heap = new Heap(m.validEdgeCount);

			int pair_count = 0;

			cout << "  Decimate:  Collecting pairs [edges]." << endl;
			for(i=0; i<m.edgeCount(); i++)
				if( m.edge(i)->isValid() )
				{
					pair_info *pair = new_pair(m.edge(i)->org(), m.edge(i)->dest());
					compute_pair_info(pair);
					pair_count++;
				}

				if( limit<0 )
				{
					limit = m.bounds.radius * 0.05;
					cout << "  Decimate:  Auto-limiting at 5% of model radius." << endl;
				}
				proximity_limit = limit * limit;
				if( proximity_limit > 0 )
				{
					cout << "  Decimate:  Collecting pairs [limit="<<limit<<"]." << endl;
					ProxGrid grid(m.bounds.min, m.bounds.max, limit);
					for(i=0; i<m.vertCount(); i++)
						grid.addPoint(m.vertex(i));

					buffer<Vec3 *> nearby(32);
					for(i=0; i<m.vertCount(); i++)
					{
						nearby.reset();
						grid.proximalPoints(m.vertex(i), nearby);

						for(j=0; j<nearby.length(); j++)
						{
							Vertex *v1 = m.vertex(i);
							Vertex *v2 = (Vertex *)nearby(j);

							if( v1->isValid() && v2->isValid() )
							{
								if( !check_for_pair(v1,v2) )
								{
									pair_info *pair = new_pair(v1,v2);
									compute_pair_info(pair);
									pair_count++;
								}
							}

						}
					}
				}
				else
					cout << "  Decimate:  Ignoring non-edge pairs [limit=0]." << endl;

				cout << "  Decimate:  Designated " << pair_count << " pairs." << endl;
}