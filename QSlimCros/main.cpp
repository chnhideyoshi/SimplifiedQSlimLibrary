#include "Mesh.h"
#include "QSlim.h"

static void qslim_init()
{
    int i;
    cerr << "Reading input ..." << endl;
    cerr << "Cleaning up initial input ..." << endl;
	int initialVertCount = M0.vertCount();
	int initialEdgeCount = M0.edgeCount();
	int initialFaceCount = M0.faceCount();
    for(i=0; i<M0.faceCount(); i++)
	if( !M0.face(i)->plane().isValid() )
	    M0.killFace(M0.face(i));
    M0.removeDegeneracy(M0.allFaces());
    for(i=0; i<M0.vertCount(); i++)
    {
	if( M0.vertex(i)->edgeUses().length() == 0 )
	    M0.vertex(i)->kill();
    }
    cerr << "Input model summary:" << endl;
    cerr << "    Vertices    : " << initialVertCount << endl;
    cerr << "    Edges       : " << initialEdgeCount << endl;
    int man=0, non=0, bndry=0, bogus=0;
    for(i=0; i<M0.edgeCount(); i++)
        switch( M0.edge(i)->faceUses().length() )
        {
        case 0:
            bogus++;
            break;
        case 1:
            bndry++;
            break;
        case 2:
            man++;
            break;
        default:
            non++;
            break;
        }
    if( bogus )
        cerr << "        Bogus       : " << bogus << endl;
    cerr << "        Boundary    : " << bndry << endl;
    cerr << "        Manifold    : " << man << endl;
    cerr << "        Nonmanifold : " << non << endl;

    cerr << "    Faces       : " << initialFaceCount << endl;
}
static void qslim_run()
{
    decimate_init(M0, pair_selection_tolerance);
    while( M0.validFaceCount > face_target&& decimate_min_error() < error_tolerance )
		decimate_contract(M0);
}
static void InitM0(Mesh& m)
{
	for(size_t i=0;i<m.Vertices.size();i++)
	{
		Point3d p=m.Vertices[i];
		Vec3 v(p.X,p.Y,p.Z);
		M0.in_Vertex(v);
	}
	for(size_t i=0;i<m.Faces.size();i++)
	{
		Triangle t=m.Faces[i];
		M0.in_Face(t.P0Index,t.P1Index,t.P2Index);
	}
}
static void ReplaceM(Mesh& m)
{
	m.Vertices.swap(std::vector<Point3d>());
	m.Faces.swap(std::vector<Triangle>());
	m.Vertices.reserve(M0.vertCount());
	m.Faces.reserve(M0.faceCount());
	int* map=new int[M0.vertCount()];
	for(int i=0;i<M0.vertCount();i++)
		map[i]=-1;
	for(int i=0;i<M0.vertCount();i++)
	{
		if(M0.vertex(i)->isValid())
		{
			real* data=M0.vertex(i)->raw();
			Point3d p((float)data[0],(float)data[1],(float)data[2]);
			map[i]=m.AddVertex(p);
		}
	}
	for(int i=0;i<M0.faceCount();i++)
	{
		if(M0.face(i)->isValid())
		{
			Vertex* v0= M0.face(i)->vertex(0);
			Vertex* v1= M0.face(i)->vertex(1);
			Vertex* v2= M0.face(i)->vertex(2);
			Triangle t(map[v0->uniqID],map[v1->uniqID],map[v2->uniqID]);
			m.AddFace(t);
		}
	}
	delete[] map;
}

int main()
{
	Mesh m;
	PlyManager::ReadFile(m,"D:\\VTKproj\\sample.ply");
	InitM0(m);
	qslim_init();
	face_target = 2*m.Faces.size()/3;
	error_tolerance = HUGE;
	will_use_plane_constraint = true;
	will_use_vertex_constraint = false;
	will_preserve_boundaries = true;
	will_preserve_mesh_quality = true;
	will_constrain_boundaries = true;
	boundary_constraint_weight = 1.0;
	will_weight_by_area = false;
	placement_policy = 1;
	pair_selection_tolerance = 0.0;
	qslim_run();
	ReplaceM(m);
	PlyManager::Output(m,"D:\\VTKproj\\sample_deci.ply");

	return 0;
}