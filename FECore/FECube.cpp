#include "stdafx.h"
#include "FECube.h"
#include "FESurface.h"

FECube::FECube() : m_mesh(0)
{
}

FECube::~FECube()
{
	for (int i = 0; i<6; ++i) 
	{
		delete m_surf[i];
		m_surf[i] = 0;
	}
}

// Get a surface
FESurface* FECube::GetSurface(int i)
{
	return m_surf[i];
}

// get the node set of the corner nodes
const FENodeSet& FECube::GetCornerNodes() const
{
	return m_corners;
}

// get the node set of boundary nodes
const FENodeSet& FECube::GetBoundaryNodes() const
{
	return m_boundary;
}

// get the mesh of this cube
FEMesh* FECube::GetMesh()
{
	return m_mesh;
}

bool FECube::Build(FEMesh* mesh)
{
	// make sure we have a mesh
	m_mesh = mesh;
	if (mesh == 0) return false;

	FEMesh& m = *mesh;
	int NN = m.Nodes();
	
	// first, get the outside surface
	FESurface* boundary = m.ElementBoundarySurface();

	// get the boundary node set
	m_boundary = boundary->GetNodeSet();

	// Next, split it up in 6 surfaces
	// We divide the surface by comparing normals to the 6 surface normals of a cube
	vec3d fn[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	for (int n = 0; n<6; ++n)
	{
		// create the surface
		m_surf[n] = new FESurface(m_mesh->GetFEModel());

		// get the normal for this face
		vec3d N = fn[n];

		int faces = 0;
		for (int i = 0; i<boundary->Elements(); ++i)
		{
			FESurfaceElement& face = boundary->Element(i);
			vec3d Ni = boundary->SurfaceNormal(face, 0, 0);
			if (Ni*N > 0.9999) faces++;
		}
		m_surf[n]->Create(faces);

		faces = 0;
		for (int i = 0; i<boundary->Elements(); ++i)
		{
			FESurfaceElement& face = boundary->Element(i);
			vec3d Ni = boundary->SurfaceNormal(face, 0, 0);
			if (Ni*N > 0.9999)
			{
				FESurfaceElement& newFace = m_surf[n]->Element(faces++);
				newFace = face;
			}
		}

		m_surf[n]->Init();
	}

	// we also need to find the 8 corner nodes
	vector<int> tag(NN, 0);
	for (int n = 0; n<6; ++n)
	{
		FESurface& sn = *m_surf[n];
		FENodeSet ns = sn.GetNodeSet();
		for (int i = 0; i<ns.size(); ++i) tag[ns[i]]++;
	}

	for (int i = 0; i<NN; ++i) if (tag[i] == 3) m_corners.add(i);
	assert(m_corners.size() == 8);

	// don't forget to cleanup
	delete boundary;

	return true;
}
