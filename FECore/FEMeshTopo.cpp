#include "stdafx.h"
#include "FEMeshTopo.h"
#include "FEElementList.h"

bool FEMeshTopo::Create(FEMesh* mesh)
{
	m_mesh = mesh;
	FEElementList elemList(*mesh);

	// create the element neighbor list
	if (m_ENL.Create(mesh) == false) return false;

	// create a face list
	if (m_faceList.Create(*mesh, m_ENL) == false) return false;

	// extract the surface facets
	m_surface = m_faceList.GetSurface();

	// create the element-face list
	if (m_EFL.Create(elemList, m_faceList) == false) return false;

	// create the element-surface facet list
	if (m_ESL.Create(elemList, m_surface) == false) return false;
	m_surface.BuildNeighbors();

	// create the edge list (from the face list)
	if (m_edgeList.Create(mesh) == false) return false;

	// create the element-edge list
	if (m_EEL.Create(elemList, m_edgeList) == false) return false;

	// create the face-edge list
	if (m_FEL.Create(m_faceList, m_edgeList) == false) return false;

	return true;
}

int FEMeshTopo::Faces()
{
	return m_faceList.Faces();
}

// return a face
const FEFaceList::FACE& FEMeshTopo::Face(int i)
{
	return m_faceList[i];
}

// return the number of edges in the mesh
int FEMeshTopo::Edges()
{
	return m_edgeList.Edges();
}

// return a edge
const FEEdgeList::EDGE& FEMeshTopo::Edge(int i)
{
	return m_edgeList[i];
}
