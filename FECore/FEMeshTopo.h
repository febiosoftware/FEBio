#pragma once
#include "FEEdgeList.h"
#include "FEFaceList.h"
#include "FEElemElemList.h"

class FEMesh;

class FECORE_API FEMeshTopo
{
public:
	FEMesh*				m_mesh;			// the mesh
	FEEdgeList			m_edgeList;		// the edge list
	FEElementEdgeList	m_EEL;			// the element-edge list
	FEFaceList			m_faceList;		// the face list (all faces)
	FEElementFaceList	m_EFL;			// the element-face list
	FEElemElemList		m_ENL;			// the element neighbor list
	FEFaceList			m_surface;		// only surface facets
	FEElementFaceList	m_ESL;			// element-surface facet list
	FEFaceEdgeList		m_FEL;			// face-edge list

public:
	bool Create(FEMesh* mesh);

	// return the number of faces in the mesh
	int Faces();

	// return a face
	const FEFaceList::FACE& Face(int i);

	// return the number of edges in the mesh
	int Edges();

	// return an edge
	const FEEdgeList::EDGE& Edge(int i);
};
