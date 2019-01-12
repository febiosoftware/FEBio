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

public:
	bool Create(FEMesh* mesh);
};
