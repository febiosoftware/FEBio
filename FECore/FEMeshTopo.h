#pragma once
#include "FEEdgeList.h"
#include "FEFaceList.h"

class FEMesh;
class FEElement;

//-----------------------------------------------------------------------------
//! This is a helper class for looping over elements, faces (internal and external), 
//! and edges of a FEMesh.
class FECORE_API FEMeshTopo
{
	class MeshTopoImp;

public:
	FEMeshTopo();
	~FEMeshTopo();

	// Create the FEMeshTopo from a mesh
	bool Create(FEMesh* mesh);

	// return number of elements
	int Elements();

	// return an element
	FEElement* Element(int i);

	// return the number of faces in the mesh
	int Faces();

	// return a face
	const FEFaceList::FACE& Face(int i);

	// return the element-face list
	const std::vector<int>& ElementFaceList(int nelem);

	// return the number of edges in the mesh
	int Edges();

	// return an edge
	const FEEdgeList::EDGE& Edge(int i);

	// return the face-edge list
	const std::vector<int>& FaceEdgeList(int nface);

	// return the element-edge list
	const std::vector<int>& ElementEdgeList(int nelem);

private:
	MeshTopoImp*	imp;
};
