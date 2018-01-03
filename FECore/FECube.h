#pragma once
#include "FEMesh.h"

//-----------------------------------------------------------------------------
// This class tries to identify surfaces, edges, and corner nodes on a mesh, 
// assuming that it is a cube.
// The surfaces are ordered as follows:
// 1: +X, 2: -X, 3: +Y, 4: -Y, 5: +Z, 6: -Z
class FECube
{
public:
	// constructor
	FECube();

	// destructor
	~FECube();

	// build the cube data
	bool Build(FEMesh* mesh);

	// get the mesh of this cube
	FEMesh* GetMesh();

public:
	// Get a surface
	FESurface* GetSurface(int i);

	// get the node set of the corner nodes
	const FENodeSet& GetCornerNodes() const;

	// get the node set of boundary nodes
	const FENodeSet& GetBoundaryNodes() const;

private:
	FEMesh*	m_mesh;

	FESurface*	m_surf[6];		// the six boundary surfaces
	FENodeSet	m_corners;		// the eight corner nodes
	FENodeSet	m_boundary;		// set of all boundary nodes
};
