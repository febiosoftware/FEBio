/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "FEEdgeList.h"
#include "FEFaceList.h"

class FEMesh;
class FEElement;
class FESurface;

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

	// get the mesh
	FEMesh* GetMesh();

	// return number of elements
	int Elements();

	// return an element
	FEElement* Element(int i);

	// get the element index (into global element array) from an element ID
	int GetElementIndexFromID(int elemId);

	// return the number of faces in the mesh
	int Faces();

	// return a face
	const FEFaceList::FACE& Face(int i);

	// return the number of surface faces
	int SurfaceFaces() const;

	// return the number of surface faces
	const FEFaceList::FACE& SurfaceFace(int i) const;

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

	// return the list of face indices of a surface
	std::vector<int> FaceIndexList(FESurface& s);

	// return the list of face indices of a surface
	std::vector<int> SurfaceFaceIndexList(FESurface& s);

	// return the element neighbor list
	std::vector<FEElement*> ElementNeighborList(int i);

	// return the element neighbor index list
	std::vector<int> ElementNeighborIndexList(int i);

private:
	MeshTopoImp*	imp;
};
