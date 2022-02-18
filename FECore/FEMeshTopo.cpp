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



#include "stdafx.h"
#include "FEMeshTopo.h"
#include "FEElementList.h"
#include "FEMesh.h"
#include "FEDomain.h"
#include "FEElemElemList.h"
#include "FESurface.h"

class FEMeshTopo::MeshTopoImp
{
public:
	MeshTopoImp()
	{
		m_minId = -1;
	}

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

	std::vector<FEElement*>	m_elem;		// element list

	std::vector<int>	m_lut;
	int					m_minId;
};

FEMeshTopo::FEMeshTopo() : imp(new FEMeshTopo::MeshTopoImp)
{
}

FEMeshTopo::~FEMeshTopo()
{
	delete imp;
}

// get the mesh
FEMesh* FEMeshTopo::GetMesh()
{
	return imp->m_mesh;
}

bool FEMeshTopo::Create(FEMesh* mesh)
{
	imp->m_mesh = mesh;
	FEElementList elemList(*mesh);

	// create a vector of all elements
	int NEL = mesh->Elements();
	imp->m_elem.resize(NEL);
	NEL = 0;
	for (int i = 0; i < mesh->Domains(); ++i)
	{
		FEDomain& dom = mesh->Domain(i);
		int nel = dom.Elements();
		for (int j = 0; j < nel; ++j)
		{
			imp->m_elem[NEL++] = &dom.ElementRef(j);
		}
	}

	// build the index lookup table
	FEElementIterator it(mesh);
	imp->m_minId = -1;
	int minId = -1, maxId = -1;
	for (; it.isValid(); ++it)
	{
		int nid = (*it).GetID(); assert(nid != -1);
		if (minId == -1)
		{
			minId = nid;
			maxId = nid;
		}
		else
		{
			if (nid < minId) minId = nid;
			if (nid > maxId) maxId = nid;
		}
	}
	int lutSize = maxId - minId + 1;
	imp->m_lut.assign(lutSize, -1);
	imp->m_minId = minId;
	int n = 0;
	for (it.reset(); it.isValid(); ++it, ++n)
	{
		int nid = (*it).GetID() - minId;
		imp->m_lut[nid] = n;
	}

	// create the element neighbor list
	if (imp->m_ENL.Create(mesh) == false) return false;

	// create a face list
	if (imp->m_faceList.Create(*mesh, imp->m_ENL) == false) return false;

	// extract the surface facets
	imp->m_surface = imp->m_faceList.GetSurface();

	// create the element-face list
	if (imp->m_EFL.Create(elemList, imp->m_faceList) == false) return false;

	// create the element-surface facet list
	if (imp->m_ESL.Create(elemList, imp->m_surface) == false) return false;
	imp->m_surface.BuildNeighbors();

	// create the edge list (from the face list)
	if (imp->m_edgeList.Create(mesh) == false) return false;

	// create the element-edge list
	if (imp->m_EEL.Create(elemList, imp->m_edgeList) == false) return false;

	// create the face-edge list
	if (imp->m_FEL.Create(imp->m_faceList, imp->m_edgeList) == false) return false;

	return true;
}

// return elements
int FEMeshTopo::Elements()
{
	return (int)imp->m_elem.size();
}

// return an element
FEElement* FEMeshTopo::Element(int i)
{
	if (i < 0) return nullptr;
	else return imp->m_elem[i];
}

int FEMeshTopo::GetElementIndexFromID(int elemId)
{
	int eid = elemId - imp->m_minId;
	if ((eid < 0) || (eid >= imp->m_lut.size())) { assert(false); return -1; }
	int lid = imp->m_lut[eid]; 
	assert(lid != -1);
	return lid;
}

int FEMeshTopo::Faces()
{
	return imp->m_faceList.Faces();
}

// return the number of surface faces
int FEMeshTopo::SurfaceFaces() const
{
	return imp->m_surface.Faces();
}

// return a face
const FEFaceList::FACE& FEMeshTopo::Face(int i)
{
	return imp->m_faceList[i];
}

// return a surface facet
const FEFaceList::FACE& FEMeshTopo::SurfaceFace(int i) const
{
	return imp->m_surface[i];
}

// return the element-face list
const std::vector<int>& FEMeshTopo::ElementFaceList(int nelem)
{
	return imp->m_EFL.FaceList(nelem);
}

// return the number of edges in the mesh
int FEMeshTopo::Edges()
{
	return imp->m_edgeList.Edges();
}

// return a edge
const FEEdgeList::EDGE& FEMeshTopo::Edge(int i)
{
	return imp->m_edgeList[i];
}

// return the face-edge list
const std::vector<int>& FEMeshTopo::FaceEdgeList(int nface)
{
	return imp->m_FEL.EdgeList(nface);
}

// return the element-edge list
const std::vector<int>& FEMeshTopo::ElementEdgeList(int nelem)
{
	return imp->m_EEL.EdgeList(nelem);
}

// return the list of face indices of a surface
std::vector<int> FEMeshTopo::FaceIndexList(FESurface& s)
{
	FENodeFaceList NFL;
	NFL.Create(imp->m_faceList);

	int NF = s.Elements();
	std::vector<int> fil(NF, -1);
	for (int i = 0; i < NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);

		int nval = NFL.Faces(el.m_node[0]);
		const std::vector<int>& nfl = NFL.FaceList(el.m_node[0]);
		for (int j = 0; j < nval; ++j)
		{
			const FEFaceList::FACE& face = imp->m_faceList[nfl[j]];
			if (face.IsEqual(&el.m_node[0]))
			{
				fil[i] = nfl[j];
				break;
			}
		}
		assert(fil[i] != -1);
	}

	return fil;
}

// return the list of face indices of a surface
std::vector<int> FEMeshTopo::SurfaceFaceIndexList(FESurface& s)
{
	FENodeFaceList NFL;
	NFL.Create(imp->m_surface);

	int NF = s.Elements();
	std::vector<int> fil(NF, -1);
	for (int i = 0; i < NF; ++i)
	{
		FESurfaceElement& el = s.Element(i);

		int nval = NFL.Faces(el.m_node[0]);
		const std::vector<int>& nfl = NFL.FaceList(el.m_node[0]);
		for (int j = 0; j < nval; ++j)
		{
			const FEFaceList::FACE& face = imp->m_surface[nfl[j]];
			if (face.IsEqual(&el.m_node[0]))
			{
				fil[i] = nfl[j];
				break;
			}
		}
		assert(fil[i] != -1);
	}

	return fil;
}

// return the element neighbor list
std::vector<FEElement*> FEMeshTopo::ElementNeighborList(int n)
{
	FEElement* el = Element(n);
	int nbrs = 0;
	switch (el->Shape())
	{
	case ET_HEX8: nbrs = 8; break;
	case ET_TET4: nbrs = 4; break;
	case ET_TET5: nbrs = 4; break;
	default:
		assert(false);
	}

	vector<FEElement*> elemList;
	for (int i = 0; i < nbrs; ++i)
	{
		elemList.push_back(imp->m_ENL.Neighbor(n, i));
	}

	return elemList;
}

// return the element neighbor list
std::vector<int> FEMeshTopo::ElementNeighborIndexList(int n)
{
	FEElement* el = Element(n);
	int nbrs = 0;
	switch (el->Shape())
	{
	case ET_HEX8: nbrs = 8; break;
	case ET_TET4: nbrs = 4; break;
	case ET_TET5: nbrs = 4; break;
	default:
		assert(false);
	}

	vector<int> elemList;
	for (int i = 0; i < nbrs; ++i)
	{
		elemList.push_back(imp->m_ENL.NeighborIndex(n, i));
	}

	return elemList;
}
