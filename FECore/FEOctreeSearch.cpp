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
#include "FEOctreeSearch.h"
#include "FEMesh.h"
#include "FEElementList.h"
#include "FESolidDomain.h"

//-----------------------------------------------------------------------------
FEOctreeSearch::Block::Block(FEMesh* mesh, int level)
{
	m_parent = nullptr;
	m_mesh = mesh;
	m_level = level;
}

//-----------------------------------------------------------------------------
FEOctreeSearch::Block::Block(FEOctreeSearch::Block* parent, int level)
{
	m_parent = parent;
	m_mesh = parent->m_mesh;
	m_level = level;
}

//-----------------------------------------------------------------------------
FEOctreeSearch::Block::~Block()
{ 
	Clear(); 
}

//-----------------------------------------------------------------------------
void FEOctreeSearch::Block::Clear()
{
	for (size_t i = 0; i < m_children.size(); ++i) delete m_children[i];
	m_children.clear();
}

//-----------------------------------------------------------------------------
// Create the eight children of an octree node and find their contents

void FEOctreeSearch::Block::CreateChildren(const int max_level, const int max_elem)
{
	vec3d dc = (m_cmax - m_cmin) / 2;
	for (int i = 0; i <= 1; ++i) {
		for (int j = 0; j <= 1; ++j) {
			for (int k = 0; k <= 1; ++k) {
				Block* node = new Block(this, m_level + 1);

				// evaluate bounding box by subdividing parent node box
				node->m_cmin = vec3d(m_cmin.x + i*dc.x,
					m_cmin.y + j*dc.y,
					m_cmin.z + k*dc.z);
				node->m_cmax = vec3d(m_cmax.x - (1 - i)*dc.x,
					m_cmax.y - (1 - j)*dc.y,
					m_cmax.z - (1 - k)*dc.z);

				// find all elements in this child node
				node->FillBlock();

				if (node->m_selist.size()) {
					// use recursion to create children of this node
					// as long as node contains too many elements
					// and max octree levels not exceeded
					if ((node->m_level < max_level) &&
						(node->m_selist.size() > max_elem))
						node->CreateChildren(max_level, max_elem);
				}

				// store this node
				m_children.push_back(node);
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Find all elements that fall inside a block

void FEOctreeSearch::Block::FillBlock()
{
	// Loop over all elements in the parent node
	for (int i = 0; i<(int)m_parent->m_selist.size(); ++i) {
		FEElement* pe = m_parent->m_selist[i];
		if (ElementIntersectsNode(pe)) {
			// add this surface element to the current node
			m_selist.push_back(pe);
		}
	}
}

//-----------------------------------------------------------------------------
// Determine whether a surface element intersects a node

bool FEOctreeSearch::Block::ElementIntersectsNode(FEElement* pel)
{
	// Extract FE node coordinates from element
	// and determine bounding box of element
	FEMesh& mesh = *m_mesh;
	FEElement& el = *pel;
	int N = el.Nodes();
	vector<vec3d> fenode(N);
	fenode[0] = mesh.Node(el.m_node[0]).m_r0;
	vec3d fmin = fenode[0];
	vec3d fmax = fenode[0];
	for (int i = 1; i<N; ++i) {
		fenode[i] = mesh.Node(el.m_node[i]).m_r0;
		if (fenode[i].x < fmin.x) fmin.x = fenode[i].x;
		if (fenode[i].x > fmax.x) fmax.x = fenode[i].x;
		if (fenode[i].y < fmin.y) fmin.y = fenode[i].y;
		if (fenode[i].y > fmax.y) fmax.y = fenode[i].y;
		if (fenode[i].z < fmin.z) fmin.z = fenode[i].z;
		if (fenode[i].z > fmax.z) fmax.z = fenode[i].z;
	}

	// Check if bounding boxes of OT node and surface element overlap
	if ((fmax.x < m_cmin.x) || (fmin.x > m_cmax.x)) return false;
	if ((fmax.y < m_cmin.y) || (fmin.y > m_cmax.y)) return false;
	if ((fmax.z < m_cmin.z) || (fmin.z > m_cmax.z)) return false;

	// At this point we find that bounding boxes overlap.
	// Technically that does not prove that the surface element is
	// inside the octree node, but any additional check would be
	// more expensive.

	return true;
}

//-----------------------------------------------------------------------------
bool FEOctreeSearch::Block::IsInside(const vec3d& r) const
{
	if ((r.x < m_cmin.x) || (r.x > m_cmax.x)) return false;
	if ((r.y < m_cmin.y) || (r.y > m_cmax.y)) return false;
	if ((r.z < m_cmin.z) || (r.z > m_cmax.z)) return false;
	return true;
}

//-----------------------------------------------------------------------------
FEElement* FEOctreeSearch::Block::FindElement(const vec3d& y, double r[3])
{
	if (IsInside(y))
	{
		if (m_children.empty() == false)
		{
			for (int i = 0; i < m_children.size(); ++i)
			{
				FEElement* pe = m_children[i]->FindElement(y, r);
				if (pe) return pe;
			}
		}
		else
		{
			vec3d x[FEElement::MAX_NODES];
			int NE = (int)m_selist.size();
			for (int i = 0; i<NE; ++i)
			{
				// get the next element
				FESolidElement& e = *((FESolidElement*)m_selist[i]);

				// get the element nodal coordinates
				int neln = e.Nodes();
				for (int j = 0; j<neln; ++j) x[j] = m_mesh->Node(e.m_node[j]).m_r0;

				// first, as a quick check, we see if y lies in the bounding box defined by x
				FEBoundingBox box(x[0]);
				for (int j = 1; j<neln; ++j) box.add(x[j]);

				// inflate a little for round-off
				double dr = box.radius()*1e-6;
				box.inflate(dr, dr, dr);

				if (box.IsInside(y))
				{
					FESolidDomain* dom = dynamic_cast<FESolidDomain*>(e.GetMeshPartition());

					// If the point y lies inside the box, we apply a Newton method to find
					// the isoparametric coordinates r
					if (dom->ProjectToReferenceElement(e, y, r)) return &e;
				}
			}
		}
	}
	
	return nullptr;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEOctreeSearch::FEOctreeSearch(FEMesh* mesh)
{
	m_root = nullptr;
	m_mesh = mesh;
	m_dom = nullptr;
	m_max_level = 6;
	m_max_elem = 9;
}

FEOctreeSearch::FEOctreeSearch(FEDomain* domain)
{
	m_root = nullptr;
	m_mesh = domain->GetMesh();
	m_dom = domain;
	m_max_level = 6;
	m_max_elem = 9;
}

//-----------------------------------------------------------------------------
FEOctreeSearch::~FEOctreeSearch()
{

}

//-----------------------------------------------------------------------------
bool FEOctreeSearch::Init(double inflate)
{
	if (m_mesh == nullptr) return false;

	// Set up the root node in the octree
	if (m_root) delete m_root;
	m_root = new Block(m_mesh, 0);

	Block& root = *m_root;

	// Create the list of all elements in the root node
	int nel = 0;
	if (m_dom == nullptr)
	{
		FEElementList EL(*m_mesh);
		nel = m_mesh->Elements();
		for (FEElementList::iterator it = EL.begin(); it != EL.end(); ++it)
			root.m_selist.push_back(it);
	}
	else
	{
		nel = m_dom->Elements();
		for (int i = 0; i < nel; ++i)
		{
			root.m_selist.push_back(&m_dom->ElementRef(i));
		}
	}

	// Find the bounding box of the mesh
	vec3d r0 = (m_mesh->Node(0)).m_r0;
	root.m_cmin = r0;
	root.m_cmax = r0;
	for (int i = 1; i<m_mesh->Nodes(); ++i) {
		vec3d r = (m_mesh->Node(i)).m_r0;
		if (r.x < root.m_cmin.x) root.m_cmin.x = r.x;
		if (r.x > root.m_cmax.x) root.m_cmax.x = r.x;
		if (r.y < root.m_cmin.y) root.m_cmin.y = r.y;
		if (r.y > root.m_cmax.y) root.m_cmax.y = r.y;
		if (r.z < root.m_cmin.z) root.m_cmin.z = r.z;
		if (r.z > root.m_cmax.z) root.m_cmax.z = r.z;
	}

	// expand bounding box by search tolerance stol
	double d = (root.m_cmax - root.m_cmin).norm()*inflate;
	root.m_cmin -= vec3d(d, d, d);
	root.m_cmax += vec3d(d, d, d);

	// figure out the levels
	int level = (int) (log((double)nel) / log(8.0));
	if (level < 1) level = 1;
	if (level > m_max_level) level = m_max_level - 1;

	// Recursively create children of this root
	if (root.m_selist.size()) {
		if ((root.m_level < m_max_level) &&
			(root.m_selist.size() > m_max_elem))
			root.CreateChildren(level, m_max_elem);
	}

	return true;
}

//-----------------------------------------------------------------------------
FEElement* FEOctreeSearch::FindElement(const vec3d& x, double r[3])
{
	return m_root->FindElement(x, r);
}
