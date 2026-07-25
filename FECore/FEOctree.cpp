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
#include "FEOctree.h"
#include "FESurface.h"
#include "FEMesh.h"
#include <stack>

OTnode::OTnode()
{
	m_po = nullptr;
	level = 0;
}

OTnode::~OTnode()
{
}

void OTnode::Clear()
{ 
	children.clear(); 
}

//-----------------------------------------------------------------------------
// Create the eight children of an octree node and find their contents

void OTnode::CreateChildren(const int max_level, const int max_elem)
{
	children.reserve(8);
	vec3d dc = (cmax - cmin)/2;
	for (int i=0; i<=1; ++i) {
		for (int j=0; j<=1; ++j) {
			for (int k=0; k<=1; ++k) {
				OTnode node;
				node.m_po = m_po;
				// evaluate bounding box by subdividing parent node box
				node.cmin = vec3d(cmin.x+i*dc.x,
								  cmin.y+j*dc.y,
								  cmin.z+k*dc.z);
				node.cmax = vec3d(cmax.x-(1-i)*dc.x,
								  cmax.y-(1-j)*dc.y,
								  cmax.z-(1-k)*dc.z);
				// update octree level
				node.level = level + 1;
				// find all surface elements in this child node
				node.FillNode(selist);
				if (node.selist.size()) {
					// use recursion to create children of this node
					// as long as node contains too many elements
					// and max octree levels not exceeded
					if ((node.level < max_level) &&
						(node.selist.size() > max_elem))
						node.CreateChildren(max_level, max_elem);
				}
				// store this node
				children.push_back(node);
			}
		}
	}
}

// determine whether two boxes intersect (cheap version)
inline bool CheckBoxIntersection(const FEOctree::Box& src, const FEOctree::Box& dst)
{
	// get element's bounding box
	vec3d fmin = src.r0;
	vec3d fmax = src.r1;

	// Check if bounding boxes of OT node and surface element overlap
	if ((fmax.x < dst.r0.x) || (fmin.x > dst.r1.x)) return false;
	if ((fmax.y < dst.r0.y) || (fmin.y > dst.r1.y)) return false;
	if ((fmax.z < dst.r0.z) || (fmin.z > dst.r1.z)) return false;

	// At this point we find that bounding boxes overlap.
	// Technically that does not prove that the surface element is
	// inside the octree node, but any additional check would be
	// more expensive.

	return true;
}

// Find all surface elements that fall inside a node
void OTnode::FillNode(const vector<int>& parent_selist)
{
	FEOctree::Box self = { cmin, cmax };
	selist.reserve(parent_selist.size());
	// Loop over all surface elements in the parent node
	int nsize = (int)parent_selist.size();
	for (int i=0; i<nsize; ++i) {
		int j = parent_selist[i];

		if (CheckBoxIntersection(m_po->m_boxes[j], self)) {
			// add this surface element to the current node
			selist.push_back(j);
		}
	}
}

//-----------------------------------------------------------------------------
// Determine if a ray intersects any of the faces of this node.
// The ray originates at p and is directed along the unit vector n

bool OTnode::RayIntersectsNode(const vec3d& p, const vec3d& n)
{
	// check intersection with x-faces
	if (n.x) {
		// face passing through cmin
		double t = (cmin.x - p.x)/n.x;
		double y = p.y + t*n.y;
		double z = p.z + t*n.z;
		if ((y >= cmin.y) && (y <= cmax.y)
			&& (z >= cmin.z) && (z <= cmax.z))
			return true;
		// face passing through cmax
		t = (cmax.x - p.x)/n.x;
		y = p.y + t*n.y;
		z = p.z + t*n.z;
		if ((y >= cmin.y) && (y <= cmax.y)
			&& (z >= cmin.z) && (z <= cmax.z))
			return true;
	}
	// check intersection with y-faces
	if (n.y) {
		// face passing through cmin
		double t = (cmin.y - p.y)/n.y;
		double x = p.x + t*n.x;
		double z = p.z + t*n.z;
		if ((x >= cmin.x) && (x <= cmax.x)
			&& (z >= cmin.z) && (z <= cmax.z))
			return true;
		// face passing through cmax
		t = (cmax.y - p.y)/n.y;
		x = p.x + t*n.x;
		z = p.z + t*n.z;
		if ((x >= cmin.x) && (x <= cmax.x)
			&& (z >= cmin.z) && (z <= cmax.z))
			return true;
	}
	// check intersection with z-faces
	if (n.z) {
		// face passing through cmin
		double t = (cmin.z - p.z)/n.z;
		double x = p.x + t*n.x;
		double y = p.y + t*n.y;
		if ((x >= cmin.x) && (x <= cmax.x)
			&& (y >= cmin.y) && (y <= cmax.y))
			return true;
		// face passing through cmax
		t = (cmax.z - p.z)/n.z;
		x = p.x + t*n.x;
		y = p.y + t*n.y;
		if ((x >= cmin.x) && (x <= cmax.x)
			&& (y >= cmin.y) && (y <= cmax.y))
			return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Find intersected octree leaves and return a set of their surface elements
void OTnode::FindIntersectedLeaves(const vec3d& p, const vec3d& n, set<int>& sel, double srad)
{
	// Check if octree node is within search radius from p.
	bool bNodeWithinSRad = ( (cmin.x - srad <= p.x) && (cmax.x + srad >= p.x) &&
	                         (cmin.y - srad <= p.y) && (cmax.y + srad >= p.y) &&
	                         (cmin.z - srad <= p.z) && (cmax.z + srad >= p.z) );

	if (bNodeWithinSRad && RayIntersectsNode(p, n)) {
		int nc = (int)children.size();
		// if this node has children, search them for intersections
		if (nc) {
			for (int ic=0; ic<nc; ++ic) {
				children[ic].FindIntersectedLeaves(p, n, sel, srad);
			}
		}
		// otherwise we have reached the smallest intersected node in this
		// branch, return its surface element list
		else {
			// using a 'set' container avoids duplication of surface
			// elements shared by multiple octree nodes
			sel.insert(selist.begin(), selist.end());
		}
	}
}

//-----------------------------------------------------------------------------
// Print node content (for debugging purposes)

void OTnode::PrintNodeContent()
{
	int nel = (int)selist.size();
	printf("Level = %d\n", level);
	for (int i=0; i<nel; ++i)
		printf("%d\n",selist[i]);
	printf("-----------------------------------------------------\n");
	
	int nc = (int)children.size();
	for (int i=0; i<nc; ++i) {
		printf("Child = %d\n", i);
		children[i].PrintNodeContent();
	}
}

//-----------------------------------------------------------------------------
// Count nodes (for debugging purposes)

void OTnode::CountNodes(int& nnode, int& nlevel)
{
	int nc = (int)children.size();
	nnode += nc;
	nlevel = (level > nlevel) ? level : nlevel;
	for (int i=0; i<nc; ++i) {
		children[i].CountNodes(nnode, nlevel);
	}
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEOctree::FEOctree(FESurface* ps)
{
	m_ps = ps;
	max_level = 5;
	max_elem = 32;
	assert(max_level && max_elem);
}

FEOctree::~FEOctree()
{
	
}

//-----------------------------------------------------------------------------

void FEOctree::Init(const double stol)
{
	assert(m_ps);
	root.Clear();

	// calculate bounding boxes for all elements
	int NE = m_ps->Elements();
	FEMesh& mesh = *(m_ps->GetMesh());
	m_boxes.resize(NE);
#pragma omp parallel for
	for (int i = 0; i < NE; ++i)
	{
		FESurfaceElement& el = m_ps->Element(i);
		vec3d rn = mesh.Node(el.m_node[0]).m_rt;
		vec3d fmin = rn;
		vec3d fmax = rn;
		int N = el.Nodes();
		for (int i = 1; i < N; ++i) {
			rn = mesh.Node(el.m_node[i]).m_rt;
			if (rn.x < fmin.x) fmin.x = rn.x;
			if (rn.x > fmax.x) fmax.x = rn.x;
			if (rn.y < fmin.y) fmin.y = rn.y;
			if (rn.y > fmax.y) fmax.y = rn.y;
			if (rn.z < fmin.z) fmin.z = rn.z;
			if (rn.z > fmax.z) fmax.z = rn.z;
		}
		m_boxes[i] = { fmin, fmax };
	}
	
	// Set up the root node in the octree
	root.m_po = this;
	root.level = 0;
	
	// Create the list of all surface elements in the root node
	int nel = m_ps->Elements();
	root.selist.resize(nel);
	for (int i=0; i<nel; ++i)
		root.selist[i] = i;
	
	// Find the bounding box of the surface
	vec3d fenode = (m_ps->Node(0)).m_rt;
	root.cmin = fenode;
	root.cmax = fenode;
	for (int i=1; i<m_ps->Nodes(); ++i) {
		fenode = (m_ps->Node(i)).m_rt;
		if (fenode.x < root.cmin.x) root.cmin.x = fenode.x;
		else if (fenode.x > root.cmax.x) root.cmax.x = fenode.x;
		
		if (fenode.y < root.cmin.y) root.cmin.y = fenode.y;
		else if (fenode.y > root.cmax.y) root.cmax.y = fenode.y;
		
		if (fenode.z < root.cmin.z) root.cmin.z = fenode.z;
		else if (fenode.z > root.cmax.z) root.cmax.z = fenode.z;
	}
    
    // expand bounding box by search tolerance stol
    double d = (root.cmax - root.cmin).norm()*stol;
    root.cmin -= vec3d(d, d, d);
    root.cmax += vec3d(d, d, d);
	
	// Recursively create children of this root
	if (root.selist.size()) {
		if ((root.level < max_level) &&
			(root.selist.size() > max_elem))
			root.CreateChildren(max_level, max_elem);
	}

	int nodes = 0, levels = 0;
	root.CountNodes(nodes, levels);
	
	return;
}

void FEOctree::FindCandidateSurfaceElements(vec3d p, vec3d n, set<int>& sel, double srad)
{
	root.FindIntersectedLeaves(p, n, sel, srad);
}

void FEOctree::VisitIntersectedLeaves(const vec3d& p, const vec3d& n, double srad, const std::function<void(int elem)>& callback)
{
	std::stack<OTnode*> S;

	S.push(&root);

	while (!S.empty())
	{
		OTnode* node = S.top(); S.pop();

		vec3d cmin = node->cmin;
		vec3d cmax = node->cmax;

		// Check if octree node is within search radius from p.
		bool bNodeWithinSRad = ((cmin.x - srad <= p.x) && (cmax.x + srad >= p.x) &&
			(cmin.y - srad <= p.y) && (cmax.y + srad >= p.y) &&
			(cmin.z - srad <= p.z) && (cmax.z + srad >= p.z));

		if (bNodeWithinSRad && node->RayIntersectsNode(p, n)) {
			int nc = (int)node->children.size();
			// if this node has children, search them for intersections
			if (nc) {
				for (int ic = 0; ic < nc; ++ic) {
					S.push(&node->children[ic]);
				}
			}
			// otherwise we have reached the smallest intersected node in this
			// branch, return its surface element list
			else {
				// using a 'set' container avoids duplication of surface
				// elements shared by multiple octree nodes
				for (int i = 0; i < (int)node->selist.size(); ++i) {
					callback(node->selist[i]);
				}
			}
		}
	}
}
