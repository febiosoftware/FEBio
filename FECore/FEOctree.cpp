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
#include "FEModel.h"

OTnode::OTnode()
{
	m_ps = nullptr;
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
	vec3d dc = (cmax - cmin)/2;
	for (int i=0; i<=1; ++i) {
		for (int j=0; j<=1; ++j) {
			for (int k=0; k<=1; ++k) {
				OTnode node;
				node.m_ps = m_ps;
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

//-----------------------------------------------------------------------------
// Find all surface elements that fall inside a node

void OTnode::FillNode(const vector<int>& parent_selist)
{
	// Loop over all surface elements in the parent node
	int nsize = (int)parent_selist.size();
	for (int i=0; i<nsize; ++i) {
		int j = parent_selist[i];
		if (ElementIntersectsNode(j)) {
			// add this surface element to the current node
			selist.push_back(j);
		}
	}
}

//-----------------------------------------------------------------------------
// Determine whether a surface element intersects a node

bool OTnode::ElementIntersectsNode(const int iel)
{
	// Extract FE node coordinates from surface element
	// and determine bounding box of surface element
	FEMesh& mesh = *(m_ps->GetMesh());
	FESurfaceElement& el = m_ps->Element(iel);
    vec3d rn = m_ps->GetNodalCoordinate(mesh.Node(el.m_node[0]));
	vec3d fmin = rn;
	vec3d fmax = rn;
	int N = el.Nodes();
	for (int i=1; i<N; ++i) {
        rn = m_ps->GetNodalCoordinate(mesh.Node(el.m_node[i]));
		if (rn.x < fmin.x) fmin.x = rn.x;
		if (rn.x > fmax.x) fmax.x = rn.x;
		if (rn.y < fmin.y) fmin.y = rn.y;
		if (rn.y > fmax.y) fmax.y = rn.y;
		if (rn.z < fmin.z) fmin.z = rn.z;
		if (rn.z > fmax.z) fmax.z = rn.z;
	}
	
	// Check if bounding boxes of OT node and surface element overlap
	if ((fmax.x < cmin.x) || (fmin.x > cmax.x)) return false;
	if ((fmax.y < cmin.y) || (fmin.y > cmax.y)) return false;
	if ((fmax.z < cmin.z) || (fmin.z > cmax.z)) return false;
	
	// At this point we find that bounding boxes overlap.
	// Technically that does not prove that the surface element is
	// inside the octree node, but any additional check would be
	// more expensive.
	
	return true;
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
	max_level = 6;
	max_elem = 9;
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
	
	// Set up the root node in the octree
	root.m_ps = m_ps;
	root.level = 0;
	
	// Create the list of all surface elements in the root node
	int nel = m_ps->Elements();
	root.selist.resize(nel);
	for (int i=0; i<nel; ++i)
		root.selist[i] = i;
	
	// Find the bounding box of the surface
    vec3d fenode = m_ps->GetNodalCoordinate(m_ps->Node(0));
	root.cmin = fenode;
	root.cmax = fenode;
	for (int i=1; i<m_ps->Nodes(); ++i) {
        fenode = m_ps->GetNodalCoordinate(m_ps->Node(i));
		if (fenode.x < root.cmin.x) root.cmin.x = fenode.x;
		if (fenode.x > root.cmax.x) root.cmax.x = fenode.x;
		if (fenode.y < root.cmin.y) root.cmin.y = fenode.y;
		if (fenode.y > root.cmax.y) root.cmax.y = fenode.y;
		if (fenode.z < root.cmin.z) root.cmin.z = fenode.z;
		if (fenode.z > root.cmax.z) root.cmax.z = fenode.z;
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
	
	return;
}

void FEOctree::FindCandidateSurfaceElements(vec3d p, vec3d n, set<int>& sel, double srad)
{
	root.FindIntersectedLeaves(p, n, sel, srad);
}
