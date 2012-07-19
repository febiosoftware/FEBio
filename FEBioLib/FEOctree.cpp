#include "stdafx.h"
#include "FEOctree.h"
#include "FECore/FESurface.h"

//-----------------------------------------------------------------------------
// Create the eight children of an octree node and find their contents

void OTnode::CreateChildren(const int max_level, const int max_elem)
{
	int i,j,k;
	vec3d dc = (cmax - cmin)/2;
	for (i=0; i<2; ++i) {
		for (j=0; j<2; ++j) {
			for (k=0; k<2; ++k) {
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
					if ((int)node.selist.size() > max_elem) {
						if (max_level) {
							if (node.level < max_level) {
								node.CreateChildren(max_level, max_elem);
							}
						} else {
							node.CreateChildren(max_level, max_elem);
						}

					}
					// store this node
					children.push_back(node);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Find all surface elements that fall inside a node

void OTnode::FillNode(vector<int> parent_selist)
{
	int i, j;
	// Loop over all surface elements in the parent node
	for (i=0; i<(int)parent_selist.size(); ++i) {
		j = parent_selist[i];
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
	int i;

	// Extract FE node coordinates from surface element
	// and determine bounding box of surface element
	FEMesh& mesh = *(m_ps->GetMesh());
	FESurfaceElement& el = m_ps->Element(iel);
	int N = el.Nodes();
	vector<vec3d> fenode(N);
	fenode[0] = mesh.Node(el.m_node[0]).m_rt;
	vec3d fmin = fenode[0];
	vec3d fmax = fenode[0];
	for (i=1; i<N; ++i) {
		fenode[i] = mesh.Node(el.m_node[i]).m_rt;
		if (fenode[i].x < fmin.x) fmin.x = fenode[i].x;
		if (fenode[i].x > fmax.x) fmax.x = fenode[i].x;
		if (fenode[i].y < fmin.y) fmin.y = fenode[i].y;
		if (fenode[i].y > fmax.y) fmax.y = fenode[i].y;
		if (fenode[i].z < fmin.z) fmin.z = fenode[i].z;
		if (fenode[i].z > fmax.z) fmax.z = fenode[i].z;
	}
	
	// Check if bounding boxes of OT node and surface element overlap
	if ((fmin.x < cmin.x) && (fmax.x < cmin.x)) return false;
	if ((fmin.y < cmin.y) && (fmax.y < cmin.y)) return false;
	if ((fmin.z < cmin.z) && (fmax.z < cmin.z)) return false;
	if ((fmin.x > cmax.x) && (fmax.x > cmax.x)) return false;
	if ((fmin.y > cmax.y) && (fmax.y > cmax.y)) return false;
	if ((fmin.z > cmax.z) && (fmax.z > cmax.z)) return false;
	
	// At this point we find that bounding boxes overlap.
	// Technically that does not prove that the surface element is
	// inside the octree node, but any additional check would be
	// more expensive.
	
	return true;
}

//-----------------------------------------------------------------------------
// Determine if a ray intersects any of the faces of this node.
// The ray originates at p and is directed along the unit vector n

bool OTnode::RayIntersectsNode(vec3d p, vec3d n)
{
	double t, x, y, z;
	
	// check intersection with x-faces
	if (n.x) {
		// face passing through cmin
		t = (cmin.x - p.x)/n.x;
		y = p.y + t*n.y;
		z = p.z + t*n.z;
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
		t = (cmin.y - p.y)/n.y;
		x = p.x + t*n.x;
		z = p.z + t*n.z;
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
		t = (cmin.z - p.z)/n.z;
		x = p.x + t*n.x;
		y = p.y + t*n.y;
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

void OTnode::FindIntersectedLeaves(vec3d p, vec3d n, set<int>& sel)
{
	if (RayIntersectsNode(p, n)) {
		int nc = children.size();
		// if this node has children, search them for intersections
		if (nc) {
			for (int ic=0; ic<nc; ++ic) {
				children[ic].FindIntersectedLeaves(p, n, sel);
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
	int i;
	int nel = selist.size();
	printf("Level = %d\n", level);
	for (i=0; i<nel; ++i)
		printf("%d\n",selist[i]);
	printf("-----------------------------------------------------\n");
	
	int nc = children.size();
	for (i=0; i<nc; ++i) {
		printf("Child = %d\n", i);
		children[i].PrintNodeContent();
	}
}

//-----------------------------------------------------------------------------
// Count nodes (for debugging purposes)

void OTnode::CountNodes(int& nnode)
{
	int nc = children.size();
	nnode += nc;
	for (int i=0; i<nc; ++i) {
		children[i].CountNodes(nnode);
	}
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEOctree::FEOctree(FESurface* ps)
{
	m_ps = ps;
	max_level = 0;	// by default there is no limit to number of levels
	max_elem = 9;
}

FEOctree::~FEOctree()
{
	
}

//-----------------------------------------------------------------------------

void FEOctree::Init()
{
	assert(m_ps);
	int i;
	root.Clear();
	
	// Set up the root node in the octree
	root.m_ps = m_ps;
	root.level = 0;
	
	// Create the list of all surface elements in the root node
	int nel = m_ps->Elements();
	root.selist.resize(nel);
	for (i=0; i<nel; ++i)
		root.selist[i] = i;
	
	// Find the bounding box of the surface
	FEMesh& mesh = *(m_ps->GetMesh());
	vec3d fenode = mesh.Node(m_ps->m_node[0]).m_rt;
	root.cmin = fenode;
	root.cmax = fenode;
	for (i=1; i<(int)m_ps->m_node.size(); ++i) {
		fenode = mesh.Node(m_ps->m_node[i]).m_rt;
		if (fenode.x < root.cmin.x) root.cmin.x = fenode.x;
		if (fenode.x > root.cmax.x) root.cmax.x = fenode.x;
		if (fenode.y < root.cmin.y) root.cmin.y = fenode.y;
		if (fenode.y > root.cmax.y) root.cmax.y = fenode.y;
		if (fenode.z < root.cmin.z) root.cmin.z = fenode.z;
		if (fenode.z > root.cmax.z) root.cmax.z = fenode.z;
	}
	
	// Recursively create children of this root
	if ((int)root.selist.size() > max_elem) {
		if (max_level) {
			if (root.level < max_level) {
				root.CreateChildren(max_level, max_elem);
			}
		} else {
			root.CreateChildren(max_level, max_elem);
		}
	}
	
// For debugging purposes...
//	root.PrintNodeContent();
//	int nnode = 0;
//	root.CountNodes(nnode);
//	printf("Number of nodes in octree = %d\n",nnode);
	
	return;
}

void FEOctree::FindCandidateSurfaceElements(vec3d p, vec3d n, set<int>& sel)
{
	root.FindIntersectedLeaves(p, n, sel);
}