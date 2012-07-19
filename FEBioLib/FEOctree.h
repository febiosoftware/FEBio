#pragma once

#include "FECore/vec3d.h"
#include "NumCore/vector.h"
#include <set>

class FESurface;

//-----------------------------------------------------------------------------
//! This is a class for an octree node

class OTnode
{
public:
	OTnode() {}
	~OTnode() {}
	void Clear() {children.clear(); }
	void CreateChildren(const int max_level, const int max_elem);
	void FillNode(vector<int> parent_selist);
	bool ElementIntersectsNode(const int j);
	void PrintNodeContent();
	bool RayIntersectsNode(vec3d p, vec3d n);
	void FindIntersectedLeaves(vec3d p, vec3d n, std::set<int>& sel);
	void CountNodes(int& nnode);
	
public:
	int				level;		//!< node level
	vec3d			cmin, cmax;	//!< node bounding box
	vector<int>		selist;		//!< list of surface elements inside this node
	vector<OTnode>	children;	//!< children of this node
	FESurface*		m_ps;		//!< the surface to search
};

//-----------------------------------------------------------------------------
//! This class is a helper class to find ray intersection with a surface

class FEOctree  
{
	
public:
	FEOctree(FESurface* ps = 0);
	~FEOctree();
	
	//! attach to a surface
	void Attach(FESurface* ps) { m_ps = ps; }
	
	//! initialize search structures
	void Init();
	
	//! find all candidate surface elements intersected by ray
	void FindCandidateSurfaceElements(vec3d p, vec3d n, std::set<int>& sel);
	
protected:
	FESurface*	m_ps;	//!< the surface to search
	OTnode root;		//!< root node in octree
	int max_level;		//!< maximum allowable number of levels in octree
	int max_elem;		//!< maximum allowable number of elements in any node
};
