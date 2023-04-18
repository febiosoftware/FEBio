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

#include "vec3d.h"
#include "vector.h"
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
	void FillNode(std::vector<int> parent_selist);
	bool ElementIntersectsNode(const int j);
	void PrintNodeContent();
	bool RayIntersectsNode(vec3d p, vec3d n);
	void FindIntersectedLeaves(vec3d p, vec3d n, std::set<int>& sel, double srad);
	void CountNodes(int& nnode, int& nlevel);
	
public:
	int				level;		//!< node level
	vec3d			cmin, cmax;	//!< node bounding box
	std::vector<int>		selist;		//!< list of surface elements inside this node
	std::vector<OTnode>	children;	//!< children of this node
	FESurface*		m_ps;		//!< the surface to search
};

//-----------------------------------------------------------------------------
//! This class is a helper class to find ray intersection with a surface

class FECORE_API FEOctree
{
	
public:
	FEOctree(FESurface* ps = 0);
	~FEOctree();
	
	//! attach to a surface
	void Attach(FESurface* ps) { m_ps = ps; }
	
	//! initialize search structures
	void Init(const double stol);
	
	//! find all candidate surface elements intersected by ray
	void FindCandidateSurfaceElements(vec3d p, vec3d n, std::set<int>& sel, double srad);
	
protected:
	FESurface*	m_ps;	//!< the surface to search
	OTnode root;		//!< root node in octree
	int max_level;		//!< maximum allowable number of levels in octree
	int max_elem;		//!< maximum allowable number of elements in any node
};
