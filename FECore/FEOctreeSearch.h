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

class FEElement;
class FEMesh;
class FEDomain;

//-----------------------------------------------------------------------------
//! This class is a helper class to find ray intersection with a surface

class FECORE_API FEOctreeSearch
{
	class Block
	{
	public:
		Block(FEMesh* mesh, int level);
		Block(Block* parent, int level);
		~Block();
		void Clear();
		void CreateChildren(int max_level, int max_elem);
		void FillBlock();
		bool ElementIntersectsNode(FEElement* pe);

		FEElement* FindElement(const vec3d& y, double r[3]);

		bool IsInside(const vec3d& r) const;

	public:
		int					m_level;		//!< node level
		vec3d				m_cmin, m_cmax;	//!< node bounding box
		std::vector<FEElement*>	m_selist;		//!< list of surface elements inside this node
		std::vector<Block*>		m_children;		//!< children of this node
		FEMesh*				m_mesh;			//!< the mesh to search
		Block*				m_parent;
	};

public:
	FEOctreeSearch(FEMesh* mesh);
	FEOctreeSearch(FEDomain* domain);
	~FEOctreeSearch();

	//! initialize search structures
	bool Init(double inflate = 0.005);

	FEElement* FindElement(const vec3d& x, double r[3]);

protected:
	FEMesh*		m_mesh;			//!< the mesh
	FEDomain*	m_dom;			//!< the domain to search (if null whole mesh will be searched)
	Block*		m_root;			//!< root node in octree
	int			m_max_level;	//!< maximum allowable number of levels in octree
	int			m_max_elem;		//!< maximum allowable number of elements in any node
};
