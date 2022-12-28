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
#include "FENodeNodeList.h"
#include <vector>

//-----------------------------------------------------------------------------
//! This class implements the idea of a level structure

//! A level structures groups the nodes per level. The basic notion is that 
//! (1) each node in level 1 is adjacent to a node in level 1 and/or level 2. 
//! (2) each node in level k is adjacent to a node in level k and/or level k-1
//! (3) each node in level 1<i<k is adjacent to a node in level i-1, i, or i+1
//! The width of the level structure is the largest width of any level.
//!
//! The level structure is used by the node - reorder algorithm implemented in
//! the FENodeReorder class.

class FECORE_API FELevelStructure
{
public:
	//! default constructor
	FELevelStructure();

	//! destructor
	virtual ~FELevelStructure();

	//! copy constructor
	FELevelStructure(FELevelStructure& L)
	{
		m_lval = L.m_lval;
		m_nref = L.m_nref;
		m_pl = L.m_pl;
		m_node = L.m_node;
		m_nwidth = L.m_nwidth;
		m_pNL = L.m_pNL;
	}

	//! assignment operator
	FELevelStructure& operator = (FELevelStructure& L)
	{
		m_lval = L.m_lval;
		m_nref = L.m_nref;
		m_pl = L.m_pl;
		m_node = L.m_node;
		m_nwidth = L.m_nwidth;
		m_pNL = L.m_pNL;

		return (*this);
	}

	//! Combine level structures L1 and L2 into one level structure
	void Merge(FELevelStructure& L1, FELevelStructure& L2, bool& bswap);

	//! Create a rooted level structure, starting at node nroot
	void Create(FENodeNodeList& L, int nroot);

	//! return the depth, ie. the number of levels.
	int Depth() const { return (int) m_lval.size(); }

	//! return the width of the level structure
	int Width() { return m_nwidth; }

	//! return the number of nodes in level l
	int Valence(int l) { return m_lval[l]; }

	//! return a list of nodes in level l
	int* NodeList(int l) { return &m_nref[0] + m_pl[l]; }

	//! return the level that node n is in
	int NodeLevel(int n) { return m_node[n]; }

	//! sort all nodes in levels l0 to l1 in order of increasing degree
	void SortLevels(int l0, int l1);

protected:
	std::vector<int>	m_lval;	//!< the level valence
	std::vector<int>	m_nref;	//!< the nodes in the level
	std::vector<int>	m_pl;	//!< start of each level
	std::vector<int>	m_node;	//!< the levels to which a nodes belongs

	FENodeNodeList*	m_pNL;	//!< The nodelist that generated the level structure

	int	m_nwidth;	//!< width of level structure

	static FELevelStructure*	m_pthis;	//!< this pointer used in static compare function

	//!< used in sorting the nodes of a level
	static int compare(const void*, const void*);
};
