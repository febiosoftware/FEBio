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
#include "fecore_api.h"
#include <vector>

class FESurface;
class FEMesh;
class FEElement;
class FEDomain;
class DumpStream;

//-----------------------------------------------------------------------------
//! The FENodeElemList class is a utility class that determines for each node 
//! to which element it belongs.

//! This class analyzes a mesh and finds for each node all elements that have
//! this node

class FECORE_API FENodeElemList
{
public:
	FENodeElemList(){}
	virtual ~FENodeElemList(){}

	//! build the node-element list for a surface
	void Create(const FESurface& s);

	//! build the node-selement list for a mesh
	void Create(FEMesh& mesh);

	//! build the node-element list for a domain
	void Create(FEDomain& dom);

	//! serialize data to/from dump file
	void Serialize(DumpStream& ar);

	//! Clear the list
	void Clear();

	int MaxValence();
	int Valence(int n) { return m_nval[n]; }
	FEElement** ElementList(int n) { return &m_eref[0] + m_pn[n]; }
	int* ElementIndexList(int n) { return &m_iref[0] + m_pn[n]; }

	int Size() { return (int) m_nval.size(); }

protected:
	std::vector<int>			m_nval;	// nodal valences
	std::vector<FEElement*>		m_eref;	// element pointers
	std::vector<int>			m_iref;	// element indices
	std::vector<int>			m_pn;	// start index into the eref array
};

//-----------------------------------------------------------------------------
//! Like the FEElemElemList, but can create multiple levels
class FENodeElemTree
{
public:
	FENodeElemTree() {}
	virtual ~FENodeElemTree() {}

	void Create(FESurface* ps, int k = 0);

	int Valence(int n) { return (int) m_nel[n].size(); }

	FEElement** ElementList(int n) { return &(m_nel[n][0]);}

	bool empty() { return m_nel.empty(); }

protected:
	std::vector< std::vector<FEElement*> >	m_nel;
};
