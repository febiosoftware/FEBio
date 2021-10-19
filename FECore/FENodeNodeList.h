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
#include <vector>
#include "fecore_api.h"

class FEMesh;
class FEDomain;

//-----------------------------------------------------------------------------
//! The FENodeNodeList class is a utility class that determines for each node 
//! the adjacent nodes

//! This class analyzes a mesh and finds for each node all nodes that are 
//! adjacent to this node

class FECORE_API FENodeNodeList
{
public:
	//! default constructor
	FENodeNodeList();

	//! desctructor
	virtual ~FENodeNodeList();

	//! create the node-node list for a mesh
	void Create(FEMesh& mesh);

	//! create the node-node list for a domain
	void Create(FEDomain& dom);

	int Size() const { return (int) m_nval.size(); }

	int Valence(int i) { return m_nval[i]; }
	int* NodeList(int i) { return &m_nref[0] + m_pn[i]; }

	void Sort();

protected:
	std::vector<int>	m_nval;	// nodal valences
	std::vector<int>	m_nref;	// adjacent nodes indices
	std::vector<int>	m_pn;	// start index into the nref array

	static FENodeNodeList*	m_pthis;
	static int compare(const void* e1, const void* e2);
};
