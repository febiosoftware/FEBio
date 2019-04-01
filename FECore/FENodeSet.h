/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEItemList.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
// Forward declarations
class FEMesh;
class FENode;
class DumpStream;

//-----------------------------------------------------------------------------
//! Defines a node set
//
class FECORE_API FENodeSet : public FEItemList
{
public:
	FENodeSet();
	FENodeSet(FEMesh* pm);
	FENodeSet(const FENodeSet& ns);

	FEMesh* GetMesh() { return m_pmesh; }

	void operator = (const FENodeSet& ns);

	void create(int n);

	void add(int n);

	void add(const std::vector<int>& ns);

	void add(const FENodeSet& ns);

	int size() const { return (int)m_Node.size(); }

	int& operator [] (int i) { return m_Node[i]; }

	const int& operator [] (int i) const { return m_Node[i]; }

	void SetID(int n) { m_nID = n; }
	int GetID() const { return m_nID; }

	std::vector<int>& GetNodeList() { return m_Node; }

	FENode* Node(int i);
	const FENode* Node(int i) const;

public:
	void Serialize(DumpStream& ar);

	static void SaveClass(DumpStream& ar, FENodeSet* p);
	static FENodeSet* LoadClass(DumpStream& ar, FENodeSet* p);

protected:
	int					m_nID;		//!< ID of nodeset
	std::vector<int>	m_Node;		//!< list of nodes

protected:
	FEMesh*	m_pmesh;	//!< pointer to parent mesh
};
