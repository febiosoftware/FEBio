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
#include "FEModelComponent.h"
#include "FENodeDataMap.h"

//-----------------------------------------------------------------------------
class FENodeSet;

//-----------------------------------------------------------------------------
//! Nodal load boundary condition
class FECORE_API FENodalLoad : public FEModelComponent
{
	FECORE_SUPER_CLASS

public:
	//! constructor
	FENodalLoad(FEModel* pfem);

	//! initialization
	bool Init() override;

	//! serialiation
	void Serialize(DumpStream& ar) override;

	//! Add a node to the node set
	void AddNode(int nid, double scale = 1.0);

	//! add a node set
	void AddNodes(const FENodeSet& ns, double scale = 1.0);

	//! number of nodes
	int Nodes() const { return (int)m_item.size(); }

	//! Node ID
	int NodeID(int n) const { return m_item[n]; }

	//! get nodal value
	double NodeValue(int n) const;

	//! get/set load 
	void SetLoad(double s);
	double GetLoad() const { return m_scale; }

	//! get/set degree of freedom
	void SetDOF(int ndof) { m_dof = ndof; }
	int GetDOF() const { return m_dof; }

private:
	int		m_dof;		// degree of freedom index

	double			m_scale;	// applied load scale factor
	vector<int>		m_item;		// item list
	FENodeDataMap	m_data;		// nodal data

	DECLARE_FECORE_CLASS();
};

