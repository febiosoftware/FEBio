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
#include "FEBoundaryCondition.h"

//-----------------------------------------------------------------------------
class FENodeSet;

//-----------------------------------------------------------------------------
//! This class represents a fixed degree of freedom
//! This boundary conditions sets the BC attribute of the nodes in the nodeset
//! to DOF_FIXED when activated.
class FECORE_API FEFixedBC : public FEBoundaryCondition
{
public:
	//! constructors
	FEFixedBC(FEModel* pfem);
	FEFixedBC(FEModel* pfem, int node, int dof);

	//! add a node to the node set
	void AddNode(int node);

	//! add a node set
	void AddNodes(const FENodeSet& ns);

	//! set the degree of freedom that will be fixed
	void SetDOF(int dof);

	//! get the node list
	std::vector<int> GetNodeList();

	//! set the node list
	void SetNodeList(const std::vector<int>& nodeList);

public:
	//! serialization
	void Serialize(DumpStream& ar);

	//! activation
	void Activate();

	//! deactivations
	void Deactivate();

public:
	std::vector<int>	m_node;		//!< node set
	int					m_dof;		//!< fixed degree of freedom
};
