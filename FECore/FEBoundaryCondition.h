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
#include "FENodeSet.h"

//-----------------------------------------------------------------------------
class FEFacetSet;

//-----------------------------------------------------------------------------
//! This class is the base class of boundary conditions.

//! Boundary conditions set the "bc" state of nodes. The bc-state determines
//! whether or not the dofs of the node will be assigned an equation number. 
//! Currently, there are two boundary conditions: a fixed (FEFixedBC) and a
//! prescribed (FEPrescribedBC) boundary condition. 
class FECORE_API FEBoundaryCondition : public FEModelComponent
{
	FECORE_SUPER_CLASS

public:
	//! constructor
	FEBoundaryCondition(FEModel* pfem);

	//! desctructor
	~FEBoundaryCondition();

	// set the dof list
	void SetDOFList(const std::vector<int>& dofs);

	// get the dof list
	const std::vector<int> GetDOFList();

	// get the node set
	const FENodeSet& GetNodeSet() const;

	//! Add a node to the node set of the bc
	virtual void AddNode(int node);

	//! will be overridden by derived classes
	virtual void AddNodes(const FENodeSet& nodeSet);

	// assign a surface to the BC
	// By default, the nodes of the surface are assigned to the BC
	virtual void AddNodes(const FEFacetSet& surf);

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! deactivate 
	void Deactivate() override;

	//! fill the prescribed values
	virtual void PrepStep(std::vector<double>& u, bool brel = true);

	// copy data from another class
	virtual void CopyFrom(FEBoundaryCondition* pbc) = 0;

protected:
	std::vector<int>	m_dofs;		//!< dof list
	FENodeSet			m_nodeSet;	//!< the node set for which this BC is defined.
};
