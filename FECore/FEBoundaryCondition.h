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
#include "FEStepComponent.h"
#include "FENodeSet.h"
#include "FEDofList.h"

//-----------------------------------------------------------------------------
class FEFacetSet;

//-----------------------------------------------------------------------------
//! This class is the base class of boundary conditions.

//! Boundary conditions set the "bc" state of nodes. The bc-state determines
//! whether or not the dofs of the node will be assigned an equation number. 
//! Currently, there are two boundary conditions: a fixed (FEFixedBC) and a
//! prescribed (FEPrescribedBC) boundary condition. 
class FECORE_API FEBoundaryCondition : public FEStepComponent
{
	FECORE_SUPER_CLASS(FEBC_ID)
	FECORE_BASE_CLASS(FEBoundaryCondition);

public:
	//! constructor
	FEBoundaryCondition(FEModel* pfem);

	//! desctructor
	~FEBoundaryCondition();

	//! fill the prescribed values
	virtual void PrepStep(std::vector<double>& u, bool brel = true);

	// copy data from another class
	virtual void CopyFrom(FEBoundaryCondition* pbc) = 0;
    
    // repair BC if needed
    virtual void Repair() {}

	void Serialize(DumpStream& ar) override;

	// TODO: Temporary construction to update some special boundary conditions in FEModel::Update
	//       Will likely remove this at some point.
	virtual void UpdateModel() {}


public:
	// set the dof list
	void SetDOFList(int ndof);
	void SetDOFList(const std::vector<int>& dofs);
	void SetDOFList(const FEDofList& dofs);

	const FEDofList& GetDofList() const { return m_dof; }

protected:
	FEDofList	m_dof;	// the dof list for the BC
};
