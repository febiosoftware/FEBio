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
#include <vector>
#include "FENodeDataMap.h"
#include "FENodeSet.h"
#include "FEDofList.h"
#include "FEModelParam.h"

//-----------------------------------------------------------------------------
//! Base class for defining initial conditions.
//! Initial conditions can be used to set the initial state of the model in an analysis. 
class FECORE_API FEInitialCondition : public FEStepComponent
{
	FECORE_SUPER_CLASS(FEIC_ID)
	FECORE_BASE_CLASS(FEInitialCondition);

public:
	FEInitialCondition(FEModel* pfem);
};

//-----------------------------------------------------------------------------
// Base class for initial conditions applied to node sets
class FECORE_API FENodalIC : public FEInitialCondition
{
public:
	FENodalIC(FEModel* fem);

	// set the nodeset for this component
	void SetNodeSet(FENodeSet* nset);

	// get the node set
	FENodeSet* GetNodeSet();

	// set the list of degrees of freedom
	void SetDOFList(const FEDofList& dofList);

	// serialization
	void Serialize(DumpStream& ar) override;

	// return the values for node i
	virtual void GetNodalValues(int inode, std::vector<double>& values) = 0;

public:
	void Activate() override;

	bool Init() override;

protected:
	FEDofList		m_dofs;
	FENodeSet*		m_nodeSet;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Class representing an initial condition on a degree of freedom
class FECORE_API FEInitialDOF : public FENodalIC
{
public:
	FEInitialDOF(FEModel* pfem);
	FEInitialDOF(FEModel* fem, int ndof, FENodeSet* nset);

	void SetDOF(int ndof);
	bool SetDOF(const char* szdof);

	void Serialize(DumpStream& ar) override;

	bool Init() override;

	// return the values for node i
	void GetNodalValues(int inode, std::vector<double>& values) override;

	void SetValue(double v);

protected:
	int				m_dof;		//!< degree of freedom
	FEParamDouble	m_data;		//!< nodal values

	DECLARE_FECORE_CLASS();
};
