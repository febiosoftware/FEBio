/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEModelComponent.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! linear constraint
class FECORE_API FELinearConstraint : public FEModelComponent
{
public:
	class DOF
	{
	public:
		DOF()
		{
			node = dof = -1;
			val = 0.0;
		}

	public:
		int		node;	// node number
		int		dof;	// degree of freedom
		double	val;	// coefficient value (ignored for parent dof)
	};

public:
	// constructors
	FELinearConstraint();
	FELinearConstraint(FEModel* pfem);
	FELinearConstraint(const FELinearConstraint& LC);

	// copy data
	void CopyFrom(const FELinearConstraint& LC);

	// serialization
	void Serialize(DumpStream& ar);

	// initialize the linear constraint
	bool Init();

	// make the constraint active
	void Activate();
	void Deactivate();

	// set the parent degree of freedom
	void SetParentDof(int dof, int node);

	// add a child degree of freedom
	void AddChildDof(int dof, int node, double v);

	// set the linear constraint offset
	void SetOffset(double d) { m_off = d; }

public:
	DOF			m_parentDof;	// parent degree of freedom
	vector<DOF>	m_childDof;		// list of child dofs
	double		m_off;			// offset value
};
