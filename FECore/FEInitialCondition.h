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
#include <vector>
#include "FENodeDataMap.h"
using namespace std;

class FENodeSet;

//-----------------------------------------------------------------------------
//! Base class for defining initial conditions.
//! Initial conditions can be used to set the initial state of the model in an analysis. 
class FECORE_API FEInitialCondition : public FEModelComponent
{
	FECORE_SUPER_CLASS

public:
	FEInitialCondition(FEModel* pfem);
};

//-----------------------------------------------------------------------------
// Class representing an initial condition on a degree of freedom
class FECORE_API FEInitialBC : public FEInitialCondition
{
public:
	FEInitialBC(FEModel* pfem);

	void SetDOF(int ndof) { m_dof = ndof; }

	void Serialize(DumpStream& ar) override;

	void Activate() override;

	void SetNodes(const FENodeSet& set);

	void Add(int node, double value);

public:
	int				m_dof;		//!< degree of freedom
	vector<int>		m_item;		//!< node IDs
	FENodeDataMap	m_data;		//!< nodal values

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Class for initializing degrees of freedom using a vec3d (useful for e.g. velocity)
class FECORE_API FEInitialBCVec3D : public FEInitialCondition
{
	struct ITEM
	{
		int		nid;	//!< node ID
		vec3d	v0;		//!< initial value

		void Serialize(DumpStream& ar);
	};

public:
	FEInitialBCVec3D(FEModel* pfem) : FEInitialCondition(pfem) { m_dof[0] = m_dof[1] = m_dof[2] = -1; }

	void SetDOF(int d0, int d1, int d2) { m_dof[0] = d0; m_dof[1] = d1; m_dof[2] = d2; }

	void Serialize(DumpStream& ar);

	void Activate();

	void Add(int nid, vec3d v) { ITEM it = {nid, v}; m_item.push_back(it); }

public:
	vector<ITEM>	m_item;
	int				m_dof[3];
};
