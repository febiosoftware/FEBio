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
#include "FENodeSet.h"
using namespace std;

//-----------------------------------------------------------------------------
//! Base class for defining initial conditions.
//! Initial conditions can be used to set the initial state of the model in an analysis. 
class FECORE_API FEInitialCondition : public FEModelComponent
{
	FECORE_SUPER_CLASS

public:
	FEInitialCondition(FEModel* pfem);

	// set the nodeset for this component
	virtual void AddNodes(const FENodeSet& nset);

	// set the list of degrees of freedom
	void SetDOFList(const std::vector<int>& dofList);

	// serialization
	void Serialize(DumpStream& ar) override;

	// return the values for node i
	virtual void NodalValues(int inode, std::vector<double>& values) = 0;

public:
	void Activate() override;

protected:
	std::vector<int>	m_dofs;
	FENodeSet			m_nodeSet;
};

//-----------------------------------------------------------------------------
// Class representing an initial condition on a degree of freedom
class FECORE_API FEInitialDOF : public FEInitialCondition
{
public:
	FEInitialDOF(FEModel* pfem);

	void SetDOF(int ndof);

	void Serialize(DumpStream& ar) override;

	bool Init() override;

	void AddNodes(const FENodeSet& set) override;

	void Add(int node, double value);

	// return the values for node i
	void NodalValues(int inode, std::vector<double>& values) override;

private:
	int				m_dof;		//!< degree of freedom
	FENodeDataMap	m_data;		//!< nodal values

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Class for initializing degrees of freedom using a vec3d (useful for e.g. velocity)
class FECORE_API FEInitialBCVec3D : public FEInitialCondition
{
public:
	FEInitialBCVec3D(FEModel* pfem);

	void SetDOF(int d0, int d1, int d2);

	bool Init() override;

	void Serialize(DumpStream& ar) override;

	void AddNodes(const FENodeSet& nset) override;

	void Add(int nid, const vec3d& v);

	// return the values for node i
	void NodalValues(int inode, std::vector<double>& values) override;

private:
	int					m_dof[3];
	std::vector<vec3d>	m_data;
};
