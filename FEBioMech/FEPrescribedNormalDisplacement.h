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
#include <FECore/FEPrescribedBC.h>

class FEPrescribedNormalDisplacement : public FEPrescribedBC
{
	struct NODE
	{
		int		nodeId;		// node ID
		vec3d	normal;		// initial normal at node
	};

public:
	// constructor
	FEPrescribedNormalDisplacement(FEModel* fem);

	// initialization
	bool Init() override;

	// activation
	void Activate() override;

	// deactivation
	void Deactivate() override;

public:
	// assign a node set to the prescribed BC
	void AddNodes(const FEFacetSet& surf) override;

	// This function is called when the solver needs to know the 
	// prescribed dof values. The brel flag indicates wheter the total 
	// value is needed or the value with respect to the current nodal dof value
	void PrepStep(std::vector<double>& ui, bool brel = true) override;

	// This is called during nodal update and should be used to enforce the 
	// nodal degrees of freedoms
	void Update() override;

	// copy data from another class
	void CopyFrom(FEPrescribedBC* pbc) override;

private:
	vector<NODE>	m_node;

	double	m_scale;

	// hint parameter helps to identify the surface geometry
	// 0 : no hint (default)
	// 1 : sphere with center at origin
	int		m_hint;	//!< hint parameter helps to identify the surface geometry

	DECLARE_FECORE_CLASS();
};
