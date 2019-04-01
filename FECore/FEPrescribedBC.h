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
class FEFacetSet;

//-----------------------------------------------------------------------------
// base class for prescribed boundary conditions
class FECORE_API FEPrescribedBC : public FEBoundaryCondition
{
public:
	FEPrescribedBC(FEModel* pfem);

public:
	// implement these functions

	// assign a node set to the prescribed BC
	virtual void AddNodes(const FENodeSet& set) {};

	// assign a surface to the BC
	// By default, the nodes of the surface are assigned to the BC
	virtual void AddNodes(const FEFacetSet& surf);

	// This function is called when the solver needs to know the 
	// prescribed dof values. The brel flag indicates wheter the total 
	// value is needed or the value with respect to the current nodal dof value
	virtual void PrepStep(std::vector<double>& ui, bool brel = true) = 0;

	// copy data from another class
	virtual void CopyFrom(FEPrescribedBC* pbc) = 0;
};
