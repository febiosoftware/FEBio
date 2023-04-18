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
#include "fecore_api.h"
#include <vector>

//-----------------------------------------------------------------------------
class FEModel;
class DumpStream;

//-----------------------------------------------------------------------------
// Convenience class for creating a list of degrees of freedom, without the 
// need to go through the FEModel class;

class FECORE_API FEDofList
{
public:
	FEDofList(FEModel* fem);
	FEDofList(const FEDofList& dofs);

	// assignment operator
	void operator = (const FEDofList& dofs);
	void operator = (const std::vector<int>& dofs);

	// clear the list
	void Clear();

	// Add a degree of freedom
	bool AddDof(const char* szdof);

	// Add a degree of freedom
	bool AddDof(int ndof);

	// Add all the dofs of a variable
	bool AddVariable(const char* szvar);

	// Add all the dofs a variable
	bool AddVariable(int nvar);

	// Add degrees of freedom
	bool AddDofs(const FEDofList& dofs);

	// is the list empty
	bool IsEmpty() const;

	// number of dofs in the list
	int Size() const;

	// access a dof
	int operator [] (int n) const;

	// serialization
	void Serialize(DumpStream& ar);

	// see if this dof list contains all the dofs of a FEDofList
	bool Contains(int dof);
	bool Contains(const FEDofList& dof);

	// Get the interpolation order 
	int InterpolationOrder(int index) const;

private:
	FEModel*			m_fem;
	std::vector<int>	m_dofList;
};
