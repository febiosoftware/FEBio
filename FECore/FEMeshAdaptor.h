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
#include <functional>

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;
class FEElement;
class FEMaterialPoint;
class FEElementSet;
class FEMeshAdaptorCriterion;

//-----------------------------------------------------------------------------
// Base class for all mesh adaptors
class FECORE_API FEMeshAdaptor : public FEModelComponent
{
	FECORE_SUPER_CLASS

public:
	FEMeshAdaptor(FEModel* fem);

	void SetElementSet(FEElementSet* elemSet);
	FEElementSet* GetElementSet();

	// The mesh adaptor should return true if the mesh remained unchanged
	// otherwise, it should return false.
	// iteration is the iteration number of the mesh adaptation loop
	virtual bool Apply(int iteration) = 0;

private:
	FEElementSet*	m_elemSet;
};

// helper function for projecting integration point data to nodes
void FECORE_API projectToNodes(FEMesh& mesh, std::vector<double>& nodeVals, std::function<double (FEMaterialPoint& mp)> f);
void FECORE_API projectToNodes(FEDomain& dom, std::vector<double>& nodeVals, std::function<double(FEMaterialPoint& mp)> f);
