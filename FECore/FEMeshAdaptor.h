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
#include "FECoreBase.h"

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;
class FEElement;

//-----------------------------------------------------------------------------
// Base class for all mesh adaptors
class FECORE_API FEMeshAdaptor : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	FEMeshAdaptor(FEModel* fem);

	// The mesh adaptor should return true if the mesh remained unchanged
	// otherwise, it should return false.
	// iteration is the iteration number of the mesh adaptation loop
	virtual bool Apply(int iteration) = 0;
};

//-----------------------------------------------------------------------------
// This class is a helper class for use in the mesh adaptors. Its purpose is to select
// elements based on some criterion. This element list is then usually passed to the 
// mesh adaptor.
class FECORE_API FEMeshAdaptorCriterion : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	FEMeshAdaptorCriterion(FEModel* fem);

	void SetSort(bool b);

	void SetMaxElements(int m);

public:

	// return a list of elements that satisfy the criterion
	virtual std::vector<int> GetElementList();

	// This function needs to be overridden in order to select some elements
	// that satisfy the selection criterion
	// return true if the element satisfies the criterion, otherwise false
	// If this function returns true, the elemVal parameter should be set
	// This is used to sort the element list
	virtual bool Check(FEElement& el, double& elemVal);

private:
	bool	m_sortList;		// sort the list
	int		m_maxelem;		// the max nr of elements to return (or 0 if don't care)

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FECORE_API FEMaxVolumeCriterion : public FEMeshAdaptorCriterion
{
public:
	FEMaxVolumeCriterion(FEModel* fem);
	bool Check(FEElement& el, double& elemVal) override;
private:
	double	m_maxVolume;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FECORE_API FEMaxVariableCriterion : public FEMeshAdaptorCriterion
{
public:
	FEMaxVariableCriterion(FEModel* fem);
	bool Check(FEElement& el, double& elemVal) override;

private:
	double	m_maxValue;
	int		m_dof;

	DECLARE_FECORE_CLASS();
};
