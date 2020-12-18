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
#include "FECoreBase.h"
#include <functional>

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;
class FEElement;
class FEMaterialPoint;
class FEElementSet;

//-----------------------------------------------------------------------------
// Base class for all mesh adaptors
class FECORE_API FEMeshAdaptor : public FECoreBase
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

//-----------------------------------------------------------------------------
class FECORE_API FEMeshAdaptorSelection
{
public:
	struct Item {
		int		m_elementIndex;
		double	m_scaleFactor;
	};

public:
	FEMeshAdaptorSelection() {}
	FEMeshAdaptorSelection(size_t size) : m_itemList(size) {}

	void resize(size_t newSize) { m_itemList.resize(newSize); }
	Item& operator [] (size_t item) { return m_itemList[item]; }
	const Item& operator [] (size_t item) const { return m_itemList[item]; }
	bool empty() const { return m_itemList.empty(); }
	size_t size() const { return m_itemList.size(); }
	void push_back(int elemIndex, double scale) { m_itemList.push_back(Item{ elemIndex, scale }); }

private:
	std::vector<Item>	m_itemList;
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
	// The elements will be selected from the element set. If nullptr is passed
	// for the element set, the entire mesh will be processed
	virtual FEMeshAdaptorSelection GetElementSelection(FEElementSet* elset);

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

//-----------------------------------------------------------------------------
class FECORE_API FEElementSelectionCriterion : public FEMeshAdaptorCriterion
{
public:
	FEElementSelectionCriterion(FEModel* fem);
	FEMeshAdaptorSelection GetElementSelection(FEElementSet* elset) override;

private:
	vector<int>	m_elemList;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FECORE_API FEDomainErrorCriterion : public FEMeshAdaptorCriterion
{
public:
	FEDomainErrorCriterion(FEModel* fem);

	FEMeshAdaptorSelection GetElementSelection(FEElementSet* elset) override;

	// derived classes must implement this function
	virtual double GetMaterialPointValue(FEMaterialPoint& mp) = 0;

private:
	double	m_pct;

	DECLARE_FECORE_CLASS();
};

// helper function for projecting integration point data to nodes
void projectToNodes(FEMesh& mesh, std::vector<double>& nodeVals, std::function<double (FEMaterialPoint& mp)> f);
